import xarray as xr
import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon
from scipy.interpolate import griddata
from geopandas.tools import sjoin
from typing import Union
from tools import timing
import dask.array as da
import glob
from dask.diagnostics import ProgressBar

def compute_scale_and_offset(
        array,
        n=16,
        num_reserveed_tail=0
    ):
    try:
        #DaskArray
        vmin = float(array.min(skipna=True).compute())
        vmax = float(array.max(skipna=True).compute())
    except:
        #Numpy Array
        vmin = array.min()
        vmax = array.max()


    # stretch/compress data to the available packed range
    nbins = 2 ** n - 1- num_reserveed_tail
    scale_factor = (vmax - vmin) / nbins

    # translate the range to be symmetric about zero
    last_bit_value = (2**(n - 1) - num_reserveed_tail )
    add_offset = vmin + last_bit_value * scale_factor

    return scale_factor, add_offset

@xr.register_dataset_accessor('d3dfm')
class D3DFM_Dataset_Accessor():
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        self._mesh2d_face_verts = None
        self._mesh2d_face_pt = None
        self._decoded_time = None
        self._land_masking = {}

        self.ugrid = {}
        self.ugrid['node_coordinates'] = None
        self.ugrid['face_node_connectivity'] = None
        self.ugrid['face_node_start_index'] = None
        self.ugrid['face_coordinates'] = None
        self._parse_ugrid_definition()

    def _parse_ugrid_definition(self):
        data_vars = self._obj.data_vars
        #Dummy variable to store topology information
        def_var = None
        if 'mesh2d' in data_vars:
            def_var = 'mesh2d'
        elif 'mesh2' in data_vars:
            def_var = 'mesh2'

        if def_var is not None:
            #Required attributes by UGRID 1.0
            self.ugrid['node_coordinates'] = self.get_var(def_var).attrs['node_coordinates'].split()
            self.ugrid['face_node_connectivity'] = self.get_var(def_var).attrs['face_node_connectivity']
            self.ugrid['face_node_start_index'] = self.get_var(self.ugrid['face_node_connectivity']).attrs.get('start_index')

            #Optional attributes by UGRID 1.0
            self.ugrid['face_coordinates'] = self.get_var(def_var).attrs.get('face_coordinates', None)
            if self.ugrid['face_coordinates'] is not None:
                self.ugrid['face_coordinates'] = self.ugrid['face_coordinates'].split()

    def get_var(
            self,
            var: str,
        ):
        """
        Return the data of `varname` of a dataset.

        :param var:
        :return:
        """

        return self._obj[var]

    @property
    def mesh2d_face_verts(self):
        """
        Get the vertices of every face, defined by the attributes `node_coordinates` and `face_node_connectivity`.
        For example, triangular mesh will have three vertices.

        :return:
        """
        if self._mesh2d_face_verts is None:
            node_x_var, node_y_var = self.ugrid['node_coordinates']
            face_node_con_var = self.ugrid['face_node_connectivity']
            face_node_start_index = self.ugrid['face_node_start_index']
            node_x = self.get_var(node_x_var).data
            node_y = self.get_var(node_y_var).data
            face_nodes = (self.get_var(face_node_con_var).data - face_node_start_index)
            face_nodes[np.isnan(face_nodes)] = -999
            face_nodes = face_nodes.astype(int)
            node_xy = np.c_[node_x, node_y]
            verts = node_xy[face_nodes]
            verts[face_nodes==-999, :] = np.nan
            self._mesh2d_face_verts = verts
        return self._mesh2d_face_verts

    @property
    def mesh2d_face_pt(self):
        """
        Get the representaitive point (usually the midpoint) of every face, defined by an optional attribute
        `face_coordinates`.  If not defined, return the calculated midpoint of every face.

        :return:
        """
        if self._mesh2d_face_pt is None:
            if self.ugrid['face_coordinates'] is not None:
                x_var, y_var = self.ugrid['face_coordinates']
                x = self.get_var(x_var).data
                y = self.get_var(y_var).data
                self._mesh2d_face_pt = np.stack([x, y]).T
            else:
                self._mesh2d_face_pt = np.nanmean(self.mesh2d_face_verts, axis=1)

        return self._mesh2d_face_pt


    @property
    def decoded_time(self):
        if self._decoded_time is None:
            ds = xr.Dataset(
                coords=dict(
                    time=self.get_var('time')
                )
            )
            ds=xr.decode_cf(ds)
            self._decoded_time = ds.time.data
        return self._decoded_time


    def __mesh2d_face_to_gdf(self):
        polygons = []
        for face in self.mesh2d_face_verts:
            nan_mask = ~np.isnan(face).any(axis=1)
            polygons.append(Polygon(face[nan_mask, :]))
        gdf = gpd.GeoDataFrame(geometry=polygons)

        return gdf

    def __mesh2d_face_union(self):
        return self.__mesh2d_face_to_gdf().geometry.unary_union

    def __mesh2d_face_union_to_gdf(self):
        union = self.__mesh2d_face_union()
        if isinstance(union, Polygon):
            union_gdf = gpd.GeoDataFrame(geometry=[union])
        else:
            union_gdf = gpd.GeoDataFrame(geometry=list(union.geoms))
        return union_gdf


    def mesh2d_face_to_file(
            self,
            path: str,
            driver: str="ESRI Shapefile"
        ):
        """
        Export the polygons of every face of the flexible mesh to a file (shapefile by default).

        :param path:
        :param driver:
        :return:
        """
        gdf = self.__mesh2d_face_to_gdf()
        gdf.to_file(path, driver=driver)

    def mesh2d_face_union_to_file(
            self,
            path: str,
            driver: str="ESRI Shapefile"
        ):
        """
        Export the union polygon of all faces of the flexible mesh to a file (shapefile by default).

        :param path:
        :param driver:
        :return:
        """
        union = self.__mesh2d_face_union()
        union_gdf = gpd.GeoDataFrame(geometry=list(union.geoms))
        union_gdf.to_file(path, driver=driver)

    def __get_land_masking_indices(
            self,
            bbox: list[float],
            ncellx: int,
            ncelly: int,
    ):
        bbox = tuple(bbox)
        if self._land_masking.get((bbox, ncellx, ncelly)) is None:
            left, bottom, right, top = bbox
            xc = np.linspace(left, right, ncellx)
            yc = np.linspace(bottom, top, ncelly)
            xpxl, ypxl = np.meshgrid(xc, yc)
            union_poly = self.__mesh2d_face_union_to_gdf()
            raster_points = gpd.GeoDataFrame(geometry=gpd.points_from_xy(xpxl.ravel(), ypxl.ravel()))
            pointInPolys = sjoin(raster_points, union_poly, how='left')
            nan_ind = pointInPolys.loc[pointInPolys.index_right.isna(), :].index
            jpxl, ipxl = np.meshgrid(np.arange(ncellx), np.arange(ncelly))
            ij = np.c_[ipxl.ravel(), jpxl.ravel()]
            nan_ij = ij[nan_ind, :]
            nan_i, nan_j = nan_ij[:, 0], nan_ij[:, 1]
            self._land_masking[(bbox, ncellx, ncelly)] = [nan_i, nan_j]
            print('get_land_masking_indices executed')
        return self._land_masking.get((bbox, ncellx, ncelly))

    def rasterize_variable(
            self,
            varname:str,
            timesteps: Union[int,list[int]]=None,
            layers: Union[int,list[int]]=None,
            bbox: list[float]=None,
            ncellx: int=1000,
            ncelly: int=1000,
            method: str = "nearest",
            maskland: bool = True,
            missing_value: int = np.nan,
    ):

        if type(timesteps) is int:
            timesteps = [timesteps]
        if type(layers) is int:
            layers = [layers]
        var = self.get_var(varname)
        var_coordinates = var.d3dfm.coordinates
        if (var_coordinates.get('lon') is None) or (var_coordinates.get('lat') is None):
            raise ValueError(f"'{varname}' has no latitude or longitude.")
        if var_coordinates.get('time') is None:
            if timesteps is not None:
                raise ValueError(f"'{varname}' does not have time coordinates but timesteps is not None.")
        if var_coordinates.get('layer') is None:
            if layers is not None:
                raise ValueError(f"'{varname}' does not have layer coordinates but layers is not None.")

        # #Check whether the `ds.some_var.data` is dask array. Fast. Without loading data into memory.
        isDaskArray = var.chunks is not None

        #Get the x- and y- coordinates of the flexible mesh, each in 1D array.
        flexMesh_x_varname, flexMesh_y_varname = var_coordinates['lon'], var_coordinates['lat']
        flexMesh_x = var[flexMesh_x_varname].data
        flexMesh_y = var[flexMesh_y_varname].data

        #Generate x- and y- coordiantes of a regular grid, each in 2D array.
        left, bottom, right, top = bbox
        xc = np.linspace(left, right, ncellx)
        yc = np.linspace(bottom, top, ncelly)
        regGrid_x, regGrid_y = np.meshgrid(xc, yc)  # x pixel, y pixel

        #Indices for nearest interpolation
        flexMesh_indices = np.arange(flexMesh_x.shape[0])
        regGrid_indices = griddata((flexMesh_x, flexMesh_y), flexMesh_indices, (regGrid_x, regGrid_y), method=method)
        regGrid_shape = regGrid_indices.shape
        regGrid_indices_1D = regGrid_indices.flatten()


        #Get the indices of raster xy which are on land (land masking indices)
        if maskland:
            nan_i, nan_j = self.__get_land_masking_indices(bbox, ncellx, ncelly)
            nan_1Di = nan_i * ncellx + nan_j

        ntimesteps = len(timesteps)
        nlayers = len(layers)
        timesteps = np.array(timesteps)
        layers = np.array(layers)

        regGrid_var = var.d3dfm.ifilter(times=timesteps, layers=layers, faces=regGrid_indices_1D).data

        if maskland:
            regGrid_var[:, nan_1Di, :] = missing_value
        regGrid_var = regGrid_var.reshape( (ntimesteps,) + regGrid_shape + (nlayers,) )
        regGrid_var = da.transpose(regGrid_var, (0, 3, 1, 2))



        lat = regGrid_y[:, 0]
        lon = regGrid_x[0, :]
        time = var['time'].d3dfm.ifilter(times=timesteps)
        dataArray_var = xr.DataArray(
            data=regGrid_var,
            dims=['time', 'layer', 'lat', 'lon'],
            coords={
                'lon': (
                    ['lon'],
                    lon,
                    {
                        "standard_name": "longitude",
                        "long_name": "longitude",
                        "units": "degrees_east",
                        "axis": "X",
                        "_FillValue": -999,
                    }
                ),
                'lat': (
                    ['lat'],
                    lat,
                    {
                        "standard_name": "latitude",
                        "long_name": "latitude",
                        "units": "degrees_north",
                        "axis": "Y",
                        "_FillValue": -999,
                    }
                ),
                'layer': (
                    ['layer'],
                    layers,
                    {
                        "standard_name": "layer",
                        # "units": "up", #this will cause error in Delft3d QuickPLot
                        "positive" : "up",
                        # By convention of Delft3D Flow or D-Flow FM, Layer 1 reprsents the bottom layer and Layer 1+
                        # represents the layer closer to the surface. So "positive": "up" should be used.
                        "axis": "Z",
                        "_FillValue": -999,
                    }
                ),
                'time': (
                    ['time'],
                    time.data,
                    {
                        "standard_name": "time",
                        "long_name": "time",
                        "units": time.attrs['units'],
                        "calendar": "proleptic_gregorian",
                        "axis": "T",
                        "_FillValue": -999,
                    }
                ),
            },
        )
        return dataArray_var


    def rasterize_variables(
            self,
            variables: Union[str, list[str]],
            times: Union[int, list[int]],
            layers: Union[int,list[int]],
            bbox: list[float] = [113.213111, 21.917770, 114.627601, 23.145613],
            ncellx: int = 1000,
            ncelly: int = 1000,
            method: str = "nearest",
            maskland: bool = True,
            missing_value: int = np.nan,
    ):
        if type(times) is int:
            times = [times]
        if type(layers) is int:
            layers = [layers]
        if type(variables) is str:
            variables = [variables]

        data_vars = {}
        for variable in variables:
            raster_variable = self.rasterize_variable(
                variable,
                times,
                layers,
                bbox,
                ncellx,
                ncelly,
                method=method,
                maskland=maskland,
                missing_value=missing_value,
            )
            data_vars[variable] = raster_variable

        ds = xr.Dataset(
            data_vars=data_vars,
            attrs={
                'Conventions': "CF-1.8",
                "Institution": "EPD",
                "Author": "Ken TSUI",
            },
        )
        return ds
        # ds.d3dfm.to_packed_netcdf(path, packing={variable:dtype}, mode=mode)


    def to_packed_netcdf(
            self,
            path,
            packing: dict = None,
            mode: str = "w",
            encoding = None,
            **kwargs
    ):
        """
        Export packed netcdf.
        :param path:
        :param packing:
        :param mode:
        :param encoding:
        :return:
        """
        if packing is not None:
            packcoding={}
            for variable, dtype in packing.items():
                dtype = dtype.lower().strip()

                if dtype == "none":
                    packcoding[variable] = {
                            "_FillValue": -999,
                    }
                    continue

                dtype_to_nbits = {
                    "int64": {'nbits': 64, '_FillValue': -9223372036854775808},
                    "int32": {'nbits': 32, '_FillValue': -2147483647},
                    "int16": {'nbits': 16, '_FillValue': -32768},
                    "int8": {'nbits': 8, '_FillValue': -128}
                }

                if dtype in dtype_to_nbits:
                    nbits = dtype_to_nbits[dtype]['nbits']
                    _FillValue = dtype_to_nbits[dtype]['_FillValue']

                    scale_factor, add_offset = compute_scale_and_offset(
                        self.get_var(variable),
                        n=nbits,
                        num_reserveed_tail=2
                    )


                    packcoding[variable] = {
                            "dtype": dtype,
                            "scale_factor": scale_factor,
                            "add_offset": add_offset,
                            "_FillValue": _FillValue,
                    }

                else:
                    raise ValueError(f"Data type: {dtype} is not supported. Only accept 'int8', 'int16', 'int32', and 'int64'.")

            if len(packcoding)!=0:
                encoding = packcoding

        return self._obj.to_netcdf(path, mode=mode, format='NETCDF4', encoding=encoding, **kwargs)

@xr.register_dataarray_accessor('d3dfm')
class D3DFM_DataArray_Accessor:
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        self._coordinates = None

    def __get_coord_vars(self):
        pass

    def ifilter(
            self,
            times: Union[int, list[int]] = None,
            layers: Union[int, list[int]] = None,
            faces: Union[int, list[int]] = None,
            edges: Union[int, list[int]] = None,
            nodes: Union[int, list[int]] = None,
            time_var: str = None, # 'time'
            layer_var: str = None, #'mesh2d_nLayers'
            face_var: str = None,  # 'mesh2d_nFaces'
            edge_var: str = None,  # 'mesh2d_nEdges'
            node_var: str = None,  # 'mesh2d_nNodes'
        ):
        if time_var is None:
            time_var = self.coordinates.get('time')
        if layer_var is None:
            layer_var = self.coordinates.get('layer')
        if face_var is None:
            face_var = self.coordinates.get('face')
        if edge_var is None:
            edge_var = self.coordinates.get('edge')
        if node_var is None:
            node_var = self.coordinates.get('node')

        filter = {}
        if times is not None:
            filter[time_var] = times
        if layers is not None:
            filter[layer_var] = layers
        if faces is not None:
            filter[face_var] = faces
        if edges is not None:
            filter[edge_var] = edges
        if nodes is not None:
            filter[node_var] = nodes

        return self._obj.isel(**filter)

    def bboxfilter(
            self,
            bbox: Union[list[float]] = None,
            x_var: str = None, #'mesh2d_face_x',
            y_var: str = None, #'mesh2d_face_y',
            face_var: str = None, # 'mesh2d_nFaces',
    ):
        if x_var is None:
            x_var = self.coordinates.get('lon')
        if y_var is None:
            y_var = self.coordinates.get('lat')
        if face_var is None:
            face_var = self.coordinates.get('face')

        filter = {}
        if bbox is not None:
            left, bottom, right, top = bbox
            x = np.logical_and(
                self._obj[face_var][x_var]<=right,
                self._obj[face_var][x_var]>=left,
            )
            y = np.logical_and(
                self._obj[face_var][y_var]<=top,
                self._obj[face_var][y_var]>=bottom,
            )
            xy = np.logical_and(x, y)
            filter[face_var] = self._obj[face_var][xy]

        return self._obj.isel(**filter)

    def mean(
            self,
            times: Union[int, list[int]]=None,
            layers: Union[int, list[int]]=None,
            bbox: Union[list[int]] = None,
            time_var: str = None, #'time',
            layer_var: str = None, #'mesh2d_nLayers',
            x_var: str = None, #'mesh2d_face_x',
            y_var: str = None, #'mesh2d_face_y',
            face_var: str = None, #'mesh2d_nFaces',
    ):
        if time_var is None:
            time_var = self.coordinates.get('time')
        if layer_var is None:
            layer_var = self.coordinates.get('layer')
        if x_var is None:
            x_var = self.coordinates.get('lon')
        if y_var is None:
            y_var = self.coordinates.get('lat')
        if face_var is None:
            face_var = self.coordinates.get('face')

        dim = []
        ifilter = {}
        bboxfilter = {}
        if times is not None:
            if type(times) in [int]:
                times = [times]
            ifilter['times'] = times
            dim.append(time_var)
        if layers is not None:
            if type(layers) in [int]:
                layers = [layers]
            ifilter['layers'] = layers
            dim.append(layer_var)
        if bbox is not None:
            bboxfilter['bbox'] = bbox
            dim.append(face_var)

        datarray_var = self._obj
        datarray_var = datarray_var.d3dfm.ifilter(times, layers, time_var, layer_var)
        datarray_var = datarray_var.d3dfm.bboxfilter(bbox, x_var, y_var, face_var)
        datarray_var = datarray_var.mean(dim=dim)
        return datarray_var

    @staticmethod
    def _coordinate_type(
            variable: xr.DataArray
    ):
        """
        Return the coordinate types (lat, lon, time, layer, face) of a variable.

        :param variable:
        :return:
        """
        coord_type = None
        if variable.attrs.get('units') is not None:
            if variable.units in ['degrees_north', 'degree_north', 'degree_N', 'degrees_N', 'degreeN', 'degreesN']:
                coord_type = 'lat'
            elif variable.units in ['degrees_east', 'degree_east', 'degree_E', 'degrees_E', 'degreeE', 'degreesE']:
                coord_type = 'lon'
            else:
                for time_unit in ['second', 'seconds', 'hour', 'hours', 'day', 'days', 'month', 'months', 'year',
                                  'years']:
                    if time_unit in variable.units:
                        coord_type = 'time'
                        break

        if coord_type is None:
            variable_name = variable.name.lower()
            if 'nlayer' in variable_name:
                coord_type = 'layer'
            elif 'nface' in variable_name:
                coord_type = 'face'
            elif 'nedge' in variable_name:
                coord_type = 'edge'
            elif 'nnodes' in variable_name:
                coord_type = 'node'

        return coord_type

    def _get_coordinates(self):
        daData = self._obj
        coords = {}

        # Scan for Coordinate Variables
        for dim in daData.dims:
            daDim = daData[dim]
            coord_type = self._coordinate_type(daDim)
            if coord_type is not None:
                coords[coord_type] = daDim.name

        # Scan for Auxiliary Coordinate Variables
        if daData.encoding.get("coordinates") is not None:
            for coord in daData.encoding.get("coordinates").split():
                coord_type = self._coordinate_type(daData[coord])
                if coord_type is not None:
                    coords[coord_type] = coord

        # # Scan for Auxiliary Coordinate Variables 2
        # for coord in daData.coords:
        #     daCoord = daData[coord]
        #     coord_type = self.coordinate_type(daCoord)
        #     if coord_type is not None:
        #         coords[coord_type] = daCoord.name

        return coords

    @property
    def coordinates(self):
        """
        Return a dictionary of the name of coordinate variable and auxiliary coordinate variable for latitude, longitude,
        layer and time.

        :return:
        """
        if self._coordinates is None:
            self._coordinates = self._get_coordinates()
        return self._coordinates


@timing
def rasterize(
        input_nc: str,
        output_nc: str,
        variables: list[str],
        times,
        layers,
        bbox,
        ncellx,
        ncelly,
        dtypes,
        timeshift=0):
    dataset = xr.open_dataset(
        input_nc,
        decode_times=False,
        chunks={'mesh2d_nLayers': 1, 'time': 1, 'mesh2d_nFaces': -1 }
    )

    try:
        decoded_time = dataset.d3dfm.decoded_time
        ds = xr.Dataset(
            data_vars=dict(
                time_index=(['time'], np.arange(len(decoded_time)) )
            ),
            coords=dict(
                time=decoded_time
            )
        )
        times = ds.time_index.sel(time=slice(times[0], times[1], int(times[2]))).data
    except:
        times = range(*[int(num) for num in times])

    raster_ds = dataset.d3dfm.rasterize_variables(
        variables=variables,
        times=times,
        layers=layers,
        bbox=bbox,
        ncellx=ncellx,
        ncelly=ncelly,
    )


    time_units = raster_ds.time.attrs['units']
    if 'second' in time_units:
        raster_ds = raster_ds.assign_coords(time=('time', raster_ds.time.data +timeshift, raster_ds.time.attrs))
    else:
        raise TypeError('Time unit is not in seconds.')

    if len(variables) == len(dtypes):
        packing = {var: dtype for var, dtype in zip(variables, dtypes) }
    elif len(dtypes) == 1:
        packing = {var: dtypes[0] for var in variables}

    write_job = raster_ds.d3dfm.to_packed_netcdf(
        output_nc,
        packing=packing,
        compute=False
    )
    with ProgressBar():
        print(f"Writing to {output_nc}")
        write_job.compute()


