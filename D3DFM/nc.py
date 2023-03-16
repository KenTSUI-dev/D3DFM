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
                        "units": "layer",
                        "positive" : "down",
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

                dtype_to_nbits = {
                    "int16": {'nbits': 16, '_FillValue': -32768},
                    "int8": {'nbits': 8, '_FillValue': -128}
                }

                if dtype in dtype_to_nbits:
                    nbits = dtype_to_nbits[dtype]['nbits']
                    _FillValue = dtype_to_nbits[dtype]['_FillValue']
                    try:
                        scale_factor, add_offset = compute_scale_and_offset(
                            self.get_var(variable).d3dfm.ifilter(times=-1),
                            n=nbits,
                            num_reserveed_tail=1
                        )
                    except:
                        scale_factor, add_offset = compute_scale_and_offset(
                            self.get_var(variable),
                            n=nbits,
                            num_reserveed_tail=1
                        )

                    packcoding[variable] = {
                            "dtype": dtype,
                            "scale_factor": scale_factor,
                            "add_offset": add_offset,
                            "_FillValue": _FillValue,
                    }

                else:
                    raise ValueError(f"Data type: {dtype} is not supported. Only accept 'int8' or 'int16'.")

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
def rasterize(input_nc, output_nc, variables, times, layers, bbox, ncellx, ncelly, dtype):
    dataset = xr.open_dataset(
        input_nc,
        decode_times=False,
        chunks={'mesh2d_nLayers': 1, 'time': 1, 'mesh2d_nFaces': -1 }
    )

    raster_ds = dataset.d3dfm.rasterize_variables(
        variables=variables,
        times=times,
        layers=layers,
        bbox=bbox,
        ncellx=ncellx,
        ncelly=ncelly,
    )
    write_job = raster_ds.d3dfm.to_packed_netcdf(
        output_nc,
        packing={var: dtype for var in variables},
        compute=False
    )
    with ProgressBar():
        print(f"Writing to {output_nc}")
        write_job.compute()


@timing
def sample_rasterize():
    # client = Client(n_workers=2, threads_per_worker=2, memory_limit='1GB')
    file_nc = r"D:\temp\D3DFM_2022Sep\HK-FM_merged_mapSep.nc" #r"D:\temp\D3DFM_2022Sep\HK-FM_merged_mapSep_duplicate.nc"
    out_file = r"D:\temp\D3DFM_2022Sep\HK-FM_merged_mapSep_raster4.nc"
    dataset = xr.open_dataset(
        file_nc,
        decode_times=False,
        chunks={'mesh2d_nLayers': 1, 'time': 1, 'mesh2d_nFaces': -1 }
    )

    raster_ds = dataset.d3dfm.rasterize_variables(
        variables="mesh2d_sa1",
        times=range(30 * 24),
        layers=[0, 10, 19],
        bbox=[113.213111, 21.917770, 114.627601, 23.145613],
        ncellx=800,
        ncelly=800,
    )
    write_job = raster_ds.d3dfm.to_packed_netcdf(
        out_file,
        packing={"mesh2d_sa1": "int8"},
        compute=False
    )
    with ProgressBar():
        print(f"Writing to {out_file}")
        write_job.compute()

def sample_rasterize_ucxucy():
    # client = Client(n_workers=2, threads_per_worker=2, memory_limit='1GB')
    file_nc = r"D:\temp\D3DFM_2022Sep\HK-FM_merged_mapSep.nc" #r"D:\temp\D3DFM_2022Sep\HK-FM_merged_mapSep_duplicate.nc"
    out_file = r"D:\temp\D3DFM_2022Sep\HK-FM_merged_mapSep_ucxucy.nc"
    dataset = xr.open_dataset(
        file_nc,
        decode_times=False,
        chunks={'mesh2d_nLayers': 1, 'time': 1, 'mesh2d_nFaces': -1 }
    )

    raster_ds = dataset.d3dfm.rasterize_variables(
        variables=["mesh2d_ucx", "mesh2d_ucy"],
        times=range(30 * 24),
        layers=range(1),
        bbox=[113.213111, 21.917770, 114.627601, 23.145613],
        ncellx=800,
        ncelly=800,
    )
    write_job = raster_ds.d3dfm.to_packed_netcdf(
        out_file,
        packing={"mesh2d_ucx": "int8", "mesh2d_ucy": "int8"},
        compute=False
    )
    with ProgressBar():
        print(f"Writing to {out_file}")
        write_job.compute()


@timing
def sample_mean():
    file_nc = r"D:\temp\D3DFM_2022Sep\HK-FM_merged_mapSep.nc" #r"D:\temp\D3DFM_2022Sep\HK-FM_merged_mapSep_duplicate.nc"
    dataset = xr.open_dataset(
        file_nc,
        decode_times=False,
        chunks={'mesh2d_nLayers': 1, 'time': 1, 'mesh2d_nFaces': -1}
    )

    daData = dataset.mesh2d_sa1
    print(daData.d3dfm.coordinates)
    daData_mean = daData.d3dfm.mean(layers=range(20))
    # daData_mean.to_netcdf(r"D:\temp\D3DFM_2022Sep\HK-FM_merged_mesh2d_sa1_mean.nc", mode='w', )



    # print(dataset.mesh2d_sa1.d3dfm.mean(range(720), range(10)))
@timing
def sample_merge_map_files():
    # Rules:
    # All datasets should contain the concatenating dimensions
    # Each dataset may not need to contain the same variables.
    # The same variable in all datasets should contain the same dimensions.
    datasets = [xr.open_dataset(filepath, chunks={'mesh2d_nLayers': -1, 'time': -1, 'mesh2d_nFaces': -1})
                 for filepath in glob.glob(r"D:\temp\D3DFM_2022Sep\HK-FM-Partitioned\HK-FM_00*_map.nc")]
    dims = ['mesh2d_nFaces', 'mesh2d_nEdges', 'mesh2d_nNodes']
    concat_dims = set(dims)
    datasets_contain_concat_dims = [concat_dims.issubset(set(ds.dims)) for ds in datasets]
    if not all(datasets_contain_concat_dims):
        raise ValueError(f'Not all datasets contain the dims - {dims}.')
    all_datasets_vars = {var for ds in datasets for var in ds.data_vars}
    var_to_datasets = {var: [ds for ds in datasets if var in ds.data_vars ] for var in all_datasets_vars }

    def all_equal(iterator):
        iterator = iter(iterator)
        try:
            first = next(iterator)
        except StopIteration:
            return True
        return all(first == x for x in iterator)

    results={}
    for var, var_datasets in var_to_datasets.items():
        datasets_var_sizes = [ds[var].sizes for ds in datasets]
        datasets_var_sizes_no_concat_dims = [
            {dim:length for dim, length in size.items() if dim not in concat_dims}
            for size in datasets_var_sizes
        ]
        datasets_var_dims = [ds[var].dims for ds in datasets]
        if not all_equal(datasets_var_sizes_no_concat_dims) or not all_equal(datasets_var_dims):
            continue

        first_dataset = var_datasets[0]
        var_dims = set(first_dataset[var].dims)
        common_concat_dims = concat_dims-(concat_dims-var_dims)
        if len(common_concat_dims) == 0:
            results[var] = first_dataset[var]
        elif len(common_concat_dims) == 1:
            if var == "mesh2d_face_nodes":
                offset = 0
                data_vars = []
                for var_dataset in var_datasets:
                    data_var = var_dataset[var].copy()
                    data_var.data += offset
                    data_vars.append(data_var)
                    offset += var_dataset.sizes['mesh2d_nNodes']
                results[var] = xr.concat(data_vars, dim=list(common_concat_dims)[0])
            else:
                results[var] = xr.concat([var_dataset[var] for var_dataset in var_datasets], dim=list(common_concat_dims)[0])
        elif len(common_concat_dims) >= 2:
            print('unexpected')

    combined_ds = xr.Dataset(
        data_vars=results
    )
    i=0
    for key, var in results.items():
        mode = "w" if i==0 else "a"
        print(f"Start to write {key}. Mode: {mode}.")
        write_job = var.to_netcdf(
            r"D:\temp\D3DFM_2022Sep\HK-FM_merged_mapSep_test.nc",
            mode=mode,
            format='NETCDF4',
            # compute=False
        )
        i+=1
        # with ProgressBar():
        #     print(f"Start to write {key}")
        #     write_job.compute()
    # print(combined_ds['mesh2d_sa1'].chunksizes)

    # combined_ds.to_netcdf(r"D:\temp\D3DFM_2022Sep\HK-FM_merged_mapSep_test.nc", format='NETCDF4')
    # print(combined_ds.d3dfm.get_var("mesh2d_sa1").attrs)
    # combined_ds.d3dfm.mesh2d_face_to_file(r"D:\temp\D3DFM_2022Sep\xr_combined_d3dfm_nc.shp")
    out_file = r"D:\temp\D3DFM_2022Sep\HK-FM_merged_mapSep_ucxucy.nc"

    # raster_ds = combined_ds.d3dfm.rasterize_variables(
    #     variables=["mesh2d_ucx", "mesh2d_ucy"],
    #     times=range(30 * 24),
    #     layers=range(1),
    #     bbox=[113.213111, 21.917770, 114.627601, 23.145613],
    #     ncellx=800,
    #     ncelly=800,
    # )
    # write_job = raster_ds.d3dfm.to_packed_netcdf(
    #     out_file,
    #     packing={"mesh2d_ucx": "int8", "mesh2d_ucy": "int8"},
    #     compute=False
    # )
    # with ProgressBar():
    #     print(f"Writing to {out_file}")
    #     write_job.compute()



if __name__ == "__main__":
    sample_rasterize()
    # sample_rasterize_ucxucy()
    # sample_mean()
    # sample_merge_map_files()

