from arcgis.gis import GIS
from arcgis.raster.analytics import copy_raster
try:
    from cred import AGOLusername, AGOLpassword
except:
    print('cred.py is missing.')
from arcgis.raster.functions import vector_field, multidimensional_filter
from tools import retry

@retry(tries=10)
def upload_nc(username, password, input_nc, output_name):
    gis = GIS(username=username, password=password)
    item = copy_raster(
        input_raster=input_nc,
        output_name=output_name,
        raster_type_name="NetCDF",
        gis=gis
    )

    return item

@retry(tries=10)
def _vectorize_item(item, var_ucx, var_ucy, output_vec_name):
    layer = item.layers[0]
    raster_ucx  = multidimensional_filter(layer, variables=[var_ucx])
    raster_ucy  = multidimensional_filter(layer, variables=[var_ucy])
    raster_vecfield = vector_field( raster_ucx, raster_ucy, input_data_type='Vector-UV', output_data_type='Vector-UV')
    raster_vecfield.save(output_vec_name)

def upload_vectorize_nc(username, password, input_nc, output_name, var_ucx, var_ucy, output_vec_name):
    item = upload_nc(username, password, input_nc, output_name)
    _vectorize_item(item, var_ucx, var_ucy, output_vec_name)

@retry(tries=10)
def vectorize_item(username, password, item_id, var_ucx, var_ucy, output_vec_name):
    gis = GIS(username=username, password=password)
    item = gis.content.get(item_id)
    layer = item.layers[0]
    raster_ucx  = multidimensional_filter(layer, variables=[var_ucx])
    raster_ucy  = multidimensional_filter(layer, variables=[var_ucy])
    raster_vecfield = vector_field( raster_ucx, raster_ucy, input_data_type='Vector-UV', output_data_type='Vector-UV')
    raster_vecfield.save(output_vec_name)

if __name__ == "__main__":
    # item = upload_nc(AGOLusername, AGOLpassword, r"D:\temp\20230118_Tree\HK-FM_merged_map_sa_wet1.nc", 'Test_Wet_Salinity')
    upload_vectorize_nc(
        AGOLusername,
        AGOLpassword,
        r"D:\temp\20230118_Tree\HK-FM_merged_map_ucxucy_dry1.nc",
        'Test_Dry_ucxucy',
        'mesh2d_ucx',
        'mesh2d_ucx',
        'Test_Dry_VecField'
    )


