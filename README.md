# D3DFM
The Delft3D Flexible Mesh Processing Tool offers two key functions to enhance the use and visualization of modelling results from Delft3D FM:
* Conversion of modelling results from Delft3D FM, which are in netCDF format using the [UGRID convention](http://ugrid-conventions.github.io/ugrid-conventions/), into a netCDF format that adheres to the [Climate & Forecast (CF) Metadata Conventions](http://cfconventions.org/). The CF conversion enables direct support and visualisation for the resulting data in various CF compliant software applications, particularly [ArcGIS Pro](https://pro.arcgis.com/en/pro-app/latest/help/data/multidimensional/essential-netcdf-vocabulary.htm).

* Uploading the converted netCDF file into ArcGIS Online / ArcGIS Server as an Imagery Layer, allowing users to easily visualise the modelling results via a web browser.

## Requirements

Create a conda virtual environment of Python 3.9:

    conda create -n D3DFM python=3.9

Set up the environment with environment.yml:

    conda env update --file environment.yml --prune


## Gallery
|                                                                                                       |                                                                                                          |
|-------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------|
| ![Photo 1](https://raw.githubusercontent.com/KenTSUI-dev/D3DFM/main/resources/ArcGISPro_D3DFM_nc.png) | ![Photo 2](https://raw.githubusercontent.com/KenTSUI-dev/D3DFM/main/resources/ArcGISOnline_D3DFM_nc.png) |
| *Rasterized D3DFM NetCDF in ArcGIS Pro*                                                               | *Rasterized D3DFM NetCDF  in a web browser.*                                                       |


## Command Line Interface (`cli.py`)
This script is designed to rasterize and upload D3DFM NetCDF data to ArcGIS Online, and optionally generate a vector field from the uploaded data. The script includes two options:

* --rasterize: Rasterizes (i.e. convert from UGRID to CF) an input NetCDF file using information from an input file, and saves the output rasterized netCDF file to a specified location.
* --uploadagol: Uploads an input rasterized NetCDF to ArcGIS Online using information from an input file, and optionally generates a vector field from the uploaded data using the --vectorize option.

### Usage

Example 1: Rasterize a NetCDF file

    python cli.py --rasterize input_rasterize.txt

In this example, the script will rasterize the D3DFM NetCDF file based on the input parameters in the "input_rasterize.txt" file.

Example 2: Upload a rasterized NetCDF file

    python cli.py --uploadagol input_upload.txt

In this example, the script will upload the rasterized NetCDF file to ArcGIS Online based on the input parameters in the "input_upload.txt" file.

Example 3: Upload a rasterized NetCDF file and generate a vector field service from the uploaded data

    python cli.py  --uploadagol input_upload.txt --vectorize
In this example, the script will upload the rasterized NetCDF file to ArcGIS Online based on the input parameters in the "input_upload.txt" file. The uploaded data will be used to generate a vector field service.

### Input File Format
The input file for the `--rasterize` option should contain nine lines of data:

1. Path to input NetCDF file
2. Path to output NetCDF file
3. Comma-separated list of variable names to rasterize
4. Range of times to rasterize (in the format start,end,step)
5. Comma-separated list of layer names to rasterize
6. Comma-separated list of bounding box coordinates (in the format xmin,ymin,xmax,ymax)
7. Number of cells in x direction
8. Number of cells in y direction
9. Data type for the output NetCDF file (int8 or int16)

Each line should be separated by a newline character.  A sample input file is as follows:

        D:\temp\D3DFM_2022Sep\HK-FM_merged_mapSep.nc
        D:\temp\D3DFM_2022Sep\HK-FM_merged_mapSep_raster4.nc
        mesh2d_sa1
        0, 720, 1
        0, 10, 19
        113.213111, 21.917770, 114.627601, 23.145613
        800 
        800
        int8



The input file for the `--uploadagol` option should contain four or seven lines of data, depending on whether the `--vectorize` option is specified:

1. ArcGIS Online username
2. ArcGIS Online password
3. Path to input rasterized NetCDF file
4. Name of output service
5. Name of x-component of vector field variable (if --vectorize is specified)
6. Name of y-component of vector field variable (if --vectorize is specified)
7. Name of output vector field service (if --vectorize is specified)

Each line should be separated by a newline character.  A sample input file is as follows:

        username
        password
        D:\temp\20230118_Tree\HK-FM_merged_map_ucxucy_dry1.nc
        Test_Dry_ucxucy
        mesh2d_ucx
        mesh2d_ucx
        Test_Dry_VecField