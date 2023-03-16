import argparse
from nc import rasterize
from agol import upload_nc, upload_vectorize_nc

# Define the command line arguments
parser = argparse.ArgumentParser(description='Rasterize and upload data to ArcGIS Online')
parser.add_argument('--rasterize', type=str, help='Rasterize the netcdf of out from D3DFM')
parser.add_argument('--uploadagol', type=str, help='Upload the rasterized netcdf to ArcGIS Online')
parser.add_argument('--vectorize', action='store_true', help='Vectorize the uploaded data')

# Parse the command line arguments
args = parser.parse_args()

# Parse each line of the input file into a list if the --rasterize option is specified
if args.rasterize:
    with open(args.rasterize, 'r') as file:
        # Split the file into lines and store each line as a separate item in a list
        file_list = [line.strip() for line in file]
        file_list[2] = [string.strip() for string in file_list[2].split(',')]
        file_list[3] = range(*[int(num) for num in file_list[3].split(',')])
        file_list[4] = [int(num) for num in file_list[4].split(',')]
        file_list[5] = [float(num) for num in file_list[5].split(',')]
        file_list[6] = int(file_list[6])
        file_list[7] = int(file_list[7])

    # Rasterize the data
    print('Rasterizing data...')
    # rasterize takes input_nc, output_nc, variables, times, layers, bbox, ncellx, ncelly and dtype as inputs
    rasterize(*file_list)

# Parse each line of the input file into a list if the --uploadagol option is specified
if args.uploadagol:
    with open(args.uploadagol, 'r') as file:
        # Split the file into lines and store each line as a separate item in a list
        file_list = [line.strip() for line in file]

    # Check if the --vectorize option is specified
    if args.vectorize:
        # Vectorize the uploaded data
        print('Uploading data and generating vector field...')
        #upload_vectorize_nc takes username, password, input_nc, output_name, var_ucx, var_ucy and output_vec_name as input parameters
        upload_vectorize_nc(*file_list)
    else:
        # Upload data
        print('Uploading data...')
        # upload_nc takes username, password, input_nc, output_name as input
        upload_nc(*file_list)

print("Completed!")