import argparse
from nc import rasterize
from agol import upload_nc, upload_vectorize_nc, vectorize_item

# Define the command line arguments
parser = argparse.ArgumentParser(prog="CLI", description='Rasterize and upload data to ArcGIS Online')
subparsers = parser.add_subparsers(dest="subcommand")

# Rasterize subcommand
rasterize_parser = subparsers.add_parser("rasterize", help="Rasterize the netcdf of out from D3DFM")
rasterize_parser.add_argument("inputfile", type=str, help="Path to the input file", metavar="INPUTFILE1.TXT")

# Uploadagol subcommand
uploadagol_parser = subparsers.add_parser("uploadagol", help="Upload the rasterized netcdf to ArcGIS Online")
uploadagol_parser.add_argument("inputfile", type=str, help="Path to the input file", metavar="INPUTFILE2.TXT")
uploadagol_parser.add_argument('--vectorize', action='store_true', help='Vectorize the uploaded data')

# Vectorize subcommand
vectorize_parser = subparsers.add_parser("vectorize", help='Vectorize the uploaded imagery data')
vectorize_parser.add_argument("inputfile", type=str, help="Path to the input file", metavar="INPUTFILE3.TXT")

# Parse the command line arguments
args = parser.parse_args()
#args = parser.parse_args(["rasterize", r"inputfile_rasterise.txt" ]) #for testing

if args.subcommand == "rasterize":
    with open(args.inputfile, 'r') as file:
        file_list = [line.strip() for line in file]
        file_list[2] = [string.strip() for string in file_list[2].split(',')]
        file_list[3] = file_list[3].split(',')
        file_list[4] = [int(num) for num in file_list[4].split(',')]
        file_list[5] = [float(num) for num in file_list[5].split(',')]
        file_list[6] = int(file_list[6])
        file_list[7] = int(file_list[7])
        file_list[8] = [s for s in file_list[8].split(',')]
        try: #for backward compatible
            file_list[9] = int(file_list[9])
        except:
            pass

    print('Rasterizing data...')
    rasterize(*file_list)

elif args.subcommand == "uploadagol":
    with open(args.inputfile, 'r') as file:
        file_list = [line.strip() for line in file]

    if args.vectorize:
        print('Uploading data and generating vector field...')
        upload_vectorize_nc(*file_list)
    else:
        print('Uploading data...')
        upload_nc(*file_list)

elif args.subcommand == "vectorize":
    with open(args.inputfile, 'r') as file:
        file_list = [line.strip() for line in file]
    vectorize_item(*file_list)
else:
    print("No subcommand specified. Use 'rasterize' or 'uploadagol' or 'vectorize'.")

print("Completed!")