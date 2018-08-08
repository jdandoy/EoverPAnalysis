import rucio
import subprocess
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='A script to get the remote directories of datasets from a specific grid site')

    parser.add_argument('--grid_site', required=True, type=str, dest='grid_site', help='Grid Site Name')
    parser.add_argument('--dataset_name', required=True, type=str, dest='dataset_name', help='Dataset Name')
    parser.add_argument('--output_filename', required=True, type=str, dest='output_filename', help='Name of output txt file')

    args = parser.parse_args()

    grid_site = args.grid_site
    dataset_name = args.dataset_name
    output_filename = args.output_filename

    files = subprocess.check_output(["rucio", "list-file-replicas", "--rse", grid_site, "--protocol", "root", dataset_name])

    list_of_files = files.split("\n")
    dataset_file_locations = []

    for line in list_of_files:
        if not "root" in line:
            continue
        dataset_filename = line.split("root://")[1]
        dataset_filename = dataset_filename.split(" | ")[0]
        while dataset_filename[-1] == " " or dataset_filename[-1] == "|" or dataset_filename[-1] == "\n":
            dataset_filename = dataset_filename[0:-1]
        while dataset_filename[0] == " " or dataset_filename[0] == "|":
            dataset_filename = dataset_filename[1:]
        dataset_filename = "root://" + dataset_filename

        dataset_file_locations.append(dataset_filename)

    thefile = file(output_filename, 'w')
    for item in dataset_file_locations:
      thefile.write("%s\n" % item)
