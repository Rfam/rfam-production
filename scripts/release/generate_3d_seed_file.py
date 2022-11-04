import os

import argparse


def create_3d_seed_file(file_with_list):
    """
    Take the IDs of all Rfam families that have been updated with 3D info
    and merge their seed files to one file
    """
    out_file = "nfs/production/agb/rfam/RELEASES/14.9/ftp/seed/Rfam.3d.seed"
    seed_dir = "nfs/production/agb/rfam/RELEASES/14.9/ftp/seed"
    with open(file_with_list, 'r') as f:
        rfam_accs = f.read()

    for acc in rfam_accs:
        seed_path = os.path.join(seed_dir, acc)
        with open(seed_path, 'r') as i:
            with open(out_file, 'w') as o:
                for line in i:
                    o.write(line)


def parse_arguments():
    """
    Performs some basic argument parsing

    return: parser object
    """

    parser = argparse.ArgumentParser(description='Generates seed file for 3D families')

    parser.add_argument("-f", help="File with list of names of families that have been updated with 3D info",
                        action="store")

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    create_3d_seed_file(args.f)
