import argparse
import os


def create_3d_seed_file(file_with_list, seed_dir):
    """
    Take the IDs of all Rfam families that have been updated with 3D info
    and merge their seed files to one file
    """
    out_file = os.path.join(seed_dir, "Rfam.3d.seed")

    with open(file_with_list, "r") as f:
        rfam_accs = f.read().splitlines()

    for acc in rfam_accs:
        seed_path = os.path.join(seed_dir, acc + ".seed")
        with open(seed_path, "r") as i:
            with open(out_file, "a") as o:
                for line in i:
                    o.write(line)


def parse_arguments():
    """
    Performs some basic argument parsing

    return: parser object
    """

    parser = argparse.ArgumentParser(description="Generates seed file for 3D families")
    required_arguments = parser.add_argument_group("required arguments")

    required_arguments.add_argument(
        "-f",
        help="File with list of names of families that have been updated with 3D info",
        action="store",
    )
    required_arguments.add_argument(
        "--seed-dir",
        help="Directory of seed files, where to write file to",
        action="store",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    create_3d_seed_file(args.f, args.seed_dir)
