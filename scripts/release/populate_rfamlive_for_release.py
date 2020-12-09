import argparse
import utils.db_utils as db


# --------------------------------------------------------------------------------------


def parse_arguments():
    """
    Basic argument parsing using
    :return:
    """

    parser = argparse.ArgumentParser("Runs table updates for new Rfam release")

    parser.add_argument("--family-ncbi",
                        help="Populates family_ncbi table with distinct family species",
                        action="store_true", default=False)
    parser.add_argument("--num-species",
                        help="Sets number of distinct species in family table",
                        action="store_true", default=False)
    parser.add_argument("--num-full", help="Sets number of full regions hits per family",
                        action="store_true", default=False)

    parser.add_argument("--num-genome-families", help="Sets number of families per genome",
                        action="store_true", default=False)

    parser.add_argument("--num-genome-families", help="Sets number of families per genome",
                        action="store_true", default=False)

    parser.add_argument("--genome-full", help="Sets number of full region hits per genome",
                        action="store_true", default=False)

    parser.add_argument("--all", help="Update all tables required for release",
                        action="store_true", default=False)

    return parser

# --------------------------------------------------------------------------------------


if __name__ == '__main__':

    parser = parse_arguments()
    args = parser.parse_args()

    if args.all is True:
        db.update_family_ncbi()
        db.set_number_of_species()
        db.set_num_full_sig_seqs()
        db.set_number_of_distinct_families_in_genome(None)
        db.set_number_of_genomic_significant_hits(None)


