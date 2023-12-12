import os
import shutil
import argparse

# -----------------------------------------------------------------------------------------

DATABASE_FILE_NAMES = ["alignment_and_tree", "clan", "clan_database_link",
                       "clan_literature_reference", "clan_membership", "database_link",
                       "db_version", "dead_clan", "dead_family", "family",
                       "family_literature_reference", "family_ncbi", "features",
                       "full_region", "genome", "genseq", "html_alignment", "keywords", "literature_reference",
                       "matches_and_fasta", "motif", "motif_database_link",
                       "motif_family_stats", "motif_file", "motif_literature",
                       "motif_matches", "motif_pdb", "motif_ss_image", "pdb_full_region",
                       "rfamseq", "secondary_structure_image", "seed_region",
                       "sunburst", "taxonomy", "taxonomy_websearch", "version",
                       "wikitext"]

# -----------------------------------------------------------------------------------------

def create_ftp_database_files_folder(source, dest_dir):
    """

    :param source:
    :param dest_dir:
    :return:
    """

    file_types = [".txt.gz", ".sql"]

    if not os.path.exists(dest_dir):
        os.mkdir(dest_dir)

    for filename in DATABASE_FILE_NAMES:
        for file_type in file_types:
            source_file = os.path.join(source, filename + file_type)
            dest_file = os.path.join(dest_dir, filename + file_type)

            shutil.copyfile(source_file, dest_file)

# -----------------------------------------------------------------------------------------


def parse_arguments():
    """
    Basic Argument parsing using python's argparse

    return: Argparse parser object
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("--source-dir", help="Source directory containing a mysqldump database dump",
                        action="store")
    parser.add_argument("--dest-dir", help="Destination directory to create ftp database_files",
                        action="store")

    return parser

# -----------------------------------------------------------------------------------------#


if __name__ == "__main__":

    #TODO: Use argpase for proper argument parsing

    parser = parse_arguments()
    args = parser.parse_args()

    source_dir = args.source_dir
    dest_dir = args.dest_dir

    create_ftp_database_files_folder(source_dir, dest_dir)
