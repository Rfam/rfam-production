import os
import sys
import shutil

# -----------------------------------------------------------------------------------------

DATABASE_FILE_NAMES = ["alignment_and_tree", "clan", "clan_database_link",
                       "clan_literature_reference", "clan_membership", "database_link",
                       "db_version", "dead_clan", "dead_family", "family",
                       "family_literature_reference", "family_ncbi", "features",
                       "full_region", "html_alignment", "keywords", "literature_reference",
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

    file_types = [".txt", ".sql"]

    if not os.path.exists(dest_dir):
        os.mkdir(dest_dir)

    for filename in DATABASE_FILE_NAMES:
        for file_type in file_types:
            source_file = os.path.join(source, filename + file_type)
            dest_file = os.path.join(dest_dir, filename + file_type)

            shutil.copyfile(source_file, dest_file)

# -----------------------------------------------------------------------------------------

if __name__ == "__main__":

    #TODO: Use argpase for proper argument parsing

    source_dir = sys.argv[1]
    dest_dir = sys.argv[2]

    create_ftp_database_files_folder(source_dir, dest_dir)
