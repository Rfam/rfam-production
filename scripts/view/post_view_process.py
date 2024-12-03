from utils import db_utils as RfamDB

"""
To be used for updating all fields prior to indexing
"""

# -----------------------------------------------------------------

if __name__ == "__main__":

    """
    parser = argparse.ArgumentParser(description='Update table fields for a new release')
    parser.add_argument('--all', action='store_const',
                        const=sum
                        help='Update fields for all accessions)
    """

    # FAMILIES
    RfamDB.set_num_full_sig_seqs()
    RfamDB.set_number_of_species()
    # truncate family_ncbi table before executing this
    RfamDB.update_family_ncbi()

    # GENOMES
    RfamDB.set_number_of_distinct_families_in_genome(upid=None)
    RfamDB.set_genome_size(upid=None)
    RfamDB.set_number_of_distinct_families_in_genome(upid=None)
    RfamDB.set_number_of_genomic_significant_hits(upid=None)
