MIRNAS_CSV = "rfam-production/mirnas_sample.csv"
SEARCH_DIRS = ["/nfs/production/agb/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk1_searches",
               "/nfs/production/agb/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk2_searches",
               "/nfs/production/agb/rfam//RELEASES/14.3/miRNA_relabelled/batch2/searches",
               "/nfs/production/agb/rfam/RELEASES/14.8/microrna/batch3_chunk1_searches"
               "/nfs/production/agb/rfam/RELEASES/14.9/microrna/batch3_chunk2_searches",
               "/nfs/production/agb/rfam/RELEASES/14.9/microrna/batch4_searches"]

STK_DIRS = ["/nfs/production/agb/rfam/microrna/batch4"]
# "/nfs/production/agb/rfam/microrna/batch3/fixed", "/nfs/production/agb/rfam/microrna/batch3/one_seed"

COPY_DIR = "/nfs/production/agb/rfam/RELEASES/14.9/microrna/batch4_searches"

UPDATE_DIR = "/nfs/production/agb/rfam/RELEASES/14.9/microrna/update_families"
NEW_DIR = "/nfs/production/agb/rfam/RELEASES/14.9/microrna/new_families"
ENV_PATH = "/nfs/production/agb/users/rfamprod/code/env2/bin/activate"
DESC_GEN_PATH = "/nfs/production/agb/users/rfamprod/code/rfam-production/scripts/preprocessing/desc_generator.py"
MEMORY = 8000
CPU = 4
LSF_GROUP = "/rfam_srch"
QUEUE = "short"