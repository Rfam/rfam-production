HTML_REPORTS = '/nfs/public/rw/xfam/rfam/' + 'searches/mirbase'
OUTPUT_FILENAME = 'mirbase-dashboard.tsv'
BLACK_LIST = [
    'MIPF0000338__mir-680'
    'MIPF0000419__mir-574',
    'MIPF0000901__mir-3470',
    'MIPF0001768__mir-7407',  # rfmake running
]  # giant families that cause crashes
COMPARE_OUTPUT_FILENAME = 'rfam-non-mirna-families-matching-mirbase.tsv'
SEARCH_DIR = '/nfs/production/agb/rfam/microrna/searches'
