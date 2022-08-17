
process create_files {
    publishDir "$params.release/metadata", mode: "copy"

    output:
    path('*_info_for_import.txt')

    """
    python $params.rfamprod/scripts/release/prepare_metadata_for_import.py -d $params.metadata_files > prepare_output.txt
    """
}

process import_to_db {
    input:
    path(query)

    """
    python $params.rfamprod/scripts/release/generate_metadata_tables.py
    """
}

workflow {
    create_files \
    | import_to_db
}
