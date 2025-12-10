nextflow.enable.dsl=2

// Validate required S3 parameters at workflow startup
if (!params.s3?.s3_host) {
    error "Missing required parameter: params.s3.s3_host"
}
if (!params.s3?.s3_key) {
    error "Missing required parameter: params.s3.s3_key"
}
if (!params.s3?.s3_secret) {
    error "Missing required parameter: params.s3.s3_secret"
}
if (!params.s3?.bucket_name) {
    error "Missing required parameter: params.s3.bucket_name"
}
if (!params.s3?.environment) {
    error "Missing required parameter: params.s3.environment"
}
if (!params.s3?.s3_base_path) {
    error "Missing required parameter: params.s3.s3_base_path"
}

process run_svn_to_s3 {
    memory '4 GB'
    cpus 2
    errorStrategy 'terminate'
    
    input:
    val flag
    
    output:
    val true
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    # Set environment variables for the script
    export S3_HOST="${params.s3.s3_host}"
    export S3_KEY="${params.s3.s3_key}"
    export S3_SECRET="${params.s3.s3_secret}"
    export BUCKET_NAME="${params.s3.bucket_name}"
    export ENVIRONMENT="${params.s3.environment}"
    export S3_BASE_PATH="${params.s3.s3_base_path}"
    
    # Make script executable if needed
    chmod +x ${projectDir}/rna-alignment-api/scripts/svn_to_s3.sh
    
    # Run the SVN to S3 script
    bash ${projectDir}/rna-alignment-api/scripts/svn_to_s3.sh
    
    echo "SVN to S3 sync completed successfully"
    """
}

workflow update_stockholm_s3 {
    take: 
    start
    
    main:
    run_svn_to_s3(start)

    emit:
    done = run_svn_to_s3.out
}

workflow {
    update_stockholm_s3(Channel.of('start'))
}