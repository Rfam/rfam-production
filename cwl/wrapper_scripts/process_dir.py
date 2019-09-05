import argparse
import shutil
import os
from subprocess import check_output, CalledProcessError, Popen, PIPE
from pathlib import Path
from yaml import dump
import logging

REFERENCE_DB_NAME = 'RF02567.cm'
CWL_EXECUTOR = 'toil-cwl-runner'
CWL_WORKFLOW_OPTIONS = []

current_dir = os.path.dirname(__file__)
WORKFLOW_PATH = os.path.realpath(os.path.join(current_dir, os.pardir, 'workflows', 'nc_rna_workflow.cwl'))

REFERENCE_CM_DB = os.path.realpath(os.path.join(current_dir, os.pardir, 'db', REFERENCE_DB_NAME))

GENOME_FILE_FORMATS = ('.fa', '.fasta')

CONFIG_FILENAME = 'config.yaml'


def validate_workflow():
    logging.info('Validating CWL workflow')
    cmd = [CWL_EXECUTOR, '--validate', WORKFLOW_PATH]
    try:
        check_output(cmd)
    except CalledProcessError:
        cmd = ' '.join(cmd)
        raise EnvironmentError(f'Failed to validate workflow, debug using "{cmd}".')


def check_setup():
    logging.info('Checking project setup')
    if not shutil.which(CWL_EXECUTOR):
        raise EnvironmentError(f'CWL executor {CWL_EXECUTOR} not found in path.')
    if not os.path.exists(WORKFLOW_PATH):
        raise EnvironmentError(f'Workflow not found. {WORKFLOW_PATH}')
    if not os.path.exists(REFERENCE_CM_DB):
        raise EnvironmentError(f'CM db not found at path {REFERENCE_CM_DB}.')
    validate_workflow()


def parse_args():
    parser = argparse.ArgumentParser(description='Tool to find and execute nc_rna_pipeline for all genome fasta files')
    parser.add_argument('directory', help='Directory to search for fasta files')
    parser.add_argument('--dry-mode', action='store_true', help='logging.info pipeline commands without launching them')
    parser.add_argument('--skip-setup-check', action='store_true', help='Skips workflow validation')
    parser.add_argument('-v', '--verbose', action='store_true')
    return parser.parse_args()


def find_seq_datasets(directory):
    files = Path(directory).glob('**/*')
    fasta_files = filter(lambda fp: any([str(fp).endswith(fmt) for fmt in GENOME_FILE_FORMATS]), files)
    fasta_files = map(lambda fp: fp.resolve(), fasta_files)
    return list(fasta_files)


def write_cwl_config_file(config_dir, fasta_file):
    logging.info('\twriting config')
    config = {
        'sequences': {
            'class': 'File',
            'format': 'http://edamontology.org/format_1929',
            'path': fasta_file
        },
        'covariance_model_database': {
            'class': 'File',
            'path': REFERENCE_CM_DB
        }
    }
    outfile = os.path.join(config_dir, CONFIG_FILENAME)
    with open(outfile, 'w') as f:
        dump(config, f)
    return outfile


def launch_pipeline(outdir, config):
    logging.info('\tlaunching pipeline')
    cmd = [CWL_EXECUTOR, *CWL_WORKFLOW_OPTIONS, '--outdir', outdir, WORKFLOW_PATH, config]
    p = Popen(cmd, stdout=PIPE, stderr=PIPE)
    logging.info('\tpipeline running')
    logging.debug('\t' + ' '.join(cmd))
    output, err = p.communicate()

    rc = p.returncode
    if rc:
        logging.info('\tAn error occurred in workflow')
        logging.info(err)
    else:
        logging.info('\tWorkflow completed successfully')


def main():
    args = parse_args()
    logging.basicConfig(format='%(message)s', level=logging.DEBUG if args.verbose else logging.INFO)

    if not args.skip_setup_check:
        check_setup()
    logging.info(f'Searching for mags in {args.directory}')
    seq_datasets = find_seq_datasets(args.directory)
    logging.info(f'Found {len(seq_datasets)} mags...')

    for sd in seq_datasets:
        logging.info(f'Processing {os.path.basename(sd)}')
        outdir = os.path.join(os.path.dirname(sd), 'out')

        os.makedirs(outdir, exist_ok=True)

        if args.dry_mode:
            logging.info('Skipping config + pipeline (dry mode)')
        else:
            config = write_cwl_config_file(os.path.dirname(sd), str(sd))
            launch_pipeline(outdir, config)
        logging.info('\n')


if __name__ == '__main__':
    main()
