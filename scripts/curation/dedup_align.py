#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "loguru>=0.7.0",
#     "click>=8.0.0",
# ]
# ///

"""
Stockholm file processor that removes duplicates and filters sequences.
Uses esl-reformat, seqkit, and esl-alimanip to process Stockholm format files.
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import click
from loguru import logger


def setup_logging(verbose: bool = False) -> None:
    """Configure logging with loguru."""
    # Remove default handler
    logger.remove()

    # Add stderr handler with appropriate level
    log_level = "DEBUG" if verbose else "INFO"
    logger.add(
        sys.stderr,
        level=log_level,
        format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | <level>{message}</level>",
    )


def process_stockholm_file(input_data: str, input_filename: str = "SEED") -> str:
    """
    Process Stockholm file data through the pipeline:
    1. Convert to FASTA format
    2. Remove duplicates based on sequence
    3. Extract sequence IDs to keep
    4. Filter original Stockholm file to keep only non-duplicate sequences
    """
    logger.info("Starting Stockholm file processing pipeline")

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)

        # Create temporary input file
        input_file = temp_path / "input.sto"
        keep_file = temp_path / "keep"

        # Write input data to temporary file
        input_file.write_text(input_data)
        logger.debug(f"Created temporary input file: {input_file}")

        # Step 1-3: Pipe esl-reformat -> seqkit rmdup -> seqkit seq to create keep file
        logger.info("Creating pipeline: esl-reformat | seqkit rmdup | seqkit seq")

        # Start the pipeline processes
        proc1 = subprocess.Popen(
            ["esl-reformat", "fasta", str(input_file)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        proc2 = subprocess.Popen(
            ["seqkit", "rmdup", "-P", "-s"],
            stdin=proc1.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        proc3 = subprocess.Popen(
            ["seqkit", "seq", "--only-id", "-n"],
            stdin=proc2.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        # Close upstream stdout pipes so processes can receive SIGPIPE if downstream exits
        proc1.stdout.close()
        proc2.stdout.close()

        # Wait for the pipeline to complete and get the final output
        keep_ids, stderr3 = proc3.communicate()

        # Check return codes
        proc1_ret = proc1.wait()
        proc2_ret = proc2.wait()
        proc3_ret = proc3.returncode

        if proc1_ret != 0:
            logger.error(f"esl-reformat failed with return code {proc1_ret}")
            stderr1 = proc1.stderr.read() if proc1.stderr else ""
            if stderr1:
                logger.error(f"esl-reformat stderr: {stderr1}")
            raise subprocess.CalledProcessError(proc1_ret, "esl-reformat")

        if proc2_ret != 0:
            logger.error(f"seqkit rmdup failed with return code {proc2_ret}")
            stderr2 = proc2.stderr.read() if proc2.stderr else ""
            if stderr2:
                logger.error(f"seqkit rmdup stderr: {stderr2}")
            raise subprocess.CalledProcessError(proc2_ret, "seqkit rmdup")

        if proc3_ret != 0:
            logger.error(f"seqkit seq failed with return code {proc3_ret}")
            if stderr3:
                logger.error(f"seqkit seq stderr: {stderr3}")
            raise subprocess.CalledProcessError(proc3_ret, "seqkit seq")

        # Write keep file
        keep_file.write_text(keep_ids)
        num_sequences = len(keep_ids.strip().split()) if keep_ids.strip() else 0
        logger.debug(f"Created keep file with {num_sequences} sequences")

        # Step 4: Filter original Stockholm file using keep list
        logger.info("Filtering original Stockholm file")
        result = subprocess.run(
            ["esl-alimanip", "--seq-k", str(keep_file), str(input_file)],
            capture_output=True,
            text=True,
            check=True,
        )

        logger.info("Processing completed successfully")
        return result.stdout


@click.command()
@click.argument(
    "input_file", type=click.Path(path_type=Path), required=False, default="-"
)
@click.argument(
    "output_file", type=click.Path(path_type=Path), required=False, default="-"
)
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose logging")
@click.version_option(version="1.0.0")
def main(input_file: Path, output_file: Path, verbose: bool) -> None:
    """
    Process Stockholm format files to remove duplicate sequences.

    This script processes Stockholm alignment files by:
    1. Converting to FASTA format
    2. Removing duplicate sequences based on sequence content
    3. Filtering the original Stockholm file to keep only unique sequences

    Arguments:
      INPUT_FILE   Input Stockholm file (use '-' for stdin, default: stdin)
      OUTPUT_FILE  Output file (use '-' for stdout, default: stdout)

    Examples:
      stockholm_processor.py input.sto output.sto    # File to file
      stockholm_processor.py - output.sto            # Stdin to file
      stockholm_processor.py - -                     # Stdin to stdout
      stockholm_processor.py                         # Stdin to stdout (default)
    """
    setup_logging(verbose)

    logger.info("Stockholm file processor starting")

    try:
        # Read input data
        if input_file == Path("-") or str(input_file) == "-":
            logger.info("Reading input from stdin")
            input_data = sys.stdin.read()
            input_filename = "SEED"

            if not input_data.strip():
                logger.error("No input data provided")
                sys.exit(1)
        else:
            if not input_file.exists():
                logger.error(f"Input file does not exist: {input_file}")
                sys.exit(1)
            logger.info(f"Reading input from file: {input_file}")
            input_data = input_file.read_text()
            input_filename = input_file.name

        # Process the data
        output_data = process_stockholm_file(input_data, input_filename)

        # Write output
        if output_file == Path("-") or str(output_file) == "-":
            logger.info("Writing output to stdout")
            sys.stdout.write(output_data)
        else:
            logger.info(f"Writing output to file: {output_file}")
            output_file.write_text(output_data)

        logger.info("Processing completed successfully")

    except Exception as e:
        logger.error(f"Error during processing: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
