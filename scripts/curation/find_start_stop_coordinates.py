#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "biopython>=1.80",
#     "requests>=2.25.0",
#     "loguru>=0.7.0",
#     "click>=8.0.0",
#     "attrs",
# ]
# ///

from __future__ import annotations

import io
import os
import re
import shutil
import subprocess as sp
import sys
import tempfile
import typing as ty
from pathlib import Path

import click
import requests
from attrs import frozen
from Bio import Align, AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from loguru import logger


@frozen
class AccessionFetcher:
    """Handles fetching and caching of FASTA sequences from NCBI."""

    cache: dict[str, SeqRecord]
    timeout: int
    temp_dir: str

    @classmethod
    def build(cls, timeout: int = 30) -> AccessionFetcher:
        return cls(
            timeout=timeout,
            cache={},
            temp_dir=tempfile.mkdtemp(prefix="ncbi_sequences_"),
        )

    def __del__(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def __fetch_fasta__(self, accession: str) -> Path | None:
        try:
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            params = {
                "db": "nuccore",
                "id": accession,
                "rettype": "fasta",
                "retmode": "text",
            }

            response = requests.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()

            if not response.text.strip() or "ERROR" in response.text.upper():
                logger.warning("Could not fetch sequence for {}", accession)
                return None

            temp_file = os.path.join(self.temp_dir, f"{accession}.fasta")
            with open(temp_file, "w") as f:
                f.write(response.text)
            return Path(temp_file)
        except requests.RequestException as e:
            logger.warning("Network error fetching {}: {}", accession, e)
            return None

    def get_sequence_and_accession(self, accession: str) -> None | SeqRecord:
        """
        Fetch sequence and version for an accession.
        Returns SeqRecord or None if not found.
        """

        if accession in self.cache:
            return self.cache[accession]

        logger.debug("Fetching {} from NCBI", accession)

        path = self.__fetch_fasta__(accession)
        if not path:
            logger.warning("No fasta for {}, so skipping", accession)
            return None

        try:
            with path.open("r") as handle:
                fasta_record = SeqIO.read(handle, "fasta")

            sequence_accession = fasta_record.id
            if "." not in sequence_accession:
                logger.debug(
                    "Assuming version 1 of of accession {}, since it was missing",
                    sequence_accession,
                )
                sequence_accession = f"{sequence_accession}.1"
            logger.debug("Found complete accession {}", sequence_accession)

            result = SeqRecord(fasta_record.seq, id=sequence_accession, description="")
            self.cache[accession] = result
            return result

        except Exception as e:
            logger.warning("Could not parse FASTA file for {}: {}", accession, e)
            return None


def coordinates_within(sequence: Seq, target: Seq) -> None | tuple[int, int]:
    start = sequence.find(target)
    assert isinstance(start, int)
    if start == -1:
        return None
    assert sequence[start : start + len(target)] == target
    start = start + 1
    end = start + len(target) - 1
    return (start, end)


def normalize_sequence(seq: Seq) -> Seq:
    return seq.replace("-", "").replace(".", "").upper().replace("U", "T")


def read_clean_alignment(path: Path, format="stockholm") -> Align.MultipleSeqAlignment:
    with tempfile.NamedTemporaryFile() as out:
        try:
            sp.check_call(
                [
                    "esl-reformat",
                    "-d",
                    "-u",
                    "--informat",
                    format,
                    "--mingap",
                    "stockholm",
                    str(path),
                ],
                stdout=out,
            )
        except Exception as e:
            logger.error(
                "Cannot format input alignment with `esl-reformat --mingap -d -u`"
            )
            logger.error(e)
        out.flush()
        return AlignIO.read(out.name, "stockholm")


def forward_name(record: SeqRecord, start: int, stop: int) -> str:
    return f"{record.id}/{start}-{stop}"


def reverse_complement_name(record: SeqRecord, start: int, stop: int) -> str:
    sequence_length = len(record.seq)
    forward_start = sequence_length - stop + 1
    forward_stop = sequence_length - start + 1
    return f"{record.id}/{forward_stop}-{forward_start}"


def search_sequences(
    sequence: Seq,
) -> ty.Iterable[tuple[str, Seq, ty.Callable[[SeqRecord, int, int], str]]]:
    yield ("forward", sequence, forward_name)
    yield ("reverse_complement", sequence.reverse_complement(), reverse_complement_name)


def validate_accession(record: SeqRecord, query: str, expected: Seq) -> bool:
    with tempfile.NamedTemporaryFile("a") as fasta:
        SeqIO.write(record, fasta, "fasta")
        fasta.flush()
        expected_accession, endpoints = query.split("/", 1)
        sp.check_call(
            ["esl-sfetch", "--index", fasta.name], stdout=sp.DEVNULL, stderr=sp.DEVNULL
        )
        output = sp.check_output(
            ["esl-sfetch", "-c", endpoints, fasta.name, expected_accession],
            text=True,
        )
        found = SeqIO.read(io.StringIO(output), "fasta")
        if found.seq != expected:
            logger.error("Found {}, but expected {}", found.seq, expected)
            return False
        return True


def map_accessions(
    filename: Path, fetcher: AccessionFetcher, start: str | None, delimiter: str
) -> None | tuple[int, dict[str, str]]:
    # TODO: Add tracking of seen accessions

    try:
        align = read_clean_alignment(filename)
    except Exception as e:
        logger.error("Error reading Stockholm file: {}", e)
        return None

    logger.info("Will try to map {}", len(align))
    mapping = {}
    for record in align:
        raw_accession = record.id
        accession = accession_splitter(start, delimiter, raw_accession)
        logger.info(
            "Fixing raw accession {}, extracted parent accession: {}",
            raw_accession,
            accession,
        )

        result = fetcher.get_sequence_and_accession(accession)
        if result is None:
            logger.warning(
                "Accession {} not downloaded and will be skipped", raw_accession
            )
            continue

        fasta_seq = result.seq
        sequence_accession = result.id
        assert sequence_accession, "Failed to find a sequence accession"
        fasta_seq = normalize_sequence(fasta_seq)
        dna_query = normalize_sequence(record.seq)

        for target_name, target, accession_builder in search_sequences(fasta_seq):
            logger.debug(
                "Looking for target {} in {} of {}",
                raw_accession,
                target_name,
                accession,
            )
            coordinates = coordinates_within(target, dna_query)
            if not coordinates:
                logger.debug(
                    "Target {} is not found in {} of {}",
                    raw_accession,
                    target_name,
                    accession,
                )
                continue

            new_accession = accession_builder(result, *coordinates)
            logger.info("{} corresponds to {}", raw_accession, new_accession)
            if not validate_accession(
                result,
                new_accession,
                dna_query,
            ):
                logger.warning("Failed to validate {}, skipping", new_accession)
                continue
            mapping[raw_accession] = new_accession
            break
        else:
            logger.warning("Could not correct {}, will skip", raw_accession)

    return (len(align), mapping)


def write_new_stockholm(
    alignment_file: Path, output: ty.IO, id_mapping: dict[str, str]
):
    with alignment_file.open("r") as handle:
        for line in handle:
            # Handle standard lines
            if (
                line.startswith("# STOCKHOLM 1.0")
                or line.startswith("//")
                or not line.strip()
            ):
                output.write(line)

            elif line.startswith("#="):
                match line[0:4]:
                    case "#=GF" | "#=GC":
                        logger.debug("Writing {} annotation line unchanged", line[0:4])
                        output.write(line)
                    case "#=GS" | "#=GR":
                        _, seqname, feature, value = re.split(r"\s+", line, 3)
                        if seqname not in id_mapping:
                            logger.warning(
                                "Did not find new accession for {} {}, skipping line",
                                line[0:4],
                                seqname,
                            )
                            continue
                        if feature == "AC":
                            logger.debug("Skipping accession line")
                            continue
                        converted = id_mapping[seqname]
                        logger.debug("Converting: {} -> {}", seqname, converted)
                        output.write(f"{line[0:4]} {converted: <30} {feature} {value}")
                    case a:
                        logger.warning("Not trying to convert annotation: {}", a)
                        output.write(line)

            else:
                line_acc, rest = re.split(r"\s+", line, 1)
                if line_acc not in id_mapping:
                    logger.warning(
                        "Did not find new accession for {}, skipping the line", line_acc
                    )
                    continue

                new_acc = id_mapping[line_acc]
                output.write(f"{new_acc: <30}")
                output.write(rest)


def accession_splitter(start: None | str, delimiter: str, accession: str) -> str:
    """
    Split the accession using the start and delimiter to find the expected
    delimiter.

    >>> accession_splitter(None, "/", "NR_00000.1/1-2")
    "NR_00000.1"
    >>> accession_splitter(None, "/", "NR_00000.1/1-2/2")
    "NR_00000.1"
    >>> accession_splitter(None, "/", "gen_NR_00000.1/1-2/2")
    "gen_NR_00000.1"
    >>> accession_splitter("_", "/", "gen_NR_00000.1/1-2/2")
    "NR_00000.1"
    >>> accession_splitter("_", ".+.", "gen_NR_00000.1/1-2/2")
    "NR_00000.1/1-2/2"
    """

    stripped = accession
    if start and start in stripped:
        _, stripped = stripped.split(start, 1)
    if delimiter in stripped:
        stripped, _ = stripped.split(delimiter, 1)
    return stripped


@click.command()
@click.argument(
    "input_file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    metavar="<stockholm_file>",
)
@click.argument(
    "output_file", type=click.File("w"), default="-", metavar="[output_file]"
)
@click.option(
    "-s",
    "--start",
    help="The string to indicate the start of the accession",
)
@click.option(
    "-d", "--delimiter", default="/", help="The accession endpoint delimiter to use"
)
@click.option(
    "-v", "--verbose", is_flag=True, help="Enable verbose logging (DEBUG level)"
)
@click.option("-q", "--quiet", is_flag=True, help="Suppress all output except errors")
@click.option(
    "--timeout",
    default=30,
    type=int,
    help="Timeout in seconds for NCBI requests (default: 30)",
)
@click.version_option(version="2.2.0", prog_name="find_start_stop_coordinates")
def main(input_file, output_file, verbose, quiet, timeout, delimiter, start):
    """
    Find start/stop coordinates for sequences in a Stockholm alignment file.

    This tool takes a Stockholm format alignment file containing INSDC or RefSeq
    accessions and maps each sequence to its genomic coordinates by fetching the
    complete sequences from NCBI and finding the alignment positions. This
    assumes that each sequence has an id in the file which contains a known
    INSDC or RefSeq accession. How the script finds it is is controlled by the
    `--delimiter` and `--start` options, see below for details.

    The alignment is passed through `esl-reformat --mingap -d -u` prior to
    being read so the resulting output should be clean, well-formatted, usable
    DNA alignment. But this requires that the input is a valid stockholm
    formatted file first.

    \b
    Arguments:
      INPUT_FILE    Stockholm format alignment file with valid accessions
      OUTPUT_FILE   Output file path, or "-" for stdout (default: stdout)

    \b
    Options:
      --delimiter, -d     The string that separates the NCBI accession from the rest of the id
      --start, -s         The string that indicates the prefix to the accession
      --quiet             Minimal logging information
      --timeout           Timeout for NCBI requests
      --verbose           Verbose logging information

    \b
    Examples:
      find_start_stop_coordinates.py input.sto                    # Output to stdout
      find_start_stop_coordinates.py input.sto output.sto         # Output to file
      find_start_stop_coordinates.py input.sto - --verbose        # Verbose stdout
      find_start_stop_coordinates.py input.sto - --quiet          # Quiet stdout

    \b
    Finding accessions:
      This assumes that each the id for each sequence starts with an NCBI
      accession. The accession is found by splitting on the `--delimiter`
      string (default `/`) and using everything before it. Some examples are
      (note, that empty string is indicated no entry in the start column):

    \b
      ID                    Start  Delimiter        Parent Accession
      --------------------  -----  ---------        ----------------
      NR_00000.1/1-2               /                NR_00000.1
      NR_00000.1/1-2               .                NR_00000
      NR_00000.1/1-2/1-2           /                NR_00000.1
      NR_00000.1/1-2/1-2           .                NR_00000
      NR_00000.1.+.1-2/1-2         .+.              NR_00000.1
      NR_00000.1.+.1-2/1-2         .                NR_00000
      NR_00000.1                   /                NR_00000.1
      NR_00000.1                   .                NR_00000
      gen_NR_00000.1/1-2    _      /                NR_00000.1
      gen_NR_00000.1/1-2    _      .                NR_00000
      gen_NR_00000.1/1-2    .      /                1
    """

    logger.remove()
    log_level = "INFO"
    if quiet:
        log_level = "ERROR"
    elif verbose:
        log_level = "DEBUG"

    logger.add(
        sys.stderr,
        format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | <level>{message}</level>",
        level=log_level,
    )

    fetcher = AccessionFetcher.build(timeout=timeout)
    input_file = Path(input_file)
    try:
        result = map_accessions(input_file, fetcher, start, delimiter)
        if not result:
            logger.error("No sequences could be mapped")
            raise ValueError("No sequences could be mapped")

        total, id_mapping = result
        buffer = io.StringIO()
        write_new_stockholm(input_file, buffer, id_mapping)
        sp.run(
            ["esl-reformat", "stockholm", "-"],
            check=True,
            stdout=output_file,
            input=buffer.getvalue(),
            text=True,
        )
        logger.info("Mapped {} sequences of {} total sequences", len(id_mapping), total)

    except KeyboardInterrupt:
        logger.info("Operation cancelled by user")
        raise click.Abort()
    except Exception as e:
        logger.error("Unexpected error: {}", e)
        if verbose:
            logger.exception("Full traceback:")
        raise click.ClickException(str(e))


if __name__ == "__main__":
    main()
