"""
Usage:

# navigate to the folder with the SEED and CM files
cd /path/to/folder/with/SEED

# identify invalid accessions that need to be replaced
rfblast.py validate

# manually upload the file `invalid.fa` to NCBI BLAST
# download results in the `Single-file JSON` format

# replace invalid accessions in the SEED file
rfblast.py replace XXXXXXXXXXX-Alignment.json

# where XXXXXXXXXXX-Alignment.json is the NCBI BLAST result.
"""

import json
import os
import re
import xmltodict

import click
import requests

from Bio.Blast import NCBIWWW
from Bio import SeqIO


IDENTITY = 90
QUERY_COVERAGE = 70
TEMPDIR = "temp"


def get_accession(gid):
    """
    Get versioned accession, for example:
    Input:
    gi|2047803076|gb|CP061286.1|
    Output:
    CP061286.1
    """
    parts = gid.split("|")
    return parts[-2]


def get_blast_data(filename):
    """
    Load BLAST JSON output.
    """
    with open(filename, "r") as f:
        return json.load(f)


def choose_replacement(data, min_identity, min_query_coverage):
    """
    Loop over BLAST results and pick best replacement for each hit.
    """
    # do not pick replacement from the same accession if already seen
    fasta = ""
    seen_accessions = set()
    for query_num, search in enumerate(data["BlastOutput2"]):
        query_title = search["report"]["results"]["search"]["query_title"]
        query_len = search["report"]["results"]["search"]["query_len"]
        replacement_found = False
        for entry in search["report"]["results"]["search"]["hits"]:
            acc = get_accession(entry["description"][0]["id"])
            if acc not in seen_accessions:
                seen_accessions.add(acc)
                replacement_found = True
            else:
                continue
            sequence = entry["hsps"][0]["hseq"]
            start = entry["hsps"][0]["hit_from"]
            end = entry["hsps"][0]["hit_to"]
            align_len = entry["hsps"][0]["align_len"]
            gaps = entry["hsps"][0]["gaps"]
            exact_matches = entry["hsps"][0]["identity"]
            identity = float(exact_matches) / align_len * 100
            query_coverage = float(align_len - gaps) / query_len * 100
            target_coverage = float(align_len - gaps) / len(sequence) * 100
            if identity >= min_identity and query_coverage >= min_query_coverage:
                warning = False
            else:
                warning = True
            summary = (
                "#{query_num} {message} {query_title} "
                "with {acc}/{start}-{end} at {identity}% identity; "
                "{gaps} gaps; query coverage {query_coverage}"
            ).format(
                acc=acc,
                start=start,
                end=end,
                query_title=query_title,
                identity=round(identity),
                query_coverage=round(query_coverage),
                target_coverage=round(target_coverage, 2),
                gaps=gaps,
                message="Replace"
                if not warning
                else "      WARNING: No replacement found for",
                query_num=query_num + 1,
            )
            print(summary)
            if not warning:
                fasta += ">{acc}/{start}-{end}\n{sequence}\n".format(
                    acc=acc,
                    start=start,
                    end=end,
                    sequence=sequence.replace("-", "").replace("T", "U"),
                )
            if replacement_found:
                break
    return fasta


def generate_new_seed(fasta):
    filename = "replacement.fasta"
    with open(filename, "w") as f:
        f.write(fasta)
    if not os.path.exists("CM"):
        cmd = "rfsearch.pl -t 30 -nodesc -relax"
        os.system(cmd)
    if not os.path.exists("SEED"):
        raise Exception("Error: SEED file not found")
    cmd = (
        "cmalign --mapali SEED --noprob CM {} > tempseed && "
        "esl-reformat pfam tempseed > NEWSEEDtemp && "
        "rm tempseed"
    ).format(filename)
    os.system(cmd)
    invalid = set()
    with open("invalid.txt", "r") as f:
        for line in f:
            invalid.add(line.strip())
    newseed = open("NEWSEED", "w")
    with open("NEWSEEDtemp", "r") as f:
        for line in f:
            skip = False
            for invalid_accession in invalid:
                if invalid_accession in line:
                    skip = True
                    break
            if not skip:
                newseed.write(line)
    newseed.close()
    os.remove("NEWSEEDtemp")
    cmd = (
        'echo "Old SEED info:" && esl-alistat SEED && '
        'echo "New SEED info:" && esl-alistat NEWSEED'
    )
    os.system(cmd)


def is_valid_accession(accession):
    """
    TB03JUN2009E__Contig_2000/988-772
    NZ_CP007501.1/771730-771924

    Found:
    {"header":{"type":"esearch","version":"0.3"},"esearchresult":{"count":"1","retmax":"1","retstart":"0","idlist":["1119664412"],"translationset":[],"querytranslation":""}}
    Found but a different ID:
    {"header":{"type":"esearch","version":"0.3"},"esearchresult":{"count":"1","retmax":"1","retstart":"0","idlist":["EP994606.1"],"translationset":[],"translationstack":[{"term":"JCVI_SCAF_1096627298421[All Fields]","field":"All Fields","count":"1","explode":"N"},"GROUP"],"querytranslation":"JCVI_SCAF_1096627298421[All Fields]"}}
    Not found:
    {"header":{"type":"esearch","version":"0.3"},"esearchresult":{"count":"0","retmax":"0","retstart":"0","idlist":[],"translationset":[],"querytranslation":"(TB03JUN2009E__Contig_2000[All Fields])","errorlist":{"phrasesnotfound":["TB03JUN2009E__Contig_2000"],"fieldsnotfound":[]},"warninglist":{"phrasesignored":[],"quotedphrasesnotfound":[],"outputmessages":["No items found."]}}}
    """
    parts = accession.split("/")
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term={}&retmode=json&idtype=acc"
    r = requests.get(url.format(parts[0]))
    data = r.json()
    if (
        int(data["esearchresult"]["count"]) == 1
        and len(data["esearchresult"]["idlist"]) > 0
        and data["esearchresult"]["idlist"][0] == parts[0]
    ):
        return True
    else:
        return False


def fetch_seqs(filename, accessions, label):
    """ """
    f_txt = "{}.txt".format(label)
    f_fa = "{}.fa".format(label)
    with open(f_txt, "w") as f:
        for accession in accessions:
            f.write(accession + "\n")
    cmd = "esl-sfetch -f {} {} > {}".format(filename, f_txt, f_fa)
    os.system(cmd)
    click.echo("Saved {} {} accessions in {}".format(len(accessions), label, f_fa))


def parse_fasta(filename):
    """ """
    accessions = set()
    with open(filename, "r") as f:
        for line in f:
            if line.startswith(">"):
                parts = re.split(r"\s+", line)
                accession = parts[0].replace(">", "")
                accessions.add(accession)
    valid = set()
    invalid = set()
    for accession in accessions:
        if is_valid_accession(accession):
            click.echo("{} is valid".format(accession))
            valid.add(accession)
        else:
            click.echo("{} is invalid".format(accession))
            invalid.add(accession)
    os.system("esl-sfetch --index {}".format(filename))
    if valid:
        fetch_seqs(filename, valid, "valid")
    else:
        click.echo("No valid accessions found")
    if invalid:
        fetch_seqs(filename, invalid, "invalid")
        click.echo("=========================\nGenerated file invalid.fa")
        os.system("esl-seqstat {}".format("invalid.fa"))
        click.echo("Upload file invalid.fa to NCBI BLAST")
    else:
        click.echo("No invalid accessions found")


@click.group()
def cli():
    pass


def validate(seed):
    """
    Convert SEED to a fasta file containing sequences with unknown IDs.
    """
    if not os.path.exists(seed):
        raise Exception("Error: SEED does not exist")
    fasta = "seed.fasta"
    cmd = "esl-reformat fasta {} > {}".format(seed, fasta)
    os.system(cmd)
    print(fasta)
    parse_fasta(fasta)


@cli.command()
@click.argument("invalid", type=click.Path(exists=True))
def blast_invalid_sequences(invalid, blast_program="blastn"):
    """
    Upload the file `invalid.fa` to NCBI BLAST and download results in the `Single-file JSON` format
    """

    blast_results = []
    fasta_file = "invalid.fa"
    sequences = SeqIO.parse(fasta_file, "fasta")

    for seq in sequences:
        result_handle = NCBIWWW.qblast("blastn", "nt", seq.seq)

        blast_results.append(result_handle.read())

    output_file = "blast_results.xml"
    with open(output_file, "w") as output:
        for res in blast_results:
            output.write(res)

    with open(output_file, "r") as xml_file:
        xml_data = xml_file.read()

    xml_dict = xmltodict.parse(xml_data)
    json_data = json.dumps(xml_dict, indent=4)

    with open("output.json", "w") as json_file:
        json_file.write(json_data)

    print(f"BLAST search completed. Results saved to {json_file}")


@cli.command()
@click.argument("blast_json", type=click.Path(exists=True))
@click.option(
    "--identity",
    default=IDENTITY,
    type=click.FLOAT,
    help="Minimum % identity between query and target",
)
@click.option(
    "--query_coverage",
    default=QUERY_COVERAGE,
    type=click.FLOAT,
    help="Minimum coverage of the seed sequence",
)
def replace(blast_json, identity, query_coverage):
    """
    Replace unknown accessions in SEED alignment using NCBI BLAST results
    """
    blast_data = get_blast_data(blast_json)
    fasta = choose_replacement(blast_data, identity, query_coverage)
    generate_new_seed(fasta)
    click.echo("Done")


@cli.command()
@click.argument("seed", type=click.Path(exists=True))
def all(seed):
    """
    Carry out all steps
    """
    validate(seed)
    blast_result = blast_invalid_sequences("invalid.fa")
    print(blast_result)
    replace(blast_json=blast_result)


if __name__ == "__main__":
    cli()
