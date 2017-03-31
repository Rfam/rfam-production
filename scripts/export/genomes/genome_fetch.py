# !/usr/bin/python
"""
Copyright [2009-2016] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
     http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

"""
TO DO:    - logging
"""

# ---------------------------------IMPORTS-------------------------------------


import httplib
import json
import os
import string
import urllib
import urllib2
import xml.etree.ElementTree as ET

import requests
from rdflib import Graph

from config import gen_config as gc

# ----------------------------------GLOBALS------------------------------------

# URLs
# This returns the list of reference proteomes from Uniprot
REF_PROT_LIST_URL = gc.REF_PROT_LIST_URL

# Retrieve the proteome rdf file
PROTEOME_URL = gc.PROTEOME_URL

# Retrieve the genome's xml file
ENA_XML_URL = gc.ENA_XML_URL

# ENA url for file download
ENA_DATA_URL = gc.ENA_DATA_URL

# ENA url for assembly data retrieval via taxon id
ENA_TAX_URL = gc.ENA_TAX_URL

# ENA GCA report file label
GCA_REP_LABEL = gc.GCA_REP_LBL

# NCBI URL for sequence download
NCBI_SEQ_URL = gc.NCBI_SEQ_URL

# ENA file formats
FORMATS = {"xml": ".xml", "fasta": ".fa"}

# ---------------------------------------------------------------------- #STEP1


def fetch_ref_proteomes():
    """
    This method returns a list of all reference proteome accessions available
    from Uniprot
    """

    ref_prot_list = []
    response = urllib2.urlopen(REF_PROT_LIST_URL)

    for ref_prot in response:
        ref_prot_list.append(ref_prot.strip())

    return ref_prot_list


# -----------------------------------------------------------------------------


def export_gca_accessions(upid_gca):
    """
    Retrieves reference proteomes ids and their associated gca accessions
    as well as the taxonomic rank/domain (eukaryotes, bacteria etc)

    upid_gca: Uniprot's tab separated file (UPID_GCA.tsv)
    """

    # need to check if the path provided is a valid file
    upid_gca_fp = open(upid_gca, 'r')

    prot_gca_pairs = {}

    for prot in upid_gca_fp:
        prot = prot.strip().split('\t')
        if prot[1] != '':
            prot_gca_pairs[prot[0]] = prot[1]
        else:
            prot_gca_pairs[prot[0]] = -1

    upid_gca_fp.close()

    return prot_gca_pairs


# -----------------------------------------------------------------------------


def extract_genome_acc(prot_rdf):
    """
    Extracts and returns the assembly accession from the proteome rdf
    which provided as input. Returns -1 if not available

    prot_rdf: A Uniprot's proteome rdf url or file path
    """
    g = Graph()

    response = requests.get(prot_rdf).status_code

    if response == httplib.OK:

        g.load(prot_rdf)

        for s, p, o in g:
            if string.find(o, "GCA") != -1:
                return os.path.split(o)[1]

    return -1


# -----------------------------------------------------------------------------


def proteome_rdf_scanner(proteome):
    """
    Scans a Uniprot's reference proteome rdf file and looks for all
    available accessions. Returns a dictionary with GCA and WGS accessions
    where applicable

    prot_rdf: Uniprot's proteome rdf url or file path
    """

    # need to do some http error handling here and if resource is unavailable
    # return None and or http error code

    prot_rdf = PROTEOME_URL % proteome

    accessions = {"GCA": -1, "WGS": -1}
    wgs_flag = False
    g = Graph()

    if requests.get(prot_rdf).status_code == httplib.OK:
        g.load(prot_rdf)

        # scan for accessions
        for s, p, o in g:
            # look for ENA accessions
            if string.find(o, "/embl/") != -1:
                if string.find(o, "GCA") != -1:
                    accessions["GCA"] = os.path.split(o)[1]

                elif wgs_flag is True:
                    accessions["WGS"] = os.path.split(o)[1]

            # if WGS keyword found, set flag to true
            elif string.find(o, "WGS") != -1:
                wgs_flag = True
    else:
        pass

    return accessions


# ---------------------------------------------------------------------- #STEP2


def fetch_genome_acc(prot):
    """
    Returns a proteome's corresponding assembly accession (ENA) in a
    dictionary format {proteome_acc:gca_acc}


    prot: One of (file|list|acc)
          - file: Uniprot's proteome list file
          - list: a list of reference proteome accessions (fetch_ref_proteomes)
          - acc: A single Uniprot ref. proteome accession
    """

    gen_acc = None
    ref_prot_list = None
    gens = {}

    # ref proteome accession list
    if isinstance(prot, list):
        ref_prot_list = prot

    # ref proteome list file
    elif os.path.isfile(prot):
        prot_file = open(prot, 'r')
        ref_prot_list = prot_file.readlines()
        prot_file.close()

    # single ref proteome accession
    else:
        proteome = string.strip(prot)
        rdf_url = gc.PROTEOME_URL % (proteome)
        res_handle = urllib.urlopen(rdf_url)
        if res_handle.getcode() == httplib.OK:
            gen_acc = extract_genome_acc(rdf_url)
        gens[proteome] = gen_acc

        return gens

    # do this for files and lists
    for proteome in ref_prot_list:
        proteome = string.strip(proteome)
        rdf_url = gc.PROTEOME_URL % (proteome)
        res_handle = urllib.urlopen(rdf_url)
        if res_handle.getcode() == httplib.OK:
            gen_acc = extract_genome_acc(rdf_url)

        gens[proteome] = gen_acc

    return gens


# -----------------------------------------------------------------------------


def fetch_ena_file(acc, file_format, dest_dir):
    """
    Retrieves a file given a valid ENA accession and stores it in the
    indicated destination in the selected format

    acc: A valid ENA entry accession
    format: A valid ENA file format (xml, fasta, txt)
    dest_dit: A valid path to destination directory
    """

    seq_url = None
    file_path = None

    if file_format.find("xml") != -1:
        seq_url = ENA_XML_URL % acc
        file_path = os.path.join(dest_dir, acc + FORMATS[file_format])

    # fetching compressed file
    else:
        seq_url = ENA_DATA_URL % (acc, file_format)
        file_path = os.path.join(dest_dir, acc + FORMATS[file_format] + ".gz")

    urllib.urlretrieve(seq_url, file_path)

    if os.path.exists(file_path):
        return True

    return False


# ---------------------------------------------------------------------- #STEP3


def extract_assembly_accs(accession):
    """
    Loads an xml tree from a file or a string (usually an http response),
    and returns a list with the genome assembly's chromosomes

    accession: A valid ENA assembly accession (without the assembly version)
    """

    accessions = []
    root = None
    chroms = None
    assembly_link = None
    assembly = None

    assembly_xml = requests.get(ENA_XML_URL % accession).content

    if os.path.isfile(assembly_xml):
        # parse xml tree and return root node
        root = ET.parse(assembly_xml).getroot()
    else:
        # fromstring returns the xml root directly
        root = ET.fromstring(assembly_xml)

    assembly = root.find("ASSEMBLY")

    if assembly is not None:

        # either parse the assembly report file or get the WGS range
        assembly_link = assembly.find(
            "ASSEMBLY_LINKS")

        if assembly_link is not None:
            # export url link and fetch all relevant assembly accessions
            url_link = assembly_link.find(
                "ASSEMBLY_LINK").find("URL_LINK").find("URL").text

            # need to check for GCA_REP_LABEL
            accessions = assembly_report_parser(url_link)

        else:
            # no assembly link provided - look for WGS element
            wgs = None
            wgs = assembly.find("WGS_SET")

            # get wgs accession
            wgs_acc = get_wgs_set_accession(
                wgs.find("PREFIX").text, wgs.find("VERSION").text)

            # get wgs range and return as single accession
            accessions.append(get_wgs_range(wgs_acc))

    # move this outside this function to where accessions are requested ??
    else:
        accessions.append(get_wgs_range(accession))

    return accessions


# ---------------------------------------------------------------------- #STEP4


def download_genomes(gen, dest_dir):
    """
    Downloads all chromosome files of a given assembly accession (ENA) in
    dest_dir

    gen: Single accession (string) or a list of genome accessions (GC*)
    dest_dir: The path of the destination directory to export the fasta
    files
    """

    # need to add logging

    accessions = None

    if os.path.isfile(gen):
        fp = open(gen, 'r')

        for gen_acc in fp:
            gen_acc = string.strip(gen_acc)
            if string.find(gen_acc, '.') != -1:
                gen_acc = gen_acc.partition('.')
                gen_acc = gen_acc[0]

            gen_dir = os.path.join(dest_dir, gen_acc)

            try:
                os.mkdir(gen_dir)
            except:
                print "Unable to generate directory for accession: ", gen

            accessions = extract_assembly_accs(gen_acc)

            if len(accessions) > 0:
                for acc in accessions:
                    fetch_ena_file(acc, "fasta", gen_dir)

            gen_acc = None
            gen_dir = None
            accessions = None

    else:
        if string.find(gen, '.') != -1:
            gen = gen.partition('.')
            gen = gen[0]

        gen_dir = os.path.join(dest_dir, gen)

        try:
            os.mkdir(gen_dir)
        except:
            print "Unable to generate directory for accession: ", gen

        accessions = extract_assembly_accs(gen)

        if len(accessions) > 0:
            for acc in accessions:
                fetch_ena_file(acc, "fasta", gen_dir)

        # if no accessions found, write to log file
        else:
            return None

    return True


# -----------------------------------------------------------------------------

def fetch_genome(gen, dest_dir):
    """
    Downloads and parses xml file of the given genome accession (gen), and
    downloads all chromosome files in fasta format in destination directory
    (dest_dir). The xml file is deleted after completion.

    gen: ENA assembly accession (GCA*)
    dest_dir: Destination of the output directory
    """

    gen_dir = os.path.join(dest_dir, gen.partition('.')[0])

    try:
        os.mkdir(gen_dir)
    except:
        pass

    fetch_ena_file(gen, "xml", gen_dir)

    gen_xml = os.path.join(gen_dir, gen + ".xml")
    chroms = extract_assembly_accs(gen_xml)

    for chrom in chroms:
        fetch_ena_file(chrom, "fasta", gen_dir)

    os.remove(gen_xml)  # remove xml file when done


# -----------------------------------------------------------------------------


def rdf_accession_search(ref_prot_acc, sub_str):
    """
    Parses rdf url and returns a list of ENA accessions

    rdf_url: The url to a Uniprot's reference proteome rdf url
    sub_str: A sub string to look for in the rdf file
    """

    accessions = []
    rdf_graph = Graph()
    rdf_url = PROTEOME_URL % ref_prot_acc

    response = requests.get(rdf_url).status_code

    if response == httplib.OK:

        rdf_graph.load(rdf_url)

        for s, p, o in rdf_graph:
            if string.find(o, sub_str) != -1:
                accessions.append(os.path.split(o)[1])
    else:
        # return http status code
        # return response.status_code
        pass

    return accessions


# -----------------------------------------------------------------------------


def assembly_report_parser(report_url):
    """
    Parses an assembly report file and returns a list of all available
    accessions (scaffolds, contigs etc)

    report_url: A url provided within an ENA assembly xml file. This is the
    text of URL tag under ASSEMBLY/ASSEMBLY_LINKS/ASSEMBLY_LINK.By default this
    is an ftp request url. Converting to http to fetch assembly accessions.
    """

    accessions = []

    # switch from ftp to http to fetch assembly accessions on the go

    report_url = report_url.replace("ftp://", "http://")

    # fetch assembly report file contents and store in a list, omitting header
    ass_rep_file = requests.get(report_url).content.split('\n')[1:]

    # if empty line, remove it
    if ass_rep_file[len(ass_rep_file) - 1] == '':
        ass_rep_file.pop(len(ass_rep_file) - 1)

    # parse list and export assembly accessions
    for line in ass_rep_file:
        line = line.strip().split('\t')
        if line[0].find('.') != -1:
            accessions.append(line[0])

    return accessions


# -----------------------------------------------------------------------------


def get_wgs_set_accession(prefix, version):
    """
    Generates ENA WGS accession using the WGS accession pattern
    (4 letter prefix) and 2-digit build version. The WGS accession is
    generated by appending prefix and version with a postfix of 6-zeros.
    For more information please visit ENA service-news: http://goo.gl/LnIyQ3

    prefix: A 4-char string representing the ENA WGS accession
    version: A 2-digit representing the WGS build version
    """

    postfix = "000000"
    wgs_accession = None

    if int(version) % 10 != 0:
        wgs_accession = prefix + '0' + version + postfix
    else:
        wgs_accession = prefix + version + postfix

    return wgs_accession


# -----------------------------------------------------------------------------


def get_wgs_range(wgs_acc):
    """
    Fetches the wgs related xml file from ENA and exports the wgs range

    wgs_acc: A valid ENA wgs accession
    """

    wgs_range = None
    wgs_xml_str = requests.get(ENA_XML_URL % wgs_acc).content
    wgs_xml_root = ET.fromstring(wgs_xml_str)

    if wgs_xml_root.find("entry") is not None:

        wgs_xrefs = wgs_xml_root.find("entry").findall("xref")

        for xref_el in wgs_xrefs:
            if xref_el.get("db") == "ENA-WGS":
                wgs_range = xref_el.get("id")

    return wgs_range


# -----------------------------------------------------------------------------


def lsf_cmd_generator(upid, gca_acc, domain, exec_path, proj_dir):
    """
    Generates an lsf job command for downloading a new genome. Returns an
    LSF specific bsub command

    upid: Uniprot's reference proteome id
    gca_acc: ENA's genome accession. -1 if there's no available id
    domain: Proteome's taxonomic domain
    exec_path: The path to the pipeline executable
    proj_dir: The path to the project directory
    """

    prot_dir = os.path.join(os.path.join(proj_dir, domain), upid)

    cmd = ("bsub -M %s "
           "-R \"rusage[mem=%s,tmp=%s]\" "
           "-o \"%s\" "
           "-e \"%s\" "
           "-u \"%s\" "
           "-Ep \"rm -rf luigi\" "
           "-g %s/%s "
           "python %s DownloadGenome --upid %s --gca-acc %s --project-dir %s --domain %s") % (
               gc.MEM, gc.MEM, gc.TMP_MEM,
               os.path.join(prot_dir, "download.out"),
               os.path.join(prot_dir, "download.err"),
               gc.USER_EMAIL, gc.LSF_GEN_GROUP,
               domain, exec_path,
               upid, gca_acc,
               proj_dir, domain)

    return cmd

# -----------------------------------------------------------------------------


def genome_script_generator(upid, domain, gen_size, out_dir):
    """
    Generates a shell script for a proteome with ip upid under out_dir.
    Memory is reserved according to genome size

    upid: Uniprot's unique proteome id
    domain: The domain under which a proteome has been classified
    gen_size: The genome's size
    out_dir: Destination directory
    """

    # create shell script for upid within
    fp = open(os.path.join(out_dir, upid + ".sh"), 'w')

    # mem = gen_size * something might not need these for downloads
    mem_size = 8000
    tmp_size = gen_size * 2
    tmp_dir = "/tmp/%s_$LSB_JOBID" % (upid)

    # generate proteome destination directory
    prot_dest_dir = os.path.join(
        os.path.join(os.path.split(out_dir)[0], domain), upid)

    fp.write("#!/bin/csh\n")
    fp.write("#BSUB -M %s\n" % mem_size)
    fp.write("#BSUB -R \"rusage[mem=%s,tmp=%s]\"\n" % (mem_size, tmp_size))

    # create a directory using the proteomes unique id extended by jobid
    fp.write("#BSUB -E \"mkdir -m 777 -p %s\"\n" % tmp_dir)
    fp.write("#BSUB -o \"%s/%sJ.out\"\n" % (tmp_dir, chr(37)))
    fp.write("#BSUB -e \"%s/%sJ.err\"\n" %
             (tmp_dir, chr(37)))  # just the error output

    fp.write("#BSUB -u \"%s\"\n" % gc.USER_EMAIL)  # email this user

    # need to write files back to genome dir prot_dest_dir
    fp.write(
        "#BSUB -f \"%s/download.out < /tmp/%sJ/%sJ.out\"\n" % (prot_dest_dir,
                                                               chr(37),
                                                               chr(37)))

    fp.write(
        "#BSUB -f \"%s/download.err < /tmp/%sJ/%sJ.err\"\n" % (prot_dest_dir,
                                                               chr(37),
                                                               chr(37)))

    # delete everything on termination or completion of job
    fp.write("#BSUB -Ep \"rm -rf %s\"\n" % tmp_dir)

    fp.write("#BSUB -g %s/%s \n\n" % (gc.LSF_GEN_GROUP % domain))

    # call executable
    fp.write("python %s %s %s \n\n" %
             (gc.GEN_DWLD_EXEC, os.path.join(prot_dest_dir, upid + ".json"),
              prot_dest_dir))

    # copy files to destination
    fp.write("cp %s/*.gz %s/.\n" % (tmp_dir, prot_dest_dir))


# -----------------------------------------------------------------------------


def load_upid_gca_file(upid_gca_file):
    """
    Parses Uniprot's upid tsv file and exports all important information in json
    format
    """

    upid_gca_dict = {}
    upid_fp = open(upid_gca_file, 'r')

    for upid_line in upid_fp:
        upid_line = upid_line.strip().split('\t')

        # add GCA accession
        if upid_line[1] != '':
            upid_gca_dict[upid_line[0]] = {"GCA": upid_line[1]}
        else:
            upid_gca_dict[upid_line[0]] = {"GCA": -1}

        upid_gca_dict[upid_line[0]]["DOM"] = upid_line[2]

    upid_fp.close()

    return upid_gca_dict


# -----------------------------------------------------------------------------


def load_upid_gca_pairs():
    """
    This is an alternative version to load_upid_gca_file loading the pairs from
    Uniprot's REST API. Returns a dictionary of upid, gca accession pairs,
    including the species kingdom
    """

    id_pairs = {}
    response = requests.get(gc.REF_PROT_REST_URL)

    if response.status_code == 200:
        content = response.content
        prot_lines = content.split('\n')

        # remove header line
        prot_lines.pop(0)

        for prot_line in prot_lines:
            if prot_line != '':
                prot_line = prot_line.strip().split('\t')
                tax_str = prot_line[2].split(',')[0].lower().strip()
                gca_acc = prot_line[1]
                if gca_acc != '':
                    id_pairs[prot_line[0]] = {
                        "GCA": prot_line[1], "DOM": tax_str}
                else:
                    id_pairs[prot_line[0]] = {"GCA": -1, "DOM": tax_str}
    else:
        # raise an error here
        pass

    return id_pairs


# -----------------------------------------------------------------------------


def fetch_genome_accessions(upid, gca_acc):
    """
    Fetches and returns a list of all accessions for a specific ref. proteome

    upid: Uniprot's ref. proteome id
    gca_acc: An ENA GCA accession associated with the upid (if available or -1)
    """

    gen_accs = []
    gca_acc = str(gca_acc)

    # there's a GCA accession
    if gca_acc != "-1":
        gca_acc = gca_acc.split('.')[0]
        gen_accs = extract_assembly_accs(gca_acc)

    else:
        prot_accs = proteome_rdf_scanner(upid)

        # no GCA or WGS, get any accessions from proteome rdf
        if prot_accs["GCA"] == -1 and prot_accs["WGS"] == -1:
            gen_accs = rdf_accession_search(upid, "/embl/")

        # found a GCA accession in the rdf file
        elif prot_accs["GCA"] != -1 and prot_accs["WGS"] == -1:
            gen_accs = extract_assembly_accs(prot_accs["GCA"])

        # WGS found
        elif prot_accs["GCA"] == -1 and prot_accs["WGS"] != -1:
            # call get_wgs_range directly here
            gen_accs = extract_assembly_accs(prot_accs["WGS"])

    # for all other cases this function will return an empty list

    return gen_accs


# -----------------------------------------------------------------------------


def fetch_wgs_range_accs(wgs_range):
    """
    Splits the WGS range into distinct accessions for metadata retrieval

    wgs_range: A valid ENA-WGS set range
    """

    wgs_accs = []

    wgs_end_points = wgs_range.strip().split('-')
    wgs_prefix = wgs_end_points[0][0:5]

    wgs_start = int(wgs_end_points[0][5:])
    wgs_end = int(wgs_end_points[1][5:])

    wgs_acc = ''

    while wgs_start < wgs_end:
        wgs_acc = wgs_prefix + str(wgs_start)
        wgs_accs.append(wgs_acc)
        wgs_start += 1
        wgs_acc = ''

    # include the last accession
    wgs_accs.append(wgs_end_points[1])

    return wgs_accs


# -----------------------------------------------------------------------------


def genome_download_validator(genome_dir):
    """
    Loop over Genome Download output directory and report any upids with
    erroneous output

    genome_dir: The path to Genome Download output directory
    """

    erroneous_genomes = {}

    # list all kingdom dirs under genome output directory
    project_dirs = os.listdir(genome_dir)
    # filter out items that are not directories
    kingdom_dirs = [x for x in project_dirs if os.path.isdir(os.path.join(genome_dir, x))]

    for kingdom in kingdom_dirs:
        erroneous_genomes[kingdom] = []
        # list all genome directories per kindom
        kingdom_dir_loc = os.path.join(genome_dir, kingdom)
        kingdom_dirs = os.listdir(kingdom_dir_loc)

        for genome in kingdom_dirs:
            genome_dir_loc = os.path.join(kingdom_dir_loc, genome)

            if os.path.exists(os.path.join(genome_dir_loc, "download.out")):
                fp = open(os.path.join(genome_dir_loc, "download.out"))
                for line in fp.readlines():
                    if line.find("Success") != -1:
                        success = 1
                        break
                fp.close()

                if success == 0:
                    erroneous_genomes[kingdom].append(genome)

            success = 0

    fp_out = open(os.path.join(genome_dir, "download_report.txt"), 'w')

    # print out erroneous genomes
    for kingdom in erroneous_genomes.keys():
        if len(erroneous_genomes[kingdom]) > 0:
            fp_out.write(kingdom + '\n')
            for proteome in erroneous_genomes[kingdom]:
                fp_out.write(proteome + '\n')

            fp_out.write('\n')

    fp_out.close()


# -----------------------------------------------------------------------------

def download_fasta_from_ncbi(accession, dest_dir):
    """
    Download fasta sequences from NCBI. In case of ENA obsolete sequences use
    this function to download the relevant files

    accession: A genome accession to download
    dest_dir: Destination directory to save the file to

    return: True on success, otherwise False
    """
    seq_url = None
    file_path = None

    seq_url = NCBI_SEQ_URL % (accession)
    file_path = os.path.join(dest_dir, accession + '.fa')

    urllib.urlretrieve(seq_url, file_path)

    if os.path.exists(file_path):
        return True

    return False


# -----------------------------------------------------------------------------

def download_sequence_report_files(project_dir, upid_gca_file):
    """
    Loads upid_gca_file json file and downloads from ENA all sequence report
    files per GCA accession. Skips if no GCA accession available

    project_dir: The path to a project directory as generated by Genome
    Download pipeline (genome_downloader.py)
    upid_gca_file: upid_gca file in json format as generated by the Genome
    download pipeline (genome_downloader.py)

    returns: void
    """

    err_seq_rep_files = {}

    fp = open(upid_gca_file, 'r')
    acc_pairs = json.load(fp)
    fp.close()

    for upid in acc_pairs.keys():

        if acc_pairs[upid]["GCA"] != -1:
            domain_dir = os.path.join(project_dir, acc_pairs[upid]["DOM"])
            if not os.path.exists(domain_dir):
                os.mkdir(domain_dir)

            updir = os.path.join(domain_dir, upid)

            if not os.path.exists(updir):
                os.mkdir(updir)

            seq_rep_url = gc.SEQ_REP_URL_TEMPLATE % (acc_pairs[upid]["GCA"][0:7],
                                                     acc_pairs[upid]["GCA"][0:10],
                                                     acc_pairs[upid]["GCA"])

            filename = "%s_sequence_report.txt" % acc_pairs[upid]["GCA"]

            urllib.urlretrieve(seq_rep_url, os.path.join(updir, filename))

            # check file exists or if it is empty
            if not os.path.exists(filename):
                err_seq_rep_files[upid] = {"GCA": acc_pairs[upid]["GCA"],
                                           "DOM": acc_pairs[upid]["DOM"]}
                continue

            elif os.path.getsize(filename) == 0:
                err_seq_rep_files[upid] = {"GCA": acc_pairs[upid]["GCA"],
                                           "DOM": acc_pairs[upid]["DOM"]}

    if len(err_seq_rep_files.keys()) > 0:
        fp_out = open(project_dir, "err_seq_rep_files.json", 'w')
        json.dump(err_seq_rep_files, fp_out)
        fp_out.close()

# -----------------------------------------------------------------------------

def sequence_report_to_json(seq_report_file, dest_dir=None):
    """
    Convert a GCA sequence report file (ENA) from .txt to .json format

    seq_report_file: The path to a valid GCA related sequence report file
    dest_dir: The path to destination directory. If None use the directory
    of the input file

    return: Accession dictionary
    """

    acc_dict = {}
    fp = open(seq_report_file, 'r')
    # discard header line
    fp.readline()

    for line in fp:
        line = line.strip().split('\t')

        acc_dict[line[0]] = {"sequence_name": line[1], "sequence-length": line[2],
                             "sequence-role": line[3], "replicon-name": line[4],
                             "replicon-type": line[5], "assembly-unit": line[6]}
    fp.close()

    if dest_dir is None:
        dest_dir = os.path.split(seq_report_file)[0]

    filename = os.path.basename(seq_report_file).partition(".")[0]
    fp_out = open(os.path.join(dest_dir, filename+".json"), 'w')
    json.dump(acc_dict, fp_out)
    fp_out.close()

    return acc_dict

# -----------------------------------------------------------------------------


if __name__ == '__main__':
    pass
