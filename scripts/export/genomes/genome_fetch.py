# !/usr/bin/python
"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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

# TO DO:    - logging

# ---------------------------------IMPORTS-------------------------------------

import httplib
import json
import os
import sys
import string
import urllib
import urllib2
import shutil
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
PROTEOME_XML_URL = gc.PROTEOME_XML_URL

# Retrieve the genome's xml file
ENA_XML_URL = gc.ENA_XML_URL

# ENA url for file download
ENA_DATA_URL = gc.ENA_DATA_URL
ENA_DATA_URL_GZIP = gc.ENA_DATA_URL_GZIP

# ENA url for assembly data retrieval via taxon id
ENA_TAX_URL = gc.ENA_TAX_URL

# ENA GCA report file label
GCA_REP_LABEL = gc.GCA_REP_LBL

# NCBI URL for sequence download
NCBI_SEQ_URL = gc.NCBI_SEQ_URL

# ENA file formats
FORMATS = {"xml": ".xml", "fasta": ".fa"}

# Maximum sequences per file
MAX_SEQS = 100000


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

def proteome_xml_scanner(proteome):
    """
    Scans a Uniprot's reference proteome rdf file and looks for all
    available accessions. Returns a dictionary with GCA and WGS accessions
    where applicable

    prot_xml: Uniprot's proteome rdf url or file path
    """

    # need to do some http error handling here and if resource is unavailable
    # return None and or http error code

    prot_xml = PROTEOME_XML_URL % proteome
    prefix = "{http://uniprot.org/uniprot}%s"

    accessions = {"GCA": -1, "WGS": -1}

    if requests.get(prot_xml).status_code == httplib.OK:
        xml_root = ET.fromstring(requests.get(prot_xml).content)
        proteome = xml_root.find(prefix % "proteome")

        # look for a GCA accession
        gen_assembly = None
        gen_assembly = proteome.find(prefix % "genome_assembly")
        if gen_assembly is not None:
            accessions["GCA"] = gen_assembly.text

        prot_components = proteome.findall(prefix % "component")

        # scan WGS accession
        for component in prot_components:
            component_name = component.get("name")
            if component_name.find("WGS") != -1:
                accessions["WGS"] = component.find(prefix % "genome_accession").text
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

def fetch_ena_file(acc, file_format,dest_dir, compressed=True):
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
        if compressed is True:
            seq_url = ENA_DATA_URL_GZIP % (acc, file_format)
            file_path = os.path.join(dest_dir, acc + FORMATS[file_format] + ".gz")
        else:
            seq_url = ENA_DATA_URL % (acc, file_format)
            file_path = os.path.join(dest_dir, acc + FORMATS[file_format])

    urllib.urlretrieve(seq_url, file_path)

    if os.path.exists(file_path):
        return True

    return False


# ---------------------------------------------------------------------- #STEP3

def extract_assembly_accs(accession):
    """
    Loads an xml tree from a file or a string (usually an http response),
    and returns a list with the genome assembly's chromosomes

    accession: A valid ENA GCA accession (without the assembly version)
    """

    accessions = []
    root = None
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
            accessions = assembly_report_parser(url_link, url=True)

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
        gen_fp = open(gen, 'r')

        for gen_acc in gen_fp:
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

        gen_fp.close()

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

def assembly_report_parser(assembly_report, url=True):
    """
    Parses an assembly report file and returns a list of all available
    accessions (scaffolds, contigs etc)

    report_url: A url provided within an ENA assembly xml file. This is the
    text of URL tag under ASSEMBLY/ASSEMBLY_LINKS/ASSEMBLY_LINK.By default this
    is an ftp request url. Converting to http to fetch assembly accessions.
    """

    accessions = []

    # switch from ftp to http to fetch assembly accessions on the go
    report_url = None
    ass_rep_file = None

    if url is True:
        report_url = assembly_report.replace("ftp://", "http://")
        # fetch assembly report file contents and store in a list, omitting header
        ass_rep_file = requests.get(report_url).content.split('\n')[1:]

        # if empty line, remove it
        if ass_rep_file[len(ass_rep_file) - 1] == '':
            ass_rep_file.pop(len(ass_rep_file) - 1)
    else:
        ass_rep_file = open(assembly_report, 'r')

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

    response = requests.get(ENA_XML_URL % wgs_acc)

    if response.status_code == httplib.OK:
        wgs_xml_str = response.content
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

    subdir_idx = upid[8:]
    prot_dir = os.path.join(os.path.join(proj_dir, subdir_idx), upid)

    cmd = ("bsub -M %s "
           "-R \"rusage[mem=%s,tmp=%s]\" "
           "-o \"%s\" "
           "-e \"%s\" "
           "-u \"%s\" "
           "-Ep \"rm -rf luigi\" "
           "-g %s "
           "python %s DownloadGenome --upid %s --gca-acc %s --project-dir %s --domain %s") % (
               gc.MEM, gc.MEM, gc.TMP_MEM,
               os.path.join(prot_dir, "download.out"),
               os.path.join(prot_dir, "download.err"),
               gc.USER_EMAIL, gc.LSF_GEN_GROUP,
               exec_path, upid, gca_acc,
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
    shell_fp = open(os.path.join(out_dir, upid + ".sh"), 'w')

    # mem = gen_size * something might not need these for downloads
    mem_size = 8000
    tmp_size = gen_size * 2
    tmp_dir = "/tmp/%s_$LSB_JOBID" % (upid)

    # generate proteome destination directory
    prot_dest_dir = os.path.join(
        os.path.join(os.path.split(out_dir)[0], domain), upid)

    shell_fp.write("#!/bin/csh\n")
    shell_fp.write("#BSUB -M %s\n" % mem_size)
    shell_fp.write("#BSUB -R \"rusage[mem=%s,tmp=%s]\"\n" % (mem_size, tmp_size))

    # create a directory using the proteomes unique id extended by jobid
    shell_fp.write("#BSUB -E \"mkdir -m 777 -p %s\"\n" % tmp_dir)
    shell_fp.write("#BSUB -o \"%s/%sJ.out\"\n" % (tmp_dir, chr(37)))
    shell_fp.write("#BSUB -e \"%s/%sJ.err\"\n" % (tmp_dir, chr(37)))

    shell_fp.write("#BSUB -u \"%s\"\n" % gc.USER_EMAIL)  # email this user

    # need to write files back to genome dir prot_dest_dir
    shell_fp.write(
        "#BSUB -f \"%s/download.out < /tmp/%sJ/%sJ.out\"\n" % (prot_dest_dir,
                                                               chr(37),
                                                               chr(37)))

    shell_fp.write(
        "#BSUB -f \"%s/download.err < /tmp/%sJ/%sJ.err\"\n" % (prot_dest_dir,
                                                               chr(37),
                                                               chr(37)))

    # delete everything on termination or completion of job
    shell_fp.write("#BSUB -Ep \"rm -rf %s\"\n" % tmp_dir)

    shell_fp.write("#BSUB -g %s/%s \n\n" % (gc.LSF_GEN_GROUP % domain))

    # call executable
    shell_fp.write("python %s %s %s \n\n" %
                   (gc.GEN_DWLD_EXEC, os.path.join(prot_dest_dir, upid + ".json"),
                    prot_dest_dir))

    # copy files to destination
    shell_fp.write("cp %s/*.gz %s/.\n" % (tmp_dir, prot_dest_dir))


# -----------------------------------------------------------------------------

def load_upid_gca_file(upid_gca_file):
    """
    Parses Uniprot's upid tsv file and exports all important information
    in json format

    upid_gca_file: UPID_GCA file provided by trembl

    returns: A dictionary of upid, gca and domain mappings {upid: {"GCA": , "DOM": }}
    """

    upid_gca_dict = {}
    upid_fp = open(upid_gca_file, 'r')

    try:
        for upid_line in upid_fp:
            upid_line = upid_line.strip().split('\t')

            # add GCA accession
            if upid_line[1] != '':
                upid_gca_dict[upid_line[0]] = {"GCA": upid_line[1]}
            else:
                upid_gca_dict[upid_line[0]] = {"GCA": -1}

            upid_gca_dict[upid_line[0]]["DOM"] = upid_line[2]

        upid_fp.close()

    except:
        raise IOError

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
                download_fp = open(os.path.join(genome_dir_loc, "download.out"))

                for line in download_fp.readlines():
                    if line.find("Success") != -1:
                        success = 1
                        break

                download_fp.close()

                if success == 0:
                    erroneous_genomes[kingdom].append(genome)

            success = 0

    fp_out = open(os.path.join(genome_dir, "download_report.txt"), 'w')

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

    upid_gca_fp = open(upid_gca_file, 'r')
    acc_pairs = json.load(upid_gca_fp)
    upid_gca_fp.close()

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
        fp_out = open(os.path.join(project_dir, "err_seq_rep_files.json"), 'w')
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
    seq_rep_fp = open(seq_report_file, 'r')

    # discard header line
    seq_rep_fp.readline()

    for line in seq_rep_fp:
        line = line.strip().split('\t')

        acc_dict[line[0]] = {"sequence_name": line[1], "sequence-length": line[2],
                             "sequence-role": line[3], "replicon-name": line[4],
                             "replicon-type": line[5], "assembly-unit": line[6]}
    seq_rep_fp.close()

    if dest_dir is None:
        dest_dir = os.path.split(seq_report_file)[0]

    filename = os.path.basename(seq_report_file).partition(".")[0]
    fp_out = open(os.path.join(dest_dir, filename + ".json"), 'w')
    json.dump(acc_dict, fp_out)
    fp_out.close()

    return acc_dict

# -----------------------------------------------------------------------------


def split_and_download(wgs_range, dest_dir):
    """
    Function to split and download smaller segments of large genome assemblies

    wgs_range: A WGS assembly sequence accession range from ENA
    (e.g. CBTL0100000001-CBTL0111673940)
    dest_dir: The path to the destination directory

    returns:  void
    """

    # split the range into separate accessions
    accessions = fetch_wgs_range_accs(wgs_range)

    file_no = len(accessions) / MAX_SEQS
    remainder = len(accessions) % MAX_SEQS

    count = 0
    # indexes
    idx1 = 0
    idx2 = MAX_SEQS
    while count < file_no:
        accession = accessions[idx1] + '-' + accessions[idx2]
        urllib.urlretrieve(ENA_DATA_URL % accession,
                           os.path.join(dest_dir, accession + '.fa'))
        idx1 = idx2 + 1
        idx2 = idx2 + MAX_SEQS

        count += 1

    # check if sequences are not split evenly, and do and extra step for the
    # remaining seqs

    if remainder != 0:
        idx1 = idx1 = idx2 + 1
        # get the last accession
        idx2 = accessions[-1]
        accession = accessions[idx1] + '-' + accessions[idx2]
        urllib.urlretrieve(ENA_DATA_URL % accession,
                           os.path.join(dest_dir, accession + '.fa'))


# -----------------------------------------------------------------------------


def fetch_accessions_from_proteome_xml(proteome):
    """
    Parses Uniprot's proteome xml and extracts all available ENA accessions

    proteome: A valid Uniprot's proteome accession

    returns: A list of genome accessions
    """
    prot_accessions = []

    # namespace prefix # or register a namespace in the ET
    prefix = "{http://uniprot.org/uniprot}%s"

    response = requests.get(gc.PROTEOME_XML_URL % proteome)

    if response.status_code == 200:
        # convert from string to xml format
        prot_tree_root = ET.fromstring(response.content)

        # get proteome node
        proteome = prot_tree_root.find(prefix % "proteome")

        # get proteome's component nodes/genome accession nodes
        component_nodes = proteome.findall(prefix % "component")

        # loop over all component nodes and extract genome accessions
        for node in component_nodes:
            gen_acc_nodes = node.findall(prefix % "genome_accession")

            for gen_acc_node in gen_acc_nodes:
                prot_accessions.append(gen_acc_node.text)

    return prot_accessions


# -----------------------------------------------------------------------------


def check_accession_availability(accession):
    """
    Check whether a specific accession is available from ENA

    accession: sequence accession
    return: True if accession is available, False otherwise
    """

    # we can expand this by adding a db option (e.g. ena, uniprot, ncbi)
    response = requests.get(ENA_XML_URL % accession)

    if response.status_code == httplib.OK:
        xml_root = ET.fromstring(response.content)

        # If the entry exists, there should be an entry node in the xml file
        entry_node = None
        entry_node = xml_root.find("entry")

        if entry_node is None:
            return False

    return True

# -----------------------------------------------------------------------------


def copy_wgs_set_from_ftp(wgs_acc, dest_dir):
    """
    Copy wgs set sequences from physical location on cluster

    wsg_acc: A valid WGS set accession (e.g. AAVU01000000)
    dest_dir: Destination directory where the sequences will be copied to

    return: void
    """

    # build path
    wgs_subdir = os.path.join(gc.ENA_FTP_WGS_PUB, wgs_acc[0:2].lower()) #AA
    wgs_filename = wgs_acc[0:6] + ".fasta.gz"

    # check if wgs sequences are in public dir and copy to destination
    if os.path.exists(os.path.join(wgs_subdir, wgs_filename)):
        shutil.copyfile(os.path.join(wgs_subdir, wgs_filename),
                        os.path.join(dest_dir, wgs_filename))

    # look for wgs set in suppressed sequences
    else:
        wgs_subdir = os.path.join(gc.ENA_FTP_WGS_SUP, wgs_acc[0:2].lower())
        if os.path.exists(os.path.join(wgs_subdir, wgs_filename)):
            shutil.copyfile(os.path.join(wgs_subdir, wgs_filename),
                            os.path.join(dest_dir, wgs_filename))

        else:
            sys.exit("WGS set %s requested does not exist." % wgs_acc)

# -----------------------------------------------------------------------------


def proteome_xml_accessions_to_dict(upid):
    """
    Parses a valid proteome xml file and returns all accessions in the form of
    a dictionary. Component names from proteome xml are used as dictionary keys

    upid: A valid Uniprot proteome upid

    returns: A dictionary with all proteome associated accessions.
    """

    proteome_accs = {"GCA": -1, "WGS": -1}
    other = {}

    # namespace prefix # or register a namespace in the ET
    prefix = "{http://uniprot.org/uniprot}%s"

    response = requests.get(gc.PROTEOME_XML_URL % upid)

    if response.status_code == 200:
        # convert from string to xml format
        prot_tree_root = ET.fromstring(response.content)

        # get proteome node
        proteome = prot_tree_root.find(prefix % "proteome")

        # get proteome's component nodes/genome accession nodes
        gca_acc = None
        gca_acc = proteome.find(prefix % "genomeAssembly").find(prefix % "genomeAssembly")
        if gca_acc is not None:
            proteome_accs["GCA"] = gca_acc.text

        component_nodes = proteome.findall(prefix % "component")

        # loop over all component nodes and extract genome accessions
        for node in component_nodes:
            name = node.get("name")

            if name.find("WGS") != -1:
                accession = node.find(prefix % "genome_accession").text
                proteome_accs["WGS"] = accession

            else:
                accession = node.find(prefix % "genome_accession")
                # if there is an accession available
                if accession is not None:
                    accession = accession.text
                    other[name] = accession

        proteome_accs["OTHER"] = other

    return proteome_accs

# -----------------------------------------------------------------------------


def copy_gca_report_file_from_ftp(gca_accession, dest_dir):
    """
    Copies the corresponding GCA report file from the ftp

    gca_accession: A valid GCA accession

    return: True if the file was found, False otherwise
    """

    seq_report_file = gca_accession + "_sequence_report.txt"
    genomic_regions_file = gca_accession + "_regions.txt"

    # 1st layer subdir GCA_XXX
    gca_dir = os.path.join(gc.ENA_GCA_SEQ_REPORT, gca_accession[0:7])

    # 2nd layer subdir GCA_XXXXXX
    gca_dir = os.path.join(gca_dir, gca_accession[0:10])

    report_file_path = os.path.join(gca_dir, seq_report_file)
    region_file_path = os.path.join(gca_dir, genomic_regions_file)

    # sanity check if the file actually exists
    if os.path.exists(report_file_path):
        shutil.copyfile(report_file_path, os.path.join(dest_dir,
                                                       seq_report_file))
        if os.path.exists(region_file_path):
            shutil.copyfile(region_file_path, os.path.join(dest_dir,
                                                       genomic_regions_file))
        return True

    return False

# -----------------------------------------------------------------------------


def get_genome_unique_accessions(upid, to_file=False, output_dir=None):
    """
    This function will extract all available accessions from the relevant
    proteome xml file and return a list of unique accessions that represent a
    complete genome. This will be a combination of assembly accessions provided
    by ENA and any additional accessions found in the proteome xml file

    upid: A valid Uniprot Proteome id
    output_dir: The path to the output dir

    return: A list with all unique genome accessions
    """

    # GCA NA - Set to 1 when GCA accession is available, but GCA report file is not available from ENA
    complete_genome_accs = {"GCA": -1, "WGS": -1, "OTHER": [], "GCA_NA": 0}

    proteome_acc_dict = proteome_xml_accessions_to_dict(upid)

    complete_genome_accs["GCA"] = proteome_acc_dict["GCA"]
    complete_genome_accs["WGS"] = proteome_acc_dict["WGS"]

    if proteome_acc_dict["GCA"] != -1:
        # create a temporary copy of the assembly report file

        if output_dir is None:
            output_dir = "/tmp"

        check_exists = copy_gca_report_file_from_ftp(proteome_acc_dict["GCA"], output_dir)

        # try downloading the files from the URL if unsuccessful
        if check_exists is False:
            url_check = download_gca_report_file_from_url(proteome_acc_dict["GCA"], output_dir)

        # get assembly report file path
        gca_report_filename = proteome_acc_dict["GCA"] + "_sequence_report.txt"

        if os.path.exists(os.path.join(output_dir, gca_report_filename)):

            gca_accs = assembly_report_parser(os.path.join(output_dir, gca_report_filename),
                                              url=False)

            accs_no_version = [x.partition('.')[0] for x in gca_accs]

            proteome_set = set(proteome_acc_dict["OTHER"].values())
            gca_set = set(accs_no_version)

            # construct a new set with unique accessions from both sets
            unique_accs = proteome_set.union(gca_set)
                
            # add unique accessions in dictionary
            complete_genome_accs["OTHER"].extend(unique_accs)

        else:
            print "Genome Assembly report file for %s is unavailable" % upid
            complete_genome_accs["OTHER"].extend(proteome_acc_dict["OTHER"].values())

            if complete_genome_accs["WGS"] != -1:
                complete_genome_accs["GCA_NA"] = 1

    else:
        complete_genome_accs["OTHER"].extend(proteome_acc_dict["OTHER"].values())

    # write proteome accessions to json file
    if to_file is True:
        fp_out = open(os.path.join(output_dir, upid+"_accessions.json"), 'w')
        json.dump(proteome_acc_dict, fp_out)
        fp_out.close()


    return complete_genome_accs

# -----------------------------------------------------------------------------


def extract_wgs_acc_from_gca_xml(gca_accession):
    """
    Parses ENA's GCA xml file and extracts the WGS set accession if available

    gca_accession:  A valid GCA accession

    return: A WGS set accesison, None if not found
    """

    xml_root = None
    wgs_acc = None
    assembly_xml = requests.get(ENA_XML_URL % gca_accession).content

    if os.path.isfile(assembly_xml):
        # parse xml tree and return root node
        xml_root = ET.parse(assembly_xml).getroot()
    else:
        # fromstring returns the xml root directly
        xml_root = ET.fromstring(assembly_xml)

    assembly = xml_root.find("ASSEMBLY")

    if assembly is not None:
        # no assembly link provided - look for WGS element
        wgs_node = assembly.find("WGS_SET")
        # get wgs accession
        wgs_acc = get_wgs_set_accession(
            wgs_node.find("PREFIX").text, wgs_node.find("VERSION").text)

    return wgs_acc

# -----------------------------------------------------------------------------


def download_gca_report_file_from_url(gca_accession, dest_dir):
    """
    Loads an xml tree from a file or a string (usually an http response),
    and returns a list with the genome assembly's chromosomes

    accession: A valid ENA GCA accession (without the assembly version)
    """

    accessions = []
    root = None
    assembly_link = None
    assembly = None
    url_links = []

    assembly_xml = requests.get(ENA_XML_URL % gca_accession).content

    if os.path.isfile(assembly_xml):
        # parse xml tree and return root node
        root = ET.parse(assembly_xml).getroot()
    else:
        # fromstring returns the xml root directly
        root = ET.fromstring(assembly_xml)

    assembly = root.find("ASSEMBLY")

    if assembly is not None:
        # either parse the assembly report file or get the WGS range
        assembly_links = assembly.find("ASSEMBLY_LINKS")

        if assembly_links is not None:
            # export url link and fetch all relevant assembly accessions
            assembly_link_nodes = assembly_links.findall(
                "ASSEMBLY_LINK")

            for node in assembly_link_nodes:
                assembly_link = node.find("URL_LINK").find("URL").text
                url_links.append(assembly_link.replace("ftp:", "http:"))

            for url in url_links:
                filename = url.split('/')[-1]
                urllib.urlretrieve(url, os.path.join(dest_dir, filename))

            return True

    return False

# -----------------------------------------------------------------------------


def get_genome_subdirectory_ranges(genome_acc_list):
    """
    This function generates a list of subdir ranges that can be used to
    organize genome files into multiple subdirs in a way that the location
    of a specific fasta file is easily detectable for the last 3 digits of the
    accession (e.g. JJRO01080032, KK558359). It takes into account LSF cluster
    limitations

    genome_acc_list: A list with all accessions in a particular assembly

    return: A list of indexes that will be used as subdirectory names
    """

    max_index = 999 # largest 3 digit number

    subdir_ranges = []
    file_indexes = []
    # construct a list with all the last 3 digits from the assembly accessions
    for acc in genome_acc_list:
        file_indexes.append(acc[-3:])

    # sort the list to devise the ranges
    file_idx_sorted = sorted(file_indexes)

    no_files = len(file_idx_sorted)

    index = gc.MAX_ALLOWED_FILES

    while index < no_files:
        subdir_ranges.append(file_idx_sorted.pop(index))
        index = index + gc.MAX_ALLOWED_FILES # increase by max allowed files

    # append the right most index
    if file_idx_sorted[-1] < max_index:
        subdir_ranges.append(file_idx_sorted[-1])
    else:
        subdir_ranges.append(max_index)

    return subdir_ranges

# -----------------------------------------------------------------------------

if __name__ == '__main__':

    pass
