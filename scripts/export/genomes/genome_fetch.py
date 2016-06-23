#!/usr/bin/python
'''
Created on 19 Nov 2015

@author: ikalvari

TO DO:    - logging
          - http request error handling (DEBUG.log)
          - move URLs in a config file
          - need to split WGS range to single accessions
'''

# ---------------------------------IMPORTS-------------------------------------


import os
import string
import xml.etree.ElementTree as ET
import urllib2
import urllib
import requests
import httplib
from rdflib import Graph

# ----------------------------------GLOBALS------------------------------------

# URLs
# This returns the list of reference proteomes from Uniprot
REF_PROT_LIST_URL = "http://www.uniprot.org/proteomes/?query=*&fil=reference%3Ayes&format=list"

# Retrieve the proteome rdf file
PROTEOME_URL = "http://www.uniprot.org/proteomes/%s.rdf"

# Retrieve the genome's xml file
ENA_XML_URL = "http://www.ebi.ac.uk/ena/data/view/%s&display=xml"

# ENA url for file download
ENA_DATA_URL = "http://www.ebi.ac.uk/ena/data/view/%s&display=%s&download=gzip"

# ENA url for assembly data retrieval via taxon id
ENA_TAX_URL = "http://www.ebi.ac.uk/ena/data/warehouse/search?query=\"tax_eq(%s)\"&result=assembly&display=xml"

# ENA file formats
FORMATS = {"xml": ".xml", "fasta": ".fa"}

GCA_REP_LABEL = 'Sequence Report'

# ---------------------------------------------------------------------- #STEP1


def fetch_ref_proteomes():
    '''
        This method returns a list of all reference proteome accessions
        available through Uniprot
    '''

    ref_prot_list = []
    response = urllib2.urlopen(REF_PROT_LIST_URL)

    for ref_prot in response:
        ref_prot_list.append(ref_prot.strip())

    return ref_prot_list

# -----------------------------------------------------------------------------


def export_gca_accessions(upid_gca):
    '''
        Retrieves reference proteomes ids and their associated gca accessions
        as well as the taxonomic rank/domain (eukaryotes, bacteria etc)

        upid_gca: Uniprot's tab separated file (UPID_GCA.tsv )
    '''

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
    '''
        Extracts and returns the assembly accession from the proteome rdf
        which provided as input. Returns -1 if not available.

        prot_rdf: A Uniprot's proteome rdf url or file path
    '''
    g = Graph()

    response = requests.get(prot_rdf).status_code

    if response == httplib.OK:

        g.load(prot_rdf)

        for s, p, o in g:
            if(string.find(o, "GCA") != -1):
                return os.path.split(o)[1]

    return -1

# -----------------------------------------------------------------------------


def proteome_rdf_scanner(proteome):
    '''
        Scans a Uniprot's reference proteome rdf file and looks for all
        available accessions. Returns a dictionary with GCA and WGS accessions
        where applicable

        prot_rdf: A Uniprot's proteome rdf url or file path
    '''

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
            if(string.find(o, "/embl/") != -1):
                if (string.find(o, "GCA") != -1):
                    accessions["GCA"] = os.path.split(o)[1]

                elif (wgs_flag is True):
                    accessions["WGS"] = os.path.split(o)[1]

            # if WGS keyword found, set flag to true
            elif (string.find(o, "WGS") != -1):
                wgs_flag = True
    else:
        pass

    return accessions

# ---------------------------------------------------------------------- #STEP2


def fetch_genome_acc(prot):
    '''
        Returns a proteome's corresponding assembly accession (ENA) in a
        dictionary format {proteome_acc:gca_acc}


        prot: One of (file|list|acc)
              - file: Uniprot's proteome list file
              - list: a list of reference proteome accessions (fetch_ref_proteomes)
              - acc: A single Uniprot ref. proteome accession
    '''

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
        rdf_url = "http://www.uniprot.org/proteomes/%s.rdf" % (proteome)
        res_handle = urllib.urlopen(rdf_url)
        if res_handle.getcode() == httplib.OK:
            gen_acc = extract_genome_acc(res_handle)
        gens[proteome] = gen_acc

        return gens

    # do this for files and lists
    for proteome in ref_prot_list:
        proteome = string.strip(proteome)
        rdf_url = "http://www.uniprot.org/proteomes/%s.rdf" % (proteome)
        res_handle = urllib.urlopen(rdf_url)
        if res_handle.getcode() == httplib.OK:
            gen_acc = extract_genome_acc(res_handle)

        gens[proteome] = gen_acc

    return gens

# -----------------------------------------------------------------------------


def fetch_ena_file(acc, file_format, dest_dir):
    '''
        Retrieves a file given a valid ENA accession and stores it in the
        indicated destination in the selected format

        acc: A valid ENA entry accession
        format: A valid ENA file format
        dest_dit: A valid path to destination directory
    '''

    seq_url = None

    seq_url = ENA_DATA_URL % (acc, file_format)

    # fetching compressed file
    ena_file = open(
        os.path.join(dest_dir, acc + FORMATS[file_format] + '.gz'), 'w')

    file_content = requests.get(seq_url).content

    ena_file.write(file_content)

    ena_file.close()


# ---------------------------------------------------------------------- #STEP3


def extract_assembly_accs(accession):
    '''
        Loads an xml tree from a file or a string (usually an http response),
        and returns a list with the genome assembly's chromosomes

        accession: A valid ENA assembly accession (without the assembly version)
    '''

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

    assembly = root.find('ASSEMBLY')

    if assembly is not None:

        # either parse the assembly report file or get the WGS range
        assembly_link = assembly.find(
            'ASSEMBLY_LINKS')

        if assembly_link is not None:
            # export url link and fetch all relevant assembly accessions
            url_link = assembly_link.find(
                'ASSEMBLY_LINK').find('URL_LINK').find('URL').text

            # need to check for GCA_REP_LABEL
            accessions = assembly_report_parser(url_link)

        else:
            # no assembly link provided - look for WGS element
            wgs = None
            wgs = assembly.find('WGS_SET')

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
    '''
        Downloads all chromosome files of a given assembly accession (ENA) in
        dest_dir

        gen: Single or a list of genome accessions (GC*)
        dest_dir: The path of the destination directory to export the fasta
        files
    '''

    # need to add logging

    accessions = None

    if(os.path.isfile(gen)):
        fp = open(gen, 'r')

        for gen_acc in fp:
            gen_acc = string.strip(gen_acc)
            if(string.find(gen_acc, '.') != -1):
                gen_acc = gen_acc.partition('.')
                gen_acc = gen_acc[0]

            gen_dir = os.path.join(dest_dir, gen_acc)

            try:
                os.mkdir(gen_dir)
            except:
                print 'Unable to generate directory for accession: ', gen

            accessions = extract_assembly_accs(gen_acc)

            if len(accessions) > 0:
                for acc in accessions:
                    fetch_ena_file(acc, "fasta", gen_dir)

            gen_acc = None
            gen_dir = None
            accessions = None

    else:
        if(string.find(gen, '.') != -1):
            gen = gen.partition('.')
            gen = gen[0]

        gen_dir = os.path.join(dest_dir, gen)

        try:
            os.mkdir(gen_dir)
        except:
            print 'Unable to generate directory for accession: ', gen

        accessions = extract_assembly_accs(gen)

        if len(accessions) > 0:
            for acc in accessions:
                fetch_ena_file(acc, "fasta", gen_dir)

        # if no accessions found, write to log file
        else:
            return


# -----------------------------------------------------------------------------


def fetch_genome(gen, dest_dir):
    '''
        Downloads and parses xml file of the given genome accession (gen), and
        downloads all chromosome files in fasta format in destination directory 
        (dest_dir)
        The xml file is removed after completion

        gen: ENA assembly accession (GCA*)
        dest_dir: destination of the output directory

    '''

    gen_dir = os.path.join(dest_dir, gen)

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
    '''
        Parses rdf url and returns a list of ENA accessions.

        rdf_url: The url to a Uniprot's reference proteome rdf url
        sub_str: A sub string to look for in the rdf file
    '''

    accessions = []
    rdf_graph = Graph()
    rdf_url = PROTEOME_URL % ref_prot_acc

    response = requests.get(rdf_url).status_code

    if response == httplib.OK:

        rdf_graph.load(rdf_url)

        for s, p, o in rdf_graph:
            if(string.find(o, sub_str) != -1):
                accessions.append(os.path.split(o)[1])
    else:
        # return http status code
        # return response.status_code
        pass

    return accessions


# -----------------------------------------------------------------------------


def assembly_report_parser(report_url):
    '''
        Parses an assembly report file and returns a list of all available
        accessions (scaffolds, contigs etc).
        To be called within extract_assembly_accs

        report_url: A url provided within an ENA assembly xml file. This is the
                    text of URL tag under ASSEMBLY/ASSEMBLY_LINKS/ASSEMBLY_LINK.
                    By default this is an ftp request url. Converting to http to
                    fetch assembly accessions.

        NOTE: accessions are exported including version

    '''

    accessions = []

    # switch from ftp to http to fetch assembly report on the go
    link_parts = ['http']
    link_parts.extend(report_url.partition(':')[1:])

    http_link = ''.join(link_parts)

    # fetch assembly report file contents and store in a list, omitting header
    ass_rep_file = requests.get(http_link).content.split('\n')[1:]

    # if empty line, remove it
    if ass_rep_file[len(ass_rep_file) - 1] == '':
        ass_rep_file.pop(len(ass_rep_file) - 1)

    # parse list and export assembly accessions
    for line in ass_rep_file:
        line = line.strip().split('\t')
        accessions.append(line[0])

    return accessions

# -----------------------------------------------------------------------------


def get_wgs_set_accession(prefix, version):
    '''
        Generates ENA WGS accession using the WGS accession pattern
        (4 letter prefix) and 2-digit build version. The WGS accession is
        generated by appending prefix and version with a postfix of 6-zeros.
        For more information please visit ENA service-news: http://goo.gl/LnIyQ3

        prefix: A 4-char string representing the ENA WGS accession
        version: A 2-digit representing the WGS build version

    '''

    postfix = '000000'
    wgs_accession = None

    if int(version) % 10 != 0:
        wgs_accession = prefix + '0' + version + postfix
    else:
        wgs_accession = prefix + version + postfix

    return wgs_accession

# -----------------------------------------------------------------------------


def get_wgs_range(wgs_acc):
    '''
        Fetches the wgs xml file from ENA and exports the wgs range

        wgs_acc: A valid ENA wgs accession
    '''

    wgs_range = None
    wgs_xml_str = requests.get(ENA_XML_URL % wgs_acc).content
    wgs_xml_root = ET.fromstring(wgs_xml_str)

    if wgs_xml_root.find('entry') is not None:

        wgs_xrefs = wgs_xml_root.find('entry').findall('xref')

        for xref_el in wgs_xrefs:
            if xref_el.get('db') == "ENA-WGS":
                wgs_range = xref_el.get('id')

    return wgs_range

# -----------------------------------------------------------------------------

if __name__ == '__main__':
    pass
