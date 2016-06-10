#!/usr/bin/python
'''
Created on 19 Nov 2015

@author: ikalvari

TO DO:    - logging
          - assembly report parser
          - handle proteomes with missing GCAs
'''
# ---------------------------------IMPORTS-------------------------------------

import os
import string
import xml.etree.ElementTree as ET
import urllib2
import requests

from rdflib import Graph

# ----------------------------------GLOBALS------------------------------------

# URLs
# This returns the list of reference proteomes from Uniprot
REF_PROT_LIST_URL = "http://www.uniprot.org/proteomes/?query=*&fil=reference%3Ayes&format=list"

# Retrieve the proteome rdf file
PROTEOME_URL = "http://www.uniprot.org/proteomes/%s&display=rdf"

# Retrieve the genome's xml file
GCA_XML_URL = "http://www.ebi.ac.uk/ena/data/view/%s&display=xml"

# ENA url for file download
ENA_DATA_URL = "http://www.ebi.ac.uk/ena/data/view/%s&display=%s&download=gzip"

# ENA file formats
FORMATS = {"xml": ".xml", "fasta": ".fa"}

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


def extract_genome_acc(prot_rdf):
    '''
        Extracts and returns the assembly accession from the proteome rdf
        which provided as input. Returns -1 if not available

        prot_rdf: A Uniprot's proteome rdf url or file path
    '''
    g = Graph()
    g.load(prot_rdf)

    for s, p, o in g:

        if(string.find(o, "GCA") != -1):
            return os.path.split(o)[1]

    return -1

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
        gen_acc = extract_genome_acc(rdf_url)
        gens[proteome] = gen_acc

        return gens

    # do this for files and lists
    for proteome in ref_prot_list:
        proteome = string.strip(proteome)
        rdf_url = "http://www.uniprot.org/proteomes/%s.rdf" % (proteome)
        gen_acc = extract_genome_acc(rdf_url)

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


def extract_chrom_accs(gen_acc):
    '''
        Loads an xml tree from a file or a string (usually an http response),
        and returns a list with the genome assembly's chromosomes

        gen_acc: A valid ENA assembly accession (without the assembly version)
    '''

    chrom_accs = []
    root = None
    chroms = None

    assembly_xml = requests.get(GCA_XML_URL % gen_acc).content

    if os.path.isfile(assembly_xml):
        # parse xml tree and return root node
        root = ET.parse(assembly_xml).getroot()
    else:
        # fromstring returns the xml root directly
        root = ET.fromstring(assembly_xml)

    # find ASSEMBLY node within genome xml tree
    chroms = root.find('ASSEMBLY').find('CHROMOSOMES')

    # if genome assembly is at chromosomal level
    if chroms is not None:
        # fetch chromosome tree
        for chrom in chroms.findall('CHROMOSOME'):
            chrom_accs.append(chrom.get('accession').partition('.')[0])

    # lower level of assembly - need to look for the url and fetch the data
    # in the report file
    else:
        # do something - need to read from the xml url here
        pass

    return chrom_accs

# ---------------------------------------------------------------------- #STEP4


def download_genomes(gen, dest_dir):
    '''
        Downloads all chromosome files of a given assembly accession (ENA) in
        dest_dir

        gen: Single or a list of genome accessions (GC*)
        dest_dir: The path of the destination directory to export the fasta
        files
    '''

    chroms = None

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
                pass

            chroms = extract_chrom_accs(gen_acc)

            if len(chroms) > 0:
                for chrom in chroms:
                    fetch_ena_file(chrom, "fasta", gen_dir)

            gen_acc = None
            gen_dir = None
            chroms = None

    else:
        if(string.find(gen, '.') != -1):
            gen = gen.partition('.')
            gen = gen[0]

        gen_dir = os.path.join(dest_dir, gen)

        try:
            os.mkdir(gen_dir)
        except:
            pass

        chroms = extract_chrom_accs(gen)

        if len(chroms) > 0:
            for chrom in chroms:
                fetch_ena_file(chrom, "fasta", gen_dir)
        # do nothing if no chroms - should write a log file
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
    chroms = extract_chrom_accs(gen_xml)

    for chrom in chroms:
        fetch_ena_file(chrom, "fasta", gen_dir)

    os.remove(gen_xml)  # remove xml file when done

# -----------------------------------------------------------------------------


def assembly_report_parser(assembly_report_url):
    '''
        Parses an assembly report file and returns a list of available
        accessions. To be called within extract_chrom_accs
    '''

    # TODO

    pass


# -----------------------------------------------------------------------------

if __name__ == '__main__':

    pass
