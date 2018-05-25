import os
import sys
import shutil
import subprocess
import json

from config import rfam_config as rc

# --------------------------------------------------------------------------------------------

def map_rfam_ensembl_accession_mapping(rfam_acc_file, ensembl_acc_file, dest_dir):
    """

    :param rfam_acc_file:
    :param ensembl_acc_file:
    :param dest_dir:
    :return:
    """

    gbh_dict = {}
    ensembl_accs = {}
    rfam_accs = {}

    # initialization
    if not os.path.exists(dest_dir):
        os.mkdir(dest_dir)

    # load ensembl accessions
    fp_ensembl = open(ensembl_acc_file, 'r')

    for line in fp_ensembl:
        line = line.strip().split('\t')
        ensembl_accs[line[0].partition('.')[0]] = {"assembly_name": line[1],
                                                   "tax_id": line[2], "species_name": line[3]}

    fp_ensembl.close()

    # load rfam accessions
    fp = open(rfam_acc_file, 'r')

    for line in fp:
        line = line.strip().split('\t')
        rfam_accs[line[1].partition('.')[0]] = line[0]
    fp.close()

    for gca_acc in ensembl_accs.keys():
        if gca_acc in rfam_accs:
            gbh_dict[rfam_accs[gca_acc]] = {}
            gbh_dict[rfam_accs[gca_acc]]["assembly_name"] = ensembl_accs[gca_acc]["assembly_name"]
            gbh_dict[rfam_accs[gca_acc]]["species_name"] = ensembl_accs[gca_acc]["species_name"]

    fp_out = open(os.path.join(dest_dir, "rfam_gbh_ids.json"), 'w')
    json.dump(gbh_dict, fp_out)

# --------------------------------------------------------------------------------------------


def generate_genome_text_file_from_file(names_list, release_version, dest_dir):
    """
    Generates the genome.txt file for the genome_browser_hup given a list
    of assembly-scientific name mappings

    names_list: A tab delimited file containing mappings of the
    assembly name and the genome scientific name e.g hg38\thomo_sapiens
    version: The version of the release
    dest_dir: The destination directory where the genome.txt
    file will be generated

    return:
    """

    trackdb_url = "ftp://ftp.ebi.ac.uk/pub/databases/Rfam/%s/genome_browser_hub/%s/trackDb.txt"

    fp_in = open(names_list, 'r')

    name_mappings = {}

    # parse input file and store information in a dictionary
    for line in fp_in:
        line = line.strip().split('\t')

        if line[0] not in name_mappings:
            name_mappings[line[0]] = line[1]

    fp_in.close()

    fp_out = open(os.path.join(dest_dir, "genomes.txt"), 'w')

    for assembly_acc in name_mappings.keys():
        fp_out.write("genome %s\n" % assembly_acc)
        fp_out.write(trackdb_url + "\n\n" % (str(release_version),
                                             name_mappings[assembly_acc]))

    fp_out.close()

# --------------------------------------------------------------------------------------------


def generate_genome_text_file_from_dict(accession_dict, release_version, dest_dir):
    """
    Generates the genome.txt file for the genome_browser_hup given a list
    of assembly-scientific name mappings

    names_list: A tab delimited file containing mappings of the
    assembly name and the genome scientific name e.g hg38\thomo_sapiens
    version: The version of the release
    dest_dir: The destination directory where the genome.txt
    file will be generated

    return:
    """

    trackdb_url = "ftp://ftp.ebi.ac.uk/pub/databases/Rfam/%s/genome_browser_hub/%s/trackDb.txt"

    fp_out = open(os.path.join(dest_dir, "genomes.txt"), 'w')

    for genome in accession_dict.keys():
        fp_out.write("genome %s\n" % accession_dict[genome]["assembly_name"])
        fp_out.write(trackdb_url % (str(release_version),
                                             accession_dict[genome]["assembly_name"]) + "\n\n")

    fp_out.close()

# --------------------------------------------------------------------------------------------


def generate_hub_txt_file(release_version, dest_dir):
    """
    Generates the hub.txt file for the new Rfam release

    release_version: The Rfam release version
    dest_dir:

    return: void
    """

    # convert to string and chop off the decimals
    rel_version_int = str(release_version).partition(".")[0]
    fp_out = open(os.path.join(dest_dir, "hub.txt"), 'w')
    fp_out.write("hub rfam%s\n" % rel_version_int)
    fp_out.write("shortLabel rfam%s_ncRNA\n" % rel_version_int)
    fp_out.write("longLabel Rfam %s non-coding RNA annotation\n" % release_version)
    fp_out.write("genomesFile genomes.txt\n")
    fp_out.write("email %s\n" % rc.RFAM_EMAIL)
    fp_out.write("descriptionUrl " + rc.BROWSER_HUB_DESC_URL % str(release_version) + '\n')
    fp_out.write("hub_description.html\n")
    fp_out.close()

# --------------------------------------------------------------------------------------------


def generate_hub_description_html_file(release_version, dest_dir):
    """

    release_version:

    return: void
    """

    fp_out = open(os.path.join(dest_dir, "hub_description.html"), 'w')

    # description section
    fp_out.write("<P><B>Description</B></P>\n\n")
    fp_out.write("<P class=\"indent1\">\n")
    fp_out.write("These tracks show non-coding RNAs from the")
    fp_out.write("<a href=\"http://rfam.org/\">Rfam database</a>.")
    fp_out.write("These tracks were created from the Rfam %s" % str(release_version))
    fp_out.write("covariance models (CMs), using the cmsearch program from the"
                 "<a href=\"http://infernal.janelia.org/\">Infernal</a> package."
                 "Further information about the methods used to build Rfam families"
                 "can be found <a href=\"http://rfam.xfam.org/help\">on our web page</a>.\n"
                 "</P>\n\n")

    # write methods
    fp_out.write("<P><B>Methods</B></P>\n\n")
    fp_out.write("<P class=\"indent1\">\n")
    fp_out.write("The cmsearch program was used with the --rfam, --cut_ga and --nohmmonly"
                 "command line options to calculate the positions of non-coding RNAs on each"
                 "genome assembly. This reproduces the method used by Rfam curators when creating"
                 "each family, including the use of expertly curated thresholds for each family.\n")
    fp_out.write("</P>\n\n")

    # credits
    fp_out.write("<P><B>Credits</B></P>\n\n")
    fp_out.write("<P class=\"indent1\">\n")
    fp_out.write("This data was generated by the <a href=\"http://rfam.org/\">Rfam database</a>."
                 " Please <a href=\"mailto:rfam-help@ebi.ac.uk\">email</a>"
                 " us with any questions or suggestions.\n</P>\n\n")

    # references
    fp_out.write("<P><B>References</B></P>\n\n<P class=\"indent1\">\n")
    fp_out.write("<a href=\"https://academic.oup.com/nar/article/46/D1/D335/4588106\">Kalvari I, Argasinska J,"
                 "Quinones-Olvera N, Nawrocki EP, Rivas E, Eddy SR, Bateman A, Finn RD, Petrov AI."
                 "Rfam 13.0: shifting to a genome-centric resource for non-coding RNA families. Nucleic"
                 "Acids Res. 2017. PMID: 29112718.</a>\n")
    fp_out.write("</P>\n\n")

    fp_out.close()

# --------------------------------------------------------------------------------------------


def genome_browser_hub_id_list_parser(genome_id_list):
    """
    Parses the input file of genome browser hub and returns
    a dictionary with the ids to be used to generate all
    related sub directories and files

    genome_id_list: A tab delimited file containing the genome upids, assembly names and
    scientific names for each genome

    return: A dictionary in the form of {upid: {"assembly_name": "assembly name",
                                                "species": "scientific name"}
    """

    accession_mapings = {}

    fp_in = open(genome_id_list, 'r')

    if genome_id_list.endswith(".json"):
        accession_mapings = json.load(fp_in)

    else:
        # parse file and generate the dictionary
        for line in fp_in:
            line = line.strip().split('\t')

            if line[0] not in accession_mapings:
                accession_mapings[line[0]] = {}
                accession_mapings[line[0]]["assembly_name"] = line[1]
                accession_mapings[line[0]]["species_name"] = line[2]

    # close input file and return dictionary
    fp_in.close()
    return accession_mapings

# --------------------------------------------------------------------------------------------


def generate_trackdb_file(species, release_version, dest_dir):
    """
    Creates a new species trackDb.txt file for a given Rfam release version

    species: The name of the species directory
    release_version: The version of the Rfam release
    dest_dir: The path to the species directory where the trackDb file will
    be generated

    return: void
    """

    trackDb_fp = open(os.path.join(dest_dir, "trackDb.txt"), 'w')

    trackDb_fp.write("track Rfam\n")
    trackDb_fp.write("bigDataUrl ftp://ftp.ebi.ac.uk/pub/databases/Rfam/%s/genome_browser_hub/%s/bigbed\n" %
                     (release_version, species))
    trackDb_fp.write("shortLabel Rfam ncRNA\n"
                     "longLabel Rfam ncRNA annotations\n"
                     "type bigBed 4\n"
                     "url http://rfam.org/family/$$\n"
                     "visibility 3\n"
                     "color 102,0,0\n")
    trackDb_fp.write("html ftp://ftp.ebi.ac.uk/pub/databases/Rfam/%s/genome_browser_hub/hub_description.html\n" %
                     release_version)

    trackDb_fp.close()

# --------------------------------------------------------------------------------------------


def generate_new_genome_browser_hub_directories(genome_id_file, release_version, dest_dir, genome_project_dir):
    """
    Generates a new genome_browser_hub directory for an upcoming Rfam release

    genome_id_file: A tab delimited file containing the genome upids, assembly names and
    scientific names for each genome
    release_version: The Rfam version the genome browser hub derives from
    dest_dir: A valid path where to generate the directories. Preferably a new Rfam release
    directory

    return: void
    """

    accession_mapings = genome_browser_hub_id_list_parser(genome_id_file)

    # creating basic directory
    if not os.path.exists(dest_dir):
        os.mkdir(dest_dir)

    # create an output directory where to generate the bidbed files
    bed_files_dir = os.path.join(dest_dir, "bed_files")
    os.mkdir(bed_files_dir)

    # create a new genome_browser_hub directory
    new_gbh_dir = os.path.join(dest_dir, "genome_browser_hub")
    os.mkdir(new_gbh_dir)

    # generate genome.txt file
    generate_genome_text_file_from_dict(accession_mapings, release_version, new_gbh_dir)

    # generate hub.txt file
    generate_hub_txt_file(release_version, new_gbh_dir)

    # generate hub_description.html file
    generate_hub_description_html_file(release_version, new_gbh_dir)

    # change to bed_files dir
    os.chdir(bed_files_dir)
    # now generate genome directories, bigbed and trackDB files
    for genome in accession_mapings.keys():
        # replace spaces in species name with underscores
        # to be used as the genome directory name. If very long replace with common name
        genome_name = accession_mapings[genome]["species_name"].replace(' ', '_')
        genome_dir = os.path.join(new_gbh_dir, genome_name)
        tbl_file_path = os.path.join(os.path.join(os.path.join(genome_project_dir,
                                                               genome[-3:]), genome),
                                     genome+'.tbl')
        # call perl script to generate the bedfiles need upid.tbl file here
        cmd = "/path/to/tblout2bigBedGenomes.pl %s"

        subprocess.call(cmd % tbl_file_path, shell=True)
        shutil.copyfile(os.path.join(bed_files_dir, genome+'.bigBed'), os.path.join(genome_dir, "bigBed"))

        generate_trackdb_file(genome_name, release_version, genome_dir)

    print "Done!"

# --------------------------------------------------------------------------------------------


if __name__ == '__main__':

    genome_id_list = sys.argv[1]
    release_version = sys.argv[2]
    genome_project_dir = sys.argv[3]
    dest_dir = sys.argv[4]

    generate_new_genome_browser_hub_directories(genome_id_list, release_version, dest_dir, genome_project_dir)