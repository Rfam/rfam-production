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

"""
Created on 27 Jan 2016

@author: ikalvari

Description: A set of database functions to ease processing and data
             retrieval from rfam_live

 Todo  - modify reset_is_significant() to enable single clan reset
       - set_num_sig_seqs() we need a new field in the family table to hold the
         number of significant sequences per family. This can be performed after
         clan competition
       - merge functions set_num_sig_seqs and set_number_of_species to one and
         pass query and fields as parameters

       *Perhaps include this module within an RfamLive class as all of the
        functions are database specific.
        Convert the code to use python mysql transactions
        """

# ---------------------------------IMPORTS---------------------------------

import json
import os
import string
import sys

from scripts.export.genomes import fetch_gen_metadata as fgm
from utils import RfamDB
from config.rfam_config import RFAMREL, RFAMLIVE, PG, FB1

# -------------------------------------------------------------------------

RFAM_ACC = 0  # full region rfam_acc
SEQ_ACC = 1  # full region rfamseq_acc
START = 2  # full region seq_start
END = 3  # full region seq_end
EVAL = 4  # full region evalue
version = '14.0'


def fetch_all_rfam_accs():
    """
    Fetch all available Rfam accessions
    """
    cnx = RfamDB.connect()
    cursor = cnx.cursor(buffered=True)
    query = "SELECT rfam_acc FROM family ORDER BY rfam_acc ASC"
    cursor.execute(query)
    rfam_accs = [str(x[0]) for x in cursor.fetchall()]
    cursor.close()
    RfamDB.disconnect(cnx)
    return rfam_accs


def set_is_significant_to_zero(rfam_acc, rfamseq_acc):
    """
    Fetch the correct db entry from full_region table according to
    rfam_acc and rfamseq_acc and set is_significant field to zero (0)

    rfam_acc: RNA family accession
    rfamseq_acc: Family specific sequence accession
    """

    # maybe have this working out of the list which will be returned from

    # connect to db
    cnx = RfamDB.connect()

    # get a new buffered cursor
    cursor = cnx.cursor(buffered=True)

    # update is_significant field to 0
    query = ("UPDATE full_region SET is_significant=0 "
             "WHERE rfam_acc=\'%s\' AND rfamseq_acc=\'%s\'") % (rfam_acc, rfamseq_acc)

    cursor.execute(query)

    cnx.commit()

    cursor.close()
    RfamDB.disconnect(cnx)


# -------------------------------------------------------------------------


def set_is_significant_to_zero_adv(rfam_acc, rfamseq_acc, region):
    """
    Fetch the correct db entry from full_region table according to
    rfam_acc and rfamseq_acc and set is_significant field to zero (0)

    rfam_acc: RNA family accession
    rfamseq_acc: Family specific sequence accession
    """

    # maybe have this working out of the list which will be returned from

    # connect to db
    cnx = RfamDB.connect()

    # get a new buffered cursor
    cursor = cnx.cursor(buffered=True)

    # update is_significant field to 0
    query = ("UPDATE full_region SET is_significant=0 "
             "WHERE rfam_acc=\'%s\' AND rfamseq_acc=\'%s\' AND seq_start=%d") % (rfam_acc, rfamseq_acc, region)

    cursor.execute(query)

    cnx.commit()

    cursor.close()
    RfamDB.disconnect(cnx)


# -------------------------------------------------------------------------


def load_clan_seqs_from_db(clan_acc):  # tested
    """
    Loads specific clan family sequences from full_region table and returns
    a dictionary structure as {Rfam_acc:{Rfseq_acc:[start, end, evalue]}}
    for clan competition.

    This has been modified to accommodate sequence duplicates

    clan_acc: Clan accession as in Rfam
    """

    fam_seqs = {}

    # connect to db
    cnx = RfamDB.connect()

    # get a new buffered cursor
    cursor = cnx.cursor(raw=True)

    # Fetch clan specific family full_region data
    query = ("SELECT full_region.rfam_acc, full_region.rfamseq_acc, \
      full_region.seq_start, full_region.seq_end, full_region.evalue_score\n"
             "FROM full_region\n"
             "JOIN (SELECT rfam_acc FROM clan_membership WHERE clan_acc=\'%s\') as CLAN_FAMS\n"
             "ON CLAN_FAMS.rfam_acc=full_region.rfam_acc") % (clan_acc)

    # execute the query
    cursor.execute(query)

    # build family dictionary of sequences
    for row in cursor:

        if str(row[RFAM_ACC]) in fam_seqs.keys():

            if str(row[SEQ_ACC]) in fam_seqs[str(row[RFAM_ACC])].keys():

                fam_seqs[str(row[RFAM_ACC])][str(row[SEQ_ACC])].append(
                    (int(row[START]), int(row[END]), float(row[EVAL])))
            else:
                fam_seqs[str(row[RFAM_ACC])][str(row[SEQ_ACC])] = [(int(row[START]),
                                                                    int(row[END]), float(row[EVAL]))]
        else:
            fam_seqs[str(row[RFAM_ACC])] = {
                str(row[SEQ_ACC]): [(int(row[START]), int(row[END]), float(row[EVAL]))]}

    # close cursor and DB connection
    cursor.close()
    RfamDB.disconnect(cnx)

    return fam_seqs


# -------------------------------------------------------------------------


def load_clan_members_from_db(clan_acc):
    """
    Retrieves all clan family members from DB and returns a list of the family
    accessions.

    clan_acc: Clan accession as in Rfam
    """

    clan_members = []

    cnx = RfamDB.connect()
    cursor = cnx.cursor(raw=True)

    query = ("SELECT rfam_acc FROM clan_membership "
             "WHERE clan_acc=\'%s\'") % (clan_acc)

    cursor.execute(query)

    rows = cursor.fetchall()

    cursor.close()
    RfamDB.disconnect(cnx)

    for fam in rows:
        clan_members.append(str(fam[0]))

        return clan_members


# -------------------------------------------------------------------------


def load_clans_from_db():
    """
    Retrieves all clan family members from DB and returns a dictionary in the
    in the form of {clan_id : [FAM1, FAM2, ... ], ... }

    clan_acc: Clan accession as in Rfam
    """

    clans = {}

    cnx = RfamDB.connect()
    cursor = cnx.cursor(raw=True)

    query = "SELECT * FROM clan_membership"

    # execute query
    cursor.execute(query)

    # fetch the data
    rows = cursor.fetchall()

    cursor.close()
    RfamDB.disconnect(cnx)

    # create the dictionary
    for row in rows:
        if str(row[0]) not in clans.keys():
            clans[str(row[0])] = [str(row[1])]
    else:
        clans[str(row[0])].append(str(row[1]))

        return clans


# -------------------------------------------------------------------------


def set_is_singificant_to_zero_multi(non_sig_seqs):
    """
    A function for batching the process of updating full_region tables upon
    clan competition. Updates the full_region table setting is_significant
    field to zero (0) for the list of non significant sequences passed in
    the form of (rfam_acc, rfamseq_acc, seq_start) tuples.

    non_sig_seqs: A list of the non significant regions to be set to zero.
                  The list is product of clan competition.

                  """

    # connect to db
    cnx = RfamDB.connect()

    # get a new buffered cursor
    cursor = cnx.cursor(raw=True)

    # query to update is_significant field to 0
    query = ("UPDATE full_region SET is_significant=0 "
             "WHERE rfam_acc=%s AND rfamseq_acc=%s AND seq_start=%s")

    try:
        # execute query batched
        cursor.executemany(query, non_sig_seqs)
        cnx.commit()

    except:
        print("MySQL Update Error. Rolling back...")
        cnx.rollback()
        cursor.close()
        RfamDB.disconnect(cnx)

        cursor.close()
        RfamDB.disconnect(cnx)


# -------------------------------------------------------------------------


def reset_is_significant(clan_comp_type='FULL'):
    """
    This function resets full_region's is_singificant field's back to 1.
    This should be able to update all or part of the table for clan
    competition initialization and restoration.
    """
    seq_regs = []

    cnx = RfamDB.connect()

    # cursor to fetch data
    d_cursor = cnx.cursor(buffered=True)

    # query to fetch all non significant sequences
    if clan_comp_type.upper() == 'FULL':
        select_query = ("SELECT rfam_acc, rfamseq_acc, seq_start FROM full_region "
                        "WHERE is_significant=0")

        # query to update 0 fields from s_query
        update_query = ("UPDATE full_region SET is_significant=1 "
                        "WHERE rfam_acc=%s AND rfamseq_acc=%s AND seq_start=%s")

    elif clan_comp_type.upper() == 'PDB':
        select_query = ("SELECT rfam_acc, pdb_id, chain, pdb_start from pdb_full_region "
                        "WHERE is_significant=0")

        update_query = ("UPDATE pdb_full_region SET is_significant=1 "
                        "WHERE rfam_acc=%s AND pdb_id=%s AND chain=%s AND pdb_start=%s")

        d_cursor.execute(select_query)

    # construct region list here
    for row in d_cursor:
        if clan_comp_type.upper() == 'FULL':
            seq_regs.append((str(row[0]), str(row[1]), int(row[2])))

        elif clan_comp_type.upper() == 'PDB':
            seq_regs.append((str(row[0]), str(row[1]), str(row[2]), int(row[3])))

    d_cursor.close()

    # get a new cursor for db updates
    u_cursor = cnx.cursor(raw=True)

    # update db
    try:
        u_cursor.executemany(update_query, seq_regs)
        cnx.commit()
    except:
        print("MySQL Update Error. Rolling back...")
        cnx.rollback()
        u_cursor.close()
        RfamDB.disconnect(cnx)

    u_cursor.close()
    RfamDB.disconnect(cnx)


# -------------------------------------------------------------------------


def update_post_process(jobs_file):
    """
    Updates _post_process table with the job_ids per family assigned by lsf

    jobs_file: This is a tab separated txt file generated from running the
    job_dequeuer.py script that submits the rfam_view_process for each
    family.
    (rfam_acc uuid job_id ...)
    """

    job_ids = []

    jobs_file_fp = open(jobs_file, 'r')

    query = ("UPDATE _post_process SET lsf_id=%s "
             "WHERE rfam_acc=%s AND uuid=%s")

    # get lsf ids from file
    for line in jobs_file_fp:
        line = line.strip()
        line = string.split(line, '\t')
        job_ids.append((line[2], line[0], line[1]))

        jobs_file_fp.close()

    # connect to db
    cnx = RfamDB.connect()
    cursor = cnx.cursor(raw=True)

    # update db
    try:
        cursor.executemany(query, job_ids)
        cnx.commit()  # move this after except statement??

    except:
        # rollback to previous state
        print("MySQL Update Error. Rollback...")
        cnx.rollback()
        cursor.close()
        RfamDB.disconnect(cnx)

    cursor.close()
    RfamDB.disconnect(cnx)


# -------------------------------------------------------------------------


def set_number_of_species():
    """
    Updates number_of_species in family table
    """

    cnx = RfamDB.connect()

    cursor = cnx.cursor(buffered=True)
    c_cursor = cnx.cursor(buffered=True)

    cursor.execute("Select rfam_acc from family")

    rfam_accs = cursor.fetchall()

    cursor.close()

    count_query = ("select count(distinct ncbi_id)\n"
                   "from full_region f, rfamseq r\n"
                   "where r.rfamseq_acc=f.rfamseq_acc\n"
                   "and is_significant=1 and rfam_acc=\'%s\'")

    count_query = """
        SELECT count(*) AS combined_count FROM (
            SELECT DISTINCT ncbi_id
            FROM full_region f, rfamseq r
            WHERE r.rfamseq_acc = f.rfamseq_acc
            AND is_significant=1 AND rfam_acc='%s'
        UNION
            SELECT DISTINCT ncbi_id
            FROM seed_region f, rfamseq r
            WHERE r.rfamseq_acc = f.rfamseq_acc
            AND rfam_acc='%s'
        ) table1"""

    # counts list
    counts = []
    for acc in rfam_accs:
        rfam_acc = str(acc[0])
        c_cursor.execute(count_query % (rfam_acc, rfam_acc))
        count = c_cursor.fetchall()
        counts.append((count[0][0], rfam_acc))
        count = 0

    c_cursor.close()
    c_cursor = cnx.cursor(buffered=True)

    # query to update number_of_species in the family table
    update_query = ("update family set number_of_species=%s where rfam_acc=%s")

    try:
        c_cursor.executemany(update_query, counts)
        cnx.commit()

    except:
        cnx.rollback()

    c_cursor.close()
    RfamDB.disconnect(cnx)

    print("Done updating number_of_species column in the family table")


# ----------------------------------------------------------------------------


def set_num_full_sig_seqs():
    """
    Updates num_full in family table to hold the number of significant
    sequences rather than the number of sequences in the full alignment
    """

    cnx = RfamDB.connect()

    cursor = cnx.cursor(buffered=True)
    c_cursor = cnx.cursor(buffered=True)

    cursor.execute("Select rfam_acc from family")

    rfam_accs = cursor.fetchall()

    cursor.close()

    # query to count all significant sequences of a family

    count_query = ("select count(*)\n"
                   "from full_region f\n"
                   "where is_significant=1\n"
                   "and type=\'full\'\n"
                   "and rfam_acc=\'%s\'")

    # counts list
    counts = []
    for acc in rfam_accs:
        c_cursor.execute(count_query % str(acc[0]))
        count = c_cursor.fetchall()[0][0]

        counts.append((count, str(acc[0])))

        count = 0

    c_cursor.close()
    c_cursor = cnx.cursor(buffered=True)

    update_query = ("update family set num_full=%s where rfam_acc=%s")

    try:
        c_cursor.executemany(update_query, counts)
        cnx.commit()

    except:
        cnx.rollback()

    c_cursor.close()
    RfamDB.disconnect(cnx)

    print("Done updating num_full column in the family table")


# ----------------------------------------------------------------------------


def update_family_ncbi():
    """
    Updates table family_ncbi by adding all distinct taxonomic ids per family.

    :return: void
    """
    cnx = RfamDB.connect()

    cursor = cnx.cursor(buffered=True)
    c_cursor = cnx.cursor(buffered=True)

    cursor.execute("select rfam_acc from family")

    rfam_accs = cursor.fetchall()

    cursor.close()

    # family_ncbi query
    get_ncbi_ids = ("select distinct rs.ncbi_id, f.rfam_id, "
                    "f.rfam_acc from full_region fr, rfamseq rs, family f "
                    "where fr.rfamseq_acc=rs.rfamseq_acc "
                    "and f.rfam_acc=fr.rfam_acc "
                    "and fr.rfam_acc=\'%s\' "
                    "and fr.is_significant=1")

    insert_query = "insert into family_ncbi (ncbi_id, rfam_id, rfam_acc) values (%s,%s,%s)"

    family_ncbi_entries = []
    cursor = cnx.cursor(buffered=True)
    cursor.execute("truncate table family_ncbi")
    for rfam_acc in rfam_accs:
        c_cursor.execute(get_ncbi_ids % rfam_acc[0])
        family_ncbi_entries = list(c_cursor.fetchall())
        entries_reformatted = [(str(x[0]), str(x[1]), str(x[2])) for x in family_ncbi_entries]

        try:
            cursor.executemany(insert_query, entries_reformatted)
            cnx.commit()

        except:
            cnx.rollback()
            sys.exit("\nError updating family_ncbi table for family %s." % rfam_acc[0])

        family_ncbi_entries = []
        entries_reformatted = []

    cursor.close()
    c_cursor.close()
    RfamDB.disconnect(cnx)

    print("Done updating family_ncbi table")


# ----------------------------------------------------------------------------


def fetch_clanin_data():
    """
    Fetches all rfam_ids per clan. To be used for clanin file generation

    :return: void
    """

    clan_members = {}
    cnx = RfamDB.connect()

    cursor = cnx.cursor(buffered=True)

    cursor.execute("select cm.clan_acc, f.rfam_id from clan_membership cm, family f "
                   "where f.rfam_acc=cm.rfam_acc "
                   "order by cm.clan_acc")

    clan_pairs = cursor.fetchall()

    cursor.close()

    # build clan membership dictionary
    for clan_pair in clan_pairs:
        clan_acc = clan_pair[0]
        rfam_id = clan_pair[1]

        if clan_acc not in clan_members.keys():
            clan_members[clan_acc] = [rfam_id]
        else:
            clan_members[clan_acc].append(rfam_id)

    cursor.close()
    RfamDB.disconnect(cnx)

    return clan_members


# ----------------------------------------------------------------------------


def set_pdb_is_significant_to_zero(non_sig_seqs, sync=False):
    """
    Sets pdb_full_region is_significant to 0 for non significant regions in non_sig_seqs
    :type non_sig_seqs: list
    :param non_sig_seqs: A list of the non significant regions to be set to zero.
    :param sync: Flag to sync release and web databases with the live db by updating the is_significant column
    The list is product of clan competition.
    """

    # reformat list by splitting pdb_id and chain
    pdb_reformatted_regions = []

    for competed_region in non_sig_seqs:
        # split pdb_id chain pairs by '_' used in concatenation for clan competition
        # pdb_id: pdb_id_chain_pairs[0] and chain: pdb_id_chain_pairs[2]
        pdb_id_chain_pairs = competed_region[1].partition('_')
        pdb_reformatted_regions.append((str(competed_region[0]), str(pdb_id_chain_pairs[0]),
                                        str(pdb_id_chain_pairs[2]), int(competed_region[2])))

    # connect to db
    if sync:
        config_list = [RFAMREL, PG, FB1]
        for db_conf in config_list:
            _update_is_significant(db_conf, pdb_reformatted_regions)

    else:
        _update_is_significant(RFAMLIVE, pdb_reformatted_regions)


def _update_is_significant(db_conf, pdb_reformatted_regions):
    """
    Update the is_significant column in the pdb_full_region_table of the given database, as specified by the db_conf
    :param db_conf: database config for the db to use
    :param pdb_reformatted_regions: formatted values to update in the table
    """
    cnx = RfamDB.connect(db_config=db_conf)

    # get a new buffered cursor
    cursor = cnx.cursor(raw=True)

    # query to update is_significant field to 0
    query = ("update pdb_full_region set is_significant=0 "
             "where rfam_acc=%s and pdb_id=%s and chain=%s and pdb_start=%s")

    try:
        # execute query batched
        cursor.executemany(query, pdb_reformatted_regions)
        cnx.commit()

    except Exception as e:
        print("MySQL Update Error. Rolling back...", e)
        cnx.rollback()

    cursor.close()
    RfamDB.disconnect(cnx)

    cursor.close()
    RfamDB.disconnect(cnx)


def fetch_clan_accessions():
    """
    Fetches all clan accessions from the database and returns then in the
    form of a list

    returns: A list of all clan accessions
    """
    cnx = RfamDB.connect()
    clan_cursor = cnx.cursor(buffered=True)

    clan_query = "SELECT clan_acc FROM clan"

    # fetch clans
    clan_cursor.execute(clan_query)
    clans = [str(x[0]) for x in clan_cursor.fetchall()]

    clan_cursor.close()
    RfamDB.disconnect(cnx)

    return clans


# ----------------------------------------------------------------------------


def fetch_clan_full_region_records(clan_acc):
    """
    Fetches all regions per clan

    param clan_acc: A valid Rfam clan accession

    returns: A list with all regions from full_region table for a specific  clan
    """

    cnx = RfamDB.connect()
    clan_cursor = cnx.cursor(buffered=True)

    clan_region_query = ("SELECT * FROM full_region\n"
                         "JOIN (SELECT rfam_acc FROM clan_membership WHERE clan_acc=\'%s\') as CLAN_FAMS\n"
                         "ON CLAN_FAMS.rfam_acc=full_region.rfam_acc")  # % (clan_acc)

    clan_cursor.execute(clan_region_query % clan_acc)

    clan_sequence_regions = clan_cursor.fetchall()

    clan_cursor.close()
    RfamDB.disconnect(cnx)

    return clan_sequence_regions


# ----------------------------------------------------------------------------


def fetch_clan_pdb_full_region_records(clan_acc):
    """
    Fetches all regions per clan

    param clan_acc: A valid Rfam clan accession

    returns: A list with all pdb regions per clan
    """

    cnx = RfamDB.connect()
    clan_cursor = cnx.cursor(buffered=True)

    clan_pdb_region_query = ("select pfr.rfam_acc, concat(pfr.pdb_id,'_',pfr.chain) as seq_acc, "
                             "pfr.pdb_start, pfr.pdb_end, pfr.bit_score, pfr.evalue_score "
                             "from pdb_full_region pfr, clan_membership cm "
                             "where cm.rfam_acc=pfr.rfam_acc "
                             "and cm.clan_acc=\'%s\' "
                             "order by seq_acc")

    clan_cursor.execute(clan_pdb_region_query % clan_acc)

    clan_sequence_regions = clan_cursor.fetchall()

    clan_cursor.close()
    RfamDB.disconnect(cnx)

    return clan_sequence_regions


# ----------------------------------------------------------------------------


def fetch_rfam_accs_sorted(order='DESC'):
    """
    Fetch all available Rfam accs and sort by specified order. DESC by default

    order: The order in which to sort the records (ASC, DESC)
    returns: void
    """

    # connect to db
    cnx = RfamDB.connect()

    # get a new buffered cursor
    cursor = cnx.cursor(buffered=True)

    # update is_significant field to 0
    query = ("select rfam_acc from seed_region\n"
             "group by rfam_acc\n"
             "order by count(*) %s" % order)

    cursor.execute(query)

    rfam_accs = [str(x[0]) for x in cursor.fetchall()]

    cursor.close()
    RfamDB.disconnect(cnx)

    return rfam_accs


# ----------------------------------------------------------------------------


def fetch_all_upids():
    """
    Fetch all available genome accessions from genome table

    return: A list of UP/RG ids as stored in genome
    """

    # connect to db
    cnx = RfamDB.connect()

    # get a new buffered cursor
    cursor = cnx.cursor(buffered=True)

    # update is_significant field to 0
    query = "select upid from genome"

    cursor.execute(query)

    genome_accs = [str(x[0]) for x in cursor.fetchall()]

    cursor.close()
    RfamDB.disconnect(cnx)

    return genome_accs


# ----------------------------------------------------------------------------


def set_genome_size(genome_sizes):
    """
    Updates total_length in genome table

    genome_sizes: This can be a json file for multiple updates or a tuple in
    the form of (size, upid) for single genome, where size is in nucleotides

    return: A list of UP/RG ids as stored in genome
    """

    # connect to db
    cnx = RfamDB.connect()

    # get a new buffered cursor
    cursor = cnx.cursor(buffered=True)

    # update is_significant field to 0
    query = "update genome set total_length=%s where upid=%s"

    genome_size_list = []
    if os.path.isfile(genome_sizes):
        gen_size_file = open(genome_sizes, 'r')
        genome_size_dict = json.load(gen_size_file)
        gen_size_file.close()
        genome_size_list = [(str(genome_size_dict[upid]), str(upid)) for upid in genome_size_dict.keys()]

    else:
        genome_size_list.append(genome_sizes)

    cursor.executemany(query, genome_size_list)
    cnx.commit()

    cursor.close()
    RfamDB.disconnect(cnx)


# ----------------------------------------------------------------------------


def set_number_of_distinct_families_in_genome(upid=None):
    """
    Sets the number distinct families with hits in a specific genome defined
    by its corresponding upid

    upid: A specific genome upid to update the number of distinct families

    return: void
    """

    upids = []
    # connect to db
    cnx = RfamDB.connect()

    # get a new buffered cursor
    cursor = cnx.cursor(buffered=True)

    if upid is None:
        upids = fetch_all_upids()

        for upid in upids:
            select_query = ("select count(distinct rfam_acc) from full_region fr, genseq gs\n"
                            "where fr.rfamseq_acc=gs.rfamseq_acc\n"
                            "and gs.upid=\'%s\'\n"
                            "and gs.version=\'%s\'")

            cursor.execute(select_query % (upid, version))
            count = cursor.fetchone()[0]

            # update is_significant field to 0
            update_query = "update genome set num_families=%d where upid=\'%s\'"

            # execute query
            cursor.execute(update_query % (count, upid))


    else:
        select_query = ("select count(distinct rfam_acc) from full_region fr, genseq gs\n"
                        "where fr.rfamseq_acc=gs.rfamseq_acc\n"
                        "and gs.upid=\'%s\'\n"
                        "and gs.version=\'%s\'")

        cursor.execute(select_query % (upid, version))
        count = cursor.fetchone()[0]

        # update is_significant field to 0
        update_query = "update genome set num_families=%d where upid=\'%s\'"

        # execute query
        cursor.execute(update_query % (count, upid))

    # commit changes and disconnect
    cnx.commit()
    cursor.close()
    RfamDB.disconnect(cnx)
    print("Done updating num_families column in the genome table")


# ----------------------------------------------------------------------------


def set_number_of_genomic_significant_hits(upid=None):
    """
    Sets the number of significant hits for a specific genome according to
    its corresponding upid id

    upid: A specific genome upid to update the number of significant hits

    return: void
    """

    # connect to db
    cnx = RfamDB.connect()

    # get a new buffered cursor
    cursor = cnx.cursor(buffered=True)

    if upid is None:

        upids = fetch_all_upids()

        for upid in upids:
            count_query = ("select count(fr.rfamseq_acc)\n"
                           "from full_region fr, genseq gs\n"
                           "where fr.rfamseq_acc=gs.rfamseq_acc\n"
                           "and fr.is_significant=1\n"
                           "and gs.upid=\'%s\'\n"
                           "and gs.version=\'%s\'")

            cursor.execute(count_query % (upid, version))
            count = cursor.fetchone()[0]

            # update is_significant field to 0
            update_query = "update genome set num_rfam_regions=%d where upid=\'%s\'"

            # execute query
            cursor.execute(update_query % (count, upid))
    else:

        count_query = ("select count(fr.rfamseq_acc)\n"
                       "from full_region fr, genseq gs\n"
                       "where fr.rfamseq_acc=gs.rfamseq_acc\n"
                       "and fr.is_significant=1\n"
                       "and gs.upid=\'%s\'\n"
                       "and gs.version=\'%s\'")

        cursor.execute(count_query % (upid, version))
        count = cursor.fetchone()[0]

        # update is_significant field to 0
        update_query = "update genome set num_rfam_regions=%d where upid=\'%s\'"

        # execute query
        cursor.execute(update_query % (count, upid))

    # commit changes and disconnect
    cnx.commit()
    cursor.close()
    RfamDB.disconnect(cnx)
    print("Done updating num_rfam_regions column in the genome table")


# ----------------------------------------------------------------------------


def fetch_author_orcid(author_name):
    """
    Searches for author by name and
    :param author_name:
    :return:
    """

    orcid = None
    cnx = RfamDB.connect()

    # Get a new buffered cursor
    cursor = cnx.cursor(buffered=True)

    query = """
    Select orcid from author
    where name like '%s%s%s' or synonyms like '%s%s%s'
    """

    cursor.execute(query % (chr(37), author_name, chr(37),
                            chr(37), author_name, chr(37)))

    result = cursor.fetchone()
    if result is not None:
        orcid = result[0]

        cursor.close()
        RfamDB.disconnect(cnx)

    # This will return none if there's no ORCiD available
    return orcid


# ----------------------------------------------------------------------------


def update_chromosome_info_in_genseq():
    # connect to db
    cnx = RfamDB.connect()

    # get a new buffered cursor
    cursor = cnx.cursor(buffered=True, dictionary=True)

    genome_query = "select upid, assembly_acc from genome where assembly_acc is not NULL"

    update_query = """
    update genseq set chromosome_type=\'%s\', chromosome_name=\'%s\'
    where upid=\'%s\' and rfamseq_acc=\'%s\' and version=14.0
    """

    cursor.execute(genome_query)
    accessions = cursor.fetchall()
    cursor.close()

    upid_gca_dict = {}

    cursor = cnx.cursor(buffered=True)

    for pair in accessions:
        upid_gca_dict[pair["upid"]] = pair["assembly_acc"]

        for upid in upid_gca_dict.keys():

            upid_gca_dict[upid]

            if upid_gca_dict[upid][0:3] == 'GCF' or upid_gca_dict[upid] == '':
                continue

            data = fgm.fetch_gca_data(upid, upid_gca_dict[upid], 'kingdom')

            if "fields" in data:
                fields = data["fields"]
                if "chromosomes" in fields:
                    for chromosome in fields["chromosomes"]:
                        cursor.execute(update_query % (str(chromosome["type"]), str(chromosome["name"]),
                                                       str(upid), str(chromosome["accession"])))

    cnx.commit()
    cursor.close()
    RfamDB.disconnect(cnx)


# ----------------------------------------------------------------------------


def update_assembly_names(upid_gca_file):
    """
    Loads the upid_gca json files and parses the corresponding assembly xml files
    from ENA to fetch the assembly names and update the fields in genome table

    param upid_gca_file: A json file with upid: {"GCA" : GCAxxx, "DOM": domain }

    return: void
    """

    fp = open(upid_gca_file, 'r')
    acc_pairs = json.load(fp)
    fp.close()

    # a list of tuples to
    assembly_names = []

    for upid in acc_pairs.keys():
        data = fgm.fetch_gca_data(upid, acc_pairs[upid]["GCA"], acc_pairs[upid]["DOM"])

        if "fields" in data:
            if data["fields"]["assembly_name"] is not None:
                assembly_names.append((data["fields"]["assembly_name"], upid))

    # connect to db
    cnx = RfamDB.connect()

    # get a new buffered cursor
    cursor = cnx.cursor(buffered=True, dictionary=True)

    query = "update genome set assembly_name=%s where upid=%s"

    cursor.executemany(query, assembly_names)
    cnx.commit()

    cursor.close()
    RfamDB.disconnect(cnx)


# ----------------------------------------------------------------------------


def get_number_of_seed_sequences(rfam_acc):
    """
    Gets the number of SEED sequences for a specific Rfam family from the
    database.

    rfam_acc: A valid Rfam family accession

    return (int): Number of SEED sequences
    """

    # connect to db
    cnx = RfamDB.connect()

    # get a new buffered cursor
    cursor = cnx.cursor(buffered=True)

    query = "Select count(*) from seed_region where rfam_acc=\'%s\'" % rfam_acc

    cursor.execute(query)

    number_seed_seqs = int(cursor.fetchone()[0])

    cursor.close()
    RfamDB.disconnect(cnx)

    return number_seed_seqs


# ----------------------------------------------------------------------------


def get_number_of_full_hits(rfam_acc):
    """
    Gets the number of FULL hits from the full_region table for a specific
    Rfam family.

    rfam_acc: A valid Rfam family accession

    return (int): Number of FULL hits from full_region table
    """

    # connect to db
    cnx = RfamDB.connect()

    # get a new buffered cursor
    cursor = cnx.cursor(buffered=True)

    query = "Select count(*) from full_region where rfam_acc=\'%s\' and type=\'full\' and is_significant=1" % rfam_acc

    cursor.execute(query)

    number_full_hits = int(cursor.fetchone()[0])

    cursor.close()
    RfamDB.disconnect(cnx)

    return number_full_hits


# ----------------------------------------------------------------------------


def fetch_metagenomic_regions():
    """
    Fetches all seed_region entries

    return: A list of tuples with all seed_region entries
    """

    # connect to db
    cnx = RfamDB.connect()

    # get a new buffered cursor
    cursor = cnx.cursor(buffered=True)

    # update is_significant field to 0
    query = ("Select rfam_acc, umgseq_acc, seq_start, seq_end "
             "from meta_full_region")

    cursor.execute(query)

    region_rows = cursor.fetchall()

    cursor.close()
    RfamDB.disconnect(cnx)

    return region_rows


# ----------------------------------------------------------------------------


def get_family_unique_ncbi_ids(rfam_acc):
    """
    Creates a list of unique NCBI ids per family based on the unique SEED and
    FULL NCBI ids

    rfam_acc: A valid Rfam family accession

    return (list): A list of unique NCBI ids associated with a specific Rfam
    family
    """

    seed_query = """
    select distinct rs.ncbi_id
    from seed_region sr, rfamseq rs
    where sr.rfamseq_acc=rs.rfamseq_acc
    and sr.rfam_acc = '%s'
    """

    full_query = """
    select distinct rs.ncbi_id
    from full_region fr, rfamseq rs
    where fr.rfamseq_acc=rs.rfamseq_acc
    and fr.rfam_acc = '%s'
    and fr.type = 'full'
    """

    cnx = RfamDB.connect()

    cursor_seed = cnx.cursor(buffered=True)
    cursor_full = cnx.cursor(buffered=True)

    cursor_seed.execute(seed_query % rfam_acc)

    seed_ncbi_ids = [x[0] for x in cursor_seed.fetchall()]

    cursor_full.execute(full_query % rfam_acc)
    full_ncbi_ids = [x[0] for x in cursor_full.fetchall()]

    unique_family_ncbi_ids = list(set(full_ncbi_ids).union(set(seed_ncbi_ids)))

    cursor_seed.close()
    cursor_full.close()
    RfamDB.disconnect(cnx)

    return unique_family_ncbi_ids


# ----------------------------------------------------------------------------


def fetch_type_specific_rfam_accessions(rna_type, return_type="list"):
    """
    Fetches all Rfam family accessions from the database matching the
    rna_type parameter

    rna_type: A string specifying a valid type of ncRNAs to extract
    from the database
    return_type: The python type the data will be return
    """

    query = """
    select rfam_acc from family
    where type like '%s%s%s'
    """

    cnx = RfamDB.connect()
    cursor = cnx.cursor(buffered=True)

    cursor.execute(query % (chr(37), rna_type, chr(37)))

    rfam_accs = {}

    # process accessions
    if return_type == "list":
        rfam_accs = [x[0] for x in cursor.fetchall()]
    elif return_type == "dict":
        for rfam_acc in [x[0] for x in cursor.fetchall()]:
            rfam_accs[rfam_acc] = ''

    cursor.close()
    RfamDB.disconnect(cnx)

    return rfam_accs


# ----------------------------------------------------------------------------


def fetch_taxonomy_fields(tax_id):
    """
    Fetches all fields from RfamLive taxonomy table based on
    the tax id provided

    tax_id: A valid tax id

    return: A dictionary with all taxonomy fields
    """

    cnx = RfamDB.connect()
    cursor = cnx.cursor(dictionary=True)

    query = "Select * from taxonomy where ncbi_id=%s"

    cursor.execute(query % tax_id)

    fields = cursor.fetchall()[0]

    return fields


# ----------------------------------------------------------------------------


def fetch_max_RG_accession_from_genome():
    """
    Fetches the maximum RFXXXXXXXXX accession from the RfamLive
    genome table. To be used for assigning accessions to genomes
    not found in Uniprot proteomes.

    return: Returns the maximum RGXXXXXXXXX id found in the genome
    table
    """

    cnx = RfamDB.connect()
    cursor = cnx.cursor(buffered=True)

    query = "Select max(upid) from genome where upid like \'RG%\'"

    cursor.execute(query)

    rfam_genome_id = cursor.fetchone()[0]

    return rfam_genome_id


# ----------------------------------------------------------------------------


def populate_genome_table(data):
    """
    Populates the RfamLive genome table with the data provided as input

    data: A list of tuples with the new genome table entries

    return: Void
    """

    # connect to db
    cnx = RfamDB.connect()

    # get a new buffered cursor
    cursor = cnx.cursor(raw=True)

    # query to update is_significant field to 0
    query = ("INSERT INTO genome (upid, assembly_acc, assembly_version, wgs_acc,"
             "wgs_version, assembly_name, assembly_level, study_ref, description,"
             "total_length, ungapped_length, circular, ncbi_id, scientific_name,"
             "common_name, kingdom, num_rfam_regions, num_families, is_reference,"
             "is_representative) VALUES (%s,%s,%s, %s, %s, %s, %s, %s, %s, %s,"
             "%s,%s,%s, %s, %s, %s, %s, %s, %s, %s)")

    try:
        # execute query batched
        cursor.executemany(query, data)
        cnx.commit()

    except:
        print("MySQL Update Error. Rolling back...")
        cnx.rollback()
        cursor.close()
        RfamDB.disconnect(cnx)

        cursor.close()
        RfamDB.disconnect(cnx)


# ----------------------------------------------------------------------------


def update_metagenomic_region_md5s(data):
    """
    Updates md5 fields of the seed region table

    data: A list of tuples specifying the entries to populate

    return: void
    """

    # connect to db
    cnx = RfamDB.connect()

    # get a new buffered cursor
    cursor = cnx.cursor(buffered=True)

    # update is_significant field to 0
    query = ("UPDATE meta_full_region SET md5=%s WHERE rfam_acc=%s "
             "AND rfamseq_acc=%s AND seq_start=%s AND seq_end=%s")

    cursor.executemany(query, data)

    cnx.commit()

    cursor.close()
    RfamDB.disconnect(cnx)


# ----------------------------------------------------------------------------


def fetch_family_tax_ids(rfam_acc):
    """
    Queries RfamLive and extracts all family taxonomy ids

    rfam_acc: A valid Rfam family accession

    return: A list of taxonomic ids associated with a specific Rfam family
    """

    query = """select distinct ncbi_id
    from family_ncbi
    where rfam_acc=\'%s\'"""

    cnx = RfamDB.connect()
    cursor = cnx.cursor(buffered=True)

    cursor.execute(query % rfam_acc)

    tax_ids = [x[0] for x in cursor.fetchall()]

    cursor.close()
    cnx.close()

    return tax_ids


# ----------------------------------------------------------------------------


def fetch_family_full_regions(rfam_acc, sort=True):
    """
    Fetches family regions from full_region table
    :param rfam_acc:

    :return: A dictionary with all FULL regions per accession belonging to a specific
    family
    """

    query = """select rfamseq_acc, seq_start, seq_end
        from full_region
        where rfam_acc=\'%s\'
        and is_significant=1
        and type=\'full\'"""

    cnx = RfamDB.connect()
    cursor = cnx.cursor(dictionary=True)

    cursor.execute(query % rfam_acc)

    regions = {}

    for region in cursor.fetchall():
        if region["rfamseq_acc"] not in regions:
            regions[region["rfamseq_acc"]] = [(int(region["seq_start"]), int(region["seq_end"]))]
        else:
            regions[region["rfamseq_acc"]].append((int(region["seq_start"]), int(region["seq_end"])))

    cursor.close()
    cnx.close()

    if sort is True:
        for accession in regions:
            # sorts hits by start points
            regions[accession].sort(key=lambda tup: tup[1])

    return regions


# ----------------------------------------------------------------------------


def fetch_family_seed_regions(rfam_acc):
    """
    Fetches family regions from full_region table
    :param rfam_acc:

    :return: A dictionary with all SEED regions per accession belonging to a specific
    family
    """

    query = """select rfamseq_acc, seq_start, seq_end
        from seed_region
        where rfam_acc=\'%s\'"""

    cnx = RfamDB.connect()
    cursor = cnx.cursor(dictionary=True)

    cursor.execute(query % rfam_acc)

    regions = cursor.fetchall()

    cursor.close()
    cnx.close()

    return regions


# ----------------------------------------------------------------------------


def fetch_family_metadata(rfam_acc):
    """
    Fetches family metadata from family table
    :param rfam_acc:

    :return: A dictionary with metadata describing an Rfam family family
    """

    query = """select rfam_id, description, type
        from family
        where rfam_acc=\'%s\'"""

    cnx = RfamDB.connect()
    cursor = cnx.cursor(dictionary=True)

    cursor.execute(query % rfam_acc)

    metadata = cursor.fetchone()

    cursor.close()
    cnx.close()

    return metadata


# ----------------------------------------------------------------------------


def fetch_mirna_families():
    """
    Fetches a list of all microRNA families from family table

    :return: A dictionary with metadata describing Rfam microRNA families
    """

    query = """select rfam_acc, rfam_id, description, gathering_cutoff
        from family
        where type like '%mirna%'"""

    cnx = RfamDB.connect()
    cursor = cnx.cursor(dictionary=True)

    cursor.execute(query)

    data = cursor.fetchall()

    cursor.close()
    cnx.close()

    return data

# ----------------------------------------------------------------------------
