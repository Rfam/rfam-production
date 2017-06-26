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
Created on 27 Jan 2016

@author: ikalvari

Description: A set of database functions to ease processing and data
             retrieval from rfam_live

TO DO: - modify reset_is_significant() to enable single clan reset
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

import string
import sys
from utils import RfamDB

# -------------------------------------------------------------------------
RFAM_ACC = 0  # full region rfam_acc
SEQ_ACC = 1  # full region rfamseq_acc
START = 2  # full region seq_start
END = 3  # full region seq_end
EVAL = 4  # full region evalue

# -------------------------------------------------------------------------


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

        if(str(row[RFAM_ACC]) in fam_seqs.keys()):

            if (str(row[SEQ_ACC]) in fam_seqs[str(row[RFAM_ACC])].keys()):

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

    query = ("SELECT * FROM clan_membership")

    # execute query
    cursor.execute(query)

    # fetch the data
    rows = cursor.fetchall()

    cursor.close()
    RfamDB.disconnect(cnx)

    # create the dictionary
    for row in rows:
        if (str(row[0]) not in clans.keys()):
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
        print "MySQL Update Error. Rolling back..."
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
    s_query = ''
    query = ''
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
        print "MySQL Update Error. Rolling back..."
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

    fp = open(jobs_file, 'r')

    query = ("UPDATE _post_process SET lsf_id=%s "
             "WHERE rfam_acc=%s AND uuid=%s")

    # get lsf ids from file
    for line in fp:
        line = line.strip()
        line = string.split(line, '\t')
        job_ids.append((line[2], line[0], line[1]))

    # connect to db
    cnx = RfamDB.connect()
    cursor = cnx.cursor(raw=True)

    # update db
    try:
        cursor.executemany(query, job_ids)
        cnx.commit()  # move this after except statement??

    except:
        # rollback to previous state
        print "MySQL Update Error. Rollback..."
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

    # counts list
    counts = []
    for acc in rfam_accs:
        c_cursor.execute(count_query % str(acc[0]))
        count = c_cursor.fetchall()

        counts.append((count[0][0], str(acc[0])))

        count = 0

    c_cursor.close()
    c_cursor = cnx.cursor(buffered=True)

    # query to update number_of_species in the family table
    update_query = (
        "update family set number_of_species=%s where rfam_acc=%s")

    try:
        c_cursor.executemany(update_query, counts)
        cnx.commit()
    except:
        cnx.rollback()

    c_cursor.close()
    RfamDB.disconnect(cnx)

    print "Done"
# ----------------------------------------------------------------------------


def set_num_sig_seqs():
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
                   "and rfam_acc=\'%s\'")

    # counts list
    counts = []
    for acc in rfam_accs:
        c_cursor.execute(count_query % str(acc[0]))
        count = c_cursor.fetchall()

        counts.append((count[0][0], str(acc[0])))

        count = 0

    c_cursor.close()
    c_cursor = cnx.cursor(buffered=True)

    update_query = (
        "update family set num_full=%s where rfam_acc=%s")

    try:
        c_cursor.executemany(update_query, counts)
        cnx.commit()
    except:
        cnx.rollback()

    c_cursor.close()
    RfamDB.disconnect(cnx)

    print "Done"

# ----------------------------------------------------------------------------


def update_family_ncbi():
    """
    Updates table family ncbi by adding all distinct taxonomic ids per family

    :return: void
    """

    cnx = RfamDB.connect()

    cursor = cnx.cursor(buffered=True)
    c_cursor = cnx.cursor(buffered=True)

    cursor.execute("Select rfam_acc from family")

    rfam_accs = cursor.fetchall()

    cursor.close()

    # family_ncbi query
    get_ncbi_ids = ("select distinct rs.ncbi_id, f.rfam_id, f.rfam_acc from full_region fr, rfamseq rs, family f "
                    "where fr.rfamseq_acc=rs.rfamseq_acc "
                    "and f.rfam_acc=fr.rfam_acc "
                    "and fr.rfam_acc=\'%s\' "
                    "and fr.is_significant=1")

    insert_query = "insert into family_ncbi (ncbi_id, rfam_id, rfam_acc) values (%s,%s,%s)"

    family_ncbi_entries = []
    cursor = cnx.cursor(buffered=True)
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

    print "Done updating family_ncbi."

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


def set_pdb_is_significant_to_zero(non_sig_seqs):
    """
    Sets pdb_full_region is_significant to 0 for non significant regions in
    non_sig_seqs list

    non_sig_seqs: A list of the non significant regions to be set to zero.
    The list is product of clan competition.

    returns: void
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
    cnx = RfamDB.connect()

    # get a new buffered cursor
    cursor = cnx.cursor(raw=True)

    # query to update is_significant field to 0
    query = ("update pdb_full_region set is_significant=0 "
             "where rfam_acc=%s and pdb_id=%s and chain=%s and pdb_start=%s")

    try:
        # execute query batched
        cursor.executemany(query, pdb_reformatted_regions)
        cnx.commit()

    except:
        print "MySQL Update Error. Rolling back..."
        cnx.rollback()
        cursor.close()
        RfamDB.disconnect(cnx)

    cursor.close()
    RfamDB.disconnect(cnx)

# ----------------------------------------------------------------------------

if __name__ == '__main__':

    pass


