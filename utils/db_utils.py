'''
Created on 27 Jan 2016

@author: ikalvari

Description: A set of database functions to ease processing and data
             retrieval from rfam_live

TO DO: modify reset_is_significant() to enable single clan reset
'''
# ---------------------------------IMPORTS---------------------------------

import string
from utils import RfamDB

# -------------------------------------------------------------------------
RFAM_ACC = 0  # full region rfam_acc
SEQ_ACC = 1  # full region rfamseq_acc
START = 2  # full region seq_start
END = 3  # full region seq_end
EVAL = 4  # full region evalue

# -------------------------------------------------------------------------


def set_is_significant_to_zero(rfam_acc, rfamseq_acc):
    '''
        Fetch the correct db entry from full_region table according to
        rfam_acc and rfamseq_acc and set is_significant field to zero (0)

        rfam_acc: RNA family accession
        rfamseq_acc: Family specific sequence accession

    '''

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
    '''
        Fetch the correct db entry from full_region table according to
        rfam_acc and rfamseq_acc and set is_significant field to zero (0)

        rfam_acc: RNA family accession
        rfamseq_acc: Family specific sequence accession
    '''

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
    '''
        Loads specific clan family sequences from full_region table and returns 
        a dictionary structure as {Rfam_acc:{Rfseq_acc:[start, end, evalue]}}
        for clan competition.

        This has been modified to accommodate sequence duplicates

        clan_acc: Clan accession as in Rfam
    '''

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
    '''
        Retrieves all clan family members from DB and returns a list of the family
        accessions.

        clan_acc: Clan accession as in Rfam
    '''

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
    '''
        Retrieves all clan family members from DB and returns a dictionary in the
        in the form of {clan_id : [FAM1, FAM2, ... ], ... }

        clan_acc: Clan accession as in Rfam
    '''

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
    '''
        A function for batching the process of updating full_region tables upon
        clan competition. Updates the full_region table setting is_significant
        field to zero (0) for the list of non significant sequences passed in
        the form of (rfam_acc, rfamseq_acc, seq_start) tuples.

        non_sig_seqs: A list of the non significant regions to be set to zero.
                      The list is product of clan competition.

    '''

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


def reset_is_significant():
    '''
        This function resets full_region's is_singificant field's back to 1.
        This should be able to update all or part of the table for clan
        competition initialization and restoration.
    '''

    seq_regs = []

    cnx = RfamDB.connect()

    # cursor to fetch data
    d_cursor = cnx.cursor(buffered=True)

    # query to fetch all non significant sequences
    s_query = ("SELECT rfam_acc, rfamseq_acc, seq_start FROM full_region "
               "WHERE is_significant=0")

    # query to update 0 fields from s_query
    query = ("UPDATE full_region SET is_significant=1 "
             "WHERE rfam_acc=%s AND rfamseq_acc=%s AND seq_start=%s")

    d_cursor.execute(s_query)

    # construct region list here
    for row in d_cursor:
        seq_regs.append((str(row[0]), str(row[1]), int(row[2])))

    d_cursor.close()

    # get a new cursor for db updates
    u_cursor = cnx.cursor(raw=True)

    # update db
    try:
        u_cursor.executemany(query, seq_regs)
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
    '''
        Updates _post_process table with the job_ids per family assigned by lsf

        jobs_file: This is a tab separated txt file generated from running the
        job_dequeuer.py script that submits the rfam_view_process for each
        family.
        (rfam_acc uuid job_id ...)
    '''

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

if __name__ == '__main__':

    # reset_is_significant()

    jobs_file = "/Users/ikalvari/Desktop/Rfam12.1/fam_view_process/view_out/famview_jobs.txt"
    update_post_process(jobs_file)
