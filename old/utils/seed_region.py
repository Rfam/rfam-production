import RfamDB

# -------------------------------------------------------------------------


def update_seed_region_md5s(data):

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
    query = (
        "UPDATE seed_region SET md5=%s WHERE rfam_acc=%s "
        "AND rfamseq_acc=%s AND seq_start=%s AND seq_end=%s"
    )

    cursor.executemany(query, data)

    cnx.commit()

    cursor.close()
    RfamDB.disconnect(cnx)


# -------------------------------------------------------------------------


def fetch_seed_regions():
    """
    Fetches all seed_region entries

    return: A list of tuples with all seed_region entries
    """

    # connect to db
    cnx = RfamDB.connect()

    # get a new buffered cursor
    cursor = cnx.cursor(buffered=True)

    # update is_significant field to 0
    query = "Select rfam_acc, rfamseq_acc, seq_start, seq_end " "from seed_region"

    cursor.execute(query)

    seed_region_rows = cursor.fetchall()

    cursor.close()
    RfamDB.disconnect(cnx)

    return seed_region_rows


# -------------------------------------------------------------------------
