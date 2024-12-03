"""
Update all RNAcentral descriptions to the latest version using RNAcentral API.

Usage:

python update_rnacentral_descriptions.py
"""

import requests
from utils import RfamDB


def update_description(cursor, cnx, rfamseq_acc, description):
    """
    Update rfamseq description for a sequence ID.
    """
    sql = """UPDATE rfamseq
             SET description=%s
             WHERE
             rfamseq_acc=%s
             """
    cursor.execute(sql, (description, rfamseq_acc))
    cnx.commit()


def get_rnacentral_ids(cursor):
    """
    Get a list of RNAcentral IDs from the rfamseq table.
    """
    data = []
    sql = """SELECT rfamseq_acc
             FROM rfamseq
             WHERE
             rfamseq_acc LIKE 'URS00%'"""
    cursor.execute(sql)
    for result in cursor.fetchall():
        data.append(result)
    print 'Found {} RNAcentral IDs'.format(len(data))
    return data


def update_descriptions(cursor, cnx):
    """
    Update all RNAcentral descriptions to the latest version.
    """
    found = 0
    not_found = 0
    errored = 0
    for entry in get_rnacentral_ids(cursor):
        url = 'http://www.ebi.ac.uk/ebisearch/ws/rest/rnacentral?query={} AND entry_type:sequence&fields=description&format=json'
        rnacentral_id = entry[0]
        data = requests.get(url.format(rnacentral_id))
        try:
            if data.json()['hitCount'] == 1:
                description = data.json()['entries'][0]['fields']['description'][0]
                print('{}: {}'.format(rnacentral_id, description))
                update_description(cursor, cnx, rnacentral_id, description)
                found += 1
            else:
                print('No description found for {}'.format(rnacentral_id))
                not_found += 1
        except (ValueError, KeyError) as e:
            print("Server Error: {err} for URL: {url}".format(err=e, url=url.format(rnacentral_id)))
            errored += 1
    print('Updated {} descriptions, not found {} descriptions'.format(found, not_found))


def main():
    """
    Main entry point.
    """
    cnx = RfamDB.connect()
    cursor = cnx.cursor(buffered=True)
    update_descriptions(cursor, cnx)
    cursor.close()
    RfamDB.disconnect(cnx)


if __name__ == '__main__':
    main()
