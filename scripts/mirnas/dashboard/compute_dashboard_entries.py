import os
import re

import mysql.connector
import requests

from utils import RfamDB
from microrna_progress import updated_families, new_commits
from dashboard_exceptions import SpeciesFileNotFound
from getters import get_family_location, get_report_url, get_rfam_id, get_mirbase_id

SEARCH_DIR = '/nfs/production/agb/rfam/microrna/searches'


def run_rfmake(location, score):
    """
    Run rfmake with a manually selected threshold if not done already.
    """
    if not score or '?' in score:
        return
    score = float(score)
    species_file = os.path.join(location, 'species')
    if not os.path.exists(species_file):
        raise SpeciesFileNotFound('Species file not found in {}'.format(location))
    ga_threshold = None
    with open(species_file, 'r') as f_in:
        for line in f_in:
            match = re.search(r'CURRENT GA THRESHOLD: (.+) BITS', line)
            if match:
                ga_threshold = float(match.group(1))
                break
    if ga_threshold and ga_threshold == score:
        pass
    else:
        cmd = 'cd {} && rfmake.pl -t {} && cd - > /dev/null'.format(location, score)
        os.system(cmd)


def get_overlaps(identifier, nocache):
    """
    Run rqc-overlaps and get overlapping families.
    External overlap [SS] of CM001599.2/64499845-64499913 with RF00600:CM001599.2/64499850-64499995 by 64
    """
    location = get_family_location(identifier)
    overlap_file = os.path.join(location, 'overlap')
    rfam_accs = set()
    if not os.path.exists(overlap_file) or os.stat(overlap_file).st_size == 0 or nocache:
        cmd = 'cd {} && touch SEED CM DESC TBLOUT SCORES && cd - > /dev/null'.format(location)
        os.system(cmd)
        cmd = 'cd {} && rqc-overlap.pl {} && cd - > /dev/null'.format(SEARCH_DIR, os.path.basename(location))
        os.system(cmd)
    with open(overlap_file, 'r') as f_overlap:
        for line in f_overlap:
            match = re.search(r'with (RF\d{5})', line)
            if match:
                rfam_accs.add(match.group(1))
    return list(rfam_accs)


def is_single_seed(location):
    """
    Check if the seed contains only 1 sequence.
    """
    seed_file = os.path.join(location, 'SEED')
    alistat_file = seed_file + '_alistat.txt'
    cmd = 'esl-alistat {} > {}'.format(seed_file, alistat_file)
    os.system(cmd)
    with open(alistat_file, 'r') as f_alistat:
        for line in f_alistat:
            if line.startswith('Number of sequences:'):
                fields = re.split(r'\:\s+', line)
                if len(fields) == 2:
                    num_seed = int(fields[1])
                    if num_seed == 1:
                        return True
                    else:
                        return False
    return None


def is_inconsistent_ss(location):
    """
    Check if seed alignment has an inconsistent secondary structure.
    """
    seed_file = os.path.join(location, 'SEED')
    alistat_file = seed_file + '_alistat.txt'
    cmd = 'esl-alistat --bpinfo test.bp.txt {} > {} 2>&1'.format(seed_file, alistat_file)
    os.system(cmd)
    with open(alistat_file, 'r') as f_alistat:
        for line in f_alistat:
            if line.startswith('Consensus structure string is inconsistent'):
                return True
    return False


def get_new_or_updated(overlaps):
    """
    Check if this is one of the microRNA family that has already been updated or
    newly committed.
    """
    if len(overlaps) == 1:
        rfam_acc = overlaps[0]
        if rfam_acc in updated_families:
            return 'Updated'
        elif rfam_acc in new_commits:
            return 'New'
    return 'False'


def get_id_matches(mirbase_id):
    """
    In case there are no overlaps, check for ID matches.
    For newly created families (not yet in the search index), run a query against RfamLive.
    """
    url = 'http://www.ebi.ac.uk/ebisearch/ws/rest/rfam?query="{}"%20AND%20entry_type:%22Family%22&fields=name&format=json'
    data = requests.get(url.format(mirbase_id))
    if data.json()['hitCount'] == 1:
        rfam_id = data.json()['entries'][0]['id']
        return [rfam_id]
    if data.json()['hitCount'] == 0:
        query = "SELECT rfam_acc FROM family WHERE rfam_id = '{}'".format(mirbase_id)
        conn = RfamDB.connect()
        cursor = conn.cursor()
        try:
            cursor.execute(query)
            rfam_id = cursor.fetchone()
            return [rfam_id]
        except mysql.connector.Error as e:
            print("MySQL error has occurred: {0}".format(e))
            raise e
        except Exception as e:
            raise e
        finally:
            cursor.close()
    return []


def get_action(identifier, location, overlaps, overlaps_by_id, score):
    """
    Identify what should be done with a family.
    """
    action = ''
    if overlaps:
        overlap_status = get_new_or_updated(overlaps)
    elif overlaps_by_id:
        overlap_status = get_new_or_updated(overlaps_by_id)
    else:
        overlap_status = None

    mirbase_id = get_mirbase_id(identifier)
    if len(overlaps_by_id) == 1 and len(overlaps) == 1 and overlaps[0] != overlaps_by_id[0]:
        return 'Fix ID mismatch'
    elif len(overlaps_by_id) == 1 and len(overlaps) == 1 and overlaps[0] == overlaps_by_id[0]:
        rfam_id = get_rfam_id(overlaps[0])
    elif len(overlaps) == 1 and not overlaps_by_id:
        rfam_id = get_rfam_id(overlaps[0])
    elif len(overlaps_by_id) == 1 and not overlaps:
        rfam_id = get_rfam_id(overlaps_by_id[0])

    if 'HYPERLINK' not in get_report_url(identifier):
        action = 'Generate report'
    elif is_single_seed(location):
        action = '1_SEED'
    elif not score or score == '?':
        action = 'Choose threshold'
    elif overlap_status == 'New' and rfam_id.lower() == mirbase_id.lower():
        action = 'Done (new family)'
    elif overlap_status == 'Updated' and rfam_id.lower() == mirbase_id.lower():
        action = 'Done (updated family)'
    elif score > 0 and not overlaps and not overlaps_by_id:
        action = 'New family'
    elif score > 0 and (overlaps or overlaps_by_id):
        action = 'Update seed'
    return action
