import glob
import os

import requests

from config.rfam_local import SLACK_WEBHOOK, PDB_FILES


def send_notification():
    """
    Send notification to Slack channel using incoming webhook
    """
    webhook_url = SLACK_WEBHOOK
    # Use the newest report rather than today's, so a run that ends after midnight still notifies.
    reports = sorted(glob.glob(os.path.join(PDB_FILES, 'pdb_families_*.txt')))
    if not reports:
        raise SystemExit('No pdb_families_*.txt report found in {0}'.format(PDB_FILES))
    report = reports[-1]
    report_date = os.path.basename(report)[len('pdb_families_'):-len('.txt')]
    slack_message = 'Report date: {0}\n'.format(report_date)
    with open(report, 'r') as f:
        for line in f:
            slack_message += line
    slack_json = {
        "text": "PDB Mapping",
        "blocks": [
            {
                "type": "section",
                "text": {
                    "type": "mrkdwn",
                    "text": slack_message
                },
            },
        ]
    }
    try:
        response = requests.post(webhook_url, json=slack_json, headers={'Content-Type': 'application/json'})
        response.raise_for_status()
    except requests.exceptions.HTTPError as e:
        raise SystemExit(e)
    except requests.exceptions.RequestException as e:
        raise SystemExit(e)
    except Exception as e:
        raise SystemExit(e)


if __name__ == '__main__':
    send_notification()
