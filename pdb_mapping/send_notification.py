import requests
import datetime

from config.rfam_local import SLACK_WEBHOOK


def send_notification():
    """
    Send notification to Slack channel using incoming webhook
    """
    slack_message = ""
    webhook_url = SLACK_WEBHOOK
    today_date = str(datetime.date.today())
    with open('pdb_mapping/pdb_families_{0}.txt'.format(today_date), 'r') as f:
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
