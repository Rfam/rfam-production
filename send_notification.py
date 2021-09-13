import requests

from config.rfam_local import SLACK_WEBHOOK


def send_notification():
    """
    Send notification to Slack channel using incoming webhook
    """
    slack_message = ""
    webhook_url = SLACK_WEBHOOK

    with open('pdb_families.txt', 'r') as f:
        for line in f:
            slack_message += line
    slack_json = {
        "text": slack_message,
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
    response = requests.post(webhook_url, json=slack_json, headers={'Content-Type': 'application/json'})
    if response.status_code != 200:
        raise ValueError("Error with request {0}, the response is:\n{1}".format(response.status_code, response.text))


if __name__ == '__main__':
    send_notification()
