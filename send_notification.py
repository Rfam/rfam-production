import json
import os

import requests


def send_notification():
    """
    Send notification to Slack channel using incoming webhook
    """
    file_data = {}

    webhook_url = os.environ('SLACK_WEBHOOK')

    with open('pdb_families.txt', 'r') as f:
        for line in f:
            description, result = line.strip().split(None, 1)
            file_data[description] = result.strip()

    data = json.dumps(file_data)

    response = requests.post(webhook_url, json={"text": data}, headers={'Content-Type': 'application/json'})
    if response.status_code != 200:
        raise ValueError('Error with request {}, the response is:\n{}'.format(response.status_code, response.text))


if __name__ == '__main__':
    send_notification()
