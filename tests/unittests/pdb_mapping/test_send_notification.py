import unittest
from unittest.mock import patch, mock_open, Mock

from requests import HTTPError

from pdb_mapping.send_notification import send_notification


class SendNotificationTest(unittest.TestCase):

    @patch("pdb_mapping.send_notification.requests.post")
    @patch("pdb_mapping.send_notification.SLACK_WEBHOOK")
    def test_send_notification_is_successful(self, mock_hook, mock_post):
        mock_hook.return_value = "/mock/url/webhook"
        with patch("builtins.open", mock_open(read_data="data")) as mock_file:
            send_notification()
            assert open("rfam-production/pdb_families.txt").read() == "data"
            mock_file.assert_called_with("rfam-production/pdb_families.txt")
            mock_post.assert_called_once()

    @patch("pdb_mapping.send_notification.requests.post")
    @patch("pdb_mapping.send_notification.SLACK_WEBHOOK")
    def test_send_notification_raises_system_exit_on_http_error(self, mock_hook, mock_post):
        mock_hook.return_value = "/mock/url/webhook"
        resp = mock_post.return_value
        resp.raise_for_status. side_effect = HTTPError()
        with patch("builtins.open", mock_open(read_data="data")):
            with self.assertRaises(SystemExit):
                send_notification()

    @patch("pdb_mapping.send_notification.requests.post")
    @patch("pdb_mapping.send_notification.SLACK_WEBHOOK")
    def test_send_notification_raises_system_exit_on_exception(self, mock_hook, mock_post):
        mock_hook.return_value = "/mock/url/webhook"
        resp = mock_post.return_value
        resp.raise_for_status. side_effect = Exception()
        with patch("builtins.open", mock_open(read_data="data")):
            with self.assertRaises(SystemExit):
                send_notification()



