"""
Tests for bulk_upload.py. No network or svn access; nothing is uploaded.

Run from the repository root:
    python3 -m unittest scripts.apicuron.tests.test_bulk_upload -v
"""

import json
import logging
import os
import sys
import tempfile
import types
import unittest
from unittest import mock

import requests

# Stand-in for config.rfam_config, so the tests need no rfam_local.py.
fake_config = types.ModuleType("config.rfam_config")
fake_config.APICURON_TOKEN = "test-token"
fake_config.BULK_REPORT_URL = "https://apicuron.invalid/api/reports/bulk"
sys.modules.setdefault("config", types.ModuleType("config"))
sys.modules["config"].rfam_config = fake_config
sys.modules["config.rfam_config"] = fake_config

from scripts.apicuron import bulk_upload as bu  # noqa: E402

NANCY = "0000-0001-8457-4455"
BLAKE = "0000-0002-6497-2883"
ANTON = "0000-0001-7279-2682"
BAD_ORCID = "0000-0002-6497-2884"  # BLAKE with a wrong check digit
FAMILY = "/trunk/Families/RF00001/DESC"
URI = "https://rfam.org/family/RF00001"


def report(timestamp, orcid=NANCY, term="update_family", uri=URI):
    return {"activity_term": term, "curator_orcid": orcid,
            "entity_uri": uri, "timestamp": timestamp}


def revision(rev, author, date, message, paths=FAMILY):
    return (rev, author, date, message, paths)


def args(**kwargs):
    values = dict(start=None, end=None, svn="/svn", rfamprod_orcid=NANCY)
    values.update(kwargs)
    return types.SimpleNamespace(**values)


class InTempDir(unittest.TestCase):
    """Run each test in an empty folder, as the script writes to the cwd."""

    def setUp(self):
        self.old_cwd = os.getcwd()
        self.tmp = tempfile.TemporaryDirectory()
        os.chdir(self.tmp.name)
        logging.disable(logging.CRITICAL)

    def tearDown(self):
        logging.disable(logging.NOTSET)
        os.chdir(self.old_cwd)
        self.tmp.cleanup()

    def load(self, name):
        with open(name) as handle:
            return json.load(handle)["reports"]


class TestOrcids(unittest.TestCase):

    def test_checksum(self):
        self.assertTrue(bu.is_valid_orcid(NANCY))
        self.assertTrue(bu.is_valid_orcid("0000-0001-6043-807X"))
        self.assertFalse(bu.is_valid_orcid(BAD_ORCID))

    def test_suggested_format(self):
        message = "CI:fix SS_cons by %s and %s" % (NANCY, BLAKE)
        self.assertEqual(bu.get_orcids(message), ([NANCY, BLAKE], []))

    def test_three_orcids(self):
        message = "CI:x by %s and %s and %s" % (NANCY, BLAKE, ANTON)
        self.assertEqual(bu.get_orcids(message), ([NANCY, BLAKE, ANTON], []))

    def test_url_form_and_repeats(self):
        message = "CI:x by https://orcid.org/%s and %s" % (NANCY, NANCY)
        self.assertEqual(bu.get_orcids(message), ([NANCY], []))

    def test_no_orcid(self):
        self.assertEqual(bu.get_orcids("CI:update GO terms"), ([], []))

    def test_bad_checksum_separated(self):
        message = "CI:x by %s and %s" % (NANCY, BAD_ORCID)
        self.assertEqual(bu.get_orcids(message), ([NANCY], [BAD_ORCID]))

    def test_bad_correction_stops(self):
        bu.check_corrections({1: [NANCY], 2: []})
        with self.assertRaisesRegex(ValueError, "r3: invalid ORCID"):
            bu.check_corrections({3: [BAD_ORCID]})


class TestParsing(unittest.TestCase):

    def test_timestamp_is_utc(self):
        self.assertEqual(bu.get_timestamp("2024-07-24T06:36:59.225407Z"),
                         "2024-07-24T06:36:59Z")
        self.assertEqual(bu.get_timestamp(""), "")

    def test_family(self):
        self.assertEqual(bu.get_family("/trunk/Families/RF02540"),
                         "https://rfam.org/family/RF02540")
        self.assertEqual(bu.get_family("/trunk/FamiliesPending/x"), "")

    def test_svn_url(self):
        self.assertEqual(bu.svn_url("/a/repo"), "file:///a/repo")
        self.assertEqual(bu.svn_url("https://svn/repo"), "https://svn/repo")

    def test_svn_error_is_in_exception(self):
        failed = mock.Mock(returncode=1, stderr=b"svn: E170013: no repo")
        with mock.patch("subprocess.run", return_value=failed):
            with self.assertRaisesRegex(RuntimeError, "E170013"):
                bu.get_youngest("/missing")

    def test_get_revisions(self):
        xml = b"""<?xml version="1.0"?><log>
          <logentry revision="5004"><author>nawrocki</author>
            <date>2014-06-11T15:02:31.1Z</date>
            <paths><path action="A">/trunk/Families/RF02540</path></paths>
            <msg>NEWMOV:Moving from pending</msg></logentry></log>"""
        with mock.patch.object(bu, "run_svn", return_value=xml):
            self.assertEqual(bu.get_revisions(5004, 5005, "/svn"), [
                (5004, "nawrocki", "2014-06-11T15:02:31.1Z",
                 "NEWMOV:Moving from pending", "/trunk/Families/RF02540")])


class TestHeader(unittest.TestCase):

    def test_bearer_added_once(self):
        for token in ("abc", " abc\n", "Bearer abc"):
            with mock.patch.object(bu.rfam_config, "APICURON_TOKEN", token):
                self.assertEqual(bu.get_header()["authorization"], "Bearer abc")

    def test_missing_token(self):
        with mock.patch.object(bu.rfam_config, "APICURON_TOKEN", ""):
            with self.assertRaises(RuntimeError):
                bu.get_header()


class TestRemoveExisting(unittest.TestCase):

    def test_exact_and_one_hour_shift(self):
        reports = [report("2022-01-14T15:35:31Z"), report("2023-01-01T10:00:00Z")]
        existing = [report("2022-01-14T14:35:31.000Z")]  # 1 h early, old upload
        self.assertEqual(bu.remove_existing(reports, existing), ([True, False], 0, 1))

    def test_outside_tolerance_or_not_whole_hours(self):
        reports = [report("2024-01-01T10:00:00Z"), report("2024-01-02T10:00:00Z")]
        existing = [report("2024-01-01T13:00:00Z"), report("2024-01-02T10:30:00Z")]
        self.assertEqual(bu.remove_existing(reports, existing), ([False, False], 0, 0))

    def test_other_curator_family_or_term(self):
        reports = [report("2024-01-01T10:00:00Z")]
        for other in (report("2024-01-01T10:00:00Z", orcid=BLAKE),
                      report("2024-01-01T10:00:00Z", uri=URI[:-1] + "2"),
                      report("2024-01-01T10:00:00Z", term="create_family")):
            self.assertEqual(bu.remove_existing(reports, [other]), ([False], 0, 0))

    def test_each_existing_report_matches_once(self):
        # Two real commits 12 h apart, both on APICURON with exact times.
        reports = [report("2024-09-07T15:58:00Z"), report("2024-09-08T03:58:00Z")]
        self.assertEqual(bu.remove_existing(reports, list(reports)), ([True, True], 2, 0))
        # Only one on APICURON: the other must not be skipped via a shift.
        reports = [report("2024-01-01T10:00:00Z"), report("2024-01-01T11:00:00Z")]
        existing = [report("2024-01-01T10:00:00Z")]
        self.assertEqual(bu.remove_existing(reports, existing), ([True, False], 1, 0))


REVISIONS = [
    revision(1, "rfamprod", "2026-11-12T09:00:00.1Z", "CI:x by %s" % BLAKE),
    revision(2, "nancyontiveros", "2026-11-12T10:00:00.1Z",
             "CI:x by %s and %s" % (BLAKE, ANTON)),
    revision(3, "apetrov", "2026-11-12T11:00:00.1Z", "CI:x"),
    revision(4, "rfamprod", "2010-01-01T12:00:00.1Z", "CI:x"),
    revision(5, "swb", "2026-11-12T13:00:00.1Z", "CI:x"),
    revision(6, "apetrov", "2026-11-12T14:00:00.1Z", "no activity term"),
    revision(7, "apetrov", "2026-11-12T15:00:00.1Z", "CI:x", "/trunk/README"),
]


class WriteHelper(InTempDir):
    """Run write_report on fake commits and read the files back."""

    def write(self, revisions=REVISIONS, existing=None, corrections=None, **kwargs):
        with mock.patch.object(bu, "get_revisions",
                               lambda start, end, svn: [r for r in revisions if start <= r[0] < end]), \
             mock.patch.object(bu, "get_youngest", lambda svn: max(r[0] for r in revisions)), \
             mock.patch.object(bu.conf, "orcid_corrections", corrections or {}):
            name, self.warnings = bu.write_report(args(**kwargs), existing)
        self.invalid = self.load("invalid_orcid_" + name)
        return self.load(name), self.load("others_" + name)


class TestWriteReport(WriteHelper):

    def test_curator_priority(self):
        upload, others = self.write(start=1, end=8)
        self.assertEqual([r["curator_orcid"] for r in upload], [
            BLAKE,           # 1: message beats the rfamprod fallback
            BLAKE, ANTON,    # 2: message beats the svn user; one report per ORCID
            ANTON,           # 3: svn user
            NANCY,           # 4: rfamprod fallback, no date cutoff
        ])
        self.assertEqual([r["curator_orcid"] for r in others], ["swb"])  # 5; 6, 7 skipped

    def test_rfamprod_fallback_off(self):
        upload, others = self.write(start=4, end=5, rfamprod_orcid="")
        self.assertEqual((upload, [r["curator_orcid"] for r in others]), ([], ["rfamprod"]))

    def test_bad_orcid_held_back(self):
        bad = [
            revision(1, "apetrov", "2026-11-12T09:00:00.1Z", "CI:x by " + BAD_ORCID),
            revision(2, "rfamprod", "2026-11-12T10:00:00.1Z", "CI:x by %s and %s" % (NANCY, BAD_ORCID)),
        ]
        logging.disable(logging.NOTSET)  # InTempDir turns logging off
        with self.assertLogs(level="WARNING") as logs:
            upload, others = self.write(revisions=bad, start=1, end=3)
        # r1: no fallback to the svn user (Anton); r2: the valid ORCID is uploaded.
        self.assertEqual([r["curator_orcid"] for r in upload], [NANCY])
        self.assertEqual(others, [])
        self.assertEqual([(r["revision"], r["curator_orcid"]) for r in self.invalid],
                         [(1, BAD_ORCID), (2, BAD_ORCID)])
        self.assertEqual(self.warnings, 2)
        # The warning suggests the line to add, keeping the valid ORCID.
        self.assertIn("2: ['%s', '<correct ORCID>'],  # message had %s" % (NANCY, BAD_ORCID),
                      "\n".join(logs.output))
        self.assertIn("2 commits need an ORCID correction: r1, r2", "\n".join(logs.output))

    def test_correction_replaces_message_orcids(self):
        bad = [revision(1, "apetrov", "2026-11-12T09:00:00.1Z", "CI:x by %s and %s" % (NANCY, BAD_ORCID))]
        upload, _ = self.write(revisions=bad, start=1, end=2, corrections={1: [NANCY, BLAKE]})
        self.assertEqual([r["curator_orcid"] for r in upload], [NANCY, BLAKE])
        self.assertEqual((self.invalid, self.warnings), ([], 0))

    def test_empty_correction_skips_commit(self):
        bad = [revision(1, "apetrov", "2026-11-12T09:00:00.1Z", "CI:x by " + BAD_ORCID)]
        upload, others = self.write(revisions=bad, start=1, end=2, corrections={1: []})
        self.assertEqual((upload, others, self.invalid, self.warnings), ([], [], [], 0))

    def test_bad_correction_stops(self):
        with self.assertRaisesRegex(ValueError, "orcid_corrections r1"):
            self.write(start=1, end=2, corrections={1: [BAD_ORCID]})

    def test_corrected_old_commit_uploaded_with_auto_start(self):
        revisions = [
            revision(1, "swb", "2026-01-01T09:00:00.1Z", "CI:x by " + BAD_ORCID),
            revision(2, "apetrov", "2026-11-12T10:00:00.1Z", "CI:x"),
        ]
        on_apicuron = [report("2026-11-12T10:00:00Z", orcid=ANTON)]  # r2
        upload, _ = self.write(revisions=revisions, existing=on_apicuron)
        self.assertEqual((upload, self.warnings), ([], 1))  # r1 held back
        upload, _ = self.write(revisions=revisions, existing=on_apicuron,
                               corrections={1: [BLAKE]})
        # Start is r2, but corrected r1 is still uploaded.
        self.assertEqual([(r["timestamp"], r["curator_orcid"]) for r in upload],
                         [("2026-01-01T09:00:00Z", BLAKE)])
        self.assertEqual(self.warnings, 0)

    def test_auto_end(self):
        upload, _ = self.write(start=3)  # newest revision is 7
        self.assertEqual(len(upload), 2)  # r3 and r4

    def test_auto_start_after_newest_on_apicuron(self):
        on_apicuron = [report("2026-11-12T10:00:00Z", orcid=BLAKE)]  # r2, first ORCID only
        upload, others = self.write(existing=on_apicuron)
        # Starts at r2 (included, half uploaded), skips what is there; r1 is older.
        self.assertEqual([(r["timestamp"], r["curator_orcid"]) for r in upload], [
            ("2026-11-12T10:00:00Z", ANTON),
            ("2026-11-12T11:00:00Z", ANTON),
            ("2010-01-01T12:00:00Z", NANCY),
        ])
        self.assertEqual([r["curator_orcid"] for r in others], ["swb"])

    def test_auto_start_needs_a_match(self):
        with self.assertRaisesRegex(RuntimeError, "give --start"):
            self.write(existing=[])

    def test_given_start_skips_existing(self):
        on_apicuron = [report("2026-11-12T10:00:00Z", orcid=BLAKE),
                       report("2026-11-12T10:00:00Z", orcid=ANTON)]
        upload, _ = self.write(start=1, end=8, existing=on_apicuron)
        self.assertEqual([r["curator_orcid"] for r in upload], [BLAKE, ANTON, NANCY])


class TestCorrectionPushedSeparately(WriteHelper):
    """An admin uploaded the right ORCID outside this script."""

    REVISIONS = [
        revision(1, "swb", "2026-11-20T09:00:00.4Z", "CI:x by %s and %s" % (NANCY, BAD_ORCID)),
        revision(2, "apetrov", "2026-11-21T09:00:00.4Z", "CI:x"),
    ]
    ALREADY = [report("2026-11-20T09:00:00.000Z", orcid=NANCY),   # r1, valid ORCID
               report("2026-11-21T09:00:00.000Z", orcid=ANTON)]   # r2
    PUSHED = report("2026-11-20T09:00:00Z", orcid=BLAKE)          # r1, by the admin

    def check(self, existing, corrections):
        upload, _ = self.write(revisions=self.REVISIONS, existing=existing,
                               corrections=corrections)
        return [r["curator_orcid"] for r in upload], self.warnings

    def test_pushed_and_corrected(self):
        self.assertEqual(self.check(self.ALREADY + [self.PUSHED], {1: [NANCY, BLAKE]}), ([], 0))

    def test_pushed_but_not_corrected(self):
        # No duplicate, but the script cannot know it was fixed: it keeps warning.
        self.assertEqual(self.check(self.ALREADY + [self.PUSHED], {}), ([], 1))

    def test_corrected_not_pushed(self):
        self.assertEqual(self.check(self.ALREADY, {1: [NANCY, BLAKE]}), ([BLAKE], 0))

    def test_pushed_with_other_time_is_uploaded_again(self):
        other_time = dict(self.PUSHED, timestamp="2026-11-25T14:00:00Z")
        self.assertEqual(self.check(self.ALREADY + [other_time], {1: [NANCY, BLAKE]}), ([BLAKE], 0))


DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

SNAPSHOT_REVISIONS = [
    revision(5004, "nawrocki", "2014-06-11T15:02:31.229894Z",
             "NEWMOV:Moving from pending to main repository.",
             "/trunk/Families/RF02540 /trunk/FamiliesPending/LSU_rRNA_archaea"),
    revision(16269, "nancyontiveros", "2023-08-08T11:13:25.746089Z",
             "CI:update GO terms", "/trunk/Families/RF04247/DESC"),
    revision(26149, "rfamprod", "2026-11-13T09:00:00.5Z",
             "CI:fix SS_cons by %s and https://orcid.org/%s" % (BLAKE, ANTON),
             "/trunk/Families/RF00005/SEED /trunk/Families/RF00005/DESC"),
    revision(26150, "rfamprod", "2026-11-13T10:00:00Z",
             "KILL:dead family", "/trunk/Families/RF00010"),
    revision(26151, "rfamprod", "2026-11-13T11:00:00Z",
             "CLCI:update clan", "/trunk/Families/RF00020/DESC"),
    revision(26152, "joanna", "2026-11-13T12:00:00Z",
             "NEW:new family", "/trunk/Families/RF04300"),
    revision(26153, "rfamprod", "2026-11-13T13:00:00Z",
             "CI:x by %s and %s" % (BLAKE, BAD_ORCID), "/trunk/Families/RF00030/DESC"),
]


class TestSnapshot(InTempDir):
    """
    Compare both output files in full with tests/data/expected_*.json.
    If a change to the output is intended, update those files by hand.
    """

    def test_output_files_unchanged(self):
        with mock.patch.object(bu, "get_revisions", lambda start, end, svn: SNAPSHOT_REVISIONS), \
             mock.patch.object(bu.conf, "orcid_corrections", {}):
            name, warnings = bu.write_report(args(start=1, end=26154), None)
        self.assertEqual(warnings, 1)
        for produced, expected in ((name, "expected_bulk_report.json"),
                                   ("others_" + name, "expected_others_bulk_report.json"),
                                   ("invalid_orcid_" + name, "expected_invalid_orcid_bulk_report.json")):
            with open(produced) as got, open(os.path.join(DATA, expected)) as want:
                self.assertEqual(json.load(got), json.load(want), expected)


class TestUpload(InTempDir):

    def setUp(self):
        super().setUp()
        self.reports = [report("2026-11-12T%02d:00:00Z" % hour) for hour in range(5)]
        with open("r.json", "w") as handle:
            json.dump({"resource_id": "rfam", "reports": self.reports}, handle)
        self.sent = []

    def post(self, fail_at=None):
        def fake_post(url, files, headers, data):
            batch = json.loads(files["reports"][1])["reports"]
            if len(self.calls) + 1 == fail_at:
                self.calls.append(data)
                error = requests.exceptions.HTTPError("500")
                return mock.Mock(text="error", raise_for_status=mock.Mock(side_effect=error))
            self.calls.append(data)
            self.sent.extend(batch)
            return mock.Mock(raise_for_status=lambda: None)
        self.calls = []
        return mock.patch("requests.post", fake_post)

    def test_batches_by_size(self):
        with self.post(), mock.patch.object(bu, "MAX_BATCH_BYTES", 300):
            bu.upload_bulk_report("r.json")
        self.assertEqual(self.sent, self.reports)
        self.assertGreater(len(self.calls), 1)
        self.assertEqual(self.calls, [{}] * len(self.calls))  # no delete_all

    def test_delete_all_on_first_batch_only(self):
        with self.post(), mock.patch.object(bu, "MAX_BATCH_BYTES", 300):
            bu.upload_bulk_report("r.json", delete_all=True)
        self.assertEqual(self.calls[0], {"delete_all": "rfam"})
        self.assertEqual(self.calls[1:], [{}] * (len(self.calls) - 1))

    def test_rerun_after_failure_uploads_the_rest_once(self):
        with self.post(fail_at=2), mock.patch.object(bu, "MAX_BATCH_BYTES", 300):
            with self.assertRaises(requests.exceptions.HTTPError):
                bu.upload_bulk_report("r.json", existing=[])
        self.assertTrue(0 < len(self.sent) < len(self.reports))
        first_run = list(self.sent)
        with self.post(), mock.patch.object(bu, "MAX_BATCH_BYTES", 300):
            bu.upload_bulk_report("r.json", existing=first_run)
        self.assertEqual(sorted(self.sent, key=lambda r: r["timestamp"]), self.reports)

    def test_existing_reports_skipped(self):
        with self.post():
            bu.upload_bulk_report("r.json", existing=self.reports[:3])
        self.assertEqual(self.sent, self.reports[3:])


class TestExistingReports(unittest.TestCase):

    def get(self, data, total):
        content = {"data": data, "meta": {"pagination": {"total": total}}}
        response = mock.Mock(json=lambda: content, raise_for_status=lambda: None)
        return mock.patch("requests.get", return_value=response)

    def test_reads_all(self):
        with self.get([report("2026-11-12T10:00:00.000Z")], 1) as get:
            self.assertEqual(len(bu.get_existing_reports()), 1)
        self.assertEqual(get.call_args.args[0], "https://apicuron.invalid/api/reports")

    def test_incomplete_download_stops(self):
        with self.get([report("2026-11-12T10:00:00.000Z")], 2):
            with self.assertRaisesRegex(RuntimeError, "1 of 2"):
                bu.get_existing_reports()


class TestCommandLine(unittest.TestCase):

    def parse(self, *argv):
        with mock.patch.object(sys, "argv", ["bulk_upload.py"] + list(argv)):
            return bu.parse_args()

    def test_svn_only_is_enough(self):
        parsed = self.parse("--svn", "/svn")
        self.assertEqual((parsed.start, parsed.end, parsed.rfamprod_orcid), (None, None, NANCY))

    def test_needs_svn_or_upload_file(self):
        with mock.patch("sys.stderr"), self.assertRaises(SystemExit):
            self.parse()

    def test_delete_all_needs_start(self):
        with mock.patch("sys.stderr"), self.assertRaises(SystemExit):
            self.parse("--svn", "/svn", "--delete-all")
        self.assertTrue(self.parse("--svn", "/svn", "--delete-all", "--start", "5").delete_all)

    def test_start_batch_removed(self):
        with mock.patch("sys.stderr"), self.assertRaises(SystemExit):
            self.parse("--svn", "/svn", "--start-batch", "2")


class TestMain(InTempDir):
    """Exit codes: 0 fine, 1 failure, 2 ORCIDs need a correction."""

    def run_main(self, revisions):
        root = logging.getLogger()
        old_handlers = root.handlers[:]
        root.handlers = []
        argv = ["bulk_upload.py", "--svn", "/svn", "--start", "1", "--end", "3", "--dry-run"]
        try:
            with mock.patch.object(sys, "argv", argv), \
                 mock.patch.object(bu, "get_existing_reports", lambda: []), \
                 mock.patch.object(bu, "get_revisions", lambda start, end, svn: revisions), \
                 mock.patch.object(bu.conf, "orcid_corrections", {}), \
                 mock.patch("requests.post", side_effect=AssertionError("upload attempted")), \
                 mock.patch("sys.stderr"):
                try:
                    bu.main()
                except SystemExit as exit:
                    return exit.code
                return 0
        finally:
            for handler in root.handlers:
                handler.close()
            root.handlers = old_handlers

    def test_exit_codes(self):
        good = [revision(1, "apetrov", "2026-11-12T09:00:00.1Z", "CI:x")]
        bad = [revision(1, "apetrov", "2026-11-12T09:00:00.1Z", "CI:x by " + BAD_ORCID)]
        self.assertEqual(self.run_main(good), 0)
        self.assertEqual(self.run_main(bad), 2)
        with mock.patch.object(bu, "get_family", side_effect=RuntimeError("boom")):
            self.assertEqual(self.run_main(good), 1)


class TestLogging(InTempDir):

    def test_default_log_name_has_date_and_time(self):
        logging.disable(logging.NOTSET)
        root = logging.getLogger()
        old_handlers = root.handlers[:]
        root.handlers = []
        try:
            name = bu.setup_logging(None)
            logging.info("hello")
            for handler in root.handlers:
                handler.flush()
            self.assertRegex(name, r"^apicuron_upload_\d{4}-\d{2}-\d{2}_\d{6}\.log$")
            with open(name) as handle:
                self.assertIn("INFO hello", handle.read())
        finally:
            for handler in root.handlers:
                handler.close()
            root.handlers = old_handlers


if __name__ == "__main__":
    unittest.main()
