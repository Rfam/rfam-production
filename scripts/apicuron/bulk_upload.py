"""
Build APICURON reports from the Rfam SVN repository and upload them.

Usage:
    bulk_upload.py --svn /path/to/svn --dry-run
    bulk_upload.py --svn /path/to/svn --start 15938 --end 16099
    bulk_upload.py --upload-file bulk_report_svn_2026-11-11.json

Defaults: --end is the newest revision + 1; --start is the newest svn commit
already on APICURON. Reports APICURON already has are skipped (same activity,
curator and family; same time, or whole hours apart up to
TIMEZONE_TOLERANCE_HOURS), so a failed run can simply be repeated.
Each run logs to apicuron_upload_<date>_<time>.log unless --log-file is given.

Curator, in order of priority:
  1. ORCIDs in the commit message (required from 11 Nov 2026), e.g.
         CI:fix SS_cons by 0000-0001-8457-4455 and 0000-0002-6497-2883
  2. the svn (HPC) username, via `conf.curator_orcids`
  3. for rfamprod, --rfamprod-orcid (Nancy by default)
One report per ORCID. Commits without a curator go to the "others_" file.

An ORCID with a bad checksum is logged as a warning with a suggested fix; its
report goes to the "invalid_orcid_" file (not uploaded) and the run ends with
exit code 2. Fix it in `conf.orcid_corrections`; the next run uploads it.
"""

import argparse
import logging
import os
import requests
import datetime
import json
import re
import subprocess
import sys
import xml.etree.ElementTree as ET

from config import rfam_config
import scripts.apicuron.conf as conf


def svn_url(svn_path):
    """Return a file:// URL for a path; URLs are returned unchanged."""

    if "://" not in svn_path:
        return "file://" + os.path.abspath(svn_path)
    return svn_path


def run_svn(cmd):
    """Run svn and return stdout; raise with svn's error message on failure."""

    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise RuntimeError("%s failed: %s" % (" ".join(cmd), result.stderr.decode().strip()))
    return result.stdout


def get_youngest(svn_path):
    """Return the newest svn revision."""

    info_cmd = ["svn", "info", "--show-item", "revision", "--non-interactive",
                svn_url(svn_path)]
    return int(run_svn(info_cmd).decode().strip())


def get_revisions(start, end, svn_path):
    """
    Read revisions start..end-1 with a single `svn log`.
    :return: list of (revision, author, date, message, changed paths)
    """

    log_cmd = ["svn", "log", "--verbose", "--xml", "--non-interactive",
               "-r", "%d:%d" % (start, end - 1), svn_url(svn_path)]
    output = run_svn(log_cmd)
    return [
        (
            int(entry.get("revision")),
            (entry.findtext("author") or "").strip(),
            entry.findtext("date") or "",
            entry.findtext("msg") or "",
            " ".join(path.text or "" for path in entry.findall("paths/path")),
        )
        for entry in ET.fromstring(output).iter("logentry")
    ]


def get_family(output):
    """Return the family URI for the first RF accession in the changed paths, or ''."""

    if match := re.search(r"(RF\d{5})", output):
        return "https://rfam.org/family/" + match.group(1)
    return ""


def get_term(output):
    """Return the activity term for a commit message, e.g. 'create_family', or ''."""

    activity_term = ""
    for ci_term, apicuron_term in conf.checkin_terms.items():
        if ci_term in output:
            activity_term = apicuron_term
    return activity_term


ORCID_PATTERN = re.compile(r"(?<![\w-])(\d{4}-\d{4}-\d{4}-\d{3}[\dX])(?![\w-])", re.I)


def is_valid_orcid(orcid):
    """Check the ORCID checksum (ISO 7064 MOD 11-2)."""

    digits = orcid.replace("-", "")
    total = 0
    for digit in digits[:-1]:
        total = (total + int(digit)) * 2
    check = (12 - total % 11) % 11
    return digits[-1].upper() == ("X" if check == 10 else str(check))


def get_orcids(output):
    """Return (valid, invalid) ORCIDs in a commit message, in order, without repeats."""

    valid, invalid = [], []
    for orcid in (match.upper() for match in ORCID_PATTERN.findall(output)):
        found = valid if is_valid_orcid(orcid) else invalid
        if orcid not in found:
            found.append(orcid)
    return valid, invalid


def check_corrections(corrections):
    """Stop if `conf.orcid_corrections` itself holds a bad ORCID."""

    for rev, orcids in corrections.items():
        for orcid in orcids:
            if not (ORCID_PATTERN.fullmatch(orcid) and is_valid_orcid(orcid)):
                raise ValueError("conf.orcid_corrections r%s: invalid ORCID %s" % (rev, orcid))


def get_timestamp(output):
    """'2024-07-24T06:36:59.225407Z' -> '2024-07-24T06:36:59Z' (UTC)."""

    return output[:19] + "Z" if output else ""


def remove_existing(reports, existing):
    """
    Mark the reports APICURON already has.
    Exact times are matched first, then whole-hour shifts up to
    TIMEZONE_TOLERANCE_HOURS; each existing report matches at most once.
    :return: (list of booleans, number exact, number shifted)
    """

    def parse(timestamp):
        return datetime.datetime.strptime(timestamp[:19], "%Y-%m-%dT%H:%M:%S")

    def key(report):
        return (report["activity_term"], report["curator_orcid"], report["entity_uri"])

    pool = {}
    for report in existing:
        pool.setdefault(key(report), []).append(parse(report["timestamp"]))

    found = [False] * len(reports)
    exact = shifted = 0
    for index, report in enumerate(reports):
        times = pool.get(key(report), [])
        time = parse(report["timestamp"])
        if time in times:
            times.remove(time)
            found[index] = True
            exact += 1
    for index, report in enumerate(reports):
        if found[index]:
            continue
        times = pool.get(key(report), [])
        time = parse(report["timestamp"])
        for other in times:
            hours = (other - time).total_seconds() / 3600
            if hours == int(hours) and abs(hours) <= TIMEZONE_TOLERANCE_HOURS:
                times.remove(other)
                found[index] = True
                shifted += 1
                break
    return found, exact, shifted


def write_report(args, existing):
    """
    Write the upload file (new reports with a curator), the "others_" file
    (reports without one, with the svn username) and the "invalid_orcid_" file
    (held back until corrected).
    :param existing: reports on APICURON, or None to skip the check
    :return: (name of the upload file, number of commits with an invalid ORCID)
    """

    corrections = conf.orcid_corrections
    check_corrections(corrections)

    end = args.end or get_youngest(args.svn) + 1
    # Without --start, read everything to find where APICURON stops.
    start = args.start or 1
    logging.info("Reading svn r%d-r%d from %s", start, end - 1, args.svn)

    reports_current = []
    revisions = []
    reports_other = []
    other_revisions = []
    reports_invalid = []

    for rev, author, date, message, paths in get_revisions(start, end, args.svn):
        entry = {
            "activity_term": get_term(message),
            "timestamp": get_timestamp(date),
            "curator_orcid": author,
            "entity_uri": get_family(paths),
        }
        if any(value == "" for value in entry.values()):
            continue
        if rev in corrections:
            orcids, invalid = list(corrections[rev]), []
            if not orcids:
                continue  # [] = never upload this commit
        else:
            orcids, invalid = get_orcids(message)
        for orcid in invalid:
            reports_invalid.append(dict(entry, curator_orcid=orcid, revision=rev))
        if invalid:
            logging.warning(
                "r%d: invalid ORCID %s (failed checksum); not uploaded. Fix: add to "
                "orcid_corrections in conf.py:\n    %d: %r,  # message had %s",
                rev, ", ".join(invalid), rev, orcids + ["<correct ORCID>"] * len(invalid),
                ", ".join(invalid))
            if not orcids:
                continue  # hold back; no fallback to the svn user
        if not orcids and author in conf.curator_orcids:
            orcids = [conf.curator_orcids[author]]
        if not orcids and author == conf.shared_svn_account and args.rfamprod_orcid:
            orcids = [args.rfamprod_orcid]
        if orcids:
            for orcid in orcids:
                reports_current.append(dict(entry, curator_orcid=orcid))
                revisions.append(rev)
        else:
            reports_other.append(entry)
            other_revisions.append(rev)
    logging.info("%d reports with a curator, %d without (others_ file)",
                 len(reports_current), len(reports_other))

    if existing is not None:
        found, exact, shifted = remove_existing(reports_current, existing)
        logging.info("%d of them already on APICURON (%d exact, %d with a timezone shift)",
                     exact + shifted, exact, shifted)
        if not args.start:
            matched = [rev for rev, is_found in zip(revisions, found) if is_found]
            if not matched:
                raise RuntimeError("No svn commit found on APICURON; give --start")
            # Include that commit: a multi-ORCID commit may be half uploaded.
            start = max(matched)
            logging.info("Newest svn commit on APICURON is r%d; starting there", start)
        # Corrected commits are always included, even before the start.
        reports_current = [
            report for report, rev, is_found in zip(reports_current, revisions, found)
            if (rev >= start or rev in corrections) and not is_found
        ]
        reports_other = [
            entry for entry, rev in zip(reports_other, other_revisions) if rev >= start
        ]
        logging.info("%d new reports to upload, %d without a curator (r%d-r%d)",
                     len(reports_current), len(reports_other), start, end - 1)

    today_date = str(datetime.date.today())
    reports_file = "bulk_report_svn_" + today_date + ".json"

    with open(reports_file, "w") as bulk_report:
        reports = {"resource_id": "rfam", "reports": reports_current}
        json.dump(reports, bulk_report, indent=4, sort_keys=True)

    with open("others_" + reports_file, "w") as bulk_report:
        reports = {"resource_id": "rfam", "reports": reports_other}
        json.dump(reports, bulk_report, indent=4, sort_keys=True)

    with open("invalid_orcid_" + reports_file, "w") as bulk_report:
        reports = {"resource_id": "rfam", "reports": reports_invalid}
        json.dump(reports, bulk_report, indent=4, sort_keys=True)
    logging.info("Wrote %s, others_%s and invalid_orcid_%s",
                 reports_file, reports_file, reports_file)

    invalid_revisions = sorted({entry["revision"] for entry in reports_invalid})
    if invalid_revisions:
        logging.warning("%d commits need an ORCID correction: %s", len(invalid_revisions),
                        ", ".join("r%d" % rev for rev in invalid_revisions))
    return reports_file, len(invalid_revisions)


def get_header():
    """Return the APICURON request headers."""

    token = rfam_config.APICURON_TOKEN
    if not token:
        raise RuntimeError("No APICURON_TOKEN in config")
    token = token.strip()
    # APICURON rejects a token without the "Bearer " prefix.
    if not token.lower().startswith("bearer "):
        token = "Bearer " + token
    header = {
        "version": "2",
        "authorization": token,
    }
    return header


def get_existing_reports():
    """Download all rfam reports on APICURON, private ones included."""

    reports_url = rfam_config.BULK_REPORT_URL.rsplit("/bulk", 1)[0]
    response = requests.get(reports_url, params={"resource_id": "rfam"},
                            headers=get_header(), timeout=300)
    response.raise_for_status()
    content = response.json()
    reports = content["data"]
    total = content.get("meta", {}).get("pagination", {}).get("total", len(reports))
    if total != len(reports):
        raise RuntimeError("APICURON returned %d of %d reports" % (len(reports), total))
    logging.info("APICURON already has %d rfam reports", len(reports))
    return reports


# Larger requests get HTTP 413 (2 MB failed, 400 KB worked).
MAX_BATCH_BYTES = 300000

# Old uploads are 1 h off (UK time read as Italian time).
TIMEZONE_TOLERANCE_HOURS = 2


def upload_bulk_report(report, delete_all=False, existing=None):
    """
    Upload a report file in batches of at most MAX_BATCH_BYTES, skipping
    reports in `existing`. delete_all replaces all rfam reports on APICURON.
    """

    with open(report, "r") as reports_file:
        reports = json.load(reports_file)["reports"]
    if existing is not None:
        found, exact, shifted = remove_existing(reports, existing)
        reports = [entry for entry, is_found in zip(reports, found) if not is_found]
        logging.info("Skipping %d reports already on APICURON; %d to upload",
                     exact + shifted, len(reports))

    batches, batch, size = [], [], 0
    for entry in reports:
        entry_size = len(json.dumps(entry)) + 2
        if batch and size + entry_size > MAX_BATCH_BYTES - 100:
            batches.append(batch)
            batch, size = [], 0
        batch.append(entry)
        size += entry_size
    if batch:
        batches.append(batch)

    for number, batch in enumerate(batches, start=1):
        payload = json.dumps({"resource_id": "rfam", "reports": batch})
        try:
            response = requests.post(
                rfam_config.BULK_REPORT_URL,
                files={"reports": ("reports.json", payload, "application/json")},
                headers=get_header(),
                data={"delete_all": "rfam"} if delete_all and number == 1 else {},
            )
            response.raise_for_status()
        except requests.exceptions.HTTPError as e:
            logging.error("Batch %d/%d failed: %s. Run again to upload the rest.",
                          number, len(batches), response.text[:500])
            raise e
        logging.info("Batch %d/%d uploaded (%d reports)", number, len(batches), len(batch))


def parse_args():
    """Parse the CLI arguments."""

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--end", type=int, action="store",
        help="last revision + 1 (default: newest revision + 1)"
    )
    parser.add_argument(
        "--start",
        type=int,
        help="first revision (default: newest svn commit on APICURON)",
        action="store",
    )
    parser.add_argument(
        "--svn", type=str, help="svn repository path or URL", action="store"
    )
    parser.add_argument(
        "--rfamprod-orcid", default=conf.rfamprod_orcid,
        help="ORCID for rfamprod commits without one (default: %(default)s; '' for none)",
    )
    parser.add_argument(
        "--upload-file", help="upload this report file instead of reading svn"
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="write the reports, do not upload"
    )
    parser.add_argument(
        "--delete-all", action="store_true",
        help="delete all rfam reports on APICURON first",
    )
    parser.add_argument(
        "--log-file",
        help="log file (default: apicuron_upload_<date>_<time>.log)",
    )
    args = parser.parse_args()
    if not args.upload_file and not args.svn:
        parser.error("give --svn, or --upload-file")
    if args.delete_all and not (args.start or args.upload_file):
        parser.error("--delete-all needs --start (otherwise older reports are lost)")
    return args


def setup_logging(log_file):
    """Log to a file and the terminal; return the log file name."""

    if not log_file:
        log_file = datetime.datetime.now().strftime("apicuron_upload_%Y-%m-%d_%H%M%S.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()],
    )
    return log_file


def main():
    """Run; exit 0 if fine, 1 on failure, 2 if ORCIDs need a correction."""

    args = parse_args()
    log_file = setup_logging(args.log_file)
    logging.info("Started: %s", " ".join(sys.argv))
    try:
        # Nothing to skip when everything is replaced.
        existing = None if args.delete_all else get_existing_reports()
        warnings = 0
        if args.upload_file:
            report = args.upload_file
        else:
            report, warnings = write_report(args, existing)
        if args.dry_run:
            logging.info("Dry run: nothing uploaded")
        else:
            upload_bulk_report(report, args.delete_all, existing)
    except Exception:
        logging.exception("Failed; see %s", log_file)
        sys.exit(1)
    if warnings:
        logging.warning("Finished with %d ORCID warnings; see %s", warnings, log_file)
        sys.exit(2)
    logging.info("Finished")


if __name__ == "__main__":
    main()
