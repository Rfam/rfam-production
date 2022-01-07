import argparse
from datetime import datetime


def create_note(version, entries):
    now = datetime.now()
    date = now.strftime("%d-%m-%Y")
    release_text = "release={version}\n" \
                   "release_date={date}\n" \
                   "entries={entries}\n".format(version=version, date=date, entries=entries)
    with open("release_note.txt", "w") as f:
        f.write(release_text)


def parse_args():
    """
    Parse the cli arguments when calling this script to insert a text file to the PDB table in the database.
    """
    parser = argparse.ArgumentParser(description='Create release note')
    parser.add_argument('--version', required=True)
    parser.add_argument('--entries', required=True)
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    create_note(args.version, args.entries)
