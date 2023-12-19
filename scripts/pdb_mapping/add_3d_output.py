import datetime

from pdb_config import add_3d_git_output, pdb_files


def write_3d_output():
    """
    Extract relevant info from the rfam-3d-seed-alignments script
    and write to a file to send to the Slack channel
    """
    changes = 0
    no_changes = "\n No updates to be made after running the rfam-3d-seed-alignments script.\n"
    with open(add_3d_git_output, "r") as output:
        contents = output.readlines()
        for line in contents:
            if "data/output" in line:
                changes += 1
                line = line.replace("data/output/", "")
                line = line.replace(".sto", "")
                if "A" in line:
                    added = "The following families have newly added 3D information: \n"
                    line = line.replace("A", "")
                    added += line
                elif "M" in line:
                    modified = "\nThe following families have been updated with 3D information: \n"
                    line = line.replace("M", "")
                    modified += line

    today_date = str(datetime.date.today())
    pdb_txt = "{dir}/pdb_families_{date}.txt".format(dir=pdb_files, date=today_date)

    with open(pdb_txt, "a") as pdb_file:
        if changes == 0:
            pdb_file.write(no_changes + "\n")
        else:
            if modified:
                pdb_file.write(modified + "\n")
            if added:
                pdb_file.write(added + "\n")


if __name__ == '__main__':
    write_3d_output()
