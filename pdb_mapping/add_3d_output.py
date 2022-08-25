import datetime

from pdb_config import add_3d_git_output, pdb_files


def write_3d_output():
    """
    Extract relevant info from the rfam-3d-seed-alignments script
    and write to a file to send to the Slack channel
    """
    modified = "\nThe following families have been updated with 3D information: \n"
    added = "The following families have been newly added with 3D information: \n"
    with open(add_3d_git_output, "r") as output:
        contents = output.readlines()
        for line in contents:
            if "data/output" in line:
                line = line.replace("data/output/", "")
                line = line.replace(".sto", "")
                if "A" in line:
                    line = line.replace("A", "")
                    added += line
                elif "M" in line:
                    line = line.replace("M", "")
                    modified += line
            else:
                msg = "No updates to be made after running the rfam-3d-seed-alignments script.\n"

    today_date = str(datetime.date.today())
    pdb_txt = "{dir}/pdb_families_{date}.txt".format(dir=pdb_files, date=today_date)

    with open(pdb_txt, "a") as pdb_file:
        if msg:
            pdb_file.write(msg + "\n")
        else:
            pdb_file.write(modified + "\n")
            pdb_file.write(added + "\n")


if __name__ == '__main__':
    write_3d_output()
