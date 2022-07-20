import datetime

from pdb_config import add_3d_output, pdb_files


def write_3d_output():
    """
    Extract relevant info from the rfam-3d-seed-alignments script
    and write to a file to send to the Slack channel
    """
    summary = ""
    with open(add_3d_output, "r") as output:
        contents = output.readlines()
        for line in contents:
            if "Created data/output" in line:
                summary += line

    today_date = str(datetime.date.today())
    pdb_txt = "{dir}/pdb_families_{date}.txt".format(dir=pdb_files, date=today_date)

    with open(pdb_txt, "a") as pdb_file:
        pdb_file.write("\nAdded 3D information into the seed of below families: \n")
        pdb_file.write(summary + "\n")
        pdb_file.write("To see the full output of the rfam-3d-add-seed-alignments script please see: {file}"
                       .format(file=add_3d_output))


if __name__ == '__main__':
    write_3d_output()
