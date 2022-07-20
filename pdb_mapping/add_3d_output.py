import datetime

from pdb_config import add_3d_output


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
                summary.join(line + '\n')

    today_date = str(datetime.date.today())
    pdb_txt = "pdb_families_{0}.txt".format(today_date)

    with open(pdb_txt, "a") as pdb_file:
        pdb_file.write("Added 3D information into the seed of below families: \n")
        pdb_file.write(summary + "\n")


if __name__ == '__main__':
    write_3d_output()
