import os
import tempfile


def get_all_rfam_accessions():
    """
    Fetch a list of all Rfam families from the SVN repository.
    """
    rfam_accessions = []
    svn_url = 'https://xfamsvn.ebi.ac.uk/svn/data_repos/trunk/Families/'
    svn_list = tempfile.NamedTemporaryFile()
    cmd = "svn list {} > {}".format(svn_url, svn_list.name)
    os.system(cmd)
    with open(svn_list.name, 'r') as f_svn_list:
        for line in f_svn_list:
            if line.startswith('RF'):
                rfam_accessions.append(line.strip().replace('/', ''))
    print('Found {} accessions on SVN'.format(len(rfam_accessions)))
    return rfam_accessions
