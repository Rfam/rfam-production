svn_url = '/nfs/production/agb/rfam/SVN_REPOSITORY/data_repos'

doc_url = "https://docs.google.com/spreadsheets/d/{doc_id}/gviz/tq?tqx=out:csv&sheet={sheet_id}"

terms_actions = {
    'Done (new family)': 'create_family',
    'Done (updated family)': 'update_family'
}

curator_orcids = {
    'Emma': '0000-0003-2958-6837'
}

report_entry = {
    'activity_term': '{term}',
    'timestamp': '{timestamp}',
    'curator_orcid': '{curator}',
    'entity_uri': '{uri}'
}

svn_authors = ['apetrov', 'evan', 'ikalveri', 'jgt', 'joanna', 'nancyontiveros', 'nawrocki', 'rdf', 'rfamprod',
               'ruthe', 'sb30', 'swb', 'testuser']

rfci_terms = {
    'CI': 'update_family',
    'NEW': 'create_family',
    'KILL': 'delete_family',
    'CLCI': 'update_clan',
    'CLKILL': 'delete_clan'
}
