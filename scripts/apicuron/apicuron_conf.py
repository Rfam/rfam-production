svn_url = ''

doc_url = 'https://docs.google.com/spreadsheets/d/{doc_id}/gviz/tq?tqx=out:csv&sheet={sheet_id}'

report_entry = {
    'activity_term': '{term}',
    'timestamp': '{timestamp}',
    'curator_orcid': '{curator}',
    'entity_uri': '{uri}'
}

svn_authors = ['apetrov', 'evan', 'ikalveri', 'jgt', 'joanna', 'nancyontiveros', 'nawrocki', 'rdf', 'rfamprod',
               'ruthe', 'sb30', 'swb', 'testuser']

curator_orcids = {
    'Emma': '0000-0003-2958-6837',
    'Nancy': '0000-0001-8457-4455',
    'Blake': '0000-0002-6497-2883'
}


sheet_terms = {
    'Done (new family)': 'create_family',
    'Done (updated family)': 'update_family'
}

rfci_terms = {
    'CI': 'update_family',
    'NEW': 'create_family',
    'KILL': 'delete_family',
    'CLCI': 'update_clan',
    'CLKILL': 'delete_clan'
}
