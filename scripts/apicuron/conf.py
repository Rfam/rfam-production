svn_authors_current = ['nancyontiveros']

svn_authors_other = ['apetrov',
                     'evan',
                     'ikalveri',
                     'jgt',
                     'joanna',
                     'nawrocki',
                     'rdf',
                     'rfamprod',
                     'ruthe',
                     'sb30',
                     'swb',
                     'testuser']

curator_orcids = {
    'Emma': '0000-0002-7894-8112',
    'Nancy': '0000-0001-8457-4455',
    'Blake': '0000-0002-6497-2883',
    'Anton': '0000-0001-7279-2682',
    'Ioanna': '0000-0001-9424-9197',
    'Sam': '0000-0001-6043-807X',
    'nancyontiveros': '0000-0001-8457-4455',
    'bsweeney': '0000-0002-6497-2883',
    'apetrov': '0000-0001-7279-2682',
    'ecooke': '0000-0002-7894-8112',
    'nawrocki': '0000-0002-2497-3427',
    'ksong': '0009-0000-9708-9918'
}

# Shared account; its commits don't say who curated.
shared_svn_account = 'rfamprod'

# Default curator for rfamprod commits without an ORCID (--rfamprod-orcid).
rfamprod_orcid = curator_orcids['Nancy']

# ORCIDs for commits whose message has a wrong ORCID, by revision. They replace
# the message's ORCIDs; [] = never upload that commit. E.g.
#     26160: ['0000-0002-6497-2883'],  # message had 0000-0002-6497-2884
orcid_corrections = {
}

checkin_terms = {
    'CI': 'update_family',
    'NEW': 'create_family',
    'KILL': 'delete_family',
    'CLCI': 'update_clan',
    'CLKILL': 'delete_clan'
}

sheet_terms = {
    'Done (new family)': 'create_family',
    'Done (updated family)': 'update_family'
}
