#!/usr/bin/env python

import json

#TODO do not replace valid accessions
#TODO generate new SEED alignment
#TODO create tests using the H742DVSC01R-Alignment.json file

def get_accession(gid):
    """
    Get versioned accession, for example:
    Input:
    gi|2047803076|gb|CP061286.1|
    Output:
    CP061286.1
    """
    parts = gid.split('|')
    return parts[-2]


def get_data(filename):
    """
    Load BLAST JSON output.
    """
    with open(filename, 'r') as f:
        return json.load(f)


def choose_replacement(data):
    """
    Loop over BLAST results and pick best replacement for each hit.
    """
    # do not pick replacement from the same accession if already seen
    seen_accessions = set()
    for query_num, search in enumerate(data['BlastOutput2']):
        query_title = search['report']['results']['search']['query_title']
        query_len = search['report']['results']['search']['query_len']
        replacement_found = False
        for entry in search['report']['results']['search']['hits']:
            acc = get_accession(entry['description'][0]['id'])
            if acc not in seen_accessions:
                seen_accessions.add(acc)
                replacement_found = True
            else:
                continue
            sequence = entry['hsps'][0]['hseq']
            start = entry['hsps'][0]['hit_from']
            end = entry['hsps'][0]['hit_to']
            align_len = entry['hsps'][0]['align_len']
            gaps = entry['hsps'][0]['gaps']
            exact_matches = entry['hsps'][0]['identity']
            identity = float(exact_matches) / align_len * 100
            query_coverage = float(align_len - gaps) / query_len * 100
            target_coverage = float(align_len - gaps) / len(sequence) * 100
            fasta = '>{acc}/{start}-{end}\n{sequence}'.format(acc=acc, start=start, end=end, sequence=sequence.replace('-', '').replace('T', 'U'))
            summary = ('#{query_num} {warning}Replace {query_title} '
                       'with {acc}/{start}-{end} at {identity}% identity; '
                       '{gaps} gaps; query coverage {query_coverage}').format(
                acc=acc, start=start, end=end, query_title=query_title,
                identity=round(identity), query_coverage=round(query_coverage),
                target_coverage=round(target_coverage, 2), gaps=gaps,
                warning='' if (identity >= 90 and query_coverage >=70) else '      WARNING: ',
                query_num=query_num+1
            )
            print summary
            if replacement_found:
                break

def main():
    filename = 'H742DVSC01R-Alignment.json'
    data = get_data(filename)
    choose_replacement(data)


if __name__ == '__main__':
    main()
