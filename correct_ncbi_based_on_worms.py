#!/usr/bin/python3

# Purpose: given two taxonomy databases of the same format, check and correct higher taxonomy levels based on one.
# Usage: correct_ncbi_based_on_worms.py
# Author: Joseph 7e

import sys




def fillTaxonomyLookup(database, prev_database):
    """

    :param database: tab seperated, sequence ID followed by taxonomy (higher to lower level)
    :return:
    """
    taxonomy_lookup = prev_database
    header = []
    sample_lookup = {}
    ranks = ['superkingdom', 'Kingdom', 'phylum', 'subphylum', 'superclass', 'class', 'subclass', 'superorder',
             'order', 'suborder', 'infraorder', 'superfamily',
             'family', 'subfamily', 'genus', 'species', 'subspecies']

    for line in open(database):
        elements = line.rstrip().split('\t')
        sample_lookup[elements[0]] = elements[1:]
        if line.startswith('#'):
            header = elements
            for h in header:
                taxonomy_lookup.setdefault(h, {})
        else:
            if elements[3] == 'Nematoda':
                cur_index = 1
                for t in elements[1:]:
                    if cur_index <= len(header):
                        rank = header[cur_index]
                        cur_taxa = elements[cur_index]
                        kept_taxonomy = elements[1:cur_index]
                        if 'unknown' not in cur_taxa:
                            taxonomy_lookup[rank][cur_taxa] = kept_taxonomy
                            if t not in taxonomy_lookup[rank].keys():
                                taxonomy_lookup[rank][cur_taxa] = kept_taxonomy
                            else:
                                if taxonomy_lookup[rank][cur_taxa] != kept_taxonomy:

                                    # try to fix them.
                                    my_index = 0
                                    for f in kept_taxonomy:
                                        if 'unknown_' not in f:
                                            if 'unknown_' in taxonomy_lookup[rank][cur_taxa][my_index]:
                                                taxonomy_lookup[rank][cur_taxa][my_index] = f
                                        if 'unknown_' not in taxonomy_lookup[rank][cur_taxa][my_index]:
                                            if 'unknown_' in f:
                                                kept_taxonomy[my_index] = taxonomy_lookup[rank][cur_taxa][my_index]
                                        my_index += 1
                                    if taxonomy_lookup[rank][cur_taxa] != kept_taxonomy:
                                        # DO NOTHING STUPID!
                                        taxonomy_lookup[rank][cur_taxa] = kept_taxonomy
                                        print (cur_taxa, rank, "UNTRUSTED!!!!!!!!!!", "Keeping the new one.")
                                        print('#', taxonomy_lookup[rank][cur_taxa])
                                        print ('#', kept_taxonomy)
                    cur_index += 1
    return taxonomy_lookup, sample_lookup

if __name__ == '__main__':
    ncbi_database = sys.argv[1]
    worms_database = sys.argv[2]
    worms_lookup, worms_sample_lookup = fillTaxonomyLookup(worms_database, {})
    print ("Parsed worms database, correcting NCBI now")
    #worms_lookup, worms_sample_lookup = fillTaxonomyLookup(worms_database, {})
    with open('expanded_ncbi_taxonomy_worms_corrected.tsv', 'w') as outhandle:
        ranks = []
        for line in open(ncbi_database):
            elements = line.rstrip().split('\t')
            if line.startswith('#'):
                outhandle.writelines(line)
                ranks = elements[1:]
            # only correct family and up!
            else:
                seq_id = elements[0]
                original_taxonomy = elements[1:]
                cur_taxonomy = elements[1:]
                if 'Nematoda' in elements: # only fix nematodes
                    cur_index = len(ranks) - 1
                    found = False
                    while cur_index >= 2 or found == False: # start at genus and move up, stop when you find a hit
                        if cur_index < 15: # genus or more
                            cur_rank = ranks[cur_index]
                            cur_taxa = cur_taxonomy[cur_index]
                            # check if taxa exists in worms database
                            if cur_taxa in worms_lookup[cur_rank].keys():
                                fix_index = 0
                                for f in worms_lookup[cur_rank][cur_taxa]:
                                    cur_taxonomy[fix_index] = f
                                    fix_index += 1
                                found = True

                        cur_index -= 1
                outhandle.writelines(elements[0] + '\t' + '\t'.join(cur_taxonomy) + '\n')
