#!/usr/bin/python3
import sys


database_input = sys.argv[1]
genes = ["Any", '12S', '16S', 'ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'NAD1', 'NAD2', 'NAD3', 'NAD4', 'NAD4L',
         'NAD5', 'NAD6']


def database_lookup(database_file):
    database = {}
    t_start = 2
    t_end = 5
    default = [0 for x in genes]
    for line in open(database_file):
        if '>' in line:
            elements = line.rstrip().lstrip('>').split('|')
            cur_gene = elements[1]
            taxonomy = '|'.join(elements[t_start:t_end])
            database.setdefault(taxonomy,default)
            gene_index = 0
            cur_index = 0
            for g in genes:
                if cur_gene == g:
                    gene_index = cur_index
                else:
                    cur_index += 1

            database[taxonomy][gene_index] += 1




database = database_lookup(database_input)

for k, v in database:
    print (k, v)



