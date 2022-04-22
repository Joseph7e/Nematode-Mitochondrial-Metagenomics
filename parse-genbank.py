#!/usr/bin/python3
# Author: Joseph Sevigny
# Date: April 04, 2022
# Purpose: Parse and reformat genbank formatted files into a standardized database

from Bio import SeqIO
from Bio import Seq
import sys
import os

# INPUT FILES
if len(sys.argv) != 3:
    print('USAGE: python3 parse_genbank.py <Genbank file> <synonyms.txt>')


output_header = 'Accession,Species,Description,Taxid,Taxonomy,Feature_Type,Loc_start,Loc_end,Strand,Gene,Final_Gene,nucleotide_sequence,protein_sequence\n'
# type = curated gene name or full(for the full sequence).

gb_input = sys.argv[1]
synonyms_input = sys.argv[2]
output = open('mitochondrial-database-raw.csv', 'w')
output.writelines(output_header)
if os.path.exists('mitochondrial-database-raw.fasta'):
    print ('Database file(s) exists, delete or rename - rm mitochondrial-database-raw*')
    sys.exit()
fasta_output = open('mitochondrial-database-raw.fasta', 'w')

def parse_midori_synonyms(synonym_file):
    synonym_lookup = {}
    for line in open(synonym_file):
        feature, gene, product, count, synonym = line.rstrip().split('\t')
        synonym_lookup.setdefault(feature,{})
        for t in [gene, product]:
            if t:
                synonym_lookup[feature].setdefault(t, synonym)

    return synonym_lookup

def parse_synonyms(synonym_file):
    synonym_lookup = {}
    for line in open(synonym_file):
        gene, synonym = line.rstrip().split('\t')

        for t in ['CDS', 'rRNA']:
            synonym_lookup.setdefault(t, {})
            if t:
                synonym_lookup[t].setdefault(gene, synonym)

    return synonym_lookup


synonym_lookup = parse_synonyms(synonyms_input)
count_source = 0

# add a check to make sure at least something is going into the final database from each accession
# remove sequences that are straight 'NNNNNNNNNN' DONE
# use taxid to get taxonomy string instead of trusting gb record


for record in SeqIO.parse(gb_input, "genbank"):
    accession = record.id
    taxonomy = record.annotations['taxonomy']
    species = record.annotations['organism']
    taxonomy.append(species)
    description = record.annotations['source']
    full_sequence = record.seq
    tax_id = None

    # Identify taxonomic lineage
    for feature in record.features:
        if feature.type == "source":
            # Check for valid cross-reference qualifier
            if "db_xref" not in feature.qualifiers: continue
            for db in feature.qualifiers["db_xref"]:
                if db.startswith("taxon:"):
                    tax_id = db.split(":")[1]

    # Ensure tax_id was found
    if tax_id is None: continue
    # lineage = "|".join(generate_lineage(tax_id, taxonomy_lookup))
    if not "Metazoa" in taxonomy: continue


    for feature in record.features:
        write = True
        f_type = feature.type
        location = feature.location
        strand = feature.location.strand
        if feature.type == "CDS":
            # Check for valid gene qualifier
            if "gene" not in feature.qualifiers: continue
            gene = feature.qualifiers["gene"][0]
            synonym = synonym_lookup.get(gene, "unknown")
            acc_id = record.annotations['accessions'][0]
            strand = feature.location.strand
            start = feature.location.start.position
            end = feature.location.start.position
            sequence = feature.extract(record.seq)
            translation = feature.qualifiers['translation'][0]
            try:
                data = [acc_id, species, description, tax_id, '|'.join(taxonomy), feature.type, str(int(start)), str(int(end)), str(strand), gene, synonym, str(sequence), translation]
            #print(','.join([x.replace(',',' ') for x in data ]))
                output.writelines(','.join([x.replace(',',' ') for x in data]) + '\n')
                header = '>' + '|'.join([synonym, accession + ':' + synonym, '|'.join(taxonomy), species])
                if len(sequence) > 10: # catch short sequences from getting through
                    with open('mito-sequence-database-' + synonym + '.fna', 'a') as fna_handle:
                        fna_handle.writelines(header + '\n' + sequence + '\n')
                if translation != 'UNKNOWN' and translation != '-':
                    with open('mito-sequence-database-' + synonym + '.faa', 'a') as faa_handle:
                        faa_handle.writelines(header + '\n' + translation + '\n')
            except: #Bio.Seq.UndefinedSequenceError: Sequence content is undefined
                continue
    try:
        fasta_output.writelines(">" + accession + ' ' + '|'.join(taxonomy) + '\n' + full_sequence + '\n')
    except: # Bio.Seq.UndefinedSequenceError: Sequence content is undefined
        continue

