# Nematode-Mitochondrial-Metagenomics-
Scripts and data for the analyses in the manuscript - Nematode Mitochondrial Metagenomics – a New Tool for Biodiversity Analysis


# Collation of publicly available data

The genomes and gene sequences downloaded from refseq and genbank were processed and reformatted using the script parse_genbank.py, which parses the genbank formatted file and writes it to a standardized format.

## NCBI refseq, complete mitochondrial genomes
mitochondrial genomes are avaialble from NCBI's ftp site.
```bash
## download genbank files
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.genomic.gbff.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.2.genomic.gbff.gz
gunzip *.gbff.gz

python3 parse_genbank.py mitochondrion.1.genomic.gbff synonyms.tsv
python3 parse_genbank.py mitochondrion.2.genomic.gbff synonyms.tsv
```

Note: The master code for parsing genbank files and creating database format is available in the NGMMAT repo. Please utilize that for research purposes.
Note: If running the parse_genbank commmands in the same directory you need to move previous versions

## NCBI genbank, partial  genomes and individual genes sequences
Search the NCBI nucleotide database with https://www.ncbi.nlm.nih.gov/nuccore/?term=(%22Nematoda%22%5BOrganism%5D+AND+mitochondrion%5Bfilter%5D+) and downlaod the full set in genbank format, optional - set a length requirement.

```bash
python3 parse_genbank.py ncbi-download.gbff synonyms.tsv
````


# synonyms.txt
The file "synonyms.txt" is a mapping file that standardizes the nomenclature for the 13 protein coding genes, 22 tRNA, and 2 ribosomal subunits found within most animal mitochondrial genomes. This file may be updated to reflect newly updated sequences deposited in genbank that use different nomenclature. 


