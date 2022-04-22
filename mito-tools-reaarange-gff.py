#!/usr/bin/python3

# example header: 10A_30000NODE_1_length_13864_cov_238744267
# working dir: /home/genome/joseph7e/HCGS/Dorota/Project_DorotaMitoV2/mito-genome-analysis/mitochondrial_genomes/mitos_data/mitos_all/gff_files
# assumes mitos annotations

import os, sys

mito_genes = ['cox1', 'nad4', 'cox3', 'cob', 'nad2', 'atp6', 'nad1', 'nad4l', 'nad6', 'nad5', 'nad3', 'cox2', 'rrnS', 'rrnL', 'atp8']

def rearrange_gff(gff, gene_id="cox1", strand="+"):
    """
    :param gff:
    :param gene_id: gene to start rearrangement on
    :return:
    """
    past_lines = []
    should_print = False
    out_name = gff + '.rearranged'
    new_start = 1
    gene_order = []
    gene_lengths = {}
    prev_node = ''

    type_counts = {'gene':[], 'tRNA':[], 'rRNA':[]}

    for g in mito_genes:
        gene_lengths[g] = 0
    with open(out_name, 'w') as out_handle:
        for line in open(gff):
            elements = line.rstrip().split('\t')
            # Grab sample name, this will need to be changed based on input data
            sample_name = elements[0].split('NODE_')[0]
            node = elements[0].replace('NODE',':NODE')
            gene_type = elements[2]
            g_length = int(elements[4]) - int(elements[3])
            g_name = elements[-1].split('=')[-1]

            cur_node = elements[0]
            if prev_node == '':
                prev_node = cur_node
            if cur_node != prev_node:
                print ('hit a new node')
                num_pcgs = str(len(type_counts['gene']))
                num_trna = str(len(type_counts['tRNA']))
                num_rrna = str(len(type_counts['rRNA']))
                num_all = ','.join([num_pcgs, num_trna, num_rrna])

                print(sample_name + ',' + node + ',' + num_all + ',' + ':'.join(gene_order) + ',' + ','.join(
                    [str(gene_lengths[x]) for x in mito_genes]))

                prev_node = cur_node
                past_lines = []
                should_print = False
                new_start = 1
                gene_order = []
                gene_lengths = {}
                prev_node = ''

                type_counts = {'gene': [], 'tRNA': [], 'rRNA': []}



            if '-' in g_name:
                g_name = g_name.split('-')[0]
            if '_' in g_name:
                g_name = g_name.split('_')[0]

            if gene_type in ['gene', 'tRNA', 'rRNA']:
                if g_name not in type_counts[gene_type]:
                    type_counts[gene_type].append(g_name)


            # add length of gene to gene_length dictionary if its pcg
            if elements[2] == 'gene' or elements[2] == 'rRNA':
                #print(elements[2], g_name, g_length)
                gene_lengths[g_name] += g_length

            # start writing based on chosen gene id
            if g_name == gene_id and not should_print:
                new_start = int(elements[3])
                should_print = True
            if should_print:
                out_handle.writelines('\t'.join(elements) + '\n')
                if elements[2] == 'gene' or elements[2] == 'rRNA':
                    gene_order.append(elements[-1].split('=')[-1])
                #     gene_lengths[g_name] += g_length
            else:
                past_lines.append(elements)

        for prev in past_lines:
            out_handle.writelines('\t'.join(prev) + '\n')
            if prev[2] == 'gene':
                gene_order.append(prev[-1].split('=')[-1])

        num_pcgs = str(len(type_counts['gene']))
        num_trna = str(len(type_counts['tRNA']))
        num_rrna = str(len(type_counts['rRNA']))
        num_all = ','.join([num_pcgs, num_trna, num_rrna])

        print(sample_name + ',' + node + ',' + num_all + ',' + ':'.join(gene_order) + ',' + ','.join([str(gene_lengths[x]) for x in mito_genes]))



        return new_start

if __name__ == '__main__':
    print ('Sample,NODE,Num_pcg,Num_tRNA,Num_rRNA,pcg_order,' + ','.join(mito_genes))
    gff_dir = sys.argv[1] # '/home/genome/joseph7e/HCGS/Dorota/Project_DorotaMitoV2/mito-genome-analysis/mitochondrial_genomes/mitos_data/mitos_all/gff_files/'

    for gff in sorted(os.listdir(gff_dir)):
        if gff.endswith('.gff'):
            rearrange_gff(gff_dir + gff)


