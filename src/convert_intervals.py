"""
Script to convert input data for the Paralogy Map between genome assemblies.
"""

import sys
import argparse
from collections import OrderedDict
from itertools import chain

from scripts.synteny.mygenome import Genome
from color_reference_species import load_nakatani_segments


def segments_to_genes(dgenes, dseg):

    """
    Convert genomic intervals in base pairs to list of genes in intervals.

    Args:
        dgenes (Genome): genes

        dseg (dict): genomic intervals

    Returns:
        (dict): for each genomic interval (key) the list of genes it harbors (value)
    """

    dgenes_seg = OrderedDict()

    for chromosome in dgenes.genes_list:

        for gene in dgenes.genes_list[chromosome]:

            chrom = str(chromosome)
            chrom = chrom.replace("group", "")

            assert gene.beginning < gene.end, "beg > end, check"

            for seg in dseg:

                if chrom == seg[0]:

                    interval = seg[1:]

                    if gene.end <= interval[1] and gene.beginning >= interval[0]:

                        dgenes_seg[seg] = dgenes_seg.get(seg, [])
                        dgenes_seg[seg].append(gene.names[0])

    return dgenes_seg


def load_id_conversion(input_file, genes_set):

    """
    Loads an ensembl id conversion file in dict

    Args:
        input_file (str) : path to the input file

        genes_set (set) : set of ensembl ids of the loaded genome

    Return:
        dict: conversion dict

    Note:
        This works here but is a bit tricky. Since ensembl conversion id history is sorted by
        release and I took the conversion file at the release of interest, I load here the latest
        used converted gene ID. However, it may have been retired after a certain release and as
        such it will not exist in `genes_set` --> in this case I do not load the conversion.

    """

    convert = {}
    with open(input_file, 'r') as infile:

        for line in infile:

            old_id, newest_id, _, _ = line.strip().split(', ') #oldid,newid,release,score

            if newest_id != '<retired>' and newest_id.split('.')[0] in genes_set:

                convert[old_id.split('.')[0]] = newest_id.split('.')[0]

    return convert


def convert_intervals(dconv, dgenes_seg):

    """
    Converts gene ids in interval to newest version.

    Args:
        dconv (dict): conversion dict
        dgenes_seg (dict): gene intervals to convert

    Return:
        dict: converted dict
    """

    dgenes_seg_conv = OrderedDict()

    for interval in dgenes_seg:

        for gene in dgenes_seg[interval]:

            if gene in dconv:

                dgenes_seg_conv[interval] = dgenes_seg_conv.get(interval, [])

                dgenes_seg_conv[interval].append(dconv[gene])

    return dgenes_seg_conv

def get_genes_at_limits(dgenes_seg_conv, dseg):

    """
    Finds all genes at intervals ends.

    Args:
        dgenes_seg_conv (dict): converted gene intervals
        dseg (dict): input genomic intervals

    Returns:
        OrderedDict: for each interval, a 5-element list with the two first genes, two last genes +
                     input coordinates
    """

    genes_at_limits = OrderedDict()
    for interval in dgenes_seg_conv:
        if interval in dseg:
            if len(dgenes_seg_conv[interval]) >= 4:
                genes_at_limits[interval] = [dgenes_seg_conv[interval][0],
                                             dgenes_seg_conv[interval][1],
                                             dgenes_seg_conv[interval][-2],
                                             dgenes_seg_conv[interval][-1], dseg[interval]]

    return genes_at_limits

def get_genomic_coord(gene_list, dgenes):

    """
    Extracts genomic start & stop for a list of genes.

    Args:
        gene_list (list): input gene list
        dgenes (Genome): genes

    Returns:
        OrderedDict: For each gene in `gene_list` present in `dgenes` gives its chr, start and end
    """

    genes_at_limits_coord = OrderedDict()
    for chromosome in dgenes.genes_list:

        for gene in dgenes.genes_list[chromosome]:

            if gene.names[0] in gene_list:

                genes_at_limits_coord[gene.names[0]] = (chromosome, gene.end, gene.beginning)

    return genes_at_limits_coord

def write_converted_seg(dgenes, dgenes_seg_conv, out, dseg=None):

    """
    Writes converted segment to file.

    Args:
        dgenes (Genome): genes
        dgenes_seg_conv (dict): converted gene intervals
        out (str): output file
        dseg (dict): input genomic intervals
    """

    genes_at_limits = get_genes_at_limits(dgenes_seg_conv, dseg)

    #extract genomic coordinates for these genes
    #TOIMPROVE break-case if a gene is named '1a', '1b', '2a' etc... (unlikely)
    values = set(chain.from_iterable(genes_at_limits.values()))
    genes_at_limits_coord = get_genomic_coord(values, dgenes)

    #print out conversion stats
    sys.stdout.write(f"Number of input intervals: {len(dseg)}, Converted: {len(dgenes_seg_conv)}\n")

    #write res
    with open(out, 'w') as outfile:
        prev_chrom = ''
        prev_reg = ''
        for interval in dgenes_seg_conv:

            if interval in dseg and interval in genes_at_limits:

                first_gene, second_gene, beforelast_gene, last_gene, anc = genes_at_limits[interval]

                chroms, _, start = genes_at_limits_coord[first_gene]

                chrome, end, _ = genes_at_limits_coord[last_gene]

                if str(chroms) != interval[0]:
                    chroms, _, start = genes_at_limits_coord[second_gene]

                if str(chroms) != interval[0]:
                    continue

                if str(chrome) != interval[0]:
                    chroms, _, start = genes_at_limits_coord[beforelast_gene]

                if str(chrome) != interval[0]:
                    continue

                if start > end:
                    start, end = end, start

                lg_old = interval[2] - interval[1]
                if prev_reg and (end - start) > 1.20*lg_old:
                    prev_reg = ''
                    continue

                if (end - start) > 1.20*lg_old and not prev_reg:

                    if prev_chrom and prev_chrom == chroms:
                        prev_reg = (max(int(start), int(prev_end)+1), end, lg_old, chrome, anc)
                        continue
                
                if prev_reg:
                    p_end = min(prev_reg[1], end)
                    if p_end - prev_reg[0] <= 1.20 * prev_reg[2]:
                        to_write = f"{prev_reg[3]}\t{prev_reg[0]}\t{prev_reg[1]}\t{prev_reg[4]}\n"
                        outfile.write(to_write)

                to_write = f"{chrome}\t{start}\t{end}\t{anc}\n"
                outfile.write(to_write)
                prev_reg = ''
                prev_end = end
                prev_chrom = chroms


if __name__ == '__main__':

    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-g', '--genes', nargs='+', help='Genes files, first old ids then new',
                        required=True)

    PARSER.add_argument('-seg', '--ancestral_seg', type=str, help='Nakatani et al. ancestral'
                        'chromosomes prediction file', required=True)

    PARSER.add_argument('-id', '--history_ids', type=str, help='Ensembl CONVERT_ID file',
                        required=True)

    PARSER.add_argument('-o', '--outfile', type=str, help='Name for the output file',
                        required=False, default='out')

    PARSER.add_argument('-f', '--genesformat', nargs='+', required=False, default=['bed'])

    ARGS = vars(PARSER.parse_args())


    assert len(ARGS["genes"]) == 2, "Error: two genes files should be provided, please check"

    SEG = load_nakatani_segments(ARGS["ancestral_seg"], agg=False)

    if len(ARGS["genesformat"]) == 1:
        ARGS["genesformat"].append(ARGS["genesformat"][0])

    GENOMES = []
    for in_file, in_format in zip(ARGS["genes"], ARGS["genesformat"]):
        GENOMES.append(Genome(in_file, in_format))

    GENES89_SEG = segments_to_genes(GENOMES[0], SEG)

    GENES_SET = {g.names[0] for g in GENOMES[1]}

    DCONV = load_id_conversion(ARGS["history_ids"], GENES_SET)

    CONV_SEG = convert_intervals(DCONV, GENES89_SEG)

    write_converted_seg(GENOMES[1], CONV_SEG, ARGS["outfile"], SEG)
