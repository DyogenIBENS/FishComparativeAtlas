"""
Module docstring todo
"""

import argparse
from collections import OrderedDict
from scripts.synteny.mygenome import Genome
from itertools import chain

from color_reference_species import load_nakatani_segments


def segments_to_genes(dgenes, dseg):

    """
    Convert genomic intervals to list of genes in interval.

    Args:
        dgenes (Genome): genes

        dseg (dict): genomic intervals
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

    """

    convert = {}
    with open(input_file, 'r') as infile:

        for line in infile:

            old_id, newest_id, release, score = line.strip().split(', ')

            if newest_id != '<retired>' and newest_id.split('.')[0] in genes_set:

                convert[old_id.split('.')[0]] = newest_id.split('.')[0]

    return convert


def convert_intervals(dgenes, dconv, dgenes_seg):

    """
    Converts gene ids in interval to newest version. 

    Args:
        dgenes (Genome): genes
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


def write_converted_seg(dgenes, dgenes_seg_conv, out, dseg=None):

    """
    Writes converted segment file.

    Args:
        dgenes (Genome): genes
        dgenes_seg_conv (dict): converted gene intervals
        out (str): output file
        dseg (dict, optional): input genomic intervals, if present will be used to print stats
    """

    #find all genes at intervals ends
    genes_at_limits = OrderedDict()
    for interval in dgenes_seg_conv:
        if interval in dseg:
            genes_at_limits[interval] = [dgenes_seg_conv[interval][0], dgenes_seg_conv[interval][-1], dseg[interval]]

    #extract genomic coordinates for these genes 
    values = set(chain.from_iterable(genes_at_limits.values()))
    genes_at_limits_coord = OrderedDict()
    for chromosome in dgenes.genes_list:

        for gene in dgenes.genes_list[chromosome]:

            if gene.names[0] in values:

                genes_at_limits_coord[gene.names[0]] = (chromosome, gene.end, gene.beginning)

    #print out conversion stats
    if dseg:
        print(len(dseg), len(dgenes_seg_conv))

    #write res
    with open(out, 'w') as outfile:

        for interval in dgenes_seg_conv:

            if interval in dseg:

                first_gene, last_gene, anc = genes_at_limits[interval]

                chrom, _, start = genes_at_limits_coord[first_gene]

                _, end, _ = genes_at_limits_coord[last_gene]

                outfile.write(str(chrom)+'\t'+str(start)+'\t'+str(end)+'\t'+str(anc)+'\n')


if __name__ == '__main__':
    
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-g', '--genes',  nargs='+', help='Genes files, first old ids then new',
                        required=True)

    PARSER.add_argument('-seg', '--ancestral_seg', type=str, help='Nakatani et al. ancestral\
                        chromosomes', required=True)

    PARSER.add_argument('-id', '--history_ids', type=str, help='Ensembl CONVERT_ID file',
                        required=True)

    PARSER.add_argument('-o', '--outfile', type=str, help='output file', required=False,
                        default='out')

    PARSER.add_argument('-f', '--genesformat', nargs='+', required=False, default=['bed'])

    ARGS = vars(PARSER.parse_args())


    assert len(ARGS["genes"]) == 2, "Error: two genes files should be provided, please check"

    SEG = load_nakatani_segments(ARGS["ancestral_seg"], agg=False)

    if len(ARGS["genesformat"]) == 1:

        ARGS["genesformat"].append(ARGS["genesformat"][0])

    GENOMES = []
    for input_file, format in zip(ARGS["genes"], ARGS["genesformat"]):

        GENOMES.append(Genome(input_file, format))

    GENES89_SEG = segments_to_genes(GENOMES[0], SEG)

    GENES_SET = {g.names[0] for g in GENOMES[1]}

    DCONV = load_id_conversion(ARGS["history_ids"], GENES_SET)

    CONV_SEG = convert_intervals(GENOMES[0], DCONV, GENES89_SEG)

    write_converted_seg(GENOMES[1], CONV_SEG, ARGS["outfile"], SEG)