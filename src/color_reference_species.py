#!/usr/bin/env python

"""
    Script to segment modern species genomes with respect to post-TGD chromosomes, based on gene
    trees and pre-TGD ancestral reconstruction.

    Example:
        $ python color_reference_species.py -seg Tetraodon.seg.txt -ag TGD_ancgenes.tsv
                                            -g genes.Tetraodon.nigroviridis.bz2 [-o out] [-f dyogen]
"""

import argparse
from collections import Counter, OrderedDict
import random
from itertools import product

from mygenome import Genome


def load_segments(seg_file, agg=True):

    """
    Loads in a dictionary a segment file (as produced by Nakatani and McLysaght 2017).
    This file lists genomic intervals on a duplicated species genome and their predicted ancestral
    pre-TGD chromosome, the file should be sorted so that intervals are ordered with respect to
    their location on chromosomes.

    Args:

        seg_file (str) : path to the input file
        agg (bool, optional): aggregate in the same dict entry several contiguous segments from the
                              same predicted ancestral chromosome

    Returns:

        (dict) : for each defined genomic interval (key) its ancestral pre-dup chromosome (value)

    Raises:
        Assertion error if any overlapping or unsorted input intervals
    """

    dseg = OrderedDict()
    prev_anc = ''
    prev_chrom = ''
    prev = ''
    chrom_set = []
    with open(seg_file, 'r') as infile:

        for line in infile:

            chrom, start, stop, anc = line.strip().split('\t')

            if not prev_chrom or chrom != prev_chrom:
                chrom_set.append(chrom)

            if agg:
                if not prev:
                    prev = (chrom, (int(start), int(stop)))

                if (anc == prev_anc and chrom == prev_chrom) or not prev_anc:
                    prev = list(prev)
                    prev.append((int(start), int(stop)))
                    prev = tuple(prev)

                else:
                    dseg[prev] = prev_anc
                    prev = (chrom, (int(start), int(stop)))

            else:
                dseg[(chrom, int(start), int(stop))] = anc

            prev_anc = anc
            prev_chrom = chrom

        if agg and prev not in dseg:
            dseg[prev] = prev_anc

    assert len(chrom_set) == len(set(chrom_set)), "Intervals in input are unordered please "\
                                                  "check your input file"
    check_order_and_overlap(dseg, agg)
    return dseg


def check_order_and_overlap(dseg, agg=True):

    """
    Checks that loaded intervals are non-overlapping.

    Args:
        dseg (dict): loaded intervals (key:genomic location, value:interval label)

    Raises:
        Assertion error if any overlapping or unsorted input intervals
    """

    interval_list = []
    prev_chrom = ''
    for interval in dseg:
        chrom = interval[0]
        if not prev_chrom or chrom == prev_chrom:
            if agg:
                interval_list += list(sum(interval[1:], ()))
            else:
                interval_list += [interval[1], interval[2]]
        else:
            # print(interval_list)
            assert interval_list == sorted(interval_list), f"Intervals in input file are either "\
            "unordered or overlapping, please check your input file (Error raised for chr {chrom})"

def genes_to_segments(dgenes, dseg, transform=False):

    """
    Builds a dictionary giving correspondance between genes and their predicted ancestral
    chromosome.


    Args:

        dgenes (Genome): considered genome represented by a mygenome.Genome instance

        dseg (dict): for each genomic interval (key) its ancestral pre-TGD chromosome
                     (value)
        transform (bool, optional): whether the resulting dict should be transformed to directly
                                    give a gene to ancestral chromosome correspondance. By
                                    default, it gives a gene to genomic interval correspondance.

    Returns:

        dict: resulting dict, genes to genomic interval from the same ancestral chromosome or
              genes to ancestral chromosomes
    """

    dgenes_seg = {}

    for my_chrom in dgenes.genes_list:

        for gene in dgenes.genes_list[my_chrom]:
            chrom = str(my_chrom)
            chrom = chrom.replace("group", "")

            assert gene.beginning < gene.end, "Start > Stop, check"

            for seg in dseg:
                if chrom == seg[0]:
                    for interval in seg[1:]:
                        if gene.end <= interval[1] and gene.beginning >= interval[0]:
                            dgenes_seg[gene.names[0]] = seg
                            if transform:
                                dgenes_seg[gene.names[0]] = dseg[seg]
    print(len(dgenes_seg))
    return dgenes_seg


def load_ohnologs(ancg_file, genes):

    """
    Loads ohnologous genes as defined in the post-duplication ancgene file, where sister
    duplicated subtree are given the same ancgene ID but a different letter 'A' or 'B'.

    Args:

        ancg_file (str): path to the input ancgene file

        genes (set): all genes of the target species

    Returns:

        dict: for each gene (key) all its ohnologs (value)
    """

    d_ohno = {}
    d_anc = {}
    prev_anc = ""

    with open(ancg_file, 'r') as infile:

        for line in infile:

            anc, descendants = line.strip().split('\t')
            anc = '_'.join(anc.split('_')[:-1])
            descendants = descendants.split()

            target_descendants = genes.intersection(descendants)

            for gene in target_descendants:
                d_anc[gene] = anc

            if prev_anc and anc == prev_anc:

                for (ohno1, ohno2) in product(target_descendants, prev_target_descendants):
                    d_ohno[ohno1] = d_ohno.get(ohno1, [])
                    d_ohno[ohno2] = d_ohno.get(ohno2, [])
                    d_ohno[ohno2].append(ohno1)
                    d_ohno[ohno1].append(ohno2)

            prev_anc = anc
            prev_target_descendants = target_descendants

    return d_ohno, d_anc


def update_color(dcolor, ohnologs, seg_anc_chr, anc_chr, cutoff, summary, d_anc=None):

    """
    Updates paralogous segments.

    Args:

        dcolor (dict): filled in place, for each genomic segments its predicted post-duplication
                       chromosome, defined by paralogy with a pre-assigned segment

        ohnologs (dict): for each gene (key) all its ohnologs (value)

        seg_anc_chr (dict): for each ancestral pre-duplication chromosome (key), all sgements in
                            the modern species that descend from it (value)

        anc_chr (dict): for each genomic interval (key) its ancestral pre-dup chromosome (value)

        cutoff (float): minimum proportion of shared paralogous genes to assign paralogous segments
    """

    ohnoletter = {'a':'b', 'b':'a'}

    first_step = list(dcolor.keys())

    for sg in first_step:
        ohnologous_seg = {}
        genes = {k for k in seg_anc_chr if seg_anc_chr[k] == sg}

        for gene in genes:

            if gene in ohnologs:
                ohnos = ohnologs[gene]
                all_seg = []
                for ohno in ohnos:
                    if ohno in seg_anc_chr:
                        all_seg.append(seg_anc_chr[ohno])

                if len(set(all_seg)) == 1:
                    ohnologous_seg[all_seg[0]] = ohnologous_seg.get(all_seg[0], 0) + 1

        counts = Counter(ohnologous_seg)
        anc, letter, _ = dcolor[sg]
        letter = ohnoletter[letter]

        for seg in counts:
            genes_ohno = [k for k in seg_anc_chr if seg_anc_chr[k] == seg]
            min_genes = min(len(genes), len(genes_ohno))

            # if counts[seg]/min_genes >= cutoff and seg in dcolor and seg != sg:

                # if dcolor[seg][:2] != [anc, letter]:

                #     print([anc, letter, counts[seg]/min_genes], dcolor[seg], seg, sg)

            if counts[seg]/min_genes >= cutoff and seg not in dcolor and anc_chr[seg] == anc:

                dcolor[seg] = [anc, letter, counts[seg]/min_genes]

                if d_anc:
                    ancgenes = {d_anc[i] for i in genes_ohno if i in d_anc}
                    for ancg in ancgenes:
                        summary[ancg] = counts[seg]/min_genes





if __name__ == '__main__':


    # Arguments

    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-g', '--genes', type=str, help='Genes file', required=False, default='')

    PARSER.add_argument('-seg', '--ancestral_seg', type=str, required=False)

    PARSER.add_argument('-ag', '--ancestral_genes', type=str, help='post TGD ancgenes file',
                        required=True)

    PARSER.add_argument('-o', '--outfile', type=str, help='output file', required=False,
                        default='out.tsv')

    PARSER.add_argument('-f', '--genesformat', type=str, required=False,
                        default='bed')

    PARSER.add_argument('--random_start', action='store_true',
                        help='randomize the seed region instead of taking the longest')

    PARSER.add_argument('--write_segments', action='store_true',
                        help='also writes segments to post TGD chromosome assignments')

    ARGS = vars(PARSER.parse_args())

    ANC_CHR = load_segments(ARGS["ancestral_seg"])

    GENES = Genome(ARGS["genes"], ARGS["genesformat"])

    GENES_ANC_CHR = genes_to_segments(GENES, ANC_CHR)

    OHNOLOGS, D_ANC = load_ohnologs(ARGS["ancestral_genes"], set(GENES_ANC_CHR.keys()))

    COLORS = {}
    SUMMARY = {}
    USED = []
    for ANC in ANC_CHR.values():

        segments = {s for s in ANC_CHR if ANC_CHR[s] == ANC}

        #identify the largest block for each ancestral chromosome
        if not ARGS["random_start"]:

            max_value = 0

            for s in segments:
                l = len({k for k in GENES_ANC_CHR if GENES_ANC_CHR[k] == s})

                if l > max_value:
                    max_segment = s
                    max_value = l

            assert max_value
            COLORS[max_segment] = [ANC, 'a', 'init']
            USED.append((max_segment))

        else:
            random_letter = random.sample(["a", "b"], 1)[0]
            random_segment = random.sample(list(segments), 1)[0]
            COLORS[random_segment] = [ANC, random_letter, 'init']
            USED.append((random_segment))


    #we update paralogous segments step by step
    #we search for segments sharing paralogous genes, iteratively relaxing constraints
    for min_prop in [0.5, 0.4, 0.3, 0.2, 0.1, 0.05]:
        for i in range(10):
            update_color(COLORS, OHNOLOGS, GENES_ANC_CHR, ANC_CHR, min_prop, SUMMARY, D_ANC)

    # print('--')

    # for v in set(SUMMARY.values()):
    #     print(v, ','.join(list({i for i in SUMMARY if SUMMARY[i]==v})))
    # print(len(SUMMARY))

    GENESA = genes_to_segments(GENES, COLORS, transform=True)

    with open(ARGS["outfile"], 'w') as fw:
        for name in GENESA:
            col = ''.join(GENESA[name][:2])
            fw.write(name+'\t'+col+'\n')

    if ARGS["write_segments"]:
        SEGS = genes_to_segments(GENES, COLORS, transform=False)
        NEW_D = {}
        SEEN = []
        for my_gene in SEGS:
            my_seg = SEGS[my_gene]

            if my_seg in SEEN:
                continue

            col = ''.join(GENESA[my_gene][:2])
            SEEN.append(my_seg)
            chromosome = my_seg[0]
            for my_interval in my_seg[1:]:

                new_key = str(chromosome)+'\t'+'\t'.join([str(i) for i in my_interval])
                NEW_D[new_key] = col

        with open(ARGS["outfile"]+'_genomic', 'w') as fw:
            for my_seg in NEW_D:
                fw.write(my_seg+'\t'+NEW_D[my_seg]+'\n')
