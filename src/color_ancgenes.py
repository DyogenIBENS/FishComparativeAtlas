"""
    Script to assign a post-duplication chromosome to ancgenes of the post-duplication ancestor,
    through a majority vote of descedant genes.

    Example:
        $ python color_ancgenes.py -ref colors_Oyzias.Latipes.txt colors_Danio.rerio.txt
                                   -g genes.Oyzias.Latipes.bz2 genes.Danio.rerio.bz2
                                   -ag ancgenes.tsv [-o out.tsv] [-f bed]
"""

import argparse

from collections import defaultdict
import operator

from scripts.synteny.mygenome import Genome

from homogenize_refs_colors import read_reference_colors


def color_ancgenes(colors, genes, ancgenes, output):

    """
    Majority vote of genes from reference species to assign a post-duplication chromosome
    to ancgenes.


    Args:

        colors (dict): for each reference duplicated species (key level 1), gives for each gene
                       (key level 2), the predicted post-duplication chromosome.

        genes (dict): for each duplicated species (key), a set with all its genes (value)

        ancgenes (str): input ancgene file

        output (str): path to the output file to write

    """

    with open(ancgenes, 'r') as infile, open(output, 'w') as outfile:

        for line in infile:

            _, descendants = line.strip().split('\t')
            descendants = descendants.split()

            votes = {}

            for sp in genes:
                sp_genes = genes[sp].intersection(descendants)
                for gene in sp_genes:
                    if gene in colors[sp]:
                        color = colors[sp][gene]

                        votes[sp] = votes.get(sp, {})
                        votes[sp][color] = votes[sp].get(color, 0) + 1

            if votes:

                for sp in votes:
                    normalization = sum(votes[sp].values())
                    votes[sp] = {k:v/normalization for k, v in votes[sp].items()}

                total = defaultdict(float)

                for vote in votes.values():
                    for k, v in vote.items():
                        total[k] += v

                winner = max(total.items(), key=operator.itemgetter(1))[0]

                ex_aequo = [i for i in total if i != winner and total[i] == total[winner]]

                if not ex_aequo:

                    outfile.write(line.strip()+'\t'+winner+'\n')

                else:
                    outfile.write(line.strip()+'\t?\n')


if __name__ == '__main__':
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-ref', '--reference_colors', type=str, nargs='+')

    PARSER.add_argument('-g', '--genes', type=str, nargs='+')

    PARSER.add_argument('-ag', '--ancestral_genes', type=str, help='post TGD ancgenes file',
                        required=True)

    PARSER.add_argument('-o', '--output', type=str, required=False, default='out.tsv')

    PARSER.add_argument('-f', '--genesformat', type=str, required=False,
                        default='bed')

    ARGS = vars(PARSER.parse_args())

    GENES = {}
    COLORS = {}

    for (ref_col_file, genes_file) in zip(ARGS["reference_colors"], ARGS["genes"]):

        COLORS[ref_col_file] = read_reference_colors(ref_col_file)

        GENES[ref_col_file] = {g.names[0] for g in Genome(genes_file, ARGS["genesformat"])}

    color_ancgenes(COLORS, GENES, ARGS["ancestral_genes"], ARGS["output"])
