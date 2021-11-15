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

from mygenome import Genome

from homogenize_refs_colors import read_reference_colors


def color_ancgenes(colors, genes, ancgenes, output, propagate=True, verbose=False):

    """
    Majority vote of genes from reference species to assign a post-duplication chromosome
    to ancgenes.

    Args:

        colors (dict): for each reference duplicated species (key level 1), gives for each gene
                       (key level 2), the predicted post-duplication chromosome.

        genes (dict): for each duplicated species (key), a set with all its genes (value)

        ancgenes (str): input ancgene file

        output (str): path to the output file to write

        propagate (bool, optional): whether to propagate annotation to ancgene of second copy, when
                                    no reference species is in its descendants. For instance,
                                    if one-post TGD ancgene is annotated '12b' through reference
                                    species votes, its sister will be annotated '12a'.

    """

    letter = {"A":"B", "B":"A"}
    maj = 0
    voting = set()
    tot = set()
    with open(ancgenes, 'r') as infile, open(output, 'w') as outfile:

        store_votes = {}

        for line in infile:

            anc, descendants = line.strip().split('\t')
            descendants = descendants.split()

            votes = {}
            all_votes = set()

            for sp in genes:
                sp_genes = genes[sp].intersection(descendants)
                for gene in sp_genes:
                    if gene in colors[sp]:
                        color = colors[sp][gene]

                        votes[sp] = votes.get(sp, {})
                        votes[sp][color] = votes[sp].get(color, 0) + 1
                        all_votes.add(color)

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

                    store_votes[anc] = winner

                    outfile.write(line.strip()+'\t'+winner+'\n')
                    tot.add(anc[:-1])
                    if len(all_votes) != 1:
                        voting.add(anc[:-1])
                        maj += 1

                # else:
                #     store_votes[anc] = "?"
                #     outfile.write(line.strip()+'\t?\n')


        #propagate to un-annotated second ancgene copy, if its sister is
        if propagate:
            infile.seek(0)

            for line in infile:
                anc, descendants = line.strip().split('\t')

                if anc not in store_votes:
                    ohnologue = anc[:-1] + letter[anc[-1]]
                    if ohnologue in store_votes and store_votes[ohnologue] != "?":
                        tot.add(anc[:-1])
                        winner_ohno = store_votes[ohnologue]
                        winner = winner_ohno[:-1] + letter[winner_ohno[-1].upper()].lower()
                        outfile.write(line.strip()+'\t'+winner+'\n')
        if verbose:
            print(f'{maj} ******{len(tot)}\n')

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

    PARSER.add_argument('-dp', '--dont_propagate', action='store_false')

    PARSER.add_argument('--verbose', action='store_true')

    ARGS = vars(PARSER.parse_args())

    GENES = {}
    COLORS = {}

    for (ref_col_file, genes_file) in zip(ARGS["reference_colors"], ARGS["genes"]):

        COLORS[ref_col_file] = read_reference_colors(ref_col_file)

        GENES[ref_col_file] = {g.names[0] for g in Genome(genes_file, ARGS["genesformat"])}

    color_ancgenes(COLORS, GENES, ARGS["ancestral_genes"], ARGS["output"], ARGS["dont_propagate"],\
                   ARGS["verbose"])
