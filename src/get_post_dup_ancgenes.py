"""
    Script to browse gene trees and extract post-duplication ancGenes.

    Example:
        $ python get_post_dup_ancgenes.py -t trees.nhx -d Clupeocephala -s sptree.nwk
                                                 [-o ancgenes_clup.tsv]
"""

import argparse
from ete3 import Tree

from scripts.synteny.duplicated_families import tag_duplicated_species
from scripts.trees.speciestree import get_species
from scripts.trees.utilities import read_multiple_objects

def write_post_dup_ancgenes(input_forest, duplicated_species, out, ancg="ancGene_TGD_"):

    """
    Browses input gene trees and writes an ancGenes file with post-duplication ancestral genes and
    their descednants. Post-duplication ancestral genes are identified as nodes under which only
    duplicated species genes are found. If this node is a duplication node, each of the two
    children subtrees represent an ancestral post-duplication genes.


    Args:
        input_forest (str): path to the input gene trees (single file)

        duplicated_species (list of str): list of the name of all duplicated species

        out (str): path to the output ancgene file

        ancg (str, optional): prefix for ancestral gene names
    """

    k = 0

    with open(input_forest, 'r') as infile, open(out, 'w') as outfile:

        for tree in read_multiple_objects(infile):

            tree = Tree(tree)

            #find all monphyletic telost groups
            tag_duplicated_species(tree.get_leaves(), duplicated_species)

            #all clades with only teleost genes in the tree
            subtrees = tree.get_monophyletic(values=["Y"], target_attr="duplicated")

            for subtree in subtrees:

                teleost_genes = sorted([i.name for i in subtree.get_leaves()])
                letter = '_A'
                if hasattr(subtree, "D"):

                    if subtree.D == 'Y':

                        for post_dup_group in subtree.children:

                            teleost_genes = sorted([i.name for i in post_dup_group])
                            outfile.write(ancg+str(k)+letter+'\t'+ ' '.join(teleost_genes)+'\n')
                            letter = '_B'

                    else:
                        outfile.write(ancg+str(k)+letter+'\t'+ ' '.join(teleost_genes)+'\n')


                else:
                    outfile.write(ancg+str(k)+letter+'\t'+ ' '.join(teleost_genes)+'\n')

                k += 1


if __name__ == '__main__':

    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-t', '--treesFile', help='Forest of trees in .nhx, with species,\
                         duplication/speciation nodes and ancestor species.',\
                         required=True)

    PARSER.add_argument('-d', '--dupSp', help='Name of the ancestor of all duplicated species.',\
                         required=True)

    PARSER.add_argument('-s', '--speciesTree', help='Species tree (newick), with ancestor names.',\
                         required=True)

    PARSER.add_argument('-o', '--outfile', help='Output file', required=False, default="out")


    ARGS = vars(PARSER.parse_args())

    #Study species
    DUP_SPECIES = get_species(ARGS["speciesTree"], ARGS["dupSp"])


    write_post_dup_ancgenes(ARGS["treesFile"], DUP_SPECIES, ARGS["outfile"])
