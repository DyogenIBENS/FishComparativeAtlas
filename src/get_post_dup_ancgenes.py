"""
    Script to browse gene trees and extract post-duplication ancGenes.

    Example:
        $ python get_post_dup_ancgenes.py -t trees.nhx -d Clupeocephala -s sptree.nwk
                                          [-o ancgenes_clup.tsv]
"""

import sys
import argparse
from collections import OrderedDict
from ete3 import Tree

from scripts.synteny.duplicated_families import tag_duplicated_species
from scripts.trees.speciestree import get_species, get_anc_order
from scripts.trees.orthologs import is_speciation
from scripts.trees.utilities import read_multiple_objects


def write_post_dup_ancgenes(input_forest, duplicated_species, out, outgr=None, ancg="ancGene_TGD_",
                            add_sp_names=False, root_check=None):

    """
    Browses input gene trees and writes an ancGenes file with post-duplication ancestral genes and
    their descednants. Post-duplication ancestral genes are identified as nodes under which only
    duplicated species genes are found. If this node is a duplication node, each of the two
    children subtrees represent an ancestral post-duplication genes.


    Args:
        input_forest (str): path to the input gene trees (single file)

        duplicated_species (list of str): list of the name of all duplicated species

        out (str): path to the output ancgene file

        outgr (list of str, optional): write orthologs in outgroup species

        ancg (str, optional): prefix for ancestral gene names
    """

    k = 0

    sys.stderr.write(f"Browsing trees...\n")
    sys.stderr.flush()

    with open(input_forest, 'r') as infile, open(out, 'w') as outfile:

        for nb, tree in enumerate(read_multiple_objects(infile)):

            if nb%1000 == 0 and nb:
                sys.stderr.write(f"Browsed {nb} trees...\n")
                sys.stderr.flush()

            tree = Tree(tree)


            if len(tree) == 1:
                # print(tree)
                continue

            if root_check and tree.S not in root_check:
                continue

            leaves = tree.get_leaves()

            #find all monphyletic telost groups
            tag_duplicated_species(leaves, duplicated_species)

            #all clades with only teleost genes in the tree
            subtrees = tree.get_monophyletic(values=["Y"], target_attr="duplicated")

            for subtree in subtrees:

                orthologs = OrderedDict()
                if outgr:

                    #if requested, write orthologs in outgroups
                    for sp in outgr:
                        orthologs[sp] = orthologs.get(sp, '')
                        homologs = {i for i in leaves if i.S == sp}
                        for homolog in homologs:
                            lca = tree.get_common_ancestor(subtree, homolog)
                            if is_speciation(lca):
                                orthologs[sp] += ','+homolog.name

                        orthologs[sp] = orthologs[sp][1:]

                if add_sp_names:
                    teleost_genes = sorted([i.name+'_'+i.S for i in subtree.get_leaves()])
                else:
                    teleost_genes = sorted([i.name for i in subtree.get_leaves()])

                letter = '_A'
                ortho = ''
                if orthologs:
                    ortho = '\t'+'\t'.join(orthologs.values())
                if hasattr(subtree, "D"):

                    if subtree.D == 'Y':

                        for post_dup_group in subtree.children:

                            if add_sp_names:
                                teleost_genes = sorted([i.name+'_'+i.S for i in\
                                                        post_dup_group.get_leaves()])
                            else:
                                teleost_genes = sorted([i.name for i in\
                                                        post_dup_group.get_leaves()])
                            outfile.write(ancg+str(k)+letter+'\t'+ ' '.join(teleost_genes)+ortho\
                                          +'\n')
                            letter = '_B'

                    else:
                        outfile.write(ancg+str(k)+letter+'\t'+ ' '.join(teleost_genes)+ortho+'\n')


                else:
                    outfile.write(ancg+str(k)+letter+'\t'+ ' '.join(teleost_genes)+ortho+'\n')

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

    PARSER.add_argument('-outgr', '--outgr_ortho', help='get orthologs for specified outgroup sp.',\
                        required=False, nargs='+')

    PARSER.add_argument('--add_sp', action='store_true',
                        help="Add '_' + species name to gene names")

    PARSER.add_argument('--check_root', action='store_true',
                        help="Check that tree is rooted above teleost to include it in ancGenes.")

    ARGS = vars(PARSER.parse_args())

    #Study species
    DUP_SPECIES = get_species(ARGS["speciesTree"], ARGS["dupSp"])

    ROOTS_OK = None
    if ARGS["check_root"]:
        ANC = get_anc_order(ARGS["speciesTree"], prune=False)
        ROOTS_OK = {i for i in ANC if ARGS["dupSp"] in ANC[i] or i == ARGS["dupSp"]}
        # print(ROOTS_OK)

    write_post_dup_ancgenes(ARGS["treesFile"], DUP_SPECIES, ARGS["outfile"], ARGS["outgr_ortho"],
                            add_sp_names=ARGS["add_sp"], root_check=ROOTS_OK)
