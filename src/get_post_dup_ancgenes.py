#!/usr/bin/env python

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

from speciestree import get_species, get_anc_order


def is_speciation(node):

    """
    Is the node a speciation node?

    Args:
        tree (ete3.TreeNode): input node, with duplications annotated with the `D` attribute.
                              D=Y if duplication, D=N otherwise. Note that dubious nodes
                              (DD=Y or DCS=0) are considered speciation nodes.

    Returns:
        bool: True if speciation, False otherwise.

    """

    speciation = False

    if (hasattr(node, "D") and node.D == 'N'):
        speciation = True

    elif (hasattr(node, "DD") and node.DD == 'Y'):
        speciation = True

    elif (hasattr(node, "DCS") and float(node.DCS) == 0.0):
        speciation = True

    return speciation


def tag_duplicated_species(leaves, duplicated):

    """
    Adds a tag to genes of duplicated species in an ete3.Tree instance, in-place.

    Args:
        leaves (list of ete3.TreeNode): leaves of the tree
        duplicated (list of str): list of the names of all duplicated species
    """

    for leaf in leaves:

        if leaf.S in duplicated:
            leaf.add_features(duplicated='Y')

        else:
            leaf.add_features(duplicated='N')


def read_multiple_objects(file_object, sep="//"):

    """
    Creates a generator to read a file with several entries (trees, alignments, or other...) one
    by one.

    Args:
        file_object (file): python file object of the input file
        sep (str, optional): the separator between entries

    Yields:
        str : the next tree (or alignment).
    """

    data = ""
    while True:

        line = file_object.readline()

        #don't forget to yield even if no separator at the end of file
        if not line:
            if data:
                yield data
            break

        #yield stored data each time we see the separator
        if line.strip() == sep:
            yield data
            data = ""

        else:
            data += line


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

    PARSER.add_argument('--add_sp', action='store_true', help="Add '_'+species to gene names")

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
