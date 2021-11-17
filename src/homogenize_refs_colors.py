#!/usr/bin/env python

"""
    Script to homogenize post-duplication chromsomes that have been assigned independently on
    reference species, so that orthologous chromsomes have the same name.

    Example:

        $ python homogeneize_refs_colors.py -i colors_Oyzias.Latipes.txt colors_Danio.rerio.txt
                                            -g genes.Oyzias.Latipes.bz2 genes.Danio.rerio.bz2
                                            -ag ancgenes.tsv -guide_sp Oyzias.Latipes [-f bed]
                                            [-os homogenized]

"""

import argparse

from mygenome import Genome


def read_reference_colors(file_ref_colors):

    """
    Loads a tab-delimited file with predicted post-duplication ancestral chromosomes for one ref
    duplicated species


    Args:
        file_ref_colors (str): path to input file


    Returns:
        dict: Correspondance between genes (key) and their post-duplication chromosomes (value)
    """

    col = {}
    with open(file_ref_colors, 'r') as infile:
        col = {line.split('\t')[0]:line.strip().split('\t')[1] for line in infile}
    return col


def homogenize_on_guide(orthologs, d_color, guide, k=13):

    """
    Homogenizes post-deuplication chromsomes names across species, so that orthologs have the same
    name, using one of the references as guide.


    Args:
        orthologs (dict): For each ref species (key level 1), correspondance between its genes
                          (key level 2) and corresponding orthologs in the guide (value)

        d_colors (dict): For each ref species (key level 1), correspondance between genes
                         (key level 2) and their predicted post-duplication chromosomes (value)

        guide (str): Name of the species to use as guide. Orthologs to the guide species will be
                     named according to names in the guide.


    Returns:
        dict: For each ref species (key level 1, excluding the guide), a correspondance between
              post-duplication chromosomes names before (key level 2) and after homogenization
              (value)

    """

    homogenized = {}

    for chrom in range(1, k+1):

        genesa_ref = {i for i in d_color[guide] if d_color[guide][i] == str(chrom)+'a'}
        genesb_ref = {i for i in d_color[guide] if d_color[guide][i] == str(chrom)+'b'}

        for species in d_color:

            homogenized[species] = homogenized.get(species, {})
            if species != guide:
                genesa = {i for i in d_color[species] if d_color[species][i] == str(chrom)+'a'}
                genesb = {i for i in d_color[species] if d_color[species][i] == str(chrom)+'b'}

                genesa_ortho = set()
                for gene in genesa:
                    if gene in orthologs[species]:
                        genesa_ortho.update(set(orthologs[species][gene]))

                genesb_ortho = set()
                for gene in genesb:
                    if gene in orthologs[species]:
                        genesb_ortho.update(set(orthologs[species][gene]))

                if len(genesa_ortho.intersection(genesb_ref)) +\
                   len(genesb_ortho.intersection(genesa_ref)) >\
                   len(genesb_ortho.intersection(genesb_ref)) +\
                   len(genesa_ortho.intersection(genesa_ref)):

                    homogenized[species][str(chrom)+'a'] = str(chrom)+'b'
                    homogenized[species][str(chrom)+'b'] = str(chrom)+'a'

                else:

                    homogenized[species][str(chrom)+'a'] = str(chrom)+'a'
                    homogenized[species][str(chrom)+'b'] = str(chrom)+'b'

    return homogenized


def orthologs_with_guide(ancgenes, genes_target, genes_guide, ortho):

    """
    Loads orthologs of genes of the guide species, with orthologs as defined by ancestral genes.
    Fills the `ortho` dict in-place.


    Args:
        ancgenes (str): path to the input ancgene file

        genes_target (set): set of all genes of a target reference species

        genes_guide (set): set of all genes of the guide species

        ortho (dict): dict to fill in-place and store orthologs in: key = gene of the target,
                      value = list of ortholog(s) in the guide

    """
    with open(ancgenes, 'r') as infile:

        for line in infile:

            _, descendants = line.strip().split('\t')
            descendants = descendants.split()

            target_descendants = genes_target.intersection(descendants)
            guide_descendants = genes_guide.intersection(descendants)

            for ortho1 in target_descendants:
                ortho[ortho1] = ortho.get(ortho1, [])
                ortho[ortho1] += list(guide_descendants)


def transform_ref_colors(input_file, d_colors, output_file):

    """
    Writes a new post-duplication chromsome prediction file, with now consistent names for
    orthologs.


    Args:
        input_file (str): path to tab-deliminted input file with predicted post-duplication
                          ancestral chromosomes for one ref species

        d_colors (dict): Correspondance between post-duplication chromosomes names before (key)
                         and after homogenization (value)

        output_file (str): path to tab-deliminted output file with predicted post-duplication
                          ancestral chromosomes for one ref species, after homogenization

    """

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            gene, anc = line.strip().split('\t')
            if d_colors:
                anc = d_colors[anc]
            outfile.write(gene+'\t'+anc+'\n')


if __name__ == '__main__':

    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-i', '--reference_colors', type=str, nargs='+')

    PARSER.add_argument('-g', '--genes', type=str, nargs='+')

    PARSER.add_argument('-ag', '--ancestral_genes', type=str, help='post TGD ancgenes file',
                        required=True)

    PARSER.add_argument('-guide_sp', '--guide_species', type=str)

    PARSER.add_argument('-os', '--output_suffix', type=str, default='homogenized', required=False)

    PARSER.add_argument('-f', '--genesformat', type=str, required=False,
                        default='bed')

    ARGS = vars(PARSER.parse_args())

    GENES = {}
    COLORS = {}
    ORTHO = {}

    for (ref_col_file, genes_file) in zip(ARGS["reference_colors"], ARGS["genes"]):

        COLORS[ref_col_file] = read_reference_colors(ref_col_file)

        GENES[ref_col_file] = {g.names[0] for g in Genome(genes_file, ARGS["genesformat"])}

    for ref_col_file in ARGS["reference_colors"]:

        if ARGS["guide_species"] in ref_col_file:
            GUIDE = ref_col_file

    for ref_col_file in ARGS["reference_colors"]:

        if ARGS["guide_species"] not in ref_col_file:

            ORTHO[ref_col_file] = {}

            orthologs_with_guide(ARGS["ancestral_genes"], GENES[ref_col_file],
                                 GENES[GUIDE], ORTHO[ref_col_file])


    HOMOGENIZED = homogenize_on_guide(ORTHO, COLORS, GUIDE)

    for ref_col_file in COLORS:

        out = ".".join(ref_col_file.split(".")[:-1]) + "_"+ ARGS['output_suffix'] + ".txt"

        transform_ref_colors(ref_col_file, HOMOGENIZED.get(ref_col_file, None), out)
