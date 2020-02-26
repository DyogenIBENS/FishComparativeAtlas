"""
    Script to draw a modern duplicated genome after assigning its genes to a post-duplication
    chromosome.

    Example:
        $ python plot_paralogy_map.py -c colored_ancgenes.tsv -g genes.Oryzias.latipes.list.bz2
                                      -o Medaka_paralogymap.svg [-s Medaka] [-f bed]
"""

import argparse

from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection

import seaborn as sns

from order_chrom import ORDER_CHROM

from scripts.synteny.mygenome import Genome


def read_ancgenes_colors(file_anc_colors, genes, anc=False, species=''):

    """
    Loads predicted colors (PREDUP_CHROM+'_'+POSTDUP_LETTER) for all genes of the considered
    modern species. If `anc` is True load only ancestral genes.


    Args:

        file_anc_colors (str): path to the tab-delimited annotated ancgenes file

        genes (set): set with names of all gene names of the considered species

        anc (bool, optional): whether only ancestral genes should be loaded
                              (if True `genes` can can be empty)


    Returns:

        dict: gives for each modern gene (key) with prediction its predicted post-duplication
              ancestral chromosome (value).
    """

    assert genes or anc, "empty genes list please check arguments"

    genes_anc = {}
    tot = len(genes)
    pred = 0

    with open(file_anc_colors, 'r') as infile:

        for line in infile:

            ancgene, descendants, anc_chr = line.strip().split('\t')
            descendants = descendants.split()

            if anc:
                genes_anc[ancgene] = anc_chr
                continue

            if anc_chr != '?':

                target_genes = genes.intersection(descendants)

                for gene in target_genes:
                    genes_anc[gene] = anc_chr
                    pred += 1

    frac = pred/float(tot) * 100

    print(species+': '+str(pred)+' annotated genes ('+str(round(frac, 2))+'%)\n')

    return genes_anc


def draw_colors(dgenes, order, genes_colors, species, out, palette, min_length=30, max_chr=30):

    """
    Uses matplotlib to draw the genome annotated by post-duplication chromosomes.


    Args:
        dgens (Genome.genes_list): genome to plot

        order (dict): pre-assigned chromosomes order based on Figures in Nakatani and McLysaght

        genes_colors (dict): for each gene (key) its predicted post-duplication chromosomes (value)

        out (str): path for output figure

        palette (dict): for each post-duplication chromosomes (key) its associated color (value)

        min_length (int, optional): minimum number of genes to plot a chromosome

        max_chr (int, optional): maximum numbre of chromosomes to plot
    """

    if species in order:
        order = order[species]

    else:
        order = [str(i) for i in sorted(dgenes.keys(), key=lambda chrom: len(dgenes[chrom]),\
                reverse=True)]
    i = 0
    chrom_draw = {}
    height = 0.9
    spacing = 0.9
    xranges, colors = [], []

    for chromosome in dgenes:
        for gene in dgenes[chromosome]:
            chrom = str(chromosome)
            chrom = chrom.replace("group", "")
            name = gene.names[0]

            col = 'whitesmoke'
            if name in genes_colors:

                col = ''.join(genes_colors[name])

                if col in palette:
                    col = palette[col]

            xranges.append((i, 1))
            colors.append(col)
            i += 1

        chrom_draw[chrom] = (xranges, colors)
        xranges, colors = [], []
        i = 0

    _, ax = plt.subplots(1, 1)
    plt.title(species +" Paralogy Map")
    yticks = []
    yticklabels = []
    ymin, nb = 0, 0
    for chrom in order:
        xranges, colors = chrom_draw[chrom]
        if len(colors) > min_length and nb < max_chr:
            ymin += height + spacing
            yrange = (ymin, height)
            coll = BrokenBarHCollection(xranges, yrange, facecolors=colors)
            ax.add_collection(coll)
            center = yrange[0] + yrange[1]/2.
            yticks.append(center)
            yticklabels.append(chrom)
            nb += 1

    ax.axis('tight')
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)
    plt.ylabel(species+ ' chromosomes')
    plt.xlabel("Genomic position (in genes unit)")
    sns.despine()
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close('all')


if __name__ == '__main__':
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-c', '--color', type=str, help='Input file', required=True)

    PARSER.add_argument('-g', '--genes', type=str, help='Genes file', required=True, default='')

    PARSER.add_argument('-o', '--output_file', type=str, help='output file', required=False,
                        default='out.svg')

    PARSER.add_argument('-s', '--species_name', type=str, required=False, default='')

    PARSER.add_argument('-f', '--genesformat', type=str, required=False,
                        default='bed')

    ARGS = vars(PARSER.parse_args())

    PALETTE = {'5a': "lime", '5b': "greenyellow", "1a": "red", "1b":"crimson", "9b":"darkorange",\
               "9a": "orangered", "13a":'#CD00CD', '13b': "#FF32FF", "3a":"darkblue",\
               "3b":"royalblue", "2a":"black", "2b":"#404040", "8a": "darkgreen",\
               "8b": "mediumseagreen", "6a":"#707070", "6b":'#B8B8B8', '4a': "brown", "4b": "peru",\
               "10a": "gold", "10b":"yellow", "11a":"mediumorchid", "11b": "plum",\
               "12a":"deeppink", "12b":"hotpink", "7a": "deepskyblue", "7b": "lightskyblue"}

    GENES = {}

    GENOME = Genome(ARGS["genes"], ARGS["genesformat"])
    GENES = {g.names[0] for g in GENOME}

    GENES_COL = read_ancgenes_colors(ARGS["color"], GENES, species=ARGS["species_name"])

    draw_colors(GENOME.genes_list, ORDER_CHROM, GENES_COL, ARGS["species_name"],\
                ARGS["output_file"], PALETTE)