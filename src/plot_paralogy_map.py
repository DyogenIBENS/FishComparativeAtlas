"""
    Script to draw a modern duplicated genome after assigning its genes to a post-duplication
    chromosome.

    Example:
        $ python plot_paralogy_map.py -c colored_ancgenes.tsv -g genes.Oryzias.latipes.list.bz2
                                      -o Medaka_paralogymap.svg [-s Medaka] [-f bed] [-maxC 30]
                                      [-minL 30]
"""

# python src/plot_paralogy_map.py -c ../SCORPiOs/ancgenes_inconsistent_trees.tsv
# -g Salmo_salar_genes_chrom_names.bz2 -o test_salmon.svg -s Salmo.salar -sort "names" -f dyogen
# -t "4R sequence/synteny conflicts --save

import argparse
import pickle


from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
from matplotlib.colors import is_color_like

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
                              (if True `genes` can be empty)

        species (str, optional): species name, optional, used to print annotation stats


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

            line = line.strip().split("\t")

            ancgene, descendants = line[:2]
            anc_chr = line[-1]
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


def draw_colors(dgenes, order, genes_colors, species, out, palette, min_length=30, max_chr=30,\
                sort_by="size", title='Paralogy Map', save=False):

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

        sort_by (str, optional): 'size' if chromosomes are to be sorted by size, 'name' by names.

        title (str, optional): title for the generated figure

        save (bool, optional): whether to pickle dump the plotted python dict
    """

    default_palette = sns.color_palette("Set2")
    assert sort_by in ["size", "names"], "Invalid `sort_by` argument, please check"

    #loaded pre-defined chrom order, if specified
    if species in order:
        order = order[species]

    #compute chrom order by size
    elif sort_by == "size":
        order = [str(i) for i in sorted(dgenes.keys(), key=lambda chrom: len(dgenes[chrom]),\
                reverse=True)]

    #compute chrom order by name
    elif sort_by == "names":
        order = sorted([i for i in dgenes.keys() if isinstance(i, int)])
        order = [str(i) for i in order]

    i = 0
    chrom_draw = {}
    height = 0.9
    spacing = 0.9
    xranges, colors = [], []
    j = 0

    #load pre-defined color palette
    for j, color in enumerate(default_palette):
        palette[str(j)] = color

    #Fill the python dict for plot
    for chromosome in dgenes:
        for gene in dgenes[chromosome]:
            chrom = str(chromosome)
            chrom = chrom.replace("group", "")
            name = gene.names[0]

            col = 'whitesmoke'
            if name in genes_colors:

                col = ''.join(genes_colors[name])

                assert col in palette or is_color_like(col), f'Cannot understand color {col}'

                if col in palette:

                    col = palette[col]

            xranges.append((i, 1))
            colors.append(col)
            i += 1

        chrom_draw[chrom] = (xranges, colors)
        xranges, colors = [], []
        i = 0

    #plot
    _, ax = plt.subplots(1, 1)
    plt.title(species +" "+title)
    yticks = []
    yticklabels = []
    ymin, nb = 0, 0
    to_save = {}
    for chrom in order:
        xranges, colors = chrom_draw[chrom]
        to_save[chrom] = colors
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

    #dump dict to file
    if save:
        with open(species+"_"+title+'.pkl', "wb") as outf:
            pickle.dump(to_save, outf)

def load_palette_from_file(input_file):

    """
    Loads a color palette from file.

    Args:
        input_file (str): Input file with chromosome name and associated color

    Returns:
        dict: for each chromosome (key) ist associated color (value)
    """

    palette = {}
    with open(input_file, 'r') as infile:
        for line in infile:
            chrom, col = line.strip().split("\t")
            palette[chrom] = col
    return palette


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

    ## plot-related options
    PARSER.add_argument('-maxC', '--max_chr', type=int, required=False,
                        default=30)

    PARSER.add_argument('-minL', '--min_length', type=int, required=False,
                        default=30)

    PARSER.add_argument('-sort', '--sort_by', type=str, required=False,
                        default="size")

    PARSER.add_argument('-t', '--title', type=str, required=False,
                        default="Paralogy Map")

    PARSER.add_argument('-pf', '--palette_from_file', type=str, required=False, default='')

    PARSER.add_argument('--save', action='store_true')

    ARGS = vars(PARSER.parse_args())

    if not ARGS["palette_from_file"]:


        PALETTE = {'5a': "lime", '5b': "greenyellow", "1a": "red", "1b":"crimson",\
                   "9b":"darkorange", "9a": "orangered", "13a":'#CD00CD', '13b': "#FF32FF",\
                   "3a":"darkblue", "3b":"royalblue", "2a":"black", "2b":"#404040",\
                   "8a": "darkgreen", "8b": "mediumseagreen", "6a":"#707070", "6b":'#B8B8B8',\
                   '4a': "brown", "4b": "peru", "10a": "gold", "10b":"yellow", "11a":"mediumorchid",
                   "11b": "plum", "12a":"deeppink", "12b":"hotpink", "7a": "deepskyblue",\
                   "7b": "lightskyblue"}
    else:

        PALETTE = load_palette_from_file(ARGS["palette_from_file"])

    GENES = {}

    GENOME = Genome(ARGS["genes"], ARGS["genesformat"])
    GENES = {g.names[0] for g in GENOME}

    GENES_COL = read_ancgenes_colors(ARGS["color"], GENES, species=ARGS["species_name"])

    draw_colors(GENOME.genes_list, ORDER_CHROM, GENES_COL, ARGS["species_name"],
                ARGS["output_file"], PALETTE, min_length=ARGS["min_length"],
                max_chr=ARGS["max_chr"], title=ARGS["title"], sort_by=ARGS["sort_by"],
                save=ARGS['save'])
