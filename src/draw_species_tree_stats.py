"""
TODO docstring

python src/draw_species_tree_stats.py -i
PM_Genofish_GENOMICUSV3_noPeri/out
PM_Genofish_GENOMICUSV3/out_newPM_Genofish_GENOMICUSV3_nocorr/out_new  -l
"GenomicusV3 + Edition + SCORPiOs" "GenomicusV3 + SCORPiOs" "Genomicus V3"
-s ../SCORPiOs/data/genofish_v3/sptree_without_gobidae.nwk
-ob PM_Genofish_GENOMICUSV3_noPeri/stats_boxplots.svg
-os PM_Genofish_GENOMICUSV3_noPeri/sptree_stats.svg -a Neopterygii -da Osteoglossocephalai
"""

import argparse

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd

from ete3 import Tree, TreeStyle, NodeStyle, TextFace



sns.set_palette("muted")

def load_annotation_stats(input_file):

    """
    Loads annotation statistics into a dict.

    Args:

    input_file (str): annotation file, as written by src/plot_paraogy_map.py

    Returns
        dict: for each species, the % of annotated genes in the paralogy Map.
    """

    stats = {}
    with open(input_file, "r") as infile:
        for line in infile:
            if line.strip():
                sp = line.split(":")[0]
                annotated = float(line.split("(")[1].split("%")[0])
                stats[sp] = annotated
    return stats


def plot_species_tree_with_stats(sptree, stats, out="out_sp.svg", dupanc=None, c_name="viridis",
                                 anc=None):
    """
    Draws a circular species tree with species colored according to a continuous variable, here
    the proportion of the species genes annotated by the paralogy_map workflow.

    Args:
        sptree (str): Newick string or Newick file name of the species tree.
        stats (dict): For a species (key) gives its associated value
        out (str, optional): output figure file name
        dupanc (str, optional): Name of an ancestor species to highlight on the tree
        c_name (str, optional): Name of the maptlotlib continuous colormap to use
        anc (str, optional): Ancestor name to only plot the subtree below it.

    """

    tree = Tree(sptree, format=1)

    if anc:
        tree = tree.search_nodes(name=anc)[0]

    cmap = plt.cm.get_cmap(c_name)
    min_a, max_a = min(stats.values()), max(stats.values())
    norm = matplotlib.colors.Normalize(vmin=min_a, vmax=max_a)

    circular_style = TreeStyle()
    circular_style.show_leaf_name = False
    circular_style.scale = 20
    circular_style.mode = 'c'
    for node in tree.traverse():
        if not node.is_leaf():
            nstyle = NodeStyle()

            if (dupanc and node.name != dupanc) or not dupanc:
                nstyle["fgcolor"] = "lightgrey"

            else:
                nstyle["fgcolor"] = "lightcoral"
                nstyle["shape"] = "square"
            nstyle["size"] = 10
            node.set_style(nstyle)

        else:
            annotated = str(stats.get(node.name, ''))
            if annotated:
                norm_a = norm(float(annotated))
                col = matplotlib.colors.to_hex(cmap(norm_a))
            else:
                col = 'black'

            node.name = node.name+" "+annotated
            nstyle = NodeStyle()
            nstyle["size"] = 0
            name_face = TextFace(node.name, fgcolor=col)
            node.add_face(name_face, column=0)
            node.set_style(nstyle)

    tree.render(out, dpi=200, tree_style=circular_style)


def boxplot_stats(datasets, out="out_box.svg", yname="Proportion of genome annotated (%)"):

    """
    Draws boxplots of the distribution of a continuous variable associated to species.

    Args:
        datasets
        out
        yname
    """

    for i, data in enumerate(datasets):
        label, stats = data
        df = pd.DataFrame(stats.items(), columns=['Species', yname])
        df["Trees"] = label
        if i > 0:
            df_all = pd.concat([df_all, df])
        else:
            df_all = df

    plt.figure(figsize=(3, 5))
    ax = sns.boxplot(data=df_all, y=yname, x="Trees", width=0.5)

    ax.set_xticklabels(ax.get_xticklabels(), fontsize=9, rotation=60, ha='right')
    sns.despine()
    plt.xlabel("")

    plt.tight_layout()
    plt.savefig(out, dpi=100)


if __name__ == '__main__':

    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    ## Required ##
    PARSER.add_argument('-i', '--input', type=str, nargs='+', help='Input annotation statistics,\
                        several files can be given',
                        required=True)

    ## Optional ##
    PARSER.add_argument('-l', '--labels', type=str, nargs='+', help='Labels for input datasets',
                        required=False, default='')

    PARSER.add_argument('-s', '--sptree', type=str, help="Species tree file, in newick. If not\
                        given, the species tree isn't plotted.", required=False, default='')

    PARSER.add_argument('-ob', '--output_box', type=str, help='output file for boxplot figure',
                        required=False, default='out_box.svg')

    PARSER.add_argument('-os', '--output_sp', type=str, help='output file for species tree figure',
                        required=False, default='out_sp.svg')

    PARSER.add_argument('-a', '--ancestor', type=str, help='Cut species at `ancestor`',
                        required=False, default='')

    PARSER.add_argument('-da', '--dupancestor', type=str, help='To draw wgd node in species tree',
                        required=False, default='')

    ARGS = vars(PARSER.parse_args())

    if not ARGS["labels"]:
        ARGS["labels"] = []
        for j, _ in enumerate(ARGS["input"]):
            ARGS["labels"].append(f"Dataset {j}")

    DATA_STATS = []
    for j, input_data in enumerate(ARGS["input"]):
        dstats = load_annotation_stats(input_data)
        lab = ARGS["labels"][j]
        DATA_STATS.append((lab, dstats))

        if j == 0 and ARGS["sptree"]:
            plot_species_tree_with_stats(ARGS["sptree"], dstats, out=ARGS["output_sp"],
                                         dupanc=ARGS["dupancestor"], anc=ARGS["ancestor"])

    boxplot_stats(DATA_STATS, out=ARGS["output_box"])
    plt.close("all")
