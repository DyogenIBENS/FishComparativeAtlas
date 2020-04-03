from ete3 import Tree, TreeStyle, NodeStyle, TextFace

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd

t = Tree("../SCORPiOs/data/genofish_v3/GenomicusV3/GenomicusV3_speciesTree.nwk", format=1)

neo = t.search_nodes(name="Neopterygii")[0]
d = {}
with open("../paralogy_map/PM_Genofish_GENOMICUSV3_nocorr/out_new", "r") as infile:
    for line in infile:
        if line.strip():
            sp = line.split(":")[0]
            a = float(line.split("(")[1].split("%")[0])
            d[sp] = a


d_cor = {}
with open("../paralogy_map/PM_Genofish_GENOMICUSV3/out_new", "r") as infile:
    for line in infile:
        if line.strip():
            sp = line.split(":")[0]
            a = float(line.split("(")[1].split("%")[0])
            d_cor[sp] = a

min_a = min(d.values())
max_a = max(d.values())

cmap = plt.cm.get_cmap('viridis')
norm = matplotlib.colors.Normalize(vmin=min_a, vmax=max_a)

circular_style = TreeStyle()
circular_style.show_leaf_name = False
circular_style.scale = 20
circular_style.mode = 'c'
for n in t.traverse():
    if not n.is_leaf():
        nstyle = NodeStyle()
        if n.name != "Osteoglossocephalai":
            nstyle["fgcolor"] = "lightgrey"
        else:
            nstyle["fgcolor"] = "lightcoral"
            nstyle["shape"] = "square"
        nstyle["size"] = 10
        n.set_style(nstyle)
    else:
        a = str(d.get(n.name, ''))
        if a:
            norm_a = norm(float(a))
            col = matplotlib.colors.to_hex(cmap(norm_a))
        else:
            col = 'black'
        n.name = n.name+" "+a
        nstyle = NodeStyle()
        nstyle["size"] = 0
        name_face = TextFace(n.name, fgcolor=col)
        n.add_face(name_face, column=0)
        n.set_style(nstyle)


# neo.render("PM_stats_species_tree.svg", tree_style=circular_style)
plt.figure(figsize=(3, 5))
data1 = pd.DataFrame(d_cor.items(), columns=['Species', "Proportion of genome annotated (%)"])
data1["Genofish trees"] = "SCORPiOs"
data2 = pd.DataFrame(d.items(), columns=['Species', "Proportion of genome annotated (%)"])
data2["Genofish trees"] = "Uncorrected"
data = pd.concat([data1, data2])
print(data)
ax = sns.boxplot(data=data, y="Proportion of genome annotated (%)", x="Genofish trees", width=0.5)
# plt.xlabel(["SCORPiOs corrected forest", "Uncorrected forest"])
# plt.ylabel("Proportion of genome annotated")
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
plt.tight_layout()
sns.despine()
plt.savefig("boxplot_stats_nocorr.svg", dpi=100)
plt.show()
