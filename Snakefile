import sys

from scripts.trees import speciestree as spt

#################################################################################################
#Runs in the conda env paralogymap + path to SCORPiOs scripts should be appended to python path #
#Dependencies = snakemake, ete3, matpotlib and seaborn                                          #
#################################################################################################

DUPLICATED_SPECIES = spt.get_species(config["species_tree"], config["ancestor"])
REF_SPECIES = ["Oryzias.latipes", "Gasterosteus.aculeatus", "Tetraodon.nigroviridis", "Danio.rerio"]
SEGMENTS = "data/MacrosyntenyTGD/Results/K=13/{ref_species}.seg.txt"
GENES = config["genes"]

#TODO: more stderr prints
#OPTIMIZE: total run 75 teleost 8 minutes on 14 cores, pyflame src to see where to gain time (there must be plenty)

rule all:
    """
    Target of the workflow:
    an .svg image for each dup species, with their genome colored by post-duplication chromosomes
    """
    input: expand(config["jobname"]+"/{dup_species}_ParalogyMap.svg", dup_species=DUPLICATED_SPECIES)

rule extract_duplicated_ancGenes:
    """
    Extracts all post-duplication ancgenes in the input gene trees.
    """
    input: trees = config.get("forest", 'test'), sptree = config["species_tree"]
    output: config["jobname"]+"/TGD_ancGenes.tsv"
    shell:"""
    python scripts/get_post_dup_ancgenes.py -t {input.trees} -d {config[ancestor]} -s {input.sptree} -o {output}
    """

rule color_each_reference:
    """
    Identifies paralogous duplicated segments within each of the 4 reference species.
    Uses paralogous genes in input gene trees.
    """
    input: segments = SEGMENTS, genes = GENES, ancGenes = config["jobname"]+"/TGD_ancGenes.tsv"
    output: config["jobname"]+"/{ref_species}_colors.txt"
    shell:"""
    python scripts/color_reference_species.py -seg {input.segments} -g {input.genes} -ag {input.ancGenes}\
                                              -o {output} -f {config[format]}
    """

def get_ref_colors(wildcards):
    return expand(config["jobname"]+"/{ref_species}_colors.txt", ref_species=REF_SPECIES)

def get_genes(wildcards):
    return expand(GENES, ref_species=REF_SPECIES)

def get_ref_colors2(wildcards):
    return expand(config["jobname"]+"/{ref_species}_colors_homogenized.txt", ref_species=REF_SPECIES)

#TODO: touching a single output was simpler to implement here but we can do better
rule homogenize_reference_ab:
    """
    Homogenizes nomenclature of paralogous segments within each references to be consistent across
    species (i.e orthologous segments should have the same name).
    """
    input: ref_colors = get_ref_colors, genes = get_genes,
           ancGenes = config["jobname"]+"/TGD_ancGenes.tsv"
    output: config["jobname"]+"/touched_file"
    params: guide = "Oryzias.latipes"
    shell:"""
    python scripts/homogenize_refs_colors.py -i {input.ref_colors} -ag {input.ancGenes} -guide_sp {params.guide}\
                                             -g {input.genes} -f {config[format]};
                                             touch {output}
    """

rule consensus_color_ancGene:
    """
    Assigns, through a majority vote of the 4 reference species descendant genes, an ancestral
    post-duplication chromosome to each ancestral gene.
    """
    input: ref_colors = config["jobname"]+"/touched_file",\
           genes = get_genes, ancGenes = config["jobname"]+"/TGD_ancGenes.tsv"
    output: config["jobname"]+"/colored_TGD_ancGenes.tsv"
    params: ref_colors = get_ref_colors2
    shell:"""
    python scripts/color_ancgenes.py -ref {params.ref_colors} -g {input.genes} -ag {input.ancGenes}\
                                     -o {output} -f {config[format]};
    """

rule draw_paralogy_map:
    """
    Draws the genome of each duplicated species colored by ancestral post-duplication chromosomes,
    using ancestral genes annotated in the previous rules.
    """
    input: colors = config["jobname"]+"/colored_TGD_ancGenes.tsv",\
           genes = GENES.replace('{ref_species}', '{dup_species}')
    output: config["jobname"]+"/{dup_species}_ParalogyMap.svg"
    shell:"""
    python scripts/plot_paralogy_map.py -c {input.colors} -g {input.genes} -o {output} -s {wildcards.dup_species}\
                                        -f {config[format]}
    """
