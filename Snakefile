import sys
import os

from scripts.trees import speciestree as spt

## TODO: PUT SEGMENTS in config

#################################################################################################
#Runs in the conda env paralogy_map + path to SCORPiOs src should be appended to python path    #
#Dependencies (see env file) = snakemake, ete3, matpotlib and seaborn                           #
#################################################################################################

DUPLICATED_SPECIES = spt.get_species(config["species_tree"], config["ancestor"])
REF_SPECIES = ["Oryzias.latipes", "Gasterosteus.aculeatus", "Tetraodon.nigroviridis", "Danio.rerio"]
SEGMENTS = "data/MacrosyntenyTGD/Results/K=13/{ref_species}.seg.txt"
SEGMENTS_OK = "data/MacrosyntenyTGD/Results/K=13/{ref_species}_ok.seg.txt"
GENES = config["genes"]

#TODO: more stderr prints
#OPTIMIZE:(?) total run 75 teleost 8 minutes on 14 cores, pyflame src to see where to gain time

rule all:
    """
    Target of the workflow:
    an .svg image for each dup species, with their genome colored by post-duplication chromosomes
    """
    input: expand(config["jobname"]+"/{dup_species}_ParalogyMap.svg", dup_species=DUPLICATED_SPECIES)


#TODO: are non-zero exit status captured here?
rule convert_intervals:
    input: s = SEGMENTS, g = GENES
    output: SEGMENTS_OK
    run:

        if wildcards.ref_species in config["genes_conv"].keys():

            cmd1 = "sed 's/Old stable ID, New stable ID, Release, Mapping score//g'\
                    data/ensembl_id_history_"+wildcards.ref_species+"_raw.csv | grep -v '^[[:space:]]*$'\
                    > data/ensembl_id_history_"+wildcards.ref_species+".csv;"

            os.system(cmd1)

            cmd2 = "python src/convert_intervals.py -g "\
                    +config["genes_conv"][wildcards.ref_species]+" "+input.g+ " -seg "+input.s\
                    +" -id data/ensembl_id_history_"+wildcards.ref_species+".csv -f "\
                    +config["gconv_format"]+" "+config["format"]+" -o "+output[0]

            os.system(cmd2)

        else:
            os.system("cp "+input[0]+" "+output[0])

rule extract_duplicated_ancGenes:
    """
    Extracts all post-duplication ancgenes in the input gene trees.
    """
    input: trees = config.get("forest", 'test'), sptree = config["species_tree"]
    output: config["jobname"]+"/TGD_ancGenes.tsv"
    shell:"""
    python src/get_post_dup_ancgenes.py -t {input.trees} -d {config[ancestor]} -s {input.sptree}\
                                        -o {output}
    """


rule color_each_reference:
    """
    Identifies paralogous duplicated segments within each of the 4 reference species.
    Uses paralogous genes in input gene trees.
    """
    input: segments = SEGMENTS_OK, genes = GENES, ancGenes = config["jobname"]+"/TGD_ancGenes.tsv"
    output: config["jobname"]+"/{ref_species}_colors.txt"
    shell:"""
    python src/color_reference_species.py -seg {input.segments} -g {input.genes}\
                                          -ag {input.ancGenes} -o {output} -f {config[format]}
    """

def get_ref_colors(wildcards):
    return expand(config["jobname"]+"/{ref_species}_colors.txt", ref_species=REF_SPECIES)

def get_genes(wildcards):
    return expand(GENES, ref_species=REF_SPECIES)

def get_ref_colors2(wildcards):
    return expand(config["jobname"]+"/{ref_species}_colors_homogenized.txt", ref_species=REF_SPECIES)


#TODO: touching a single output was simpler to implement in snakemake but we can do better
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
    python src/homogenize_refs_colors.py -i {input.ref_colors} -ag {input.ancGenes}\
                                         -guide_sp {params.guide} -g {input.genes}\
                                         -f {config[format]}; touch {output}
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
    python src/color_ancgenes.py -ref {params.ref_colors} -g {input.genes} -ag {input.ancGenes}\
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
    python src/plot_paralogy_map.py -c {input.colors} -g {input.genes} -o {output}\
                                    -s {wildcards.dup_species} -f {config[format]}
    """
