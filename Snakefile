import sys
import os

from scripts.trees import speciestree as spt

#################################################################################################
#Runs in the conda env paralogy_map + path to SCORPiOs scripts should be appended to python path #
#Dependencies (see env file) = snakemake, ete3, matpotlib and seaborn                            #
#################################################################################################

DUPLICATED_SPECIES = spt.get_species(config["species_tree"], config["ancestor"])
config["prune_ancestor"] = config.get("prune_ancestor", "Neopterygii")
REF_SPECIES = ["Oryzias.latipes", "Gasterosteus.aculeatus", "Tetraodon.nigroviridis", "Danio.rerio"]
SEGMENTS = config.get("seg", "data/MacrosyntenyTGD/Results/K=13/{ref_species}.seg.txt")
SEGMENTS_OK = config.get("seg_ok", "data/MacrosyntenyTGD/Results/K=13/{ref_species}_ok.seg.txt")
GENES = config["genes"]


rule all:
    """
    Target of the workflow:
    an .svg image for each dup species, with their genome colored by post-duplication chromosomes
    """
    input: f'{config["jobname"]}/sptree_stats.svg',
           f'{config["jobname"]}/box_stats.svg'


rule convert_intervals:
    """
    Maps intervals of current assembly to intervals of previous assembly used in Nakatani and McLysaght.
    """
    input: s = SEGMENTS, g = GENES
    output: SEGMENTS_OK
    run:

        if wildcards.ref_species in config["genes_conv"].keys():

            cmd1 = f"sed 's/Old stable ID, New stable ID, Release, Mapping score//g'\
                    data/ensembl_id_history_{wildcards.ref_species}_raw.csv | grep -v '^[[:space:]]*$'\
                    > data/ensembl_id_history_{wildcards.ref_species}.csv;"

            os.system(cmd1)

            cmd2 = f'python src/convert_intervals.py -g {config["genes_conv"][wildcards.ref_species]}\
                    {input.g} -seg {input.s} -id data/ensembl_id_history_{wildcards.ref_species}.csv\
                    -f {config["gconv_format"]} {config["format"]} -o {output[0]}'

            os.system(cmd2)

        else:
            os.system(f"cp {input[0]} {output[0]}")


rule extract_duplicated_ancGenes:
    """
    Extracts all post-duplication ancgenes in the input gene trees.
    """
    input: trees = config.get("forest", 'test'), sptree = config["species_tree"]
    output: f'{config["jobname"]}/TGD_ancGenes.tsv'
    shell:
        "python src/get_post_dup_ancgenes.py -t {input.trees} -d {config[ancestor]} "
        "-s {input.sptree} -o {output} --check_root"


rule color_each_reference:
    """
    Identifies paralogous duplicated segments within each of the 4 reference species.
    Uses paralogous genes in input gene trees and pre-TGD segments.
    """
    input: segments = SEGMENTS_OK, genes = GENES, ancGenes = f'{config["jobname"]}/TGD_ancGenes.tsv'
    output: f'{config["jobname"]}/{{ref_species}}_colors.txt'
    shell:
        "python src/color_reference_species.py -seg {input.segments} -g {input.genes} "
        "-ag {input.ancGenes} -o {output} -f {config[format]}"


#TODO: touching a single output was simpler to implement in snakemake but we can do better
rule homogenize_references_ab:
    """
    Homogenizes nomenclature of paralogous segments within each references to be consistent across
    species (i.e orthologous segments should have the same name).
    """
    input:
        ref_colors = expand(f'{config["jobname"]}/{{ref_species}}_colors.txt',
                              ref_species=REF_SPECIES),
        genes = expand(GENES, ref_species=REF_SPECIES),
        ancGenes = f'{config["jobname"]}/TGD_ancGenes.tsv'
    output: f'{config["jobname"]}/touched_file'
    params: guide = "Gasterosteus.aculeatus"
    shell:
        "python src/homogenize_refs_colors.py -i {input.ref_colors} -ag {input.ancGenes} "
        "-guide_sp {params.guide} -g {input.genes} -f {config[format]}; touch {output}"


rule consensus_color_ancGene:
    """
    Assigns, through a majority vote of the 4 reference species descendant genes, an ancestral
    post-duplication chromosome to each ancestral gene.
    """
    input: ref_colors = f'{config["jobname"]}/touched_file',
           genes = expand(GENES, ref_species=REF_SPECIES),
           ancGenes = f'{config["jobname"]}/TGD_ancGenes.tsv'

    output: f'{config["jobname"]}/colored_TGD_ancGenes.tsv'

    params:
            ref_colors = expand(f'{config["jobname"]}/{{ref_species}}_colors_homogenized.txt',
                                ref_species=REF_SPECIES)
    shell:
        "python src/color_ancgenes.py -ref {params.ref_colors} -g {input.genes} "
        "-ag {input.ancGenes} -o {output} -f {config[format]}"


rule draw_paralogy_map:
    """
    Draws the genome of each duplicated species colored by ancestral post-duplication chromosomes,
    using ancestral genes annotated in the previous rules.
    """
    input: colors = f'{config["jobname"]}/colored_TGD_ancGenes.tsv',
           genes = GENES.replace('{ref_species}', '{dup_species}')

    output: plot = report(f'{config["jobname"]}/{{dup_species}}_ParalogyMap.svg', category="Paralogy Maps"),\
                   # caption="Paralogy map for {wildcards.dup_species}",\
            stats = temp(f'{config["jobname"]}/{{dup_species}}_out_stats.txt')

    params: draw = config.get('draw', '')

    shell:
        "python src/plot_paralogy_map.py -c {input.colors} -g {input.genes} -o {output.plot} "
        "-s {wildcards.dup_species} -f {config[format]} -os {output.stats} {params.draw}"


rule stats:
    input:
        st = expand(f'{config["jobname"]}/{{dup_species}}_out_stats.txt',
                    dup_species=DUPLICATED_SPECIES),
        fig = expand(f'{config["jobname"]}/{{dup_species}}_ParalogyMap.svg',
                      dup_species=DUPLICATED_SPECIES)

    output: f'{config["jobname"]}/out_stats.txt'

    shell:
        "cat {input.st} > {output}"


if "comparisons" in config:
    rule plot_annotation_statistics1:
        """
        Plots proportion of genome annotated.
        """
        input: stats = f'{config["jobname"]}/out_stats.txt',

        output: boxplots = report(f'{config["jobname"]}/box_stats.svg', category="Annotation statistics"),\
                       # caption="Proportion of genomes annotated and comparisons with previous results",\
                sptree = report(f'{config["jobname"]}/sptree_stats.svg', category="Annotation statistics")
                       # caption="Visualizaiton of annotation statistics on the species tree",\
                       
        shell:
            "python src/draw_species_tree_stats.py -i {input.stats} {config[comparisons]} "
            "-l {config[labels]} -s {config[species_tree]} -ob {output.boxplots} "
            "-os {output.sptree} -a {config[prune_ancestor]} -da {config[ancestor]}"

else:
    rule plot_annotation_statistics2:
        """
        Plots proportion of genome annotated.
        """
        input: stats = f'{config["jobname"]}/out_stats.txt',

        output: boxplots = report(f'{config["jobname"]}/box_stats.svg', category="Annotation statistics"),\
                       # caption="Proportion of genomes annotated and comparisons with previous results",\
                sptree = report(f'{config["jobname"]}/sptree_stats.svg', category="Annotation statistics")
                       # caption="Visualizaiton of annotation statistics on the species tree",\
                       
        shell:
            "python src/draw_species_tree_stats.py -i {input.stats} -s {config[species_tree]} "
            "-ob {output.boxplots} -os {output.sptree} -a {config[prune_ancestor]} "
            "-da {config[ancestor]}"

# rule compare_zfin:
# rule plot_zfin: