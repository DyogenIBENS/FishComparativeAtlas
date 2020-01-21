*** MacrosyntenyTGD/Data/Orthologs/

This directory contains 1-to-1 orthologs (or 1-to-2 orthologs between non- and post-TGD species) in Ensembl gene trees (release 75). The position of the genes can be found in map.txt file, which was extracted from Ensembl gene tree NHX file.

Species
HUMAN	Human
MOUSE	Mouse
CANFA	Dog
MONDO	Opossum
CHICK	Chicken
MELGA	Turkey
POEGU	Zebra finch
ANOCA	Anole lizard
LEPOC	Spotted gar
BRARE	Zebrafish
ORYLA	Medaka
GASAC	Stickleback
TETNG	Tetraodon


*** MacrosyntenyTGD/Data/NonTGDSegments/

Non-TGD segments are shown in [species].txt files with the following information for each segment.
1. Non-TGD segment name
2. Chromosome name
3. Start position
4. End position


*** MacrosyntenyTGD/Data/Macrosynteny/

Each file summarizes macrosynteny information.
The first line shows the species names used in the segmentation step.
The second line shows a non-TGD segment name.
The following lines list non-TGD genes in the segment and their orthologs: individual columns show orthologs in the respective species in the first line, and each number indicates the chromosome index of an ortholog.
(Genes with no orthologs were deleted in this file.)
Then, information for the next segment follows after a blank line.

Chromosomes are arranged in the following order and indexes start from zero (e.g., an ortholog located on medaka chr1 is represented by "0").
CHICK: chr1, ..., chr28, chr32, chrW, chrZ, LGE22C19W28_E50C23, LGE64
MELGA: chr1, ..., chr30, chrW
POEGU: chr1, chr1A, chr1B, chr2, chr3, chr4, chr4A, chr5, ..., chr28, chrZ, LG2, LG5, LG22
ANOCA: chr1, ..., chr6, LGa, LGb, LGc, LGd, LGf, LGg, LGh
BRARE: chr1, ..., chr25
ORYLA: chr1, ..., chr24
GASAC: chrI, chrII, chrIII, ..., chrXXI
TETNG: chr1, ..., chr21


*** MacrosyntenyTGD/Results/K=13/

-- The inference result is written in est_Usk.txt file. This file shows estimated values E[\hat{U}_{s,k}] (k=1,...,13) for each non-TGD segment s.

-- Individual [speciesName].seg.txt files contain information for painting the respective genomes. Each line shows the following information.
1. Chromosome name
2. Start position
3. End position
4. Value of k (i.e., inferred pre-TGD chromosome)

-- The macrosynteny model is symmetric with respect to all permutations of k=1,...,13. As a result, two chromosomes with the same value of k in two different reconstructions (e.g., with K=12 and 13) are not guaranteed to be the same pre-TGD chromosome. Therefore, the values of k in the above files are not directly associated with chr1,...chr13 in the figures. Instead, the association is defined in ColorCoding.txt file, in which columns indicate
1. value of k,
2. pre-TGD chromosome name,
3. RGB color, and
4. color name.

We chose these RGB colors following "Color Universal Design" guidelines (see below for references). As a result, these colors are not consistent with the color coding in previously published papers.
References:
- http://jfly.iam.u-tokyo.ac.jp/color/ (English)
- http://jfly.iam.u-tokyo.ac.jp/colorset/ (Japanese. The RGB colors were chosen from here.)

-- Genomes were painted according to [speciesName].seg.txt and shown in [speciesName].pdf, preWGD.pdf, and genomePaint.pdf files.

-- Orthologs were plotted in orthologPlot_nonTGD.pdf file, which visualizes ortholog distributions between non-TGD segments and post-TGD genomes.


*** MacrosyntenyTGD/Results/K=12/

An alternative reconstruction with K=12 is shown in this directory.



