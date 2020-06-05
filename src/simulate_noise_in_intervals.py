"""
    Script to inject random noise in input pre-TGD segments to asses robustness of the pipeline.

    Example:
        $

"""

import argparse
import numpy as np
from scripts.synteny.mygenome import Genome


from color_reference_species import load_nakatani_segments, genes_to_segments
from convert_intervals import segments_to_genes


def apply_noise(dseg, sigma):

    """
    Add noise to labelled intervals on a chromosome: gaussian noise on borders.

    Args:
        dseg (dict): labelled intervals
        sigma (float): noise parameter

    Returns:
        dict: intervals after noise application
    """

    #sort intervals by start
    sorted_intervals = sorted(dseg.keys(), key=lambda x: x[1])
    stored_na = []

    #create a list of contiguous intervals, if hole --> store to remove it later
    intervals_to_list = []
    ancs = []
    for interval in sorted_intervals:
        ch, beg, end = interval
        anc = dseg[interval]

        #if first store start and end
        if not intervals_to_list:
            intervals_to_list += [beg, end]
            ancs.append(anc)

        #if interval directly follows previous interval, store only end
        elif beg == intervals_to_list[-1]+1:

            intervals_to_list += [end]
            ancs.append(anc)

        #otherwise, store gap as start and end: gap will be kept as gaps in results (==no label)
        else:
            intervals_to_list += [beg, end]
            stored_na.append((beg, end))
            ancs.append('na')
            ancs.append(anc)

    #add nois to interval borders
    intervals_noise = gaussian_noise(np.array(intervals_to_list[1:-1]), sigma)

    #convert limits back to integer
    intervals_noise = np.rint(intervals_noise)

    #do not change chromosome start and stop
    intervals_noise = [intervals_to_list[0]]+list(intervals_noise)+[intervals_to_list[-1]]

    #sort new limits (small intervals may have been inverted)
    intervals_noise = sorted(intervals_noise)

    #map corresponding labels back
    dseg_noise = mapping_back_anc(intervals_noise, ancs, ch)

    return dseg_noise


def mapping_back_anc(intervals, ancs, ch):

    """
    Maps predicted labels to intervals after random perturbation.

    Args:
        intervals (list): list of perturbed intervals (sorted)
        ancs (list): labels (sorted, i.e label[0] correspond to label of intervals[0] before
                     randomization)
        ch (str): name of current chrom, to store along with interval start and stop

    Returns:
        dict: Correspondence randomized interval (key) <--> label (value)
    """

    dseg = {}
    for i in range(len(intervals)-1):
        (beg, end) = intervals[i:i+2]
        if ancs[i] != 'na':
            dseg[(ch, (int(beg), int(end)))] = ancs[i]

    return dseg

def gaussian_noise(x, sigma):

    """
    Adds gaussian noise to all elements of a list.

    Args:
        x (numpy.array):
        sigma (int): gaussian noise parameter

    Returns:
        list: List after gaussian noise injection

    """

    return list(np.random.normal(x, sigma))


def write_intervals(dseg, output):

    """
    Writes intervals + labels to file.

    Args:
        dseg (dict): interval (key) and label (value)
        output (str): output name
    """

    with open(output, 'w') as outfile:

        for interval in dseg:

            anc = dseg[interval]
            to_write = '\t'.join([str(i) for i in interval])
            to_write += '\t'+ anc +'\n'
            outfile.write(to_write)


def write_genomic_intervals_from_genes(dseg, output):
    """
    """
    with open(output, 'w') as outfile:

        for interval in dseg:

            anc = dseg[interval]
            to_write = '\t'.join([str(i) for i in interval])
            to_write += '\t'+ anc +'\n'
            outfile.write(to_write)


def genes_intervals_to_dummy_coord(dgenes_seg, dseg):
    coord_g = {}
    d_genes = {}
    prev_ch = ''
    for interval in dgenes_seg:
        ch, st, end = interval


        if prev_ch and ch == prev_ch:
            start = prev_stop+1
            stop = prev_stop+len(dgenes_seg[interval])+1
            coord_g[(ch, start, stop)] = dseg[interval]
            j = 0
            for i in range(start, stop):
                d_genes[ch][i] = dgenes_seg[interval][j]
                j += 1

            prev_stop = prev_stop+len(dgenes_seg[interval])+1

        else:
            prev_stop = len(dgenes_seg[interval])-1
            coord_g[(ch, 0, prev_stop)] = dseg[interval]
            for i in range(prev_stop+1):
                d_genes[ch] = d_genes.get(ch, {})
                d_genes[ch][i] = dgenes_seg[interval][i]

        prev_ch = ch
    return coord_g, d_genes


def genes_to_segments2(genes_coord, dseg, dgenes=None, transform=True):

    dgenes_seg = {}
    for chromosome in genes_coord.keys():

        for i, gene in enumerate(genes_coord[chromosome]):
            chrom = str(chromosome)
            if "group" in chrom:
                chrom = chrom.replace("group", "")

            for seg in dseg:
                if chrom == seg[0]:
                    for interval in seg[1:]:
                        if i <= interval[1] and i >= interval[0]:
                            if transform:
                                dgenes_seg[genes_coord[chromosome][gene]] = dseg[seg]
                            else:
                                gene_name = genes_coord[chromosome][gene]
                                chrom_name = chromosome
                                # if gr:
                                # chrom_name = "group" + chromosome
                                # print(chrom_name)
                                # print([i.names[0] for i in dgenes.genes_list[chrom_name]])
                                try:
                                    genen = [i for i in dgenes.genes_list[int(chrom_name)]\
                                             if i.names[0] == gene_name]
                                except ValueError:
                                    genen = [i for i in dgenes.genes_list[chrom_name]\
                                             if i.names[0] == gene_name]

                                if not genen:
                                    chrom_name = "group" + chromosome
                                    genen = [i for i in dgenes.genes_list[chrom_name]\
                                             if i.names[0] == gene_name]

                                if not genen:
                                    continue

                                genen = genen[0]
                                dgenes_seg[(chrom, interval)] = dgenes_seg.get((chrom, interval),
                                                                               [])
                                dgenes_seg[(chrom, interval)].append((genen.beginning, genen.end,
                                                                      dseg[seg]))


    return dgenes_seg



if __name__ == '__main__':

    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-seg', '--ancestral_seg', type=str, help='Nakatani et al. ancestral\
                        chromosomes', required=True)

    PARSER.add_argument('-o', '--outfile', type=str, help='output file', required=False,
                        default='out')

    PARSER.add_argument('-g', '--genes', type=str, help='Genes file', required=False, default='')


    PARSER.add_argument('-s', '--sigma', type=float, help='gaussian noise parameter',
                        required=False, default=5)

    PARSER.add_argument('-f', '--genesformat', type=str, required=False,
                        default='bed')

    PARSER.add_argument('-ori', '--save_ori', type=str, required=False,
                        default='')

    PARSER.add_argument('--genes_intervals', action='store_true')

    ARGS = vars(PARSER.parse_args())

    SEG = load_nakatani_segments(ARGS["ancestral_seg"], agg=False)

    GENES = Genome(ARGS["genes"], ARGS["genesformat"])

    CHROMS = {i for i, _, _ in SEG}

    NOISY_SEGS = {}


    if not ARGS["genes_intervals"]:


        for chrom in CHROMS:
            #TODO: improve this
            SEG_CHROMS = {(i, j, k):SEG[(i, j, k)] for i, j, k in SEG if i == chrom}
            NOISY_SEG = apply_noise(SEG_CHROMS, ARGS["sigma"])
            NOISY_SEGS.update(NOISY_SEG)
        # write_intervals(NOISY_SEGS, ARGS["outfile"])

        genes = genes_to_segments(GENES, NOISY_SEGS, transform=True)


    else:
        # SEG = {(i, (j, k)):SEG[(i, j, k)] for i, j, k in SEG.keys()}
        genes_int = segments_to_genes(GENES, SEG)
        genes_int, genes_coord = genes_intervals_to_dummy_coord(genes_int, SEG)

        for chrom in CHROMS:
            #TODO: improve this
            SEG_CHROMS = {(i, j, k):genes_int[(i, j, k)] for i, j, k in genes_int if i == chrom}
            if SEG_CHROMS:
                NOISY_SEG = apply_noise(SEG_CHROMS, ARGS["sigma"])
                NOISY_SEGS.update(NOISY_SEG)

        genes = genes_to_segments2(genes_coord, NOISY_SEGS)

        INTERVALS = genes_to_segments2(genes_coord, NOISY_SEGS, GENES, transform=False)
        with open(ARGS["outfile"]+"_genomic", 'w') as out:
            for chrom in CHROMS:
                SEG_CHROMS = {(i, j, k):INTERVALS[(i, (j, k))] for i, (j, k) in INTERVALS
                              if i == chrom}
                i = 0
                for interval in sorted(SEG_CHROMS.keys(), key=lambda x: x[1]):
                    sorted_genes = sorted(SEG_CHROMS[interval], key=lambda x: x[0])

                    if i == 0:
                        start, _, anc = sorted_genes[0]
                        _, end, _ = sorted_genes[-1]
                        prev_st = start
                        prev_e = end

                    else:
                        start = prev_e
                        _, end, _ = sorted_genes[-1]
                        prev_e = end

                    out.write(chrom+'\t'+str(start)+'\t'+str(end)+'\t'+anc+'\n')
                    i += 1


    with open(ARGS["outfile"], 'w') as fw:
        for name in genes:

            col = ''.join(genes[name])
            # print(col, name)

            fw.write(name+'\t'+col+'\n')


    if ARGS["save_ori"]:
        SEG = {(i, (j, k)):SEG[(i, j, k)] for i, j, k in SEG} #TODO: improve this

        genes = genes_to_segments(GENES, SEG, transform=True)

        with open(ARGS["save_ori"], 'w') as fw:
            for name in genes:
                col = ''.join(genes[name])
                # print(col, name)
                fw.write(name+'\t'+col+'\n')
