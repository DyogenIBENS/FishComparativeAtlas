"""
TODO

type of noise to apply:
move interval boundaries? --> i'll start with that. careful to check for disappearing intervals
                          --> plot images to see what is going on
shuffle a few labels?
break regions?

"""
import numpy as np
import argparse
from scripts.synteny.mygenome import Genome


from color_reference_species import load_nakatani_segments, genes_to_segments

#call this on each chromosome
def apply_noise(dseg, sigma):
    #sort keys by lower bound
    sorted_intervals = sorted(dseg.keys(), key=lambda x:x[1])
    stored_na = []

    #create a list of contiguous intervals, if hole --> store to remove it later
    intervals_to_list = []
    ancs = []
    for interval in sorted_intervals:
        ch, beg, end = interval
        # print(ch, beg)
        anc = dseg[interval]
        if not intervals_to_list:
            intervals_to_list += [beg, end]
            ancs.append(anc)

        elif beg == intervals_to_list[-1]+1:

            intervals_to_list += [end]
            ancs.append(anc)

        else:
            intervals_to_list += [beg, end]
            stored_na.append((beg, end))
            ancs.append('na')
            ancs.append(anc)

    intervals_noise = gaussian_noise(np.array(intervals_to_list[1:-1]), sigma)
    intervals_noise = np.rint(intervals_noise)
    intervals_noise = [intervals_to_list[0]]+list(intervals_noise)+[intervals_to_list[-1]]
    intervals_noise = sorted(intervals_noise)
    dseg_noise = mapping_back_anc(intervals_noise, ancs, ch)
    return dseg_noise

def mapping_back_anc(intervals, ancs, ch):
    dseg = {}
    for i in range(len(intervals)-1):
        (beg, end) = intervals[i:i+2]
        if ancs[i] != 'na':
            dseg[(ch, (int(beg), int(end)))] = ancs[i]

    return dseg

def gaussian_noise(x, sigma):
    return list(np.random.normal(x, sigma))

def write_intervals(dseg, output):
    """
    """
    with open(output, 'w') as outfile:

        for interval in dseg:

            anc = dseg[interval]
            to_write = '\t'.join([str(i) for i in interval])
            to_write += '\t'+ anc +'\n'
            outfile.write(to_write)


if __name__ == '__main__':

    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-seg', '--ancestral_seg', type=str, help='Nakatani et al. ancestral\
                        chromosomes', required=True)

    PARSER.add_argument('-o', '--outfile', type=str, help='output file', required=False,
                        default='out')

    PARSER.add_argument('-g', '--genes', type=str, help='Genes file', required=False, default='')


    PARSER.add_argument('-s', '--sigma', type=float, help='gaussian noise parameter', required=False,
                        default=2000)

    PARSER.add_argument('-f', '--genesformat', type=str, required=False,
                        default='bed')

    PARSER.add_argument('-ori', '--save_ori', type=str, required=False,
                        default='')

    ARGS = vars(PARSER.parse_args())

    SEG = load_nakatani_segments(ARGS["ancestral_seg"], agg=False)

    GENES = Genome(ARGS["genes"], ARGS["genesformat"])

    CHROMS = set([i for i, _, _ in SEG.keys()])

    NOISY_SEGS = {}

    for chrom in CHROMS:
        SEG_CHROMS = {(i, j, k):SEG[(i, j, k)] for i, j, k in SEG.keys() if i == chrom} #TODO: improve this
        NOISY_SEG = apply_noise(SEG_CHROMS, ARGS["sigma"])
        NOISY_SEGS.update(NOISY_SEG)
    # write_intervals(NOISY_SEGS, ARGS["outfile"])

    genes = genes_to_segments(GENES, NOISY_SEGS, transform=True)

    with open(ARGS["outfile"], 'w') as fw:
        for name in genes:
            col = ''.join(genes[name])
            fw.write(name+'\t'+col+'\n')

    if ARGS["save_ori"] :
        SEG = {(i, (j, k)):SEG[(i, j, k)] for i, j, k in SEG.keys()} #TODO: improve this

        genes = genes_to_segments(GENES, SEG, transform=True)

        with open(ARGS["save_ori"], 'w') as fw:
            for name in genes:
                col = ''.join(genes[name])
                fw.write(name+'\t'+col+'\n')