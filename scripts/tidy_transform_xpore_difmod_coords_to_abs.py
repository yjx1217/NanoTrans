# ============================================================================ #
# info

# date: 231023
# auth: yld
# desc: Transform relative pos of motifs mapped on transcriptome to absolute coordinates
# exam: python tidy_transform_difmod_coords_to_abs.py \
#           -r results/arabidopsis/Batch_Dataset2.all_samples_combined.flair_all_collapsed.isoforms.with_productivity.bed \
#           -s results/arabidopsis/difmod-filtered-bydiffrate.tsv \
#           -o results/arabidopsis/difmod_abs_coordinates.bed

# ============================================================================ #

# ============================================================================ #
# load

import argparse
import pandas as pd
from itertools import accumulate

# ============================================================================ #
# parse args

parser = argparse.ArgumentParser(description='Transform relative pos of motifs mapped on transcriptome to absolute coordinates')
parser.add_argument('-r', '--refbed',  help='reference annotation in bed format, required', required=True)
parser.add_argument('-s', '--searchbed', help='motifs to search in bed format, required', required=True)
parser.add_argument('-o', '--outname',   help='output directory, required', default="./",  required=True)

args = parser.parse_args()

refbed = pd.read_table(args.refbed, header  = None)
motifbed = args.searchbed
outname = args.outname
# ============================================================================ #
# test

if False:
    refbed = pd.read_table("results/arabidopsis/Batch_Dataset2.all_samples_combined.flair_all_collapsed.isoforms.with_productivity.bed", header  = None)
    motifbed = "results/arabidopsis/difmod-filtered-bydiffrate.tsv"
    outname = "results/arabidopsis/difmod_abs_coordinates.bed"

# ============================================================================ #
# functions

# motifpos = 1200
# isoformId = "AT1G78930.1_AT1G78930"
def trans_relpos_to_abspos(refbed, isoformId, motifpos):
    sel_bed = refbed[refbed[3] == isoformId]
    exons_size =  list(map(int, (sel_bed.iloc[0,10].split(",")[:-1])))
    exons_start = list(map(int, (sel_bed.iloc[0,11].split(",")[:-1])))
    chrom = sel_bed.iloc[0, 0]
    strand = sel_bed.iloc[0, 5]
    tx_start = int(sel_bed.iloc[0,1])
    if strand == "-":
        motifpos = sum(exons_size) - motifpos - 1 
    exons_accum_length = list(accumulate(exons_size))
    for i, v in enumerate(exons_accum_length):
        if v > motifpos:
            idx = i
            if idx == 0:
                relpos = motifpos
            else:
                relpos = motifpos - exons_accum_length[idx - 1]
            break
    abspos = tx_start + exons_start[idx] + relpos
    return chrom, abspos, strand

# ============================================================================ #
# transform coordinates


with open(motifbed) as bed:
    next(bed)
    chroms = []
    absposs = []
    strands = []
    rawinfos = []
    for record in bed:
        record = record.split("\t")
        isoform_id = record[0] + "_" + record[1]
        motifpos = int(record[3])
        raw_info = isoform_id + "_" + record[3] + "_" + record[4]
        chrom, abspos, strand = trans_relpos_to_abspos(refbed, isoform_id, motifpos)
        chroms.append(chrom)
        absposs.append(abspos)
        strands.append(strand)
        rawinfos.append(raw_info)

difmod_abs_bed = {'chrom': chroms, 'start': absposs, 'end': absposs, 'name': rawinfos, 'score': 1, 'strand': strands}

df_difmod_abs_bed = pd.DataFrame(difmod_abs_bed)

df_difmod_abs_bed.to_csv(outname, sep = "\t", index=False, header=False)

