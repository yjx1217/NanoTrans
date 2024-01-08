
import argparse
import numpy as np
import pandas as pd

import pysam
from collections import Counter
import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42
#plt.rcParams['font.family'] = 'Arial'
# ============================================================================ #
# parse args

parser = argparse.ArgumentParser(description='Plot frequency of motifs around stop codons')
parser.add_argument('-r', '--refbed',  help='reference annotation in bed format, required', required=True)
parser.add_argument('-s', '--searchbed', help='motifs to search in bed format, required', required=True)
parser.add_argument('-o', '--outname',   help='output directory, required', default="./",  required=True)

args = parser.parse_args()

genes_bed = args.refbed

# ============================================================================ #
# functions

def get_cds_region(bed, chrom, stop_pos, gene_strand, use_score=False, use_strand=True):
    region = np.zeros(1000, dtype=np.int32)
    try:
        bed_iter = bed.fetch(chrom, stop_pos - 500, stop_pos + 500)  # extract records located in the search region from bed
    except ValueError:
        return region
    for record in bed_iter:
        record = record.split()
        pos = int(record[1])  # abs pos of m6a
        if pos > (stop_pos - 500) and pos < (stop_pos + 500):
            strand = record[5]
            if use_strand:
                if strand == gene_strand:
                    idx = pos - stop_pos + 500  # the region starts with -500 of stop codon
                    region[idx] += float(record[4]) if use_score else 1
            else:
                idx = pos - stop_pos + 500
                region[idx] += float(record[4]) if use_score else 1
    return region

def get_stop_profiles(tabix_file, use_score=False, use_strand=True):
    stop_profiles = []
    mismatches = pysam.TabixFile(tabix_file)
    n_records = 0
    with open(genes_bed) as bed:
        for record in bed:
            record = record.split()
            chrom = record[0].replace('Chr', '')
            if chrom in ['Mt', 'Pt']:
                continue
            strand = record[5]
            cds_end = int(record[7]) if strand == '+' else int(record[6])
            stop_prof = get_cds_region(mismatches, chrom, cds_end, strand, use_score=use_score, use_strand=use_strand)
            if strand == '-':
                stop_prof = stop_prof[::-1]
            n_records += 1
            stop_profiles.append(stop_prof)
    mismatches.close()
    return stop_profiles, n_records

der_site_stop_profiles, n_records_der = get_stop_profiles(args.searchbed)
# ============================================================================ #
# plot
fig, ax = plt.subplots(figsize=(5, 3))
ax.plot(np.sum(der_site_stop_profiles, axis=0), color="steelblue", lw=2, zorder=1)
#ax2 = ax.twinx()
ax.axvline(500, color='#555555', ls='--', lw=2, zorder=-1)
ax.set_xticks([250, 500, 750])
#plt.setp(ax2.get_yticklabels(), color="blue")
ax.set_xticklabels(['-250nt', 'StopCodon', '+250nt'])
ax.set_xlabel('    ')
ax.set_ylabel('Error sites (frequency)', color="steelblue")
ax.set_xlim(200, 800)
plt.tight_layout()
plt.savefig(args.outname, format = "pdf")
