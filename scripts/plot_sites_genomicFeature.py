
import numpy as np
import pandas as pd

import pysam
from collections import Counter
import matplotlib.pyplot as plt

import argparse

import os

plt.rcParams['pdf.fonttype'] = 42
#plt.rcParams['font.family'] = 'Arial'
# ============================================================================ #
# parse args

parser = argparse.ArgumentParser(description='Plot frequency of motifs around stop codons')
parser.add_argument('-r', '--refbed',  help='reference annotation in bed format, required', required=True)
parser.add_argument('-s', '--searchbed', help='motifs to search in bed format, required', required=True)
parser.add_argument('-o', '--outname',   help='output directory, required', default="./",  required=True)

args = parser.parse_args()

searchbed_idx = args.searchbed + ".tbi"
if os.path.exists(searchbed_idx):
    print("index file exists")
else:
    searchbed_idx = args.searchbed + ".csi"

# ============================================================================ #
# test

# class NamespaceMock:
#     def __init__(self, **kwargs):
#         self.__dict__.update(kwargs)
# 
# args = NamespaceMock(
#         refbed =    "results/arabidopsis/Batch_Dataset2.all_samples_combined.flair_all_collapsed.isoforms.with_productivity.bed",
#         searchbed = "results/arabidopsis/difmod_abs_coordinates.bed.gz",
#         outname =   "results/arabidopsis/plot_sites_genomic_feature.pdf"
#     )

# ============================================================================ #
# prep | keep chrs of refbed as same as searchbed

genes_bed = args.refbed
df = pd.read_csv(args.searchbed, sep='\t', header = None)
unique_chrs = set(df[df.columns[0]].astype(str))

# ============================================================================ #
# functions

def parse_exons_introns_flank(record, flanksize=200):
    start = int(record[1])
    end = int(record[2])
    exstarts = np.fromstring(record[11], sep=',') + start
    exends = exstarts + np.fromstring(record[10], sep=',')
    exons = np.dstack([exstarts, exends])[0]
    left_flank = np.array([[max(0, start - flanksize), start]])
    right_flank = np.array([[end, end + flanksize]])
    if len(exons) > 1:
        introns = np.dstack([exons[:-1, 1], exons[1:, 0]])[0]
    else:
        introns = np.array([])
    return exons, introns, left_flank, right_flank


def split_intervals(invs, pos, side='left'):
    idx = np.searchsorted(invs.ravel(), pos)
    split = np.insert(invs.ravel(), idx, [pos, pos]).reshape(-1, 2)
    split_idx = (idx + 1) // 2
    return split[:split_idx], split[split_idx:]


def parse_cds_utr_introns_flank(record, flanksize):
    exons, introns, left_flank, right_flank = parse_exons_introns_flank(record, flanksize)
    cds_start = int(record[6])
    cds_end = int(record[7])
    utr1, cds = split_intervals(exons, cds_start)
    cds, utr2 = split_intervals(cds, cds_end)
    return utr1, cds, utr2, introns, left_flank, right_flank, exons

def parse_features(record, flanksize=200):
    features = {}
    invs = {}
    features['chrom'] = record[0].replace('Chr', '')
    features['strand'] = record[5]
    utr1, invs['cds'], utr2, invs['introns'], left_flank, right_flank, invs['exons'] = parse_cds_utr_introns_flank(record, flanksize)
    if features['strand'] == '+':
        invs['5utr'] = utr1
        invs['3utr'] = utr2
        invs['upstream'] = left_flank
        invs['downstream'] = right_flank
    else:
        invs['5utr'] = utr2
        invs['3utr'] = utr1
        invs['upstream'] = right_flank
        invs['downstream'] = left_flank
    features['invs'] = invs
    return features



def get_lengths_for_norm():
    feat_lengths = Counter()
    with open(genes_bed) as bed:
        for record in bed:
            record = parse_features(record.split())
            if record['chrom'] not in unique_chrs:
                continue
            for feat_type, invs in record['invs'].items():
                for inv in invs:
                    feat_lengths[feat_type] += (inv[1] - inv[0])
    return pd.Series(feat_lengths) / 1000

def count_mismatches_in_features(der_fn, use_strand=True):
    feature_counts = Counter()
    feat_lengths = get_lengths_for_norm()
    n_records = 0
    with open(genes_bed) as bed, pysam.TabixFile(der_fn, index = searchbed_idx) as tabix:
        for record in bed:
            record = parse_features(record.split())
            if record['chrom'] not in unique_chrs:
                continue
            n_records += 1
            for feat_type, invs in record['invs'].items():
                for inv in invs:
                    if not inv[0] == inv[1]:
                        for mm in tabix.fetch(record['chrom'], *inv):
                            mm = mm.split('\t')
                            if use_strand:
                                if mm[5] == record['strand']:
                                    feature_counts[feat_type] += 1
                            else:
                                feature_counts[feat_type] += 1
    feature_counts = pd.Series(feature_counts) / feat_lengths
    return feature_counts, n_records

mod_feature_counts, n_records_der = count_mismatches_in_features(
    args.searchbed # 'results/difmod_abs_coordinates.bed.gz'
)
# ============================================================================ #
# plot

fig, ax = plt.subplots(figsize=(3, 3))
count_by_name = mod_feature_counts.reset_index()

df_counts_features = pd.DataFrame(
   dict(
      names=count_by_name['index'],
      counts=count_by_name[0]
   )
)

df_counts_features_main = df_counts_features[df_counts_features.names.isin(['upstream','5utr', 'cds', '3utr', 'downstream'])]
df_counts_features_main["names"] = pd.Categorical(df_counts_features_main["names"], categories = ['upstream','5utr', 'cds', '3utr', 'downstream'])
df_counts_features_main_sorted = df_counts_features_main.sort_values(by = "names")

new_names = {
    "upstream": "Upstream",
    "5utr": "5'UTR",
    "cds": "CDS",
    "3utr": "3'UTR",
    "downstream": "Downstream"
}
df_counts_features_main_sorted["names"] = df_counts_features_main_sorted["names"].map(new_names)

plt.bar('names', 'counts', data=df_counts_features_main_sorted, color ='steelblue', width = 0.5)
plt.setp(ax.get_xticklabels(), rotation=35, ha='right')
#ax.set_xticklabels(['Upstream', '5\'UTR', 'CDS', '3\'UTR', 'Downstream'])
ax.set_ylabel('Error sites per kb')
plt.tight_layout()
plt.savefig(args.outname, format = "pdf")
