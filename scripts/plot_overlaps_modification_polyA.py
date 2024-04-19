# ============================================================================ #
# info

# date: 231008
# desc: overlaps between RNA modification & polyA sites
# refs: "https://github.com/bartongroup/Simpson_Barton_Nanopore_1/blob/master/notebooks/18_pas_strength_at_m6a_sites.ipynb"

# ============================================================================ #
# load

import argparse
import itertools as it
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt

import pysam
from ushuffle import shuffle as ushuffle

plt.rcParams['pdf.fonttype'] = 42
#plt.rcParams['font.family'] = 'Arial'
# ============================================================================ #
# parse args

parser = argparse.ArgumentParser(description='Plot overlaps between motifs of m6A & polyA')
parser.add_argument('-b', '--kmerbed',  help='kmers bed, required', required=True)
parser.add_argument('-r', '--reffasta', help='fasta file for searching flank sequences of motifs, required', required=True)
parser.add_argument('-o', '--outname',  help='output directory, required', default="./plot_overlap_mod_pas.pdf",  required=True)

args = parser.parse_args()

# ============================================================================ #
# functions

def main(kmerbed, reffasta, outname):
    print('the kmerbed is', kmerbed)
    print('the reference fasta is', reffasta)
    print('the outname is', outname)

RC = str.maketrans('ACGTSWRYN', 'TGCASWYRN')

def rev_comp(seq):
    return seq.translate(RC)[::-1]


def motif_count(seqs, motifs, motif_len):
    kmer_pos_count = defaultdict(lambda: np.zeros(shape=len(seqs[0])))
    for seq in seqs:
        assert len(seq) == len(seqs[0]), 'All seqs should be uniform length'
        for i in range(len(seq) - motif_len + 1):
            kmer = seq[i: i + motif_len]
            if kmer in motifs:
                kmer_pos_count[kmer][i: i + motif_len] += 1
    return np.array([kmer_pos_count[k] for k in motifs])


def permute_seqs(seqs, w):
    shuf_seqs = []
    for seq in seqs:
        shuf = ushuffle(seq.encode(), w).decode()
        shuf_seqs.append(shuf)
    return shuf_seqs

def motif_enrichment(seqs, w, n_perm, motifs, motif_len):
    obs = motif_count(seqs, motifs, motif_len)
    exp = []
    for _ in range(n_perm):
        shuf_seqs = permute_seqs(seqs, w)
        exp.append(motif_count(shuf_seqs, motifs, motif_len))
    exp = np.array(exp)
    enrichment = np.log2(obs.sum(0) + 0.5) - np.log2(exp.sum(1) + 0.5)
    return enrichment

# ============================================================================ #
# m6a

chars = [[ "G", "T", "A"], ["G", "A"], ["A"], ["C"], ["C", "T", "A"]]
combos = list(it.product(*chars))

M6A_MOTIFS = []
for r in combos:
    joined = ''.join(r)
    M6A_MOTIFS.append(joined)

# ============================================================================ #
# polyA

CANON_PAS = 'AATAAA'
PAS_KMERS = set([CANON_PAS])
for i, n in it.product(range(6), 'ACGT'):
    kmer = CANON_PAS[:i] + n + CANON_PAS[i + 1:]
    if not kmer == 'AAAAAA':
        PAS_KMERS.add(kmer)
PAS_KMERS = list(PAS_KMERS)

def get_pas_enrichment(der_sites_fn, fasta_fn, w, n_perm):
    seqs = []
    with pysam.FastaFile(fasta_fn) as fasta, open(der_sites_fn) as gtf:
        for record in gtf:
            record = record.split()
            #pos = ((int(record[3]) - 1) + int(record[4])) // 2
            chrom = record[0]
            start = int(record[1])
            end   = int(record[2])
            if start < 0:
                start = 0
            #strand = record[6]
            try:
                seq = fasta.fetch(chrom, start, end)
            except Exception as e:
                print(e)
                continue
            #if record[3] == 'start':
            seq = "{:N>101}".format(seq)  
            #if strand == '-':
            #    seq = rev_comp(seq)
            seqs.append(seq)
    seqs = list(set(seqs))
    pas_enrichment = {}
    for motif in PAS_KMERS:
        pas_enrichment[motif] = motif_enrichment(seqs, w, n_perm, [motif,], 6)
    m6a_enrichment = motif_enrichment(seqs, w, n_perm, M6A_MOTIFS, 5)
    return pas_enrichment, m6a_enrichment, seqs

# ============================================================================ #
# run

pas_enrichment, m6a_enrichment, seqs = get_pas_enrichment(
    args.kmerbed,
    args.reffasta,
    w=2, n_perm=100
)

# ============================================================================ #
# plot

pal = ['#0072b2', '#d55e00', '#009e73', '#f0e442', '#cc79a7']
fig, ax = plt.subplots(figsize=(5, 3))

for k, v in pas_enrichment.items():
    ax.plot(np.arange(-51, 50), v.mean(0), color='#cccccc', zorder=-1)
ax.plot(np.arange(-51, 50), pas_enrichment['AATAAA'].mean(0), label='AATAAA',    color=pal[0],    lw=2, zorder=5)
ax.plot(np.arange(-51, 50), pas_enrichment['TATAAA'].mean(0), label='TATAAA',    color=pal[1],    lw=2, zorder=4)
ax.plot(np.arange(-51, 50), pas_enrichment['AACAAA'].mean(0), label='AACAAA',    color=pal[2],    lw=2, zorder=2)
ax.plot(np.arange(-51, 50), pas_enrichment['AAGAAA'].mean(0), label='AAGAAA',    color=pal[3],    lw=2, zorder=1)
ax.plot(np.arange(-51, 50), m6a_enrichment.mean(0),           label='m6A motif', color='#252525', lw=2, zorder=0)
ax.legend()
ax.set_xlabel('Distance from m6A motif (nt)')
ax.set_ylabel('Motif enrichment\nover shuffled seqs (log2 scale)')
ax.set_xlim(-50, 50)
plt.tight_layout()
plt.savefig(args.outname, format = "pdf")

# ============================================================================ #
# test

if __name__ == '__main__':
    try:
        main(args.kmerbed, args.reffasta, args.outname)
    except Exception as e:
        print(e)
