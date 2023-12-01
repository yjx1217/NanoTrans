#!/bin/bash
set -e -o pipefail

#######################################
# load environment variables
source ./../../env.sh
PATH=$ucsc_dir:$PATH

#######################################
# set project-specific variables
batch_id="Batch_Example" # The batch_id used for the processing batch. Default = "Batch_Example".
master_sample_table="Master_Sample_Table.$batch_id.txt" # The master sample table for the processing batch. Default = "Master_Sample_Table.$batch_id.txt".
#threads=4 # The number of threads to use. Default = 4.
#top_n=20 # The number of genes with top differential RNA modification for plotting. Default = 20.
cutoff_p=0.05 # The threshold of p-value to filter diffmod kmers. Default = 0.05.
cutoff_modrate=0.5 # The threshold of the lowest modified rate change between the two condition. Default = 0.5.
#debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
############################################################
# Normally no need to change the following settings.
# Inputs
        ref_tx_bed="./../02.Isoform_Clustering_and_Quantification/${batch_id}/all_samples_combined/${batch_id}.all_samples_combined.flair_all_collapsed.isoforms.with_productivity.bed"
         ref_tx_fa="./../02.Isoform_Clustering_and_Quantification/${batch_id}/all_samples_combined/${batch_id}.all_samples_combined.flair_all_collapsed.isoforms.fa"
      raw_diffmods="./${batch_id}/all_samples_combined/${batch_id}.rna_modification.majority_direction_kmer_diffmod.table.tidy.txt"
# Outputs
 filtered_diffmods="./${batch_id}/all_samples_combined/${batch_id}.rna_modification.majority_direction_kmer_diffmod.table.tidy.filtered.txt"
      abs_diffmods="./${batch_id}/all_samples_combined/${batch_id}.rna_modification.majority_direction_kmer_diffmod.table.tidy.filtered.abs.txt"
diffmods_flank_bed="./${batch_id}/all_samples_combined/${batch_id}.rna_modification.majority_direction_kmer_diffmod.table.tidy.filtered.flank.bed"
plot_seqlogo_kmers="./${batch_id}/all_samples_combined/plot_seqlogo_kmers.pdf"
plot_topfreq_kmers="./${batch_id}/all_samples_combined/plot_frequencies_kmers.pdf"
plot_metaStopCodon="./${batch_id}/all_samples_combined/plot_sites_stop_codon_metaplot.pdf"
 plot_genomicFeats="./${batch_id}/all_samples_combined/plot_sites_genomic_feature.pdf"
  plot_overlap_pas="./${batch_id}/all_samples_combined/plot_overlap_mod_pas.pdf"
##########################################################

echo "Filtering differentially modifed kmers..."
Rscript --vanilla $NANOTRANS_HOME/scripts/tidy_filtering_kmers.R \
    -i $raw_diffmods \
    -f higher \
    -p ${cutoff_p} -l 0 -d ${cutoff_modrate} \
    -o $filtered_diffmods

echo "Plotting seqlogo and frequencies of kmers..."
Rscript --vanilla $NANOTRANS_HOME/scripts/plot_seqlogo_freq_kmers.R \
    -i  ${filtered_diffmods}  \
    -s ${plot_seqlogo_kmers} \
    -f ${plot_topfreq_kmers} 

# prep | transform relative coordinates in diffmod table to absolute coordinates
echo "Transforming relative coordinates of kmers to absolute..."
source $miniconda3_dir/activate $build_dir/miniconda3
python $NANOTRANS_HOME/scripts/tidy_transform_xpore_difmod_coords_to_abs.py \
    -r ${ref_tx_bed} \
    -s ${filtered_diffmods} \
    -o ${abs_diffmods} 

# prep | sort and index bed to plot freq of motifs around stop codon
cat ${abs_diffmods} | sort -k1,1 -k2,2n | bgzip > ${abs_diffmods}.gz
tabix -b 2 -e 3 -p bed ${abs_diffmods}.gz

# plot | freq of sites around stop codons
echo "Plotting sites distribution around stop codons..."
python $NANOTRANS_HOME/scripts/plot_sites_around_stopCodon.py \
    -r ${ref_tx_bed} \
    -s ${abs_diffmods}.gz \
    -o ${plot_metaStopCodon}

# plot | distribution of sites on genomic features
echo "Plotting sites distribution on distinct genomic features..."
python $NANOTRANS_HOME/scripts/plot_sites_genomicFeature.py \
    -r ${ref_tx_bed} \
    -s ${abs_diffmods}.gz \
    -o ${plot_genomicFeats}

# prep input bed
awk -F "\t" 'NR > 1 {print $1"_"$2"\t"$4-50"\t"$4+51}' ${filtered_diffmods} > ${diffmods_flank_bed}
# plot | overlap between mod sites & polyA sites
echo "Plotting overlaps between modification & polyA sites..."
python $NANOTRANS_HOME/scripts/plot_overlaps_modification_polyA.py \
    -b ${diffmods_flank_bed} \
    -r ${ref_tx_fa} \
    -o ${plot_overlap_pas}

source $miniconda3_dir/deactivate


############################
# checking bash exit status
if [[ $? -eq 0 ]]
then
    echo ""
    echo "#########################################################################"
    echo ""
    echo "NanoTrans message: This bash script has been successfully processed! :)"
    echo ""
    echo "#########################################################################"
    echo ""
    exit 0
fi
############################
