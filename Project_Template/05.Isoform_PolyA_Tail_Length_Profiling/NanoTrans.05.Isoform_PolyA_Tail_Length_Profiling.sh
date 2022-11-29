#!/bin/bash
set -e -o pipefail

#######################################
# load environment variables
source ./../../env.sh

#######################################
# set project-specific variables
batch_id="Batch_Example" # The batch_id used for the processing batch. Default = "Batch_Example".
master_sample_table="Master_Sample_Table.$batch_id.txt" # The master sample table for the processing batch. Default = "Master_Sample_Table.$batch_id.txt".
threads=12 # The number of threads to use. Default = 12.
debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
############################################################
# Normally no need to change the following settings.
transcript2gene_map="./../00.Reference_Genome/ref.transcript2gene_map.txt" # The reference transcript_id to gene_id and gene_name  mapping file. Default = "./../00.Reference_Genome/ref.transcript2gene_map.txt".
long_reads_dir="./../00.Long_Reads" # The directory for the Nanopore direct-RNA-seq reads. Default = "./../00.Long_Reads".
isoform_cq_dir="./../02.Isoform_Clustering_and_Quantification" # The directory for isoform clustering and quantification. Default = "./../02.Isoform_Clustering_and_Quantification".
min_mapping_quality="1"; # The minimal mapping quality to use for filtering long-read mapping alignment. Default = 1.
method="nanopolish" # The computational method for estimating polyA length: "nanopolish". Default = "nanopolish".

##########################################################

perl $NANOTRANS_HOME/scripts/batch_polya_tail_length_profiling.pl \
    -sample_table $master_sample_table \
    -threads $threads \
    -long_reads_dir $long_reads_dir \
    -isoform_cq_dir $isoform_cq_dir \
    -mapping_dir $mapping_dir \
    -method $method \
    -transcript2gene_map $transcript2gene_map \
    -batch $batch_id \
    -debug $debug

perl $NANOTRANS_HOME/scripts/batch_plot_polya_tail_length_results.pl \
    -sample_table $master_sample_table \
    -threads $threads \
    -batch $batch_id \
    -debug $debug


# clean up intermediate files
if [[ $debug == "no" ]]
then
    echo ""
fi


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
