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
threads=4 # The number of threads to use. Default = 4.
log2foldchange_cutoff=1 # The log2 foldchange cutoff for the differential expression volcano plot. Default = 1.
adj_p_value_cutoff=0.05 # The adjusted p-value cutoff for the differential expression volcano plot. Default = 0.05.
debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
############################################################
# Normally no need to change the following settings.
transcript2gene_map="./../00.Reference_Genome/ref.transcript2gene_map.txt" # The transcript to gene mapping file. Default = "./../00.Reference_Genome/ref.transcript2gene_map.txt".
ref_dir="./../00.Reference_Genome" # The directory for the reference genome. Default = "./../00.Reference_Genome".
long_reads_dir="./../00.Long_Reads" # The directory for the Nanopore direct-RNA-seq reads. Default = "./../00.Long_Reads".
mapping_dir="./../01.Reference_Genome_based_Read_Mapping" # The directory for reference-based long-read mapping. Default = "./../01.Reference_Genome_based_Read_Mapping". 
isoform_cq_dir="./../02.Isoform_Clustering_and_Quantification" # The directory for isoform clustering and quantification. Default = "./../02.Isoform_Clustering_and_Quantification". 
min_mapping_quality="1"; # The minimal mapping quality to use for filtering long-read mapping alignment. Default = 1.
##########################################################

R_LIBS=""
source $miniconda3_dir/activate $build_dir/flair_conda_env
perl $NANOTRANS_HOME/scripts/batch_isoform_expression_and_splicing_comparison.pl \
    -sample_table $master_sample_table \
    -threads $threads \
    -isoform_cq_dir $isoform_cq_dir \
    -batch $batch_id \
    -transcript2gene_map $transcript2gene_map \
    -debug $debug
source $miniconda3_dir/deactivate

source ./../../env.sh
perl $NANOTRANS_HOME/scripts/batch_plot_diffexp_results.pl \
    -sample_table $master_sample_table \
    -threads $threads \
    -batch $batch_id \
    -log2foldchange_cutoff $log2foldchange_cutoff \
    -adj_p_value_cutoff $adj_p_value_cutoff \
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
