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
threads=12 # The number of threads to use. Default = 12.
top_n=20 # The number of genes with top differential RNA modification for plotting. Default = 20.
debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
############################################################
# Normally no need to change the following settings.
transcript2gene_map="./../00.Reference_Genome/ref.transcript2gene_map.txt" # The reference transcript to gene mapping file. Default = "./../00.Reference_genome/ref.transcript2gene_map.txt".
isoform_cq_dir="./../02.Isoform_Clustering_and_Quantification" # The directory for isoform clustering and quantification. Default = "./../02.Isoform_Clustering_and_Quantification".
long_reads_dir="./../00.Long_Reads" # The directory for FASTQ-formatted Nanopore direct-RNA-seq reads. Default = "./../00.Long_Reads".
##########################################################

R_LIBS=""
source $miniconda3_dir/activate $build_dir/xpore_conda_env
perl $NANOTRANS_HOME/scripts/batch_rna_modification_detection.pl \
    -batch_id $batch_id \
    -sample_table $master_sample_table \
    -threads $threads \
    -long_reads_dir $long_reads_dir \
    -isoform_cq_dir $isoform_cq_dir \
    -transcript2gene_map $transcript2gene_map \
    -debug $debug
source $miniconda3_dir/deactivate

source ./../../env.sh
perl $NANOTRANS_HOME/scripts/batch_plot_rna_modification_results.pl \
    -sample_table $master_sample_table \
    -batch $batch_id \
    -top_n $top_n \
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
