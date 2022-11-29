#!/bin/bash
set -e -o pipefail

#######################################
# load environment variables
source ./../../env.sh
PATH=$ucsc_dir:$PATH

#######################################
# set project-specific variables
# Important note: Please only run this bash script after you have successfully run NanoTrans.02.Isoform_Clustering_and_Quantification.sh.
batch_id="Batch_Example" # The batch_id used for the processing batch. Default = "Batch_Example".
query_gene_id="AT1G05850" # The query gene id for plotting isoform usage. Default = "".
min_reads_support=6 # The minimal reads support required for the isoform to be plotted. Default = 6.
debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
############################################################
# Normally no need to change the following settings.
output_prefix="./$batch_id/all_samples_combined/query_gene_isoform_usage/$batch_id.$query_gene_id.isoform_usage"
isoform_bed="./$batch_id/all_samples_combined/$batch_id.all_samples_combined.flair_all_collapsed.isoforms.with_productivity.bed"
count_matrix="./$batch_id/all_samples_combined/$batch_id.all_samples_combined.counts_matrix.tsv"
##########################################################

echo "##########################"
echo "Important note: Please only run this bash script after you have successfully run NanoTrans.02.Isoform_Clustering_and_Quantification.sh."
echo "##########################"

if [[ ! -d $batch_id/all_samples_combined/query_gene_isoform_usage ]]
then
    mkdir -p $batch_id/all_samples_combined/query_gene_isoform_usage
fi

$flair_dir/plot_isoform_usage $isoform_bed $count_matrix $query_gene_id --min_reads $min_reads_support  -o $output_prefix


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
