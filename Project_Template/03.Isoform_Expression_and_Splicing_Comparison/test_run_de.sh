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
contrast="vir1,VIRc" # format: test_group_tag,control_group_tag. Must be the same as the name in the comparison_group column of the master_sample_table.
read_counts_cutoff=5 # The read count cutoff for the differentially expressed genes/isoforms table. Default = 5.
log2foldchange_cutoff=1 # The log2 foldchange cutoff for the differential expression volcano plot. Default = 1.
adj_p_value_cutoff=0.05 # The adjusted p-value cutoff for the differential expression volcano plot. Default = 0.05.
debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
############################################################
# Normally no need to change the following settings.
quant_table="./../02.Isoform_Clustering_and_Quantification/${batch_id}/all_samples_combined/${batch_id}.all_samples_combined.counts_matrix.tidy.txt"  # The directory for isoform clustering and quantification. Default = "./../02.Isoform_Clustering_and_Quantification". 
##########################################################

echo "Detecting genes/isoforms with differential expression ... \n"
Rscript --vanilla $NANOTRANS_HOME/scripts/differential_expression.R \
    --tidy_count_table ${quant_table} \
    --read_counts_cutoff ${read_counts_cutoff} \
    --log2foldchange_cutoff ${log2foldchange_cutoff} \
    --adj_p_value_cutoff ${adj_p_value_cutoff} \
    --contrast ${contrast} \
    --batch_id ${batch_id} 

source activate ../../build/flair_conda_env/
Rscript --vanilla $NANOTRANS_HOME/scripts/differential_isoform_usage.R \
    --tidy_count_table ${quant_table} \
    #--adj_p_value_cutoff ${adj_p_value_cutoff} \
    --contrast ${contrast} \
    --batch_id ${batch_id} 
source deactivate

    