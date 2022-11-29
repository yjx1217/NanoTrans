#!/bin/bash
set -e -o pipefail

#######################################
# load environment variables
source ./../../env.sh

#######################################
# set project-specific variables
batch_id="Batch_Example" # The batch_id used for the processing batch. Default = "Batch_Example".
master_sample_table="Master_Sample_Table.$batch_id.txt" # The master sample table for the processing batch. Default = "Master_Sample_Table.Batch_polyA.txt".
for_human="no" # Whether the sequenced species is human: "yes" or "no". Default = "no".
threads=4 # The number of threads to use. Default = 4.
debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
############################################################
# Normally no need to change the following settings.
ref_dir="./../00.Reference_Genome" # The directory for FASTA-formatted reference genome. Default = "./../00.Reference_Genome".
long_reads_dir="./../00.Long_Reads" # The directory for FASTQ-formatted Nanopore direct-RNA-seq reads. Default = "./../00.Long_Reads".
transcript2gene_map="./../00.Reference_Genome/ref.transcript2gene_map.txt" # The transcript to gene mapping file. Default = "./../00.Reference_Genome/ref.transcript2gene_map.txt".
ref_gtf_file="./../00.Reference_Genome/ref.genome.gtf" # The ref.genome.gtf file. Default = "./../00.Reference_Genome/ref.genome.gtf".
min_mapping_quality="1"; # The minimal mapping quality to use for filtering long-read mapping alignment. Default = 1.


##########################################################

perl $NANOTRANS_HOME/scripts/batch_gene_fusion_detection.pl \
    -sample_table $master_sample_table \
    -threads $threads \
    -for_human $for_human \
    -ref_dir $ref_dir \
    -long_reads_dir $long_reads_dir \
    -transcript2gene_map $transcript2gene_map \
    -batch $batch_id \
    -debug $debug

perl $NANOTRANS_HOME/scripts/batch_plot_gene_fusion_results.pl \
    -sample_table $master_sample_table \
    -threads $threads \
    -gtf $ref_gtf_file \
    -batch $batch_id 

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
