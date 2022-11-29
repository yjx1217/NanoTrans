#!/bin/bash
set -e -o pipefail

##########################################
# load environment variables
source ./../../env.sh

###########################################
# set project-specific variables
batch_id="Batch_Example" # The batch_id used for the processing batch. Default = "Batch_Example".
master_sample_table="Master_Sample_Table.$batch_id.txt" # The master sample table for the processing batch. Default = "Master_Sample_Table.$batch_id.txt".
threads=4 # The number of threads to use. Defualt = 4.
debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
###########################################

###########################################
# Normally no need to change the following settings.
mapping_mode="map2genome" # The mapping mode for the processing batch: "map2transcriptome" or "map2genome". Default = "map2transcriptome".
long_read_mapper="minimap2" # The long-read mapper used for mapping. Default = "minimap2".
long_read_technology="nanopore_drs" # The long-read technology used: "nanopore_drs". Default = "nanopore_drs".
ref_genome_dir="./../00.Reference_Genome" # The directory for the reference genome. Default = "./../00.Reference_Genome".
long_reads_dir="./../00.Long_Reads" # The directory for the Nanopore direct-RNA-seq reads. Default = "./../00.Long_Reads".
min_mapping_quality=0 # The minimal mapping quality to use for filtering long-read mapping alignment. Default = 0.
###########################################

# process the pipeline

echo ""
echo "Running long-read mapping for $long_read_technology reads with $long_read_mapper ..."
echo ""

test_file_existence () {
    filename=$1
    if [[ ! -e $filename ]]
    then
        echo "The file/directory $filename does not exists! Process terminated!"
        exit
    else
        echo "Test passed!"
    fi
}

echo "Testing the existence of ref_genome_dir: $ref_genome_dir .."
test_file_existence $ref_genome_dir
echo "Testing the existence of long_reads_dir: $long_reads_dir .."
test_file_existence $long_reads_dir

if [[ $mapping_mode == "map2genome" ]]
then
    perl $NANOTRANS_HOME/scripts/batch_long_read_spliced_mapping.pl \
    -i $master_sample_table \
    -t $threads \
    -mapping_mode $mapping_mode \
    -ref_dir $ref_genome_dir \
    -long_reads_dir $long_reads_dir \
    -long_read_mapper $long_read_mapper \
    -long_read_technology $long_read_technology \
    -batch_id ${batch_id} \
    -min_mapping_quality $min_mapping_quality \
    -debug $debug
fi


echo "Long-read mapping has finished ..."
echo ""
# clean up intermediate files
echo ""
echo "Clean up intermediate files ..."
echo ""
if [[ $debug = "no" ]]
then
    :
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
