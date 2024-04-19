#!/bin/bash
set -e -o pipefail

#######################################
# load environment variables
source ./../../env.sh
#######################################
# set project-specific variables
batch_id="Batch_Example" # The batch_id used for the processing batch. Default = "Batch_Example".
master_sample_table="Master_Sample_Table.$batch_id.txt" # The master sample table for the processing batch. Default = "Master_Sample_Table.Batch_polyA.txt".
############################################################
# Normally no need to change the following settings.
wkdir=${PWD}
outname=NanoTrans_Report_${batch_id}.html
##########################################################
r_dir=$(which R)
r_dir=${r_dir%/*}
source $miniconda3_dir/activate  $build_dir/quarto_conda_env
export PATH=${r_dir}:$PATH
R_LIBS=$NANOTRANS_HOME/build/R_libs
echo "Generating HTML report .."
quarto render NanoTrans_Report.qmd -P "wkdir:${wkdir%/*}" -P "dataset:${batch_id}" --output ${outname}
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