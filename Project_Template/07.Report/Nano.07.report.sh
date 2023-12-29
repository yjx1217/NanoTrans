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
##########################################################

source $miniconda3_dir/activate  $build_dir/quarto_conda_env
R_LIBS=$NANOTRANS_HOME/build/R_libs
quarto render NanoTrans_Report.qmd -P "wkdir:${wkdir%/*}" -P "dataset:${batch_id}"
source $miniconda3_dir/deactivate
