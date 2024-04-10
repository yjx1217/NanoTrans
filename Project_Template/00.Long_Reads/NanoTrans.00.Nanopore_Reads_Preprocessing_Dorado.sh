#!/bin/bash
set -e -o pipefail
#######################################
# load environment variables
source ./../../env.sh

########################
dorado_run_device="cpu" # The running device of basecalling: "cuda:0,...,N", "cuda:all", or "cpu". Default = "cpu".
if [[ $dorado_run_device != "cpu" ]]
then
    # set CUDA environment for GPU
    gpu_bin_path="/public/software/cuda-11.4/bin" # The path of CUDA's bin directory. Default = "/public/software/cuda-11.4/bin".
    gpu_lib_path="/public/software/cuda-11.4/lib64" # The path of CUDA's lib directory. Default = "/public/software/cuda-11.4/lib64".
    gpu_include_path="/public/software/cuda-11.4/include" # The path of CUDA's include directory. Default = "/public/software/cuda-11.4/include".
fi
##########################3
export PATH=$gpu_bin_path:$PATH
export LD_LIBRARY_PATH=$gpu_lib_path:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=$gpu_include_path:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$C_INCLUDE_PATH

#######################################
# set project-specific variables
sample_id="vir1.rep1" # The flowcell ID of the nanopore run. Default = "vir1.rep1".
raw_fast5_dir="./raw_fast5/$sample_id" # The directory containing the raw nanopore reads, fast5 or pod5 format, before basecalling. Default = "./raw_fast5/$sample_id".
decode_models_name="rna002_70bps_hac@v3" # The RNA model name used for basecalling. For more details, please to see: https://github.com/nanoporetech/dorado?tab=readme-ov-file#rna-models. Default = "rna002_70bps_hac@v3".

#############################
# Normally no need to change the following parameter settings.
threads=8 # The number of CPU threads to use. Default = 8.
run_basecalling="yes"   # Whether to perform basecalling:    "yes" or "no". Default = "yes". 
run_nanoplotting="yes"  # Whether to perform nanoplotting:   "yes" or "no". Default = "yes". 
run_demultiplexing="no" # Whether to perform demultiplexing: "yes" or "no". Default = "no". 
barcode_kit_version="" # The barcode kit version of the nanopore run. Default = "".
basecalled_fastq_dir="./basecalled_fastq/$sample_id" # The directory containing the basecalled fastq reads. This directory will be automatically generated when running basecalling.
basecalled_summary_dir="./basecalled_summary/$sample_id" # The directory containing nanoplot summary outputs. This directory will be automatically generated when running nanoplot.
qual=5 # read quality filter for basecalling. Default = 5.
demultiplexing_threads=$threads # The number of threads to use for demultiplexing. Default = 8.
demultiplexed_fastq_dir="$basecalled_fastq_dir/demultiplexed_fastq" # The directory containing the demultiplexed basecalled nanopore reads. This directory will be automatically generated when running demultiplexing. 
##############################

wkdir=$(pwd)
if [[ "$run_demultiplexing" == "yes" && "barcode_kit_version" == "" ]]
then
    echo "The variable run_demultiplexing has been set to \"yes\" but the barcode_kit_version variable has not been specified!"
    echo "Please specified barcode_kit_version if you want to run demultiplexing."
    exit;
fi

if [[ "$run_basecalling" == "yes" ]]
then
    echo "Check if $basecalled_fastq_dir is empty for running basecalling."

    if [[ -d $basecalled_fastq_dir && "$(ls $basecalled_fastq_dir)" ]]
    then
        echo "Warning! The basecalled fastq directory $basecalled_fastq_dir exists and it is not empty! Please empty its content if you want to run basecalling."
        echo "Exit!!!"
        exit
    else
        echo "Check passed!"
        echo "Running basecalling .."
        if [[ ! -d $basecalled_fastq_dir ]]
        then
            mkdir -p $basecalled_fastq_dir
        fi

        $dorado_dir/dorado download --model ${decode_models_name}
        $dorado_dir/dorado basecaller \
            --device $dorado_run_device \
            --min-qscore $qual \
            --emit-fastq \
            --recursive \
            --verbose \
            ${decode_models_name} \
            $raw_fast5_dir | gzip -c > ${basecalled_fastq_dir}/${sample_id}.basecalled_reads.Q${qual}.pass.fastq.gz
    fi
fi

cd $wkdir
if [[ "$run_demultiplexing" == "yes" ]]
then
    echo "Check if $basecalled_fastq_dir has basecalled reads for running demultiplexing."
    if [[ "$(ls $basecalled_fastq_dir)" ]]
    then
        echo "Running demultiplexing."
        $dorado_dir/dorado demux \
            --kit-name $barcode_kit_version \
            --emit-fastq \
            --output-dir $demultiplexed_fastq_dir \
            --threads    $demultiplexing_threads \
            ${basecalled_fastq_dir}/${sample_id}.basecalled_reads.Q${qual}.pass.fastq.gz

        cd $demultiplexed_fastq_dir
        for b in ${barcode_kit_version}_barcode*.fastq
        do 
            barcode_id=${b%.*}
            echo "for demultiplexing: barcode=$barcode_id"
            cat ./$b | gzip -c > ${sample_id}.basecalled_reads.Q${qual}.pass.${b}.gz
        done
        cat ./unclassified.fastq | gzip -c > ${sample_id}.basecalled_reads.Q${qual}.pass.unclassified.fastq.gz
    else
        echo "There is no reads in $basecalled_fastq_dir!"
        echo "Please put the basecalled reads in $basecalled_fastq_dir for demultiplexing!"
        echo "Exit!!!"
        exit
    fi
fi

set +oe pipefail 

cd $wkdir
if [[ "$run_nanoplotting" == "yes" ]]
then
    echo "Check if $basecalled_fastq_dir has basecalled reads for running nanoplotting."
    if [[ "$(ls $basecalled_fastq_dir)" ]]
    then
        echo "Running nanoplotting."
        mkdir -p $basecalled_summary_dir
        fastq_input="$basecalled_fastq_dir/$sample_id.basecalled_reads.Q${qual}.pass.fastq.gz"
        $nanoplot_dir/NanoPlot \
            --threads $threads \
            --fastq $fastq_input \
            --N50 \
            -o "$basecalled_summary_dir/${sample_id}_basecalled_reads_Q${qual}_pass_NanoPlot_out"
    fi
    if [[ "$run_demultiplexing" == "yes" ]]
    then
        cd $demultiplexed_fastq_dir
        for b in ${barcode_kit_version}_barcode*.fastq
        do
            barcode_id=${b%.*}
            echo "for nanoplotting: barcode=$barcode_id"
            fastq_input="./$sample_id.basecalled_reads.Q${qual}.pass.${b}.gz"
            $nanoplot_dir/NanoPlot \
            --threads $threads \
            --fastq $fastq_input \
            --N50 \
            -o "./../../../$basecalled_summary_dir/${sample_id}_basecalled_reads_Q${qual}_pass_${barcode_id}_NanoPlot_out"
        done
        echo "for nanoplotting: unclassified"
        fastq_input="./$sample_id.basecalled_reads.Q${qual}.pass.unclassified.fastq.gz"
        $nanoplot_dir/NanoPlot \
                --threads $threads \
                --fastq $fastq_input \
                --N50 \
                -o "./../../../$basecalled_summary_dir/${sample_id}_basecalled_reads_Q${qual}_pass_unclassified_NanoPlot_out" 
    fi
fi



############################
# checking bash exit status
if [[ $? -eq 0 ]]
then
    echo "" 
    echo "##########################################################################"
    echo "" 
    echo "NanoTrans message: This bash script has been successfully processed! :)"
    echo ""
    echo "##########################################################################"
    echo ""
    exit 0
fi
############################
