#!/bin/bash
set -e -o pipefail

#######################################
# load environment variables 
source ./../../env.sh

#######################################
# set project-specific variables

ref_genome_download_URL="http://ftp.ensemblgenomes.org/pub/plants/release-54/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz" # Reference genome assembly FASTA file download URL. We recommend using the genome assembly FASTA file from Esembl, Ensembl Metazoa, Ensembl Plants, or Ensembl Fungi. Default = "http://ftp.ensemblgenomes.org/pub/plants/release-54/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz" for arabidopsis. Example for human: "http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz".
ref_annotation_download_URL="http://ftp.ensemblgenomes.org/pub/plants/release-54/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.54.gtf.gz" # Reference genome annotation GTF file download URL. We recommend using the genome annotation GTF file from Esembl, Ensembl Metazoa, Ensembl Plants, or Ensembl Fungi. Default = "http://ftp.ensemblgenomes.org/pub/plants/release-54/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.54.gtf.gz" for arabidopsis. Example for human: "http://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/Homo_sapiens.GRCh38.107.gtf.gz".    

if [[ ! -e ref.genome.fa ]]
then
    if [[ $ref_genome_download_URL =~ \.gz$ ]]
    then
	wget -c --no-check-certificate $ref_genome_download_URL -O  ref.genome.raw.fa.gz
	gunzip -c ref.genome.raw.fa.gz > ref.genome.raw.fa
    else
	wget -c --no-check-certificate $ref_genome_download_URL -O  ref.genome.raw.fa
    fi
fi

if [[ ! -e ref.genome.raw.gtf ]]
then
    if [[ $ref_annotation_download_URL =~ \.gz$ ]]
    then
	wget -c --no-check-certificate $ref_annotation_download_URL -O ref.genome.raw.gtf.gz
	gunzip -c ref.genome.raw.gtf.gz > ref.genome.raw.gtf
    else
	wget -c --no-check-certificate $ref_annotation_download_URL -O ref.genome.raw.gtf
    fi
fi

if [[ ! -d tmp ]]
then
    mkdir tmp
fi

perl $NANOTRANS_HOME/scripts/tidy_fasta.pl -i ref.genome.raw.fa -o ref.genome.fa
perl $NANOTRANS_HOME/scripts/tidy_id_in_ensembl_gtf.pl -i ref.genome.raw.gtf -o ref.genome.gtf
perl $NANOTRANS_HOME/scripts/transcript2gene_map_by_ensembl_gtf.pl -i ref.genome.gtf -o ref.transcript2gene_map.txt -ignore_version_number yes

if [[ "$debug" == "no" ]]
then
     rm ref.genome.raw.fa
     rm ref.genome.raw.gtf
fi


# get cDNA sequence
$gffread_dir/gffread ref.genome.gtf -g ref.genome.fa -w ref.transcriptome.fa

# index genome
$samtools_dir/samtools faidx ref.genome.fa
# $minimap2_dir/minimap2 -d ref.genome.mmi ref.genome.fa
$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -jar $picard_dir/picard.jar CreateSequenceDictionary -R ref.genome.fa -O ref.genome.dict

# index transcriptome
$samtools_dir/samtools faidx ref.transcriptome.fa
# $minimap2_dir/minimap2 -d ref.transcriptome.mmi ref.transcriptome.fa
$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -jar $picard_dir/picard.jar CreateSequenceDictionary -R ref.transcriptome.fa -O ref.transcriptome.dict

# prepare reference genome and annotation files for JAFFAL
perl $NANOTRANS_HOME/scripts/prepare_ref_genome_for_JAFFAL.pl -f ref.genome.fa -g ref.genome.gtf -p genome_annotation

if [[ -d tmp ]]
then
    rm -r tmp
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
