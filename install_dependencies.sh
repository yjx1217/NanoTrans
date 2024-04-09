#!/bin/bash
# last update: 2024/04/09
set -e -o pipefail

#########################
NANOTRANS_HOME=$(pwd)
BUILD="build"
gpu_support="yes"
mainland_china_installation="no";
#########################                                                                                                       

timestamp () {
  date +"%F %T"
}

clean () {
    dir=$1
    if [ -d $dir ] 
    then
	echo "remove previously failed installation in $BUILD/$dir"
	rm -rf $dir
    fi
}

download () {
  url=$1
  download_location=$2
  echo "Downloading $url to $download_location"
  if [[ -f $1 ]];then
    # if the tool has been downloaded and deposited in the local, then copy that to destination and install it
    echo -n "The tool exists in the local: $1"
    cp $1 $2
  else
  wget -c --no-check-certificate --max-redirect=30 $url -O $download_location
  fi
}

clone () {
  url=$1
  dir=$(basename $url)
  echo "run clone for \"git clone $url\""
  git clone $url --depth 1
  cd $dir
  git fetch --unshallow
}

tidy_version () { 
    echo "$1" | awk -F. '{ printf("%d%03d%03d%03d\n", $1,$2,$3,$4); }';
}

check_installed () {
    if [ -e "$1/installed" ]; then
        echo "installed"
    else
        echo ""
    fi
}

note_installed () {
    touch "$1/installed"
}

install_r_pkg () {
    if [ -z $(check_installed "$rlib_dir/$1") ]; then
        echo "[$(timestamp)] Installing $1 (R package)..."
        clean "$rlib_dir/$1"
        R -e "install.packages(\"$1\", repos=\"http://cran.rstudio.com/\", lib=\"$build_dir/R_libs/\")"
        note_installed "$rlib_dir/$1"
    fi  
}
install_r_pkg_source () {
    if [ -z $(check_installed "$rlib_dir/$1") ]; then
        echo "[$(timestamp)] Installing $1 (R package)..."
        clean "$rlib_dir/$1"
        url_source_rpkg="https://cran.r-project.org/src/contrib/Archive/$1/${1}_${2}.tar.gz"
        R -e "install.packages(\"$url_source_rpkg\", repos=NULL, type=\"source\", lib=\"$build_dir/R_libs/\")"
        note_installed "$rlib_dir/$1"
    fi  
}
#install_r_pke_byBiocManager () {
#    if [ -z $(check_installed "$rlib_dir/$1") ]; then
#        echo "[$(timestamp)] Installing $1 (R package)..."
#        clean "$rlib_dir/$1"
#        R -e ".libPaths(\"$build_dir/R_libs/\"); BiocManager::install(\"$1\", lib=\"$build_dir/R_libs/\")"
#        note_installed "$rlib_dir/$1"
#    fi  
#}

# aim: marking that pkg (python) has beed installed successfully through conda 
install_pkg_byconda () {
    if [ ! -d $site_packages/$1 ]; then mkdir -p $site_packages/$1 ; fi
    if [ -z $(check_installed $site_packages/$1) ]; then
        echo "[$(timestamp)] Installing $1 ..."
        source $miniconda3_dir/activate $build_dir/miniconda3
        $miniconda3_dir/conda install -y -c conda-forge -c bioconda $1=$2
        source $miniconda3_dir/deactivate 
    fi
    note_installed $site_packages/$1
}

install_and_create_pkg_condaenv () {
    if [ -z $(check_installed $build_dir/${1}_conda_env/bin) ]; then
        echo "[$(timestamp)] Installing $1 ..."
        cd $build_dir
        $miniconda3_dir/conda create -y -p $build_dir/${1}_conda_env
        source $miniconda3_dir/activate    $build_dir/${1}_conda_env
        $miniconda3_dir/conda install -y -c conda-forge -c bioconda $1=$2
        source $miniconda3_dir/deactivate
        note_installed $build_dir/${1}_conda_env/bin
    fi
}


echo ""
echo ""
echo "##################################################################"
echo "###                                                            ###"
echo "###                  Welcome to NanoTrans                      ###"
echo "###                                                            ###"
echo "##################################################################"
echo ""
echo ""
echo ""
echo ""
echo "[$(timestamp)] Performing gcc version check .." 
GCC_VERSION=$(gcc -dumpversion)
if [[ $(tidy_version "$GCC_VERSION") -ge $(tidy_version "4.9.0") ]]
then 
   echo "Detected gcc version: $GCC_VERSION."
   echo "NanoTrans requires gcc version >= 4.9."
   echo "Requirement satisfied."
   echo "Installation continues .."
else
   echo "Detected gcc version: $GCC_VERSION."
   echo "NanoTrans requires gcc version >= 4.9."
   echo "Requirement unsatisfied"
   echo "Please upgrade your gcc version. Installation Terminates!"
   echo "Exit!"
   exit 1
fi
echo ""
echo ""
echo "[$(timestamp)] Installation starts ..."
echo ""


while getopts ":hgca" opt
do
    case "${opt}" in
        h)
            echo "Usage:"
            echo "bash install_dependencies.sh"
            echo "When running installation with CUDA GPU support, please run this script with the '-g' option >"
            echo "bash install_dependencies.sh -g"
            echo "When running installation without CUDA GPU support, please run this script with the '-c' option >"
            echo "bash install_dependencies.sh -c"
            echo "When running installation within mainland China and use alternative conda repository sources, please run this script with the '-a' option >"
            echo "bash install_dependencies.sh -a";;
        g)
	    echo "Detected the '-g' option >"
	    echo "Set gpu_support as 'yes'"
	    gpu_support="yes";;
        c)
	    echo "Detected the '-c' option >"
	    echo "Set gpu_support as 'no'"
	    gpu_support="no";;
        a)
	    echo "Detected the '-a' option >"
	    echo "Set installation location as 'mainland_china' and use alternative conda repository sources"
	    mainland_china_installation="yes";;
    esac
done
echo "";

# echo "mainland_china_installation=$mainland_china_installation"

if [ -z "$MAKE_JOBS" ]
then
    echo "[$(timestamp)] Defaulting to 2 concurrent jobs when executing make. Override with MAKE_JOBS=<NUM>"
    MAKE_JOBS=2
fi

if [ ! -z "$INSTALL_DEPS" ]; then
    echo "Installing ONT build dependencies for Debian/Ubuntu."
    echo "sudo privileges are required and you will be prompted to enter your password"
    sudo apt-get update
    xargs -a debiandeps sudo apt-get install -y
fi

MINICONDA3_VERSION="py39_23.11.0-2" # released on 2023.11.16
if [[ "$mainland_china_installation" == "no" ]]
then
    MINICONDA3_DOWNLOAD_URL="https://repo.anaconda.com/miniconda/Miniconda3-${MINICONDA3_VERSION}-Linux-x86_64.sh"
else
    MINICONDA3_DOWNLOAD_URL="https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-${MINICONDA3_VERSION}-Linux-x86_64.sh"
fi


# for reads preparation and preprocessing
BLAST_VERSION="2.2.31" #
BLAST_DOWNLOAD_URL="http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz"

# RMBLAST_VERSION="2.2.28" #
# RMBLAST_DOWNLOAD_URL="http://ftp.ncbi.nlm.nih.gov/blast/executables/rmblast/${RMBLAST_VERSION}/ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz"

# SRA_VERSION="3.0.0" # released on 2022.02.10
# SRA_DOWNLOAD_URL="https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz"

TRIMMOMATIC_VERSION="0.38" #
TRIMMOMATIC_DOWNLOAD_URL="http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${TRIMMOMATIC_VERSION}.zip"

BPIPE_VERSION="0.9.9.2"
BPIPE_DOWNLOAD_URL="https://github.com/ssadedin/bpipe/releases/download/${BPIPE_VERSION}/bpipe-${BPIPE_VERSION}.tar.gz"

VELVET_VERSION="1.2.10"
VELVET_DOWNLOAD_URL="https://github.com/dzerbino/velvet/archive/v${VELVET_VERSION}.tar.gz"

OASES_VERSION="0.2.09"
OASES_DOWNLOAD_URL="https://github.com/dzerbino/oases/archive/refs/tags/${OASES_VERSION}.tar.gz"

BOWTIE2_VERSION="2.4.5"
BOWTIE2_DOWNLOAD_URL="https://github.com/BenLangmead/bowtie2/releases/download/v${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip"

FASTXTK_VERSION="0.0.13_binaries_Linux_2.6_amd64"
FASTXTK_DOWNLOAD_URL="http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_${FASTXTK_VERSION}.tar.bz2"

BBMAP_VERSION="39.00"
BBMAP_DOWNLOAD_URL="https://sourceforge.net/projects/bbmap/files/BBMap_${BBMAP_VERSION}.tar.gz"


# basecalling
GUPPY_VERSION="6.2.1" # released on 2022.06.07
#if [[ "$gpu_support" == "yes" ]]
#then
GUPPY_GPU_DOWNLOAD_URL="https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_${GUPPY_VERSION}_linux64.tar.gz"
#else
GUPPY_CPU_DOWNLOAD_URL="https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_${GUPPY_VERSION}_linux64.tar.gz"
#fi

DORADO_VERSION="0.6.0" # released on 2024.04.02
DORADO_DOWNLOAD_URL="https://cdn.oxfordnanoportal.com/software/analysis/dorado-${DORADO_VERSION}-linux-x64.tar.gz"

# QC
NANOPLOT_VERSION="1.40.0" # released on 2022.04.08
# NANOFILT_VERSION="2.8.0" # released on 2021.02.26
#NUMPY_VERSION="1.16.3"
SCIPY_VERSION="1.2.1"

NANOPOLISH_VERSION="0.14.0" # released on 2022.05.25
NANOPOLISH_GITHUB_VERSION="07cb03d"

FLAIR_VERSION="2.0.0" # 

GFFREAD_VERSION="0.12.7" # released on 2021.07.23
GFFREAD_DOWNLOAD_URL="http://ccb.jhu.edu/software/stringtie/dl/gffread-${GFFREAD_VERSION}.Linux_x86_64.tar.gz"

PICARD_VERSION="2.27.4" # released on 2022.06.30
PICARD_DOWNLOAD_URL="https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar"

# SEQTK_DOWNLOAD_URL="https://github.com/lh3/seqtk"
# SEQTK_GITHUB_COMMIT_VERSION="7c04ce7" # commited on 2021.03.23

XPORE_VERSION="2.1"
JAFFAL_VERSION="2.3"
JAFFAL_DOWNLOAD_URL="https://github.com/Oshlack/JAFFA/releases/download/version-${JAFFAL_VERSION}/JAFFA-version-${JAFFAL_VERSION}.tar.gz"

# UCSC Utilities
BEDPARTITION_DOWNLOAD_URL="http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedPartition"
GTFTOGENEPRED_DOWNLOAD_URL="http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred"
BLAT_DOWNLOAD_URL="http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat"

# pkgs used in generating HTML reports
QUARTO_VERSION="1.3.450"

# Create the $BUILD directory for dependency installation

if [ -d $BUILD ]
then
    echo ""
    echo "[$(timestamp)] Detected previously generated $BUILD directory."
else
    echo "[$(timestamp)] Create the new $BUILD directory."
    mkdir $BUILD
    echo ""
fi

cd $BUILD
build_dir=$(pwd)

echo ""
echo "[$(timestamp)] Download and install all the dependencies"
echo ""

# Downloading all the dependencies

# ---------- set Perl & Python environment variables -------------
PYTHONPATH="$build_dir"
PERL5LIB="$build_dir:$PERL5LIB"
PERL5LIB="$build_dir/cpanm/perlmods/lib/perl5:$PERL5LIB"

echo ""
echo "[$(timestamp)] Installing Perl modules ..."
cpanm_dir=$build_dir/cpanm
mkdir -p  $cpanm_dir
if [ -z $(check_installed $cpanm_dir) ]; then
    cd $cpanm_dir
    # wget -c --no-check-certificate -O - https://cpanmin.us/ > cpanm
    # work around for the unstable downloading issue
    cp $NANOTRANS_HOME/misc/cpanm .
    ########
    chmod +x cpanm
    mkdir -p perlmods
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --prompt --notest --skip-installed Test::More@1.302086
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --prompt --notest --skip-installed Env@1.04
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --prompt --skip-installed Statistics::Descriptive@3.0612
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --prompt --skip-installed Statistics::Descriptive::Discrete@0.07
    # $cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Math::Random@0.72
    # $cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Math::Round@0.07
    note_installed $cpanm_dir
fi

echo ""
echo "[$(timestamp)] Installing R libraries ..."
which R || echo "Please install R first."
r_path=$(which R)
rlib_dir="$build_dir/R_libs"
mkdir -p $rlib_dir
cd $rlib_dir
R_VERSION=$(R --version |head -1 |cut -d " " -f 3)

if [ $(tidy_version "$R_VERSION") -ge $(tidy_version "3.6.0") ]; then
    if [ -z $(check_installed "$rlib_dir/BiocManager") ]; then
        clean "$rlib_dir/BiocManager"
        echo "R_VERSION=$R_VERSION, use the new bioconductor installation protocol"
        R -e ".libPaths(\"$build_dir/R_libs/\");install.packages(\"BiocManager\", repos=\"http://cran.rstudio.com/\", lib=\"$build_dir/R_libs/\")";
        note_installed "$rlib_dir/BiocManager"
    fi
else
    die "R >= v3.6.0 is needed! Exit!"
fi

#if [ -z $(check_installed "$rlib_dir/optparse") ]; then
#    clean "$rlib_dir/optparse"
#    R -e "install.packages(\"optparse\", repos=\"http://cran.rstudio.com/\", lib=\"$build_dir/R_libs/\")"
#    note_installed "$rlib_dir/optparse"
#fi

# if [ -z $(check_installed "$rlib_dir/RColorBrewer") ]; then
#     clean "$rlib_dir/RColorBrewer"
#     R -e "install.packages(\"RColorBrewer\", repos=\"http://cran.rstudio.com/\", lib=\"$build_dir/R_libs/\")"
#     note_installed "$rlib_dir/RColorBrewer"
# fi

#if [ -z $(check_installed "$rlib_dir/ggplot2") ]; then
#    clean "$rlib_dir/ggplot2"
#    R -e "install.packages(\"ggplot2\", repos=\"http://cran.rstudio.com/\", lib=\"$build_dir/R_libs/\")"
#    note_installed "$rlib_dir/ggplot2"
#fi
#
#if [ -z $(check_installed "$rlib_dir/ggrepel") ]; then
#    clean "$rlib_dir/ggrepel"
#    R -e "install.packages(\"ggrepel\", repos=\"http://cran.rstudio.com/\", lib=\"$build_dir/R_libs/\")"
#    note_installed "$rlib_dir/ggrepel"
#fi
#
#if [ -z $(check_installed "$rlib_dir/dplyr") ]; then
#    clean "$rlib_dir/dplyr"
#    R -e "install.packages(\"dplyr\", repos=\"http://cran.rstudio.com/\", lib=\"$build_dir/R_libs/\")"
#    note_installed "$rlib_dir/dplyr"
#fi

# if [ -z $(check_installed "$rlib_dir/samplesizeCMH") ]; then
#     clean "$rlib_dir/samplesizeCMH"
#     R -e "install.packages(\"samplesizeCMH\", repos=\"http://cran.rstudio.com/\", lib=\"$build_dir/R_libs/\")"
#     note_installed "$rlib_dir/samplesizeCMH"
# fi

install_r_pkg optparse
install_r_pkg ggplot2
install_r_pkg ggrepel
install_r_pkg dplyr
install_r_pkg data.table 
install_r_pkg_source ggseqlogo "0.1"
install_r_pkg DT         
install_r_pkg tidyr      
install_r_pkg tibble     
install_r_pkg generics

# ------------- Miniconda3 --------------------
echo ""
echo "[$(timestamp)] Installing miniconda3 ..."
miniconda3_dir="$build_dir/miniconda3/bin"
if [ -z $(check_installed $miniconda3_dir) ]; then
    cd $build_dir
    clean "$build_dir/miniconda3"
    download $MINICONDA3_DOWNLOAD_URL "Miniconda3-${MINICONDA3_VERSION}-Linux-x86_64.sh"
    bash Miniconda3-${MINICONDA3_VERSION}-Linux-x86_64.sh -b -p $build_dir/miniconda3
    if [[ "$mainland_china_installation" == "yes" ]]
    then

	$miniconda3_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/main
        $miniconda3_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/free
        $miniconda3_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/pro
        $miniconda3_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/msys2
        $miniconda3_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/bioconda
        $miniconda3_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge

        $miniconda3_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/main
        $miniconda3_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/free
        $miniconda3_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/pro
        $miniconda3_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/msys2
        $miniconda3_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/cloud/bioconda
        $miniconda3_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/cloud/conda-forge

    else 
	$miniconda3_dir/conda config --add channels defaults
	$miniconda3_dir/conda config --add channels bioconda
	$miniconda3_dir/conda config --add channels conda-forge
    fi
    $miniconda3_dir/conda config --set show_channel_urls yes
    $miniconda3_dir/conda config --set ssl_verify no
    $miniconda3_dir/conda config --set channel_priority flexible

    cd $build_dir
    rm Miniconda3-${MINICONDA3_VERSION}-Linux-x86_64.sh 
    note_installed $miniconda3_dir
fi

# --------------- ncbi-blast+ ------------------
echo ""
echo "[$(timestamp)] Installing ncbi-blast+ ..."
blast_dir="$build_dir/ncbi-blast-${BLAST_VERSION}+/bin"
windowmasker_dir="$build_dir/ncbi-blast-${BLAST_VERSION}+/bin"
if [ -z $(check_installed $blast_dir) ]; then
    cd $build_dir
    clean "$build_dir/ncbi-blast-${BLAST_VERSION}+"
    echo ""
    echo "Installing ncbi-blast+ ..."
    echo "Download ncbi-blast-v${BLAST_VERSION}"
    download $BLAST_DOWNLOAD_URL "ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz"
    tar xvzf ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz
    cd $build_dir
    rm ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz
    note_installed $blast_dir
fi    

# # --------------- ncbi-rmblast ------------------
# echo ""
# echo "[$(timestamp)] Installing ncbi-rmblast ..."
# rmblast_dir="$build_dir/ncbi-rmblastn-${RMBLAST_VERSION}/bin"
# if [ -z $(check_installed $rmblast_dir) ]; then
#     cd $build_dir
#     clean "$build_dir/ncbi-rmblastn-${RMBLAST_VERSION}"
#     echo ""
#     echo "Installing ncbi-rmblastn ..."
#     echo "Download ncbi-rmblastn-v${BLAST_VERSION}"
#     download $RMBLAST_DOWNLOAD_URL "ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz"
#     tar xvzf ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz
#     # copy rmblastn binary file to ncbi-blast+ directory for easy RepeatMasker configuration
#     cp $rmblast_dir/rmblastn $blast_dir
#     cd $build_dir
#     rm ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz
#     note_installed $rmblast_dir
# fi

# # ------------- SRA Toolkit -------------------
# echo ""
# echo "[$(timestamp)] Installing SRAtoolkit ..."
# sra_dir="$build_dir/sratoolkit.${SRA_VERSION}-centos_linux64/bin"
# if [ -z $(check_installed $sra_dir) ]; then
#     cd $build_dir
#     echo "Download SRAtoolkit-v${SRA_VERSION}"
#     download $SRA_DOWNLOAD_URL sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz
#     tar xvzf sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz
#     rm sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz
#     note_installed $sra_dir
# fi

# --------------- Trimmomatic -----------------
echo ""
echo "[$(timestamp)] Installing Trimmomatic ..."
trimmomatic_dir="$build_dir/Trimmomatic-${TRIMMOMATIC_VERSION}"
if [ -z $(check_installed $trimmomatic_dir) ]; then
    cd $build_dir
    clean "$build_dir/Trimmomatic-${TRIMMOMATIC_VERSION}"
    echo "Download Trimmomatic-v${TRIMMOMATIC_VERSION}"
    download $TRIMMOMATIC_DOWNLOAD_URL "Trimmomatic-${TRIMMOMATIC_VERSION}.zip"
    unzip Trimmomatic-${TRIMMOMATIC_VERSION}.zip
    cd $trimmomatic_dir
    chmod 755 trimmomatic-${TRIMMOMATIC_VERSION}.jar
    ln -s trimmomatic-${TRIMMOMATIC_VERSION}.jar trimmomatic.jar 
    echo "java -XX:ParallelGCThreads=4 -jar $trimmomatic_dir/trimmomatic-${TRIMMOMATIC_VERSION}.jar \$*"  > trimmomatic
    chmod 755 trimmomatic
    cd $build_dir
    rm Trimmomatic-${TRIMMOMATIC_VERSION}.zip
    note_installed $trimmomatic_dir
fi

#--------------- Guppy-GPU --------------------
echo ""
echo "[$(timestamp)] Installing guppy-GPU ..."
guppy_gpu_dir="$build_dir/ont-guppy-gpu/bin"
if [ -z $(check_installed $guppy_gpu_dir) ]; then
    cd $build_dir
    clean "$build_dir/ont-guppy-gpu"
    echo "Download Guppy-v${GUPPY_VERSION}"
    download $GUPPY_GPU_DOWNLOAD_URL "ont-guppy_${GUPPY_VERSION}_linux64.tar.gz"
    tar xvzf ont-guppy_${GUPPY_VERSION}_linux64.tar.gz
    mv ont-guppy ont-guppy-gpu
    rm ont-guppy_${GUPPY_VERSION}_linux64.tar.gz
    note_installed $guppy_gpu_dir
fi

# --------------- Guppy-CPU --------------------
echo ""
echo "[$(timestamp)] Installing guppy-CPU ..."
guppy_cpu_dir="$build_dir/ont-guppy-cpu/bin"
if [ -z $(check_installed $guppy_cpu_dir) ]; then
    cd $build_dir
    clean "$build_dir/ont-guppy-cpu"
    echo "Download Guppy-v${GUPPY_VERSION}"
    download $GUPPY_CPU_DOWNLOAD_URL "ont-guppy_${GUPPY_VERSION}_linux64.tar.gz"
    tar xvzf ont-guppy_${GUPPY_VERSION}_linux64.tar.gz
    rm ont-guppy_${GUPPY_VERSION}_linux64.tar.gz
    note_installed $guppy_cpu_dir
fi

# -------------------- Dorado --------------------
echo ""
echo "[$(timestamp)] Installing Dorado ..."
dorado_dir="$build_dir/dorado-${DORADO_VERSION}-linux-x64/bin"
if [ -z $(check_installed $dorado_dir) ]; then
    cd $build_dir
    clean "$build_dir/dorado-${DORADO_VERSION}-linux-x64/"
    echo "Download Dorado-v${GUPPY_VERSION}"
    download $DORADO_DOWNLOAD_URL "dorado-${DORADO_VERSION}-linux-x64.tar.gz"
    tar xvzf dorado-${DORADO_VERSION}-linux-x64.tar.gz
    rm dorado-${DORADO_VERSION}-linux-x64.tar.gz
    note_installed $dorado_dir
fi

# # ------------- seqtk -------------------
# echo ""
# echo "[$(timestamp)] Installing seqtk ..."
# seqtk_dir="$build_dir/seqtk"
# if [ -z $(check_installed $seqtk_dir) ]; then
#     cd $build_dir
#     clean "$build_dir/seqtk"
#     echo "Download SEQTK-v${SEQTK_VERSION}"
#     git clone $SEQTK_DOWNLOAD_URL
#     cd seqtk
#     git checkout -f -q $SEQTK_GITHUB_COMMIT_VERSION
#     make
#     cd $build_dir
#     note_installed $seqtk_dir
# fi

# # --------------- Nanofilt --------------------
# echo ""
# echo "[$(timestamp)] Installing NanoFilt ..."
# nanofilt_dir="$build_dir/nanofilt_conda_env/bin"
# if [ -z $(check_installed $nanofilt_dir) ]; then
#     cd $build_dir
#     $miniconda3_dir/conda create -y -p $build_dir/nanofilt_conda_env
#     source $miniconda3_dir/activate $build_dir/nanofilt_conda_env
#     $miniconda3_dir/conda install -y -c bioconda nanofilt=${NANOFILT_VERSION}
#     source $miniconda3_dir/deactivate
#     note_installed $nanofilt_dir
# fi

# --------------- Nanoplot --------------------
echo ""
echo "[$(timestamp)] Installing NanoPlot ..."
nanoplot_dir="$build_dir/nanoplot_conda_env/bin"
if [ -z $(check_installed $nanoplot_dir) ]; then
    cd $build_dir
    $miniconda3_dir/conda create -y -p $build_dir/nanoplot_conda_env python=3.9
    source $miniconda3_dir/activate $build_dir/nanoplot_conda_env
    $miniconda3_dir/conda install -y -c bioconda nanoplot=${NANOPLOT_VERSION}
    source $miniconda3_dir/deactivate
    note_installed $nanoplot_dir
fi

# --------------- Picard -----------------
echo ""
echo "[$(timestamp)] Installing picard ..."
picard_dir="$build_dir/Picard-v${PICARD_VERSION}"
if [ -z $(check_installed $picard_dir) ]; then
    cd $build_dir
    clean "$build_dir/Picard-v${PICARD_VERSION}"
    echo "Download Picard-v${PICARD_VERSION}"
    download $PICARD_DOWNLOAD_URL "picard.jar"
    mkdir Picard-v${PICARD_VERSION}
    mv picard.jar $picard_dir
    cd $picard_dir
    chmod 755 picard.jar
    note_installed $picard_dir
fi

# --------------- gffread ------------------
echo ""
echo "[$(timestamp)] Installing gffread ..."
gffread_dir="$build_dir/gffread"
if [ -z $(check_installed $gffread_dir) ]; then
    cd $build_dir
    clean "$build_dir/gffread"
    echo "Download gffread-v${GFFREAD_VERSION}"
    download $GFFREAD_DOWNLOAD_URL "gffread-${GFFREAD_VERSION}.Linux_x86_64.tar.gz"
    tar xzf gffread-${GFFREAD_VERSION}.Linux_x86_64.tar.gz
    mv gffread-${GFFREAD_VERSION}.Linux_x86_64 gffread
    cd $build_dir
    rm gffread-${GFFREAD_VERSION}.Linux_x86_64.tar.gz
    note_installed $gffread_dir
fi

# -------------- Nanopolish ----------------
echo ""
echo "[$(timestamp)] Installing Nanopolish ..."
nanopolish_dir="$build_dir/nanopolish_conda_env/bin"
if [ -z $(check_installed $nanopolish_dir) ]; then
    cd $build_dir
    clean "$build_dir/nanopolish_conda_env"
    $miniconda3_dir/conda create -y -p $build_dir/nanopolish_conda_env python=3.9
    source $miniconda3_dir/activate $build_dir/nanopolish_conda_env
    $miniconda3_dir/conda install -y -c bioconda nanopolish=${NANOPOLISH_VERSION} 
    source $miniconda3_dir/deactivate
    note_installed $nanopolish_dir
fi

# --------------- XPORE -----------------
echo ""
echo "[$(timestamp)] Installing XPORE ..."
xpore_dir="$build_dir/xpore_conda_env/bin"
if [ -z $(check_installed $xpore_dir) ]; then
    cd $build_dir
    clean "$build_dir/xpore_conda_env"
    echo "Download XPORE-v${XPORE_VERSION}"
    $miniconda3_dir/conda create -y -p $build_dir/xpore_conda_env python=3.9
    source $miniconda3_dir/activate $build_dir/xpore_conda_env
    $miniconda3_dir/conda install -y -c bioconda xpore=${XPORE_VERSION}
    source $miniconda3_dir/deactivate
    note_installed $xpore_dir
fi

# --------------- flair -----------------
echo ""
echo "[$(timestamp)] Installing flair ..."
flair_dir="$build_dir/flair_conda_env/bin"
minimap2_dir="$build_dir/flair_conda_env/bin"
samtools_dir="$build_dir/flair_conda_env/bin"
bedtools_dir="$build_dir/flair_conda_env/bin"
if [ -z $(check_installed $flair_dir) ]; then
    cd $build_dir
    clean "$build_dir/flair_conda_env"
    echo "Download flair-v${FLAIR_VERSION}"
    $miniconda3_dir/conda create -y -p $build_dir/flair_conda_env python=3.9
    source $miniconda3_dir/activate $build_dir/flair_conda_env
    $miniconda3_dir/conda install -y -c conda-forge -c bioconda flair=${FLAIR_VERSION}
    source $miniconda3_dir/deactivate
    # fix a bug of flair
    cd $flair_dir
    head -1 bam2Bed12 > bam2Bed12.header
    tail -n +2 $NANOTRANS_HOME/misc/bam2Bed12 >bam2Bed12.content
    cat bam2Bed12.header bam2Bed12.content > bam2Bed12
    chmod a+x bam2Bed12
    rm bam2Bed12.header
    rm bam2Bed12.content
    # turn on pdf output for plot_isoform_usage
    cp $NANOTRANS_HOME/misc/plot_isoform_usage_for_pdf plot_isoform_usage
    # fix a bug in subset_unassigned_reads.py
    cd $build_dir/flair_conda_env/lib/python3.9/site-packages/flair
    sed -i -e "25s/$/[:line[3].rfind(';')]/" -e "27s/$/[:line[9].rfind(';')]/" subset_unassigned_reads.py
    note_installed $flair_dir
fi

# --------------- rna modification deps -----------------

site_packages="$build_dir/site_packages"
mkdir -p $site_packages
# pkgs used in plotting characteristics of rna modifications
CYUSHUFFLE_VERSION="1.1.2"
PYSAM_VERSION="0.22.0"
TABIX_VERSION="1.11"
MATPLOTLIB_VERSION="3.8.2"
NUMPY_VERSION="1.26.0"
PANDAS_VERSION="2.1.4"

install_pkg_byconda pysam      $PYSAM_VERSION
install_pkg_byconda tabix      $TABIX_VERSION
install_pkg_byconda matplotlib $MATPLOTLIB_VERSION
install_pkg_byconda numpy      $NUMPY_VERSION 
install_pkg_byconda pandas     $PANDAS_VERSION
install_pkg_byconda cyushuffle $CYUSHUFFLE_VERSION

# --------------- bpipe ------------------
echo ""
echo "[$(timestamp)] Installing bpipe ..."
bpipe_dir="$build_dir/bpipe-${BPIPE_VERSION}/bin"
if [ -z $(check_installed $bpipe_dir) ]; then
    cd $build_dir
    clean "$build_dir/bpipe-${BPIPE_VERSION}"
    echo "Download bpipe-v${BPIPE_VERSION}"
    download $BPIPE_DOWNLOAD_URL "bpipe-${BPIPE_VERSION}.tar.gz"
    tar xzf bpipe-${BPIPE_VERSION}.tar.gz
    rm bpipe-${BPIPE_VERSION}.tar.gz
    note_installed $bpipe_dir
fi

# --------------- velvet ------------------
echo ""
echo "[$(timestamp)] Installing velvet ..."
velvet_dir="$build_dir/velvet-${VELVET_VERSION}"
if [ -z $(check_installed $velvet_dir) ]; then
    cd $build_dir
    clean "$build_dir/velvet-${VELVET_VERSION}"
    echo "Download velvet-v${VELVET_VERSION}"
    download $VELVET_DOWNLOAD_URL "velvet-${VELVET_VERSION}.tar.gz"
    tar xzf velvet-${VELVET_VERSION}.tar.gz
    cd "$build_dir/velvet-${VELVET_VERSION}"
    make MAXKMERLENGTH=37 LONGSEQUENCES=1 #OPENMP=1
    cd $build_dir
    rm velvet-${VELVET_VERSION}.tar.gz
    note_installed $velvet_dir
fi

# --------------- oases ------------------
echo ""
echo "[$(timestamp)] Installing oases ..."
oases_dir="$build_dir/oases-${OASES_VERSION}"
if [ -z $(check_installed $oases_dir) ]; then
    cd $build_dir
    clean "$build_dir/oases-${OASES_VERSION}"
    echo "Download oases-v${OASES_VERSION}"
    download $OASES_DOWNLOAD_URL "oases-${OASES_VERSION}.tar.gz"
    tar xzf oases-${OASES_VERSION}.tar.gz
    cd "$build_dir/oases-${OASES_VERSION}"
    make MAXKMERLENGTH=37 LONGSEQUENCES=1 "VELVET_DIR=$build_dir/velvet-${VELVET_VERSION}" #OPENMP=1
    cd $build_dir
    rm oases-${OASES_VERSION}.tar.gz
    note_installed $oases_dir
fi

# --------------- bowtie2 ------------------
echo ""
echo "[$(timestamp)] Installing bowtie2 ..."
bowtie2_dir="$build_dir/bowtie2-${BOWTIE2_VERSION}-linux-x86_64"
if [ -z $(check_installed $bowtie2_dir) ]; then
    cd $build_dir
    clean "$build_dir/bowtie2-${BOWTIE2_VERSION}"
    echo "Download bowtie2-v${BOWTIE2_VERSION}"
    download $BOWTIE2_DOWNLOAD_URL "bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip"
    unzip bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip
    cd $build_dir
    rm bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip
    note_installed $bowtie2_dir
fi

# --------------- fastx_toolkit ------------------
echo ""
echo "[$(timestamp)] Installing fastx_toolkit ..."
fastxtk_dir="$build_dir/fastx_toolkit/bin"
if [ -z $(check_installed $fastxtk_dir) ]; then
    cd $build_dir
    clean "$build_dir/fastx_toolkit"
    mkdir fastx_toolkit
    cd fastx_toolkit
    echo "Download fastx_toolkit-v${FASTXTK_VERSION}"
    download $FASTXTK_DOWNLOAD_URL "fastx_toolkit_${FASTXTK_VERSION}.tar.bz2"
    tar -xjf fastx_toolkit_${FASTXTK_VERSION}.tar.bz2
    rm fastx_toolkit_${FASTXTK_VERSION}.tar.bz2
    note_installed $fastxtk_dir
fi

# --------------- BBmap ------------------
echo ""
echo "[$(timestamp)] Installing BBmap ..."
bbmap_dir="$build_dir/bbmap/bin"
if [ -z $(check_installed $bbmap_dir) ]; then
    cd $build_dir
    clean "$build_dir/bbmap-${BBMAP_VERSION}"
    echo "Download bbmap-v${BBMAP_VERSION}"
    download $BBMAP_DOWNLOAD_URL "bbmap-${BBMAP_VERSION}.tar.gz"
    tar -xzf bbmap-${BBMAP_VERSION}.tar.gz
    cd bbmap
    mkdir bin
    cd bin
    for script in $(ls ./../*.sh) ; do
	s=$(basename $script)
	s_pre=$(echo $s | sed 's/.sh//g')
	ln -s $script $s_pre
	chmod 755 $s_pre
    done
    cd $build_dir
    rm bbmap-${BBMAP_VERSION}.tar.gz
    note_installed $bbmap_dir
fi

# --------------- UCSC -----------------
ucsc_dir="$build_dir/UCSC_Utilities"
echo ""
echo "[$(timestamp)] Installing UCSC utilities ..."
if [ -z $(check_installed $ucsc_dir) ]; then
    cd $build_dir
    clean "$build_dir/UCSC_Utilities"
    mkdir UCSC_Utilities
    cd $ucsc_dir
    download $BEDPARTITION_DOWNLOAD_URL "bedPartition"
    download $GTFTOGENEPRED_DOWNLOAD_URL "gtfToGenePred"
    download $BLAT_DOWNLOAD_URL "blat"
    chmod 755 $ucsc_dir/*
    note_installed $ucsc_dir
fi

# --------------- JAFFAL -----------------
jaffal_dir="$build_dir/JAFFA-version-${JAFFAL_VERSION}"
echo ""
echo "[$(timestamp)] Installing JAFFAL ..."
if [ -z $(check_installed $jaffal_dir) ]; then
    cd $build_dir
    download $JAFFAL_DOWNLOAD_URL "JAFFA-version-${JAFFAL_VERSION}.tar.gz"
    tar xzf JAFFA-version-${JAFFAL_VERSION}.tar.gz
    rm JAFFA-version-${JAFFAL_VERSION}.tar.gz
    cd JAFFA-version-${JAFFAL_VERSION}
    mkdir bin
    g++ -std=c++11 -O3 -o ./bin/make_3_gene_fusion_table ./src/make_3_gene_fusion_table.c++
    g++ -std=c++11 -O3 -o ./bin/extract_seq_from_fasta ./src/extract_seq_from_fasta.c++
    g++ -std=c++11 -O3 -o ./bin/make_simple_read_table ./src/make_simple_read_table.c++
    g++ -std=c++11 -O3 -o ./bin/process_transcriptome_align_table ./src/process_transcriptome_align_table.c++
    g++ -O3 -o ./bin/make_count_table ./src/make_count_table.c++
    echo "minimap2=\"$minimap2_dir/minimap2\"" >> tools.groovy
    echo "blastn=\"$blast_dir/blastn\"" >> tools.groovy
    echo "makeblastdb=\"$blast_dir/makeblastdb\"" >> tools.groovy
    echo "bpipe=\"$bpipe_dir/bpipe\"" >> tools.groovy
    echo "velvetg=\"$velvet_dir/velvetg\"" >> tools.groovy
    echo "velveth=\"$velvet_dir/velveth\"" >> tools.groovy
    echo "oases=\"$velvet_dir/oases\"" >> tools.groovy
    echo "trimmomatic=\"$trimmomatic_dir/trimmomatic\"" >> tools.groovy
    echo "samtools=\"$samtools_dir/samtools\"" >> tools.groovy
    echo "bowtie2=\"$bowtie2_dir/bowtie2\"" >> tools.groovy
    echo "blat=\"$ucsc_dir/blat\"" >> tools.groovy
    echo "fasta_formatter=\"$fastxtk_dir/fasta_formatter\"" >> tools.groovy
    echo "dedupe=\"$bbmap_dir/dedupe\"" >> tools.groovy
    echo "reformat=\"$bbmap_dir/reformat\"" >> tools.groovy
    echo "make_3_gene_fusion_table=\"$jaffal_dir/bin/make_3_gene_fusion_table\"" >> tools.groovy
    echo "extract_seq_from_fasta=\"$jaffal_dir/bin/extract_seq_from_fasta\"" >> tools.groovy
    echo "make_simple_read_table=\"$jaffal_dir/bin/make_simple_read_table\"" >> tools.groovy
    echo "make_count_table=\"$jaaffal_dir/bin/make_count_table\"" >> tools.groovy
    echo "process_transcriptome_align_table=\"$jaffal_dir/bin/process_transcriptome_align_table\"" >> tools.groovy
    echo "R=\"$r_path\"" >> tools.groovy
    cp $jaffal_dir/known_fusions.txt $NANOTRANS_HOME/misc/human.known_fusions.txt
    note_installed $jaffal_dir
fi

# --------------- Quarto -----------------

install_and_create_pkg_condaenv quarto ${QUARTO_VERSION}


echo ""
echo "#########################################################"
# Configure executable paths
cd $NANOTRANS_HOME
echo ""
echo "[$(timestamp)] All dependency instalations have finished." 
echo "[$(timestamp)] Configuring executable paths ..."
echo "export NANOTRANS_HOME=${NANOTRANS_HOME}" > env.sh
echo "export build_dir=${build_dir}" >> env.sh
echo "export PYTHONPATH=${PYTHONPATH}" >> env.sh
echo "export R_LIBS=${rlib_dir}:${build_dir}/flair_conda_env/lib/R/library" >> env.sh
echo "export PERL5LIB=${PERL5LIB}" >> env.sh 
echo "export cpanm_dir=${cpanm_dir}" >> env.sh
echo "export miniconda3_dir=${miniconda3_dir}" >> env.sh
#echo "export sra_dir=${sra_dir}" >> env.sh
echo "export trimmomatic_dir=${trimmomatic_dir}" >> env.sh
echo "export guppy_gpu_dir=${guppy_gpu_dir}" >> env.sh
echo "export guppy_cpu_dir=${guppy_cpu_dir}" >> env.sh
echo "export dorado_dir=${dorado_dir}" >> env.sh
echo "export nanopolish_dir=${nanopolish_dir}" >> env.sh
echo "export blast_dir=${blast_dir}" >> env.sh
# echo "export seqtk_dir=${seqtk_dir}" >> env.sh
# echo "export nanofilt_dir=${nanofilt_dir}" >> env.sh
echo "export nanoplot_dir=${nanoplot_dir}" >> env.sh
echo "export minimap2_dir=${minimap2_dir}" >> env.sh
echo "export gffread_dir=${gffread_dir}" >> env.sh
echo "export samtools_dir=${samtools_dir}" >> env.sh
echo "export bedtools_dir=${bedtools_dir}" >> env.sh
#echo "export htslib_dir=${htslib_dir}" >> env.sh
#echo "export tabix_dir=${tabix_dir}" >> env.sh
echo "export picard_dir=${picard_dir}" >> env.sh
echo "export xpore_dir=${xpore_dir}" >> env.sh
echo "export flair_dir=${flair_dir}" >> env.sh
echo "export bowtie2_dir=${bowtie2_dir}" >> env.sh
echo "export fastxtk_dir=${fastxtk_dir}" >> env.sh
echo "export bbmap_dir=${bbmap_dir}" >> env.sh
echo "export oases_dir=${oases_dir}" >> env.sh
echo "export velvet_dir=${velvet_dir}" >> env.sh
echo "export bpipe_dir=${bpipe_dir}" >> env.sh
echo "export ucsc_dir=${ucsc_dir}" >> env.sh
echo "export jaffal_dir=${jaffal_dir}" >> env.sh
echo "export quarto_dir=${build_dir}/quarto_conda_env/bin" >> env.sh


# test java configuration: requireds java 1.8 
echo ""
echo "#########################################################"
echo "[$(timestamp)] Testing java configuration ..."
echo ""
java_bin=""
if type -p java
then 
    java_bin=$(which java)
    echo "found java executable in PATH: $java_bin"
elif [[ -n "$JAVA_HOME" ]] && [[ -x "$JAVA_HOME/bin/java" ]]
then 
    java_bin="$JAVA_HOME/bin/java"
    echo "found java executable in JAVA_HOME: $java_bin" 
else 
    echo "";
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
    echo "Failed to detect Java installation in the system!"
    echo "Please install java 1.8, which is a dependency of NanoTrans!\n";
    echo "After the java installation, please manually set the directory path to java 1.8 executable on the last line of the env.sh file generated by this installation script!"
    echo "export java_dir=" >> env.sh
fi  

if [[ -n "$java_bin" ]]
then
    java_version=$("$java_bin" -version 2>&1 | awk -F '"' '/version/ {print $2}')
    echo "detected java_version: $java_version"
    if [ $(tidy_version "$java_version") -eq $(tidy_version "1.8") ]
    then
	java_dir=$(dirname $java_bin)
	echo "export java_dir=${java_dir}" >> env.sh
        echo "You have the correct java version for NanoTrans! NanoTrans will take care of the configuration."
    else
	echo "";
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
	echo "Your java version is not the version required by NanoTrans (java v1.8)!"
        echo "Please manually set the directory path to java 1.8 executable on the last line of the env.sh file generated by this installation script!"
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
	echo "export java_dir=" >> env.sh
    fi
fi

echo ""
echo "#########################################################"
# echo "[$(timestamp)] Decompress large supporting files ..."
# gunzip $NANOTRANS_HOME/data/Proteome_DB_for_annotation.CDhit_I95.fa.gz
echo "done!"
echo ""
echo ""



############################
# checking Bash exit status
if [[ $? -eq 0 ]]
then
    echo "#########################################################"
    echo ""
    echo "[$(timestamp)] NanoTrans message: This bash script has been successfully processed! :)"
    echo ""
    echo "#########################################################"
    echo ""
    exit 0
fi
############################

