# NanoTrans

<p align="center">
  <img src="https://github.com/yjx1217/NanoTrans/blob/main/NanoTrans.logo.png" alt="NanoTrans_logo" width="408" height="132"/>
</p>

**NanoTrans: An integrated computational framework for comprehensive transcriptome analyses with Nanopore direct-RNA sequencing**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Description
<div style="text-align: justify"> 
Nanopore direct-RNA sequencing (DRS) provides the direct access to native RNA strands with full-length information, shedding light on rich qualitative and quantitative properties of gene expression profiles. Here with NanoTrans, we present an integrated computational framework that comprehensively covers all major DRS-based application scopes, including isoform clustering and quantification, poly(A) tail length estimation, RNA modification profiling, and gene fusion detection. In addition to its merit in providing such a streamlined one-stop solution, NanoTrans also shines in its workflow-orientated modular design, batch processing capability, rich tabular and graphic report outputs, as well as automatic installation and configuration support. Given the rising adoption of Nanopore DRS technology in the field, we believe NanoTrans will become a highly useful tool to help researchers to fully explore the power of this exciting technology with rich biological insights obtained.
</div>

<p align="center">
  <img src="https://github.com/yjx1217/NanoTrans/blob/main/NanoTrans.overview.png" alt="NanoTrans_overview" width="800" height="851"/>
</p>

Under the hood, a series of task-specific modules are provided to carry out the full workflow of NanoTrans:

* **00.Reference_Genome**
  * donwloading and preprocessing the reference genome and annotation
* **00.Long_Reads**
  * performing basecalling and length/quality summarization of raw Nanopore DRS fast5 reads
* **01.Reference_Genome_based_Read_Mapping**
  * mapping the nanopore DRS reads against the reference genome
* **02.Isoform_Clustering_and_Quantification**
  * clustering and polishing isoforms and quantifying their expression levels
* **03.Isoform_Expression_and_Splicing_Comparison**
  * comparing isoform usages and splicing preferences among different sample groups or samples
* **04.Isoform_RNA_Modification_Identification**
  * identifying RNA modification profile of each isoform
* **05.Isoform_PolyA_Tail_Length_Profiling**
  * profiling poly(A) tail length of each isoform
* **06.Gene_Fusion_Detection**
  * identifying gene fusion based on the chimeric isoform evidence


## Citation
Fan Wang, Xinxin Zhang, Li Zhang, Jing Li, Jia-Xing Yue. (2022) NanoTrans: an integrated computational framework for comprehensive transcriptome analyses with Nanopore direct-RNA sequencing. BioRxiv. (doi: https://doi.org/10.1101/2022.11.29.518309)


## License
NanoTrans itself is distributed under the MIT license but some of its dependencies might have more strict license for commercial use. Please check the licensing details of those dependencies.

## Release history
* v0.0.1 Released on 2022/11/30

## Installation
```sh
git clone https://github.com/yjx1217/NanoTrans.git
cd NanoTrans
bash ./install_dependencies.sh
```
If the installation succeeds, you should see the following massage:
“NanoTrans message: This bash script has been successfully processed! :)”
This signifies the success of the installation process. 

Upon the success of the installation, a subdirectory named build and a file named env.sh will be generated. The build subdirectory holds all the installed dependencies, while the env.sh file contains the execution paths of these dependencies. This file will be automatically loaded to set up the working environment for various modules of NanoTrans. The base directory of NanoTrans is defined as $NANOTRANS_HOME in this file.

If unexpected error occurs during installation, normally you can just re-do “bash ./install_dependencies.sh” step and the installation should be able to automatically resume from the previous interruption point. 


## Requirements
### Hardware, operating system and network requirements
NanoTrans is designed for a desktop or computing server running an x86-64-bit Linux operating system. Multithreaded processors are preferred to speed up the process since many modules support multithreaded processing. A stable internet connection is required for its installation. 

### Software requirements
* bash (https://www.gnu.org/software/bash/)
* bzip2 and libbz2-dev (http://www.bzip.org/)
* curl (https://curl.haxx.se/)
* gcc and g++ (https://gcc.gnu.org/)
* git (https://git-scm.com/)
* GNU make (https://www.gnu.org/software/make/)
* gzip (https://www.gnu.org/software/gzip/)
* libopenssl-devel
* libcurl-devel
* java runtime environment (JRE) v1.8.0 (https://www.java.com)
* perl v5.12 or newer (https://www.perl.org/)
* tar (https://www.gnu.org/software/tar/)
* unzip (http://infozip.sourceforge.net/UnZip.html)
* wget (https://www.gnu.org/software/wget/)
* zlib and zlib-devel (https://zlib.net/)
* xz and xz-devel (https://tukaani.org/xz/)
