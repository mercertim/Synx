# Synx
Analysis scripts for Synx pipeline

## Synx analysis pipeline
These scripts a required for the analysis of ONT or Illumina data using Synx. These can be used with DNA or RNA (cDNA) libraries using the aligned bam files as input. A demo dataset (ONT fasta file) and reference files for use with these scripts are contained in this repository. Input and output directories need to be specified by the user along with reference files according to the genome of interest.

## System requirements and installation guide
The Synx analysis is performed by third party software and does not require installation beyond the below dependencies. Running the demo should take <1min. While the below demo can be run on a standard laptop, we recommend using a high performance computing (HPC) cluster for data preprocessing and calculating per base error stats. The below demonstration has been tested on a HPC 64-bit running CentOS Linux release 7.9.2009 with python (3.6.7) and R (4.0.2). Downstream analysis in R can be run on a standard laptop computer and was tested on an Apple Macbook 16G RAM, 500GB memory.

Dependencies command-line:

minimap2 (v2.24) https://github.com/lh3/minimap2 (for alignment of ONT data, not required if starting with an aligned bam file)

samtools (using htslib v1.9) http://www.htslib.org/

Dependencies Python:

pip install pysamstats (v1.1.2)

## Demo
A sample Oxford nanopore dataset containing a Synx library is contained in the file ont_demo_data.fastq for demonstration. To run:

#### Clone repository scripts and reference files
git clone https://github.com/mercertim/Synx.git

cd Synx

#### Get demo datasets
curl -OL https://github.com/mercertim/Synx/releases/download/v0.1/ont_demo_data.fastq

curl -OL https://github.com/mercertim/Synx/releases/download/v0.1/illumina_demo_data.bam

#### Align to Synx sequence with minimap2 and sort and index bam
minimap2 -ax map-ont -t 8 synx.fa ont_demo_data.fastq | samtools sort - > ont_demo_data_synx.bam

samtools index ont_demo_data_synx.bam

#### Get pileup stats per base for Synx genome for Ont and Illumina (starting with aligned bam)
pysamstats --fasta synx.fa --type variation ont_demo_data_synx.bam > ont_demo_data_synx.bam.bed

pysamstats --fasta synx.fa --type variation illumina_demo_data.bam > illumina_demo_data.bam.bed
 
#### Collate pileup stats 
python3 analyzePile.py ont_demo_data_synx.bam.bed > ont_demo_data_synx.bam.bed.tsv

python3 analyzePile.py illumina_demo_data.bam.bed > illumina_demo_data.bam.bed.tsv

The expected output for Illumina and ONT libraries can be found in the files test_ont_demo_data_synx.bam.bed.tsv and test_illumina_demo_data_synx.bam.bed.tsv
