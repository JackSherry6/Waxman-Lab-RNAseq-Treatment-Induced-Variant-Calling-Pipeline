# Waxman Lab RNA-seq SNP Calling Pipeline for CD1 Mouse Reference Genome

This project was prompted by my desire to improve the quality of SNP calling from RNA-seq data by using a CD1-specific reference genome (Jung et al., 2023). I was tasked with developing a method to accurately detect mutations driving tumorigenesis across the genome from RNA-seq data generated from our outbred CD1 mouse strain. Previously, our lab aligned reads and performed differential expression analyses using the mm9 reference genome; however, this approach was inadequate for identifying polymorphisms in an outbred population. The mm9 reference represents a single, inbred laboratory strain and therefore lacks the genetic diversity present in CD1 mice. As a result, alignment accuracy was reduced, variant calls were unreliable, and many true CD1-specific SNPs were missed. Additionally, there was no reliable gtf file or known snps vcf for our CD1 mouse strain. Since our lab's main focus is mutagen-induced liver cancer, I wanted the pipeline to handle tumor samples accurately, while also being able to process non-tumor samples in the same run and compare both. Therefore I decided to switch the tool to varscan somatic for control vs tumor variant calling (although I am currently looking to switch to Mutect2). These improvements allowed for more accurate mapping, improved SNP detection, and better representation of the genetic variation inherent to this outbred strain. 

**Add recent progress

## Table of Contents
1. Features
2. Requirements
3. Installation
4. Usage
5. Configuration
6. Input and Output
7. Contributing
8. References

## Features
- Modular Nextflow pipeline with clearly separated steps:
  - Read preprocessing and quality control
  - Alignment to CD1 reference genome
  - Custom experimental group combining and merging 
  - SNP calling and variant filtering
  - Supports BU HPC cluster execution.
- Docker/Singularity container support for reproducibility.
- Automatic logging and error handling.
- Scalable to large RNA-seq datasets.

## Requirements
- Must have a conda environment with nextflow in order to run nextflow
- Modules already installed on BU Shared Computing cluster (SCC)
- If using aws, see envs file for all packages to install
- If not using BU SCC, see envs directory for software and version information
 
## Installation
  - Clone this repository in the SCC
  - git clone 'https://github.com/JackSherry6/cd1-varscan-snp-pipeline.git'
 
## Usage
Basic execution: 
- ```module load miniconda```
- ```conda activate <name_of_your_nexflow_conda_env>```
- See configuration in order to set sample paths, variables and names
- ```nextflow run main.nf -profile conda,cluster``` (for waxman lab you should always run on the cluster, but if using aws, substitute ```aws``` for ```cluster```)
- NOTE: other nextflow commands will work but these are the ones I use and recommend.

## Configuration
- Lines for configuration in config file:
  - Fill in the example samplesheet with your fastq paths (name sure the sample names end in an alphabetical character and then a number ex. _d5).
  - Mandatory: Set ref_genome and gtf to their respective paths.
  - Mandatory: Set control_cnt to the number of control samples multiplied by 10.
  - Mandatory: Define your group names as the same names you used in the samplesheet minus the trailing integer. Add or remove groups depending on how many groups                you are analyzing.
  - Optional: set paths for fa_dict, star_index, lncRNAs_ref, and known_snps_vcf or set them to 'null' (fa_dict and star_index will help speed the pipeline up                    while lncRNAs_ref and known_snps_vcfgrant will additional output data).

## Inputs and Outputs
- Input:
  - Sample Fastq files
  - Reference Fasta file
  - Matching reference GTF file (preferably with full transcriptomic features, not exonic only)
  - Optional: fasta dictionary, star index, known lncRNA vcf, and known snps vcf
- Output:
  - Comprehensive Multiqc report
  - Bam and bam index files for genome browser viewing
  - Per-sample vcf files (snp and indel)
  - Annotated, processed, group specific csv files (snp and indel)

## Contributing 
- Email me at jgsherry@bu.edu for additional information or contributing information

## References

Jung, Y. H., Wang, H.-L. V., Ali, S., Corces, V. G., … & Kremsky, I. (2023). Characterization of a strain-specific CD-1 reference genome reveals potential inter- and intra-strain functional variability. BMC Genomics, 24, Article 437. https://doi.org/10.1186/s12864-023-09523-x
