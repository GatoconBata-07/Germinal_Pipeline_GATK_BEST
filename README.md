# Germline WES Automated Pipeline

## Overview
This repository contains a robust, end-to-end automated bioinformatics pipeline (`germinal_auto.py`) designed for the identification of germline variants from Whole Exome Sequencing (WES) data. The pipeline strictly adheres to the **GATK Best Practices (https://gatk.broadinstitute.org/hc/en-us/articles/360035894711-About-the-GATK-Best-Practices)** for variant discovery, ensuring high-quality, reproducible results suitable for clinical genomic research.

## Pipeline Steps
1. **Alignment:** Fastq files are aligned to the reference genome using `BWA-MEM`.
2. **Sorting & Indexing:** Alignments are sorted and indexed using `Samtools`.
3. **Duplicate Marking:** PCR duplicates are flagged using `GATK MarkDuplicatesSpark`.
4. **Base Quality Score Recalibration (BQSR):** Systematic errors are corrected using `GATK BaseRecalibrator` and `ApplyBQSR`.
5. **Variant Calling:** Germline variants are called using `GATK HaplotypeCaller`.
6. **Hard Filtering:** SNPs and INDELs are separated and filtered using strict threshold parameters via `GATK VariantFiltration` to minimize false positives.
7. **Merging:** Filtered SNPs and INDELs are merged back into a final, high-quality VCF file.

## Prerequisites
Ensure the following tools are installed and available in your system's `$PATH` (or use the Docker image):
* Python 3.7+
* BWA (version: 0.7.17-r1188)
* GATK4 (version 4.6)

## Required Files
1. Paired-end FASTQ files (`_R1.fastq.gz`, `_R2.fastq.gz`)
2. Reference Genome in FASTA format (e.g., `hg38.fasta`). **Note:** The reference genome must be indexed (`.fai`, `.dict`, and BWA index files).
3. Known Sites VCF file for BQSR (e.g., `dbsnp_146.hg38.vcf.gz`), also indexed (`.tbi`).

## Usage
Run the pipeline via the command line:

./germinal_auto.py
