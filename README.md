# Germline WES Automated Pipeline (Python + Docker Orchestration)

## Overview
This repository contains a robust, end-to-end automated bioinformatics pipeline (`germinal_auto.py`) designed for the identification of germline variants from Whole Exome Sequencing (WES) data. The pipeline strictly adheres to the **[GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035894711-About-the-GATK-Best-Practices)** for variant discovery, ensuring high-quality, reproducible results suitable for clinical genomic research.

To maximize both performance and absolute clinical reproducibility, this pipeline utilizes a **Hybrid Orchestration Architecture**:
1. **Native Execution:** High-throughput alignment and sorting steps are executed locally on the host Linux system using `BWA-MEM` and `Samtools` for optimal I/O performance.
2. **Docker Orchestration:** Complex variant calling and annotation steps are containerized. The Python script dynamically pulls and runs the **official, unaltered Docker images** from the Broad Institute (GATK4) and Ensembl (VEP). This prevents dependency conflicts and ensures 100% adherence to GATK Best Practices.

## Pipeline Steps
1. **Alignment:** Fastq files are aligned to the reference genome using `BWA-MEM`.
2. **Sorting & Indexing:** Alignments are sorted and indexed using `Samtools`.
3. **Duplicate Marking:** PCR duplicates are flagged using `GATK MarkDuplicates`.
4. **Base Quality Score Recalibration (BQSR):** Systematic errors are corrected using `GATK BaseRecalibrator` and `ApplyBQSR`.
5. **Variant Calling:** Germline variants are called using `GATK HaplotypeCaller`.
6. **Hard Filtering:** SNPs and INDELs are separated and filtered using strict threshold parameters via `GATK VariantFiltration` to minimize false positives.
7. **Merging:** Filtered SNPs and INDELs are merged back into a final, high-quality VCF file.

## System Requirements & Dependencies
To run this pipeline, the host system must be a Linux-based environment with the following installed:

### 1. Local Tools (Host Environment)
* **Python** (v3.7 or higher, version 3.10.12 recommended)
* **BWA** (version 0.7.17 recommended)
* **Samtools** (version 1.13 recommended)

### 2. Docker Engine
* **Docker** must be installed and the user running the Python script must have permissions to execute Docker commands (e.g., added to the `docker` group).
* **Internet Connection:** Required during the first run to pull the official images.

### 3. Official Docker Images utilized (Auto-pulled by the script)
* `broadinstitute/gatk:4.6.0.0` (For MarkDuplicates, BQSR, HaplotypeCaller, VariantFiltration)
* `ensemblorg/ensembl-vep:release_112.0` (For variant functional annotation) 

## Required Files
1. Paired-end FASTQ files (`_R1.fastq.gz`, `_R2.fastq.gz`)
2. Reference Genome in FASTA format (e.g., `hg38.fasta`). **Note:** The reference genome must be indexed (`.fai`, `.dict`, and BWA index files).
3. Known Sites VCF file for BQSR (`Homo_sapiens_assembly38.dbsnp138.vcf.gz`), also indexed (`.tbi`).

## Hard Filtering Quality Parameters
To ensure maximum clinical accuracy, the script applies the following strict filters using `GATK VariantFiltration`:
* **SNPs:** `QD < 2.0 || QUAL < 30.0 || SOR > 3.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0`
* **INDELs:** `QD < 2.0 || QUAL < 30.0 || FS > 200.0 || ReadPosRankSum < -20.0`

## Usage

Unlike tools that rely on complex command-line arguments, this pipeline is designed with a **hardcoded configuration block** for strict execution tracking (highly recommended in clinical research settings). 

### Step 1: Configure the Script
Open the `germinal_auto.py` file in your preferred text editor and locate the configuration variables at the top of the script. Modify these absolute paths and variables according to your specific patient data and system architecture:

```python
# --- CONFIGURATION BLOCK ---
sample_name = "Patient_01"  # Define your output prefix
fastq1 = "/path/to/sample_R1.fastq.gz"
fastq2 = "/path/to/sample_R2.fastq.gz"
reference_genome = "/path/to/hg38.fasta"
known_sites = "/path/to/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
output_dir = "/path/to/output_folder"
```

### Step 2: Execute the Pipeline
Once the variables are set and saved, execute the script directly from your terminal:

```bash
chmod +x germinal_auto.py
./germinal_auto.py
```
