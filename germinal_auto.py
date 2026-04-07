#!/usr/bin/env python3

import os
import subprocess
import sys
import logging
import shutil

# First, in your terminal where the script is located, run > chmod +x germinal_auto.py

# TO RUN THIS PIPELINE: SAVE THE SCRIPT AND EXECUTE > ./germinal_auto.py

# =============================
# Global Configuration
# =============================

threads = 16
sample_name = "my_sample_name"

# Input .fastq files 

input_fastq1 = "my_sample_name_R1.fastq.gz"
input_fastq2 = "my_sample_name_R2.fastq.gz"

" required reference files > reference_genome = hg38.fa and known_sites = Homo_sapiens_assembly38.dbsnp138.vcf "

# Reference files (using the absolute path of the host and the same path within the container)

# Replace the path: /home/gato/RefGen/Tumoral_Exomes_Ref/ with the /folder_path_to_my_files 

reference_genome = "/home/gato/RefGen/Tumoral_Exomes_Ref/hg38.fa"

known_sites = "/home/gato/RefGen/Tumoral_Exomes_Ref/Homo_sapiens_assembly38.dbsnp138.vcf"

# Output Directory

output_dir = "output"

# Docker Images

docker_gatk = "broadinstitute/gatk:4.6.0.0"
docker_vep = "ensemblorg/ensembl-vep:release_112.0"

# Route where the VEP caches are located

# Replace the path: /media/gato/22D187F822421089/vep_data with the /folder_path_to_vep_data

vep_data_dir = "/media/gato/22D187F822421089/vep_data"

# Log
log_file = os.path.join(output_dir, "pipeline.log")

# =============================
# Output files
# =============================
output_sam = os.path.join(output_dir, f"{sample_name}.aligned.sam")
output_bam = os.path.join(output_dir, f"{sample_name}.sorted_dedup.bam")
recal_table = os.path.join(output_dir, f"{sample_name}.recal_data.table")
recalibrated_bam = os.path.join(output_dir, f"{sample_name}.recalibrated.bam")
alignment_metrics = os.path.join(output_dir, "alignment_metrics.txt")
insert_metrics = os.path.join(output_dir, "insert_size_metrics.txt")
insert_histogram = os.path.join(output_dir, "insert_size_histogram.pdf")
raw_variants_vcf = os.path.join(output_dir, "raw_variants.vcf")
snps_vcf = os.path.join(output_dir, "raw_snps.vcf")
indels_vcf = os.path.join(output_dir, "raw_indels.vcf")
filtered_snps_vcf = os.path.join(output_dir, "filtered_snps.vcf")
filtered_indels_vcf = os.path.join(output_dir, "filtered_indels.vcf")
analysis_ready_snps_vcf = os.path.join(output_dir, "analysis-ready-snps.vcf")
analysis_ready_indels_vcf = os.path.join(output_dir, "analysis-ready-indels.vcf")
filtered_snps_gt_vcf = os.path.join(output_dir, "analysis-ready-snps-filteredGT.vcf")
filtered_indels_gt_vcf = os.path.join(output_dir, "analysis-ready-indels-filteredGT.vcf")
merged_vcf = os.path.join(output_dir, "final_merge.filtered.vcf")
annotated_vcf = os.path.join(output_dir, f"{sample_name}.final.annotated.vcf")

# =============================
# Logging Configuration
# =============================

os.makedirs(output_dir, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    filename=log_file,
    filemode='w'
)
logging.info("Pipeline starting")

# =============================
# Utility functions
# =============================

def validate_input_files():
    """Verify that the input files exist before continuing."""
    required_files = [input_fastq1, input_fastq2, reference_genome, known_sites]
    for file in required_files:
        if not os.path.exists(file):
            logging.error(f"Input file not found: {file}")
            raise FileNotFoundError(f"Required file not found: {file}")
    logging.info("All input files are present.")

def run_command(cmd, docker_image=None, volumes=None, workdir="/gatk/data"):
    """
    Execute a command in the shell or within Docker.
      - docker_image: str with the image name (ej. 'broadinstitute/gatk:4.6.0.0')
      - volumes: list of tuples (host_dir, container_dir) for mounting additional volumes.
      - workdir: working directory within the container.
    """
    try:
        if docker_image:
            # We built the string of assemblies (-v)
            #
            # 1) Mount the pwd /gatk/data for (BAMs, VCFs, logs, etc.)
            volume_str = "-v $(pwd):/gatk/data"

            # 2) Mount the references folder in the same absolute path
            #    from the host inside the container, so as not to change the paths:
            volume_str += " -v /home/gato/RefGen/Tumoral_Exomes_Ref:/home/gato/RefGen/Tumoral_Exomes_Ref"

            # 3) Mount additional volumes if specified
            if volumes:
                for (host_path, container_path) in volumes:
                    volume_str += f" -v {host_path}:{container_path}"

            docker_cmd = (
                f"docker run --rm "
                f"{volume_str} "
                f"-w {workdir} "
                f"{docker_image} "
                f"{cmd}"
            )
            logging.info(f"Running Docker command: {docker_cmd}")
            subprocess.run(docker_cmd, shell=True, check=True)
        else:
            # Local execution (outside of container)
            logging.info(f"Running local command: {cmd}")
            subprocess.run(cmd, shell=True, check=True)

        logging.info(f"Command completed: {cmd}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error executing the command: {cmd}")
        raise e

# =============================
# STEPS
# =============================
def align_bwa_mem():
    """STEP 1: BWA-MEM alignment (local run)."""
    logging.info("Initiating alignment with BWA-MEM...")
    read_group = f"@RG\\tID:{sample_name}\\tPL:ILLUMINA\\tSM:{sample_name}"
    cmd = (
        f"bwa mem -t {threads} "
        f"-R \"{read_group}\" "
        f"{reference_genome} {input_fastq1} {input_fastq2} > {output_sam}"
    )
    run_command(cmd, docker_image=None)  # BWA local
    logging.info("Alignment completed")

def mark_duplicates():
    """STEP 2: Mark duplicates with GATK (within Docker GATK)."""
    logging.info("Marking duplicates...")
    cmd = (
        f"gatk MarkDuplicatesSpark "
        f"-I {output_sam} "
        f"-O {output_bam}"
    )
    run_command(
        cmd,
        docker_image=docker_gatk,
        volumes=None,  # Use the pwd mount by default -> /gatk/data + /home/gato/RefGen/Tumoral_Exomes_Ref
        workdir="/gatk/data"
    )
    logging.info("Duplicates marked and BAM file generated.")
    os.remove(output_sam)

def base_quality_recalibration():
    """STEP 3: Base quality recalibration with GATK."""
    logging.info("Recalibrating baseline quality...")

    # BaseRecalibrator
    cmd_recal = (
        f"gatk BaseRecalibrator "
        f"-I {output_bam} "
        f"-R {reference_genome} "
        f"--known-sites {known_sites} "
        f"-O {recal_table}"
    )
    run_command(cmd_recal, docker_image=docker_gatk)

    # ApplyBQSR
    cmd_apply = (
        f"gatk ApplyBQSR "
        f"-I {output_bam} "
        f"-R {reference_genome} "
        f"--bqsr-recal-file {recal_table} "
        f"-O {recalibrated_bam}"
    )
    run_command(cmd_apply, docker_image=docker_gatk)
    logging.info("Recalibration completed.")

def collect_metrics():
    """STEP 4: Collect alignment and insert size metrics."""
    logging.info("Collecting alignment metrics and inserts...")

    # CollectAlignmentSummaryMetrics
    cmd_alignment = (
        f"gatk CollectAlignmentSummaryMetrics "
        f"-R {reference_genome} "
        f"-I {recalibrated_bam} "
        f"-O {alignment_metrics}"
    )
    run_command(cmd_alignment, docker_image=docker_gatk)

    # CollectInsertSizeMetrics
    cmd_insert = (
        f"gatk CollectInsertSizeMetrics "
        f"-I {recalibrated_bam} "
        f"-O {insert_metrics} "
        f"-H {insert_histogram}"
    )
    run_command(cmd_insert, docker_image=docker_gatk)
    logging.info("Metrics collected.")

def call_variants():
    """STEP 5: Calling variants with GATK HaplotypeCaller."""
    logging.info("Calling variants...")
    cmd = (
        f"gatk HaplotypeCaller "
        f"-R {reference_genome} "
        f"-I {recalibrated_bam} "
        f"-O {raw_variants_vcf}"
    )
    run_command(cmd, docker_image=docker_gatk)
    logging.info("Variant calling completed.")

def extract_snps_indels():
    """STEP 6: Split SNPs and INDELs."""
    logging.info("Extracting SNPs and INDELs...")

    cmd_snps = (
        f"gatk SelectVariants "
        f"-R {reference_genome} "
        f"-V {raw_variants_vcf} "
        f"--select-type SNP "
        f"-O {snps_vcf}"
    )
    cmd_indels = (
        f"gatk SelectVariants "
        f"-R {reference_genome} "
        f"-V {raw_variants_vcf} "
        f"--select-type INDEL "
        f"-O {indels_vcf}"
    )
    run_command(cmd_snps, docker_image=docker_gatk)
    run_command(cmd_indels, docker_image=docker_gatk)
    logging.info("Extraction completed.")

def filter_variants():
    """STEP 7: Filtering variants with VariantFiltration."""
    logging.info("Filtering variants...")

    # -- SNPs: Filtering parameters
    cmd_snps = (
        f"gatk VariantFiltration "
        f"-R {reference_genome} "
        f"-V {snps_vcf} "
        f"-O {filtered_snps_vcf} "
        f"-filter \"QD < 2.0\" --filter-name FILTER_QD2 "
        f"-filter \"QUAL < 30.0\" --filter-name FILTER_QUAL30 "
        f"-filter \"SOR > 3.0\" --filter-name FILTER_SOR3 "
        f"-filter \"FS > 60.0\" --filter-name FILTER_FS60 "
        f"-filter \"MQ < 40.0\" --filter-name FILTER_MQ40 "
        f"-filter \"MQRankSum < -12.5\" --filter-name FILTER_MQRankSum12.5 "
        f"-filter \"ReadPosRankSum < -8.0\" --filter-name FILTER_ReadPosRankSum8 "
        f"--genotype-filter-expression \"DP < 10\" --genotype-filter-name GENOTYPE_FILTER_DP10 "
        f"--genotype-filter-expression \"GQ < 10\" --genotype-filter-name GENOTYPE_FILTER_GQ10"
    )
    run_command(cmd_snps, docker_image=docker_gatk)

    # -- INDELs: Filtering parameters
    cmd_indels = (
        f"gatk VariantFiltration "
        f"-R {reference_genome} "
        f"-V {indels_vcf} "
        f"-O {filtered_indels_vcf} "
        f"-filter \"QD < 2.0\" --filter-name FILTER_QD2 "
        f"-filter \"QUAL < 30.0\" --filter-name FILTER_QUAL30 "
        f"-filter \"FS > 200.0\" --filter-name FILTER_FS200 "
        f"-filter \"ReadPosRankSum < -20.0\" --filter-name FILTER_ReadPosRankSum20 "
        f"--genotype-filter-expression \"DP < 10\" --genotype-filter-name GENOTYPE_FILTER_DP10 "
        f"--genotype-filter-expression \"GQ < 10\" --genotype-filter-name GENOTYPE_FILTER_GQ10"
    )
    run_command(cmd_indels, docker_image=docker_gatk)
    logging.info("Filtering completed.")

def select_pass_variants():
    """STEP 8: Select the variants with PASS status (Variant-level)."""
    logging.info("Selecting variants that passed the filters (PASS) at the variant level...")

    cmd_snps = (
        f"gatk SelectVariants "
        f"--exclude-filtered "
        f"-V {filtered_snps_vcf} "
        f"-O {analysis_ready_snps_vcf}"
    )
    cmd_indels = (
        f"gatk SelectVariants "
        f"--exclude-filtered "
        f"-V {filtered_indels_vcf} "
        f"-O {analysis_ready_indels_vcf}"
    )
    run_command(cmd_snps, docker_image=docker_gatk)
    run_command(cmd_indels, docker_image=docker_gatk)
    logging.info("Selection of variants completed.")

def exclude_failed_genotypes():
    """
    STEP 9: Exclude genotypes that failed the filters.
    """
    logging.info("Excluding variants that failed genotype filters...")

    # First copy headers (#), then discard unwanted GENOTYPE FILTER rows.
    cmd_snps = (
        f"grep '^#' {analysis_ready_snps_vcf} > {filtered_snps_gt_vcf} && "
        f"grep -v '^#' {analysis_ready_snps_vcf} | "
        f"grep -v 'GENOTYPE_FILTER_DP10' | "
        f"grep -v 'GENOTYPE_FILTER_GQ10' >> {filtered_snps_gt_vcf}"
    )
    run_command(cmd_snps, docker_image=None)

    cmd_indels = (
        f"grep '^#' {analysis_ready_indels_vcf} > {filtered_indels_gt_vcf} && "
        f"grep -v '^#' {analysis_ready_indels_vcf} | "
        f"grep -v 'GENOTYPE_FILTER_DP10' | "
        f"grep -v 'GENOTYPE_FILTER_GQ10' >> {filtered_indels_gt_vcf}"
    )
    run_command(cmd_indels, docker_image=None)

    logging.info("Exclusion completed.")

def merge_variants():
    """STEP 10: Merge SNP and INDEL into a single VCF."""
    logging.info("Merge SNPs and INDELs...")
    cmd = (
        f"gatk MergeVcfs "
        f"-I {filtered_snps_gt_vcf} "
        f"-I {filtered_indels_gt_vcf} "
        f"-O {merged_vcf}"
    )
    run_command(cmd, docker_image=docker_gatk)
    logging.info("Merge Completed.")

def annotate_with_vep():
    """STEP 11: Functional annotation with VEP within Docker (ensembl-vep)."""
    logging.info("Annotating variants with VEP...")

    # filename final_merge.filtered.vcf" (previous pipeline output)
    # It's located in your local output folder
    final_vcf_local = merged_vcf  # "output/final_merge.filtered.vcf"
    
    # 1) Copy the final_merge.filtered.vcf file to the VEP/cache folder
    # (on the host) without changing the name.
    vep_input_local = os.path.join(vep_data_dir, "final_merge.filtered.vcf")
    shutil.copy(final_vcf_local, vep_input_local)
    logging.info(f"Copied {final_vcf_local} -> {vep_input_local}")

    # 2) Prepare the command to run in Docker
    #    - Mount "/media/gato/22D187F822421089/vep_data" into "/data"
    #    - We will use "/data" as the working directory
    #    - The input and output are in the /data folder inside the container.
    cmd = (
        "vep "
        "--input_file final_merge.filtered.vcf "          # within /data
        "--output_file final_merge.filtered.annotated.vcf "  # generate an annotated file
        "--vcf "
        "--offline "
        "--cache "
        "--force_overwrite "
        "--dir_cache /data "  # <-- We link the cache to /data, where we mount vep data_dir
    )

    # 3) Run VEP with Docker by mounting the VEP folder (vep_data_dir)
    #    in /data to match the --dir_cache /data
    volumes = [
        (os.path.abspath(vep_data_dir), "/data")
    ]
    run_command(cmd, docker_image=docker_vep, volumes=volumes, workdir="/data")

    # 4) Copy the annotated file back to the output folder (output_dir)
    annotated_in_vepdata = os.path.join(vep_data_dir, "final_merge.filtered.annotated.vcf")
    final_annotated_local = os.path.join(output_dir, "final_merge.filtered.annotated.vcf")

    if os.path.exists(annotated_in_vepdata):
        shutil.copy(annotated_in_vepdata, final_annotated_local)
        logging.info(f"Annotated file copied back to {final_annotated_local}")
    else:
        logging.warning("The annotated file was not found in the VEP folder.")

    logging.info("Finished Annotation.")

# =============================
# MAIN
# =============================
def main():
    validate_input_files()
    align_bwa_mem()
    mark_duplicates()
    base_quality_recalibration()
    collect_metrics()
    call_variants()
    extract_snps_indels()
    filter_variants()
    select_pass_variants()
    exclude_failed_genotypes()
    merge_variants()
    annotate_with_vep()
    logging.info("Pipeline successfully completed!")

if __name__ == "__main__":
    main()

