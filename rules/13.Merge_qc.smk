#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
"""
ATACFlow Pipeline - Merged Quality Control Report Module

This module aggregates all preprocessing quality control reports into a single
comprehensive MultiQC report, providing a unified view of data quality across
all preprocessing steps from raw data through cleaning and initial processing.

Key Components:
- merge_qc_report: Aggregates all preprocessing QC reports using MultiQC

This module generates a final preprocessing quality control report that combines
results from FastQC, Fastp, and other quality assessment tools, enabling quick
evaluation of overall data quality and identification of any systematic issues
or outlier samples.
"""

import os

# ----- rule ----- #
rule merge_qc_report:
    """
    Aggregate all preprocessing quality control reports into a comprehensive summary.

    This rule collects all quality control outputs from the preprocessing pipeline
    and synthesizes them into a single interactive HTML report using MultiQC.
    This comprehensive report provides a unified view of data quality across all
    samples and all preprocessing steps, from raw data through cleaning and
    initial quality assessment.

    Key quality control metrics aggregated in this report include:
    - Raw read quality statistics from FastQC (both R1 and R2 reads)
    - Trimming statistics and quality improvements from Fastp
    - Adapter content and trimming efficiency
    - Read length distributions after quality filtering
    - Overall quality metrics before and after preprocessing
    - Sample-to-sample comparison of all quality metrics

    The aggregated report enables rapid identification of:
    - Outlier samples with quality issues at any preprocessing step
    - Systematic biases affecting multiple samples
    - Batch effects or plate-based quality patterns
    - Overall success of the preprocessing pipeline
    - Trends in quality improvements through the pipeline

    This comprehensive preprocessing QC summary is essential for quality assurance,
    providing a clear audit trail of data quality improvements and ensuring that
    only high-quality, properly processed data proceeds to downstream analysis
    steps including read alignment and peak calling. The report also serves as
    valuable documentation for publication and reproducibility.
    """
    input:
        md5_check = "01.qc/md5_check.tsv",
        fastqc_files_r1 = expand("01.qc/short_read_qc_r1/{sample}_R1_fastqc.zip", sample=samples.keys()),
        fastqc_files_r2 = expand("01.qc/short_read_qc_r2/{sample}_R2_fastqc.zip", sample=samples.keys()),
        r1_trimmed = expand("01.qc/short_read_trim/{sample}.R1.trimed.fq.gz", sample=samples.keys()),
        r2_trimmed = expand("01.qc/short_read_trim/{sample}.R2.trimed.fq.gz", sample=samples.keys()),
        fastp_report = expand("01.qc/short_read_trim/{sample}.trimed.html", sample=samples.keys()),
    output:
        report = '01.qc/multiqc_merge_qc/multiqc_merge_qc_report.html',
    resources:
        **rule_resource(config, 'low_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/multiqc.yaml"),
    params:
        search_dir = "01.qc",
        out_dir = "01.qc/multiqc_merge_qc",
        report_name = "multiqc_merge_qc_report.html",
        title = "merge qc report",
    log:
        "logs/01.qc/multiqc_merge_qc.log",
    threads:
        config['parameter']['threads']['multiqc'],
    shell:
        """
        mkdir -p {params.out_dir}
        multiqc ./{params.search_dir} \
                --force \
                --outdir {params.out_dir} \
                --title "{params.title}" \
                --filename {params.report_name} \
                &> {log}
        """
# ----- rule ----- #