#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
"""
ATACFlow Pipeline - ATAC-seq Quality Control Module

This module provides comprehensive quality control assessment specifically tailored
for ATAC-seq experiments, evaluating key metrics that indicate the success of
the transposase accessibility assay and the quality of the resulting sequencing data.

Key Components:
- get_organelle_filter_expr: Helper function to retrieve organellar chromosome names
- ataqv_qc: Performs comprehensive ATAC-seq QC using the Ataqv tool
- multiqc_ATAC_QC: Aggregates all ATAC-specific QC reports using MultiQC
- multiqc_macs2_samples: Aggregates MACS2 reports from individual samples
- multiqc_macs2_group: Aggregates MACS2 reports including group-level analyses

This module generates quality control metrics that are specific to ATAC-seq,
including TSS enrichment scores, fragment length distributions, nucleosome
positioning patterns, and library complexity estimates, providing a comprehensive
assessment of experiment quality.
"""

def get_organelle_filter_expr(wildcards):
    """
    Retrieve organellar chromosome names (mitochondrial and plastid) from configuration.

    This helper function extracts the names of organellar chromosomes (mitochondria
    and chloroplasts/plastids) from the genome configuration, which are used for
    specialized QC metrics and filtering in ATAC-seq analysis. Organellar reads
    are often highly abundant in ATAC-seq data and require special handling.

    Args:
        wildcards: Snakemake wildcards object containing sample information

    Returns:
        str: Name of the mitochondrial chromosome, or empty string if not configured
    """
    build = config.get("Genome_Version")
    chrMID = config.get("genome_info", {}).get(build, {}).get("chrMID", {})

    if isinstance(chrMID, list):
        # Return the first element if it's a list
        return chrMID[0] if chrMID else ""
    else:
        # Return the value directly if it's not a list
        return chrMID

rule ataqv_qc:
    """
    Perform comprehensive ATAC-seq quality control assessment using Ataqv.

    This rule runs Ataqv (ATAC-seq Quality Validation), a specialized tool designed
    specifically for evaluating the quality of ATAC-seq experiments. Ataqv computes
    a comprehensive set of quality metrics that assess the success of the Tn5
    transposase assay and the overall quality of the sequencing data.

    Key ATAC-seq quality metrics evaluated include:
    - TSS enrichment score: Measures enrichment of reads at transcription start sites
    - Fragment length distribution: Assesses nucleosome positioning patterns
    - Library complexity: Estimates unique molecules and sequencing saturation
    - Mitochondrial read proportion: Evaluates organellar contamination levels
    - Peak calling quality metrics: Assesses the quality of identified peaks
    - Read alignment quality: Evaluates mapping statistics and quality

    Ataqv generates detailed JSON and text output that can be visualized with
    specialized tools, providing a comprehensive assessment of ATAC-seq experiment
    quality and identifying potential issues that may affect downstream analysis.
    This QC step is critical for ensuring that only high-quality ATAC-seq data
    proceeds to peak calling and differential accessibility analysis.
    """
    input:
        shifted_sort_bam = '02.mapping/shifted/{sample}.shifted.sorted.bam',
        shifted_sort_bam_bai = '02.mapping/shifted/{sample}.shifted.sorted.bam.bai',
        narrow_peak = "03.peak_calling/MACS2/{sample}/{sample}_peaks.narrowPeak",
    output:
        json = "02.mapping/ataqv/{sample}.ataqv.json",
        log_out = "02.mapping/ataqv/{sample}.ataqv.out"
    conda:
        workflow.source_path("../envs/ataqv.yaml"),
    resources:
        **rule_resource(config, 'high_resource', skip_queue_on_local=True, logger=logger),
    log:
        "logs/05.qc/ataqv/{sample}.ataqv.log",
    message:
        "Running ataqv QC on {wildcards.sample}",
    params:
        tss = config['Bowtie2_index'][config['Genome_Version']]['tss_bed'],
        autosomes = config['Bowtie2_index'][config['Genome_Version']]['autosomes'],
        mito_name = lambda wildcards:get_organelle_filter_expr(wildcards),
        organism = config['species']
    shell:
        """
        ataqv \
            --peak-file {input.narrow_peak} \
            --name {wildcards.sample} \
            --metrics-file {output.json} \
            --tss-file {params.tss} \
            --autosomal-reference-file {params.autosomes} \
            --mitochondrial-reference-name {params.mito_name} \
            --ignore-read-groups \
            {params.organism} \
            {input.shifted_sort_bam} > {output.log_out} 2> {log}
        """

rule multiqc_ATAC_QC:
    """
    Aggregate all ATAC-seq quality control reports into a unified summary using MultiQC.

    This rule collects all ATAC-seq specific quality control outputs and synthesizes
    them into a single interactive HTML report that provides a comprehensive overview
    of experiment quality across all samples. MultiQC parses outputs from Ataqv,
    Preseq, Samtools, and GATK to present a unified view of ATAC-seq data quality.

    Key QC metrics aggregated in this report include:
    - Ataqv ATAC-seq specific metrics (TSS enrichment, fragment lengths, etc.)
    - Library complexity estimates from Preseq
    - Mapping statistics from Samtools flagstat and stats
    - Duplication metrics from GATK MarkDuplicates
    - Sample-to-sample comparison of all quality metrics

    The aggregated report enables rapid identification of:
    - Outlier samples with quality issues
    - Systematic biases affecting multiple samples
    - Batch effects or plate-based patterns
    - Overall experiment quality and success

    This comprehensive QC summary is essential for quality assurance, providing a
    clear audit trail of data quality for publication and ensuring that only
    high-quality samples proceed to downstream analysis.
    """
    input:
        json = expand("02.mapping/ataqv/{sample}.ataqv.json",sample=samples.keys()),
        log_out = expand("02.mapping/ataqv/{sample}.ataqv.out",sample=samples.keys()),
        preseq = expand('02.mapping/preseq/{sample}.lc_extrap.txt',sample=samples.keys()),
        c_curve = expand('02.mapping/preseq/{sample}.c_curve.txt',sample=samples.keys()),
        samtools_flagstat = expand('02.mapping/samtools_flagstat/{sample}_bam_flagstat.tsv',sample=samples.keys()),
        samtools_stats = expand('02.mapping/samtools_stats/{sample}_bam_stats.tsv',sample=samples.keys()),
        metrics = expand('02.mapping/gatk/{sample}/{sample}.rg.dedup.metrics.txt',sample=samples.keys()),
    output:
        report = '05.ATAC_QC/multiqc_ATAC_report.html',
    resources:
        **rule_resource(config, 'low_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/multiqc.yaml"),
    message:
        "Running MultiQC to aggregate fastp reports",
    benchmark:
        "benchmarks/05.ATAC_QC/multiqc_ATAC_report.txt",
    params:
        origin_reports = "02.mapping/",
        report_dir = "05.ATAC_QC/",
        report = "multiqc_ATAC_report.html",
        title = "ATAC_report",
    log:
        "logs/05.ATAC_QC/multiqc_ATAC_report.log",
    threads:
        config['parameter']['threads']['multiqc'],
    shell:
        """
        multiqc {params.origin_reports} \
                --force \
                --outdir {params.report_dir} \
                -i {params.title} \
                -n {params.report} &> {log}
        """

rule multiqc_macs2_samples:
    """
    Aggregate MACS2 peak calling reports from individual samples using MultiQC.

    This rule collects the MACS2 output from each individual sample and synthesizes
    them into a single interactive HTML report that enables cross-sample comparison
    of peak calling statistics. MultiQC parses the MACS2 XLS output files to
    extract key peak calling metrics and present them in a visual format.

    Key MACS2 metrics aggregated include:
    - Number of peaks called per sample
    - Peak width distributions across samples
    - Fragment length estimates from MACS2 modeling
    - Peak calling statistics (q-values, fold changes, etc.)
    - Comparison of peak calling results across all samples

    The aggregated report enables rapid identification of:
    - Samples with unusually high or low numbers of peaks
    - Systematic differences in peak calling characteristics
    - Outlier samples that may require reanalysis
    - Overall consistency of peak calling across the experiment

    This summary of peak calling results provides valuable quality control information
    about the peak calling step and helps identify any issues that may affect
    downstream differential accessibility analysis.
    """
    input:
        xls = expand("03.peak_calling/MACS2/{sample}/{sample}_peaks.xls",sample=samples.keys())
    output:
        report = '05.ATAC_QC/multiqc_MACS2_Samples_report.html',
    resources:
        **rule_resource(config, 'low_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/multiqc.yaml"),
    message:
        "Running MultiQC to aggregate MACS2",
    benchmark:
        "benchmarks/05.ATAC_QC/multiqc_MACS2_Samples_report.txt",
    params:
        origin_reports = "03.peak_calling/",
        report_dir = "05.ATAC_QC/",
        report = "multiqc_MACS2_Samples_report.html",
        title = "MACS2_Samples_report",
    log:
        "logs/05.ATAC_QC/multiqc_MACS2_Samples_report.log",
    threads:
        config['parameter']['threads']['multiqc'],
    shell:
        """
        multiqc {params.origin_reports} \
                --force \
                --outdir {params.report_dir} \
                -i {params.title} \
                -n {params.report} &> {log}
        """

rule multiqc_macs2_group:
    """
    Aggregate MACS2 peak calling reports including both individual and group-level analyses.

    This rule collects MACS2 outputs from both individual samples and merged group
    analyses, synthesizing them into a comprehensive report that enables comparison
    of peak calling results at both levels. This is particularly valuable for
    experiments with biological replicates where both sample-level and group-level
    peak calling are performed.

    Key features of this comprehensive report:
    - Includes peak calling statistics from all individual samples
    - Includes peak calling statistics from all merged group analyses
    - Enables comparison of individual vs. group-level peak calling results
    - Shows how merging replicates affects peak calling statistics

    The aggregated report helps assess:
    - The benefit of merging replicates for peak calling sensitivity
    - Consistency between individual and group-level results
    - Quality of both sample-level and group-level peak calls
    - Overall robustness of the peak calling strategy

    This comprehensive view of peak calling results at multiple levels provides
    valuable insights into the experiment and helps validate the analytical approach.
    """
    input:
        xls =  expand("03.peak_calling/MACS2/{sample}/{sample}_peaks.xls",sample=samples.keys()),
        group_xls = expand("03.peak_calling/MERGE_MACS2/{group}/{group}_peaks.xls",group = groups.keys()),
    output:
        report = '05.ATAC_QC/multiqc_MACS2_Merge_report.html',
    resources:
        **rule_resource(config, 'low_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/multiqc.yaml"),
    message:
        "Running MultiQC to aggregate MACS2",
    benchmark:
        "benchmarks/05.ATAC_QC/multiqc_MACS2_Merge_report.txt",
    params:
        origin_reports = "03.peak_calling/",
        report_dir = "05.ATAC_QC/",
        report = "multiqc_MACS2_Merge_report.html",
        title = "MACS2_Merge_report",
    log:
        "logs/05.ATAC_QC/multiqc_MACS2_Merge_report.log",
    threads:
        config['parameter']['threads']['multiqc'],
    shell:
        """
        multiqc {params.origin_reports} \
                --force \
                --outdir {params.report_dir} \
                -i {params.title} \
                -n {params.report} &> {log}
        """
# ----- end of rules ----- #