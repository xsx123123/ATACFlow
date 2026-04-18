#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
"""
ATACFlow Pipeline - Group-level Differential Accessibility Analysis Module

This module performs differential accessibility analysis using the consensus peakset
derived from merged group-level BAM files, providing a more robust approach to
identifying chromatin regions that show significant changes in accessibility between
experimental conditions. By using group-level peaks as the consensus set, this
approach can capture more consistent regulatory changes.

Key Components:
- merge_DEG: Performs differential accessibility analysis using DESeq2 on group-level peaks
- merge_Enrichments: Performs gene ontology enrichment analysis on group-level differential peaks

This module is particularly valuable for experiments with multiple biological replicates,
where the group-level peak calling approach provides a more robust foundation for
differential accessibility analysis.
"""

rule merge_DEG:
    """
    Perform differential accessibility analysis using DESeq2 on the group-level consensus peakset.

    This rule identifies chromatin regions that show statistically significant differences
    in accessibility between experimental conditions using DESeq2, but uses the consensus
    peakset derived from merging group-level MACS2 results rather than individual sample
    peaks. This approach provides a more robust set of peaks for differential analysis
    by focusing on consistent accessible regions within each experimental condition.

    Key advantages of the group-level approach:
    - Reduces noise by focusing on peaks consistently called within groups
    - Provides a more comprehensive view of the regulatory landscape
    - Can capture weak but consistent accessibility differences
    - Reduces multiple testing burden by using a smaller, more robust peakset

    Key steps in the differential accessibility analysis:
    - Normalization of read counts across all samples
    - Estimation of dispersion parameters using the group-level peakset
    - Statistical testing for differential accessibility between conditions
    - Multiple testing correction using FDR (False Discovery Rate)
    - Generation of visualizations including PCA plots for quality assessment

    Important analysis parameters:
    - LFC (log2 fold change) threshold for calling significant differences
    - P-value (adjusted) threshold for statistical significance
    - Sample information for defining the experimental design
    - Contrast information for specifying which comparisons to perform

    Outputs from this analysis include:
    - Global PCA plot showing overall sample relationships and clustering
    - Detailed results tables for each contrast with comprehensive statistics
    - Summary statistics of differential peaks across all comparisons
    - Various diagnostic plots and quality control visualizations

    Using the group-level consensus peakset for differential analysis often provides
    more biologically meaningful results by focusing on the consistent accessible
    regions that are most likely to represent genuine regulatory elements.
    """
    input:
        counts = lambda wildcards: get_diff_analysis_input(config, config.get('_merge_group', False)),
    output:
        output = '06.deg_enrich/DEG_merge/Global_PCA.pdf',
        summary = '06.deg_enrich/DEG_merge/All_Contrast_Differential_Peaks_Statistics.csv',
    resources:
        **rule_resource(config, 'low_resource',  skip_queue_on_local=True,logger = logger),
    message:
        "Running merge_DEG",
    conda:
        workflow.source_path("../envs/deg_deseq2.yaml"),
    log:
        "logs/06.DEG/merge_DEG.log",
    benchmark:
        "benchmarks/06.DEG/merge_DEG.txt",
    params:
        deg_dir = '06.deg_enrich/DEG_merge',
        samples = config['sample_csv'],
        paired = config['paired_csv'],
        PATH = workflow.source_path(config['parameter']['DEG']['PATH']),
        LFC = config['parameter']['DEG']['LFC'],
        PVAL = config['parameter']['DEG']['PVAL'],
    threads: 
        1
    shell:
        """
        chmod +x {params.PATH} && \
        Rscript {params.PATH} -c {input.counts} \
                -m {params.samples} \
                -p {params.paired} \
                -o {params.deg_dir} \
                --lfc={params.LFC} \
                --pval={params.PVAL} \
                --label_co=gene_id &> {log}
        """

rule merge_Enrichments:
    """
    Perform gene ontology enrichment analysis on group-level differentially accessible peaks.

    This rule performs functional enrichment analysis to identify biological processes,
    molecular functions, and cellular components that are significantly associated with
    the genes near differentially accessible chromatin regions identified using the
    group-level consensus peakset. This analysis links changes in chromatin accessibility
    to potential functional consequences at the gene and pathway level, using the more
    robust set of peaks derived from group-level analysis.

    Key enrichment analysis steps using group-level peaks:
    - Assignment of group-identified peaks to putative target genes based on genomic proximity
    - Statistical testing for enrichment of gene ontology terms in the associated genes
    - Multiple testing correction to control the false discovery rate
    - Generation of comprehensive enrichment results tables and visualizations

    Important parameters for enrichment analysis:
    - Gene ontology (GO) database structure (OBO file)
    - Gene ontology annotations specific to the organism being studied
    - Significance cutoff for enrichment (default: 0.05 adjusted p-value)
    - Method for assigning peaks to potential target genes
    - Regular expression pattern for extracting gene identifiers from annotations

    The group-level enrichment analysis provides several advantages:
    - Higher confidence in peak-to-gene assignments due to more robust peak calls
    - More reliable enrichment results due to reduced noise in the peakset
    - Better ability to detect subtle but consistent regulatory changes
    - More biologically meaningful functional annotations

    The results from this analysis, together with the differential accessibility results
    from the group-level analysis, provide a comprehensive view of the regulatory changes
    occurring between experimental conditions, enabling deeper and more reliable biological
    interpretation of the ATAC-seq data.
    """
    input:
        DEG_info = "06.deg_enrich/DEG_merge/All_Contrast_Differential_Peaks_Statistics.csv",
    output:
        Enrichments_dir = directory("06.deg_enrich/merge_enrich/"),
    resources:
        **rule_resource(config, 'low_resource',  skip_queue_on_local=True,logger = logger),
    benchmark:
        "benchmarks/06.DEG/merge_Enrichments.txt",
    message:
        "Running merge_Enrichments",
    conda:
        workflow.source_path("../envs/go_enrich_r.yaml"),
    log:
        "logs/06.DEG/merge_Enrichments.log",
    params:
        obo = config['Bowtie2_index']['GO']['obo'],
        go_annotation = config['Bowtie2_index'][config['Genome_Version']]['go_annotation'],
        gene_col = config['parameter']['Enrichments']['gene_col'],
        r_script = workflow.source_path(config['parameter']['Enrichments']['PATH']),
        wrapper = workflow.source_path(config['parameter']['Enrichments']['PATH_py']),
        lib_type = config['parameter']['Enrichments']['lib_type'],
        gene_regex = config['parameter']['Enrichments']['gene_regex'],
        deg_dir = "06.deg_enrich/DEG_merge",
        cutoff = config['parameter']['Enrichments'].get('cutoff', 0.05),
    threads:
        1
    shell:
        """
        python {params.wrapper} \
            --rscript {params.r_script} \
            --deg_info {input.DEG_info} \
            --deg_dir {params.deg_dir} \
            --lib_type {params.lib_type} \
            -o {params.obo} \
            -a {params.go_annotation} \
            -d {output.Enrichments_dir} \
            --gene_col {params.gene_col} \
            --gene_regex '{params.gene_regex}' \
            --cutoff {params.cutoff} > {log} 2>&1
        """
# ------- rule ------- #