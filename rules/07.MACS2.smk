#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
"""
ATACFlow Pipeline - Single Sample Peak Calling Module

This module handles individual sample-level peak calling and consensus peak
generation from all samples. It provides the foundational peak calling analysis
for ATAC-seq experiments.

Key Components:
- macs2_callpeak: Single sample peak calling using MACS2
- homer_annotate_peaks: Annotate single sample peaks using HOMER
- create_consensus_peakset: Create consensus peaks from ALL samples

Output Structure:
- 03.peak_calling/single/{sample}/ - Individual sample MACS2 results
- 03.peak_calling/single_HOMER/ - Individual sample HOMER annotations
- 04.consensus/single/ - Consensus peaks from all samples
"""

rule macs2_callpeak:
    """
    Identify open chromatin regions (peaks) using MACS2 on the shifted ATAC-seq BAM file.

    This rule performs peak calling on each individual sample to identify genomic
    regions with significant chromatin accessibility. MACS2 is used with parameters
    optimized for ATAC-seq data to detect regions where Tn5 transposase insertion
    is significantly enriched.

    Key parameters:
    - -f BAMPE: Uses paired-end information for improved peak modeling
    - -g {{genome_size}}: Effective genome size for p-value calculation
    - -q {{qvalue}}: Q-value threshold (FDR control)
    - --keep-dup all: Retains all duplicates (already filtered)
    - -B --SPMR: Generates normalized bedGraph files
    """
    input:
        bam = '02.mapping/shifted/{sample}.shifted.sorted.bam',
        bai = '02.mapping/shifted/{sample}.shifted.sorted.bam.bai'
    output:
        narrow_peak = "03.peak_calling/single/{sample}/{sample}_peaks.narrowPeak",
        xls = "03.peak_calling/single/{sample}/{sample}_peaks.xls",
        summits = "03.peak_calling/single/{sample}/{sample}_summits.bed",
        bdg = "03.peak_calling/single/{sample}/{sample}_treat_pileup.bdg"
    resources:
        **rule_resource(config, 'high_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/macs2.yaml"),
    log:
        "logs/03.peak_calling/single/macs2_{sample}.log",
    message:
        "Running MACS2 peak calling for {wildcards.sample}",
    benchmark:
        "benchmarks/03.peak_calling/single/macs2_{sample}.txt",
    threads:
        config['parameter']['threads']['macs2'],
    params:
        gsize = config['genome_info'][config['Genome_Version']]['effectiveGenomeSize'],
        qvalue = config['parameter']['macs2']['qvalue'],
        outdir = "03.peak_calling/single/{sample}"
    shell:
        """
        macs2 callpeak \
            -t {input.bam} \
            -f BAMPE \
            -g {params.gsize} \
            --name {wildcards.sample} \
            --outdir {params.outdir} \
            -q {params.qvalue} \
            --keep-dup all \
            -B --SPMR 2> {log}
        """

rule homer_annotate_peaks:
    """
    Annotate peaks relative to gene features using HOMER.
    """
    input:
        peak = "03.peak_calling/single/{sample}/{sample}_peaks.narrowPeak"
    output:
        annotation = "03.peak_calling/single_HOMER/{sample}_annotation.txt",
        stats = "03.peak_calling/single_HOMER/{sample}_stats.txt"
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/homer.yaml"),
    log:
        "logs/03.peak_calling/single/homer_{sample}.log",
    message:
        "Running HOMER annotation for {wildcards.sample}",
    benchmark:
        "benchmarks/03.peak_calling/single/homer_{sample}.txt",
    threads: config['parameter']['threads']['homer']
    params:
        gtf = config['Bowtie2_index'][config['Genome_Version']]['genome_gtf'],
        genome_fasta = config['Bowtie2_index'][config['Genome_Version']]['genome_fa'],
    shell:
        """
        annotatePeaks.pl {input.peak} {params.genome_fasta} \
            -gtf {params.gtf} \
            -annStats {output.stats} \
            -p {threads} > {output.annotation} 2> {log}
        """

rule create_consensus_peakset:
    """
    Merge peaks from ALL samples to create a single consensus BED file.
    
    This rule combines narrowPeak files from all individual samples, sorts them,
    and merges overlapping regions to create a comprehensive consensus peakset
    representing the union of all accessible chromatin regions across the dataset.
    """
    input:
        peaks = expand("03.peak_calling/single/{sample}/{sample}_peaks.narrowPeak", sample=samples.keys())
    output:
        consensus = "04.consensus/single/all_samples_consensus_peaks.bed"
    resources:
        **rule_resource(config, 'high_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/bedtools.yaml"),
    log:
        "logs/04.consensus/single/merge_peaks.log",
    message:
        "Creating consensus peakset from all samples",
    benchmark:
        "benchmarks/04.consensus/single/merge_peaks.txt",
    threads: 
        config['parameter']['threads']['bedtools']
    shell:
        """
        cat {input.peaks} | \
        sort -k1,1 -k2,2n | \
        bedtools merge -i stdin > {output.consensus} 2> {log}
        """

rule homer_annotate_consensus_peaks:
    """
    Annotate consensus peaks relative to gene features using HOMER.
    """
    input:
        consensus = "04.consensus/single/all_samples_consensus_peaks.bed"
    output:
        annotation = "04.consensus/single/consensus_peaks_annotation.txt",
        stats = "04.consensus/single/consensus_peaks_stats.txt"
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/homer.yaml"),
    log:
        "logs/04.consensus/single/consensus_annotate.log",
    message:
        "Running HOMER annotation for consensus peaks",
    benchmark:
        "benchmarks/04.consensus/single/consensus_annotate.txt",
    threads: config['parameter']['threads']['homer']
    params:
        gtf = config['Bowtie2_index'][config['Genome_Version']]['genome_gtf'],
        genome_fasta = config['Bowtie2_index'][config['Genome_Version']]['genome_fa'],
    shell:
        """
        annotatePeaks.pl {input.consensus} {params.genome_fasta} \
            -gtf {params.gtf} \
            -annStats {output.stats} \
            -p {threads} > {output.annotation} 2> {log}
        """

rule generate_count_matrix_by_featureCounts:
    """
    Generate a read count matrix using featureCounts for the consensus peakset.

    This rule counts ATAC-seq paired-end fragments overlapping each consensus peak
    across all samples. The -p flag ensures proper fragment-level counting.
    """
    input:
        consensus = "04.consensus/single/all_samples_consensus_peaks.bed",
        bams = expand("02.mapping/shifted/{sample}.shifted.sorted.bam", sample=samples.keys())
    output:
        counts_matrix = "04.consensus/single/consensus_counts_matrix.txt",
        description = "04.consensus/single/matrix_description.txt",
        summary = "04.consensus/single/consensus_counts_matrix.txt.summary",
        saf = temp("04.consensus/single/consensus_peaks.saf") 
    conda:
        workflow.source_path("../envs/subread.yaml"),
    log:
        "logs/04.consensus/single/featureCounts.log"
    benchmark:
        "benchmarks/04.consensus/single/featureCounts.txt"
    resources:
        **rule_resource(config, 'high_resource', skip_queue_on_local=True, logger=logger),
    threads: 
        config['parameter']['threads'].get('featurecounts', 16)
    shell:
        """
        echo "1. Converting consensus BED to SAF format..." > {log}
        awk 'BEGIN{{OFS="\\t"; print "GeneID\\tChr\\tStart\\tEnd\\tStrand"}} \
            {{print $1":"$2"-"$3, $1, $2, $3, "+"}}' {input.consensus} > {output.saf}

        echo "2. Running featureCounts for ATAC-seq PE fragments..." >> {log}
        featureCounts \
            -p \
            -B \
            -C \
            -T {threads} \
            -F SAF \
            -a {output.saf} \
            -o {output.counts_matrix} \
            {input.bams} >> {log} 2>&1
            
        echo "3. Cleaning up matrix header for downstream R compatibility..." >> {log}
        sed -i 's|02.mapping/shifted/||g; s|.shifted.sorted.bam||g' {output.counts_matrix}

        echo "4. Generating description file..." >> {log}
        echo "File Name: $(basename {output.counts_matrix})" > {output.description}
        echo "Generated Date: $(date +'%Y-%m-%d %H:%M:%S')" >> {output.description}
        echo "--------------------------------------------------" >> {output.description}
        echo "Total Samples: $(echo "{input.bams}" | wc -w)" >> {output.description}
        echo "Total Consensus Peaks: $(wc -l < {input.consensus})" >> {output.description}
        """

rule generate_count_matrix_ann:
    """
    Add gene annotation to the count matrix.
    """
    input:
        annotation = "04.consensus/single/consensus_peaks_annotation.txt",
        counts_matrix = "04.consensus/single/consensus_counts_matrix.txt",
    output:
        counts_matrix_ann = "04.consensus/single/consensus_counts_matrix_ann.txt",
    conda:
        workflow.source_path("../envs/bedtools.yaml"),
    resources:
        **rule_resource(config, 'low_resource', skip_queue_on_local=True, logger=logger),
    log:
        "logs/04.consensus/single/consensus_counts_matrix_ann.log",
    benchmark:
        "benchmarks/04.consensus/single/consensus_counts_matrix_ann.txt",
    params:
        path = workflow.source_path(config['parameter']['merge_peaks']['path'])
    threads: 
        1
    shell:
        """
        chmod +x {params.path}
        python3 {params.path} \
            -a {input.annotation} \
            -c {input.counts_matrix} \
            -o {output.counts_matrix_ann} &> {log}
        """