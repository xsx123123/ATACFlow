#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
import os

rule macs2_callpeak:
    """
    Call peaks using MACS2 on the shifted BAM file.
    Note: Since the BAM is already shifted, we use BAMPE mode without additional shift parameters.
    """
    input:
        bam = '02.mapping/shifted/{sample}.shifted.sorted.bam',
        bai = '02.mapping/shifted/{sample}.shifted.sorted.bam.bai'
    output:
        narrow_peak = "03.peak_calling/MACS2/{sample}/{sample}_peaks.narrowPeak",
        xls = "03.peak_calling/MACS2/{sample}/{sample}_peaks.xls",
        summits = "03.peak_calling/MACS2/{sample}/{sample}_summits.bed",
        bdg = "03.peak_calling/MACS2/{sample}/{sample}_treat_pileup.bdg"
    resources:
        **rule_resource(config, 'high_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/macs2.yaml"),
    log:
        "logs/03.peak_calling/macs2_{sample}.log",
    message:
        "Running MACS2 peak calling for {wildcards.sample}",
    benchmark:
        "benchmarks/03.peak_calling/macs2_{sample}.txt",
    threads:
        1
    params:
        gsize = config['genome_info'][config['Genome_Version']]['effectiveGenomeSize'],
        qvalue = config['parameter']['macs2']['qvalue'],
        outdir = "03.peak_calling/MACS2/{sample}"
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
        peak = "03.peak_calling/MACS2/{sample}/{sample}_peaks.narrowPeak"
    output:
        annotation = "03.peak_calling/HOMER/{sample}_annotation.txt",
        stats = "03.peak_calling/HOMER/{sample}_stats.txt"
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/homer.yaml"),
    log:
        "logs/03.peak_calling/homer_{sample}.log",
    message:
        "Running HOMER annotation for {wildcards.sample}",
    benchmark:
        "benchmarks/03.peak_calling/homer_{sample}.txt",
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
    Filter: strict sorting and merging overlapping regions.
    """
    input:
        peaks = expand("03.peak_calling/MACS2/{sample}/{sample}_peaks.narrowPeak", sample=samples.keys())
    output:
        consensus = "04.consensus/consensus_peaks.bed"
    resources:
        **rule_resource(config, 'high_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/bedtools.yaml"),
    log:
        "logs/04.consensus/merge_peaks.log",
    message:
        "Creating consensus peakset from all samples",
    benchmark:
        "benchmarks/04.consensus/merge_peaks.txt",
    threads: config['parameter']['threads']['bedtools']
    shell:
        """
        cat {input.peaks} | \
        sort -k1,1 -k2,2n | \
        bedtools merge -i stdin > {output.consensus} 2> {log}
        """

rule generate_count_matrix:
    """
    Create a tabular file (Count Matrix) using BEDTools multicov.
    Rows = Consensus Peaks, Columns = Samples.
    This file is ready for filtering and DE analysis (DESeq2/edgeR).
    """
    input:
        consensus = "04.consensus/consensus_peaks.bed",
        bams = expand("02.mapping/shifted/{sample}.shifted.sorted.bam", sample=samples.keys()),
        bais = expand("02.mapping/shifted/{sample}.shifted.sorted.bam.bai", sample=samples.keys())
    output:
        counts = "04.consensus/raw_counts.txt",
        counts_with_header = "04.consensus/consensus_counts_matrix.txt"
    resources:
        **rule_resource(config, 'high_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/bedtools.yaml"),
    log:
        "logs/04.consensus/multicov.log",
    message:
        "Generating count matrix from consensus peaks and BAM files",
    benchmark:
        "benchmarks/04.consensus/multicov.txt",
    threads: config['parameter']['threads']['bedtools']
    params:
        sample_names = " ".join(samples.keys())
    shell:
        """
        echo "Running multicov..." > {log}
        bedtools multicov -bams {input.bams} -bed {input.consensus} > {output.counts}.tmp 2>> {log}
        HEADER="chrom\tstart\tend\t"{params.sample_names}
        echo -e "${{HEADER}}" | tr ' ' '\t' | cat - {output.counts}.tmp > {output.counts}
        echo -e "${{HEADER}}" | tr ' ' '\t' | cat - {output.counts}.tmp > {output.counts_with_header}
        rm {output.counts}.tmp
        """
# ----- end of rules ----- #