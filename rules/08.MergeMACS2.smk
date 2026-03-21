#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
"""
ATACFlow Pipeline - Group-Level Analysis Module

This module handles group-level peak calling and IDR analysis for ATAC-seq data.
It provides two distinct analysis workflows:

1. POOLED ANALYSIS: Merge replicates within each group and call peaks
2. IDR ANALYSIS: Identify reproducible peaks across replicates using IDR

Output Structure:
- 03.peak_calling/pooled/{group}/ - Group-level pooled MACS2 results
- 03.peak_calling/pooled_HOMER/ - Group-level HOMER annotations
- 03.peak_calling/idr/{group}/ - IDR analysis results
- 03.peak_calling/idr_HOMER/ - IDR peak annotations
- 04.consensus/pooled/ - Group consensus peaks
- 04.consensus/idr/ - IDR-filtered high-confidence peaks
"""

import os

rule merge_shifted_bams:
    """
    Merge shifted ATAC-seq BAM files by experimental group.
    
    This rule combines BAM files from all samples belonging to the same
    experimental group, creating a single consolidated BAM file with increased
    sequencing depth for more robust peak calling.
    """
    input:
        lambda wildcards: expand("02.mapping/shifted/{group}.shifted.sorted.bam",group=groups[wildcards.group])
    output:
        bam = "02.mapping/merged/{group}.merged.bam",
        bai = "02.mapping/merged/{group}.merged.bam.bai"
    log:
        "logs/02.mapping/merge_{group}.log"
    threads:
        8
    conda:
        workflow.source_path("../envs/samtools.yaml")
    shell:
        """
        echo "Merging group {wildcards.group}..." > {log}
        samtools merge -@ {threads} {output.bam} {input} 2>> {log}
        samtools index -@ {threads} {output.bam} 2>> {log}
        """

rule macs2_pooled_callpeak:
    """
    Identify peaks using MACS2 on merged group BAM files (Pooled Analysis).
    
    This performs peak calling on merged BAM files from each experimental group,
    leveraging increased sequencing depth to identify more robust accessible
    chromatin regions.
    """
    input:
        bam = "02.mapping/merged/{group}.merged.bam",
        bai = "02.mapping/merged/{group}.merged.bam.bai"
    output:
        narrow_peak = "03.peak_calling/pooled/{group}/{group}_peaks.narrowPeak",
        xls = "03.peak_calling/pooled/{group}/{group}_peaks.xls",
        summits = "03.peak_calling/pooled/{group}/{group}_summits.bed",
        bdg = "03.peak_calling/pooled/{group}/{group}_treat_pileup.bdg"
    resources:
        **rule_resource(config, 'high_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/macs2.yaml"),
    log:
        "logs/03.peak_calling/pooled/macs2_{group}.log",
    message:
        "Running pooled MACS2 peak calling for {wildcards.group}",
    benchmark:
        "benchmarks/03.peak_calling/pooled/macs2_{group}.txt",
    threads:
        config['parameter']['threads']['macs2'],
    params:
        gsize = config['genome_info'][config['Genome_Version']]['effectiveGenomeSize'],
        qvalue = config['parameter']['macs2']['qvalue'],
        outdir = "03.peak_calling/pooled/{group}"
    shell:
        """
        macs2 callpeak \
            -t {input.bam} \
            -f BAMPE \
            -g {params.gsize} \
            --name {wildcards.group} \
            --outdir {params.outdir} \
            -q {params.qvalue} \
            --keep-dup all \
            -B --SPMR 2> {log}
        """

rule homer_annotate_pooled_peaks:
    """
    Annotate group-level pooled peaks using HOMER.
    """
    input:
        peak = "03.peak_calling/pooled/{group}/{group}_peaks.narrowPeak",
    output:
        annotation = "03.peak_calling/pooled_HOMER/{group}_annotation.txt",
        stats = "03.peak_calling/pooled_HOMER/{group}_stats.txt"
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/homer.yaml"),
    log:
        "logs/03.peak_calling/pooled/homer_{group}.log",
    message:
        "Running HOMER annotation for {wildcards.group} pooled peaks",
    benchmark:
        "benchmarks/03.peak_calling/pooled/homer_{group}.txt",
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

rule create_pooled_consensus_peakset:
    """
    Create consensus peaks from group-level pooled peak calls.
    
    Merges peaks from all experimental groups to create a comprehensive
    consensus peakset for differential accessibility analysis.
    """
    input:
        peaks = expand("03.peak_calling/pooled/{group}/{group}_peaks.narrowPeak", group=groups.keys())
    output:
        consensus = "04.consensus/pooled/all_groups_consensus_peaks.bed"
    resources:
        **rule_resource(config, 'high_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/bedtools.yaml"),
    log:
        "logs/04.consensus/pooled/merge_peaks.log",
    message:
        "Creating consensus peakset from pooled group peaks",
    benchmark:
        "benchmarks/04.consensus/pooled/merge_peaks.txt",
    threads: config['parameter']['threads']['bedtools']
    shell:
        """
        cat {input.peaks} | \
        sort -k1,1 -k2,2n | \
        bedtools merge -i stdin > {output.consensus} 2> {log}
        """

rule homer_annotate_pooled_consensus_peaks:
    """
    Annotate pooled consensus peaks using HOMER.
    """
    input:
        consensus = "04.consensus/pooled/all_groups_consensus_peaks.bed"
    output:
        annotation = "04.consensus/pooled/consensus_peaks_annotation.txt",
        stats = "04.consensus/pooled/consensus_peaks_stats.txt"
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/homer.yaml"),
    log:
        "logs/04.consensus/pooled/consensus_annotate.log",
    message:
        "Running HOMER annotation for pooled consensus peaks",
    benchmark:
        "benchmarks/04.consensus/pooled/consensus_annotate.txt",
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

rule idr_analysis:
    """
    Perform Irreproducible Discovery Rate (IDR) analysis on group replicates.
    
    Identifies reproducible peaks across biological replicates within a group
    using the IDR framework. For groups with only one replicate, returns the
    original peaks as consensus.
    
    NOTE: IDR-filtered peaks are used for QC reporting, but differential analysis
    uses the pooled consensus peaks for sufficient coverage.
    """
    input:
        peaks = lambda wildcards: expand("03.peak_calling/single/{sample}/{sample}_peaks.narrowPeak",
            sample=groups[wildcards.group]
        )
    output:
        idr_bed = "03.peak_calling/idr/{group}/Final_Consensus_Peaks.bed",
        idr_log = "03.peak_calling/idr/{group}/idr_pipeline.log",
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/idr.yaml"),
    log:
        "logs/03.peak_calling/idr/{group}.log"
    benchmark:
        "benchmarks/03.peak_calling/idr/{group}.txt"
    threads:
        config['parameter']['threads'].get('idr', 4)
    params:
        script = workflow.source_path(config["parameter"]["idr_scripts"]["path"]),
        outdir = "03.peak_calling/idr/{group}"
    shell:
        """
        mkdir -p {params.outdir}

        peak_count=$(echo "{input.peaks}" | wc -w)

        if [ "$peak_count" -ge 2 ]; then
            echo "[$(date)] Running IDR for group {wildcards.group} with $peak_count replicates..." > {log}
            python3 {params.script} \
                -i {input.peaks} \
                -o {params.outdir} \
                -t {threads} >> {log} 2>&1
        else
            echo "[$(date)] Group {wildcards.group} has only one replicate. Using original peaks..." >> {log}
            awk 'BEGIN{{OFS="\\t"}} {{print $1, $2, $3}}' {input.peaks} | \
            sort -k1,1 -k2,2n | \
            bedtools merge > {output.idr_bed} 2>> {log}
        fi
        """

rule homer_annotate_idr_peaks:
    """
    Annotate IDR-filtered consensus peaks using HOMER.
    """
    input:
        idr_peaks = "03.peak_calling/idr/{group}/Final_Consensus_Peaks.bed",
    output:
        annotation = "03.peak_calling/idr_HOMER/{group}_annotation.txt",
        stats = "03.peak_calling/idr_HOMER/{group}_stats.txt",
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/homer.yaml"),
    log:
        "logs/03.peak_calling/idr/homer_{group}.log",
    message:
        "Running HOMER annotation for {wildcards.group} IDR peaks",
    benchmark:
        "benchmarks/03.peak_calling/idr/homer_{group}.txt",
    params:
        gtf = config['Bowtie2_index'][config['Genome_Version']]['genome_gtf'],
        genome_fasta = config['Bowtie2_index'][config['Genome_Version']]['genome_fa'],
    threads: 
        config['parameter']['threads']['homer'],
    shell:
        """
        annotatePeaks.pl {input.idr_peaks} {params.genome_fasta} \
            -gtf {params.gtf} \
            -annStats {output.stats} \
            -p {threads} > {output.annotation} 2> {log}
        """

rule calculate_frip_score:
    """
    Calculate FRiP (Fraction of Reads in Peaks) score for each sample.
    
    FRiP measures the proportion of reads that fall within called peaks,
    which is a quality metric for peak calling performance.
    """
    input:
        bam = "02.mapping/shifted/{sample}.shifted.sorted.bam",
        peaks = "03.peak_calling/single/{sample}/{sample}_peaks.narrowPeak"
    output:
        frip = "03.peak_calling/single/{sample}/{sample}_frip.txt"
    conda:
        workflow.source_path("../envs/subread.yaml"),
    resources:
        **rule_resource(config, 'low_resource', skip_queue_on_local=True, logger=logger),
    log:
        "logs/03.peak_calling/single/frip_{sample}.log"
    threads: 4
    shell:
        """
        totalReads=$(samtools view -c {input.bam})
        
        awk 'BEGIN{{OFS="\\t"; print "GeneID\\tChr\\tStart\\tEnd\\tStrand"}}{{print $4, $1, $2+1, $3, "."}}' {input.peaks} > {output.frip}.saf
        
        peakReads=$(featureCounts -p -a {output.frip}.saf -F SAF -o {output.frip}.counts {input.bam} 2>> {log} | grep -A1 "TotalReads" | tail -n1 | awk '{{print $1}}')
        
        frip=$(echo "scale=4; $peakReads / $totalReads" | bc)
        
        echo "TotalReads: $totalReads" > {output.frip}
        echo "PeakReads: $peakReads" >> {output.frip}
        echo "FRiP: $frip" >> {output.frip}
        
        rm -f {output.frip}.saf {output.frip}.counts {output.frip}.counts.summary
        """

rule generate_pooled_count_matrix:
    """
    Generate count matrix using pooled consensus peaks for differential analysis.
    
    This is the PRIMARY count matrix used for downstream DEG analysis.
    Uses group-level pooled consensus peaks to ensure sufficient coverage
    for statistical analysis.
    """
    input:
        consensus = "04.consensus/pooled/all_groups_consensus_peaks.bed",
        bams = expand("02.mapping/shifted/{sample}.shifted.sorted.bam", sample=samples.keys())
    output:
        counts_matrix = "04.consensus/pooled/consensus_counts_matrix.txt",
        description = "04.consensus/pooled/matrix_description.txt",
        summary = "04.consensus/pooled/consensus_counts_matrix.txt.summary",
        saf = temp("04.consensus/pooled/consensus_peaks.saf")
    conda:
        workflow.source_path("../envs/subread.yaml"),
    resources:
        **rule_resource(config, 'high_resource', skip_queue_on_local=True, logger=logger),
    log:
        "logs/04.consensus/pooled/featureCounts.log"
    benchmark:
        "benchmarks/04.consensus/pooled/featureCounts.txt"
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
            
        echo "3. Cleaning up matrix header..." >> {log}
        sed -i 's|02.mapping/shifted/||g; s|.shifted.sorted.bam||g' {output.counts_matrix}

        echo "4. Generating description file..." >> {log}
        echo "File Name: $(basename {output.counts_matrix})" > {output.description}
        echo "Generated Date: $(date +'%Y-%m-%d %H:%M:%S')" >> {output.description}
        echo "--------------------------------------------------" >> {output.description}
        echo "Total Samples: $(echo "{input.bams}" | wc -w)" >> {output.description}
        echo "Total Consensus Peaks: $(wc -l < {input.consensus})" >> {output.description}
        """

rule generate_pooled_count_matrix_ann:
    """
    Add gene annotation to the pooled consensus count matrix.
    """
    input:
        annotation = "04.consensus/pooled/consensus_peaks_annotation.txt",
        counts_matrix = "04.consensus/pooled/consensus_counts_matrix.txt",
    output:
        counts_matrix_ann = "04.consensus/pooled/consensus_counts_matrix_ann.txt",
    conda:
        workflow.source_path("../envs/bedtools.yaml"),
    resources:
        **rule_resource(config, 'low_resource', skip_queue_on_local=True, logger=logger),
    log:
        "logs/04.consensus/pooled/consensus_counts_matrix_ann.log",
    benchmark:
        "benchmarks/04.consensus/pooled/consensus_counts_matrix_ann.txt",
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