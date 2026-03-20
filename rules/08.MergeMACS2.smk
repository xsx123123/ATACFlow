#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
"""
ATACFlow Pipeline - Group-level Peak Calling Module

This module handles the merging of BAM files by experimental group and performs
group-level peak calling to identify consensus accessible chromatin regions across
biological replicates. Group-level analysis provides more robust peak calls by
increasing sequencing depth and reducing technical noise.

Key Components:
- merge_shifted_bams: Merges shifted ATAC-seq BAM files by experimental group
- macs2_merge_callpeak: Performs MACS2 peak calling on merged group BAMs
- merge_homer_annotate_peaks: Annotates group-level peaks using HOMER
- merge_create_consensus_peakset: Creates consensus peaks from group-level calls
- merge_homer_annotate_consensus_peaks: Annotates group consensus peaks
- merge_generate_count_matrix: Generates count matrix using group consensus peaks
- merge_generate_count_matrix_ann: Adds gene annotation to group count matrix

This module is particularly useful for experimental designs with multiple biological
replicates, where group-level analysis can reveal more consistent and robust
accessible chromatin regions that may be missed in individual sample analysis.
"""

import os

rule merge_shifted_bams:
    """
    Merge shifted ATAC-seq BAM files by experimental group.

    This rule combines the shifted BAM files from all samples belonging to the same
    experimental group, creating a single consolidated BAM file with increased
    sequencing depth. Merging biological replicates prior to peak calling improves
    the sensitivity and robustness of peak detection by reducing technical noise
    and increasing coverage of accessible chromatin regions.

    Key features of this merging process:
    - Uses samtools merge for efficient BAM combination
    - Maintains proper pairing information for paired-end reads
    - Creates an index for the merged BAM file for rapid access
    - Preserves all quality information and alignment tags

    The merged BAM file serves as input for group-level peak calling, which can
    identify accessible regions that may be too weak to detect in individual samples
    but become apparent when replicates are combined. This approach is particularly
    valuable for experiments with limited sequencing depth per sample.
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

rule macs2_merge_callpeak:
    """
    Identify open chromatin regions (peaks) using MACS2 on merged group BAM files.

    This rule performs peak calling on the merged BAM files from each experimental group,
    leveraging the increased sequencing depth to identify more robust and sensitive
    accessible chromatin regions. MACS2 is used with parameters optimized for ATAC-seq
    data to detect regions of open chromatin where Tn5 transposase insertion is
    significantly enriched relative to background.

    Key parameters for group-level ATAC-seq peak calling:
    - -f BAMPE: Uses paired-end information directly for improved peak modeling
    - -g {genome_size}: Specifies effective genome size for proper p-value calculation
    - -q {qvalue}: Sets q-value threshold for statistical significance (FDR control)
    - --keep-dup all: Retains all duplicates (already marked and filtered in prior steps)
    - -B --SPMR: Generates bedGraph files normalized by reads per million mapped reads

    Since the input BAM files have already been shifted to correct for Tn5 transposase bias,
    no additional shifting parameters are required in MACS2. The BAMPE mode leverages
    the paired-end nature of ATAC-seq data to build a more accurate model of fragment
    size distribution and peak shape.

    Group-level peak calling generates several important output files:
    - _peaks.narrowPeak: Standard narrowPeak format with peak coordinates and statistics
    - _peaks.xls: Detailed peak information in Excel-compatible format
    - _summits.bed: Precise summit coordinates for each peak (single base resolution)
    - _treat_pileup.bdg: Normalized read coverage signal track in bedGraph format

    These group-level peak calls represent the consensus accessible chromatin landscape
    of each experimental condition and serve as the foundation for downstream analyses
    including differential accessibility analysis between groups.
    """
    input:
        bam = "02.mapping/merged/{group}.merged.bam",
        bai = "02.mapping/merged/{group}.merged.bam.bai"
    output:
        narrow_peak = "03.peak_calling/MERGE_MACS2/{group}/{group}_peaks.narrowPeak",
        xls = "03.peak_calling/MERGE_MACS2/{group}/{group}_peaks.xls",
        summits = "03.peak_calling/MERGE_MACS2/{group}/{group}_summits.bed",
        bdg = "03.peak_calling/MERGE_MACS2/{group}/{group}_treat_pileup.bdg"
    resources:
        **rule_resource(config, 'high_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/macs2.yaml"),
    log:
        "logs/03.peak_calling/macs2_{group}.log",
    message:
        "Running MACS2 peak calling for {wildcards.group}",
    benchmark:
        "benchmarks/03.peak_calling/MERGE_macs2_{group}.txt",
    threads:
        config['parameter']['threads']['macs2'],
    params:
        gsize = config['genome_info'][config['Genome_Version']]['effectiveGenomeSize'],
        qvalue = config['parameter']['macs2']['qvalue'],
        outdir = "03.peak_calling/MERGE_MACS2/{group}"
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

rule merge_homer_annotate_peaks:
    """
    Annotate group-level peaks relative to gene features using HOMER.

    This rule annotates the consensus peaks identified from merged group BAM files
    with respect to known gene structures, providing functional context to the
    accessible chromatin regions. HOMER (Hypergeometric Optimization of Motif
    EnRichment) is used to determine the genomic context of each peak, including
    its distance to the nearest transcription start site (TSS) and the gene it
    potentially regulates.

    Key annotation information provided includes:
    - Peak location relative to gene features (promoter, exon, intron, intergenic)
    - Distance to the nearest transcription start site (TSS)
    - Nearest gene identifier and name
    - Peak statistics and enrichment metrics
    - Summary statistics of peak distribution across genomic features

    This annotation is crucial for interpreting the biological significance of
    accessible chromatin regions, as it links peaks to potential target genes and
    provides insights into the regulatory landscape of the genome under different
    experimental conditions. The annotated peaks facilitate downstream analyses
    such as gene ontology enrichment and pathway analysis of potential regulatory regions.
    """
    input:
        peak = "03.peak_calling/MERGE_MACS2/{group}/{group}_peaks.narrowPeak",
    output:
        annotation = "03.peak_calling/MERGE_HOMER/{group}_annotation.txt",
        stats = "03.peak_calling/MERGE_HOMER/{group}_stats.txt"
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/homer.yaml"),
    log:
        "logs/03.peak_calling/MERGE_homer_{group}.log",
    message:
        "Running HOMER annotation for {wildcards.group}",
    benchmark:
        "benchmarks/03.peak_calling/MERGE_homer_{group}.txt",
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

rule merge_create_consensus_peakset:
    """
    Create a consensus peakset by merging peaks from group-level MACS2 results.

    This rule generates a unified set of consensus accessible chromatin regions by
    merging peaks from all experimental groups, rather than from individual samples.
    This approach creates a more robust consensus set that captures the complete
    regulatory landscape across all experimental conditions in the study.

    Key features of this group-based consensus approach:
    - Merges peaks from group-level calls rather than individual samples
    - Uses strict sorting and overlap-based merging to create non-overlapping regions
    - Maintains the genomic coordinates of all accessible regions across conditions
    - Provides a comprehensive reference for differential accessibility analysis

    By using group-level peaks as input, this approach reduces noise from individual
    sample variation and focuses on the consistent accessible regions within each
    experimental condition. The resulting consensus peakset serves as a common
    coordinate system for quantifying accessibility across all samples and groups,
    enabling robust statistical comparisons between different experimental conditions.
    """
    input:
        peaks = expand("03.peak_calling/MERGE_MACS2/{group}/{group}_peaks.narrowPeak", group=groups.keys())
    output:
        consensus = "04.consensus/consensus_peak/{group}_consensus_peaks.bed"
    resources:
        **rule_resource(config, 'high_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/bedtools.yaml"),
    log:
        "logs/04.consensus/merge_group_peaks.log",
    message:
        "Creating consensus peakset from GROUP merged peaks",
    threads: config['parameter']['threads']['bedtools']
    shell:
        """
        cat {input.peaks} | \
        sort -k1,1 -k2,2n | \
        bedtools merge -i stdin > {output.consensus} 2> {log}
        """

rule merge_idr_peaks:
    """
    Perform Irreproducible Discovery Rate (IDR) analysis on group replicates.

    This rule identifies reproducible peaks across biological replicates within a group
    using the IDR framework. It compares all pairs of replicates, identifies peaks 
    that are consistent (IDR < 0.05), and merges them into a final consensus peakset.
    For groups with only one replicate, it provides the original peaks (converted to 
    3-column BED) as the consensus set for consistency.
    """
    input:
        peaks = lambda wildcards: expand("03.peak_calling/MACS2/{sample}/{sample}_peaks.narrowPeak",
            sample=groups[wildcards.group]
        )
    output:
        idr_bed = "03.peak_calling/MERGE_IDR/{group}/Final_Consensus_Peaks.bed",
        idr_batch_log = "03.peak_calling/MERGE_IDR/{group}/idr_batch_run.log",
        idr_peaks = lambda wildcards: expand(
            "03.peak_calling/MERGE_IDR/{sample}/{sample}_peaks.idr.narrowPeak",
            sample=groups[wildcards.group]
        )
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/idr.yaml"),
    log:
        "logs/03.peak_calling/idr_{group}.log"
    benchmark:
        "benchmarks/03.peak_calling/idr_{group}.txt"
    threads:
        config['parameter']['threads'].get('idr', 4)
    params:
        script = workflow.source_path(config["parameter"]["idr_scripts"]["path"]),
        outdir = "03.peak_calling/MERGE_IDR/{group}"
    shell:
        """
        # Ensure output directory exists
        mkdir -p {params.outdir}

        # Count number of replicates
        peak_count=$(echo "{input.peaks}" | wc -w)

        if [ "$peak_count" -ge 2 ]; then
            echo "[$(date)] Running IDR for group {wildcards.group} with $peak_count replicates..." > {log}
            python3 {params.script} \
                -i {input.peaks} \
                -o {params.outdir} \
                -t {threads} >> {log} 2>&1
        else
            echo "[$(date)] Group {wildcards.group} has only one replicate. Skipping IDR and creating fallback consensus BED..." >> {log}
            # Convert narrowPeak to 3-column BED (chr, start, end), sort and merge to ensure consistency
            awk 'BEGIN{{OFS="\\t"}} {{print $1, $2, $3}}' {input.peaks} | \
            sort -k1,1 -k2,2n | \
            bedtools merge > {output.idr_bed} 2>> {log}
        fi
        """


rule merge_homer_annotate_consensus_peaks:
    """
    Annotate the group-level consensus peaks relative to gene features using HOMER.

    This rule provides functional annotation for the consensus peakset derived from
    merging group-level peaks, linking accessible chromatin regions to potential
    target genes and regulatory elements. The annotation process uses HOMER to
    determine the genomic context of each consensus peak and its relationship to
    known gene structures.

    Key annotation information for consensus peaks includes:
    - Genomic feature classification (promoter, exon, intron, intergenic)
    - Distance to nearest transcription start site (TSS)
    - Nearest gene identifiers and functional information
    - Summary statistics of peak distribution across the genome

    Annotating the consensus peakset ensures that all regions used for differential
    accessibility analysis have consistent functional annotation, facilitating the
    interpretation of changes in chromatin accessibility between experimental conditions
    and enabling downstream functional enrichment analyses.
    """
    input:
        consensus = "04.consensus/consensus_peak/{group}_consensus_peaks.bed"
    output:
        annotation = "04.consensus/consensus_peak/{group}_consensus_peaks_annotation.txt",
        stats = "04.consensus/consensus_peak/{group}_consensus_peaks_stats.txt"
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/homer.yaml"),
    log:
        "logs/03.peak_calling/{group}_consensus_peaks.log",
    message:
        "Running HOMER annotation for {wildcards.group} consensus peaks",
    benchmark:
        "benchmarks/03.peak_calling/{group}_consensus_peaks_annotation.txt",
    threads: 
        config['parameter']['threads']['homer']
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


rule merge_homer_annotate_idr_peaks:
    """
    Annotate the group-level consensus peaks relative to gene features using HOMER.

    This rule provides functional annotation for the consensus peakset derived from
    merging group-level peaks, linking accessible chromatin regions to potential
    target genes and regulatory elements. The annotation process uses HOMER to
    determine the genomic context of each consensus peak and its relationship to
    known gene structures.

    Key annotation information for consensus peaks includes:
    - Genomic feature classification (promoter, exon, intron, intergenic)
    - Distance to nearest transcription start site (TSS)
    - Nearest gene identifiers and functional information
    - Summary statistics of peak distribution across the genome

    Annotating the consensus peakset ensures that all regions used for differential
    accessibility analysis have consistent functional annotation, facilitating the
    interpretation of changes in chromatin accessibility between experimental conditions
    and enabling downstream functional enrichment analyses.
    """
    input:
        idr_peaks = "03.peak_calling/MERGE_IDR/{sample}/{sample}_peaks.idr.narrowPeak",
    output:
        annotation = "04.consensus/idr_peak/{sample}_idr_peaks_annotation.txt",
        stats = "04.consensus/idr_peak/{sample}_idr_peaks_stats.txt",
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/homer.yaml"),
    log:
        "logs/03.idr_calling/{sample}_consensus_peaks.log",
    message:
        "Running HOMER annotation for {sample} consensus peaks",
    benchmark:
        "benchmarks/03.peak_calling/{sample}_consensus_peaks_annotation.txt",
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

#rule merge_generate_count_matrix:
#    """
#    Generate a read count matrix for the group-level consensus peakset across all samples.
#
#    This rule quantifies chromatin accessibility by counting the number of sequencing
#    reads overlapping each region in the group-level consensus peakset for every
#    sample in the experiment. The resulting count matrix provides a quantitative
#    measure of accessibility for each peak across all experimental conditions,
#    serving as the foundation for differential accessibility analysis.
#
#    Key features of the count matrix generation:
#    - Uses the consensus peakset as a common coordinate system
#    - Counts reads in each shifted BAM file for every sample
#    - Maintains sample identity for downstream statistical comparisons
#    - Generates a matrix description file for documentation
#
#    The count matrix is structured with peaks as rows and samples as columns, with
#    each cell containing the number of reads from that sample overlapping that peak.
#    This format is compatible with popular differential accessibility tools like DESeq2
#    and edgeR, enabling robust statistical analysis of chromatin accessibility changes
#    between experimental conditions.
#    """
#    input:
#        consensus = "04.consensus/group_consensus_peaks.bed",
#        bams = expand("02.mapping/shifted/{sample}.shifted.sorted.bam", sample=samples.keys())
#    output:
#        counts_matrix = "04.consensus/merge_consensus_counts_matrix.txt",
#        description = "04.consensus/merge_matrix_description.txt"
#    conda:
#        workflow.source_path("../envs/bedtools.yaml"),
#    log:
#        "logs/04.consensus/multicov.log"
#    params:
#        sample_names = list(samples.keys()),
#        path = workflow.source_path(config['parameter']['generate_atac_matrix']['path'])
#    threads: 
#        config['parameter']['threads']['bedtools']
#    shell:
#        """
#        chmod +x {params.path}
#        python3 {params.path} \
#            --bed {input.consensus} \
#            --inputs {input.bams} \
#            --samples {params.sample_names} \
#            --output {output.counts_matrix} \
#            --desc {output.description} \
#            --log {log}
#        """

rule merge_generate_count_matrix_by_featureCounts:
    """
    Generate a read count matrix using featureCounts for the group-level consensus peakset.

    This rule quantifies chromatin accessibility by counting the number of ATAC-seq
    paired-end fragments overlapping each region in the group-level consensus peakset.
    By replacing bedtools multicov with featureCounts (-p flag), we ensure that the 
    entire DNA fragment (insert size) is correctly evaluated against the peak boundaries, 
    providing a much more accurate representation of the open chromatin state.
    """
    input:
        consensus = "04.consensus/consensus_peak/{group}_consensus_peaks.bed",
        bams = expand("02.mapping/shifted/{sample}.shifted.sorted.bam", sample=samples.keys())
    output:
        counts_matrix = "04.consensus/consensus_peak/merge_consensus_counts_matrix.txt",
        description = "04.consensus/consensus_peak/merge_matrix_description.txt",
        summary = "04.consensus/consensus_peak/merge_consensus_counts_matrix.txt.summary",
        saf = temp("04.consensus/consensus_peak/{group}_consensus_peaks.saf")
    conda:
        workflow.source_path("../envs/subread.yaml"), 
    resources:
        **rule_resource(config, 'high_resource',  skip_queue_on_local=True,logger = logger),
    log:
        "logs/04.consensus/merge_featureCounts.log"
    benchmark:
        "benchmarks/04.consensus/merge_featureCounts.txt"
    threads: 
        config['parameter']['threads'].get('featurecounts', 16)
    shell:
        """
        echo "1. Converting GROUP consensus BED to SAF format..." > {log}
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
            
        echo "3. Cleaning up matrix header for downstream compatibility..." >> {log}
        sed -i 's|02.mapping/shifted/||g; s|.shifted.sorted.bam||g' {output.counts_matrix}

        echo "4. Generating description file..." >> {log}
        echo "File Name: $(basename {output.counts_matrix})" > {output.description}
        echo "Generated Date: $(date +'%Y-%m-%d %H:%M:%S')" >> {output.description}
        echo "--------------------------------------------------" >> {output.description}
        echo "Total Samples Quantified: $(echo "{input.bams}" | wc -w)" >> {output.description}
        echo "Total Group Consensus Peaks: $(wc -l < {input.consensus})" >> {output.description}
        """

rule count_matrix_idr_featureCounts:
    """
    Generate a high-confidence read count matrix using featureCounts 
    based on the IDR-filtered consensus peakset.
    """
    input:
        consensus = "idr_results/Final_Consensus_Peaks.bed",
        bams = expand("02.mapping/shifted/{sample}.shifted.sorted.bam", sample=samples.keys())
    output:
        counts_matrix = "04.consensus/idr_peak/idr_consensus_counts_matrix.txt",
        description = "04.consensus/idr_peak/idr_matrix_description.txt",
        summary = "04.consensus/idr_peak/idr_consensus_counts_matrix.txt.summary",
        saf = temp("04.consensus/idr_peak/idr_consensus_peaks.saf")
    conda:
        workflow.source_path("../envs/subread.yaml"), 
    resources:
        **rule_resource(config, 'high_resource',  skip_queue_on_local=True, logger=logger),
    log:
        "logs/04.consensus/idr_featureCounts.log"
    benchmark:
        "benchmarks/04.consensus/idr_featureCounts.txt"
    threads: 
        config['parameter']['threads'].get('featurecounts', 16)
    shell:
        """
        echo "1. Converting IDR consensus BED to SAF format..." > {log}
        # 注意这里 awk 把原来的名称列改成了坐标拼接，保证 featureCounts 出来的 GeneID 是独一无二的
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
            
        echo "3. Cleaning up matrix header for downstream compatibility..." >> {log}
        sed -i 's|02.mapping/shifted/||g; s|.shifted.sorted.bam||g' {output.counts_matrix}

        echo "4. Generating description file..." >> {log}
        echo "File Name: $(basename {output.counts_matrix})" > {output.description}
        echo "Generated Date: $(date +'%Y-%m-%d %H:%M:%S')" >> {output.description}
        echo "--------------------------------------------------" >> {output.description}
        echo "Total Samples Quantified: $(echo "{input.bams}" | wc -w)" >> {output.description}
        echo "Total IDR Consensus Peaks: $(wc -l < {input.consensus})" >> {output.description}
        """

rule merge_generate_count_matrix_ann:
    """
    Add gene annotation information to the group-level count matrix.

    This rule combines the quantitative accessibility data from the count matrix with
    the functional annotation from HOMER, creating an integrated matrix that includes
    both accessibility measurements and gene context for each consensus peak. This
    annotated matrix facilitates the interpretation of differential accessibility results
    by immediately linking accessibility changes to potential target genes.

    Key integration features:
    - Merges peak coordinates with HOMER annotation data
    - Preserves all count data for differential analysis
    - Adds gene identifiers, names, and functional information
    - Maintains genomic location and peak statistics

    The resulting annotated count matrix serves as a comprehensive resource for both
    computational analysis and biological interpretation, enabling researchers to
    quickly connect changes in chromatin accessibility to potential gene regulatory
    effects and functional pathways.
    """
    input:
        annotation = "04.consensus/consensus_peak/{group}_consensus_peaks_annotation.txt",
        counts_matrix = "04.consensus/consensus_peak/merge_consensus_counts_matrix..txt",
    output:
        counts_matrix_ann = "04.consensus/consensus_peak/merge_consensus_counts_matrix_ann.txt",
    conda:
        workflow.source_path("../envs/bedtools.yaml"),
    resources:
        **rule_resource(config, 'low_resource',  skip_queue_on_local=True,logger = logger),
    log:
        "logs/04.consensus/consensus_counts_matrix_ann.log",
    benchmark:
        "benchmarks/04.consensus/consensus_counts_matrix_ann.txt",
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
# ----- end of rules ----- #