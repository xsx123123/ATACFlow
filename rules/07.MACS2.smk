#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
import os

rule macs2_callpeak:
    """
    Identify open chromatin regions (peaks) using MACS2 on the shifted ATAC-seq BAM file.

    This rule performs peak calling, the core analysis step in ATAC-seq experiments,
    to identify genomic regions with significant chromatin accessibility. MACS2 (Model-based
    Analysis for ChIP-Seq) is used here in a configuration optimized for ATAC-seq data
    to detect regions of open chromatin where Tn5 transposase insertion is significantly
    enriched relative to background.

    Key parameters for ATAC-seq peak calling:
    - -f BAMPE: Uses paired-end information directly for improved peak modeling
    - -g {genome_size}: Specifies effective genome size for proper p-value calculation
    - -q {qvalue}: Sets q-value threshold for statistical significance (FDR control)
    - --keep-dup all: Retains all duplicates (already marked and filtered in prior steps)
    - -B --SPMR: Generates bedGraph files normalized by reads per million mapped reads

    Since the input BAM file has already been shifted to correct for Tn5 transposase bias,
    no additional shifting parameters are required in MACS2. The BAMPE mode leverages
    the paired-end nature of ATAC-seq data to build a more accurate model of fragment
    size distribution and peak shape.

    Peak calling with MACS2 generates several important output files:
    - _peaks.narrowPeak: Standard narrowPeak format with peak coordinates and statistics
    - _peaks.xls: Detailed peak information in Excel-compatible format
    - _summits.bed: Precise summit coordinates for each peak (single base resolution)
    - _treat_pileup.bdg: Normalized read coverage signal track in bedGraph format

    These peak calls represent the accessible chromatin landscape of the sample and serve
    as the foundation for all downstream analyses including annotation, differential
    accessibility analysis, transcription factor motif analysis, and pathway enrichment.
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
    Annotate peaks relative to gene features using HOMER.
    """
    input:
        consensus = "04.consensus/consensus_peaks.bed"
    output:
        annotation = "04.consensus/consensus_peaks_annotation.txt",
        stats = "04.consensus/consensus_peaks_stats.txt"
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/homer.yaml"),
    log:
        "logs/03.peak_calling/consensus_peaks.log",
    message:
        "Running HOMER annotation for consensus peaks",
    benchmark:
        "benchmarks/03.peak_calling/consensus_peaks_annotation.txt",
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

# rule generate_count_matrix:
#    input:
#        consensus = "04.consensus/consensus_peaks.bed",
#        bams = expand("02.mapping/shifted/{sample}.shifted.sorted.bam", sample=samples.keys())
#    output:
#        counts_matrix = "04.consensus/consensus_counts_matrix.txt",
#        description = "04.consensus/matrix_description.txt"
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

rule generate_count_matrix_by_featureCounts:
    """
    Generate a read count matrix using featureCounts for the consensus peakset across all samples.

    This rule replaces bedtools multicov with featureCounts to correctly count paired-end 
    ATAC-seq fragments. The -p flag ensures that the entire DNA fragment (insert size) 
    between Read 1 and Read 2 is evaluated against the peak boundaries, which is 
    the standard and most accurate quantification method for ATAC-seq.
    """
    input:
        consensus = "04.consensus/consensus_peaks.bed",
        bams = expand("02.mapping/shifted/{sample}.shifted.sorted.bam", sample=samples.keys())
    output:
        counts_matrix = "04.consensus/consensus_counts_matrix.txt",
        # 保留原有的 description 输出，防止 Snakemake 报错
        description = "04.consensus/matrix_description.txt",
        # 捕获 featureCounts 默认生成的比对统计报告
        summary = "04.consensus/consensus_counts_matrix.txt.summary",
        # 临时生成的 SAF 文件，用完自动删除
        saf = temp("04.consensus/consensus_peaks.saf") 
    conda:
        workflow.source_path("../envs/subread.yaml"), # 确保你的环境中安装了 subread 包
    log:
        "logs/04.consensus/featureCounts.log"
    benchmark:
        "benchmarks/04.consensus/featureCounts.txt"
    threads: 
        config['parameter']['threads'].get('featurecounts', 16)
    shell:
        """
        echo "1. Converting consensus BED to SAF format..." > {log}
        # 使用 awk 快速将 BED 转换为 SAF (ID, Chr, Start, End, Strand)
        awk 'BEGIN{{OFS="\\t"; print "GeneID\\tChr\\tStart\\tEnd\\tStrand"}} \
            {{print $1":"$2"-"$3, $1, $2, $3, "+"}}' {input.consensus} > {output.saf}

        echo "2. Running featureCounts for ATAC-seq PE fragments..." >> {log}
        # 核心计数命令
        # -p: 识别为双端测序并计算整个 fragment
        # -B: 仅保留两端都比对上的 fragments
        # -C: 过滤掉跨染色体的嵌合体噪音
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
        # 裁掉列名中的 BAM 路径和后缀，只保留纯净的 Sample Name，防止下游合并报错
        sed -i 's|02.mapping/shifted/||g; s|.shifted.sorted.bam||g' {output.counts_matrix}

        echo "4. Generating description file..." >> {log}
        # 用 Bash 原生命令生成轻量级的描述文件，满足原流程的 output 检查
        echo "File Name: $(basename {output.counts_matrix})" > {output.description}
        echo "Generated Date: $(date +'%Y-%m-%d %H:%M:%S')" >> {output.description}
        echo "--------------------------------------------------" >> {output.description}
        echo "Total Samples: $(echo "{input.bams}" | wc -w)" >> {output.description}
        echo "Total Consensus Peaks: $(wc -l < {input.consensus})" >> {output.description}
        """



rule generate_count_matrix_ann:
    input:
        annotation = "04.consensus/consensus_peaks_annotation.txt",
        counts_matrix = "04.consensus/consensus_counts_matrix.txt",
    output:
        counts_matrix_ann = "04.consensus/consensus_counts_matrix_ann.txt",
    conda:
        workflow.source_path("../envs/bedtools.yaml"),
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