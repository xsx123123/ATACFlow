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
    log:
        "logs/03.peak_calling/macs2_{sample}.log"
    benchmark:
        "benchmarks/03.peak_calling/macs2_{sample}.txt"
    conda:
        workflow.source_path("../envs/macs2.yaml")
    threads: 1
    params:
        gsize = config.get('parameter', {}).get('macs2', {}).get('gsize', 'hs'),
        qvalue = config.get('parameter', {}).get('macs2', {}).get('qvalue', 0.05),
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
    conda:
        workflow.source_path("../envs/homer.yaml")
    log:
        "logs/03.peak_calling/homer_{sample}.log"
    threads: 1
    params:
        gtf = config['genome']['gtf'],
        genome_fasta = config['genome']['fasta']
    shell:
        """        
        annotatePeaks.pl {input.peak} {params.genome_fasta} \
            -gtf {params.gtf} \
            -annStats {output.stats} \
            > {output.annotation} 2> {log}
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
    conda:
        workflow.source_path("../envs/bedtools.yaml")
    log:
        "logs/04.consensus/merge_peaks.log"
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
        bams = expand("02.mapping/shifted/{sample}.shifted.sorted.bam", sample=SAMPLES),
        bais = expand("02.mapping/shifted/{sample}.shifted.sorted.bam.bai", sample=SAMPLES)
    output:
        counts = "04.consensus/raw_counts.txt",
        counts_with_header = "04.consensus/consensus_counts_matrix.txt"
    conda:
        workflow.source_path("../envs/bedtools.yaml")
    log:
        "logs/04.consensus/multicov.log"
    params:
        sample_names = samples.keys()
    shell:
        """
        # 1. 使用 bedtools multicov 计算
        # -bams 后面跟所有的 bam 文件
        # -bed 输入 consensus bed
        echo "Running multicov..." > {log}
        bedtools multicov -bams {input.bams} -bed {input.consensus} > {output.counts} 2>> {log}
        
        # 2. 添加表头 (Header) 让文件更易读
        # multicov 的输出前三列是 chrom, start, end，后面跟着每个 BAM 的计数
        
        HEADER="chrom\tstart\tend\t"$(echo "{params.sample_names}" | sed 's/ /\t/g')
        echo -e "$HEADER" | cat - {output.counts} > {output.counts_with_header}
        """