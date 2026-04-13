#!/usr/bin/snakemake
# -*- coding: utf-8 -*-

rule Bowtie2_mapping:
    """
    Align paired-end ATAC-seq reads to the reference genome using Bowtie2.

    This rule performs the critical step of mapping cleaned ATAC-seq reads to the
    reference genome using Bowtie2, a fast and memory-efficient aligner optimized
    for short reads. Proper alignment is essential for all downstream ATAC-seq
    analyses including peak calling and transcription factor footprinting.

    Bowtie2 is configured with parameters specifically optimized for ATAC-seq:
    - --very-sensitive: Enables high sensitivity alignment for better mapping rates
    - -X 2000: Sets maximum fragment length to 2000bp to accommodate nucleosome-free
                regions and mono-/di-nucleosome fragments typical in ATAC-seq
    - --no-mixed: Disables unpaired alignments, ensuring only proper pairs are considered
    - --no-discordant: Excludes discordant read pairs from alignment results

    The alignment process uses the cleaned, adapter-trimmed reads from the previous
    step and generates an unsorted BAM file as output. This BAM file serves as input
    for subsequent processing steps including coordinate sorting, duplicate marking,
    and quality filtering.

    For ATAC-seq experiments, accurate alignment is crucial because:
    - It determines the genomic location of Tn5 transposase insertions
    - Mapping quality directly impacts peak calling sensitivity and specificity
    - Proper alignment is required for nucleosome positioning analysis
    - It enables identification of accessible chromatin regions

    The resulting alignments undergo extensive quality control and filtering in
    subsequent rules to ensure only high-quality, properly paired reads are used
    for downstream analysis.
    """
    input:
        r1 = "01.qc/short_read_trim/{sample}.R1.trimed.fq.gz",
        r2 = "01.qc/short_read_trim/{sample}.R2.trimed.fq.gz",
    output:
        bam = temp('02.mapping/Aligner/{sample}/{sample}.bam'),
    resources:
        **rule_resource(config, 'high_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../../../envs/bowtie2.yaml"),
    log:
        stats = "logs/02.mapping/Bowtie2_{sample}.stats.log",
        err = "logs/02.mapping/Bowtie2_{sample}.err.log",
    message:
        "Running Bowtie2 mapping on {wildcards.sample} R1 and R2",
    benchmark:
        "benchmarks/{sample}_Bowtie2_benchmark.txt",
    params:
        index = config['Bowtie2_index'][config['Genome_Version']]['index'],
    threads:
        config['parameter']['threads']['bowtie2'],
    shell:
        """
        ulimit -n 65535 2>/dev/null || true
        bowtie2 \
            -p {threads} \
            -X 2000 \
            --very-sensitive \
            --no-mixed \
            --no-discordant \
            -x {params.index} \
            -1 {input.r1} \
            -2 {input.r2} 2> {log.stats} | \
        samtools view -b - > {output.bam} 2> {log.err}
        """
