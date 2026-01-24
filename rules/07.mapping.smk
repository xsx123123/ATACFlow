#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
import os

def get_java_opts(wildcards, input, resources):
    mem_gb = max(int(resources.mem_mb / 1024) - 4, 2)
    return f"-Xmx{mem_gb}g -XX:+UseParallelGC -XX:ParallelGCThreads=4"

def get_blacklist_path(wildcards):
    """
    Get the blacklist path based on genome build
    """
    build = config.get("Genome_Version")
    blacklist_dict = config.get('Bowtie2_index',{}).get(build,{}).get("blacklists", {})

    if not blacklist_dict:
        return ""
    else:
        blacklist_dir = os.path.join(config.get("reference_path"),blacklist_dict)
    return blacklist_dir

def get_organelle_names(wildcards):
    build = config.get("Genome_Version")
    genome_cfg = config.get("genome_info", {}).get(build, {})

    org_list = []
    for key in ["chrMID", "plast"]:
        val = genome_cfg.get(key)
        if val:
            if isinstance(val, list): org_list.extend(val)
            else: org_list.append(val)

    # Return as space-separated string for easier shell processing
    return " ".join(org_list) if org_list else ""

rule Bowtie2_mapping:
    """
    Map paired-end reads to reference genome using Bowtie2
    """
    input:
        r1 = "01.qc/short_read_trim/{sample}.R1.trimed.fq.gz",
        r2 = "01.qc/short_read_trim/{sample}.R2.trimed.fq.gz",
    output:
        bam = temp('02.mapping/Bowtie2/{sample}/{sample}.bam'),
    resources:
        **rule_resource(config, 'high_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/bowtie2.yaml"),
    log:
        "logs/02.mapping/Bowtie2_{sample}.log",
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
        ( ulimit -n 65535 && bowtie2 \
                -p {threads} \
                -X 2000 \
                --very-sensitive \
                --no-mixed \
                --no-discordant \
                -x {params.index} \
                -1 {input.r1} \
                -2 {input.r2} | samtools view -bS - > {output.bam} ) &> {log}
        """

rule sam_to_sorted_bam:
    """
    Convert SAM to sorted BAM and create index
    """
    input:
        bam = '02.mapping/Bowtie2/{sample}/{sample}.bam',
    output:
        sort_bam = '02.mapping/Bowtie2/{sample}/{sample}.sorted.bam',
        sort_bam_bai = '02.mapping/Bowtie2/{sample}/{sample}.sorted.bam.bai',
    resources:
        **rule_resource(config, 'high_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/samtools.yaml"),  # Changed from bwa2.yaml to samtools.yaml
    message:
        "Converting SAM to BAM, sorting and indexing for {wildcards.sample}",
    log:
        "logs/02.mapping/sam_to_sorted_bam_{sample}.log",
    benchmark:
        "benchmarks/{sample}_sam_to_sorted_bam_benchmark.txt",
    threads:
        config['parameter']['threads']['samtools'],
    shell:
        """
        ( samtools sort -@ {threads} -o {output.sort_bam} {input.bam} &&
        samtools index -@ {threads} {output.sort_bam} ) &>{log}
        """

rule estimate_library_complexity:
    """
    Estimate library complexity using Preseq
    """
    input:
        sort_bam = '02.mapping/Bowtie2/{sample}/{sample}.sorted.bam',
        sort_bai = '02.mapping/Bowtie2/{sample}/{sample}.sorted.bam.bai',
    output:
        preseq = '02.mapping/preseq/{sample}.lc_extrap.txt',
        c_curve = '02.mapping/preseq/{sample}.c_curve.txt',
    resources:
        **rule_resource(config, 'low_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/Preseq.yaml"),
    message:
        "Running Preseq for {wildcards.sample}",
    log:
        "logs/02.mapping/preseq_{sample}.log",
    benchmark:
        "benchmarks/{sample}_preseq_benchmark.txt",
    threads:
        1
    shell:
        """
        exec 2> {log}
        set -x
        preseq lc_extrap -pe -v -output {output.preseq} -B {input.sort_bam}
        preseq c_curve -pe -v -output {output.c_curve} -B  {input.sort_bam}
        """

rule samtools_flagst:
    input:
        bam = '02.mapping/Bowtie2/{sample}/{sample}.sorted.bam',
        bai = '02.mapping/Bowtie2/{sample}/{sample}.sorted.bam.bai',
    output:
        samtools_flagstat = '02.mapping/samtools_flagstat/{sample}_bam_flagstat.tsv',
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/bwa2.yaml"),
    message:
        "Running flagst for MarkDuplicates of BAM : {input.bam}",
    log:
        "logs/02.mapping/bam_dup_lagstat_{sample}.log",
    benchmark:
        "benchmarks/{sample}_Dup_bam_lagstat_benchmark.txt",
    threads:
        config['parameter']["threads"]["samtools_flagstat"],
    shell:
        """
        samtools flagstat \
                 -@ {threads} \
                 -O tsv \
                 {input.bam} > {output.samtools_flagstat} 2>{log}
        """
s
rule samtools_stats:
    input:
        bam = '02.mapping/Bowtie2/{sample}/{sample}.sorted.bam',
        bai = '02.mapping/Bowtie2/{sample}/{sample}.sorted.bam.bai',
    output:
        samtools_stats = '02.mapping/samtools_stats/{sample}_bam_stats.tsv',
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/bwa2.yaml"),
    message:
        "Running stats for MarkDuplicates of BAM : {input.bam}",
    log:
        "logs/02.mapping/bam_dup_stats_{sample}.log",
    benchmark:
        "benchmarks/{sample}_Dup_bam_stats_benchmark.txt",
    threads:
        config['parameter']['threads']['samtools_stats'],
    params:
        reference = config['Bowtie2_index'][config['Genome_Version']]['genome_fa'],
    shell:
        """
        samtools stats \
                 -@ {threads} \
                 --reference {params.reference} \
                 {input.bam} > {output.samtools_stats}  2>{log}
        """


rule add_read_groups:
    """
    Add read groups to BAM file using GATK
    """
    input:
        bam = '02.mapping/Bowtie2/{sample}/{sample}.sorted.bam',
        sort_bai = '02.mapping/Bowtie2/{sample}/{sample}.sorted.bam.bai',
    output:
        bam = temp('02.mapping/gatk/{sample}/{sample}.rg.bam'),
    conda:
        workflow.source_path("../envs/gatk.yaml")
    log:
        "logs/02.mapping/gatk/AddRG/{sample}.log"
    benchmark:
        "benchmarks/02.mapping/gatk/AddRG/{sample}.txt"
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    threads: 1
    params:
        java_opts = get_java_opts
    shell:
        """
        gatk --java-options "{params.java_opts}" AddOrReplaceReadGroups \
             -I {input.bam} \
             -O {output.bam} \
             -SO coordinate \
             -ID 1 -LB lib1 \
             -PL illumina -PU unit1 \
             -SM {wildcards.sample} \
             --CREATE_INDEX true 2> {log}
        """

rule mark_duplicates:
    """
    Mark PCR duplicates using GATK MarkDuplicates
    """
    input:
        bam = '02.mapping/gatk/{sample}/{sample}.rg.bam',
    output:
        bam = temp('02.mapping/gatk/{sample}/{sample}.rg.dedup.bam'),
        metrics = '02.mapping/gatk/{sample}/{sample}.rg.dedup.metrics.txt',
    resources:
        **rule_resource(config, 'high_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/gatk.yaml")
    log:
        "logs/02.mapping/gatk/mark_dup_{sample}.log"
    benchmark:
        "benchmarks/02.mapping/gatk/mark_dup_{sample}.txt"
    threads: 2
    params:
        java_opts = get_java_opts
    shell:
        """
        gatk --java-options "{params.java_opts}" MarkDuplicates \
             -I {input.bam} \
             -O {output.bam} \
             -M {output.metrics} \
             --CREATE_INDEX true \
             --MAX_RECORDS_IN_RAM 5000000 \
             --SORTING_COLLECTION_SIZE_RATIO 0.5 2> {log}
        """

rule filter_blacklist_and_mito:
    input:
        bam = '02.mapping/gatk/{sample}/{sample}.rg.dedup.bam',
    output:
        bam = '02.mapping/filtered/{sample}.clean.bam',
        bai = '02.mapping/filtered/{sample}.clean.bam.bai'
    log:
        "logs/02.mapping/filter_blacklist_mito_{sample}.log"
    threads: 4
    conda:
        workflow.source_path("../envs/samtools.yaml")
    params:
        mapq = 30,
        flag_filter = 1804,
        blacklist = lambda wildcards: get_blacklist_path(wildcards),
        organelle_filter = lambda wildcards: " && ".join([f'rname != \\"{n}\\"' for n in get_organelle_names(wildcards).split()]) if get_organelle_names(wildcards) else "1"
    shell:
        """
        # 设置 blacklist 过滤命令
        if [ -n "{params.blacklist}" ]; then
            FILTER_CMD="bedtools intersect -v -a stdin -b {params.blacklist}"
            echo "Using blacklist: {params.blacklist}" > {log}
        else
            FILTER_CMD="cat"
            echo "No blacklist used." > {log}
        fi

        # 直接在命令中使用 params.organelle_filter
        echo "Filter expression: {params.organelle_filter}" >> {log}

        (samtools view -h -b -F {params.flag_filter} -q {params.mapq} -e "{params.organelle_filter}" {input.bam} | \
         eval "$FILTER_CMD" > {output.bam}) 2>> {log}

        samtools index {output.bam} >> {log} 2>&1
        """

# rule refine_bam_strict:
#    """
#    Apply strict filtering for ATAC-seq: insert size, mismatch count, soft clips
#    """
#    input:
#        bam = '02.mapping/filtered/{sample}.clean.bam',
#        bai = '02.mapping/filtered/{sample}.clean.bam.bai'
#   output:
#        bam = '02.mapping/refined/{sample}.strict.bam',
#        bai = '02.mapping/refined/{sample}.strict.bam.bai'
#    log:
#        "logs/02.mapping/refine_strict_{sample}.log"
#    benchmark:
#        "benchmarks/02.mapping/refine_strict_{sample}.txt"
#    conda:
#        workflow.source_path("../envs/samtools.yaml")
#    threads: 4
#    params:
#        max_insert = 2000,   # Maximum insert size
#        max_mismatch = 4     # Maximum allowed mismatches
#    shell:
#        """
#        echo "Starting Strict Filtering..." > {log}
#
#        samtools view -h -b -F 4 -@ {threads} \
#        -e 'abs(tlen) <= {params.max_insert} && [NM] <= {params.max_mismatch} && ! (cigar =~ "S")' \
#        {input.bam} \
#        > {output.bam} 2>> {log}
#
#        echo "------------------------------------------------" >> {log}
#        echo "Reads before refine:" >> {log}
#        samtools view -c {input.bam} >> {log}
#
#        echo "Reads after refine:" >> {log}
#        samtools view -c {output.bam} >> {log}
#
#        # Index the output
#        samtools index -@ {threads} {output.bam} >> {log} 2>&1
#        """

rule filter_proper_pairs:
    """
    Filter for properly paired reads and sort by name then by position
    """
    input:
        bam = '02.mapping/filtered/{sample}.clean.bam',
        bai = '02.mapping/filtered/{sample}.clean.bam.bai'
    output:
        sort_name_bam = temp('02.mapping/filter_pe/{sample}.sort_name.bam'),
        bam = temp('02.mapping/filter_pe/{sample}.filter_pe.bam'),
        sort_bam = '02.mapping/filter_pe/{sample}.filter_pe.sorted.bam',
        sort_bam_bai = '02.mapping/filter_pe/{sample}.filter_pe.sorted.bam.bai'
    log:
        "logs/02.mapping/filter_proper_pairs_{sample}.log"
    benchmark:
        "benchmarks/02.mapping/filter_proper_pairs_{sample}.txt"
    conda:
        workflow.source_path("../envs/samtools.yaml")
    threads: 10
    params:
        filter_pe = workflow.source_path(config['parameter']['filter_pe']['path']),
    shell:
        """
         ( chmod +x {params.filter_pe} && \
        samtools sort -n -@ {threads} -o {output.sort_name_bam} {input.bam} && \
        {params.filter_pe} -t {threads} -i {output.sort_name_bam} -o {output.bam} && \
        samtools sort -@ {threads} {output.bam} -o {output.sort_bam} && \
        samtools index -@ {threads} {output.sort_bam} ) > {log} 2>&1
        """

rule atac_seq_shift:
    """
    Perform ATAC-seq specific shift to account for Tn5 transposase binding bias
    """
    input:
        bam = '02.mapping/filter_pe/{sample}.filter_pe.sorted.bam',
    output:
        shifted_bam = '02.mapping/shifted/{sample}.shifted.bam',
        shifted_sort_bam = '02.mapping/shifted/{sample}.shifted.sorted.bam',
        shifted_sort_bam_bai = '02.mapping/shifted/{sample}.shifted.sorted.bam.bai'
    log:
        "logs/02.mapping/atac_seq_shift_{sample}.log"
    conda:
        workflow.source_path("../envs/deeptools.yaml")
    threads:
        20
    shell:
        """
        alignmentSieve -b {input.bam} -o {output.shifted_bam} --ATACshift -p {threads} 2>> {log} && \
        samtools sort -@ {threads} {output.shifted_bam} -o {output.shifted_sort_bam} && \
        samtools index {output.shifted_sort_bam}
        """

rule generate_bigwig_coverage:
    """
    Generate normalized BigWig coverage files from BAM
    """
    input:
        shifted_sort_bam = '02.mapping/shifted/{sample}.shifted.sorted.bam',
        shifted_sort_bam_bai = '02.mapping/shifted/{sample}.shifted.sorted.bam.bai'
    output:
        bw = f"02.mapping/bamCoverage/{{sample}}_{config['parameter']['bamCoverage']['normalizeUsing']}.bw"
    resources:
        **rule_resource(config, 'high_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/deeptools.yaml"),
    message:
        "Running bamCoverage (bigwig generation) for {input.shifted_sort_bam}"
    log:
        "logs/02.mapping/bamCoverage_{sample}.log",
    benchmark:
        "benchmarks/{sample}_bamCoverage_benchmark.txt",
    threads:
        config['parameter']['threads']['bamCoverage'],
    params:
        binSize = config['parameter']['bamCoverage']['binSize'],
        smoothLength = config['parameter']['bamCoverage']['smoothLength'],
        normalizeUsing = config['parameter']['bamCoverage']['normalizeUsing'],
        effectiveGenomeSize = config['genome_info'][config['Genome_Version']]['effectiveGenomeSize'],
    shell:
        """
        bamCoverage --bam {input.shifted_sort_bam} -o {output.bw} \
            --binSize {params.binSize} \
            --normalizeUsing {params.normalizeUsing} \
            --effectiveGenomeSize {params.effectiveGenomeSize} \
            --ignoreForNormalization chrX chrY chrM \
            -p {threads} &> {log}
        """

rule tss_enrichment_analysis:
    """
    Compute TSS enrichment profile and generate plot
    """
    input:
        bw = f"02.mapping/bamCoverage/{{sample}}_{config['parameter']['bamCoverage']['normalizeUsing']}.bw"
    output:
        matrix = "02.mapping/computeMatrix/{sample}_TSS_matrix.gz",
        plot = "02.mapping/plots/{sample}_TSS_enrichment.png"
    params:
        gene_bed = config['Bowtie2_index'][config['Genome_Version']]['bed'],
        referencePoint = config['parameter']['draw_tss_plot']['referencePoint'],
        range_up_down = config['parameter']['draw_tss_plot']['range'],
    log:
        "logs/02.mapping/tss_enrichment_{sample}.log"
    conda:
        workflow.source_path("../envs/deeptools.yaml"),
    threads: 20
    shell:
        """
        # Compute matrix around TSS
        computeMatrix reference-point --referencePoint {params.referencePoint} \
                      -b {params.range_up_down} -a {params.range_up_down} \
                      -R {params.gene_bed} \
                      -S {input.bw} \
                      --skipZeros \
                      -o {output.matrix} \
                      -p {threads} &> {log}

        # Generate TSS enrichment plot
        plotProfile -m {output.matrix} \
                    -out {output.plot} \
                    --plotTitle "ATAC-seq TSS Enrichment for {wildcards.sample}" \
                    --dpi 1000 \
                    --perGroup >> {log} 2>&1
        """
# ----- end of rules ----- #