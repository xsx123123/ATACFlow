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

def get_organelle_filter_expr(wildcards):
    """
    Generate filter expression for organelle chromosomes
    """
    build = config.get("Genome_Version")
    org_names = config.get("genome_info", {}).get(build, []).get("chrMID", [])

    if not org_names:
        return ""

    expr_list = [f'rname != "{name}"' for name in org_names]
    full_expr = " && ".join(expr_list)

    return f"-e '{full_expr}'"

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
        ( preseq lc_extrap -pe -output {output.preseq} {input.sort_bam} && \
        preseq c_curve -pe -output {output.c_curve} {input.sort_bam} ) 2> {log}
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
        bai = temp('02.mapping/gatk/{sample}/{sample}.rg.bam.bai')
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
        bai = '02.mapping/gatk/{sample}/{sample}.rg.bam.bai',
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
    """
    Filter out blacklist regions and mitochondrial/organelle reads
    """
    input:
        bam = '02.mapping/gatk/{sample}/{sample}.rg.dedup.bam',  # Fixed path
    output:
        bam = '02.mapping/filtered/{sample}.clean.bam',
        bai = '02.mapping/filtered/{sample}.clean.bam.bai'
    log:
        "logs/02.mapping/filter_blacklist_mito_{sample}.log"
    conda:
        workflow.source_path("../envs/samtools.yaml")
    threads: 4
    params:
        mapq = 30,
        flag_filter = 1804,  # Remove unmapped, secondary, QC failed, duplicates
        blacklist = lambda wildcards: get_blacklist_path(wildcards),
        mito_filter_arg = lambda wildcards: get_organelle_filter_expr(wildcards)
    shell:
        """
        # 1. 处理 Blacklist 逻辑
        if [ -n "{params.blacklist}" ]; then
            FILTER_CMD="bedtools intersect -v -a stdin -b {params.blacklist}"
            echo "Using blacklist: {params.blacklist}" > {log}
        else
            FILTER_CMD="cat"
            echo "No blacklist used." > {log}
        fi

        # 2. 打印细胞器过滤信息到日志
        if [ -n "{params.mito_filter_arg}" ]; then
            echo "Filtering organelles with: {params.mito_filter_arg}" >> {log}
        else
            echo "No organelle filtering applied (names not found in config)." >> {log}
        fi

        # 3. 执行管道 - first apply filters, then remove blacklist regions
        (samtools view -h -b -F {params.flag_filter} -q {params.mapq} {params.mito_filter_arg} {input.bam} | \\
         $FILTER_CMD > {output.bam}) 2>> {log}

        samtools index {output.bam} >> {log} 2>&1
        """

rule refine_bam_strict:
    """
    Apply strict filtering for ATAC-seq: insert size, mismatch count, soft clips
    """
    input:
        bam = '02.mapping/filtered/{sample}.clean.bam',
        bai = '02.mapping/filtered/{sample}.clean.bam.bai'
    output:
        bam = '02.mapping/refined/{sample}.strict.bam',
        bai = '02.mapping/refined/{sample}.strict.bam.bai'
    log:
        "logs/02.mapping/refine_strict_{sample}.log"
    benchmark:
        "benchmarks/02.mapping/refine_strict_{sample}.txt"
    conda:
        workflow.source_path("../envs/samtools.yaml")
    threads: 4
    params:
        max_insert = 2000,   # Maximum insert size
        max_mismatch = 4     # Maximum allowed mismatches
    shell:
        """
        echo "Starting Strict Filtering..." > {log}

        # 1. abs(tlen) <= {params.max_insert}: Filter for insert size <= 2000bp
        # 2. [NM] <= {params.max_mismatch}: Filter for mismatch count <= 4
        # 3. ! (cigar =~ "S"): Exclude reads with soft clipping

        samtools view -h -b -e 'abs(tlen) <= {params.max_insert} && [NM] <= {params.max_mismatch}' \
        {input.bam} \
        > {output.bam} 2>> {log}

        echo "------------------------------------------------" >> {log}
        echo "Reads before refine:" >> {log}
        samtools view -c {input.bam} >> {log}
        echo "Reads after refine:" >> {log}
        samtools view -c {output.bam} >> {log}

        # Index the output
        samtools index {output.bam} >> {log} 2>&1
        """

rule filter_proper_pairs:
    """
    Filter for properly paired reads and sort by name then by position
    """
    input:
        bam = '02.mapping/refined/{sample}.strict.bam',
        bai = '02.mapping/refined/{sample}.strict.bam.bai'
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
        chmod +x {params.filter_pe} && \
        # Sort by name for PE filtering
        samtools sort -n -@ {threads} -o {output.sort_name_bam} {input.bam} && \
        # Filter properly paired reads
        {params.filter_pe} -t {threads} -i {output.sort_name_bam} -o {output.bam} && \
        # Sort by coordinate for downstream analysis
        samtools sort -@ {threads} {output.bam} -o {output.sort_bam} && \
        # Index the final BAM
        samtools index -@ {threads} {output.sort_bam} > {log} 2>&1
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
        workflow.source_path("../envs/deeptools.yaml")  # alignmentSieve is often in deeptools env
    threads:
        20
    shell:
        """
        # Apply ATAC-seq shift correction
        alignmentSieve -b {input.bam} -o {output.shifted_bam} --ATACshift -p {threads} 2>> {log} && \
        # Sort and index the shifted BAM
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
        "Running bamCoverage (bigwig generation) for {input.bam}"
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