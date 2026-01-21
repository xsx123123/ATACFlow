#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
import os

def get_blacklist_path(wildcards):
    # 1. 获取当前分析的基因组版本
    # 假设你的 config 结构是 config['genome_build']
    # 如果你是每个 sample 对应不同物种，这里逻辑要根据 sample_table 来改
    build = config.get("genome_build") 
    
    # 2. 获取对应的 blacklist 字典
    blacklist_dict = config.get("blacklists", {})
    
    # 3. 返回路径，如果没有找到则返回 None 或 空字符串
    return blacklist_dict.get(build, "")

def get_organelle_filter_expr(wildcards):
    build = config.get("genome_build")
    org_names = config.get("organelle_names", {}).get(build, [])
    
    if not org_names:
        return "" 

    expr_list = [f'rname != "{name}"' for name in org_names]
    full_expr = " && ".join(expr_list)
    
    return f"-e '{full_expr}'"

rule Bowtie_mapping:
    input:
        r1 = "01.qc/short_read_trim/{sample}.R1.trimed.fq.gz",
        r2 = "01.qc/short_read_trim/{sample}.R2.trimed.fq.gz",
    output:
        Aligned_sam =  '02.mapping/Bowtie2/{sample}/{sample}.sam',
    resources:
        **rule_resource(config, 'high_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/Bowtie.yml"),
    log:
        "logs/02.mapping/Bowtie_{sample}.log",
    message:
        "Running Bowtie mapping on {wildcards.sample} R1 and R2",
    benchmark:
        "benchmarks/{sample}_Bowtie_benchmark.txt",
    params:
        index = config['Bowtie2_index'][config['Genome_Version']]['index'],
    threads:
        config['parameter']['threads']['bowtie2'],
    shell:
        """
        ulimit -n 65535 && bowtie2 \
                -p {threads} \
                -X 2000 \
                --very-sensitive \
                --no-mixed \
                --no-discordant \
                -x {params.index} \
                -1 {input.r1} \
                -2 {input.r2} \
                -S {output.Aligned_sam} &> {log}
        """

rule sort_index:
    input:
        Aligned_sam =  '02.mapping/Bowtie2/{sample}/{sample}.sam',
    output:
        bam = '02.mapping/Bowtie/sort_index/{sample}.bam',
        sort_bam = '02.mapping/Bowtie/sort_index/{sample}.sort.bam',
        sort_bam_bai = '02.mapping/Bowtie/sort_index/{sample}.sort.bam.bai',
    resources:
        **rule_resource(config, 'high_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/bwa2.yaml"),
    message:
        "Running samtools sort & index for {wildcards.sample}",
    log:
        "logs/02.mapping/bwa_sort_index_{sample}.log",
    benchmark:
            "benchmarks/{sample}_bam_sort_index_benchmark.txt",
    threads:
        config['parameter']['threads']['samtools'],
    shell:
        """
        ( samtools view -@ {threads} -b -S {input.Aligned_sam}  -o  {output.bam} &&
        samtools sort -@ {threads} -o {output.sort_bam} {output.bam} &&
        samtools index -@ {threads} {output.sort_bam})  &>{log}
        """


rule AddOrReplaceReadGroups:
    input:
        bam = '02.mapping/Bowtie/sort_index/{sample}.sort.bam',
        bai = '02.mapping/Bowtie/sort_index/{sample}.sort.bam.bai',
    output:
        bam = temp('02.mapping/gatk/{sample}/{sample}.rg.bam'),
        bai = temp('02.mapping/gatk/{sample}/{sample}.rg.bai')
    conda:
        workflow.source_path("../envs/gatk.yaml")
    log:
        "logs/02.mapping/gatk/AddRG/{sample}.log"
    benchmark:
        "benchmarks/02.mapping/gatk/AddRG/{sample}.txt"
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
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

rule MarkDuplicates:
    input:
        bam = '02.mapping/gatk/{sample}/{sample}.rg.bam',
    output:
        bam = temp('02.mapping/gatk/{sample}/{sample}.rg.dedup.bam'),
        metrics = '02.mapping/gatk_MarkDuplicates/{sample}.rg.dedup.metrics.txt',
    resources:
        **rule_resource(config, 'high_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/gatk.yaml")
    log:
        "logs/02.mapping/gatk/MarkDup/{sample}.log"
    benchmark:
        "benchmarks/02.mapping/gatk/MarkDup/{sample}.txt"
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

rule FilterBam:
    input:
        bam = '04.variant/gatk/{sample}/{sample}.rg.dedup.bam',
        bai = '04.variant/gatk/{sample}/{sample}.rg.dedup.bam.bai'
        # blacklist 在 params 里处理
    output:
        bam = '05.filter/{sample}.clean.bam',
        bai = '05.filter/{sample}.clean.bam.bai'
    log:
        "logs/05.filter/{sample}.filter.log"
    conda:
        workflow.source_path("../envs/samtools.yaml")
    threads: 4
    params:
        mapq = 30,
        flag_filter = 1804,
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

        # 3. 执行管道
        
        (samtools view -h -b -F {params.flag_filter} -q {params.mapq} \
            {params.mito_filter_arg} \
            {input.bam} \
        | $FILTER_CMD \
        > {output.bam}) 2>> {log}

        samtools index {output.bam} >> {log} 2>&1
        """

rule RefineBam:
    input:
        # 承接上一步 FilterBam 的输出
        bam = '05.filter/{sample}.clean.bam',
        bai = '05.filter/{sample}.clean.bam.bai'
    output:
        # 生成最终用于 Peak Calling 的严格过滤文件
        bam = '06.refine/{sample}.strict.bam',
        bai = '06.refine/{sample}.strict.bam.bai'
    log:
        "logs/06.refine/{sample}.refine.log"
    benchmark:
        "benchmarks/06.refine/{sample}.txt"
    conda:
        # 依然使用包含 samtools 的环境
        workflow.source_path("../envs/samtools.yaml")
    threads: 4
    params:
        max_insert = 2000,   # 最大插入片段长度
        max_mismatch = 4     # 最大允许错配数
    shell:
        """
        echo "Starting Strict Filtering..." > {log}
        
        # 1. abs(tlen) <= {params.max_insert}: 筛选 Insert Size <= 2000 的片段
        # 2. [NM] <= {params.max_mismatch}: 筛选错配数 (NM tag) <= 4 的 reads
        # 3. ! (cigar =~ "S"): 【高风险】剔除 CIGAR 字符串中包含 "S" (Soft-clip) 的 reads
        
        samtools view -h -b -e 'abs(tlen) <= {params.max_insert} && [NM] <= {params.max_mismatch} && ! (cigar =~ "S")' \
        {input.bam} \
        > {output.bam} 2>> {log}

        echo "------------------------------------------------" >> {log}
        echo "Reads before refine:" >> {log}
        samtools view -c {input.bam} >> {log}
        echo "Reads after refine:" >> {log}
        samtools view -c {output.bam} >> {log}

        # 建索引
        samtools index {output.bam} >> {log} 2>&1
        """





rule Preseq:
    shell:
        """
        preseq lc_extrap -pe -output SampleID.lc_extrap.txt input.bam
        """



rule bam_shift:
    shell:
        """
        alignmentSieve --bam input.bam --ATACshift --outFile shifted.bam
        """

rule samtools_flagst:
    input:
        bam = '02.mapping/Bowtie/sort_index/{sample}.sort.bam',
        bai = '02.mapping/Bowtie/sort_index/{sample}.sort.bam.bai'
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

rule samtools_stats:
    input:
        bam = '02.mapping/Bowtie/sort_index/{sample}.sort.bam',
        bai = '02.mapping/Bowtie/sort_index/{sample}.sort.bam.bai'
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
        reference = config['Bowtie_index'][config['Genome_Version']]['genome_fa'],
    shell:
        """
        samtools stats \
                 -@ {threads} \
                 --reference {params.reference} \
                 {input.bam} > {output.samtools_stats}  2>{log}
        """

rule bamCoverage:
    input:
        bam = '02.mapping/Bowtie/sort_index/{sample}.sort.bam',
        bai = '02.mapping/Bowtie/sort_index/{sample}.sort.bam.bai'
    output:
        bw = f"02.mapping/bamCoverage/{{sample}}_{config['parameter']['bamCoverage']['normalizeUsing']}.bw"
    resources:
        **rule_resource(config, 'high_resource',  skip_queue_on_local=True,logger = logger),
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
        effectiveGenomeSize = 
    shell:
        """
        bamCoverage --bam {input.bam} -o {output.bw} \
            --binSize 10 \
            --normalizeUsing RPGC \
            --effectiveGenomeSize {params.effectiveGenomeSize} \
            --extendReads &> {log}           
        """

rule mapping_report:
    input:
        log_final = expand('02.mapping/Bowtie/{sample}/{sample}.Log.final.out',sample=samples.keys()),
        qualimap_report_txt = expand('02.mapping/qualimap_report/{sample}/genome_results.txt',sample=samples.keys()),
        samtools_flagstat = expand('02.mapping/samtools_flagstat/{sample}_bam_flagstat.tsv',sample=samples.keys()),
        samtools_stats = expand('02.mapping/samtools_stats/{sample}_bam_stats.tsv',sample=samples.keys()),
    output:
        report = "02.mapping/mapping_report/multiqc_mapping_report.html",
    resources:
        **rule_resource(config, 'low_resource',skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/multiqc.yaml"),
    message:
        "Running MultiQC to aggregate mapping reports",
    params:
        fastqc_reports = "02.mapping/",
        report_dir = '02.mapping/mapping_report',
        report = "multiqc_mapping_report.html",
        title = "mapping_report",
    log:
        "logs/02.mapping/multiqc_mapping_report.log",
    benchmark:
        "benchmarks/multiqc_mapping_report_benchmark.txt",
    threads:
        config['parameter']['threads']['multiqc'],
    shell:
        """
        multiqc {params.fastqc_reports} \
                --force \
                --outdir {params.report_dir} \
                -i {params.title} \
                -n {params.report} &> {log}
        """
# ----- rule ----- #