#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
import os

def get_java_opts(wildcards, input, resources):
    mem_gb = max(int(resources.mem_mb / 1024) - 4, 2)
    if mem_gb > 50:
        mem_gb = 20
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

# --------------- Mapping Rules --------------- #
if config.get("mapping_tools","chromap") == "chromap":
    include: "./subrules/mapping/chromap.smk"
    logger.info("ATAC mapping powered by Chromap") 
else:
    include: "./subrules/mapping/bowtie2.smk"
    logger.info("ATAC mapping powered by bowtie2") 
# --------------- Mapping Rules --------------- #
rule sorted_bam:
    """
    Convert SAM to sorted BAM format and create index for efficient access.

    This rule converts the SAM format output from the aligner to the compressed BAM format,
    which is the standard binary format for storing aligned sequencing data. The BAM file
    is coordinate-sorted to enable efficient random access by genomic position, and an
    index is generated to support rapid data retrieval for downstream analysis tools.

    Key processing steps:
    - Convert SAM to BAM format using samtools (compression reduces file size)
    - Sort BAM by chromosomal coordinates for indexed access
    - Generate BAM index (.bai) for fast random access to genomic regions

    The sorted BAM and its index serve as the foundation for downstream analyses
    including peak calling, coverage visualization, and quality assessment.
    """
    input:
        bam = '02.mapping/Aligner/{sample}/{sample}.bam',
    output:
        sort_bam = '02.mapping/Aligner/{sample}/{sample}.sorted.bam',
        sort_bam_bai = '02.mapping/Aligner/{sample}/{sample}.sorted.bam.bai',
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

rule bam2cram:
    """
    Convert BAM files to CRAM format for storage efficiency.

    This rule compresses the sorted BAM files into CRAM format, which typically
    achieves 40-60% smaller file sizes compared to BAM while maintaining full
    compatibility with most bioinformatics tools. The CRAM format uses reference-
    based compression, making it ideal for large-scale ATAC-seq projects where
    storage costs are a concern.
    """
    input:
        sort_bam = '02.mapping/Aligner/{sample}/{sample}.sorted.bam',
        sort_bam_bai = '02.mapping/Aligner/{sample}/{sample}.sorted.bam.bai',
    output:
        cram = '02.mapping/cram/{sample}.cram',
        cram_index = '02.mapping/cram/{sample}.cram.crai',
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/bwa2.yaml"),
    message:
        "Converting BAM to CRAM format : {input.sort_bam}",
    log:
        "logs/02.mapping/bam_dup_stats_{sample}.log",
    benchmark:
        "benchmarks/{sample}_Dup_bam_stats_benchmark.txt",
    threads:
        config['parameter']['threads']['bam2cram'],
    params:
        reference = config['Bowtie2_index'][config['Genome_Version']]['genome_fa'],
    shell:
        """
        samtools view -@ {threads} -C -T {params.reference} -o {output.cram} {input.sort_bam}
        samtools index  -@ {threads} {output.cram}
        """

rule samtools_flagst:
    """
    Generate flag statistics for BAM files using samtools flagstat.

    This rule produces comprehensive statistics about the alignment flags in the BAM file,
    providing a breakdown of how reads are categorized based on their SAM flag values.
    Flagstat reports essential quality metrics that help assess the success of the
    alignment step and identify potential issues with the sequencing data.

    Key statistics reported include:
    - Total number of reads in the BAM file
    - Number of reads mapped to the reference genome
    - Properly paired reads (reads mapped in correct orientation and distance)
    - Reads mapped as singletons (only one end of pair mapped)
    - Reads mapped to different chromosomes or with unexpected orientations
    - Duplicate reads and supplementary alignments

    These metrics are crucial for quality control, enabling the identification of
    alignment problems such as high unmapped rates, poor pairing efficiency, or
    contamination with adapter dimers. The tab-separated format facilitates downstream
    parsing and integration into quality control reports.
    """
    input:
        bam = '02.mapping/Aligner/{sample}/{sample}.sorted.bam',
        bai = '02.mapping/Aligner/{sample}/{sample}.sorted.bam.bai',
    output:
        samtools_flagstat = '02.mapping/samtools_flagstat/{sample}_bam_flagstat.tsv',
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/samtools.yaml"),
    message:
        "Running samtools flagstat for BAM : {input.bam}",
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
    """
    Generate comprehensive statistics for BAM files using samtools stats.

    This rule produces detailed statistical summaries of the BAM file, including
    alignment quality metrics, coverage statistics, insert size distributions,
    and base composition analyses. The samtools stats tool provides comprehensive
    quality control data that helps assess the overall success of the sequencing
    experiment and identify potential technical issues.

    Key statistics and metrics reported include:
    - Summary statistics: total reads, mapped reads, duplicate rates
    - Quality metrics: average base quality scores, quality distributions
    - Insert size statistics: mean, median, and standard deviation of fragment lengths
    - Coverage statistics: mean coverage, coverage histograms
    - Base composition: nucleotide frequencies and GC content
    - Read length distributions and mapping quality scores
    - Chromosome-specific statistics and coverage summaries

    The reference genome is used to calculate coverage statistics and evaluate
    the quality of alignments against the expected genomic sequence. This detailed
    statistical report is essential for quality control, enabling the assessment
    of library quality, sequencing depth, and the identification of any technical
    artifacts that might affect downstream analyses. The tab-separated output format
    facilitates easy parsing and integration into quality control summaries and
    multi-sample comparison reports.
    """
    input:
        bam = '02.mapping/Aligner/{sample}/{sample}.sorted.bam',
        bai = '02.mapping/Aligner/{sample}/{sample}.sorted.bam.bai',
    output:
        samtools_stats = '02.mapping/samtools_stats/{sample}_bam_stats.tsv',
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/samtools.yaml"),
    message:
        "Running samtools stats for BAM : {input.bam}",
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
    Add read groups to BAM file using samtools addreplacerg.

    This rule assigns read group (RG) information to each alignment in the BAM file,
    which is essential for downstream GATK tools and other analysis pipelines that
    require read group metadata. Read groups allow tracking of sequencing data
    from different libraries, lanes, or flow cells, enabling proper handling of
    technical artifacts and batch effects.

    Key read group fields assigned:
    - ID: Read group identifier (set to sample name for unique identification)
    - LB: Library identifier (set to "lib1" for tracking library preparation)
    - PL: Platform/technology used (configured via config, typically "illumina")
    - PU: Platform unit/flow cell barcode and lane (set to "unit1")
    - SM: Sample name (derived from wildcards.sample for sample tracking)

    The read group information is crucial for:
    - Distinguishing reads from different sequencing runs or lanes
    - Identifying and correcting batch effects in downstream analyses
    - Enabling duplicate marking algorithms to work correctly across read groups
    - Facilitating variant calling and other analyses that require read group awareness

    The output BAM file maintains all original alignment information while adding
    the structured read group metadata required by GATK and other analysis tools.
    The coordinate sorting and indexing ensure compatibility with downstream
    processing steps including duplicate marking and variant analysis.
    """
    input:
        bam = '02.mapping/Aligner/{sample}/{sample}.sorted.bam',
        sort_bai = '02.mapping/Aligner/{sample}/{sample}.sorted.bam.bai',
    output:
        bam = '02.mapping/gatk/{sample}/{sample}.rg.bam',
        bai = '02.mapping/gatk/{sample}/{sample}.rg.bam.bai',
    conda:
        workflow.source_path("../envs/samtools.yaml")
    log:
        "logs/02.mapping/gatk/AddRG/{sample}.log"
    benchmark:
        "benchmarks/02.mapping/gatk/AddRG/{sample}.txt"
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    threads: 
        4
    params:
        java_opts = get_java_opts,
        PL = config["parameter"]["AddOrReplaceReadGroups"]["PL"],
    shell:
        """
        # 使用 samtools addreplacerg 为 BAM 文件添加 read group 信息
        # -@ {threads}: 指定使用的线程数
        # -r: 指定 read group 字符串，格式为 @RG\tID=xxx\tLB=xxx\tPL=xxx\tPU=xxx\tSM=xxx
        #     ID: Read group 唯一标识符，使用样本名
        #     LB: Library 标识符，标识文库
        #     PL: Platform 平台 (如 illumina)
        #     PU: Platform Unit 平台单元 (flow cell + lane)
        #     SM: Sample 样本名
        # -o: 输出 BAM 文件路径
        # -O BAM: 指定输出格式为 BAM
        ( samtools addreplacerg \
                    -@ {threads} \
                    -r "@RG\tID={wildcards.sample}\tLB=lib1\tPL={params.PL}\tPU=unit1\tSM={wildcards.sample}" \
                    -o {output.bam} \
                    -O BAM {input.bam} && \
        # 为输出的 BAM 文件创建索引 (.bai)
        samtools index -@ {threads} {output.bam} ) 2> {log}
        """

# --------------- mark dup Rules --------------- #
if config.get("mapping_tools","chromap") == "chromap":
    include: "./subrules/mapping/chromap_mark_duplicates.smk"
    logger.info("skiping ATAC mark dup for Chromap") 
else:
    include: "./subrules/mapping/bowtie2_mark_duplicates.smk"
    logger.info("ATAC mark dup for bowtie2") 
# --------------- mark dup Rules --------------- #


rule filter_blacklist_and_mito:
    """
    Filter BAM files to remove blacklist regions, organellar reads, low mapping quality reads, and unwanted flags
    """
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
    Filter for properly paired reads and sort by name then by position.

    This rule processes the filtered BAM files to retain only properly paired reads
    and performs coordinate-based sorting to prepare the data for downstream analyses.
    Properly paired reads are those where both ends of the DNA fragment map to the
    reference genome in the expected orientation and within a reasonable distance
    from each other.

    Key processing steps:
    - Sort the input BAM file by read name to group paired reads together
    - Filter for properly paired reads using specialized filtering tools
    - Sort the resulting BAM file by genomic coordinates for downstream compatibility
    - Generate index files for rapid random access to genomic regions

    This filtering step is critical for ATAC-seq analysis as it ensures that only
    high-quality, properly paired fragments are used for peak calling and other
    downstream analyses. Removing improperly paired reads improves the accuracy
    of fragment length estimation, insertion site analysis, and peak detection.

    The coordinate-sorted output is compatible with downstream tools that require
    position-sorted input, such as peak callers, coverage analysis tools, and
    visualization software.
    """
    input:
        bam = '02.mapping/filtered/{sample}.clean.bam',
        bai = '02.mapping/filtered/{sample}.clean.bam.bai'
    output:
        sort_name_bam = '02.mapping/filter_pe/{sample}.sort_name.bam',
        bam = '02.mapping/filter_pe/{sample}.filter_pe.bam',
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


# --------------- tn5 shift Rules --------------- #
if config.get("mapping_tools","chromap") == "chromap":
    include: "./subrules/mapping/chromap_shift.smk"
    logger.info("skiping bam tn5 shift for Chromap") 
else:
    include: "./subrules/mapping/bowtie2_shift.smk "
    logger.info("ATAC bam tn5 shift for bowtie2") 
# ---------------tn5 shift Rules --------------- #


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