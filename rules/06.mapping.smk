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
        sort_bam = temp('02.mapping/Aligner/{sample}/{sample}.sorted.bam'),
        sort_bam_bai = temp('02.mapping/Aligner/{sample}/{sample}.sorted.bam.bai'),
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

rule estimate_library_complexity:
    """
    Estimate library complexity and predict sequencing depth using Preseq.

    This rule analyzes the complexity of the sequencing library by estimating the
    number of unique molecules present and predicting how many additional unique
    reads could be obtained with increased sequencing depth. Library complexity
    is a critical quality metric that indicates whether the library was sufficiently
    complex or if PCR amplification bias has reduced diversity.

    Key analyses performed:
    - lc_extrap: Extrapolates the complexity curve to predict unique read yield
    - c_curve: Generates the complexity curve showing observed vs. expected unique reads

    High complexity libraries show a curve that plateaus slowly, indicating that
    additional sequencing would yield many new unique reads. Low complexity libraries
    show rapid saturation, suggesting PCR duplicates dominate and deeper sequencing
    would be inefficient. This information guides decisions about whether additional
    sequencing would be beneficial for the experiment.
    """
    input:
        sort_bam = '02.mapping/Aligner/{sample}/{sample}.sorted.bam',
        sort_bai = '02.mapping/Aligner/{sample}/{sample}.sorted.bam.bai',
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
        preseq lc_extrap -pe -v -output {output.preseq} -B {input.sort_bam} || echo "Preseq lc_extrap failed for {wildcards.sample}" > {output.preseq}
        preseq c_curve -pe -v -output {output.c_curve} -B {input.sort_bam} || echo "Preseq c_curve failed for {wildcards.sample}" > {output.c_curve}
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
    Add read groups to BAM file using GATK AddOrReplaceReadGroups.

    This rule assigns read group information to each alignment in the BAM file,
    which is essential for downstream GATK tools and other analysis pipelines that
    require read group metadata. Read groups allow tracking of sequencing data
    from different libraries, lanes, or flow cells, enabling proper handling of
    technical artifacts and batch effects.

    Key read group fields assigned:
    - ID: Read group identifier (set to "1" for single sample processing)
    - LB: Library identifier (set to "lib1" for tracking library preparation)
    - PL: Platform/technology used (set to "illumina" for Illumina sequencing)
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
    Mark PCR duplicates using GATK MarkDuplicates and configure aligner-specific duplicate handling.

    This rule identifies and marks PCR duplicates in the BAM file, which are reads that
    originate from the same DNA fragment due to PCR amplification during library preparation.
    Duplicate marking is essential for accurate downstream analyses as PCR duplicates can
    skew coverage estimates, affect variant calling, and introduce bias in peak detection
    for ATAC-seq experiments.

    For aligner-specific duplicate handling:
    - When using Chromap as the aligner, this rule skips GATK MarkDuplicates because
      Chromap performs duplicate marking internally during the alignment process.
      The rule creates symbolic links to the input BAM file and generates a metrics
      file indicating that deduplication was handled by Chromap.
    - When using Bowtie2 or other aligners, this rule runs GATK MarkDuplicates
      to identify and mark duplicate reads based on alignment coordinates and
      read sequences.

    Key MarkDuplicates parameters:
    - CREATE_INDEX: Generates index for the output BAM file for rapid access
    - MAX_RECORDS_IN_RAM: Sets memory buffer size for sorting (5 million records)
    - SORTING_COLLECTION_SIZE_RATIO: Controls temporary file usage during sorting

    The duplicate marking process examines read pairs to identify those with
    identical alignment coordinates, which likely represent PCR duplicates of
    the same original DNA fragment. One read from each duplicate set is
    retained as the primary alignment, while others are marked with the
    DUPLICATE flag for downstream filtering.

    The metrics file generated by MarkDuplicates provides statistics about
    duplication rates, library complexity, and optical duplicates, which
    are valuable quality control metrics for assessing library preparation
    success and determining whether additional sequencing would be beneficial.
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
        java_opts = get_java_opts,  
        aligner = config.get("mapping_tools", "bowtie2"),
    shell:
        """
        if [ "{params.aligner}" = "chromap" ]; then
            echo "Using Chromap: Skipping GATK MarkDuplicates..." > {log}
        
            INPUT_ABS=$(readlink -f {input.bam})
            OUTPUT_ABS=$(readlink -f {output.bam})
            
            ln -sf "$INPUT_ABS" {output.bam}
            
            echo "MarkDuplicates skipped. Deduplication was automatically done by Chromap." > {output.metrics}

            if [ -f "{input.bam}.bai" ]; then
                ln -sf "$INPUT_ABS.bai" {output.bam}.bai
            else
                INPUT_PREFIX="${{INPUT_ABS%.bam}}"
                OUTPUT_PREFIX="${{OUTPUT_ABS%.bam}}"
                if [ -f "$INPUT_PREFIX.bai" ]; then
                    ln -sf "$INPUT_PREFIX.bai" "$OUTPUT_PREFIX.bai"
                fi
            fi
        else
            echo "Running GATK MarkDuplicates..." > {log}
            gatk --java-options "{params.java_opts}" MarkDuplicates \
                 -I {input.bam} \
                 -O {output.bam} \
                 -M {output.metrics} \
                 --CREATE_INDEX true \
                 --MAX_RECORDS_IN_RAM 5000000 \
                 --SORTING_COLLECTION_SIZE_RATIO 0.5 2>> {log}
        fi
        """

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
    Perform ATAC-seq specific read shifting to correct for Tn5 transposase binding bias.

    This rule applies a critical correction to ATAC-seq alignments to account for the
    binding properties of the Tn5 transposase, which creates a 9bp duplication during
    the tagmentation reaction. This shift is essential for accurate determination of
    transcription factor binding sites and open chromatin region boundaries.

    The Tn5 transposase binds as a dimer and inserts adapters separated by 9 base pairs.
    Consequently:
    - Reads mapping to the positive strand are shifted +4 bp
    - Reads mapping to the negative strand are shifted -5 bp
    - This centers the read pileup over the actual Tn5 insertion site

    This rule uses deepTools' alignmentSieve with the --ATACshift parameter, which
    automatically applies this standard ATAC-seq correction. The shifted BAM file
    is then re-sorted and indexed to maintain compatibility with downstream tools.

    Proper shifting is critical for ATAC-seq analysis because:
    - It improves the resolution of transcription factor footprinting analysis
    - It ensures accurate calling of peak boundaries
    - It enables proper TSS enrichment calculations
    - It provides more precise nucleosome positioning information
    - It is required for TOBIAS footprinting and other high-resolution analyses

    The shifted BAM file serves as the primary input for peak calling, bigwig
    generation, TSS enrichment analysis, and all subsequent downstream analyses
    in the ATAC-seq workflow.
    """
    input:
        bam = '02.mapping/filter_pe/{sample}.filter_pe.sorted.bam',
        bai = '02.mapping/filter_pe/{sample}.filter_pe.sorted.bam.bai' 
    output:
        shifted_bam = temp('02.mapping/shifted/{sample}.shifted.bam'),
        shifted_sort_bam = '02.mapping/shifted/{sample}.shifted.sorted.bam',
        shifted_sort_bam_bai = '02.mapping/shifted/{sample}.shifted.sorted.bam.bai'
    log:
        "logs/02.mapping/atac_seq_shift_{sample}.log"
    conda:
        workflow.source_path("../envs/deeptools.yaml")
    params:
        whether_shift = config.get("mapping_tools", "bowtie2") 
    threads:
        20
    shell:
        """
        # 使用 {{params.whether_shift}} 将 Snakemake 变量传给 Bash
        if [ "{params.whether_shift}" = "chromap" ]; then
            echo "Using Chromap: Skipping deepTools shift and creating symlinks..." > {log}
            ln -sf $(readlink -f {input.bam}) {output.shifted_bam}
            ln -sf $(readlink -f {input.bam}) {output.shifted_sort_bam}
            ln -sf $(readlink -f {input.bai}) {output.shifted_sort_bam_bai}
        else
            echo "Running deepTools alignmentSieve for Tn5 shifting..." > {log}
            alignmentSieve -b {input.bam} -o {output.shifted_bam} --ATACshift -p {threads} 2>> {log} && \
            samtools sort -@ {threads} {output.shifted_bam} -o {output.shifted_sort_bam} && \
            samtools index {output.shifted_sort_bam}
        fi
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