
def get_organelle_filter_expr(wildcards):
    """
    获取线粒体染色体名称 (chrMID)
    """
    build = config.get("Genome_Version")
    chrMID = config.get("genome_info", {}).get(build, {}).get("chrMID", {})

    if isinstance(chrMID, list):
        # Return the first element if it's a list
        return chrMID[0] if chrMID else ""
    else:
        # Return the value directly if it's not a list
        return chrMID

rule ataqv_qc:
    input:
        shifted_sort_bam = '02.mapping/shifted/{sample}.shifted.sorted.bam',
        shifted_sort_bam_bai = '02.mapping/shifted/{sample}.shifted.sorted.bam.bai',
        narrow_peak = "03.peak_calling/MACS2/{sample}/{sample}_peaks.narrowPeak",
    output:
        json = "02.mapping/ataqv/{sample}.ataqv.json",
        log_out = "02.mapping/ataqv/{sample}.ataqv.out"
    conda:
        workflow.source_path("../envs/ataqv.yaml"),
    resources:
        **rule_resource(config, 'high_resource', skip_queue_on_local=True, logger=logger),
    log:
        "logs/05.qc/ataqv/{sample}.ataqv.log",
    message:
        "Running ataqv QC on {wildcards.sample}",
    params:
        tss = config['Bowtie2_index'][config['Genome_Version']]['tss_bed'],
        autosomes = config['Bowtie2_index'][config['Genome_Version']]['autosomes'],
        mito_name = lambda wildcards:get_organelle_filter_expr(wildcards),
        organism = config['species']
    shell:
        """
        ataqv \
            --peak-file {input.narrow_peak} \
            --name {wildcards.sample} \
            --metrics-file {output.json} \
            --tss-file {params.tss} \
            --autosomal-reference-file {params.autosomes} \
            --mitochondrial-reference-name {params.mito_name} \
            --ignore-read-groups \
            {params.organism} \
            {input.shifted_sort_bam} > {output.log_out} 2> {log}
        """
rule multiqc_ATAC_QC:
    """
    Run MultiQC to aggregate ATAC QC reports
    """
    input:
        json = expand("02.mapping/ataqv/{sample}.ataqv.json",sample=samples.keys()),
        log_out = expand("02.mapping/ataqv/{sample}.ataqv.out",sample=samples.keys()),
        preseq = expand('02.mapping/preseq/{sample}.lc_extrap.txt',sample=samples.keys()),
        c_curve = expand('02.mapping/preseq/{sample}.c_curve.txt',sample=samples.keys()),
        samtools_flagstat = expand('02.mapping/samtools_flagstat/{sample}_bam_flagstat.tsv',sample=samples.keys()),
        samtools_stats = expand('02.mapping/samtools_stats/{sample}_bam_stats.tsv',sample=samples.keys()),
        metrics = expand('02.mapping/gatk/{sample}/{sample}.rg.dedup.metrics.txt',sample=samples.keys()),
    output:
        report = '05.ATAC_QC/multiqc_ATAC_report.html',
    resources:
        **rule_resource(config, 'low_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/multiqc.yaml"),
    message:
        "Running MultiQC to aggregate fastp reports",
    benchmark:
        "benchmarks/05.ATAC_QC/multiqc_ATAC_report.txt",
    params:
        origin_reports = "02.mapping/",
        report_dir = "05.ATAC_QC/",
        report = "multiqc_ATAC_report.html",
        title = "ATAC_report",
    log:
        "logs/05.ATAC_QC/multiqc_ATAC_report.log",
    threads:
        config['parameter']['threads']['multiqc'],
    shell:
        """
        multiqc {params.origin_reports} \
                --force \
                --outdir {params.report_dir} \
                -i {params.title} \
                -n {params.report} &> {log}
        """
# ----- end of rules ----- #