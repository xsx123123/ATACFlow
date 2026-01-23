
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
        json = "05.qc/ataqv/{sample}.ataqv.json",
        log_out = "05.qc/ataqv/{sample}.ataqv.out"
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

rule mkarv_report:
    """
    Merge all JSONs into one HTML report.
    """
    input:
        jsons = expand("05.qc/ataqv/{sample}.ataqv.json.gz", sample=samples.keys())
    output:
        directory("05.qc/ataqv_report")
    conda:
        workflow.source_path("../envs/ataqv.yaml")
    log:
        "logs/05.qc/mkarv.log"
    shell:
        """
        mkarv 05.qc/ataqv_report {input.jsons} > {log} 2>&1
        """
# ----- end of rules ----- #