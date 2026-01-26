#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
import os
# ----- rule ----- #
rule merge_qc_report:
    input:
        md5_check = "01.qc/md5_check.tsv",
        fastqc_files_r1 = expand("01.qc/short_read_qc_r1/{sample}_R1_fastqc.zip", sample=samples.keys()),
        fastqc_files_r2 = expand("01.qc/short_read_qc_r2/{sample}_R2_fastqc.zip", sample=samples.keys()),
        r1_trimmed = expand("01.qc/short_read_trim/{sample}.R1.trimed.fq.gz", sample=samples.keys()),
        r2_trimmed = expand("01.qc/short_read_trim/{sample}.R2.trimed.fq.gz", sample=samples.keys()),
        fastp_report = expand("01.qc/short_read_trim/{sample}.trimed.html", sample=samples.keys()),
    output:
        report = '01.qc/multiqc_merge_qc/multiqc_merge_qc_report.html',
    resources:
        **rule_resource(config, 'low_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/multiqc.yaml"),
    params:
        search_dir = "01.qc",
        out_dir = "01.qc/multiqc_merge_qc",
        report_name = "multiqc_merge_qc_report.html",
        title = "merge qc report",
    log:
        "logs/01.qc/multiqc_merge_qc.log",
    threads:
        config['parameter']['threads']['multiqc'],
    shell:
        """
        mkdir -p {params.out_dir}
        multiqc ./{params.search_dir} \
                --force \
                --outdir {params.out_dir} \
                --title "{params.title}" \
                --filename {params.report_name} \
                &> {log}
        """
# ----- rule ----- #