#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
import os

rule generate_fastq_screen_conf:
    """
    Generate fastq_screen configuration file with database path
    """
    input:
        template = workflow.source_path(config['parameter']['validate_fastq_screen']['path_conf']),
    output:
        conf = "01.qc/fastq_screen.conf"
    params:
        db_path = config.get("fastq_screen_db_path", "/data/jzhang/reference/"),
    log:
        "logs/01.qc/generate_fastq_screen_conf.log",
    benchmark:
        "benchmarks/01.qc/generate_fastq_screen_conf.txt",
    localrule: True
    shell:
        """
        sed "s|__FASTQ_SCREEN_DB_PATH__|{params.db_path}|g" {input.template} > {output.conf} 2> {log}
        """

rule check_fastq_screen_conf:
    """
    Validate fastq_screen configuration file
    """
    input:
        conf = "01.qc/fastq_screen.conf",
    output:
        log = "01.qc/fastq_screen_config_check.log",
    resources:
        **rule_resource(config, 'low_resource',  skip_queue_on_local=True,logger = logger),
    params:
        validate_fastq_screen = workflow.source_path(config['parameter']['validate_fastq_screen']['path']),
    log:
        "logs/01.qc/check_fastq_screen_conf.log",
    benchmark:
        "benchmarks/01.qc/check_fastq_screen_conf.txt",
    threads:
        1
    shell:
        """
        chmod +x {params.validate_fastq_screen} && \
        {params.validate_fastq_screen} {input.conf} --log {output.log} &> {log}
        """

rule short_read_fastq_screen_r1:
    """
    Check contamination in R1 reads using fastq_screen
    """
    input:
        md5_check = "01.qc/md5_check.tsv",
        log = "01.qc/fastq_screen_config_check.log",
        conf = "01.qc/fastq_screen.conf"
    output:
        fastq_screen_result = "01.qc/fastq_screen_r1/{sample}_R1_screen.txt",
    resources:
        **rule_resource(config, 'high_resource',  skip_queue_on_local=True,logger = logger),
    log:
        "logs/01.qc/fastq_screen_r1_{sample}.log",
    conda:
        workflow.source_path('../envs/fastq_screen.yaml'),
    params:
        out_dir = "01.qc/fastq_screen_r1/",
        link_r1_dir = os.path.join("00.raw_data",
                                      config['convert_md5'],
                                      "{sample}/{sample}_R1.fq.gz"),
        subset = config['parameter'][ 'fastq_screen']['subset'],
        aligner = config['parameter']['fastq_screen']['aligner'],
    message:
        "Running fastq_screen on {wildcards.sample} r1",
    benchmark:
        "benchmarks/01.qc/fastq_screen_r1_{sample}.txt",
    threads:
        config['parameter']['threads']['fastq_screen'],
    shell:
        """
        fastq_screen --threads  {threads} \
                     --force \
                     --subset  {params.subset} \
                     --aligner  {params.aligner} \
                     --conf {input.conf} \
                     --outdir {params.out_dir} \
                     {params.link_r1_dir} &> {log}
        """

rule short_read_fastq_screen_r2:
    """
    Check contamination in R2 reads using fastq_screen
    """
    input:
        md5_check = "01.qc/md5_check.tsv",
        log = "01.qc/fastq_screen_config_check.log",
        conf = "01.qc/fastq_screen.conf"
    output:
        fastq_screen_result = "01.qc/fastq_screen_r2/{sample}_R2_screen.txt",
    resources:
        **rule_resource(config, 'high_resource',  skip_queue_on_local=True,logger = logger),
    log:
        "logs/01.qc/fastq_screen_r2_{sample}.log",
    conda:
        workflow.source_path('../envs/fastq_screen.yaml'),
    params:
        out_dir = "01.qc/fastq_screen_r2/",
        subset = config['parameter'][ 'fastq_screen']['subset'],
        aligner = config['parameter']['fastq_screen']['aligner'],
        link_r2_dir = os.path.join("00.raw_data",
                                      config['convert_md5'],
                                      "{sample}/{sample}_R2.fq.gz"),
    message:
        "Running fastq_screen on {wildcards.sample} r2",
    benchmark:
        "benchmarks/01.qc/fastq_screen_r2_{sample}.txt",
    threads:
        config['parameter']['threads']['fastq_screen'],
    shell:
        """
        fastq_screen --threads  {threads} \
                     --force \
                     --subset  {params.subset} \
                     --aligner  {params.aligner} \
                     --conf {input.conf} \
                     --outdir {params.out_dir} \
                     {params.link_r2_dir} &> {log}
        """

rule fastq_screen_multiqc_r1:
    """
    Run MultiQC to aggregate R1 fastq screen reports
    """
    input:
        fastqc_files_r1 = expand("01.qc/fastq_screen_r1/{sample}_R1_screen.txt",\
                                  sample=samples.keys()),
    output:
        report = "01.qc/fastq_screen_multiqc_r1/multiqc_r1_fastq_screen_report.html",
    resources:
        **rule_resource(config, 'low_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/multiqc.yaml"),
    message:
        "Running MultiQC to aggregate R1 fastq screen reports",
    params:
        fastqc_reports = "01.qc/fastq_screen_r1",
        report_dir = "01.qc/fastq_screen_multiqc_r1/",
        report = "multiqc_r1_fastq_screen_report.html",
        title = "r1-fastq-screen-multiqc-report",
    log:
        "logs/01.qc/fastq_screen_multiqc_r1.log",
    benchmark:
        "benchmarks/01.qc/fastq_screen_multiqc_r1.txt",
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

rule fastq_screen_multiqc_r2:
    """
    Run MultiQC to aggregate R2 fastq screen reports
    """
    input:
        fastqc_files_r1 = expand("01.qc/fastq_screen_r2/{sample}_R2_screen.txt",
                                sample=samples.keys()),
    output:
        report = "01.qc/fastq_screen_multiqc_r2/multiqc_r2_fastq_screen_report.html",
    resources:
        **rule_resource(config, 'low_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/multiqc.yaml"),
    message:
        "Running MultiQC to aggregate R2 fastq screen reports",
    params:
        fastqc_reports = "01.qc/fastq_screen_r2",
        report_dir = "01.qc/fastq_screen_multiqc_r2/",
        report = "multiqc_r2_fastq_screen_report.html",
        title = "r2-fastq-screen-multiqc-report",
    log:
        "logs/01.qc/fastq_screen_multiqc_r2.log",
    benchmark:
        "benchmarks/01.qc/fastq_screen_multiqc_r2.txt",
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
# ----- end of rules ----- #
