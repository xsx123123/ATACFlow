#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
rule merge_DEG:
    input:
        counts = "04.consensus/merge_consensus_counts_matrix_ann.txt",
    output:
        output = '06.deg_enrich/DEG_merge/Global_PCA.pdf',
        summary = '06.deg_enrich/DEG_merge/All_Contrast_Differential_Peaks_Statistics.csv',
    resources:
        **rule_resource(config, 'low_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/deg_deseq2.yaml"),
    log:
        "logs/06.DEG/deseq2_peak.log",
    benchmark:
        "benchmarks/deseq2_peak_benchmark.txt",
    params:
        deg_dir = '06.deg_enrich/DEG_merge',
        samples = config['sample_csv'],
        paired = config['paired_csv'],
        PATH = workflow.source_path(config['parameter']['DEG']['PATH']),
        LFC = config['parameter']['DEG']['LFC'],
        PVAL = config['parameter']['DEG']['PVAL'],
    threads: 
        1
    shell:
        """
        chmod +x {params.PATH} && \
        Rscript {params.PATH} -c {input.counts} \
                -m {params.samples} \
                -p {params.paired} \
                -o {params.deg_dir} \
                --lfc={params.LFC} \
                --pval={params.PVAL} \
                --label_co=gene_id &> {log}
        """

rule merge_Enrichments:
    input:
        DEG_info = "06.deg_enrich/DEG_merge/All_Contrast_Differential_Peaks_Statistics.csv",
    output:
        Enrichments_dir = directory("06.deg_enrich/merge_enrich/"),
    resources:
        **rule_resource(config, 'low_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/go_enrich_r.yaml"),
    log:
        "logs/06.deg_enrich/go_enrich.log",
    params:
        obo = config['Bowtie2_index']['GO']['obo'],
        go_annotation = config['Bowtie2_index'][config['Genome_Version']]['go_annotation'],
        gene_col = config['parameter']['Enrichments']['gene_col'],
        r_script = workflow.source_path(config['parameter']['Enrichments']['PATH']),
        wrapper = workflow.source_path(config['parameter']['Enrichments']['PATH_py']),
        lib_type = config['parameter']['Enrichments']['lib_type'],
        gene_regex = config['parameter']['Enrichments']['gene_regex'],
        deg_dir = "06.deg_enrich/DEG_merge",
        cutoff = config['parameter']['Enrichments'].get('cutoff', 0.05)
    shell:
        """
        python {params.wrapper} \
            --rscript {params.r_script} \
            --deg_info {input.DEG_info} \
            --deg_dir {params.deg_dir} \
            --lib_type {params.lib_type} \
            -o {params.obo} \
            -a {params.go_annotation} \
            -d {output.Enrichments_dir} \
            --gene_col {params.gene_col} \
            --gene_regex '{params.gene_regex}' \
            --cutoff {params.cutoff} > {log} 2>&1
        """
# ------- rule ------- #
