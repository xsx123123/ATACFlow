#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
"""
ATACFlow Pipeline - Data Delivery and Organization Module

This module handles the final organization, packaging, and delivery of ATAC-seq
analysis results using a Rust-accelerated delivery tool. It ensures that all
analysis outputs are properly organized, documented, and transferred to the
final delivery directory in a structured and reproducible manner.

Key Components:
- delivery: Organizes and transfers all analysis results to the final delivery directory
- delivery_report: Prepares a subset of results specifically for report generation

This module ensures that all analysis outputs are properly packaged for delivery,
including comprehensive manifest files, MD5 checksums for data integrity, and
detailed logs documenting the delivery process. This is essential for ensuring
reproducibility and providing clear documentation of all analysis results.
"""

import os
import pandas as pd

rule delivery:
    input:
        DataDeliver(config)
    output:
        manifest_json = os.path.join(config['data_deliver'],'delivery_manifest.json'),
        manifest_md5 = os.path.join(config['data_deliver'],'delivery_manifest.md5'),
        manifest_log = os.path.join(config['data_deliver'],'delivery_details.log'),
    resources:
        **rule_resource(config, 'high_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/py3.12.yaml"),
    params:
        out_dir = config['data_deliver'],
        config_path = workflow.source_path(config['parameter']['ATACFlow_Deliver_Tool']['config_path']),
        source_dir = config['workflow'],
    log:
        "logs/delivery.log",
    benchmark:
        "benchmark/delivery.txt",
    threads:
        config['parameter']['threads']['rnaflow-cli'],
    shell:
        """
        ( rnaflow-cli deliver \
                    -d {params.source_dir} \
                    -o {params.out_dir} \
                    -c {params.config_path} ) &>{log}
        """

rule delivery_report:
    input:
        ReportData(config)
    output:
        manifest_json = os.path.join(config['data_deliver'],'report_data','delivery_manifest.json'),
        manifest_md5 = os.path.join(config['data_deliver'],'report_data','delivery_manifest.md5'),
        manifest_log = os.path.join(config['data_deliver'],'report_data','delivery_details.log'),
    resources:
        **rule_resource(config, 'high_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/py3.12.yaml"),
    params:
        out_dir =  os.path.join(config['data_deliver'],'report_data'),
        config_path = workflow.source_path(config['parameter']['ATACFlow_Deliver_Tool']['config_path_report']),
        source_dir = config['workflow'],
    log:
        "logs/delivery_report.log",
    benchmark:
        "benchmark/delivery_report.txt",
    threads:
        config['parameter']['threads']['rnaflow-cli'],
    shell:
        """
        ( rnaflow-cli deliver \
                    -d {params.source_dir} \
                    -o {params.out_dir} \
                    -c {params.config_path}  ) &>{log}
        """