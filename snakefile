#!/usr/bin/env python3
# *---utf-8---*
# Version: RNAFlow_v0.1
# Author : JZHANG

import sys
import os
from snakemake.utils import min_version, validate

# ------- Import Custom Modules ------- #
from rules.utils.id_convert import load_samples, load_contrasts,parse_groups
from rules.utils.validate import check_reference_paths,load_user_config,validate_genome_version
from rules.utils.reference_update import resolve_reference_paths
from rules.utils.resource_manager import rule_resource

# Lock Snakemake Version
min_version("9.9.0")

# --------- 1. Config Loading --------- #
# Load default configs
configfile: "config/config.yaml"
configfile: "config/reference.yaml"
configfile: "config/run_parameter.yaml"
configfile: "config/cluster_config.yaml" 
# configfile: "config.yaml"

# Load CLI argument config (Highest Priority)
load_user_config(config, cmd_arg_name='analysisyaml')

# --------- 2. Processing & Validation --------- #
# Update absolute paths for references
resolve_reference_paths(config,
                        config.get('can_use_genome_version', []),
                        base_path=config.get('reference_path'))
# Validate schema and file existence
validate(config, "schema/config.schema.yaml")
check_reference_paths(config.get("STAR_index", {}))
# Get logger instance for validation
from snakemake_logger_plugin_rich_loguru import get_analysis_logger
logger = get_analysis_logger()
validate_genome_version(config=config, logger=logger)
# --------- 3. Workspaces & Samples --------- #
workdir: config["workflow"]
merge_group,samples = load_samples(config["sample_csv"], required_cols=["sample", "sample_name", "group"])
groups = parse_groups(samples)
ALL_CONTRASTS, CONTRAST_MAP = load_contrasts(config["paired_csv"], samples)
# --------- 4. Rules Import --------- #
include: 'rules/00.log.smk'
include: 'rules/01.common.smk'
include: 'rules/03.file_convert_md5.smk'
include: 'rules/04.short_read_qc.smk'
include: 'rules/05.Contamination_check.smk'
include: 'rules/06.short_read_clean.smk'
include: 'rules/07.mapping.smk'
include: 'rules/08.MACS2.smk'
include: 'rules/09.MergeMACS2.smk'

# 传递对比组信息给motif规则
exec("""
try:
    wildcard_constraints(
        contrast="|".join(ALL_CONTRASTS) if ALL_CONTRASTS else ".*",
        treatment="|".join(CONTRAST_MAP) if CONTRAST_MAP else ".*"
    )
except:
    pass
""")

include: 'rules/11.motifs.smk'
include: 'rules/14.Merge_qc.smk'
# --------- 5. Target Rule --------- #
rule all:
    input:
        DataDeliver(config = config,merge_group = merge_group,groups = groups)

# 新增：专门用于运行差异motif分析的目标规则
rule differential_motif_analysis:
    input:
        expand("06.motif_analysis/06.final_report/differential_motif_analysis_report.txt",
               zip, contrast=ALL_CONTRASTS, treatment=CONTRAST_MAP) if ALL_CONTRASTS and CONTRAST_MAP else []