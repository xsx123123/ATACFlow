#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import sys
import time
from pathlib import Path
from typing import Dict, Union, List, Callable
from rich import print as rich_print
from snakemake.io import expand
from loguru import logger


from utils.datadeliver import (
    qc_clean,
    mapping,
    peak_calling,
    motif_analysis,
    consensus_peaks,
    diff_peaks,
    merge_group_analysis,
    atac_qc,
    get_diff_analysis_input,
    should_run_pooled_analysis,
)


def DataDeliver(
    config: Dict = None,
    samples: Dict = None,
    merge_group: bool = False,
    groups: Dict = None,
    run_pooled: bool = False,
) -> List[str]:
    """
    Main data delivery orchestrator for ATAC-seq pipeline.
    Controls the flow of the pipeline based on 'only_qc' and specific module flags.
    """
    config = config or {}
    samples = samples or {}
    groups = groups or {}

    convert_md5_path = config.get("convert_md5", "md5_check")
    data_deliver = [
        "01.qc/md5_check.tsv",
        os.path.join("00.raw_data", convert_md5_path),
        os.path.join("00.raw_data", convert_md5_path, "raw_data_md5.json"),
    ]

    if config.get("parameter", {}).get("fastq_screen", {}).get("run"):
        data_deliver.extend(
            expand(
                "01.qc/fastq_screen_r1/{sample}_R1_screen.txt", sample=samples.keys()
            )
        )
        data_deliver.extend(
            expand(
                "01.qc/fastq_screen_r2/{sample}_R2_screen.txt", sample=samples.keys()
            )
        )
        data_deliver.append(
            "01.qc/fastq_screen_multiqc_r1/multiqc_r1_fastq_screen_report.html"
        )
        data_deliver.append(
            "01.qc/fastq_screen_multiqc_r2/multiqc_r2_fastq_screen_report.html"
        )

    basic_modules = [
        "qc_clean",
        "mapping",
        "peak_calling",
        "motif_analysis",
        "consensus_peaks",
        "atac_qc",
    ]
    deep_qc_flags = ["bamCoverage"]
    downstream_modules = ["diff_peaks"]

    def execute_qc_clean(samples, data_deliver):
        return qc_clean(samples, data_deliver)

    def execute_mapping(samples, data_deliver, config):
        return mapping(samples, data_deliver, config)

    def execute_peak_calling(samples, data_deliver):
        return peak_calling(samples, data_deliver)

    def execute_motif_analysis(samples, data_deliver):
        return motif_analysis(samples, data_deliver)

    def execute_consensus_peaks(samples, data_deliver):
        return consensus_peaks(samples, data_deliver)

    def execute_diff_peaks(samples, data_deliver):
        return diff_peaks(samples, data_deliver, run_pooled)

    def execute_atac_qc(samples, data_deliver):
        return atac_qc(samples, data_deliver, run_pooled)

    def execute_merge_group_analysis(groups, data_deliver):
        return merge_group_analysis(groups, data_deliver)

    module_functions: Dict[str, Callable] = {
        "qc_clean": execute_qc_clean,
        "mapping": execute_mapping,
        "peak_calling": execute_peak_calling,
        "motif_analysis": execute_motif_analysis,
        "consensus_peaks": execute_consensus_peaks,
        "diff_peaks": execute_diff_peaks,
        "atac_qc": execute_atac_qc,
    }

    global _qc_warning_logged

    # Use a local dictionary to track enabled modules to avoid side effects on global config
    run_modules = {}

    for module in basic_modules:
        run_modules[module] = config.get(module) is not False

    for flag in deep_qc_flags:
        run_modules[flag] = config.get(flag) is not False

    if config.get("only_qc"):
        if not _qc_warning_logged:
            logger.warning(
                "**********************************************************************"
            )
            logger.warning(
                "   [MODE] ONLY QC ENABLED                                            "
            )
            logger.warning(
                "   - Running: Raw QC, Mapping, Peak Calling, Motif Analysis, QC      "
            )
            logger.warning(
                "   - Skipping: Differential Peak Analysis                              "
            )
            logger.warning(
                "**********************************************************************"
            )
            time.sleep(1)
            _qc_warning_logged = True

        for module in downstream_modules:
            run_modules[module] = False

    else:
        for module in downstream_modules:
            run_modules[module] = config.get(module) is not False

    for module, func in module_functions.items():
        if run_modules.get(module):
            if module == "mapping":
                data_deliver = func(samples, data_deliver, config)
            else:
                data_deliver = func(samples, data_deliver)

    if run_pooled:
        data_deliver = execute_merge_group_analysis(groups, data_deliver)

    if config.get("print_target"):
        rich_print("[bold green]Generated Target Files:[/bold green]")
        rich_print(data_deliver)

    return data_deliver


def ReportData(config: dict = None) -> List[str]:
    """
    Collects all files required for generating the final Quarto ATAC-seq analysis report.
    """
    if config.get("report"):
        return [
            os.path.join(config["data_deliver"], "delivery_manifest.json"),
            os.path.join(config["data_deliver"], "delivery_manifest.md5"),
            os.path.join(config["data_deliver"], "delivery_details.log"),
            os.path.join(config["data_deliver"], "report_data/project_summary.json"),
            os.path.join(
                config["data_deliver"], "report_data", "delivery_manifest.json"
            ),
            os.path.join(
                config["data_deliver"], "report_data", "delivery_manifest.md5"
            ),
            os.path.join(config["data_deliver"], "report_data", "delivery_details.log"),
            os.path.join(config["data_deliver"], "Analysis_Report/index.html"),
        ]
    else:
        return []


def get_sample_data_dir(sample_id: str = None, config: dict = None) -> str:
    """
    根据 sample_id 查找包含 fastq 文件的目录。

    逻辑更新：
    1. 优先查找是否存在以 sample_id 命名的【子目录】。
    2. 如果子目录不存在，则查找该目录下是否存在以 sample_id 开头的【文件】。
       如果存在文件，则返回该 base_dir。
    """

    # 确保 config 里有这个 key，防止报错
    if "raw_data_path" not in config:
        raise ValueError("Config dictionary missing 'raw_data_path' key.")

    # 遍历配置中的所有原始数据路径
    for base_dir in config["raw_data_path"]:
        # --- 情况 A: 也就是你之前的逻辑 (raw_data/SampleID/xxx.fq) ---
        sample_subdir = os.path.join(base_dir, sample_id)
        if os.path.isdir(sample_subdir):
            return sample_subdir

        # --- 情况 B: 也就是你现在的 ls 结果 (raw_data/SampleID.R1.fq) ---
        # 我们使用 glob 模糊匹配：查看该目录下是否有以 sample_id 开头的文件
        # pattern 类似于: /data/.../00.raw_data/L1MKK1806607-a1*
        pattern = os.path.join(base_dir, f"{sample_id}*")

        # 获取匹配的文件列表
        matching_files = glob.glob(pattern)

        # 只要找到了匹配的文件（并且是文件而不是文件夹），就说明数据在 base_dir 这一层
        if matching_files:
            # 简单的过滤：确保找到的是文件 (防止碰巧有一个叫 SampleID_tmp 的文件夹)
            # 只要有一个是文件，我们就认为找到了
            if any(os.path.isfile(f) for f in matching_files):
                return base_dir

    # 如果循环结束还没找到
    raise FileNotFoundError(
        f"无法在 {config['raw_data_path']} 中找到 {sample_id} 的数据目录或文件"
    )


def get_all_input_dirs(sample_keys: str = None, config: dict = None) -> list:
    """
    遍历所有样本 ID，调用 get_sample_data_dir，
    返回一个包含所有数据目录的列表。
    """
    dir_list = []
    for sample_id in sample_keys:
        dir_list.append(get_sample_data_dir(sample_id, config=config))

    return list(set(dir_list))


def judge_bwa_index(config: dict = None) -> bool:
    """
    判断是否需要重新构建bwa索引
    """
    bwa_index = config["bwa_mem2"]["index"]
    bwa_index_files = [
        bwa_index + suffix
        for suffix in [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac", ".alt"]
    ]

    return not all(os.path.exists(f) for f in bwa_index_files)


def judge_star_index(config: dict, Genome_Version: str) -> bool:
    """
    判断是否需要重新构建 STAR 索引
    Returns:
        True: 文件缺失，需要构建
        False: 文件完整，不需要构建
    """

    try:
        star_config = config["STAR_index"][Genome_Version]
        index_dir = star_config["index"]
    except KeyError:
        print(
            f"Error: Genome Version '{Genome_Version}' not found in config or structure incorrect."
        )
        sys.exit(1)

    if not os.path.isdir(index_dir):
        return True

    required_files = [
        "chrLength.txt",
        "exonGeTrInfo.tab",
        "genomeParameters.txt",
        "sjdbInfo.txt",
        "chrNameLength.txt",
        "exonInfo.tab",
        "Log.out",
        "sjdbList.fromGTF.out.tab",
        "chrName.txt",
        "geneInfo.tab",
        "SA",
        "sjdbList.out.tab",
        "chrStart.txt",
        " Genome",
        "SAindex",
        "transcriptInfo.tab",
    ]

    full_paths = [os.path.join(index_dir, f) for f in required_files]

    missing_files = [f for f in full_paths if not os.path.exists(f)]

    if missing_files:
        return True

    return False


def check_gene_version(config: dict = None, logger=None) -> None:
    """
    Check if the gene version in config matches allowed list.
    """
    # Use the provided logger or get the unified logger
    if logger is None:
        from snakemake_logger_plugin_rich_loguru import get_analysis_logger

        logger = get_analysis_logger()

    try:
        version = config["Genome_Version"]
        allowed = config["can_use_genome_version"]

        if version not in allowed:
            logger.error(f"Version mismatch! '{version}' is not in {allowed}")
            raise ValueError(f"Unsupported genome version: {version}")

        logger.info(f"Config check passed: Genome_Version '{version}' is supported.")

    except KeyError as e:
        logger.error(f"Config structure error: Missing key {e}")
        raise
    except TypeError:
        logger.error("Config must be a valid dictionary.")
        raise
