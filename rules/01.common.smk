#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# loading packages
import os
import glob
from pathlib import Path
from typing import Dict, Union
from rich import print as rich_print
# Target rule function
def DataDeliver(config:dict = None,merge_group = False,groups:dict = None) -> list:
    """
    This function performs Bioinformation analysis on the input configuration
    and returns a list of results.
    """
    # short-read raw-data qc result
    data_deliver = ["01.qc/md5_check.tsv",
                    os.path.join('00.raw_data',config['convert_md5']),
                    os.path.join('00.raw_data',config['convert_md5'],"raw_data_md5.json")
                    ]
    # fastq-screen
    if config['parameter']['fastq_screen']['run']:
        data_deliver.extend(expand("01.qc/fastq_screen_r1/{sample}_R1_screen.txt",
                                          sample=samples.keys()))
        data_deliver.extend(expand("01.qc/fastq_screen_r2/{sample}_R2_screen.txt",
                                          sample=samples.keys()))
        data_deliver.append("01.qc/fastq_screen_multiqc_r1/multiqc_r1_fastq_screen_report.html")
        data_deliver.append("01.qc/fastq_screen_multiqc_r2/multiqc_r2_fastq_screen_report.html")
    # fastqc & multiqc
    data_deliver.append("01.qc/short_read_r1_multiqc/multiqc_r1_raw-data_report.html")
    data_deliver.append("01.qc/short_read_r2_multiqc/multiqc_r2_raw-data_report.html")        
    # short-read trim & clean result
    data_deliver.append("01.qc/multiqc_short_read_trim/multiqc_short_read_trim_report.html")
    # merge qc
    data_deliver.append("01.qc/multiqc_merge_qc/multiqc_merge_qc_report.html")
    # mapping
    data_deliver.extend(expand('02.mapping/Bowtie2/{sample}/{sample}.sorted.bam',sample=samples.keys()))
    data_deliver.extend(expand('02.mapping/Bowtie2/{sample}/{sample}.sorted.bam.bai',sample=samples.keys()))
    data_deliver.extend(expand('02.mapping/preseq/{sample}.lc_extrap.txt',sample=samples.keys())),
    data_deliver.extend(expand('02.mapping/preseq/{sample}.c_curve.txt',sample=samples.keys())),
    data_deliver.extend(expand('02.mapping/filter_pe/{sample}.filter_pe.sorted.bam',sample=samples.keys())),
    data_deliver.extend(expand('02.mapping/filter_pe/{sample}.filter_pe.sorted.bam.bai',sample=samples.keys())),
    data_deliver.extend(expand('02.mapping/shifted/{sample}.shifted.sorted.bam',sample=samples.keys())),
    data_deliver.extend(expand('02.mapping/shifted/{sample}.shifted.sorted.bam.bai',sample=samples.keys())),
    data_deliver.extend(expand(f"02.mapping/bamCoverage/{{sample}}_{config['parameter']['bamCoverage']['normalizeUsing']}.bw",
                                          sample=samples.keys()))
    data_deliver.extend(expand('02.mapping/computeMatrix/{sample}_TSS_matrix.gz',sample=samples.keys())),
    data_deliver.extend(expand('02.mapping/plots/{sample}_TSS_enrichment.png',sample=samples.keys())),
    data_deliver.append(expand('02.mapping/samtools_flagstat/{sample}_bam_flagstat.tsv',sample=samples.keys()))
    data_deliver.append(expand('02.mapping/samtools_stats/{sample}_bam_stats.tsv',sample=samples.keys()))
    # macs2
    data_deliver.extend(expand('03.peak_calling/MACS2/{sample}/{sample}_peaks.narrowPeak',sample=samples.keys()))
    data_deliver.extend(expand('03.peak_calling/MACS2/{sample}/{sample}_peaks.xls',sample=samples.keys()))
    data_deliver.extend(expand('03.peak_calling/MACS2/{sample}/{sample}_summits.bed',sample=samples.keys()))
    data_deliver.extend(expand('03.peak_calling/MACS2/{sample}/{sample}_treat_pileup.bdg',sample=samples.keys()))
    # HOMER
    data_deliver.append(expand("03.peak_calling/HOMER/{sample}_annotation.txt",sample=samples.keys()))
    data_deliver.append(expand("03.peak_calling/HOMER/{sample}_stats.txt",sample=samples.keys()))
    # count_matrix
    data_deliver.append("04.consensus/consensus_counts_matrix.txt")
    data_deliver.append("04.consensus/matrix_description.txt")
    data_deliver.append("04.consensus/consensus_counts_matrix_ann.txt")
    # ataqv
    data_deliver.append("05.qc/ataqv_report")


    # ------------- merge group ---------------- #
    if merge_group:
        data_deliver.extend(expand("02.mapping/merged/{group}.merged.bam",group = groups.keys()))
        data_deliver.extend(expand("02.mapping/merged/{group}.merged.bam.bai",group = groups.keys()))
        data_deliver.extend(expand("03.peak_calling/MERGE_MACS2/{group}/{group}_peaks.narrowPeak",group = groups.keys()))
        data_deliver.extend(expand("03.peak_calling/MERGE_MACS2/{group}/{group}_peaks.xls",group = groups.keys()))
        data_deliver.extend(expand("03.peak_calling/MERGE_MACS2/{group}/{group}_summits.bed",group = groups.keys()))
        data_deliver.extend(expand("03.peak_calling/MERGE_MACS2/{group}/{group}_treat_pileup.bdg",group = groups.keys()))
        data_deliver.extend(expand("03.peak_calling/MERGE_HOMER/{group}_annotation.txt",group = groups.keys()))
        data_deliver.extend(expand("03.peak_calling/MERGE_HOMER/{group}_stats.txt",group = groups.keys()))
        data_deliver.append("04.consensus/merge_matrix_description.txt")
        data_deliver.append("04.consensus/merge_consensus_counts_matrix.txt")
        data_deliver.append("04.consensus/merge_consensus_counts_matrix_ann.txt")
    # deg
    if merge_group:
        data_deliver.append('06.deg_enrich/DEG_merge/Global_PCA.pdf')
        data_deliver.append('06.deg_enrich/DEG_merge/All_Contrast_Differential_Peaks_Statistics.csv')
        data_deliver.append('06.deg_enric/merge_enrich/')
    else:
        data_deliver.append('06.deg_enrich/DEG/Global_PCA.pdf')
        data_deliver.append('06.deg_enrich/DEG/All_Contrast_Differential_Peaks_Statistics.csv')
        data_deliver.append('06.deg_enric/enrich/')
    
    # ATAC Report
    data_deliver.append('05.ATAC_QC/multiqc_ATAC_report.html')

    if config['print_target']:
       rich_print(data_deliver)
    return  data_deliver

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
    raise FileNotFoundError(f"无法在 {config['raw_data_path']} 中找到 {sample_id} 的数据目录或文件")

def get_all_input_dirs(sample_keys:str = None,
                       config:dict = config) -> list:
    """
    遍历所有样本 ID，调用 get_sample_data_dir，
    返回一个包含所有数据目录的列表。
    """
    dir_list = []
    for sample_id in sample_keys:
        dir_list.append(get_sample_data_dir(sample_id,config = config))

    return list(set(dir_list))

def judge_bwa_index(config:dict = None) -> bool:
    """
    判断是否需要重新构建bwa索引
    """
    bwa_index = config['bwa_mem2']['index']
    bwa_index_files = [bwa_index + suffix for suffix in ['.0123', '.amb', '.ann', '.bwt.2bit.64', '.pac', '.alt']]
    
    return not all(os.path.exists(f) for f in bwa_index_files)

def judge_star_index(config: dict, Genome_Version: str) -> bool:
    """
    判断是否需要重新构建 STAR 索引
    Returns:
        True: 文件缺失，需要构建
        False: 文件完整，不需要构建
    """

    try:
        star_config = config['STAR_index'][Genome_Version]
        index_dir = star_config['index']
    except KeyError:
        print(f"Error: Genome Version '{Genome_Version}' not found in config or structure incorrect.")
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
        "transcriptInfo.tab"
    ]
    
    full_paths = [os.path.join(index_dir, f) for f in required_files]

    missing_files = [f for f in full_paths if not os.path.exists(f)]
    
    if missing_files:
        return True
    
    return False

def check_gene_version(config: dict = None, logger = None) -> None:
    """
    Check if the gene version in config matches allowed list.
    """
    # Use the provided logger or get the unified logger
    if logger is None:
        from snakemake_logger_plugin_rich_loguru import get_analysis_logger
        logger = get_analysis_logger()

    try:
        version = config['Genome_Version']
        allowed = config['can_use_genome_version']

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
        
# --------------------- #
