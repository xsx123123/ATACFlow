#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
from pathlib import Path
from typing import Dict, List, Union
from snakemake.io import expand


def qc_clean(samples: Dict = None, data_deliver: List = None) -> List:
    """
    Handle quality control and cleaning steps.

    Args:
        samples: Dictionary of sample information
        data_deliver: List of deliverable files to extend

    Returns:
        Updated list of deliverable files
    """
    if samples is None:
        samples = {}
    if data_deliver is None:
        data_deliver = []

    data_deliver.append("01.qc/short_read_r1_multiqc/multiqc_r1_raw-data_report.html")
    data_deliver.append("01.qc/short_read_r2_multiqc/multiqc_r2_raw-data_report.html")
    data_deliver.append(
        "01.qc/multiqc_short_read_trim/multiqc_short_read_trim_report.html"
    )
    data_deliver.append("01.qc/multiqc_merge_qc/multiqc_merge_qc_report.html")
    return data_deliver


def mapping(
    samples: Dict = None, data_deliver: List = None, config: Dict = None
) -> List:
    """
    Handle ATAC-seq sequence alignment and related outputs.

    Args:
        samples: Dictionary of sample information
        data_deliver: List of deliverable files to extend
        config: Configuration dictionary

    Returns:
        Updated list of deliverable files
    """
    if samples is None:
        samples = {}
    if data_deliver is None:
        data_deliver = []
    if config is None:
        config = {}

    data_deliver.extend(
        expand("02.mapping/cram/{sample}.cram", sample=samples.keys())
    )
    data_deliver.extend(
        expand(
            "02.mapping/cram/{sample}.cram.crai", sample=samples.keys()
        )
    )
    data_deliver.extend(
        expand("02.mapping/preseq/{sample}.lc_extrap.txt", sample=samples.keys())
    )
    data_deliver.extend(
        expand("02.mapping/preseq/{sample}.c_curve.txt", sample=samples.keys())
    )
    data_deliver.extend(
        expand(
            "02.mapping/filter_pe/{sample}.filter_pe.sorted.bam", sample=samples.keys()
        )
    )
    data_deliver.extend(
        expand(
            "02.mapping/filter_pe/{sample}.filter_pe.sorted.bam.bai",
            sample=samples.keys(),
        )
    )
    data_deliver.extend(
        expand("02.mapping/shifted/{sample}.shifted.sorted.bam", sample=samples.keys())
    )
    data_deliver.extend(
        expand(
            "02.mapping/shifted/{sample}.shifted.sorted.bam.bai", sample=samples.keys()
        )
    )

    if config.get("bamCoverage"):
        normalize_method = (
            config.get("parameter", {})
            .get("bamCoverage", {})
            .get("normalizeUsing", "RPKM")
        )
        for sample in samples.keys():
            data_deliver.append(
                f"02.mapping/bamCoverage/{sample}_{normalize_method}.bw"
            )

    data_deliver.extend(
        expand("02.mapping/computeMatrix/{sample}_TSS_matrix.gz", sample=samples.keys())
    )
    data_deliver.extend(
        expand("02.mapping/plots/{sample}_TSS_enrichment.png", sample=samples.keys())
    )
    data_deliver.extend(
        expand(
            "02.mapping/samtools_flagstat/{sample}_bam_flagstat.tsv",
            sample=samples.keys(),
        )
    )
    data_deliver.extend(
        expand(
            "02.mapping/samtools_stats/{sample}_bam_stats.tsv", sample=samples.keys()
        )
    )

    return data_deliver


def peak_calling(samples: Dict = None, data_deliver: List = None) -> List:
    """
    Handle ATAC-seq peak calling with MACS2.

    Args:
        samples: Dictionary of sample information
        data_deliver: List of deliverable files to extend

    Returns:
        Updated list of deliverable files
    """
    if samples is None:
        samples = {}
    if data_deliver is None:
        data_deliver = []

    data_deliver.extend(
        expand(
            "03.peak_calling/MACS2/{sample}/{sample}_peaks.narrowPeak",
            sample=samples.keys(),
        )
    )
    data_deliver.extend(
        expand(
            "03.peak_calling/MACS2/{sample}/{sample}_peaks.xls", sample=samples.keys()
        )
    )
    data_deliver.extend(
        expand(
            "03.peak_calling/MACS2/{sample}/{sample}_summits.bed", sample=samples.keys()
        )
    )
    data_deliver.extend(
        expand(
            "03.peak_calling/MACS2/{sample}/{sample}_treat_pileup.bdg",
            sample=samples.keys(),
        )
    )

    return data_deliver


def motif_analysis(samples: Dict = None, data_deliver: List = None) -> List:
    """
    Handle ATAC-seq motif analysis with HOMER.

    Args:
        samples: Dictionary of sample information
        data_deliver: List of deliverable files to extend

    Returns:
        Updated list of deliverable files
    """
    if samples is None:
        samples = {}
    if data_deliver is None:
        data_deliver = []

    data_deliver.extend(
        expand("03.peak_calling/HOMER/{sample}_annotation.txt", sample=samples.keys())
    )
    data_deliver.extend(
        expand("03.peak_calling/HOMER/{sample}_stats.txt", sample=samples.keys())
    )

    return data_deliver


def consensus_peaks(samples: Dict = None, data_deliver: List = None) -> List:
    """
    Handle consensus peak calling and count matrix generation.

    Args:
        samples: Dictionary of sample information
        data_deliver: List of deliverable files to extend

    Returns:
        Updated list of deliverable files
    """
    if samples is None:
        samples = {}
    if data_deliver is None:
        data_deliver = []

    data_deliver.append("04.consensus/consensus_counts_matrix.txt")
    data_deliver.append("04.consensus/matrix_description.txt")
    data_deliver.append("04.consensus/consensus_counts_matrix_ann.txt")

    return data_deliver


def diff_peaks(
    samples: Dict = None, data_deliver: List = None, merge_group: bool = False
) -> List:
    """
    Handle differential peak analysis.

    Args:
        samples: Dictionary of sample information
        data_deliver: List of deliverable files to extend
        merge_group: Whether to use merge group mode

    Returns:
        Updated list of deliverable files
    """
    if samples is None:
        samples = {}
    if data_deliver is None:
        data_deliver = []

    if merge_group:
        data_deliver.append("06.deg_enrich/DEG_merge/Global_PCA.pdf")
        data_deliver.append(
            "06.deg_enrich/DEG_merge/All_Contrast_Differential_Peaks_Statistics.csv"
        )
        data_deliver.append("06.deg_enrich/merge_enrich/")
    else:
        data_deliver.append("06.deg_enrich/DEG/Global_PCA.pdf")
        data_deliver.append(
            "06.deg_enrich/DEG/All_Contrast_Differential_Peaks_Statistics.csv"
        )
        data_deliver.append("06.deg_enrich/enrich/")

    return data_deliver


def merge_group_analysis(groups: Dict = None, data_deliver: List = None) -> List:
    """
    Handle merge group specific analysis.

    Args:
        groups: Dictionary of group information
        data_deliver: List of deliverable files to extend

    Returns:
        Updated list of deliverable files
    """
    if groups is None:
        groups = {}
    if data_deliver is None:
        data_deliver = []

    data_deliver.extend(
        expand("02.mapping/merged/{group}.merged.bam", group=groups.keys())
    )
    data_deliver.extend(
        expand("02.mapping/merged/{group}.merged.bam.bai", group=groups.keys())
    )
    data_deliver.extend(
        expand(
            "03.peak_calling/MERGE_MACS2/{group}/{group}_peaks.narrowPeak",
            group=groups.keys(),
        )
    )
    data_deliver.extend(
        expand(
            "03.peak_calling/MERGE_MACS2/{group}/{group}_peaks.xls", group=groups.keys()
        )
    )
    data_deliver.extend(
        expand(
            "03.peak_calling/MERGE_MACS2/{group}/{group}_summits.bed",
            group=groups.keys(),
        )
    )
    data_deliver.extend(
        expand(
            "03.peak_calling/MERGE_MACS2/{group}/{group}_treat_pileup.bdg",
            group=groups.keys(),
        )
    )
    data_deliver.extend(
        expand(
            "03.peak_calling/MERGE_HOMER/{group}_annotation.txt", group=groups.keys()
        )
    )
    data_deliver.extend(
        expand("03.peak_calling/MERGE_HOMER/{group}_stats.txt", group=groups.keys())
    )
    data_deliver.append("04.consensus/merge_matrix_description.txt")
    data_deliver.append("04.consensus/merge_consensus_counts_matrix.txt")
    data_deliver.append("04.consensus/merge_consensus_counts_matrix_ann.txt")
    data_deliver.append("05.ATAC_QC/multiqc_MACS2_Merge_report.html")

    return data_deliver


def atac_qc(
    samples: Dict = None, data_deliver: List = None, merge_group: bool = False
) -> List:
    """
    Handle ATAC-seq specific QC reports.

    Args:
        samples: Dictionary of sample information
        data_deliver: List of deliverable files to extend
        merge_group: Whether to use merge group mode

    Returns:
        Updated list of deliverable files
    """
    if samples is None:
        samples = {}
    if data_deliver is None:
        data_deliver = []

    data_deliver.append("05.ATAC_QC/multiqc_ATAC_report.html")

    if not merge_group:
        data_deliver.append("05.ATAC_QC/multiqc_MACS2_Samples_report.html")

    return data_deliver
