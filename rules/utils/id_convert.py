#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
from pathlib import Path
from rich import print as rprint
from typing import List, Dict, Tuple

def _validate_df(df: pd.DataFrame, required_cols: List[str], index_col: str) -> None:
    """[内部函数] 校验 DataFrame 的完整性和唯一性 (保持不变)"""
    try:
        from snakemake_logger_plugin_rich_loguru import get_analysis_logger
        logger = get_analysis_logger()
    except ImportError:
        import logging
        logger = logging.getLogger("Analysis")

    # 1. 校验必填列
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        rprint(f"[bold red]❌ 样本表格式错误！缺失列: {missing_cols}[/bold red]")
        sys.exit(1)

    # 2. 校验 ID 唯一性
    if df[index_col].duplicated().any():
        duplicated_ids = df[df[index_col].duplicated()][index_col].unique().tolist()
        rprint(f"[bold red]❌ 样本ID不唯一！重复ID: {duplicated_ids}[/bold red]")
        sys.exit(1)

    # 3. 校验空值
    if df[required_cols].isnull().any().any():
        nan_rows = df[df[required_cols].isnull().any(axis=1)][index_col].tolist()
        rprint(f"[yellow]⚠️ 警告: 样本存在空值: {nan_rows}[/yellow]")

def load_samples(csv_path, required_cols=None, index_col="sample") -> Tuple[bool, Dict]:
    """
    读取 CSV，自动生成 BAM 路径。
    
    Returns:
        merge_group (bool): 只有当【所有组】的样本数都 > 1 时，才返回 True。
        samples_dict (dict): 样本信息字典。
    """
    if required_cols is None:
        required_cols = [index_col, "group"]

    file_path = Path(csv_path)
    if not file_path.exists():
        rprint(f"[bold red]❌ Error: 找不到样本表文件: {file_path}[/bold red]")
        sys.exit(1)

    try:
        # 读取并清洗
        df = pd.read_csv(file_path, dtype=str, comment='#')
        df.columns = df.columns.str.strip()
        df = df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

        _validate_df(df, required_cols, index_col)

        # =========================================================
        # 【核心修改】 检查逻辑：必须所有组样本数 > 1
        # =========================================================
        group_counts = df['group'].value_counts()
        
        # .all() : 只有当 Series 中所有值都为 True 时，结果才为 True
        merge_group = (group_counts > 1).all()
        
        # 详细的 Debug 信息
        if merge_group:
            rprint(f"[bold green]✅ 合并条件满足:[/bold green] 所有组均包含生物学重复 (All groups have >1 samples).")
            rprint(f"   Merge Mode -> [bold green]ON[/bold green]")
        else:
            # 找出哪些组是单样本，导致了 False
            single_sample_groups = group_counts[group_counts <= 1].index.tolist()
            rprint(f"[bold yellow]⚠️ 合并条件未满足:[/bold yellow] 存在单样本组 (Singletons detected).")
            rprint(f"   导致无法完全合并的组: [bold red]{single_sample_groups}[/bold red]")
            rprint(f"   Merge Mode -> [bold red]OFF[/bold red]")

        # 自动构建 BAM 路径
        df['bam'] = df[index_col].apply(
            lambda x: f"02.mapping/STAR/sort_index/{x}.sort.bam"
        )

        samples_dict = df.set_index(index_col, drop=False).to_dict(orient="index")
        return merge_group, samples_dict

    except Exception as e:
        rprint(f"[bold red]❌ Error: load_samples 解析失败: {e}[/bold red]")
        sys.exit(1)


def parse_groups(samples_dict: Dict) -> Dict[str, List[str]]:
    """
    将 SampleID -> Info 的字典反转为 Group -> [SampleID_1, SampleID_2] 的字典。
    
    Args:
        samples_dict: load_samples 返回的字典
        
    Returns:
        dict: {'Control': ['s1', 's2'], 'Treat': ['s3', 's4']}
    """
    # 使用 defaultdict(list) 可以省去 "if key not in dict" 的判断，代码更简洁
    groups = defaultdict(list)
    
    for sample_id, info in samples_dict.items():
        # 获取该样本的组名
        group_name = info.get('group')
        
        if group_name:
            groups[group_name].append(sample_id)
        else:
            # 防御性编程：万一没有 group 字段 (虽然 load_samples 校验过)
            rprint(f"[red]⚠️ Warning: Sample {sample_id} has no group info![/red]")
            
    return dict(groups) # 转回普通字典返回


if __name__ == "__main__":
    import tempfile
    
    # 场景 1: 完美情况 (所有组都有重复) -> 期望 True
    csv_perfect = """sample, sample_name, group
    s1, A_rep1, GroupA
    s2, A_rep2, GroupA
    s3, B_rep1, GroupB
    s4, B_rep2, GroupB
    """
    
    # 场景 2: 混合情况 (GroupA 有重复，GroupB 只有一个) -> 期望 False
    csv_mixed = """sample, sample_name, group
    s1, A_rep1, GroupA
    s2, A_rep2, GroupA
    s3, B_rep1, GroupB
    """

    print("\n--- 测试场景 1: 所有组都有重复 (Expect: True) ---")
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.csv', delete=False) as tmp:
        tmp.write(csv_perfect)
        path1 = tmp.name
    load_samples(path1)
    os.remove(path1)

    print("\n--- 测试场景 2: 混合情况 (Expect: False) ---")
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.csv', delete=False) as tmp:
        tmp.write(csv_mixed)
        path2 = tmp.name
    load_samples(path2)
    os.remove(path2)