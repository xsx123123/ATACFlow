# 🧬 单样本 ATAC-seq 核心质控与预处理流程 (QC & Filtering)

本流程采用多级过滤策略，旨在从原始测序数据中提取高纯度的 Open Chromatin 信号。针对植物基因组特性（如高细胞器污染）及配对末端测序（Paired-End）的结构完整性进行了深度优化。

流程核心包含以下五个关键阶段：

## 1. 预处理与重复标记 (Pre-processing)

### 目标
标准化 BAM 文件格式并计算文库复杂度，而非盲目删除数据。

### 工具
- GATK AddOrReplaceReadGroups
- GATK MarkDuplicates (Picard)

### 关键操作
- **RG 标准化**: 为 BAM 添加 Read Group 信息，兼容 GATK 下游标准。
- **Soft-Marking**: 识别 PCR 重复，但仅修改 Flag (标记为 1024) 而不物理删除。
- **目的**: 保留完整数据以准确计算重复率 (Duplication Rate)，作为评估文库构建质量 (Library Complexity) 的核心指标。

## 2. 基础过滤 (Basic Filtering / "粗筛")

### 目标
剔除低质量比对及非基因组信号干扰。

### 工具
- samtools
- bedtools

### 过滤标准
- **去重**: 物理剔除 Flag 包含 1024 的 Reads。
- **比对质量**: 剔除未比对 (Flag 4)、次要比对 (Flag 256, Secondary) 及低质量比对 (MAPQ < 30)。
- **🌿 植物学特别优化**: 动态剔除比对到线粒体 (chrM/Mt) 和叶绿体 (Pt/Ct) 的 Reads，有效降低细胞器 DNA 对核基因组信号的干扰。
- **黑名单清洗 (可选)**: 使用 bedtools intersect -v 剔除 Blacklist 区域（高信号假阳性区域）。

## 3. 进阶精修 (Refinement / "细筛")

### 目标
基于物理特征进一步清洗数据。

### 工具
- samtools view -e (Expression Filter)

### 过滤标准
- **插入片段长度**: abs(TLEN) <= 2000bp。ATAC-seq 主要富集核小体间区及单/双核小体片段，剔除过长片段有助于降低背景噪音。
- **错配容忍度**: [NM] <= 4。剔除错配过多的 Reads，确保序列比对的高度置信。

## 4. 结构性深度过滤 (Structural PE Filtering / "核心技术")

### 目标
基于配对关系 (Pair-Awareness) 的严格质控，解决传统工具无法处理的单端残留问题。

### 前置要求
数据必须先转换为 名称排序 (Name Sorted)。

### 工具
自研高性能 Rust 工具 filter_pe (基于 rust-htslib 开发，多线程加速)。

### 核心逻辑
- **嵌合体去除 (Chimeric Removal)**: 严格剔除 R1 和 R2 比对到不同染色体的异常 Reads。
- **方向校正 (Orientation Check)**: 仅保留标准的 FR (Forward-Reverse) 配对方向，剔除 FF/RR 等倒位或串联异常。
- **孤儿去除 (Orphan Removal)**: 采用 "All-or-Nothing" 策略。如果一对 Reads 中有一条因质量问题被前序步骤丢弃，其配对 Read 将在此步被强制丢弃，严禁"单端 Reads"混入双端分析流程。
- **输出**: 生成详细的 JSON 统计报告及被丢弃 Reads 的 BAM 文件以供回溯。

## 5. 最终标准化 (Final Standardization)

### 目标
生成适配下游分析的最终文件。

### 工具
- samtools sort
- samtools index

### 操作
将 Rust 工具输出的名称排序 BAM，重新转换为 坐标排序 (Coordinate Sorted) 并建立索引 (.bai)。

### 目的
产出符合 MACS2 (Peak Calling)、deepTools (Visualization)、IGV 等软件标准的最终文件。

## 流程设计哲学

通过 "Sort-by-Name -> Structural Filter -> Sort-by-Coordinate" 的双重排序策略，我们在保证下游软件兼容性的同时，实现了对 Paired-End 数据的绝对控制，最大限度提升了数据的信噪比。