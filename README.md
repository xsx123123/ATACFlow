# 🧬 ATACFlow: 完整的ATAC-seq数据分析流程

**版本**: v0.0.5

ATACFlow是一个全面的ATAC-seq数据分析流程，涵盖了从原始数据质控到最终报告生成的完整分析过程。该流程针对植物基因组特性进行了优化，能够有效处理细胞器污染等问题。

## 📜 版本更新日志

### v0.0.5 (2024-03)
- **新增**: 灵活的 Peak Calling 策略配置 (`peak_calling.use_pooled_peaks`)
- **优化**: 重组了 peak calling 模块的输出目录结构，清晰分离单样本、pooled 和 IDR 分析
- **优化**: 将 FRiP (Fraction of Reads in Peaks) 质控集成到流程中
- **修复**: 修复了单样本组情况下 DEG 分析无法运行的 bug
- **改进**: 添加了自动回退逻辑，当无法使用 pooled peaks 时自动使用单样本共识peaks

## 🌟 核心特色 (Core Highlights)

ATACFlow 不仅仅是一个基础的比对和 Peak Calling 流程，它集成了多项前沿特性：

*   **高级足迹分析 (TOBIAS Footprinting)**: 集成了完整的 TOBIAS 流程（ATACorrect -> EstimateFootprints -> BINDetect），能够以单碱基分辨率推断转录因子的动态结合。
*   **AI 驱动的自动化报告**: 利用大语言模型 (LLM) 驱动的 AI 引擎，结合容器化技术 (Apptainer/Docker)，自动生成包含生物学解读的交互式分析报告。
*   **植物学深度优化**: 针对植物基因组中高比例的线粒体/叶绿体污染，实现了动态剔除和结构性过滤算法，最大化保留有效读段。
*   **灵活的 Peak Calling 策略** *(v0.0.5)*: 默认采用 "Group Merge -> Call Peak" 策略，显著提升了生物学重复样本中弱信号 Peak 的检测灵敏度。支持配置化选择使用 pooled 或单样本共识peaks。

## 📋 流程概览

ATACFlow包含以下主要分析阶段：

1. **数据预处理与质控**
2. **序列比对与过滤**
3. **Peak识别与注释**
4. **合并样本分析**
5. **转录因子结合位点分析**
6. **质量控制与报告生成**

```mermaid
%% 初始化配置：使用 base 主题，强制让背景透明，线条用中性色 %%
%%{
  init: {
    'theme': 'base',
    'themeVariables': {
      'primaryColor': '#E3F2FD',
      'primaryTextColor': '#2c3e50',
      'lineColor': '#7f8c8d',
      'fontFamily': 'Helvetica'
    }
  }
}%%

graph LR

    %% -------------------- 样式库 -------------------- %%
    %% 核心技巧：颜色不要太深，也不要太亮，保持中间调 %%

    %% 普通节点：圆角 + 中性边框 %%
    classDef base fill:#fff,stroke:#7f8c8d,stroke-width:1px,rx:5,ry:5,color:#333;

    %% 1. 数据流 (蓝色系) - 在黑夜模式下会显得很亮眼 %%
    classDef raw fill:#e3f2fd,stroke:#1565c0,stroke-width:2px,rx:5,ry:5,color:#0d47a1;

    %% 2. 比对流 (绿色系) %%
    classDef map fill:#e8f5e9,stroke:#2e7d32,stroke-width:2px,rx:5,ry:5,color:#1b5e20;

    %% 3. 分析流 (橙色系) %%
    classDef core fill:#fff3e0,stroke:#ef6c00,stroke-width:2px,rx:5,ry:5,color:#e65100;

    %% 4. 高级流 (紫色系) %%
    classDef adv fill:#f3e5f5,stroke:#7b1fa2,stroke-width:2px,rx:5,ry:5,color:#4a148c;

    %% 5. 终点 (深灰底白字) - 注意：fill 用深灰而不是纯黑，stroke 用浅灰，这样在黑底也能看清边界 %%
    classDef endNode fill:#37474f,stroke:#cfd8dc,stroke-width:2px,rx:15,ry:15,color:#fff;

    %% 决策节点 %%
    classDef decision fill:#fffde7,stroke:#f57f17,stroke-width:2px,rx:5,ry:5,color:#ff6f00;


    %% -------------------- 流程图内容 -------------------- %%

    %% 1. 数据清洗 %%
    subgraph S1 ["Step 1: Data Cleaning"]
        direction TB
        Raw[Raw Data]:::raw --> MD5{MD5 Check}:::base
        MD5 --> QC1[FastQC & Screen]:::base
        QC1 --> Trim[Fastp Trimming]:::base
        Trim --> Clean[Clean Data]:::raw
    end

    %% 2. 比对与过滤 %%
    subgraph S2 ["Step 2: Mapping & Filtering"]
        direction TB
        Clean --> Bowtie2[Bowtie2 Mapping]:::map
        Bowtie2 --> Filter[Structural Filtering]:::base
        Filter --> Shift[Fragment Shifting]:::base
        Shift --> BAM[Processed BAM]:::map
    end

    %% 3. Peak识别 - 三路径 (v0.0.5更新) %%
    subgraph S3 ["Step 3: Peak Calling"]
        direction TB
        BAM --> SinglePeaks[Single Sample MACS2]:::core
        SinglePeaks --> IDR[IDR Analysis (if >=2 reps)]:::core
        SinglePeaks --> SingleConsensus[Single Consensus]:::core
        
        BAM --> PooledPeaks{Pooled?}:::decision
        PooledPeaks -->|Yes| MergeBAM[Merge Group BAMs]:::core
        MergeBAM --> MergeMACS2[Merged Sample MACS2]:::core
        MergeMACS2 --> PooledConsensus[Pooled Consensus]:::core
        PooledPeaks -->|No| SingleConsensus
        
        SingleConsensus --> DEG[DEG Analysis]:::adv
        PooledConsensus --> DEG
    end

    %% 4. 高级分析 %%
    subgraph S4 ["Step 4: Advanced Analysis"]
        direction TB
        BAM -.-> ATACv[ATACv QC]:::adv
        MergeMACS2 -.-> TOBIAS[TOBIAS Motifs]:::adv
        Consensus -.-> DEG_Merge["DESeq2 Merged"]:::adv
        IndividualPeaks -.-> DEG_Single["DESeq2 Single"]:::adv
        DEG_Merge --> Enrich_Merge["GO/KEGG Enrichment Merged"]:::adv
        DEG_Single --> Enrich_Single["GO/KEGG Enrichment Single"]:::adv
    end

    %% 5. 交付 %%
    Report(Final Report):::endNode

    %% -------------------- 连线逻辑 -------------------- %%
    ATACv --> Report
    TOBIAS --> Report
    Enrich_Merge --> Report
    Enrich_Single --> Report

    %% -------------------- 关键美化：透明化 Subgraph -------------------- %%
    %% 这一步把那块黑色的背景去掉了！ %%
    %% fill:none = 透明 %%
    %% stroke:#7f8c8d = 中性灰边框 (黑白背景都可见) %%
    %% stroke-dasharray = 虚线，看起来更轻盈 %%

    style S1 fill:none,stroke:#7f8c8d,stroke-width:2px,stroke-dasharray: 5 5,color:#7f8c8d
    style S2 fill:none,stroke:#7f8c8d,stroke-width:2px,stroke-dasharray: 5 5,color:#7f8c8d
    style S3 fill:none,stroke:#7f8c8d,stroke-width:2px,stroke-dasharray: 5 5,color:#7f8c8d
    style S4 fill:none,stroke:#7f8c8d,stroke-width:2px,stroke-dasharray: 5 5,color:#7f8c8d

    %% 统一线条颜色为中性灰 %%
    linkStyle default stroke:#7f8c8d,stroke-width:1px,fill:none;
```

## 🔬 详细分析流程

### 1. 数据预处理与质控 (Quality Control & Preprocessing)

#### 目标
对原始测序数据进行质量控制和预处理，确保后续分析的数据质量。

#### 工具
- FastQC: 序列质量评估
- Fastp: 序列修剪和过滤
- MultiQC: 质控结果汇总
- FastQ Screen: 污染检测

#### 处理步骤
- **MD5校验**: 验证原始数据完整性
- **FastQC分析**: 评估原始数据质量
- **污染检测**: 使用FastQ Screen检测外源序列污染
- **序列修剪**: 使用Fastp去除低质量序列和接头
- **质控汇总**: 使用MultiQC生成综合质控报告

### 2. 序列比对与过滤 (Mapping & Filtering)

#### 目标
将高质量的reads比对到参考基因组，并进行严格的过滤以获得高质量的比对结果。

#### 工具
- Bowtie2: 序列比对
- Samtools: BAM文件处理
- 自研Rust工具: 结构性过滤

#### 处理步骤
- **序列比对**: 使用Bowtie2将reads比对到参考基因组
- **重复标记**: 使用GATK标记PCR重复
- **基础过滤**: 去除未比对、次要比对和低质量比对
- **植物学优化**: 动态剔除线粒体和叶绿体比对reads
- **结构性过滤**: 使用自研工具进行配对关系严格质控
- **坐标排序**: 生成最终的坐标排序BAM文件

### 3. Peak识别与注释 (Peak Calling & Annotation)

#### 目标
识别开放染色质区域并进行功能注释。

#### 工具
- MACS2: Peak识别
- HOMER: Peak注释

#### 处理步骤
- **Peak识别**: 使用MACS2识别开放染色质区域
- **Peak注释**: 使用HOMER对peaks进行基因功能注释
- **TSS富集分析**: 评估ATAC-seq数据质量

### 4. 合并样本分析 (Merged Sample Analysis)

#### 目标
对同一组的生物学重复样本进行合并分析，提高检测效力。

#### 处理步骤
- **BAM文件合并**: 合并同一组的生物学重复样本
- **合并Peak识别**: 在合并样本上进行peak识别
- **Peak注释**: 对合并样本的peaks进行注释

### 5. 转录因子结合位点分析 (Motif Analysis)

#### 目标
识别开放染色质区域中的转录因子结合位点，并比较不同实验组间的差异。

#### 工具
- TOBIAS (Transcription factor Occupancy prediction By Investigation of ATAC-seq Signal)

#### 分析流程
- **足迹分析 (Footprinting)**: 识别蛋白质-DNA结合足迹，推断转录因子结合位点
- **结合偏好分析 (BINDetect)**: 基于motif数据库检测转录因子结合偏好
- **组间差异分析**: 比较预定义对比组间的转录因子结合差异
- **可视化**: 生成火山图、MA图等可视化结果

#### 输出文件
- 格式化的peak文件
- 校正后的信号bigwig文件
- 足迹分析结果
- 转录因子结合检测结果和可视化图表
- 组间差异分析报告

### 6. 质量控制与报告生成 (QC & Reporting)

#### 目标
生成全面的质量控制报告和分析结果报告。

#### 工具
- ataqv: ATAC-seq特异性质控
- MultiQC: 结果汇总
- 自定义报告工具

#### 处理步骤
- **ataqv质控**: 生成ATAC-seq特异性质量控制报告
- **结果汇总**: 使用MultiQC汇总所有分析结果
- **数据交付**: 组织和打包所有分析结果
- **报告生成**: 生成最终分析报告

## 📁 输出目录结构

```
├── 00.raw_data/           # 原始数据链接和MD5校验
├── 01.qc/                 # 质控结果
│   ├── short_read_qc_r1/  # R1 reads质控
│   ├── short_read_qc_r2/  # R2 reads质控
│   ├── short_read_trim/   # 修剪后数据和质控报告
│   └── ...
├── 02.mapping/            # 比对结果
│   ├── Bowtie2/           # Bowtie2比对结果
│   ├── filter_pe/         # 过滤后结果
│   ├── shifted/           # 位移后结果
│   ├── merged/            # 合并样本结果 (run_pooled=True时)
│   └── ataqv/             # ATAC-seq质控结果
├── 03.peak_calling/       # Peak识别结果
│   ├── single/            # 单样本peak calling (始终运行)
│   ├── single_HOMER/      # 单样本peak注释
│   ├── pooled/            # 组级别pooled peak calling (run_pooled=True时)
│   ├── pooled_HOMER/      # 组级别peak注释
│   ├── idr/               # IDR分析 (组内>=2样本时)
│   └── idr_HOMER/         # IDR peak注释
├── 04.consensus/          # 共识peak集
│   ├── single/            # 单样本共识 (始终运行)
│   │   ├── all_samples_consensus_peaks.bed
│   │   ├── consensus_peaks_annotation.txt
│   │   ├── consensus_counts_matrix.txt
│   │   └── consensus_counts_matrix_ann.txt
│   └── pooled/            # 组级别共识 (run_pooled=True时)
│       ├── all_groups_consensus_peaks.bed
│       ├── consensus_peaks_annotation.txt
│       ├── consensus_counts_matrix.txt
│       └── consensus_counts_matrix_ann.txt
├── 05.ATAC_QC/            # ATAC-seq特异性质控报告
├── 06.deg_enrich/         # 差异分析和富集分析
│   ├── DEG/               # 差异peak分析结果
│   └── enrich/            # GO/KEGG富集分析
├── 06.motif_analysis/     # 转录因子结合位点分析
│   ├── 01.formatted_peaks/
│   ├── 02.signal_corrected/
│   ├── 03.footprints/
│   ├── 04.bindetect/
│   ├── 05.differential_motifs/
│   └── 06.final_report/
└── report/                # 最终报告
```

## ⚙️ Peak Calling 策略配置 (v0.0.5 新增)

ATACFlow v0.0.5 引入了灵活的 Peak Calling 策略配置，用户可以通过配置文件控制使用哪种共识peaks进行下游差异分析。

### 配置选项

在 `config.yaml` 中添加以下配置：

```yaml
# Peak Calling Configuration
# 控制peaks的调用和合并方式
peak_calling:
  # 使用pooled (合并组) 共识peaks进行DEG分析
  # true: 在组内合并生物学重复后调用peaks (更robust，推荐)
  # false: 使用单样本共识peaks
  use_pooled_peaks: true
```

### 运行逻辑

| 配置 | merge_group | run_pooled | DEG分析输入 |
|------|-------------|------------|-------------|
| use_pooled_peaks: true | True | ✅ 运行 | `04.consensus/pooled/consensus_counts_matrix_ann.txt` |
| use_pooled_peaks: true | False (有单样本组) | ❌ 跳过 | `04.consensus/single/consensus_counts_matrix_ann.txt` (自动回退) |
| use_pooled_peaks: false | 任意 | ❌ 跳过 | `04.consensus/single/consensus_counts_matrix_ann.txt` |

> **注意**: `merge_group` 是自动检测的变量，仅当所有实验组都有多于1个生物学重复时为 True。如果存在单样本组，系统会自动回退到单样本共识peaks。

### 三种分析模式说明

1. **单样本分析 (Single Sample)**: 对每个样本独立进行 peak calling，然后合并所有样本的peaks生成共识集

2. **Pooled分析 (合并样本)**: 先将同一组的生物学重复样本的BAM文件合并，然后在合并后的BAM上进行peak calling
   - 优势: 增加测序深度，提高弱信号peaks的检测灵敏度
   - 适用于: 有多个生物学重复的实验组

3. **IDR分析**: 在组内样本间进行 Irreproducible Discovery Rate 分析，识别可重复的peaks
   - 用途: 质量控制，评估生物学重复间的一致性
   - 注意: IDR过滤后的peaks数量可能较少，不建议用于差异分析

### 推荐的配置

```yaml
# 推荐配置: 使用pooled peaks
# 适用于: 有2个或以上生物学重复的实验
peak_calling:
  use_pooled_peaks: true

# 替代配置: 使用单样本共识peaks
# 适用于: 想保留更多peaks，或生物学重复较少的情况
peak_calling:
  use_pooled_peaks: false
```

## ⚙️ 配置文件

流程使用多个配置文件来控制分析参数：

- `config.yaml`: 基本配置参数
- `reference.yaml`: 参考基因组路径配置
- `run_parameter.yaml`: 软件运行参数
- `cluster_config.yaml`: 集群资源配置

## 🧠 流程设计哲学

ATACFlow遵循以下设计理念：

1. **模块化设计**: 每个分析步骤都是独立的模块，便于维护和扩展
2. **质量优先**: 严格的质控和过滤步骤确保分析结果的可靠性
3. **植物学优化**: 针对植物基因组特点进行特殊优化处理
4. **自动化交付**: 自动生成完整的分析报告和交付清单
5. **可重现性**: 使用Snakemake确保分析流程的可重现性

## 🗺️ 路线图与未来改进 (Roadmap)

为了进一步提升流程的稳健性和科学性，计划在后续版本中引入以下改进：

### 已有计划
1.  **IDR (Irreproducible Discovery Rate) 支持**: 引入 ENCODE 推荐的 IDR 框架，定量评估生物学重复之间的一致性，为 Peak 过滤提供更严谨的统计学依据。
2.  **自动 Blacklist 过滤**: 针对常见物种（如人、小鼠、拟南芥等）集成 ENCODE Blacklist 自动过滤步骤，剔除已知的信号伪影区域。
3.  **动态 QC 阈值预警**: 在 MultiQC 报告中集成基于 `ataqv` 指标（如 TSS Enrichment Score, Fragment Length Distribution）的自动判定系统，对低质量样本进行实时报警。
4.  **调控网络增强**: 深化 TOBIAS 结果与 DEG 结果的关联分析，构建更精细的 TF-Target 基因调控网络。

### 与 nf-core/atacseq 对比后的新增优化方向

#### 🔴 高优先级
1.  **差异可及性分析模块**: 集成 DESeq2/edgeR 进行组间差异可及性分析，包括 PCA 聚类、差异 peak 鉴定和可视化。
2.  **容器化支持**: 增加 Docker/Singularity 容器支持，与 Conda 环境并存，提高流程的可重复性和跨平台部署能力。
3.  **精细化 BAM 过滤策略**: 参考 nf-core 的详尽过滤策略，增加：
    - 错配数限制（>4 mismatches）过滤
    - soft-clipped reads 过滤
    - 插入大小限制（>2kb）过滤
    - 染色体间配对 reads 过滤
    - FR 方向验证

#### 🟡 中优先级
4.  **IGV 会话文件自动生成**: 自动创建包含 bigWig  tracks、peaks 和差异位点的 IGV 会话文件，方便用户直接可视化。
5.  **测试数据集与 CI/CD**: 添加完整的测试数据集和自动化测试流程，确保流程更新的稳健性。
6.  **扩展比对工具选择**: 增加 BWA 和 STAR 作为可选比对工具，提供更多选择。

#### 🟢 低优先级
7.  **模块化重构**: 将常用分析步骤抽象为可复用模块，便于跨项目共享和维护。
8.  **多工作流系统支持**: 评估 Nextflow 版本的可行性，与现有 Snakemake 版本并存。

## 🚀 使用方法

```bash
# 运行完整分析流程
snakemake --cores=60 -p --conda-frontend mamba --use-conda --rerun-triggers mtime  --logger rich-loguru --config analysisyaml=/data/jzhang/project/Temp/atac_human_PRJNA427322/01.workflow/config.yaml
```