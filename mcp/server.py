#!/usr/bin/env python3
"""
ATACFlow MCP Server
A Model Context Protocol server for ATAC-seq analysis using ATACFlow
"""

import os
import sys
import yaml
import csv       # 优化：移至顶部全局导入
import asyncio   # 优化：引入异步处理
import subprocess
from pathlib import Path
from typing import Optional, List, Dict, Any
from fastmcp import FastMCP
from pydantic import BaseModel, Field

# Add the parent directory to path to access ATACFlow
sys.path.insert(0, str(Path(__file__).parent.parent))

# Initialize MCP server
mcp = FastMCP("ATACFlow")

# Configuration
ATACFLOW_ROOT = Path(__file__).parent.parent
SKILLS_DIR = ATACFLOW_ROOT / "skills"
CONFIG_DIR = ATACFLOW_ROOT / "config"
EXAMPLES_DIR = SKILLS_DIR / "examples"
MCP_CONFIG_FILE = Path(__file__).parent / "mcp_config.yaml"

# Load local MCP configuration for tool paths
def load_mcp_config():
    config = {
        "conda_path": "conda",
        "snakemake_path": "snakemake",
        "default_env": "atacflow"
    }
    if MCP_CONFIG_FILE.exists():
        try:
            with open(MCP_CONFIG_FILE, "r") as f:
                user_config = yaml.safe_load(f)
                if user_config:
                    config.update(user_config)
        except Exception:
            pass
    return config

MCP_PATHS = load_mcp_config()

class ProjectConfig(BaseModel):
    """Project configuration model"""
    project_name: str = Field(..., description="Name of the project")
    genome_version: str = Field(..., description="Genome version (e.g., hg38, TAIR10.1)")
    species: str = Field(..., description="Species name (e.g., Homo_sapiens)")
    raw_data_path: List[str] = Field(..., description="Paths to raw data directories")
    workflow_dir: str = Field(..., description="Working directory for workflow")
    output_dir: str = Field(..., description="Output directory for results")
    mapping_tool: str = Field(default="chromap", description="Mapping tool (chromap or bowtie2)")
    use_pooled_peaks: bool = Field(default=True, description="Use pooled peaks for DEG analysis")
    only_qc: bool = Field(default=False, description="Only run QC analysis")


@mcp.tool()
def list_supported_genomes() -> List[Dict[str, str]]:
    """List all supported genome versions in ATACFlow"""
    reference_yaml = CONFIG_DIR / "reference.yaml"
    if not reference_yaml.exists():
        return [{"error": "reference.yaml not found"}]

    with open(reference_yaml, "r", encoding="utf-8") as f:
        ref_config = yaml.safe_load(f)

    genomes = []
    if "Bowtie2_index" in ref_config:
        for genome_name, genome_info in ref_config["Bowtie2_index"].items():
            genomes.append(
                {"name": genome_name, "description": f"Reference genome: {genome_name}"}
            )

    common_genomes = [
        {"name": "hg38", "description": "Human reference genome (GRCh38)"},
        {"name": "GRCm39", "description": "Mouse reference genome (GRCm39)"},
        {"name": "TAIR10.1", "description": "Arabidopsis reference genome (TAIR10.1)"},
        {"name": "Lsat_Salinas_v11", "description": "Lettuce reference genome (v11)"},
        {"name": "Lsat_Salinas_v8", "description": "Lettuce reference genome (v8)"},
        {"name": "ITAG4.1", "description": "Tomato reference genome (ITAG4.1)"},
    ]

    return common_genomes


@mcp.tool()
def get_config_template(template_type: str = "standard") -> str:
    """Get an ATACFlow configuration template"""
    template_files = {
        "complete": EXAMPLES_DIR / "config_complete.yaml",
        "standard": EXAMPLES_DIR / "config_standard.yaml",
        "qc_only": EXAMPLES_DIR / "config_qc_only.yaml",
    }
    template_file = template_files.get(template_type, template_files["standard"])

    if not template_file.exists():
        return f"Error: Template {template_type} not found"

    with open(template_file, "r", encoding="utf-8") as f:
        return f.read()


@mcp.tool()
def generate_config_file(config: ProjectConfig, output_path: str) -> str:
    """
    Generate an ATACFlow config.yaml file using the structured ProjectConfig model
    """
    try:
        out_path = Path(output_path).resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        
        # 优化：利用预定义的 Pydantic 模型生成标准的 YAML 配置
        config_dict = config.model_dump()
        
        with open(out_path, "w", encoding="utf-8") as f:
            yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)
        
        return f"Successfully generated config file at {out_path}"
    except Exception as e:
        return f"Error generating config file: {str(e)}"


@mcp.tool()
def create_sample_csv(sample_data: List[Dict[str, str]], output_path: str) -> str:
    """Create a samples.csv file for ATACFlow"""
    try:
        out_path = Path(output_path).resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)

        with open(out_path, "w", newline="", encoding="utf-8") as csvfile:
            fieldnames = ["sample", "sample_name", "group"]
            # 优化：添加 extrasaction='ignore' 防止 LLM 传入多余字段导致崩溃
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, extrasaction='ignore')
            writer.writeheader()
            for sample in sample_data:
                # 优化：补全缺失字段防止 KeyError
                safe_sample = {key: sample.get(key, "unknown") for key in fieldnames}
                writer.writerow(safe_sample)
        return f"Successfully created samples.csv at {out_path}"
    except Exception as e:
        return f"Error creating samples.csv: {str(e)}"


@mcp.tool()
def create_contrasts_csv(contrasts: List[Dict[str, str]], output_path: str) -> str:
    """Create a contrasts.csv file for differential peak analysis"""
    try:
        out_path = Path(output_path).resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)

        with open(out_path, "w", newline="", encoding="utf-8") as csvfile:
            fieldnames = ["contrast", "treatment"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, extrasaction='ignore')
            writer.writeheader()
            for contrast in contrasts:
                safe_contrast = {key: contrast.get(key, "unknown") for key in fieldnames}
                writer.writerow(safe_contrast)
        return f"Successfully created contrasts.csv at {out_path}"
    except Exception as e:
        return f"Error creating contrasts.csv: {str(e)}"


@mcp.tool()
def validate_config(config_path: str) -> Dict[str, Any]:
    """Validate an ATACFlow configuration file"""
    config_file = Path(config_path).resolve()
    if not config_file.exists():
        return {"valid": False, "errors": ["Configuration file not found"]}

    try:
        with open(config_file, "r", encoding="utf-8") as f:
            config = yaml.safe_load(f)

        errors = []
        required_fields = [
            "project_name", "Genome_Version", "species",
            "raw_data_path", "workflow", "data_deliver"
        ]

        for field in required_fields:
            if field not in config:
                errors.append(f"Missing required field: {field}")

        return {
            "valid": len(errors) == 0,
            "errors": errors,
            "config": config if len(errors) == 0 else None,
        }
    except Exception as e:
        return {"valid": False, "errors": [str(e)]}


@mcp.tool()
async def run_atacflow(config_path: str, cores: int = 20, dry_run: bool = False) -> str:
    """
    Run the ATACFlow pipeline (Asynchronous / Detached)
    """
    config_file = Path(config_path).resolve()
    if not config_file.exists():
        return f"Error: Configuration file not found at {config_file}"

    # 优化：更优雅的命令拼接方式
    snakemake_bin = MCP_PATHS.get("snakemake_path", "snakemake")
    cmd = [snakemake_bin]
    
    if dry_run:
        cmd.extend(["-n", "--quiet"])
        
    cmd.extend([
        f"--cores={cores}",
        "-p",
        "--conda-frontend", "mamba",
        "--use-conda",
        "--rerun-triggers", "mtime",
        "--logger", "rich-loguru",
        "--config", f"analysisyaml={config_file}"
    ])

    try:
        if dry_run:
            # 优化：dry_run 使用 asyncio 快速返回预览结果
            process = await asyncio.create_subprocess_exec(
                *cmd,
                cwd=str(ATACFLOW_ROOT),
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE
            )
            stdout, stderr = await process.communicate()
            
            output = f"Dry Run Command: {' '.join(cmd)}\n\n"
            if stdout:
                output += f"Stdout:\n{stdout.decode()}\n"
            if stderr:
                output += f"Stderr:\n{stderr.decode()}\n"
            return output
        else:
            # 优化：实际运行使用后台游离进程 (Detached Process)，彻底释放 MCP 服务器
            log_file = config_file.parent / "atacflow_run.log"
            
            with open(log_file, "w") as f:
                subprocess.Popen(
                    cmd,
                    cwd=str(ATACFLOW_ROOT),
                    stdout=f,
                    stderr=subprocess.STDOUT,
                    start_new_session=True  # 将进程与当前服务器剥离
                )
            
            return (
                f"🚀 ATACFlow pipeline has been successfully started in the background!\n"
                f"Command: {' '.join(cmd)}\n"
                f"Execution logs are being written to: {log_file}\n"
                f"You can safely continue other tasks."
            )
            
    except FileNotFoundError:
        return (
            f"Error: '{snakemake_bin}' not found. "
            "Please update the 'snakemake_path' in mcp/mcp_config.yaml "
            "with the absolute path to your snakemake executable."
        )
    except Exception as e:
        return f"Error starting ATACFlow: {str(e)}"


@mcp.tool()
def check_conda_environment(env_name: Optional[str] = None) -> Dict[str, Any]:
    """Check if the required conda environment exists and is valid"""
    env_name = env_name or MCP_PATHS.get("default_env", "atacflow")
    conda_bin = MCP_PATHS.get("conda_path", "conda")
    
    try:
        result = subprocess.run([conda_bin, "--version"], capture_output=True, text=True)
        if result.returncode != 0:
            return {
                "available": False, 
                "error": f"Conda not found at '{conda_bin}'. Please update mcp/mcp_config.yaml."
            }

        result = subprocess.run([conda_bin, "env", "list"], capture_output=True, text=True)
        env_exists = env_name in result.stdout

        snakemake_available = False
        if env_exists:
            try:
                # 尝试在该环境下运行 snakemake
                result = subprocess.run(
                    [conda_bin, "run", "-n", env_name, "snakemake", "--version"],
                    capture_output=True, text=True, timeout=10
                )
                snakemake_available = result.returncode == 0
            except:
                pass

        return {
            "available": env_exists,
            "env_name": env_name,
            "snakemake_available": snakemake_available,
            "conda_version": result.stdout.strip().split("\n")[0] if result.stdout else "unknown",
            "message": "Success" if env_exists else f"Environment '{env_name}' not found. Please create it or update config."
        }
    except FileNotFoundError:
        return {
            "available": False, 
            "error": f"Conda executable '{conda_bin}' not found. Please specify absolute path in mcp/mcp_config.yaml"
        }
    except Exception as e:
        return {"available": False, "error": str(e)}


@mcp.tool()
def get_project_structure() -> Dict[str, Any]:
    """Get the recommended ATACFlow project structure"""
    return {
        "directories": [
            "Project_Root/00.raw_data/",
            "Project_Root/01.workflow/",
            "Project_Root/02.data_deliver/",
        ],
        "files": [
            "Project_Root/01.workflow/config.yaml",
            "Project_Root/01.workflow/samples.csv",
            "Project_Root/01.workflow/contrasts.csv",
        ],
        "description": "Recommended project structure for ATACFlow analysis",
    }

# Resources remain the same...
@mcp.resource("atacflow://config-templates/complete")
def get_complete_config_template() -> str:
    template_file = EXAMPLES_DIR / "config_complete.yaml"
    if template_file.exists():
        with open(template_file, "r", encoding="utf-8") as f:
            return f.read()
    return "Template not found"

@mcp.resource("atacflow://config-templates/standard")
def get_standard_config_template() -> str:
    template_file = EXAMPLES_DIR / "config_standard.yaml"
    if template_file.exists():
        with open(template_file, "r", encoding="utf-8") as f:
            return f.read()
    return "Template not found"

@mcp.resource("atacflow://config-templates/qc-only")
def get_qc_config_template() -> str:
    template_file = EXAMPLES_DIR / "config_qc_only.yaml"
    if template_file.exists():
        with open(template_file, "r", encoding="utf-8") as f:
            return f.read()
    return "Template not found"

@mcp.resource("atacflow://skills/main")
def get_skill_documentation() -> str:
    skill_file = SKILLS_DIR / "SKILL.md"
    if skill_file.exists():
        with open(skill_file, "r", encoding="utf-8") as f:
            return f.read()
    return "SKILL.md not found"


# --- Prompts section ---

@mcp.prompt()
def setup_new_project(project_name: str = "My_ATAC_Project") -> str:
    """Guide the AI to set up a new ATACFlow analysis project"""
    return f"""I want to start a new ATAC-seq analysis project named "{project_name}" using ATACFlow. 

Please follow these steps:
1. First, use `list_supported_genomes` to show me which reference genomes are available.
2. Check if my environment is ready using `check_conda_environment`.
3. Explain the recommended project structure using `get_project_structure`.
4. Ask me for the following details:
   - Species name
   - Absolute path to raw FASTQ data
   - Desired output directory
5. Once I provide those, use `generate_config_file` to create the configuration.
"""

@mcp.prompt()
def troubleshoot_failure(log_path: str = "01.workflow/atacflow_run.log") -> str:
    """Guide the AI to help debug a pipeline failure"""
    return f"""My ATACFlow pipeline failed. The log file is located at "{log_path}".

Please help me debug by:
1. Reading the last 50 lines of the log file to identify the specific rule that failed.
2. Checking if there are any common issues like "Command not found" or "Out of memory".
3. Using `validate_config` to ensure my config.yaml is still valid.
4. Suggesting a fix or explaining the error message in simple terms.
"""

@mcp.prompt()
def prepare_publication_report() -> str:
    """Guide the AI to summarize results for a report or publication"""
    return """I have finished my ATAC-seq analysis. Please help me summarize the results for a report.

Please:
1. Guide me to find the key QC metrics (e.g., TSS enrichment, FRiP scores).
2. Look for the peak calling results in the output directory.
3. Help me describe the differential analysis (DEG) results if they were performed.
4. Provide a structure for the "Methods" section of my paper based on the current ATACFlow workflow.
"""

if __name__ == "__main__":
    print(f"Starting ATACFlow MCP Server...")
    print(f"ATACFlow Root: {ATACFLOW_ROOT}")
    mcp.run()
