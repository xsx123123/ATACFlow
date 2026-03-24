# ATACFlow MCP Server

这是一个基于 [Model Context Protocol (MCP)](https://modelcontextprotocol.io/) 构建的服务端，使用 **[uv](https://docs.astral.sh/uv/)** 进行现代化的环境管理。

## 🌟 功能特性

- **基因组查询**：列出系统支持的参考基因组。
- **配置生成**：自动化生成 `config.yaml`, `samples.csv`, `contrasts.csv`。
- **异步运行**：后台启动 Snakemake，不阻塞 AI。
- **环境隔离**：使用 `uv` 确保依赖库与分析环境互不干扰。
- **双模式支持**：本地 stdio 模式 + 远程部署能力。

## 🛠️ 快速上手 (使用 uv)

### 1. 安装依赖
如果你还没有安装 `uv`，请先安装：
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

在 `mcp/` 目录下同步环境：
```bash
cd mcp
uv sync
```

### 2. 测试运行（使用便捷脚本）
```bash
# 方式一：快速测试（推荐）
./start.sh test

# 方式二：前台运行本地模式
./start.sh local

# 方式三：后台运行
./start.sh background
```

或者手动使用官方的 MCP Inspector（需要安装 Node.js）：
```bash
npx @modelcontextprotocol/inspector uv --directory mcp run server.py
```

## 🤖 在 AI 客户端中使用

### 本地使用（同一台机器）
编辑 `claude_desktop_config.json`：
```json
{
  "mcpServers": {
    "atacflow": {
      "command": "uv",
      "args": [
        "--directory",
        "/home/zj/pipeline/ATACFlow/mcp",
        "run",
        "server.py"
      ]
    }
  }
}
```

### 远程部署（通过 SSH 隧道）
在你的本地机器上编辑 `claude_desktop_config.json`：
```json
{
  "mcpServers": {
    "atacflow": {
      "command": "ssh",
      "args": [
        "zj@your-server-ip",
        "cd", "/home/zj/pipeline/ATACFlow/mcp", "&&",
        "uv", "run", "python", "server.py"
      ]
    }
  }
}
```

## 📂 目录结构
- `server.py`: 服务器核心代码。
- `start.sh`: 便捷启动脚本。
- `mcp_config.yaml`: 服务配置文件。
- `pyproject.toml`: `uv` 项目配置文件。
- `TEST_AND_DEPLOY.md`: 详细的测试与部署指南。

## ⚠️ 注意事项
- **项目命名**：本项目在 `pyproject.toml` 中命名为 `atacflow-mcp`。
- **环境要求**：运行 Snakemake 仍需系统路径中有 `conda` 或 `mamba`。
- **远程部署**：详见 `TEST_AND_DEPLOY.md` 获取更多远程部署方案。
