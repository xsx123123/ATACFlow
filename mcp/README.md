# ATACFlow MCP Server

这是一个基于 [Model Context Protocol (MCP)](https://modelcontextprotocol.io/) 构建的服务端，使用 **[uv](https://docs.astral.sh/uv/)** 进行现代化的环境管理。

## 🌟 功能特性

- **基因组查询**：列出系统支持的参考基因组。
- **配置生成**：自动化生成 `config.yaml`, `samples.csv`, `contrasts.csv`。
- **异步运行**：后台启动 Snakemake，不阻塞 AI。
- **环境隔离**：使用 `uv` 确保依赖库与分析环境互不干扰。

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

### 2. 测试运行
推荐使用官方的 MCP Inspector（需要安装 Node.js）来手动测试工具：
```bash
# 在 ATACFlow 根目录下
npx @modelcontextprotocol/inspector uv --directory mcp run server.py
```

## 🤖 在 AI 客户端中使用

### 配置 Claude Desktop
编辑 `claude_desktop_config.json`，添加以下配置。`uv` 会自动处理路径和虚拟环境。

```json
{
  "mcpServers": {
    "atacflow": {
      "command": "uv",
      "args": [
        "--directory",
        "/YOUR_PATH/ATACFlow/mcp",
        "run",
        "server.py"
      ],
      "env": {
        "PATH": "/usr/local/bin:/usr/bin:/bin:/YOUR_CONDA_PATH/bin"
      }
    }
  }
}
```

## 📂 目录结构
- `server.py`: 服务器核心代码。
- `pyproject.toml`: `uv` 项目配置文件（包含依赖定义）。
- `uv.lock`: 依赖锁定文件，确保环境可复现。

## ⚠️ 注意事项
- **项目命名**：本项目在 `pyproject.toml` 中命名为 `atacflow-mcp` 以避免与官方 `mcp` 库冲突。
- **环境要求**：运行 Snakemake 仍需系统路径中有 `conda` 或 `mamba`。
