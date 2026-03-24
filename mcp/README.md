# ATACFlow MCP Server

这是一个基于 [Model Context Protocol (MCP)](https://modelcontextprotocol.io/) 构建的服务端，使用 **[uv](https://docs.astral.sh/uv/)** 进行现代化的环境管理。

## 🌟 功能特性

- **基因组查询**：列出系统支持的参考基因组。
- **配置生成**：自动化生成 `config.yaml`, `samples.csv`, `contrasts.csv`。
- **异步运行**：后台启动 Snakemake，不阻塞 AI。
- **环境隔离**：使用 `uv` 确保依赖库与分析环境互不干扰。
- **双模式支持**：本地 stdio 模式 + 远程部署能力。

---

## 🚀 快速开始

### 前置要求
- Python 3.13+
- uv (包管理器)
- conda/mamba (用于运行 Snakemake)
- Node.js (可选，用于 MCP Inspector 测试)

### 1. 安装 uv (如果未安装)
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

### 2. 安装项目依赖
```bash
cd /home/zj/pipeline/ATACFlow/mcp
uv sync
```

### 3. 测试运行

#### 方式一：使用便捷脚本（推荐）
```bash
# 快速测试（使用 MCP Inspector）
./start.sh test

# 前台运行本地模式
./start.sh local

# 后台运行
./start.sh background
```

#### 方式二：手动使用 MCP Inspector
```bash
npx @modelcontextprotocol/inspector uv --directory mcp run server.py
```

---

## 🤖 在 AI 客户端中使用

### 场景一：本地使用（服务端和客户端在同一台机器）

编辑你的 `claude_desktop_config.json`：
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

### 场景二：远程部署（服务端在远程服务器，客户端在本地）

#### 方案 A：SSH 隧道方式（最简单，无需修改代码）

在你的**本地机器**上编辑 `claude_desktop_config.json`：
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

**优点**：
- 不需要额外配置网络服务
- SSH 加密传输，安全
- 配置简单，即插即用

#### 方案 B：生产环境部署（systemd + Nginx）

如果需要作为长期运行的服务，请参考 `TEST_AND_DEPLOY.md` 中的完整方案，包括：
- systemd 服务配置（自动重启、日志管理）
- Nginx 反向代理配置
- 防火墙配置建议

---

## 📂 目录结构

| 文件 | 说明 |
|------|------|
| `server.py` | MCP 服务器核心代码 |
| `start.sh` | 便捷启动脚本 |
| `mcp_config.yaml` | 服务配置文件（路径、网络等） |
| `pyproject.toml` | uv 项目配置（依赖定义） |
| `uv.lock` | 依赖锁定文件，确保环境可复现 |
| `TEST_AND_DEPLOY.md` | 详细的测试与部署指南 |
| `README.md` | 本文档 |

---

## ⚙️ 配置说明

编辑 `mcp_config.yaml` 可以自定义：

```yaml
# 路径配置
conda_path: "micromamba"
snakemake_path: "/path/to/snakemake"
default_env: "snakemake9"

# 远程服务配置（可选）
host: "0.0.0.0"    # 监听地址
port: 8000          # 监听端口
```

---

## 🔧 故障排查

### 问题 1：`uv` 命令找不到
```bash
# 确认 uv 安装位置
which uv

# 或重新安装
curl -LsSf https://astral.sh/uv/install.sh | sh
```

### 问题 2：依赖安装失败
```bash
cd /home/zj/pipeline/ATACFlow/mcp
rm -rf .venv uv.lock
uv sync
```

### 问题 3：测试时工具不显示
- 确认 `uv sync` 成功完成
- 检查 Python 版本是否 >= 3.13
- 查看 MCP Inspector 的日志输出

---

## ⚠️ 注意事项

- **项目命名**：本项目在 `pyproject.toml` 中命名为 `atacflow-mcp`，避免与官方 `mcp` 库冲突。
- **环境要求**：运行 Snakemake 仍需系统路径中有 `conda` 或 `mamba`。
- **远程部署**：生产环境建议使用 SSH 隧道或 systemd + Nginx 方案，详见 `TEST_AND_DEPLOY.md`。

---

## 📚 更多信息

- [Model Context Protocol 官方文档](https://modelcontextprotocol.io/)
- [fastmcp 项目](https://github.com/jlowin/fastmcp)
- [uv 包管理器](https://docs.astral.sh/uv/)
