# Sanger16S Web 部署说明

Sanger16S Web 是一个适合局域网部署的批量分析界面，支持：

- 上传 FASTA，直接做 16S BLAST 鉴定
- 上传 FASTQ/FASTQ.GZ，先跑共识流程，再做 16S BLAST 鉴定

## 运行前提

建议在 Linux 服务器上运行，并先进入正确的 Conda 环境：

```bash
conda activate sanger16s
```

必须满足以下条件：

- `blastn` 已在 `PATH`
- 本地 16S BLAST 数据库已准备完成
- `Flask` 已安装
- `bash` 可用
- `run_16s_consensus_pipeline.sh` 存在
- `scripts/` 中的 helper 脚本完整存在

## 启动命令

注意：`--db-dir` 指向“本地 BLAST 数据库文件所在目录”，例如目录内应包含
`16S_ribosomal_RNA.nsq`、`16S_ribosomal_RNA.nin`、`16S_ribosomal_RNA.nhr` 等文件，
并不是固定写成 `db` 这个名字。

```bash
python sanger16s.py serve \
  --host 0.0.0.0 \
  --port 8080 \
  --db-dir /data/blast/16s-db \
  --db-version 2026-06-05 \
  --consensus-pipeline-script ./run_16s_consensus_pipeline.sh \
  --work-dir ./web-work \
  --threads 8 \
  --queue-workers 1 \
  --max-target-seqs 1 \
  --max-upload-mb 64
```

浏览器访问：

```text
http://<server-ip>:8080
```

## 当前能力

- 支持多文件上传
- 支持 FASTA 模式和 FASTQ 模式
- FASTQ 模式自动先生成共识序列
- 后台任务排队执行
- 页面轮询显示任务状态
- 页面内预览结果前 20 行
- 任务完成后下载 CSV

## 工作目录

`--work-dir` 下会生成：

```text
web-work/
  jobs/
    job-xxxx/
      input/
      pipeline-output/
      result.csv
      status.json
```

说明：

- `input/`：上传的原始文件
- `pipeline-output/`：FASTQ 模式下的共识流程输出
- `result.csv`：最终 16S 鉴定结果
- `status.json`：前端任务状态数据

## 建议的生产部署方式

- 用 `systemd` 常驻运行 Python 服务
- 用 `nginx` 做反向代理
- 限制内网访问范围
- 定时清理旧任务目录

仓库附带示例文件：

- `deploy/systemd/sanger16s-web.service.example`
- `deploy/nginx/sanger16s-web.conf.example`
- `deploy/sanger16s-web.env.example`

## 参数建议

- `--queue-workers 1`
  - 避免多个任务同时占满 CPU 和 I/O
- `--threads`
  - 控制单个 BLAST 任务的线程数
- `--max-target-seqs 1`
  - 默认只保留最佳命中
- `--db-version`
  - 建议显式填写数据库发布日期或内部版本号
- `--consensus-pipeline-script`
  - 指向 `run_16s_consensus_pipeline.sh`

## 后续可增强项

- 登录或反向代理鉴权
- 自动清理历史任务
- 导出 Excel
- 增加任务详情页
- 增加 pipeline 中间结果下载
