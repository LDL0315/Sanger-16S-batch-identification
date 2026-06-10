# Sanger16S

Sanger16S 是一个基于本地 NCBI BLAST 数据库的 16S 细菌鉴定工具，支持：

- 直接批量上传 FASTA 做 16S 鉴定
- 上传 FASTQ/FASTQ.GZ，先生成共识序列，再做 16S 鉴定
- 通过局域网网页批量提交任务、查看状态并下载结果

## 快速入口

- 中文 CLI 文档：[README.zh-CN.md](./README.zh-CN.md)
- 中文网页部署文档：[WEBAPP.zh-CN.md](./WEBAPP.zh-CN.md)
- 英文首页说明：[README.md](./README.md)
- 版本变更记录：[CHANGELOG.md](./CHANGELOG.md)
- 首个 release 文案：[RELEASE_NOTES_v0.5.0.md](./RELEASE_NOTES_v0.5.0.md)
- 贡献说明：[CONTRIBUTING.md](./CONTRIBUTING.md)
- 部署示例说明：[deploy/README.md](./deploy/README.md)

## 仓库结构

```text
.
├── sanger16s.py
├── run_16s_consensus_pipeline.sh
├── scripts/
├── web/
├── deploy/
├── environment.yml
├── LICENSE
├── CHANGELOG.md
├── README.md
├── README.zh-CN.md
└── WEBAPP.zh-CN.md
```

## 环境安装

推荐在 Linux 上使用 Conda 或 Mamba：

```bash
conda env create -f environment.yml
conda activate sanger16s
```

`environment.yml` 已包含：

- `python>=3.10`
- `blast`
- `flask`
- `seqkit`
- `vsearch`
- `minimap2`
- `samtools`

FASTQ 共识流程还要求系统可用：

- `bash`
- `gzip`

## 命令行用法

检查 BLAST 和数据库：

```bash
python sanger16s.py check --db-dir /data/blast/16s-db
```

下载本地 16S BLAST 数据库：

```bash
python sanger16s.py setup-db --db-dir /data/blast/16s-db
```

这里的 `--db-dir` 指的是“本地 BLAST 数据库文件所在目录”，例如里面应包含
`16S_ribosomal_RNA.nsq`、`16S_ribosomal_RNA.nin` 等文件，并不是固定要写成 `db`。

直接从 FASTA 批量鉴定：

```bash
python sanger16s.py run \
  --input ./fasta_files \
  --db-dir /data/blast/16s-db \
  --out ./results.csv
```

## CSV 输出字段

- `sample_file`
- `query_id`
- `query_length`
- `best_accession`
- `best_title`
- `organism_name`
- `percent_identity`
- `query_coverage`
- `alignment_length`
- `mismatches`
- `gaps`
- `evalue`
- `bitscore`
- `confidence`
- `note`

`confidence` 的含义：

- `species_high`
- `genus_or_close_species`
- `low_confidence`
- `no_result`

## 网页服务

网页服务支持两种任务模式：

- `FASTA -> BLAST`
- `FASTQ -> consensus -> BLAST`

启动示例：

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

更完整的网页部署说明见 [WEBAPP.zh-CN.md](./WEBAPP.zh-CN.md)。

## 部署示例配置

仓库中附带了公开版常用部署配置示例：

- `deploy/systemd/sanger16s-web.service.example`
- `deploy/nginx/sanger16s-web.conf.example`
- `deploy/sanger16s-web.env.example`

## 注意事项

- `--work-dir` 下会保存上传文件、任务状态和结果 CSV。
- FASTQ 模式依赖 `run_16s_consensus_pipeline.sh` 及 `scripts/` 下的 helper 脚本。
- 16S 序列并不总能分辨非常接近的物种，低置信度结果应谨慎解释。

## 许可证

MIT，见 [LICENSE](./LICENSE)。
