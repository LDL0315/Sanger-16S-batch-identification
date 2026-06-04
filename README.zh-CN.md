# Sanger 16S 批量鉴定工具

这个工具用于批量处理 Sanger 测序得到的细菌 16S FASTA 文件，并用本地
NCBI 16S BLAST 数据库为每条序列找到最高一致性的命中结果，最终输出一个
CSV 汇总文件。

## 1. 在 Linux 服务器上安装环境

推荐使用 Conda 或 Mamba，不需要 `sudo` 权限。

```bash
conda env create -f environment.yml
conda activate sanger16s
```

如果使用 Mamba：

```bash
mamba env create -f environment.yml
conda activate sanger16s
```

## 2. 检查环境

```bash
python sanger16s.py check --db-dir db
```

第一次运行时，只要已经安装好 Conda 环境，通常会看到 `blastn` 和
`update_blastdb.pl` 可以找到，但数据库会提示缺失。这是正常的，因为本地
NCBI 16S 数据库还没有下载。

## 3. 下载本地 NCBI 16S 数据库

服务器需要能够连接 NCBI。

```bash
python sanger16s.py setup-db --db-dir db
```

这个命令会把 NCBI `16S_ribosomal_RNA` BLAST 数据库下载并解压到 `db/`
目录中。数据库只需要下载一次，后续批量分析可以重复使用。

## 4. 批量运行鉴定

把 FASTA 文件放到同一个文件夹中，例如：

```text
fasta_files/
  sample_001.fasta
  sample_002.fasta
  sample_003.fa
```

然后运行：

```bash
python sanger16s.py run --input fasta_files --db-dir db --out results.csv
```

脚本会递归读取 `.fasta`、`.fa`、`.fas` 和 `.fna` 文件。默认使用文件名作为
样本名。如果一个 FASTA 文件里有多条序列，每条序列会在结果 CSV 中单独占
一行。

## CSV 结果字段

输出 CSV 包含以下列：

- `sample_file`：输入 FASTA 文件名
- `query_id`：FASTA 中的序列 ID
- `query_length`：输入序列长度
- `best_accession`：最佳命中的 NCBI accession
- `best_title`：最佳命中的完整标题
- `organism_name`：最佳命中的物种或分类名称
- `percent_identity`：序列一致性百分比
- `query_coverage`：query 覆盖度
- `alignment_length`：比对长度
- `mismatches`：错配数量
- `gaps`：gap 数量
- `evalue`：BLAST E-value
- `bitscore`：BLAST bitscore
- `confidence`：鉴定置信等级
- `note`：备注，例如无命中、序列字符被替换、存在多个接近命中等

## 置信等级说明

最佳 BLAST 命中按照以下优先级选择：

1. `percent_identity` 最高
2. `query_coverage` 最高
3. `evalue` 最低
4. `bitscore` 最高

`confidence` 字段含义如下：

- `species_high`：identity >= 99.0 且 coverage >= 95，物种水平可信度较高。
- `genus_or_close_species`：identity >= 97.0 且 coverage >= 90，通常可作为属
  或近缘物种参考。
- `low_confidence`：低于上述阈值，但仍报告最高一致性的命中结果。
- `no_result`：空序列、空文件或没有 BLAST 命中。

如果多个物种水平命中非常接近，`note` 列会显示
`multiple close hits; review manually`，建议人工复核。

## 常用参数

指定 BLAST 使用的 CPU 线程数：

```bash
python sanger16s.py run --input fasta_files --db-dir db --out results.csv --threads 8
```

调整每条序列保留的候选命中数量，用于判断是否存在多个接近命中：

```bash
python sanger16s.py run --input fasta_files --db-dir db --out results.csv --max-target-seqs 20
```

## 注意事项

16S 序列并不总是能区分非常接近的细菌物种。因此，本工具报告的是“最高一致性
命中 + 置信等级”，不会把低置信度结果强行解释为确定的物种鉴定。
