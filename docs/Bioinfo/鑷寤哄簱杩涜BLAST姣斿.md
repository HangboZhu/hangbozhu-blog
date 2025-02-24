---
title: 自行建库进行Blast比对
description: 自行建库，在自己的数据上实现比对
tags: [linux, bioinfo]
sidebar_position: 1
---


# 自行建库进行Blast比对
## 1、下载 BLAST 软件
根据自身系统版本，下载对应的源代码，然后进行解压操作。


```bash
# 下载相关的源代码
wget -c https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.15.0/ncbi-blast-2.15.0+-x64-linux.tar.gz
# 解压
tar -zxvf ncbi-blast-2.15.0+-x64-linux.tar.gz
# 复制 bin 的路径，添加到环境变量中，此处省略具体操作
# 最后保存配置文件
source ~/.bashrc
```

使用 `help` 命令进行测试，如果能成功打印出相关信息，则表明安装成功。

```bash
$ blastn -h
USAGE
  blastn [-h] [-help] [-import_search_strategy filename]
    [-export_search_strategy filename] [-task task_name] [-db database_name]
    [-dbsize num_letters] [-gilist filename] [-seqidlist filename]
    [-negative_gilist filename] [-negative_seqidlist filename]
    [-taxids taxids] [-negative_taxids taxids] [-taxidlist filename]
    [-negative_taxidlist filename] [-entrez_query entrez_query]
    [-db_soft_mask filtering_algorithm] [-db_hard_mask filtering_algorithm]
    [-subject subject_input_file] [-subject_loc range] [-query input_file]
    [-out output_file] [-evalue evalue] [-word_size int_value]
    [-gapopen open_penalty] [-gapextend extend_penalty]
    [-perc_identity float_value] [-qcov_hsp_perc float_value]
    [-max_hsps int_value] [-xdrop_ungap float_value] [-xdrop_gap float_value]
    [-xdrop_gap_final float_value] [-searchsp int_value]
    [-sum_stats bool_value] [-penalty penalty] [-reward reward] [-no_greedy]
    [-min_raw_gapped_score int_value] [-template_type type]
    [-template_length int_value] [-dust DUST_options]
    [-filtering_db filtering_database]
    [-window_masker_taxid window_masker_taxid]
    [-window_masker_db window_masker_db] [-soft_masking soft_masking]
    [-ungapped] [-culling_limit int_value] [-best_hit_overhang float_value]
    [-best_hit_score_edge float_value] [-subject_besthit]
    [-window_size int_value] [-off_diagonal_range int_value]
    [-use_index boolean] [-index_name string] [-lcase_masking]
    [-query_loc range] [-strand strand] [-parse_deflines] [-outfmt format]
    [-show_gis] [-num_descriptions int_value] [-num_alignments int_value]
    [-line_length line_length] [-html] [-sorthits sort_hits]
    [-sorthsps sort_hsps] [-max_target_seqs num_sequences]
    [-num_threads int_value] [-mt_mode int_value] [-remote] [-version]

DESCRIPTION
   Nucleotide-Nucleotide BLAST 2.12.0+

Use '-help' to print detailed descriptions of command line arguments
```

## 2、构建私有数据库

### 2.1 使用 makeblastdb 工具构建数据库
安装 BLAST 软件并获取目标序列数据后，可使用 BLAST 提供的 `makeblastdb` 工具构建本地数据库。

`makeblastdb` 的基本命令格式如下：

```plaintext
makeblastdb -in input.fasta -dbtype nucl/prot -out database_name
```

参数说明：
- `-in`：指定输入的 FASTA 文件。
- `-dbtype`：指定数据库类型，`nucl` 表示核酸序列，`prot` 表示蛋白质序列。
- `-out`：指定输出数据库的名称。

例如，要构建一个 `prot` 类型的私有库，可使用以下命令：

```plaintext
makeblastdb -in chuand_version/anti-CRIPRdb-2022.fasta -dbtype prot -out anti-CRISPR_db_2022
```

### 2.2 检验构建的数据库的正确性
使用以下命令检验数据库的正确性：

```plaintext
blastdbcmd -info -db database_name
```

以之前构建的数据库为例：

```plaintext
$ blastdbcmd -info -db anti-CRISPR_db_2022

Database: chuand_version/anti-CRIPRdb-2022.fasta
	3,686 sequences; 456,340 total residues

Date: Feb 19, 2025  3:31 PM	Longest sequence: 322 residues

BLASTDB Version: 5

Volumes:/xxxx
```

查看数据库文件，将其放置在一个文件夹中：

```plaintext
$ ls -lh
总计 1.1M
-rwxrwxrwx 1 root root  20K  2月 19 15:39 anti-CRISPR_db_2022.pdb
-rwxrwxrwx 1 root root 497K  2月 19 15:39 anti-CRISPR_db_2022.phr
-rwxrwxrwx 1 root root  29K  2月 19 15:39 anti-CRISPR_db_2022.pin
-rwxrwxrwx 1 root root  44K  2月 19 15:39 anti-CRISPR_db_2022.pot
-rwxrwxrwx 1 root root 450K  2月 19 15:39 anti-CRISPR_db_2022.psq
-rwxrwxrwx 1 root root  16K  2月 19 15:39 anti-CRISPR_db_2022.ptf
-rwxrwxrwx 1 root root  15K  2月 19 15:39 anti-CRISPR_db_2022.pto
```

## 3、使用本地数据库进行比对
构建并验证数据库后，可使用 BLAST 工具进行本地搜索。例如，使用 `blastp` 进行蛋白序列比对：

```plaintext
blastp -query query_sequence.fasta -db database_name -out results.out
```

参数说明：
- `-query`：指定查询序列文件。
- `-db`：指定本地数据库名称。
- `-out`：指定输出结果文件。

直接输出的 `out` 文件信息较为冗杂，若只需要 `qseqid`、`sseqid`、`pident`、`length`、`mismatch`、`evalue`、`bitscore` 这些值，可修改输出的默认参数：

```plaintext
(
echo -e "QueryID\tSubjectID\tIdentity%\tAlignLen\tMismatches\tE-value\tBitScore"
blastp -query input.fa -db database \
       -outfmt "6 qseqid sseqid pident length mismatch evalue bitscore" \
       -max_target_seqs 29
) > blast_results.tsv
```

查看结果：

```plaintext
$ cat blast_results.tsv 
QueryID	SubjectID	Identity%	AlignLen	Mismatches	E-value	BitScore
BFH87998.1	anti_CRISPR3654_PWD66516.1_AcrIF22_I-F/I-E_Pectobacterium	37.500	32	20	1.8	25.0
BFH87998.1	anti_CRISPR0700_WP_089045469.1_AcrIIA7_II-A_Sinorhizobium	32.000	50	34	2.3	24.6
BFH87998.1	anti_CRISPR2188_AUA71579.1_AcrIF7_I-F_Pseudomonas	41.379	29	13	3.8	23.5
BFH87998.1	anti_CRISPR2187_AUA96137.1_AcrIF7_I-F_Pseudomonas	41.379	29	13	3.8	23.5
BFH87998.1	anti_CRISPR1283_WP_091963223.1_AcrIIA9_II-A_Prevotella	23.636	55	39	8.8	23.9
XKF41245.1	anti_CRISPR2788_EJU89670.1_AcrIIA5_II-A_Enterococcus	32.353	34	15	1.3	26.2
XKF41245.1	anti_CRISPR2428_CUM01408.1_AcrIIA1_II-A_Listeria	22.917	96	65	2.4	25.4
XKF41245.1	anti_CRISPR2335_OSA28332.1_AcrIIA1_II-A_Listeria	22.917	96	65	2.4	25.4
XKF41245.1	anti_CRISPR2334_PDG69377.1_AcrIIA1_II-A_Listeria	22.917	96	65	2.4	25.4
XKF41245.1	anti_CRISPR2333_OER47233.1_AcrIIA1_II-A_Listeria	22.917	96	65	2.4	25.4
XKF41245.1	anti_CRISPR2332_PCY40019.1_AcrIIA1_II-A_Listeria	22.917	96	65	2.4	25.4
XKF41245.1	anti_CRISPR2331_PCY66075.1_AcrIIA1_II-A_Listeria	22.917	96	65	2.4	25.4
XKF41245.1	anti_CRISPR2330_OFG88478.1_AcrIIA1_II-A_Listeria	22.917	96	65	2.4	25.4
XKF41245.1	anti_CRISPR2329_PDB12248.1_AcrIIA1_II-A_Listeria	22.917	96	65	2.4	25.4
XKF41245.1	anti_CRISPR2328_OEQ65283.1_AcrIIA1_II-A_Listeria	22.917	96	65	2.4	25.4
XKF41245.1	anti_CRISPR2327_MCT87412.1_AcrIIA1_II-A_Listeria	22.917	96	65	2.4	25.4
XKF41245.1	anti_CRISPR2326_PDC44312.1_AcrIIA1_II-A_Listeria	22.917	96	65	2.4	25.4
XKF41245.1	anti_CRISPR2287_MCK68445.1_AcrIIA1_II-A_Listeria	22.000	100	69	6.2	24.3
XKF41245.1	anti_CRISPR2286_MCK95540.1_AcrIIA1_II-A_Listeria	22.000	100	69	6.2	24.3
XKF41245.1	anti_CRISPR2285_MCJ57161.1_AcrIIA1_II-A_Listeria	22.000	100	69	6.2	24.3
XKF41245.1	anti_CRISPR2324_PCZ79994.1_AcrIIA1_II-A_Listeria	21.875	96	66	6.8	23.9
XKF41245.1	anti_CRISPR2340_MCV66918.1_AcrIIA1_II-A_Listeria	21.875	96	66	7.0	23.9
```

## 4、编写 Python 脚本将比对流程化
创建一个名为 `easy_custom_blastp.py` 的 Python 脚本：

```python
#!/usr/bin/env python3
import argparse
import subprocess
import os
import sys

def main():
    # 设置命令行参数解析
    parser = argparse.ArgumentParser(description="Run BLASTP and save results with headers.")
    parser.add_argument("-query", required=True, help="Input query FASTA file")
    parser.add_argument("-db", required=True, help="BLAST database name")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file name")
    parser.add_argument("-max_target_seqs", type=int, default=29,
                        help="Maximum number of target sequences (default: 29)")

    args = parser.parse_args()

    # 检查输入文件是否存在
    if not os.path.exists(args.query):
        sys.exit(f"Error: Query file '{args.query}' not found!")

    # 写入标题行
    header = "QueryID\tSubjectID\tIdentity%\tAlignLen\tMismatches\tE-value\tBitScore\n"
    with open(args.output, 'w') as f:
        f.write(header)

    # 构建 BLAST 命令
    blast_cmd = [
        "blastp",
        "-query", args.query,
        "-db", args.db,
        "-outfmt", "6 qseqid sseqid pident length mismatch evalue bitscore",
        "-max_target_seqs", str(args.max_target_seqs)
    ]

    # 执行 BLAST 并将结果追加到文件
    try:
        with open(args.output, 'a') as f:
            subprocess.run(blast_cmd, stdout=f, check=True)
        print(f"BLAST results saved to: {args.output}")
    except subprocess.CalledProcessError as e:
        os.remove(args.output)  # 删除不完整的输出文件
        sys.exit(f"BLASTP failed with error: {e}")
    except FileNotFoundError:
        sys.exit("Error: 'blastp' command not found. Ensure BLAST+ is installed and in PATH.")

if __name__ == "__main__":
    main()
```

查看脚本的参数设置：

```plaintext
$ python3 easy_custom_blastp.py -h
usage: easy_custom_blastp.py [-h] -query QUERY -db DB -o OUTPUT [-max_target_seqs MAX_TARGET_SEQS]

Run BLASTP and save results with headers.

options:
  -h, --help            show this help message and exit
  -query QUERY          Input query FASTA file
  -db DB                BLAST database name
  -o OUTPUT, --output OUTPUT
                        Output TSV file name
  -max_target_seqs MAX_TARGET_SEQS
                        Maximum number of target sequences (default: 29)
```

使用脚本进行比对：

```plaintext
$ python3 easy_custom_blastp.py -query ./test_query.fasta -db ./blast_database/anti-CRISPR_db_2022 -o ./test_results.tsv 
BLAST results saved to: ./test_results.tsv
```

至此，大功告成！🍻