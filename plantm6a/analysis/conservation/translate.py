#!/usr/bin/env python3
"""
CDS翻译模块
将CDS核苷酸序列翻译为蛋白序列
"""

import sys
import re
from pathlib import Path


def parse_fasta(fasta_path):
    """逐条解析fasta，yield (header, seq)"""
    header, seq_parts = None, []
    
    opener = open
    if str(fasta_path).endswith('.gz'):
        import gzip
        opener = gzip.open
    
    with opener(fasta_path, 'rt') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line.upper())
    if header is not None:
        yield header, "".join(seq_parts)


CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def translate(seq, frame=0):
    """
    将核苷酸序列翻译为氨基酸序列
    
    参数:
        seq: 核苷酸序列
        frame: 翻译框 (0, 1, 2)
    
    返回:
        蛋白序列
    """
    seq = seq.upper()[frame:]
    pep = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa = CODON_TABLE.get(codon, 'X')
        if aa == '*':
            break
        pep.append(aa)
    return "".join(pep)


def translate_cds(input_fa, output_fa, verbose=True):
    """
    翻译CDS文件
    
    参数:
        input_fa: 输入CDS fasta文件
        output_fa: 输出蛋白fasta文件
        verbose: 是否显示进度
    
    返回:
        (written, skipped): 写入和跳过的序列数
    """
    written = 0
    skipped = 0
    
    with open(output_fa, 'w') as out:
        for header, seq in parse_fasta(input_fa):
            gene_id = header.split()[0]
            
            # 尝试提取CDS坐标（1-based, inclusive）
            cds_match = re.search(r'CDS=(\d+)-(\d+)', header)
            if cds_match:
                start = int(cds_match.group(1)) - 1  # 转为0-based
                end = int(cds_match.group(2))
                cds_seq = seq[start:end]
            else:
                cds_seq = seq

            pep = translate(cds_seq)
            if len(pep) == 0:
                skipped += 1
                continue
            
            out.write(f">{gene_id}\n")
            # 每60个字符换行
            for i in range(0, len(pep), 60):
                out.write(pep[i:i+60] + "\n")
            written += 1
    
    if verbose:
        print(f"  {Path(input_fa).name}: 写入 {written} 条蛋白序列，跳过 {skipped} 条")
    
    return written, skipped


def main():
    """命令行接口"""
    if len(sys.argv) != 3:
        print("用法: python3 translate.py <input_cds.fa> <output_pep.fa>")
        sys.exit(1)
    
    input_fa, output_fa = sys.argv[1], sys.argv[2]
    translate_cds(input_fa, output_fa, verbose=True)


if __name__ == "__main__":
    main()


__all__ = ['translate_cds', 'translate', 'parse_fasta']
