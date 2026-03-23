#!/usr/bin/env python3
"""
统计基因组exon中A的含量和落在exon内的bed位点含量
修复版本：避免重叠exon导致的A碱基重复计数
"""

import sys
from collections import defaultdict
import gzip

def merge_overlapping_intervals(intervals):
    """
    合并重叠的区间
    输入: [(start1, end1), (start2, end2), ...]
    输出: 合并后的区间列表
    """
    if not intervals:
        return []
    
    # 按start位置排序
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    
    merged = []
    current_start, current_end = sorted_intervals[0]
    
    for start, end in sorted_intervals[1:]:
        if start <= current_end:  # 有重叠或相邻
            current_end = max(current_end, end)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = start, end
    
    merged.append((current_start, current_end))
    return merged

def parse_gtf_exons(gtf_file):
    """从GTF文件提取exon区域并合并重叠区域"""
    exons_by_chr = defaultdict(list)
    opener = gzip.open if gtf_file.endswith('.gz') else open
    
    try:
        with opener(gtf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                if fields[2] == 'exon':
                    chrom = fields[0]
                    start = int(fields[3]) - 1  # 转换为0-based
                    end = int(fields[4])
                    exons_by_chr[chrom].append((start, end))
    except Exception as e:
        print(f"Error parsing GTF file: {e}", file=sys.stderr)
        return exons_by_chr
    
    # 对每条染色体的exon进行排序和合并
    for chrom in exons_by_chr:
        exons_by_chr[chrom] = merge_overlapping_intervals(exons_by_chr[chrom])
    
    return exons_by_chr

def parse_gff_exons(gff_file):
    """从GFF文件提取exon区域并合并重叠区域"""
    exons_by_chr = defaultdict(list)
    opener = gzip.open if gff_file.endswith('.gz') else open
    
    try:
        with opener(gff_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                # GFF可能用不同的feature type
                if fields[2] in ['exon', 'CDS', 'mRNA']:
                    chrom = fields[0]
                    start = int(fields[3]) - 1
                    end = int(fields[4])
                    exons_by_chr[chrom].append((start, end))
    except Exception as e:
        print(f"Error parsing GFF file: {e}", file=sys.stderr)
        return exons_by_chr
    
    # 对每条染色体的exon进行排序和合并
    for chrom in exons_by_chr:
        exons_by_chr[chrom] = merge_overlapping_intervals(exons_by_chr[chrom])
    
    return exons_by_chr

def read_fasta(fasta_file):
    """读取FASTA文件"""
    sequences = {}
    current_chr = None
    current_seq = []
    
    opener = gzip.open if fasta_file.endswith('.gz') else open
    
    try:
        with opener(fasta_file, 'rt') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_chr:
                        sequences[current_chr] = ''.join(current_seq)
                    current_chr = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line.upper())
            if current_chr:
                sequences[current_chr] = ''.join(current_seq)
    except Exception as e:
        print(f"Error reading FASTA file: {e}", file=sys.stderr)
    
    return sequences

def is_in_exon(chrom, pos, exons_by_chr):
    """判断某个位点是否在exon区域内"""
    if chrom not in exons_by_chr:
        return False
    
    for start, end in exons_by_chr[chrom]:
        if start <= pos < end:
            return True
    return False

def count_A_in_exons(genome_seq, exons_by_chr):
    """统计exon区域中A的数量（已合并重叠区域，不会重复计数）"""
    total_A = 0
    total_bases = 0
    
    for chrom, exons in exons_by_chr.items():
        if chrom not in genome_seq:
            continue
        for start, end in exons:
            if end > len(genome_seq[chrom]):
                continue
            seq = genome_seq[chrom][start:end]
            total_A += seq.count('A')
            total_bases += len(seq)
    
    return total_A, total_bases

def count_bed_sites_in_exons(bed_file, exons_by_chr):
    """统计bed文件中落在exon区域内的位点数"""
    total_sites = 0
    exon_sites = 0
    
    opener = gzip.open if bed_file.endswith('.gz') else open
    
    try:
        with opener(bed_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                line = line.strip()
                if not line:
                    continue
                    
                total_sites += 1
                fields = line.split('\t')
                if len(fields) < 3:
                    continue
                
                chrom = fields[0]
                start = int(fields[1])  # bed文件是0-based
                pos = start  # 通常m6A位点用start位置
                
                if is_in_exon(chrom, pos, exons_by_chr):
                    exon_sites += 1
    except Exception as e:
        print(f"Error reading BED file: {e}", file=sys.stderr)
    
    return total_sites, exon_sites

def analyze(genome_file, annotation_file, bed_file, verbose=True):
    """执行分析并返回结果"""
    
    if verbose:
        print(f"正在分析...")
        print(f"基因组: {genome_file}")
        print(f"注释: {annotation_file}")
        print(f"Bed文件: {bed_file}")
        print("-" * 60)
    
    # 读取注释
    if verbose:
        print("\n[1/4] 读取注释文件并合并重叠exon...")
    
    if annotation_file.endswith('.gtf') or annotation_file.endswith('.gtf.gz'):
        exons_by_chr = parse_gtf_exons(annotation_file)
    else:
        exons_by_chr = parse_gff_exons(annotation_file)
    
    total_exons = sum(len(exons) for exons in exons_by_chr.values())
    if verbose:
        print(f"      找到 {len(exons_by_chr)} 条染色体")
        print(f"      合并后共 {total_exons:,} 个非重叠exon区域")
    
    # 读取基因组
    if verbose:
        print("\n[2/4] 读取基因组序列...")
    genome_seq = read_fasta(genome_file)
    if verbose:
        print(f"      读取了 {len(genome_seq)} 条染色体序列")
    
    # 统计exon中的A
    if verbose:
        print("\n[3/4] 统计exon中的A碱基（无重复计数）...")
    total_A, total_bases = count_A_in_exons(genome_seq, exons_by_chr)
    
    # 统计bed位点
    if verbose:
        print("\n[4/4] 统计bed文件中的m6A位点...")
    total_sites, exon_sites = count_bed_sites_in_exons(bed_file, exons_by_chr)
    
    # 计算结果
    results = {
        'exon_bases': total_bases,
        'exon_A_count': total_A,
        'A_percentage': total_A/total_bases*100 if total_bases > 0 else 0,
        'total_m6A_sites': total_sites,
        'exon_m6A_sites': exon_sites,
        'non_exon_sites': total_sites - exon_sites,
        'exon_site_percentage': exon_sites/total_sites*100 if total_sites > 0 else 0,
        'modification_rate': exon_sites/total_A*100 if total_A > 0 else 0,
        'm6A_per_1000A': exon_sites/total_A*1000 if total_A > 0 else 0,
        'm6A_per_1MA': exon_sites/total_A*1e6 if total_A > 0 else 0,
    }
    
    # 输出结果
    if verbose:
        print("\n" + "=" * 60)
        print("统计结果".center(60))
        print("=" * 60)
        
        print(f"\n【Exon区域统计】")
        print(f"  Exon总碱基数:        {results['exon_bases']:>15,} bp")
        print(f"  Exon中A的数量:       {results['exon_A_count']:>15,}")
        print(f"  A含量比例:           {results['A_percentage']:>15.2f} %")
        
        print(f"\n【m6A位点统计】")
        print(f"  Bed文件总位点数:     {results['total_m6A_sites']:>15,}")
        print(f"  落在exon内的位点:    {results['exon_m6A_sites']:>15,}")
        print(f"  不在exon内的位点:    {results['non_exon_sites']:>15,}")
        print(f"  Exon内位点比例:      {results['exon_site_percentage']:>15.2f} %")
        
        print(f"\n【修饰率计算】")
        print(f"  修饰率 (m6A/A):      {results['modification_rate']:>15.4f} %")
        print(f"  每1000个A的m6A数:    {results['m6A_per_1000A']:>15.2f}")
        print(f"  每100万个A的m6A数:   {results['m6A_per_1MA']:>15.0f}")
        
        print("\n" + "=" * 60 + "\n")
    
    return results

def main():
    if len(sys.argv) != 4:
        print("用法: python count_A_and_m6A.py <genome.fa> <annotation.gtf/gff> <sites.bed>")
        sys.exit(1)
    
    genome_file = sys.argv[1]
    annotation_file = sys.argv[2]
    bed_file = sys.argv[3]
    
    analyze(genome_file, annotation_file, bed_file, verbose=True)

if __name__ == '__main__':
    main()
