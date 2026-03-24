#!/usr/bin/env python3
"""
m6A Motif分析脚本 - 完全修复版
修复了：
1. 去重逻辑（现在考虑链方向）
2. 序列提取和验证逻辑
3. 统计计数问题
"""

import yaml
import pysam
from collections import defaultdict, Counter
import pandas as pd
from pathlib import Path
import sys


def parse_annotation_file(annotation_file):
    """
    解析注释文件（GTF或GFF3格式），提取exon区域
    返回: {chrom: [(start, end, strand), ...]}
    """
    exons = defaultdict(list)
    file_format = 'gff' if annotation_file.endswith(('.gff', '.gff3')) else 'gtf'
    
    print(f"  解析 {file_format.upper()} 文件: {annotation_file}")
    
    with open(annotation_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            chrom = fields[0]
            feature = fields[2]
            start = int(float(fields[3])) - 1  # Convert to 0-based
            end = int(float(fields[4]))
            strand = fields[6]
            
            # 提取exon特征
            if feature == 'exon' or (file_format == 'gff' and feature in ['mRNA', 'CDS']):
                exons[chrom].append((start, end, strand))
    
    # 合并重叠区域（保留链信息）
    merged_exons = {}
    for chrom in exons:
        merged_exons[chrom] = merge_overlapping_regions(exons[chrom])
    
    total_exons = sum(len(v) for v in merged_exons.values())
    print(f"    提取到 {total_exons} 个exon区域")
    
    return merged_exons


def merge_overlapping_regions(regions):
    """合并重叠的基因组区域（按链分开合并）"""
    if not regions:
        return []
    
    # 分别处理正链和负链
    plus_regions = [r for r in regions if r[2] == '+']
    minus_regions = [r for r in regions if r[2] == '-']
    
    def merge_strand(strand_regions):
        if not strand_regions:
            return []
        sorted_regions = sorted(strand_regions, key=lambda x: (x[0], x[1]))
        merged = [sorted_regions[0]]
        
        for current in sorted_regions[1:]:
            last = merged[-1]
            if current[0] <= last[1]:
                merged[-1] = (last[0], max(last[1], current[1]), last[2])
            else:
                merged.append(current)
        return merged
    
    return merge_strand(plus_regions) + merge_strand(minus_regions)


def check_overlap(pos, regions):
    """
    检查位点是否与任何区域重叠
    返回: (is_overlap, strand)
    """
    for start, end, strand in regions:
        if start <= pos < end:
            return True, strand
    return False, None


def extract_sequence_context(fasta_file, chrom, pos, strand, window=2):
    """
    提取位点周围的序列 - 完全修复版
    
    关键点：
    1. pos是m6A位点在基因组上的位置
    2. 对于正链：基因组序列 = mRNA序列，中心应该是A
    3. 对于负链：基因组序列中心是T，反向互补后得到mRNA序列，中心才是A
    
    window: 上下游窗口大小（默认2，提取5bp序列）
    返回：5bp mRNA序列，中心位置应该是A
    """
    try:
        # 提取基因组序列
        seq = fasta_file.fetch(chrom, max(0, pos - window), pos + window + 1)
        seq = seq.upper()
        
        # 验证序列长度
        if len(seq) != 2 * window + 1:
            return None
        
        center_idx = len(seq) // 2
        
        # 验证中心碱基
        if strand == '+':
            # 正链：基因组序列就是mRNA序列，中心必须是A
            if seq[center_idx] != 'A':
                return None
            return seq
        else:
            # 负链：基因组序列中心应该是T，反向互补后中心才是A
            if seq[center_idx] != 'T':
                return None
            # 反向互补得到mRNA序列
            rc_seq = reverse_complement(seq)
            # 再次验证（双重保险）
            if rc_seq[center_idx] != 'A':
                return None
            return rc_seq
            
    except Exception as e:
        return None


def reverse_complement(seq):
    """计算反向互补序列"""
    complement_map = str.maketrans('ATGCN', 'TACGN')
    return seq.translate(complement_map)[::-1]


def classify_motif_type(sequence):
    """
    分类motif类型
    
    序列格式: 5bp，中心位置(index=2)是A(m6A位点)
    提取中间3bp作为motif
    
    Motif定义:
    - RAC: R(A/G) + A + C  -> AAC 或 GAC
    - GAT: G + A + T
    - Others: 其他所有模式
    """
    if not sequence or len(sequence) < 5:
        return 'unknown', sequence
    
    # 提取中间3个碱基作为motif (index 1-3)
    motif = sequence[1:4]
    
    if len(motif) != 3:
        return 'unknown', sequence
    
    # 二次验证中心位置是A
    if motif[1] != 'A':
        return 'not_A_centered', sequence
    
    # 分类motif
    if motif in ['AAC', 'GAC']:  # RAC motif
        return 'RAC', motif
    elif motif == 'GAT':  # GAT motif
        return 'GAT', motif
    else:
        return 'others', motif


def analyze_single_species(name, genome_path, annotation_path, bed_path):
    """分析单个物种 - 完全修复版"""
    print(f"\n{'='*60}")
    print(f"分析物种: {name}")
    print(f"{'='*60}")
    
    # 展开路径
    genome_path = str(Path(genome_path).expanduser())
    annotation_path = str(Path(annotation_path).expanduser())
    bed_path = str(Path(bed_path).expanduser())
    
    # 检查文件存在性
    missing_files = []
    for path, desc in [(genome_path, "基因组"), (annotation_path, "注释"), (bed_path, "BED")]:
        if not Path(path).exists():
            missing_files.append(f"{desc}: {path}")
    
    if missing_files:
        print("  错误: 以下文件不存在:")
        for f in missing_files:
            print(f"    - {f}")
        return None
    
    # 1. 解析exon区域
    exons = parse_annotation_file(annotation_path)
    
    # 2. 打开基因组文件
    print(f"  加载基因组: {genome_path}")
    try:
        genome = pysam.FastaFile(genome_path)
    except Exception as e:
        print(f"  错误: 无法打开基因组文件 - {e}")
        return None
    
    # 3. 处理BED文件
    print(f"  处理BED文件: {bed_path}")
    
    stats = {
        'total_sites': 0,
        'in_exon': 0,
        'not_in_exon': 0,
        'RAC': 0,
        'GAT': 0,
        'others': 0,
        'unknown': 0,
        'not_A_centered': 0,
        'positive_strand': 0,
        'negative_strand': 0,
        'seq_extraction_failed': 0
    }
    
    motif_sequences = []
    processed_positions = set()  # 修复：现在包含链信息 (chrom, pos, strand)
    
    with open(bed_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#') or line.startswith('track'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            
            chrom = fields[0]
            start = int(float(fields[1]))
            end = int(float(fields[2]))
            
            # 获取BED文件的strand信息
            bed_strand = fields[5] if len(fields) > 5 else None
            
            # m6A位点位置（取中心）
            m6a_pos = (start + end) // 2
            stats['total_sites'] += 1
            
            # 检查是否在exon中
            if chrom not in exons:
                stats['not_in_exon'] += 1
                continue
            
            is_in_exon, exon_strand = check_overlap(m6a_pos, exons[chrom])
            
            if not is_in_exon:
                stats['not_in_exon'] += 1
                continue
            
            stats['in_exon'] += 1
            
            # 确定最终使用的链方向
            if bed_strand and bed_strand in ['+', '-']:
                strand = bed_strand
            elif exon_strand:
                strand = exon_strand
            else:
                strand = '+'  # 默认正链
            
            # 修复的去重逻辑：考虑链方向
            pos_key = (chrom, m6a_pos, strand)
            if pos_key in processed_positions:
                continue  # 跳过重复位点
            processed_positions.add(pos_key)
            
            # 统计链方向
            if strand == '+':
                stats['positive_strand'] += 1
            else:
                stats['negative_strand'] += 1
            
            # 提取序列上下文
            seq_context = extract_sequence_context(genome, chrom, m6a_pos, strand, window=2)
            
            if seq_context:
                # 分类motif
                motif_type, motif_seq = classify_motif_type(seq_context)
                stats[motif_type] += 1
                
                if motif_type in ['RAC', 'GAT', 'others']:
                    motif_sequences.append({
                        'chrom': chrom,
                        'pos': m6a_pos,
                        'strand': strand,
                        'motif_type': motif_type,
                        'motif_seq': motif_seq,
                        'full_seq': seq_context
                    })
            else:
                stats['seq_extraction_failed'] += 1
    
    genome.close()
    
    # 计算百分比（只统计在exon中且成功提取序列的位点）
    exon_total = stats['RAC'] + stats['GAT'] + stats['others']
    
    if exon_total == 0:
        print("  警告: 未找到有效的exon m6A位点")
        return None
    
    result = {
        'Species': name,
        'Total_m6A_sites': stats['total_sites'],
        'Unique_sites': len(processed_positions),
        'Not_in_exon': stats['not_in_exon'],
        'In_exon_total': stats['in_exon'],
        'Valid_motifs': exon_total,
        'Positive_strand': stats['positive_strand'],
        'Negative_strand': stats['negative_strand'],
        'RAC_count': stats['RAC'],
        'RAC_percentage': (stats['RAC'] / exon_total * 100) if exon_total > 0 else 0,
        'GAT_count': stats['GAT'],
        'GAT_percentage': (stats['GAT'] / exon_total * 100) if exon_total > 0 else 0,
        'Others_count': stats['others'],
        'Others_percentage': (stats['others'] / exon_total * 100) if exon_total > 0 else 0,
        'Seq_extraction_failed': stats['seq_extraction_failed'],
        'Not_A_centered_count': stats['not_A_centered']
    }
    
    # 打印统计信息
    print(f"\n  统计结果:")
    print(f"    总m6A位点: {stats['total_sites']:,}")
    print(f"    去重后位点: {len(processed_positions):,}")
    print(f"    不在exon中: {stats['not_in_exon']:,}")
    print(f"    在exon中: {stats['in_exon']:,}")
    print(f"    - 正链: {stats['positive_strand']:,}")
    print(f"    - 负链: {stats['negative_strand']:,}")
    print(f"    有效motif总数: {exon_total:,}")
    print(f"    Motif分布:")
    print(f"    - RAC: {stats['RAC']:,} ({result['RAC_percentage']:.2f}%)")
    print(f"    - GAT: {stats['GAT']:,} ({result['GAT_percentage']:.2f}%)")
    print(f"    - Others: {stats['others']:,} ({result['Others_percentage']:.2f}%)")
    if stats['seq_extraction_failed'] > 0:
        print(f"    序列提取失败: {stats['seq_extraction_failed']:,}")
    if stats['not_A_centered'] > 0:
        print(f"    非A中心: {stats['not_A_centered']:,}")
    
    # 打印前10个样例用于验证
    print(f"\n  样例motif序列（前10个）:")
    for i, m in enumerate(motif_sequences[:10], 1):
        print(f"    {i}. {m['chrom']}:{m['pos']} ({m['strand']}) -> {m['full_seq']} [{m['motif_seq']}] - {m['motif_type']}")
    
    result['motif_details'] = motif_sequences[:100]
    
    return result


def main():
    config_file = "./species_config.yaml"
    
    print("="*60)
    print("m6A Motif Distribution Analysis - Fully Fixed Version")
    print("="*60)
    
    # 读取配置
    print(f"\n读取配置文件: {config_file}")
    with open(config_file, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    
    print(f"发现 {len(config['species'])} 个物种")
    
    # 分析每个物种
    results = []
    for species_info in config['species']:
        result = analyze_single_species(
            species_info['name'],
            species_info['genome'],
            species_info['annotation'],
            species_info['bed']
        )
        
        if result:
            results.append(result)
    
    # 生成结果
    if not results:
        print("\n错误: 没有成功分析任何物种")
        return
    
    print(f"\n{'='*60}")
    print("生成结果文件")
    print(f"{'='*60}")
    
    # 创建DataFrame
    df = pd.DataFrame([{k: v for k, v in r.items() if k != 'motif_details'} 
                       for r in results])
    
    # 保存CSV
    output_csv = "./outputs/m6a_motif_distribution_fully_fixed.csv"
    df.to_csv(output_csv, index=False, float_format='%.2f')
    print(f"已保存CSV文件: {output_csv}")
    
    # 生成详细报告
    report_file = "./outputs/m6a_motif_analysis_report_fully_fixed.txt"
    with open(report_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("m6A Motif Distribution Analysis Report - Fully Fixed Version\n")
        f.write("="*80 + "\n\n")
        
        f.write("Key Fixes in This Version:\n")
        f.write("  1. Deduplication now considers strand information\n")
        f.write("  2. Fixed sequence extraction and validation logic:\n")
        f.write("     - Positive strand: genome seq = mRNA seq, center must be A\n")
        f.write("     - Negative strand: genome seq center is T, RC to get mRNA seq with A\n")
        f.write("  3. Improved statistics counting\n")
        f.write("  4. Added sample output for verification\n\n")
        
        f.write("Motif Classification Rules:\n")
        f.write("  - RAC: AAC or GAC (R = A/G)\n")
        f.write("  - GAT: GAT\n")
        f.write("  - Others: All other 3-mer patterns with A in the center\n\n")
        
        f.write("="*80 + "\n\n")
        
        for result in results:
            f.write(f"Species: {result['Species']}\n")
            f.write(f"{'-'*80}\n")
            f.write(f"  Total m6A sites: {result['Total_m6A_sites']:,}\n")
            f.write(f"  Unique sites (with strand): {result['Unique_sites']:,}\n")
            f.write(f"  Sites not in exons: {result['Not_in_exon']:,}\n")
            f.write(f"  Sites in exons: {result['In_exon_total']:,}\n")
            f.write(f"  Valid motifs extracted: {result['Valid_motifs']:,}\n")
            f.write(f"    - Positive strand: {result['Positive_strand']:,}\n")
            f.write(f"    - Negative strand: {result['Negative_strand']:,}\n\n")
            
            f.write(f"  Motif Distribution:\n")
            f.write(f"    RAC:    {result['RAC_count']:6,} ({result['RAC_percentage']:6.2f}%)\n")
            f.write(f"    GAT:    {result['GAT_count']:6,} ({result['GAT_percentage']:6.2f}%)\n")
            f.write(f"    Others: {result['Others_count']:6,} ({result['Others_percentage']:6.2f}%)\n")
            
            if result['Seq_extraction_failed'] > 0 or result['Not_A_centered_count'] > 0:
                f.write(f"\n  Quality Control:\n")
                if result['Seq_extraction_failed'] > 0:
                    f.write(f"    Sequence extraction failed: {result['Seq_extraction_failed']:,}\n")
                if result['Not_A_centered_count'] > 0:
                    f.write(f"    Not A-centered: {result['Not_A_centered_count']:,}\n")
            
            # 添加样例
            if result['motif_details']:
                f.write(f"\n  Sample motifs (first 10):\n")
                for i, m in enumerate(result['motif_details'][:10], 1):
                    f.write(f"    {i}. {m['chrom']}:{m['pos']} ({m['strand']}) -> {m['full_seq']} [{m['motif_seq']}]\n")
            
            f.write("\n" + "="*80 + "\n\n")
    
    print(f"已保存详细报告: {report_file}")
    
    # 显示汇总表格
    print("\n" + "="*60)
    print("汇总结果")
    print("="*60)
    display_columns = ['Species', 'Valid_motifs', 'Positive_strand', 'Negative_strand',
                      'RAC_count', 'RAC_percentage', 'GAT_count', 'GAT_percentage', 
                      'Others_count', 'Others_percentage']
    print(df[display_columns].to_string(index=False))
    
    print(f"\n分析完成！")
    print(f"结果文件保存在: ./outputs")


if __name__ == "__main__":
    main()
