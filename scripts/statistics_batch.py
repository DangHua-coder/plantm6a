#!/usr/bin/env python3
"""
批量分析多个物种的m6A修饰情况
"""

import os
import sys
import yaml
from pathlib import Path
import pandas as pd
from datetime import datetime

# 导入主分析函数
from count_A_and_m6A import analyze

def expand_path(path):
    """展开路径中的~符号"""
    return os.path.expanduser(path)

def load_config(config_file):
    """加载配置文件"""
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    return config['species']

def batch_analyze(config_file, output_file=None):
    """批量分析所有物种"""
    
    print("=" * 70)
    print("批量m6A修饰分析".center(70))
    print("=" * 70)
    print(f"\n开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"配置文件: {config_file}\n")
    
    # 加载配置
    species_list = load_config(config_file)
    print(f"共需要分析 {len(species_list)} 个物种/样本\n")
    
    # 存储所有结果
    all_results = []
    
    # 逐个分析
    for i, species in enumerate(species_list, 1):
        name = species['name']
        genome = expand_path(species['genome'])
        annotation = expand_path(species['annotation'])
        bed = expand_path(species['bed'])
        
        print("-" * 70)
        print(f"[{i}/{len(species_list)}] 正在分析: {name}")
        print("-" * 70)
        
        # 检查文件是否存在
        missing_files = []
        for file_path, file_type in [(genome, '基因组'), (annotation, '注释'), (bed, 'Bed')]:
            if not os.path.exists(file_path):
                missing_files.append(f"{file_type}: {file_path}")
        
        if missing_files:
            print(f"⚠️  警告: 以下文件不存在，跳过此物种:")
            for mf in missing_files:
                print(f"   - {mf}")
            print()
            continue
        
        try:
            # 执行分析（verbose=False不显示详细信息）
            results = analyze(genome, annotation, bed, verbose=False)
            
            # 添加物种名称
            results['species'] = name
            results['genome_file'] = os.path.basename(genome)
            results['annotation_file'] = os.path.basename(annotation)
            results['bed_file'] = os.path.basename(bed)
            
            all_results.append(results)
            
            # 显示关键结果
            print(f"✓ 完成!")
            print(f"  Exon中A:          {results['exon_A_count']:>12,}")
            print(f"  Exon内m6A位点:    {results['exon_m6A_sites']:>12,}")
            print(f"  修饰率:           {results['modification_rate']:>12.4f}%")
            print(f"  每1000个A的m6A:   {results['m6A_per_1000A']:>12.2f}")
            print()
            
        except Exception as e:
            print(f"❌ 错误: {e}")
            print()
            continue
    
    # 生成汇总表格
    if all_results:
        df = pd.DataFrame(all_results)
        
        # 重新排列列的顺序
        column_order = [
            'species',
            'exon_bases',
            'exon_A_count',
            'A_percentage',
            'total_m6A_sites',
            'exon_m6A_sites',
            'non_exon_sites',
            'exon_site_percentage',
            'modification_rate',
            'm6A_per_1000A',
            'm6A_per_1MA',
            'genome_file',
            'annotation_file',
            'bed_file'
        ]
        df = df[column_order]
        
        # 设置更友好的列名
        df.columns = [
            'Species',
            'Exon_Bases',
            'Exon_A_Count',
            'A_Percentage(%)',
            'Total_m6A_Sites',
            'Exon_m6A_Sites',
            'Non_Exon_Sites',
            'Exon_Site_Percentage(%)',
            'Modification_Rate(%)',
            'm6A_per_1000A',
            'm6A_per_1MA',
            'Genome_File',
            'Annotation_File',
            'Bed_File'
        ]
        
        # 输出到文件
        if output_file is None:
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            output_file = f'm6A_analysis_summary_{timestamp}.tsv'
        
        df.to_csv(output_file, sep='\t', index=False, float_format='%.4f')
        print("=" * 70)
        print(f"✓ 汇总结果已保存到: {output_file}")
        print("=" * 70)
        
        # 显示汇总表格
        print("\n汇总表格预览:\n")
        
        # 创建一个简化的显示版本
        display_df = df[[
            'Species',
            'Exon_A_Count',
            'Exon_m6A_Sites',
            'Modification_Rate(%)',
            'm6A_per_1000A'
        ]].copy()
        
        print(display_df.to_string(index=False))
        print()
        
        # 显示一些统计信息
        print("\n跨物种统计:")
        print(f"  平均修饰率: {df['Modification_Rate(%)'].mean():.4f}%")
        print(f"  修饰率范围: {df['Modification_Rate(%)'].min():.4f}% - {df['Modification_Rate(%)'].max():.4f}%")
        print(f"  平均每1000A的m6A: {df['m6A_per_1000A'].mean():.2f}")
        print()
        
    else:
        print("❌ 没有成功分析任何物种")
    
    print(f"\n结束时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

def main():
    if len(sys.argv) < 2:
        print("用法: python batch_analysis.py <config.yaml> [output.tsv]")
        print("\n示例:")
        print("  python batch_analysis.py species_config.yaml")
        print("  python batch_analysis.py species_config.yaml results.tsv")
        sys.exit(1)
    
    config_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    if not os.path.exists(config_file):
        print(f"错误: 配置文件不存在: {config_file}")
        sys.exit(1)
    
    batch_analyze(config_file, output_file)

if __name__ == '__main__':
    main()
