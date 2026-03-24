#!/usr/bin/env python3
"""
完整的保守m6A分析pipeline示例
"""

import sys
from pathlib import Path

project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from plantm6a.analysis.conservation import (
    analyze_species_pair,
    summarize_conserved_m6a,
)


def example_analyze_pair():
    """示例: 分析两个物种的保守m6A位点"""
    
    n_conserved = analyze_species_pair(
        spA='arabidopsis',
        spB='rice',
        genome_A='/path/to/arabidopsis.fa',
        genome_B='/path/to/rice.fa',
        annot_A='/path/to/arabidopsis.gtf',
        annot_B='/path/to/rice.gtf',
        fmt_A='gtf',
        fmt_B='gtf',
        m6a_A_file='/path/to/arabidopsis_m6a.bed',
        m6a_B_file='/path/to/rice_m6a.bed',
        ortholog_file='/path/to/arabidopsis_vs_rice.tsv',
        output_file='./arabidopsis_vs_rice_conserved_m6A.tsv',
        orthologs_only=True,
        verbose=True
    )
    
    print(f"\n找到 {n_conserved} 个保守m6A位点")


def example_summarize():
    """示例: 汇总所有物种对的结果"""
    
    stats = summarize_conserved_m6a(
        conserved_dir='./conserved_results/',
        output_dir='./summary/',
        verbose=True
    )
    
    print(f"\n汇总了 {len(stats.get('all_sites', []))} 个保守位点")


if __name__ == '__main__':
    print("PlantM6A 保守m6A分析完整示例")
    print("="*60)
    
    # 取消注释以运行
    # example_analyze_pair()
    # example_summarize()
    
    print("\n提示: 修改文件路径为实际数据路径后运行")
