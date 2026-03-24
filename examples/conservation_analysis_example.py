#!/usr/bin/env python3
"""
保守性分析使用示例
"""

import sys
from pathlib import Path

project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from plantm6a.analysis.conservation import translate_cds, extract_pairwise_orthologs


def example_translate():
    """示例: 翻译CDS序列"""
    print("="*60)
    print("示例 1: 翻译CDS为蛋白序列")
    print("="*60)
    
    translate_cds(
        input_fa="path/to/cds.fa",
        output_fa="path/to/protein.fa",
        verbose=True
    )


def example_extract_orthologs():
    """示例: 提取同源基因对"""
    print("\n" + "="*60)
    print("示例 2: 从OrthoFinder结果提取同源基因对")
    print("="*60)
    
    stats = extract_pairwise_orthologs(
        orthogroups_tsv="orthofinder_out/Results_*/Orthogroups/Orthogroups.tsv",
        orthologues_dir="orthofinder_out/Results_*/Orthologues/",
        output_dir="./pairwise_orthologs/",
        verbose=True
    )
    
    print(f"\n提取了 {len(stats)} 对物种的同源基因")


if __name__ == '__main__':
    print("PlantM6A 保守性分析示例")
    print("注意: 请修改文件路径为实际数据路径")
    print()
    
    # example_translate()
    # example_extract_orthologs()
