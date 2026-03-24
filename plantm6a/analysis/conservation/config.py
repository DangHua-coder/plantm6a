#!/usr/bin/env python3
"""
物种配置文件

用户需要根据实际数据路径修改此文件
"""

from pathlib import Path

# 默认路径（用户可以通过环境变量或配置文件覆盖）
BASE_DIR = Path.home() / "data_source" / "m6A_cross_species"

# 物种配置
# 每个物种需要提供：
#   - fasta: 基因组序列文件
#   - annot: 注释文件（GTF或GFF3）
#   - fmt: 注释格式 ('gtf' 或 'gff')
#   - m6a: m6A位点BED文件
SPECIES = {
    'arabidopsis': {
        'fasta': str(BASE_DIR / 'genomes' / 'arabidopsis.fa'),
        'annot': str(BASE_DIR / 'annotations' / 'arabidopsis.gtf'),
        'fmt': 'gtf',
        'm6a': str(BASE_DIR / 'm6a_sites' / 'arabidopsis_m6a.bed'),
    },
    'rice': {
        'fasta': str(BASE_DIR / 'genomes' / 'rice.fa'),
        'annot': str(BASE_DIR / 'annotations' / 'rice.gtf'),
        'fmt': 'gtf',
        'm6a': str(BASE_DIR / 'm6a_sites' / 'rice_m6a.bed'),
    },
    # 添加更多物种...
}

# OrthoFinder输出目录
PAIRWISE_DIR = str(BASE_DIR / 'pairwise_tables')

# 保守位点分析输出目录
OUTPUT_DIR = str(BASE_DIR / 'conserved_m6a_results')


def load_species_config(config_file=None):
    """
    从外部文件加载配置（可选）
    
    配置文件格式（YAML）：
    species:
      arabidopsis:
        fasta: /path/to/genome.fa
        annot: /path/to/annotation.gtf
        fmt: gtf
        m6a: /path/to/m6a_sites.bed
    """
    if config_file is None:
        return SPECIES, PAIRWISE_DIR, OUTPUT_DIR
    
    import yaml
    with open(config_file) as f:
        config = yaml.safe_load(f)
    
    return (
        config.get('species', SPECIES),
        config.get('pairwise_dir', PAIRWISE_DIR),
        config.get('output_dir', OUTPUT_DIR)
    )


__all__ = ['SPECIES', 'PAIRWISE_DIR', 'OUTPUT_DIR', 'load_species_config']
