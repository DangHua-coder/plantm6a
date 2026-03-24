"""
Conservation analysis module for plantm6a

提供保守性分析功能，包括:
- CDS翻译
- OrthoFinder同源基因对提取
"""

from .translate import translate_cds, translate, parse_fasta
from .extract_orthologs import (
    extract_pairwise_orthologs,
    parse_orthogroups,
    parse_ortholog_files
)

__all__ = [
    'translate_cds',
    'translate',
    'parse_fasta',
    'extract_pairwise_orthologs',
    'parse_orthogroups',
    'parse_ortholog_files',
]
