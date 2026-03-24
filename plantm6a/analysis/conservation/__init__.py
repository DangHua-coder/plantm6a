"""
Conservation analysis module for plantm6a（完整版）
"""

from .translate import translate_cds, translate, parse_fasta
from .extract_orthologs import (
    extract_pairwise_orthologs,
    parse_orthogroups,
    parse_ortholog_files
)

# 保守位点分析
try:
    from .conserved_sites import (
        TranscriptModel,
        parse_annotation,
        build_position_index,
        get_transcript_sequence,
        genomic_to_transcript_pos,
        load_m6a_sites,
        map_m6a_to_transcripts,
        sw_align,
        find_conserved_sites,
        analyze_species_pair,
    )
    CONSERVED_SITES_AVAILABLE = True
except ImportError as e:
    CONSERVED_SITES_AVAILABLE = False

# 汇总分析
try:
    from .summarize import summarize_conserved_m6a
    SUMMARIZE_AVAILABLE = True
except ImportError:
    SUMMARIZE_AVAILABLE = False

# 配置
try:
    from .config import SPECIES, PAIRWISE_DIR, OUTPUT_DIR, load_species_config
    CONFIG_AVAILABLE = True
except ImportError:
    CONFIG_AVAILABLE = False
# 可视化（需要matplotlib, pandas等）
try:
    from .visualize import plot_chord, SPECIES_ORDER, LABELS, SPECIES_COLORS
    VISUALIZE_AVAILABLE = True
except ImportError:
    VISUALIZE_AVAILABLE = False

__all__ = [
    'translate_cds',
    'translate',
    'parse_fasta',
    'extract_pairwise_orthologs',
    'parse_orthogroups',
    'parse_ortholog_files',
]

if CONSERVED_SITES_AVAILABLE:
    __all__.extend([
        'TranscriptModel',
        'parse_annotation',
        'build_position_index',
        'get_transcript_sequence',
        'genomic_to_transcript_pos',
        'load_m6a_sites',
        'map_m6a_to_transcripts',
        'sw_align',
        'find_conserved_sites',
        'analyze_species_pair',
    ])

if SUMMARIZE_AVAILABLE:
    __all__.append('summarize_conserved_m6a')

if CONFIG_AVAILABLE:
    __all__.extend(['SPECIES', 'PAIRWISE_DIR', 'OUTPUT_DIR', 'load_species_config'])

if VISUALIZE_AVAILABLE:
    __all__.extend(['plot_chord', 'SPECIES_ORDER', 'LABELS', 'SPECIES_COLORS'])
