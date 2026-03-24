"""
Analysis module for plantm6a
"""

from .statistics import analyze, count_A_in_exons, count_bed_sites_in_exons
from .batch import batch_analyze

# Motif分析模块（需要pysam）
try:
    from .motif import analyze_motifs, perform_clustering_analysis, generate_all_3mers_with_A
    MOTIF_AVAILABLE = True
except ImportError:
    MOTIF_AVAILABLE = False

# Conservation分析模块
try:
    from . import conservation
    CONSERVATION_AVAILABLE = True
except ImportError:
    CONSERVATION_AVAILABLE = False

if MOTIF_AVAILABLE:
    __all__ = [
        'analyze',
        'count_A_in_exons', 
        'count_bed_sites_in_exons',
        'batch_analyze',
        'analyze_motifs',
        'perform_clustering_analysis',
        'generate_all_3mers_with_A',
        'conservation',
    ]
else:
    __all__ = [
        'analyze',
        'count_A_in_exons', 
        'count_bed_sites_in_exons',
        'batch_analyze',
    ]

