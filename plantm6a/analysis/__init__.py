"""
Analysis module for plantm6a
"""

from .statistics import analyze, count_A_in_exons, count_bed_sites_in_exons
from .batch import batch_analyze

__all__ = [
    'analyze',
    'count_A_in_exons', 
    'count_bed_sites_in_exons',
    'batch_analyze',
]
