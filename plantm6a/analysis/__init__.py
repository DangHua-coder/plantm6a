"""
Analysis module for plantm6a.
"""

from .statistics import analyze, count_A_in_exons, count_bed_sites_in_exons

try:
    from .batch import batch_analyze
    BATCH_AVAILABLE = True
except ImportError:
    BATCH_AVAILABLE = False

try:
    from .motif import analyze_motifs, perform_clustering_analysis, generate_all_3mers_with_A
    from .site_annotation import annotate_m6a_sites
    MOTIF_AVAILABLE = True
except ImportError:
    MOTIF_AVAILABLE = False

try:
    from . import conservation
    CONSERVATION_AVAILABLE = True
except ImportError:
    CONSERVATION_AVAILABLE = False

try:
    from .metagene import plot_metagene, run_region2bin
    METAGENE_AVAILABLE = True
except ImportError:
    METAGENE_AVAILABLE = False

try:
    from .ejc import (
        DEFAULT_EJC_SPECIES_MAP,
        ExonJunctionTripletAnalyzer,
        parse_species_map_file,
        plot_ejc_triplet,
        run_ejc_batch,
        run_ejc_triplet,
    )
    EJC_AVAILABLE = True
except ImportError:
    EJC_AVAILABLE = False

try:
    from .differential import (
        find_differential_m6a,
        gene_level_summary,
        plot_differential_m6a,
        run_differential_m6a,
    )
    DIFFERENTIAL_AVAILABLE = True
except ImportError:
    DIFFERENTIAL_AVAILABLE = False

__all__ = [
    'analyze',
    'count_A_in_exons',
    'count_bed_sites_in_exons',
]

if BATCH_AVAILABLE:
    __all__.append('batch_analyze')

if MOTIF_AVAILABLE:
    __all__.extend([
        'analyze_motifs',
        'perform_clustering_analysis',
        'generate_all_3mers_with_A',
        'annotate_m6a_sites',
    ])

if CONSERVATION_AVAILABLE:
    __all__.append('conservation')

if METAGENE_AVAILABLE:
    __all__.extend(['plot_metagene', 'run_region2bin'])

if EJC_AVAILABLE:
    __all__.extend([
        'DEFAULT_EJC_SPECIES_MAP',
        'ExonJunctionTripletAnalyzer',
        'parse_species_map_file',
        'plot_ejc_triplet',
        'run_ejc_batch',
        'run_ejc_triplet',
    ])

if DIFFERENTIAL_AVAILABLE:
    __all__.extend([
        'find_differential_m6a',
        'gene_level_summary',
        'plot_differential_m6a',
        'run_differential_m6a',
    ])
