# PlantM6A

A comprehensive toolkit for plant m6A RNA modification analysis.

## Features

- Statistical analysis of m6A sites in plant genomes
- Exon junction analysis
- Motif analysis
- Batch processing for multiple species

## Installation
```bash
# Coming soon
pip install plantm6a
```

## Quick Start
```bash
# Analyze a single species
python scripts/batch_analysis.py config.yaml
```

## Development Status

🚧 This project is under active development.

## Author

[Your Name] - PhD student at [Your University]

## License

MIT License

## Motif Analysis

PlantM6A provides comprehensive motif analysis for m6A sites.

### Features
- **Simple mode**: RAC/GAT/Others classification
- **Complete mode**: All 16 possible 3-mer motifs (XAY patterns)
- **Clustering analysis**: PCA, t-SNE, hierarchical clustering for multi-species comparison

### Usage

#### Simple Mode
```python
from plantm6a.analysis import analyze_motifs

result = analyze_motifs(
    genome_path="genome.fa",
    annotation_path="annotation.gtf",
    bed_path="sites.bed",
    mode='simple',
    verbose=True
)

print(f"RAC: {result['RAC_percentage']:.2f}%")
print(f"GAT: {result['GAT_percentage']:.2f}%")
```

#### Complete Mode (All 16 3-mers)
```python
result = analyze_motifs(
    genome_path="genome.fa",
    annotation_path="annotation.gtf",
    bed_path="sites.bed",
    mode='complete',
    verbose=True
)

# Access specific motif
print(f"AAC: {result['AAC_count']} ({result['AAC_pct']:.2f}%)")
```

#### Multi-species Clustering
```python
from plantm6a.analysis import perform_clustering_analysis
import pandas as pd

# Analyze multiple species
results = []
for species_info in species_list:
    result = analyze_motifs(
        genome_path=species_info['genome'],
        annotation_path=species_info['gtf'],
        bed_path=species_info['bed'],
        mode='complete',
        verbose=False
    )
    result['Species'] = species_info['name']
    results.append(result)

df = pd.DataFrame(results)
clustering = perform_clustering_analysis(df, output_dir='./clustering_results')
```

### Dependencies
- Core: `pysam`
- Clustering (optional): `matplotlib`, `seaborn`, `scikit-learn`, `scipy`
