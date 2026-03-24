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
## Conservation Analysis (Complete Pipeline)

### Overview
Complete pipeline for identifying conserved m6A sites between species using orthologous genes.

### Features
- Full transcript analysis (5'UTR + CDS + 3'UTR)
- Genomic to transcript coordinate mapping
- Smith-Waterman sequence alignment
- Conserved site identification with region annotation
- Multi-species conservation summarization
- Beautiful Circos-style chord diagram visualization

### Quick Start

#### 1. Translate CDS to Proteins
```python
from plantm6a.analysis.conservation import translate_cds

translate_cds("species1_cds.fa", "species1_protein.fa")
```

#### 2. Run OrthoFinder
```bash
bash scripts/run_orthofinder.sh 40  # 40 threads
```

#### 3. Extract Pairwise Orthologs
```python
from plantm6a.analysis.conservation import extract_pairwise_orthologs

extract_pairwise_orthologs(
    orthogroups_tsv="orthofinder_out/Results_*/Orthogroups/Orthogroups.tsv",
    orthologues_dir="orthofinder_out/Results_*/Orthologues/",
    output_dir="./pairwise_orthologs/"
)
```

#### 4. Identify Conserved m6A Sites
```python
from plantm6a.analysis.conservation import analyze_species_pair

n_conserved = analyze_species_pair(
    spA='arabidopsis', spB='rice',
    genome_A='arabidopsis.fa', genome_B='rice.fa',
    annot_A='arabidopsis.gtf', annot_B='rice.gtf',
    fmt_A='gtf', fmt_B='gtf',
    m6a_A_file='arabidopsis_m6a.bed',
    m6a_B_file='rice_m6a.bed',
    ortholog_file='arabidopsis_vs_rice.tsv',
    output_file='arabidopsis_vs_rice_conserved_m6A.tsv'
)
```

#### 5. Summarize Results
```python
from plantm6a.analysis.conservation import summarize_conserved_m6a

stats = summarize_conserved_m6a(
    conserved_dir='./conserved_results/',
    output_dir='./summary/'
)
```

### Batch Processing Scripts

For analyzing multiple species pairs, use the complete pipeline scripts:
```bash
# Full conservation analysis
python scripts/conserved_m6A_pipeline.py \
    --pairwise_dir ./pairwise_tables \
    --output_dir ./conserved_results \
    --threads 8

# Comprehensive summarization  
python scripts/summarize_conserved_m6A.py \
    --conserved_dir ./conserved_results \
    --output_dir ./summary

# Chord diagram visualization
python scripts/plot_chord_v4.py \
    --input ./summary/pair_region_breakdown.tsv \
    --output ./figures/conserved_m6A_chord \
    --min_sites 100
```

### Dependencies
```bash
# Required
pip install parasail  # or biopython
conda install -c bioconda samtools orthofinder

# For visualization
pip install matplotlib seaborn pandas numpy
```
