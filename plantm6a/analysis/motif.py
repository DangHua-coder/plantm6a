#!/usr/bin/env python3
"""
m6A Motif Analysis Module

功能:
1. 简单模式: RAC/GAT/Others三类分类
2. 完整模式: 所有16种3-mer motif分析
3. 聚类模式: PCA、t-SNE、K-means、层次聚类
"""

import pysam
from collections import defaultdict, Counter
import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# 可选的机器学习库
try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE
    from sklearn.cluster import KMeans, AgglomerativeClustering
    from sklearn.preprocessing import StandardScaler
    from sklearn.ensemble import RandomForestClassifier
    from scipy.cluster.hierarchy import dendrogram, linkage
    from scipy.spatial.distance import pdist, squareform
    ML_AVAILABLE = True
except ImportError:
    ML_AVAILABLE = False


def parse_annotation_file(annotation_file):
    """解析注释文件（GTF或GFF3格式），提取exon区域"""
    exons = defaultdict(list)
    file_format = 'gff' if annotation_file.endswith(('.gff', '.gff3', '.gff.gz', '.gff3.gz')) else 'gtf'
    
    print(f"  解析 {file_format.upper()} 文件: {annotation_file}")
    
    opener = open
    if annotation_file.endswith('.gz'):
        import gzip
        opener = gzip.open
    
    with opener(annotation_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            chrom = fields[0]
            feature = fields[2]
            start = int(float(fields[3])) - 1
            end = int(float(fields[4]))
            strand = fields[6]
            
            if feature == 'exon' or (file_format == 'gff' and feature in ['mRNA', 'CDS']):
                exons[chrom].append((start, end, strand))
    
    merged_exons = {}
    for chrom in exons:
        merged_exons[chrom] = merge_overlapping_regions(exons[chrom])
    
    total_exons = sum(len(v) for v in merged_exons.values())
    print(f"    提取到 {total_exons} 个exon区域")
    
    return merged_exons


def merge_overlapping_regions(regions):
    """合并重叠的基因组区域（按链分开合并）"""
    if not regions:
        return []
    
    plus_regions = [r for r in regions if r[2] == '+']
    minus_regions = [r for r in regions if r[2] == '-']
    
    def merge_strand(strand_regions):
        if not strand_regions:
            return []
        sorted_regions = sorted(strand_regions, key=lambda x: (x[0], x[1]))
        merged = [sorted_regions[0]]
        
        for current in sorted_regions[1:]:
            last = merged[-1]
            if current[0] <= last[1]:
                merged[-1] = (last[0], max(last[1], current[1]), last[2])
            else:
                merged.append(current)
        return merged
    
    return merge_strand(plus_regions) + merge_strand(minus_regions)


def check_overlap(pos, regions):
    """检查位点是否与任何区域重叠"""
    for start, end, strand in regions:
        if start <= pos < end:
            return True, strand
    return False, None


def extract_sequence_context(fasta_file, chrom, pos, strand, window=2):
    """提取位点周围的序列"""
    try:
        seq = fasta_file.fetch(chrom, max(0, pos - window), pos + window + 1)
        seq = seq.upper()
        
        if len(seq) != 2 * window + 1:
            return None
        
        center_idx = len(seq) // 2
        
        if strand == '+':
            if seq[center_idx] != 'A':
                return None
            return seq
        else:
            if seq[center_idx] != 'T':
                return None
            rc_seq = reverse_complement(seq)
            if rc_seq[center_idx] != 'A':
                return None
            return rc_seq
            
    except Exception:
        return None


def reverse_complement(seq):
    """计算反向互补序列"""
    complement_map = str.maketrans('ATGCN', 'TACGN')
    return seq.translate(complement_map)[::-1]


def extract_3mer_motif(sequence):
    """从5bp序列中提取3mer motif"""
    if not sequence or len(sequence) < 5:
        return None, False
    
    motif = sequence[1:4]
    
    if len(motif) != 3 or motif[1] != 'A':
        return None, False
    
    if not all(base in 'ATGC' for base in motif):
        return None, False
    
    return motif, True


def classify_motif_simple(sequence):
    """简单模式: RAC/GAT/Others分类"""
    if not sequence or len(sequence) < 5:
        return 'unknown', sequence
    
    motif = sequence[1:4]
    
    if len(motif) != 3 or motif[1] != 'A':
        return 'not_A_centered', sequence
    
    if motif in ['AAC', 'GAC']:
        return 'RAC', motif
    elif motif == 'GAT':
        return 'GAT', motif
    else:
        return 'others', motif


def generate_all_3mers_with_A():
    """生成所有中心为A的3mer motif（共16种）"""
    bases = ['A', 'T', 'G', 'C']
    motifs = []
    for b1 in bases:
        for b3 in bases:
            motifs.append(f"{b1}A{b3}")
    return sorted(motifs)


def analyze_motifs(genome_path, annotation_path, bed_path, 
                   mode='simple', verbose=True):
    """
    分析m6A位点的motif分布
    
    参数:
        genome_path: 基因组FASTA文件
        annotation_path: GTF/GFF注释文件
        bed_path: m6A位点BED文件
        mode: 'simple' (RAC/GAT/Others) 或 'complete' (16种3-mer)
        verbose: 是否显示详细信息
    
    返回:
        结果字典
    """
    
    if verbose:
        print(f"\n{'='*60}")
        print(f"Motif分析 (模式: {mode})")
        print(f"{'='*60}")
    
    # 展开路径
    genome_path = str(Path(genome_path).expanduser())
    annotation_path = str(Path(annotation_path).expanduser())
    bed_path = str(Path(bed_path).expanduser())
    
    # 检查文件
    for path, desc in [(genome_path, "基因组"), (annotation_path, "注释"), (bed_path, "BED")]:
        if not Path(path).exists():
            raise FileNotFoundError(f"{desc}文件不存在: {path}")
    
    # 解析exon
    exons = parse_annotation_file(annotation_path)
    
    # 打开基因组
    if verbose:
        print(f"  加载基因组: {genome_path}")
    genome = pysam.FastaFile(genome_path)
    
    # 处理BED文件
    if verbose:
        print(f"  处理BED文件: {bed_path}")
    
    # 初始化计数器
    if mode == 'simple':
        stats = {
            'total_sites': 0,
            'in_exon': 0,
            'not_in_exon': 0,
            'RAC': 0,
            'GAT': 0,
            'others': 0,
            'unknown': 0,
            'not_A_centered': 0,
            'positive_strand': 0,
            'negative_strand': 0,
            'seq_extraction_failed': 0
        }
        motif_counter = None
    else:  # complete mode
        all_3mers = generate_all_3mers_with_A()
        motif_counter = Counter({motif: 0 for motif in all_3mers})
        stats = {
            'total_sites': 0,
            'in_exon': 0,
            'not_in_exon': 0,
            'valid_motifs': 0,
            'positive_strand': 0,
            'negative_strand': 0,
            'seq_extraction_failed': 0,
            'invalid_motif': 0
        }
    
    processed_positions = set()
    motif_sequences = []
    
    with open(bed_path, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('track'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            
            chrom = fields[0]
            start = int(float(fields[1]))
            end = int(float(fields[2]))
            bed_strand = fields[5] if len(fields) > 5 else None
            
            m6a_pos = (start + end) // 2
            stats['total_sites'] += 1
            
            # 检查exon
            if chrom not in exons:
                stats['not_in_exon'] += 1
                continue
            
            is_in_exon, exon_strand = check_overlap(m6a_pos, exons[chrom])
            
            if not is_in_exon:
                stats['not_in_exon'] += 1
                continue
            
            stats['in_exon'] += 1
            
            # 确定链方向
            if bed_strand and bed_strand in ['+', '-']:
                strand = bed_strand
            elif exon_strand:
                strand = exon_strand
            else:
                strand = '+'
            
            # 去重
            pos_key = (chrom, m6a_pos, strand)
            if pos_key in processed_positions:
                continue
            processed_positions.add(pos_key)
            
            # 统计链
            if strand == '+':
                stats['positive_strand'] += 1
            else:
                stats['negative_strand'] += 1
            
            # 提取序列
            seq_context = extract_sequence_context(genome, chrom, m6a_pos, strand, window=2)
            
            if seq_context:
                if mode == 'simple':
                    motif_type, motif_seq = classify_motif_simple(seq_context)
                    stats[motif_type] += 1
                    
                    if motif_type in ['RAC', 'GAT', 'others']:
                        motif_sequences.append({
                            'chrom': chrom,
                            'pos': m6a_pos,
                            'strand': strand,
                            'motif_type': motif_type,
                            'motif_seq': motif_seq,
                            'full_seq': seq_context
                        })
                else:  # complete mode
                    motif_3mer, is_valid = extract_3mer_motif(seq_context)
                    if is_valid:
                        stats['valid_motifs'] += 1
                        motif_counter[motif_3mer] += 1
                    else:
                        stats['invalid_motif'] += 1
            else:
                stats['seq_extraction_failed'] += 1
    
    genome.close()
    
    # 构建结果
    if mode == 'simple':
        exon_total = stats['RAC'] + stats['GAT'] + stats['others']
        
        result = {
            'total_sites': stats['total_sites'],
            'unique_sites': len(processed_positions),
            'not_in_exon': stats['not_in_exon'],
            'in_exon': stats['in_exon'],
            'valid_motifs': exon_total,
            'positive_strand': stats['positive_strand'],
            'negative_strand': stats['negative_strand'],
            'RAC_count': stats['RAC'],
            'RAC_percentage': (stats['RAC'] / exon_total * 100) if exon_total > 0 else 0,
            'GAT_count': stats['GAT'],
            'GAT_percentage': (stats['GAT'] / exon_total * 100) if exon_total > 0 else 0,
            'others_count': stats['others'],
            'others_percentage': (stats['others'] / exon_total * 100) if exon_total > 0 else 0,
            'motif_sequences': motif_sequences[:100]
        }
        
        if verbose:
            print(f"\n  统计结果:")
            print(f"    总位点: {stats['total_sites']:,}")
            print(f"    有效motif: {exon_total:,}")
            print(f"    RAC: {stats['RAC']:,} ({result['RAC_percentage']:.2f}%)")
            print(f"    GAT: {stats['GAT']:,} ({result['GAT_percentage']:.2f}%)")
            print(f"    Others: {stats['others']:,} ({result['others_percentage']:.2f}%)")
    
    else:  # complete mode
        all_3mers = generate_all_3mers_with_A()
        result = {
            'total_sites': stats['total_sites'],
            'in_exon': stats['in_exon'],
            'valid_motifs': stats['valid_motifs'],
            'positive_strand': stats['positive_strand'],
            'negative_strand': stats['negative_strand']
        }
        
        # 添加所有3mer计数和百分比
        for motif in all_3mers:
            count = motif_counter[motif]
            percentage = (count / stats['valid_motifs'] * 100) if stats['valid_motifs'] > 0 else 0
            result[f'{motif}_count'] = count
            result[f'{motif}_pct'] = percentage
        
        if verbose:
            print(f"\n  统计结果:")
            print(f"    总位点: {stats['total_sites']:,}")
            print(f"    有效motif: {stats['valid_motifs']:,}")
            print(f"    Top 5 motifs:")
            for motif, count in motif_counter.most_common(5):
                pct = (count / stats['valid_motifs'] * 100) if stats['valid_motifs'] > 0 else 0
                print(f"      {motif}: {count:,} ({pct:.2f}%)")
    
    return result


def perform_clustering_analysis(df_motifs, output_dir):
    """执行聚类分析（需要安装sklearn等库）"""
    if not ML_AVAILABLE:
        raise ImportError("聚类分析需要安装 matplotlib, seaborn, scikit-learn")
    
    print(f"\n{'='*60}")
    print("聚类分析")
    print(f"{'='*60}")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 提取特征
    all_3mers = generate_all_3mers_with_A()
    feature_cols = [f'{motif}_pct' for motif in all_3mers]
    
    X = df_motifs[feature_cols].values
    species_names = df_motifs['Species'].values
    
    # 标准化
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # 1. PCA
    print("\n1. PCA降维...")
    pca = PCA(n_components=min(len(species_names), len(feature_cols)))
    X_pca = pca.fit_transform(X_scaled)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    axes[0].scatter(X_pca[:, 0], X_pca[:, 1], s=100, alpha=0.6, 
                    c=range(len(species_names)), cmap='tab20')
    for i, species in enumerate(species_names):
        axes[0].annotate(species, (X_pca[i, 0], X_pca[i, 1]), 
                        fontsize=8, ha='center')
    axes[0].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
    axes[0].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
    axes[0].set_title('PCA: Species Distribution')
    axes[0].grid(True, alpha=0.3)
    
    axes[1].plot(range(1, len(pca.explained_variance_ratio_[:10])+1), 
                 np.cumsum(pca.explained_variance_ratio_[:10]), 'bo-')
    axes[1].set_xlabel('Number of Components')
    axes[1].set_ylabel('Cumulative Explained Variance')
    axes[1].set_title('PCA Variance Explained')
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'pca_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. t-SNE
    print("2. t-SNE降维...")
    tsne = TSNE(n_components=2, random_state=42, 
                perplexity=min(30, len(species_names)-1))
    X_tsne = tsne.fit_transform(X_scaled)
    
    plt.figure(figsize=(10, 8))
    plt.scatter(X_tsne[:, 0], X_tsne[:, 1], s=100, alpha=0.6, 
               c=range(len(species_names)), cmap='tab20')
    for i, species in enumerate(species_names):
        plt.annotate(species, (X_tsne[i, 0], X_tsne[i, 1]), 
                    fontsize=9, ha='center')
    plt.xlabel('t-SNE 1')
    plt.ylabel('t-SNE 2')
    plt.title('t-SNE: Species Distribution')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / 'tsne_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. K-means
    print("3. K-means聚类...")
    optimal_k = min(4, len(species_names) - 1)
    kmeans = KMeans(n_clusters=optimal_k, random_state=42, n_init=10)
    kmeans_labels = kmeans.fit_predict(X_scaled)
    
    # 4. 层次聚类
    print("4. 层次聚类...")
    linkage_matrix = linkage(X_scaled, method='ward')
    
    plt.figure(figsize=(12, 6))
    dendrogram(linkage_matrix, labels=species_names, 
              leaf_rotation=45, leaf_font_size=10)
    plt.xlabel('Species')
    plt.ylabel('Distance')
    plt.title('Hierarchical Clustering Dendrogram')
    plt.tight_layout()
    plt.savefig(output_dir / 'hierarchical_clustering.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 5. 距离矩阵
    print("5. 距离矩阵...")
    dist_matrix = squareform(pdist(X_scaled, metric='euclidean'))
    
    plt.figure(figsize=(10, 8))
    sns.heatmap(dist_matrix, xticklabels=species_names, yticklabels=species_names,
                cmap='YlOrRd', annot=True, fmt='.2f', square=True)
    plt.title('Species Distance Matrix (Euclidean)')
    plt.tight_layout()
    plt.savefig(output_dir / 'distance_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 保存结果
    clustering_results = pd.DataFrame({
        'Species': species_names,
        'KMeans_Cluster': kmeans_labels,
        'PC1': X_pca[:, 0],
        'PC2': X_pca[:, 1],
        'tSNE1': X_tsne[:, 0],
        'tSNE2': X_tsne[:, 1]
    })
    clustering_results.to_csv(output_dir / 'clustering_results.csv', index=False)
    
    print(f"\n聚类分析完成！结果保存在: {output_dir}")
    
    return clustering_results


__all__ = [
    'analyze_motifs',
    'perform_clustering_analysis',
    'generate_all_3mers_with_A'
]
