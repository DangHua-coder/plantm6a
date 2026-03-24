#!/usr/bin/env python3
"""
m6A 全motif分析和聚类脚本
功能：
1. 统计所有64种3碱基motif的分布
2. 使用多种方法进行聚类分析（PCA、K-means、层次聚类、t-SNE等）
3. 使用分类器评估物种间差异（Random Forest、SVM）
"""

import yaml
import pysam
from collections import defaultdict, Counter
import pandas as pd
import numpy as np
from pathlib import Path
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans, AgglomerativeClustering, DBSCAN
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.model_selection import cross_val_score
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform
import warnings
warnings.filterwarnings('ignore')


def parse_annotation_file(annotation_file):
    """解析注释文件（GTF或GFF3格式），提取exon区域"""
    exons = defaultdict(list)
    file_format = 'gff' if annotation_file.endswith(('.gff', '.gff3')) else 'gtf'
    
    print(f"  解析 {file_format.upper()} 文件: {annotation_file}")
    
    with open(annotation_file, 'r') as f:
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
            
    except Exception as e:
        return None


def reverse_complement(seq):
    """计算反向互补序列"""
    complement_map = str.maketrans('ATGCN', 'TACGN')
    return seq.translate(complement_map)[::-1]


def extract_3mer_motif(sequence):
    """
    从5bp序列中提取3mer motif
    返回: (motif_3mer, is_valid)
    """
    if not sequence or len(sequence) < 5:
        return None, False
    
    motif = sequence[1:4]
    
    if len(motif) != 3:
        return None, False
    
    # 验证中心位置是A
    if motif[1] != 'A':
        return None, False
    
    # 验证只包含ATGC
    if not all(base in 'ATGC' for base in motif):
        return None, False
    
    return motif, True


def generate_all_3mers_with_A():
    """生成所有中心为A的3mer motif（共16种：XAY，X和Y可以是A/T/G/C）"""
    bases = ['A', 'T', 'G', 'C']
    motifs = []
    for b1 in bases:
        for b3 in bases:
            motifs.append(f"{b1}A{b3}")
    return sorted(motifs)


def analyze_single_species(name, genome_path, annotation_path, bed_path):
    """分析单个物种，统计所有3mer motif"""
    print(f"\n{'='*60}")
    print(f"分析物种: {name}")
    print(f"{'='*60}")
    
    genome_path = str(Path(genome_path).expanduser())
    annotation_path = str(Path(annotation_path).expanduser())
    bed_path = str(Path(bed_path).expanduser())
    
    missing_files = []
    for path, desc in [(genome_path, "基因组"), (annotation_path, "注释"), (bed_path, "BED")]:
        if not Path(path).exists():
            missing_files.append(f"{desc}: {path}")
    
    if missing_files:
        print("  错误: 以下文件不存在:")
        for f in missing_files:
            print(f"    - {f}")
        return None
    
    exons = parse_annotation_file(annotation_path)
    
    print(f"  加载基因组: {genome_path}")
    try:
        genome = pysam.FastaFile(genome_path)
    except Exception as e:
        print(f"  错误: 无法打开基因组文件 - {e}")
        return None
    
    print(f"  处理BED文件: {bed_path}")
    
    # 初始化所有可能的3mer计数
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
            
            if chrom not in exons:
                stats['not_in_exon'] += 1
                continue
            
            is_in_exon, exon_strand = check_overlap(m6a_pos, exons[chrom])
            
            if not is_in_exon:
                stats['not_in_exon'] += 1
                continue
            
            stats['in_exon'] += 1
            
            if bed_strand and bed_strand in ['+', '-']:
                strand = bed_strand
            elif exon_strand:
                strand = exon_strand
            else:
                strand = '+'
            
            pos_key = (chrom, m6a_pos, strand)
            if pos_key in processed_positions:
                continue
            processed_positions.add(pos_key)
            
            if strand == '+':
                stats['positive_strand'] += 1
            else:
                stats['negative_strand'] += 1
            
            seq_context = extract_sequence_context(genome, chrom, m6a_pos, strand, window=2)
            
            if seq_context:
                motif_3mer, is_valid = extract_3mer_motif(seq_context)
                
                if is_valid and motif_3mer:
                    motif_counter[motif_3mer] += 1
                    stats['valid_motifs'] += 1
                else:
                    stats['invalid_motif'] += 1
            else:
                stats['seq_extraction_failed'] += 1
    
    genome.close()
    
    if stats['valid_motifs'] == 0:
        print("  警告: 未找到有效的motif")
        return None
    
    # 构建结果
    result = {
        'Species': name,
        'Total_m6A_sites': stats['total_sites'],
        'Unique_sites': len(processed_positions),
        'In_exon_total': stats['in_exon'],
        'Valid_motifs': stats['valid_motifs'],
        'Positive_strand': stats['positive_strand'],
        'Negative_strand': stats['negative_strand'],
    }
    
    # 添加所有3mer的计数和百分比
    for motif in all_3mers:
        count = motif_counter[motif]
        percentage = (count / stats['valid_motifs'] * 100) if stats['valid_motifs'] > 0 else 0
        result[f'{motif}_count'] = count
        result[f'{motif}_pct'] = percentage
    
    print(f"\n  统计结果:")
    print(f"    总m6A位点: {stats['total_sites']:,}")
    print(f"    有效motif: {stats['valid_motifs']:,}")
    print(f"    Top 5 motifs:")
    for motif, count in motif_counter.most_common(5):
        pct = (count / stats['valid_motifs'] * 100)
        print(f"      {motif}: {count:,} ({pct:.2f}%)")
    
    return result


def perform_clustering_analysis(df_motifs, output_dir):
    """执行多种聚类分析"""
    print(f"\n{'='*60}")
    print("聚类分析")
    print(f"{'='*60}")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 准备数据：提取所有3mer的百分比作为特征
    all_3mers = generate_all_3mers_with_A()
    feature_cols = [f'{motif}_pct' for motif in all_3mers]
    
    X = df_motifs[feature_cols].values
    species_names = df_motifs['Species'].values
    
    # 标准化
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # 1. PCA分析
    print("\n1. PCA降维分析...")
    pca = PCA(n_components=min(len(species_names), len(feature_cols)))
    X_pca = pca.fit_transform(X_scaled)
    
    # 绘制PCA结果
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # PCA散点图
    axes[0].scatter(X_pca[:, 0], X_pca[:, 1], s=100, alpha=0.6, c=range(len(species_names)), cmap='tab20')
    for i, species in enumerate(species_names):
        axes[0].annotate(species, (X_pca[i, 0], X_pca[i, 1]), fontsize=8, ha='center')
    axes[0].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
    axes[0].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
    axes[0].set_title('PCA: Species Distribution')
    axes[0].grid(True, alpha=0.3)
    
    # PCA方差解释
    axes[1].plot(range(1, len(pca.explained_variance_ratio_[:10])+1), 
                 np.cumsum(pca.explained_variance_ratio_[:10]), 'bo-')
    axes[1].set_xlabel('Number of Components')
    axes[1].set_ylabel('Cumulative Explained Variance')
    axes[1].set_title('PCA Variance Explained')
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'pca_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. t-SNE分析
    print("2. t-SNE降维分析...")
    tsne = TSNE(n_components=2, random_state=42, perplexity=min(30, len(species_names)-1))
    X_tsne = tsne.fit_transform(X_scaled)
    
    plt.figure(figsize=(10, 8))
    plt.scatter(X_tsne[:, 0], X_tsne[:, 1], s=100, alpha=0.6, c=range(len(species_names)), cmap='tab20')
    for i, species in enumerate(species_names):
        plt.annotate(species, (X_tsne[i, 0], X_tsne[i, 1]), fontsize=9, ha='center')
    plt.xlabel('t-SNE 1')
    plt.ylabel('t-SNE 2')
    plt.title('t-SNE: Species Distribution')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / 'tsne_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. K-means聚类
    print("3. K-means聚类...")
    optimal_k = min(4, len(species_names) - 1)
    kmeans = KMeans(n_clusters=optimal_k, random_state=42, n_init=10)
    kmeans_labels = kmeans.fit_predict(X_scaled)
    
    # 4. 层次聚类
    print("4. 层次聚类分析...")
    linkage_matrix = linkage(X_scaled, method='ward')
    
    plt.figure(figsize=(12, 6))
    dendrogram(linkage_matrix, labels=species_names, leaf_rotation=45, leaf_font_size=10)
    plt.xlabel('Species')
    plt.ylabel('Distance')
    plt.title('Hierarchical Clustering Dendrogram')
    plt.tight_layout()
    plt.savefig(output_dir / 'hierarchical_clustering.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 5. 距离矩阵热图
    print("5. 生成物种间距离热图...")
    dist_matrix = squareform(pdist(X_scaled, metric='euclidean'))
    
    plt.figure(figsize=(10, 8))
    sns.heatmap(dist_matrix, xticklabels=species_names, yticklabels=species_names,
                cmap='YlOrRd', annot=True, fmt='.2f', square=True)
    plt.title('Species Distance Matrix (Euclidean)')
    plt.tight_layout()
    plt.savefig(output_dir / 'distance_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 6. 特征重要性（使用Random Forest）
    print("6. 特征重要性分析（Random Forest）...")
    if len(species_names) >= 3:
        # 创建虚拟标签用于分类
        rf = RandomForestClassifier(n_estimators=100, random_state=42)
        rf.fit(X_scaled, range(len(species_names)))
        
        feature_importance = pd.DataFrame({
            'Motif': all_3mers,
            'Importance': rf.feature_importances_
        }).sort_values('Importance', ascending=False)
        
        plt.figure(figsize=(12, 6))
        top_features = feature_importance.head(16)
        plt.barh(range(len(top_features)), top_features['Importance'])
        plt.yticks(range(len(top_features)), top_features['Motif'])
        plt.xlabel('Feature Importance')
        plt.title('Top 16 Discriminative Motifs (Random Forest)')
        plt.tight_layout()
        plt.savefig(output_dir / 'feature_importance.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        feature_importance.to_csv(output_dir / 'feature_importance.csv', index=False)
    
    # 保存聚类结果
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
    
    return {
        'pca': pca,
        'X_pca': X_pca,
        'kmeans_labels': kmeans_labels,
        'linkage_matrix': linkage_matrix
    }


def main():
    config_file = "./species_config.yaml"
    
    print("="*60)
    print("m6A Complete 3-mer Motif Analysis with Clustering")
    print("="*60)
    
    print(f"\n读取配置文件: {config_file}")
    with open(config_file, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    
    print(f"发现 {len(config['species'])} 个物种")
    
    # 分析每个物种
    results = []
    for species_info in config['species']:
        result = analyze_single_species(
            species_info['name'],
            species_info['genome'],
            species_info['annotation'],
            species_info['bed']
        )
        
        if result:
            results.append(result)
    
    if not results:
        print("\n错误: 没有成功分析任何物种")
        return
    
    print(f"\n{'='*60}")
    print("生成结果文件")
    print(f"{'='*60}")
    
    output_dir = Path("./outputs")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 创建完整DataFrame
    df = pd.DataFrame(results)
    
    # 保存完整结果
    output_csv = output_dir / "m6a_all_3mer_motifs.csv"
    df.to_csv(output_csv, index=False, float_format='%.4f')
    print(f"已保存完整CSV文件: {output_csv}")
    
    # 提取motif百分比数据用于聚类
    all_3mers = generate_all_3mers_with_A()
    pct_cols = ['Species'] + [f'{motif}_pct' for motif in all_3mers]
    df_motifs = df[pct_cols].copy()
    
    # 执行聚类分析
    clustering_results = perform_clustering_analysis(df_motifs, output_dir)
    
    # 生成摘要报告
    print("\n生成分析报告...")
    report_file = output_dir / "analysis_report.txt"
    with open(report_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("m6A Complete 3-mer Motif Distribution and Clustering Analysis\n")
        f.write("="*80 + "\n\n")
        
        f.write("Analysis Overview:\n")
        f.write(f"  Total species analyzed: {len(results)}\n")
        f.write(f"  Total 3-mer motifs tracked: {len(all_3mers)}\n")
        f.write(f"  Clustering methods applied: PCA, t-SNE, K-means, Hierarchical\n\n")
        
        f.write("="*80 + "\n\n")
        
        for result in results:
            f.write(f"Species: {result['Species']}\n")
            f.write(f"{'-'*80}\n")
            f.write(f"  Total m6A sites: {result['Total_m6A_sites']:,}\n")
            f.write(f"  Valid motifs: {result['Valid_motifs']:,}\n")
            f.write(f"  Positive strand: {result['Positive_strand']:,}\n")
            f.write(f"  Negative strand: {result['Negative_strand']:,}\n\n")
            
            # Top 10 motifs
            motif_data = [(m, result[f'{m}_count'], result[f'{m}_pct']) 
                          for m in all_3mers]
            motif_data.sort(key=lambda x: x[1], reverse=True)
            
            f.write(f"  Top 10 Motifs:\n")
            for i, (motif, count, pct) in enumerate(motif_data[:10], 1):
                f.write(f"    {i:2d}. {motif}: {count:6,} ({pct:6.2f}%)\n")
            
            f.write("\n" + "="*80 + "\n\n")
    
    print(f"已保存分析报告: {report_file}")
    print(f"\n分析完成！所有结果保存在: {output_dir}")


if __name__ == "__main__":
    main()
