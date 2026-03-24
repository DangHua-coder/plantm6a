#!/usr/bin/env python3
"""
保守m6A位点汇总分析（完整版）
"""

import os
import sys
from pathlib import Path
from collections import defaultdict, Counter


def summarize_conserved_m6a(conserved_dir, output_dir, verbose=True):
    """汇总保守m6A分析结果"""
    conserved_dir = Path(conserved_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    all_files = list(conserved_dir.glob("*_conserved_m6A.tsv"))
    
    if verbose:
        print(f"找到 {len(all_files)} 个物种对结果文件")

    if not all_files:
        print("错误: 未找到任何结果文件")
        return {}

    # 1. 每对物种的保守位点数量
    summary = []
    for f in sorted(all_files):
        name = f.stem.replace('_conserved_m6A', '')
        parts = name.split('_vs_')
        if len(parts) != 2:
            continue
        spA, spB = parts[0], parts[1]
        with open(f) as fh:
            n = sum(1 for line in fh) - 1
        summary.append((spA, spB, n))

    summary_file = output_dir / 'pair_summary.tsv'
    with open(summary_file, 'w') as out:
        out.write("SpeciesA\tSpeciesB\tConserved_m6A_sites\n")
        for spA, spB, n in sorted(summary):
            out.write(f"{spA}\t{spB}\t{n}\n")
    
    if verbose:
        print(f"汇总统计 → {summary_file}")

    # 2. 详细位点表
    all_sites = []
    for f in sorted(all_files):
        name = f.stem.replace('_conserved_m6A', '')
        parts = name.split('_vs_')
        if len(parts) != 2:
            continue
        spA, spB = parts[0], parts[1]
        
        with open(f) as fh:
            fh.readline()
            for line in fh:
                parts_line = line.rstrip('\n').split('\t')
                if len(parts_line) < 12:
                    continue
                
                all_sites.append({
                    'spA': spA,
                    'spB': spB,
                    'og': parts_line[0],
                    'tA': parts_line[1],
                    'tB': parts_line[2],
                    'posA': int(parts_line[3]),
                    'posB': int(parts_line[4]),
                    'regionA': parts_line[5],
                    'regionB': parts_line[6],
                    'scoreA': float(parts_line[7]),
                    'scoreB': float(parts_line[8]),
                    'colA': int(parts_line[9]),
                    'colB': int(parts_line[10]),
                    'diff': int(parts_line[11]),
                    'aln_score': float(parts_line[12]) if len(parts_line) > 12 else 0.0,
                    'ctx_A': parts_line[14] if len(parts_line) > 14 else '',
                    'ctx_B': parts_line[15] if len(parts_line) > 15 else '',
                })

    if verbose:
        print(f"总保守位点数（含重复）: {len(all_sites)}")

    # 3. 区域分布统计
    region_counter = Counter()
    for s in all_sites:
        key = f"{s['regionA']}-{s['regionB']}"
        region_counter[key] += 1

    region_file = output_dir / 'region_distribution.tsv'
    with open(region_file, 'w') as out:
        out.write("regionA-regionB\tcount\tpercent\n")
        total = len(all_sites)
        for combo, cnt in sorted(region_counter.items(), key=lambda x: -x[1]):
            pct = cnt / total * 100 if total > 0 else 0
            out.write(f"{combo}\t{cnt}\t{pct:.2f}%\n")
    
    if verbose:
        print(f"区域分布统计 → {region_file}")

    # 4. 按区域分解的物种对统计
    pair_region_breakdown = []
    for f in sorted(all_files):
        name = f.stem.replace('_conserved_m6A', '')
        parts = name.split('_vs_')
        if len(parts) != 2:
            continue
        spA, spB = parts[0], parts[1]
        
        region_count = Counter()
        with open(f) as fh:
            fh.readline()
            for line in fh:
                parts_line = line.rstrip('\n').split('\t')
                if len(parts_line) < 7:
                    continue
                combo = f"{parts_line[5]}-{parts_line[6]}"
                region_count[combo] += 1
        
        for combo, cnt in region_count.items():
            pair_region_breakdown.append({
                'spA': spA,
                'spB': spB,
                'region_combo': combo,
                'count': cnt
            })
    
    pair_region_file = output_dir / 'pair_region_breakdown.tsv'
    with open(pair_region_file, 'w') as out:
        out.write("spA\tspB\tregion_combo\tcount\n")
        for item in pair_region_breakdown:
            out.write(f"{item['spA']}\t{item['spB']}\t{item['region_combo']}\t{item['count']}\n")
    
    if verbose:
        print(f"物种对区域分解 → {pair_region_file}")

    # 5. 多物种保守OG统计
    og_species_pairs = defaultdict(set)
    for s in all_sites:
        og_species_pairs[s['og']].add((s['spA'], s['spB']))

    multi_conserved = {og: pairs for og, pairs in og_species_pairs.items() if len(pairs) >= 2}
    
    if verbose:
        print(f"在≥2对物种中有保守m6A的OG数: {len(multi_conserved)}")

    multi_file = output_dir / 'multi_species_conserved_OGs.tsv'
    with open(multi_file, 'w') as out:
        out.write("OG_ID\tN_species_pairs\tSpecies_pairs\n")
        for og, pairs in sorted(multi_conserved.items(), key=lambda x: -len(x[1])):
            pair_str = ';'.join(f"{a}:{b}" for a, b in sorted(pairs))
            out.write(f"{og}\t{len(pairs)}\t{pair_str}\n")
    
    if verbose:
        print(f"多物种保守OG表 → {multi_file}")
    
    return {
        'summary': summary,
        'all_sites': all_sites,
        'region_distribution': dict(region_counter),
        'multi_conserved_ogs': multi_conserved,
        'pair_region_breakdown': pair_region_breakdown
    }


def main():
    """命令行接口"""
    import argparse
    parser = argparse.ArgumentParser(description="汇总保守m6A分析结果")
    parser.add_argument('--conserved_dir', required=True)
    parser.add_argument('--output_dir', required=True)
    args = parser.parse_args()
    
    summarize_conserved_m6a(args.conserved_dir, args.output_dir)


if __name__ == '__main__':
    main()


__all__ = ['summarize_conserved_m6a']
