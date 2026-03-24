#!/usr/bin/env python3
"""
从OrthoFinder结果中提取两两物种之间的同源基因对
"""

import os
import re
import sys
import argparse
import itertools
from pathlib import Path
from collections import defaultdict


def parse_ortholog_files(orthologues_dir, verbose=True):
    """
    解析OrthoFinder的Orthologues目录（直系同源）
    
    返回:
        dict: (spA, spB) -> list of (og_id, geneA, geneB)
    """
    pairs = defaultdict(list)
    
    ortho_path = Path(orthologues_dir)
    if not ortho_path.exists():
        if verbose:
            print(f"  警告: Orthologues目录不存在: {orthologues_dir}")
        return pairs

    tsv_files = list(ortho_path.rglob("*__v__*.tsv"))
    if verbose:
        print(f"  找到 {len(tsv_files)} 个直系同源文件")

    for tsv_file in tsv_files:
        fname = tsv_file.stem
        match = re.search(r'(.+)__v__(.+)', fname)
        if not match:
            continue
        spA = match.group(1).strip()
        spB = match.group(2).strip()

        with open(tsv_file) as f:
            header = f.readline().rstrip('\n').split('\t')
            if len(header) < 3:
                continue
            
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 3:
                    continue
                og_id = parts[0]
                genes_a = [g.strip() for g in parts[1].split(',') if g.strip()]
                genes_b = [g.strip() for g in parts[2].split(',') if g.strip()]
                
                if not genes_a or not genes_b:
                    continue
                
                for ga, gb in itertools.product(genes_a, genes_b):
                    pairs[(spA, spB)].append((og_id, ga, gb))

    return pairs


def parse_orthogroups(orthogroups_tsv, verbose=True):
    """
    解析Orthogroups.tsv
    
    返回:
        species_list: 物种列表
        og_dict: og_id -> {species: [genes]}
    """
    og_dict = {}
    species_list = []
    
    with open(orthogroups_tsv) as f:
        header = f.readline().rstrip('\n').split('\t')
        species_list = header[1:]
        
        for line in f:
            parts = line.rstrip('\n').split('\t')
            og_id = parts[0]
            genes_by_sp = {}
            for i, sp in enumerate(species_list):
                col = parts[i+1] if i+1 < len(parts) else ''
                genes = [g.strip() for g in col.split(',') if g.strip()]
                if genes:
                    genes_by_sp[sp] = genes
            og_dict[og_id] = genes_by_sp
    
    if verbose:
        print(f"  读取到 {len(og_dict)} 个Orthogroups，{len(species_list)} 个物种")
    
    return species_list, og_dict


def write_pairwise_tables(species_list, og_dict, ortholog_pairs, output_dir, verbose=True):
    """
    为每对物种写出TSV文件
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    all_species = sorted(species_list)
    stats = {}

    for spA, spB in itertools.combinations(all_species, 2):
        outfile = output_dir / f"{spA}_vs_{spB}.tsv"
        rows = []

        # 直系同源
        for key in [(spA, spB), (spB, spA)]:
            for og_id, gA, gB in ortholog_pairs.get(key, []):
                if key == (spA, spB):
                    rows.append((og_id, gA, gB, 'ortholog'))
                else:
                    rows.append((og_id, gB, gA, 'ortholog'))

        # 同OG内的跨物种对
        for og_id, genes_by_sp in og_dict.items():
            genes_a = genes_by_sp.get(spA, [])
            genes_b = genes_by_sp.get(spB, [])
            if not genes_a or not genes_b:
                continue
            for ga, gb in itertools.product(genes_a, genes_b):
                rows.append((og_id, ga, gb, 'inparalog_og'))

        # 去重
        seen = set()
        deduped = []
        for row in rows:
            key = (row[1], row[2])
            if key not in seen:
                seen.add(key)
                deduped.append(row)

        with open(outfile, 'w') as f:
            f.write(f"OG_ID\t{spA}_gene\t{spB}_gene\ttype\n")
            for row in deduped:
                f.write('\t'.join(row) + '\n')
        
        n_ortho = sum(1 for r in deduped if r[3] == 'ortholog')
        n_inpara = sum(1 for r in deduped if r[3] == 'inparalog_og')
        stats[(spA, spB)] = {'orthologs': n_ortho, 'inparalogs': n_inpara, 'total': len(deduped)}
        
        if verbose:
            print(f"  {spA} vs {spB}: {n_ortho} orthologs, {n_inpara} inparalogs → {outfile.name}")

    # 汇总统计
    summary_file = output_dir / "summary_statistics.tsv"
    with open(summary_file, 'w') as f:
        f.write("SpeciesA\tSpeciesB\tOrthologs\tInParalogs\tTotal_pairs\n")
        for (spA, spB), s in sorted(stats.items()):
            f.write(f"{spA}\t{spB}\t{s['orthologs']}\t{s['inparalogs']}\t{s['total']}\n")
    
    if verbose:
        print(f"\n  汇总统计: {summary_file}")
    
    return stats


def extract_pairwise_orthologs(orthogroups_tsv, orthologues_dir, output_dir, verbose=True):
    """
    主函数：提取两两物种同源基因对
    
    参数:
        orthogroups_tsv: Orthogroups.tsv文件路径
        orthologues_dir: Orthologues/目录路径
        output_dir: 输出目录
        verbose: 是否显示详细信息
    
    返回:
        stats: 统计信息字典
    """
    if verbose:
        print("=== 解析Orthogroups.tsv ===")
    species_list, og_dict = parse_orthogroups(orthogroups_tsv, verbose)

    if verbose:
        print("\n=== 解析直系同源文件 ===")
    ortholog_pairs = parse_ortholog_files(orthologues_dir, verbose)

    if verbose:
        print("\n=== 生成两两物种对应表 ===")
    stats = write_pairwise_tables(species_list, og_dict, ortholog_pairs, output_dir, verbose)
    
    return stats


def main():
    """命令行接口"""
    parser = argparse.ArgumentParser(description="从OrthoFinder结果提取两两物种同源基因对")
    parser.add_argument('--orthogroups', required=True,
                        help='Orthogroups.tsv 文件路径')
    parser.add_argument('--orthologues_dir', required=True,
                        help='OrthoFinder Orthologues/ 目录路径')
    parser.add_argument('--output_dir', required=True,
                        help='输出目录')
    args = parser.parse_args()

    extract_pairwise_orthologs(
        args.orthogroups,
        args.orthologues_dir,
        args.output_dir,
        verbose=True
    )
    
    print("\n✓ 完成！")


if __name__ == "__main__":
    main()


__all__ = ['extract_pairwise_orthologs', 'parse_orthogroups', 'parse_ortholog_files']
