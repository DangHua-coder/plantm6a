#!/usr/bin/env python3
"""
保守m6A位点分析模块（完整版）

在全转录本（5'UTR + CDS + 3'UTR）上做保守m6A位点分析
"""

import os
import sys
import subprocess
import logging
from pathlib import Path
from collections import defaultdict

try:
    import parasail
    USE_PARASAIL = True
except ImportError:
    try:
        from Bio import pairwise2
        USE_PARASAIL = False
    except ImportError:
        USE_PARASAIL = None

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(message)s',
    datefmt='%H:%M:%S',
)
log = logging.getLogger(__name__)

# 比对参数
TOLERANCE = 1
MATCH = 2
MISMATCH = -1
GAP_OPEN = -2
GAP_EXT = -1


class TranscriptModel:
    """存储一个转录本的完整结构"""
    __slots__ = ['tid', 'chrom', 'strand', 'exons', 'cds_intervals']

    def __init__(self, tid, chrom, strand):
        self.tid = tid
        self.chrom = chrom
        self.strand = strand
        self.exons = []
        self.cds_intervals = []

    def finalize(self):
        """排序并去重"""
        self.exons = sorted(set(self.exons), key=lambda x: x[0])
        self.cds_intervals = sorted(set(self.cds_intervals), key=lambda x: x[0])

    @property
    def tx_start(self):
        return self.exons[0][0] if self.exons else None

    @property
    def tx_end(self):
        return self.exons[-1][1] if self.exons else None

    @property
    def cds_start(self):
        return self.cds_intervals[0][0] if self.cds_intervals else None

    @property
    def cds_end(self):
        return self.cds_intervals[-1][1] if self.cds_intervals else None


def parse_annotation(annot_path, fmt='gtf'):
    """解析GTF或GFF3"""
    models = {}

    def get_tid_gtf(attrs):
        for tok in attrs.split(';'):
            tok = tok.strip()
            if tok.startswith('transcript_id'):
                return tok.split('"')[1] if '"' in tok else tok.split()[-1]
        return None

    def get_tid_gff3(attrs):
        for tok in attrs.split(';'):
            tok = tok.strip()
            if tok.startswith('Parent='):
                return tok[7:].split(',')[0]
        for tok in attrs.split(';'):
            tok = tok.strip()
            if tok.startswith('transcript_id='):
                return tok[14:]
        return None

    with open(annot_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            feat = parts[2]
            if feat not in ('exon', 'CDS'):
                continue

            chrom = parts[0]
            start = int(parts[3]) - 1
            end = int(parts[4])
            strand = parts[6]
            attrs = parts[8]

            tid = get_tid_gtf(attrs) if fmt == 'gtf' else get_tid_gff3(attrs)
            if not tid:
                continue

            if tid not in models:
                models[tid] = TranscriptModel(tid, chrom, strand)

            if feat == 'exon':
                models[tid].exons.append((start, end))
            elif feat == 'CDS':
                models[tid].cds_intervals.append((start, end))

    for m in models.values():
        m.finalize()

    models = {tid: m for tid, m in models.items() if m.exons}
    log.info(f"    解析得到 {len(models)} 个有exon注释的transcript")
    return models


def build_position_index(models):
    """建立染色体位置索引"""
    idx = defaultdict(list)
    for tid, m in models.items():
        if m.tx_start is None:
            continue
        idx[m.chrom].append((m.tx_start, m.tx_end, tid, m.strand))
    for chrom in idx:
        idx[chrom].sort(key=lambda x: x[0])
    return idx


def reverse_complement(seq):
    """反向互补"""
    comp = str.maketrans('ACGTN', 'TGCAN')
    return seq.translate(comp)[::-1]


def get_transcript_sequence(model, fasta_path):
    """提取全转录本序列"""
    exons = model.exons
    pieces = []
    for s, e in exons:
        region = f"{model.chrom}:{s+1}-{e}"
        cmd = ['samtools', 'faidx', fasta_path, region]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            seq = ''.join(result.stdout.split('\n')[1:]).upper()
            pieces.append(seq)
        except subprocess.CalledProcessError:
            return None

    tx_seq = ''.join(pieces)
    if model.strand == '-':
        tx_seq = reverse_complement(tx_seq)
    return tx_seq


def cds_to_tx_pos_plus(model):
    """正链：计算CDS在转录本中的位置"""
    if not model.cds_intervals:
        return None, None
    cds_genomic_start = model.cds_intervals[0][0]
    cds_genomic_end = model.cds_intervals[-1][1]

    offset = 0
    cds_tx_start = None
    cds_tx_end = None
    for s, e in model.exons:
        if s <= cds_genomic_start < e:
            cds_tx_start = offset + (cds_genomic_start - s)
        if s < cds_genomic_end <= e:
            cds_tx_end = offset + (cds_genomic_end - s)
            break
        offset += (e - s)
    return cds_tx_start, cds_tx_end


def cds_to_tx_pos_minus(model):
    """负链：计算CDS在转录本中的位置"""
    if not model.cds_intervals:
        return None, None
    cds_genomic_start = model.cds_intervals[0][0]
    cds_genomic_end = model.cds_intervals[-1][1]

    offset = 0
    cds_tx_start = None
    cds_tx_end = None
    for s, e in reversed(model.exons):
        if s < cds_genomic_end <= e:
            cds_tx_start = offset + (e - cds_genomic_end)
        if s <= cds_genomic_start < e:
            cds_tx_end = offset + (e - 1 - cds_genomic_start) + 1
            break
        offset += (e - s)
    return cds_tx_start, cds_tx_end


def classify_region_plus(tx_pos, model):
    """正链：判断区域"""
    cds_start, cds_end = cds_to_tx_pos_plus(model)
    if cds_start is None:
        return 'noncoding'
    if tx_pos < cds_start:
        return '5utr'
    elif tx_pos < cds_end:
        return 'cds'
    else:
        return '3utr'


def classify_region_minus(tx_pos, model):
    """负链：判断区域"""
    cds_start, cds_end = cds_to_tx_pos_minus(model)
    if cds_start is None:
        return 'noncoding'
    if tx_pos < cds_start:
        return '5utr'
    elif tx_pos < cds_end:
        return 'cds'
    else:
        return '3utr'


def genomic_to_transcript_pos(genomic_pos, model):
    """基因组坐标映射到转录本坐标"""
    exons = model.exons

    if model.strand == '+':
        offset = 0
        for s, e in exons:
            if s <= genomic_pos < e:
                tx_pos = offset + (genomic_pos - s)
                region = classify_region_plus(tx_pos, model)
                return tx_pos, region
            offset += (e - s)
    else:
        offset = 0
        for s, e in reversed(exons):
            if s <= genomic_pos < e:
                tx_pos = offset + (e - 1 - genomic_pos)
                region = classify_region_minus(tx_pos, model)
                return tx_pos, region
            offset += (e - s)

    return None, None


def load_m6a_sites(bed_file):
    """加载m6A位点"""
    sites = defaultdict(list)
    with open(bed_file) as f:
        for line in f:
            if line.startswith('#') or line.startswith('track'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            pos = (start + end) // 2
            score = float(parts[4]) if len(parts) > 4 else 1.0
            sites[chrom].append((pos, score))
    return sites


def map_m6a_to_transcripts(m6a_sites, position_index, models):
    """将m6A映射到转录本"""
    result = defaultdict(list)
    
    for chrom, sites in m6a_sites.items():
        if chrom not in position_index:
            continue
        
        tx_list = position_index[chrom]
        
        for pos, score in sites:
            for tx_start, tx_end, tid, strand in tx_list:
                if tx_start <= pos < tx_end:
                    model = models[tid]
                    tx_pos, region = genomic_to_transcript_pos(pos, model)
                    if tx_pos is not None:
                        result[tid].append((tx_pos, region, score))

    return result

def sw_align(seqA, seqB):
    """Smith-Waterman比对"""
    if USE_PARASAIL is None:
        raise ImportError("需要安装 parasail 或 biopython")
    
    if USE_PARASAIL:
        matrix = parasail.matrix_create("ACGTN", MATCH, MISMATCH)
        res = parasail.sw_trace_striped_16(
            seqA, seqB, int(abs(GAP_OPEN)), int(abs(GAP_EXT)), matrix
        )
        tb = res.traceback
        return tb.query, tb.ref, res.score
    else:
        from Bio import pairwise2
        alns = pairwise2.align.localms(
            seqA, seqB, MATCH, MISMATCH, GAP_OPEN, GAP_EXT,
            one_alignment_only=True
        )
        if not alns:
            return None, None, 0
        a = alns[0]
        return str(a.seqA), str(a.seqB), a.score


def build_alignment_map(aligned_seq):
    """构建比对映射"""
    cds_to_col = {}
    col_to_cds = {}
    pos = 0
    for col, ch in enumerate(aligned_seq):
        if ch != '-':
            cds_to_col[pos] = col
            col_to_cds[col] = pos
            pos += 1
    return cds_to_col, col_to_cds


def find_conserved_sites(tidA, tidB, seqA, seqB, m6a_A, m6a_B, tolerance=TOLERANCE):
    """找保守m6A位点"""
    if not seqA or not seqB or not m6a_A or not m6a_B:
        return []

    alnA, alnB, score = sw_align(seqA, seqB)
    if alnA is None:
        return []

    tx_to_colA, _ = build_alignment_map(alnA)
    tx_to_colB, _ = build_alignment_map(alnB)

    m6a_dict_A = {}
    for tx_pos, region, s in m6a_A:
        if tx_pos not in m6a_dict_A or s > m6a_dict_A[tx_pos][1]:
            m6a_dict_A[tx_pos] = (region, s)

    m6a_dict_B = {}
    for tx_pos, region, s in m6a_B:
        if tx_pos not in m6a_dict_B or s > m6a_dict_B[tx_pos][1]:
            m6a_dict_B[tx_pos] = (region, s)

    conserved = []
    for posA, (regionA, scoreA) in m6a_dict_A.items():
        if posA not in tx_to_colA:
            continue
        colA = tx_to_colA[posA]

        for posB, (regionB, scoreB) in m6a_dict_B.items():
            if posB not in tx_to_colB:
                continue
            colB = tx_to_colB[posB]

            if abs(colA - colB) <= tolerance:
                conserved.append({
                    'tidA': tidA, 'tidB': tidB,
                    'posA': posA, 'regionA': regionA, 'scoreA': scoreA,
                    'posB': posB, 'regionB': regionB, 'scoreB': scoreB,
                    'colA': colA, 'colB': colB,
                    'col_diff': abs(colA - colB),
                    'aln_score': score,
                    'aln_len': len(alnA),
                    'ctx_A': alnA[max(0, colA-10):colA+11],
                    'ctx_B': alnB[max(0, colB-10):colB+11],
                })
    return conserved


def analyze_species_pair(
    spA, spB,
    genome_A, genome_B,
    annot_A, annot_B,
    fmt_A, fmt_B,
    m6a_A_file, m6a_B_file,
    ortholog_file,
    output_file,
    orthologs_only=True,
    verbose=True
):
    """分析两个物种的保守m6A位点"""
    if verbose:
        log.info(f"{'='*60}")
        log.info(f"处理 {spA} vs {spB}")
    
    if verbose:
        log.info(f"  [{spA}] 解析注释...")
    modelsA = parse_annotation(annot_A, fmt_A)
    
    if verbose:
        log.info(f"  [{spB}] 解析注释...")
    modelsB = parse_annotation(annot_B, fmt_B)

    idxA = build_position_index(modelsA)
    idxB = build_position_index(modelsB)

    if verbose:
        log.info(f"  [{spA}] 映射m6A...")
    m6a_raw_A = load_m6a_sites(m6a_A_file)
    gene_m6a_A = map_m6a_to_transcripts(m6a_raw_A, idxA, modelsA)

    if verbose:
        log.info(f"  [{spB}] 映射m6A...")
    m6a_raw_B = load_m6a_sites(m6a_B_file)
    gene_m6a_B = map_m6a_to_transcripts(m6a_raw_B, idxB, modelsB)

    if verbose:
        log.info(f"  [{spA}] m6A转录本: {len(gene_m6a_A)}")
        log.info(f"  [{spB}] m6A转录本: {len(gene_m6a_B)}")

    ortho_pairs = []
    with open(ortholog_file) as f:
        f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            og, gA, gB, tp = parts[0], parts[1], parts[2], parts[3]
            if orthologs_only and tp != 'ortholog':
                continue
            ortho_pairs.append((og, gA, gB))

    if verbose:
        log.info(f"  同源基因对: {len(ortho_pairs)}")

    n_conserved = 0
    n_with_m6a = 0

    with open(output_file, 'w') as out:
        out.write('\t'.join([
            'OG_ID',
            f'{spA}_transcript', f'{spB}_transcript',
            f'{spA}_tx_pos', f'{spB}_tx_pos',
            f'{spA}_region', f'{spB}_region',
            f'{spA}_m6a_score', f'{spB}_m6a_score',
            f'{spA}_aln_col', f'{spB}_aln_col',
            'col_diff', 'aln_score', 'aln_len',
            f'{spA}_context', f'{spB}_context',
        ]) + '\n')

        for og, gA, gB in ortho_pairs:
            m6a_A = gene_m6a_A.get(gA, [])
            m6a_B = gene_m6a_B.get(gB, [])
            if not m6a_A or not m6a_B:
                continue

            n_with_m6a += 1
            seqA = get_transcript_sequence(modelsA[gA], genome_A) if gA in modelsA else None
            seqB = get_transcript_sequence(modelsB[gB], genome_B) if gB in modelsB else None

            if not seqA or not seqB:
                continue

            conserved = find_conserved_sites(gA, gB, seqA, seqB, m6a_A, m6a_B)
            for site in conserved:
                out.write('\t'.join(map(str, [
                    og,
                    site['tidA'], site['tidB'],
                    site['posA'], site['posB'],
                    site['regionA'], site['regionB'],
                    f"{site['scoreA']:.3f}", f"{site['scoreB']:.3f}",
                    site['colA'], site['colB'],
                    site['col_diff'],
                    f"{site['aln_score']:.1f}",
                    site['aln_len'],
                    site['ctx_A'], site['ctx_B'],
                ])) + '\n')
                n_conserved += 1

    if verbose:
        log.info(f"  结果: {n_with_m6a} 对含m6A，{n_conserved} 个保守位点")
        log.info(f"  输出: {output_file}")
    
    return n_conserved


__all__ = [
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
]
