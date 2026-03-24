#!/bin/bash
# =============================================================================
# OrthoFinder Pairwise Ortholog Analysis Pipeline
# 使用PlantM6A conservation模块进行保守性分析
# =============================================================================

set -euo pipefail

# ====== 参数配置 ======
INPUT_CDS_DIR="./input_cds"          # CDS核苷酸序列目录
INPUT_PEP_DIR="./input_proteins"     # 蛋白序列目录（翻译后）
OUTPUT_DIR="./orthofinder_out"       # OrthoFinder输出目录
PAIRWISE_DIR="./pairwise_tables"     # 两两物种基因对表输出目录
THREADS=${1:-40}                      # 线程数（默认40）
# =====================

mkdir -p "$INPUT_PEP_DIR" "$PAIRWISE_DIR"

# =============================================================================
# STEP 0: 翻译CDS为蛋白序列
# =============================================================================
echo "========================================"
echo "STEP 0: 翻译CDS为蛋白序列"
echo "========================================"

if ls "${INPUT_CDS_DIR}"/*.fa 1>/dev/null 2>&1; then
    for cds_fa in "${INPUT_CDS_DIR}"/*.fa; do
        species=$(basename "$cds_fa" .fa)
        pep_fa="${INPUT_PEP_DIR}/${species}.fa"
        
        if [ ! -f "$pep_fa" ]; then
            echo "  翻译: $species"
            python3 -m plantm6a.analysis.conservation.translate "$cds_fa" "$pep_fa"
        else
            echo "  跳过: $species (已存在)"
        fi
    done
else
    echo "  警告: ${INPUT_CDS_DIR} 为空或不存在"
    echo "  请将CDS fasta文件放入此目录"
    exit 1
fi

# =============================================================================
# STEP 1: 运行OrthoFinder
# =============================================================================
echo "========================================"
echo "STEP 1: 运行OrthoFinder"
echo "========================================"

if ! command -v orthofinder &>/dev/null; then
    echo "ERROR: 未找到orthofinder，请先安装："
    echo "  conda install -c bioconda orthofinder"
    exit 1
fi

orthofinder \
    -f "$INPUT_PEP_DIR" \
    -o "$OUTPUT_DIR" \
    -t "$THREADS" \
    -a "$THREADS" \
    -S diamond

echo "OrthoFinder 完成！"

# =============================================================================
# STEP 2: 定位OrthoFinder结果目录
# =============================================================================
echo "========================================"
echo "STEP 2: 定位结果文件"
echo "========================================"

RESULTS_DIR=$(find "$OUTPUT_DIR" -name "Orthogroups" -type d | head -1 | xargs dirname)
if [ -z "$RESULTS_DIR" ]; then
    echo "ERROR: 未找到OrthoFinder结果目录"
    exit 1
fi
echo "  结果目录: $RESULTS_DIR"

ORTHOGROUPS_TSV="${RESULTS_DIR}/Orthogroups/Orthogroups.tsv"
ORTHOLOGS_DIR="${RESULTS_DIR}/Orthologues"

# =============================================================================
# STEP 3: 提取两两物种基因对
# =============================================================================
echo "========================================"
echo "STEP 3: 提取两两物种基因对"
echo "========================================"

python3 -m plantm6a.analysis.conservation.extract_orthologs \
    --orthogroups "$ORTHOGROUPS_TSV" \
    --orthologues_dir "$ORTHOLOGS_DIR" \
    --output_dir "$PAIRWISE_DIR"

echo "========================================"
echo "全部完成！结果在: $PAIRWISE_DIR"
echo "========================================"
