#!/usr/bin/env python3
"""
Command line interface for plantm6a
"""

import sys
import argparse
from plantm6a.analysis.statistics import analyze
from plantm6a.analysis.batch import batch_analyze


def main():
    parser = argparse.ArgumentParser(
        description='PlantM6A: A toolkit for plant m6A analysis'
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # stats命令：单个物种统计
    stats_parser = subparsers.add_parser('stats', help='Statistical analysis for single species')
    stats_parser.add_argument('--genome', required=True, help='Genome FASTA file')
    stats_parser.add_argument('--gtf', required=True, help='Annotation GTF/GFF file')
    stats_parser.add_argument('--bed', required=True, help='m6A sites BED file')
    
    # batch命令：批量分析
    batch_parser = subparsers.add_parser('batch', help='Batch analysis for multiple species')
    batch_parser.add_argument('--config', required=True, help='Configuration YAML file')
    batch_parser.add_argument('--output', help='Output TSV file (optional)')
    
    args = parser.parse_args()
    
    if args.command == 'stats':
        analyze(args.genome, args.gtf, args.bed, verbose=True)
    elif args.command == 'batch':
        batch_analyze(args.config, args.output)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()
