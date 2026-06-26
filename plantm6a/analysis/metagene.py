#!/usr/bin/env python3
"""Metagene profile binning and plotting utilities."""

import math
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

REGION_TYPES = [
    "UTR5",
    "CDS",
    "UTR3",
    "exonFirst",
    "exonInternal",
    "exonLast",
    "mRNAIntron",
    "ncRNAIntron",
]


def mean(data):
    if not data:
        return 0
    return sum(data) / len(data)


def _parse_int_list(value):
    return [int(x) for x in value.rstrip(",").split(",") if x]


def calculate_weights(annot_file):
    mRNA_num = 0
    ncRNA_num = 0
    UTR5_sizes = []
    CDS_sizes = []
    UTR3_sizes = []
    exon_first_sizes = []
    exon_internal_sizes = []
    exon_last_sizes = []

    with open(annot_file) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue

            fields = line.split("\t")
            tx_start = int(fields[1])
            tx_strand = fields[5]
            cds_start = int(fields[6])
            cds_end = int(fields[7])
            exon_lengths = _parse_int_list(fields[10])
            exon_starts = _parse_int_list(fields[11])
            exon_len_sum = sum(exon_lengths)

            if cds_start != cds_end:
                mRNA_num += 1
                left_size_sum = 0
                left_UTR_size = 0
                for exon_len, exon_start in zip(exon_lengths, exon_starts):
                    start = tx_start + exon_start
                    end = start + exon_len
                    if cds_start >= start and cds_start < end:
                        left_UTR_size = left_size_sum + cds_start - start
                        break
                    left_size_sum += exon_len

                right_size_sum = 0
                right_UTR_size = 0
                for i in range(len(exon_lengths) - 1, -1, -1):
                    start = tx_start + exon_starts[i]
                    end = start + exon_lengths[i]
                    if cds_end > start and cds_end <= end:
                        right_UTR_size = right_size_sum + end - cds_end
                        break
                    right_size_sum += exon_lengths[i]

                if tx_strand == "+":
                    UTR5_sizes.append(left_UTR_size)
                    UTR3_sizes.append(right_UTR_size)
                elif tx_strand == "-":
                    UTR5_sizes.append(right_UTR_size)
                    UTR3_sizes.append(left_UTR_size)
                else:
                    print(f"Line {line_num} had unknown strand!", file=sys.stderr)
                CDS_sizes.append(exon_len_sum - left_UTR_size - right_UTR_size)
            else:
                ncRNA_num += 1
                exon_num = len(exon_lengths)
                if exon_num == 1:
                    exon_first_sizes.append(exon_lengths[0])
                elif exon_num == 2:
                    if tx_strand == "+":
                        exon_first_sizes.append(exon_lengths[0])
                        exon_last_sizes.append(exon_lengths[1])
                    elif tx_strand == "-":
                        exon_first_sizes.append(exon_lengths[1])
                        exon_last_sizes.append(exon_lengths[0])
                    else:
                        print(f"Line {line_num} had unknown strand!", file=sys.stderr)
                else:
                    if tx_strand == "+":
                        exon_first_sizes.append(exon_lengths[0])
                        exon_last_sizes.append(exon_lengths[-1])
                        exon_internal_sizes.append(exon_len_sum - exon_lengths[0] - exon_lengths[-1])
                    elif tx_strand == "-":
                        exon_first_sizes.append(exon_lengths[-1])
                        exon_last_sizes.append(exon_lengths[0])
                        exon_internal_sizes.append(exon_len_sum - exon_lengths[0] - exon_lengths[-1])
                    else:
                        print(f"Line {line_num} had unknown strand!", file=sys.stderr)

    mean_UTR5 = mean(UTR5_sizes)
    mean_CDS = mean(CDS_sizes)
    mean_UTR3 = mean(UTR3_sizes)
    mean_mRNA = mean_UTR5 + mean_CDS + mean_UTR3
    mean_exon_first = mean(exon_first_sizes)
    mean_exon_internal = mean(exon_internal_sizes)
    mean_exon_last = mean(exon_last_sizes)
    mean_ncRNA = mean_exon_first + mean_exon_internal + mean_exon_last
    if mean_ncRNA == 0:
        mean_ncRNA = 1

    weights = {
        "UTR5": round(3 * mean_UTR5 / mean_mRNA, 3) if mean_mRNA > 0 else 1,
        "CDS": round(3 * mean_CDS / mean_mRNA, 3) if mean_mRNA > 0 else 1,
        "UTR3": round(3 * mean_UTR3 / mean_mRNA, 3) if mean_mRNA > 0 else 1,
        "exonFirst": round(3 * mean_exon_first / mean_ncRNA, 3),
        "exonInternal": round(3 * mean_exon_internal / mean_ncRNA, 3),
        "exonLast": round(3 * mean_exon_last / mean_ncRNA, 3),
        "mRNAIntron": 1,
        "ncRNAIntron": 1,
    }

    print(
        f"{mRNA_num}({mean_mRNA:.1f}) Length-scaled 5' UTR ({mean_UTR5:.1f}, {weights['UTR5']}), "
        f"CDS ({mean_CDS:.1f}, {weights['CDS']}) and 3' UTR ({mean_UTR3:.1f}, {weights['UTR3']})",
        file=sys.stderr,
    )
    print(
        f"{ncRNA_num}({mean_ncRNA:.1f}) Length-scaled First Exon ({mean_exon_first:.1f}, {weights['exonFirst']}), "
        f"Internal Exon ({mean_exon_internal:.1f}, {weights['exonInternal']}) and "
        f"Last Exon ({mean_exon_last:.1f}, {weights['exonLast']})",
        file=sys.stderr,
    )
    return weights, mRNA_num, ncRNA_num


def default_weights():
    return {
        "UTR5": 1,
        "CDS": 1,
        "UTR3": 1,
        "exonFirst": 1,
        "exonInternal": 1,
        "exonLast": 1,
        "mRNAIntron": 1,
        "ncRNAIntron": 1,
    }


def run_intersect(in_file, annot_file, out_file, strand=False):
    tmp_file = f"{out_file}.tmp"
    awk_cmd = "awk -F '\\t' 'BEGIN {OFS=\"\\t\"} {a=int(($2+$3)/2)+1;print $1,a-1,a,$4,$5,$6}'"
    intersect_cmd = f"intersectBed -a - -b {annot_file} -wb -s" if strand else f"intersectBed -a - -b {annot_file} -wb"
    full_cmd = f"{awk_cmd} {in_file} | {intersect_cmd} > {tmp_file}"
    try:
        subprocess.run(full_cmd, shell=True, check=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"Error running intersectBed: {exc}") from exc
    return tmp_file


def count_peak_total(in_file):
    with open(in_file) as f:
        return sum(1 for _ in f)


def count_positions(tmp_file):
    pos_num = defaultdict(int)
    with open(tmp_file) as f:
        for line in f:
            fields = line.strip().split("\t")
            loc = fields[0] + fields[2] + fields[5]
            pos_num[loc] += 1
    return pos_num


class RegionBinCalculator:
    def __init__(self, bin_sum, pos_num, peak_sum=0, rpm=False, loc_count=False, loc_pct=False):
        self.bin_sum = bin_sum
        self.pos_num = pos_num
        self.peak_sum = peak_sum
        self.rpm = rpm
        self.loc_count = loc_count
        self.loc_pct = loc_pct
        self.region_bin = defaultdict(lambda: defaultdict(float))

    def add_value(self, region, bin_idx, loc, peak_value):
        if self.rpm:
            self.region_bin[region][bin_idx] += 1 * 1000000 / self.peak_sum
        elif self.loc_count:
            if self.loc_pct:
                self.region_bin[region][bin_idx] += 1 / (self.pos_num[loc] * self.peak_sum)
            else:
                self.region_bin[region][bin_idx] += 1 / self.pos_num[loc]
        elif self.loc_pct:
            self.region_bin[region][bin_idx] += peak_value / (self.pos_num[loc] * self.peak_sum)
        else:
            self.region_bin[region][bin_idx] += peak_value / self.pos_num[loc]


def process_intersection(tmp_file, calculator, bin_nums=None, bin_out_file=None):
    bin_out = open(bin_out_file, "w") if bin_out_file else None
    try:
        with open(tmp_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("chrM"):
                    continue
                fields = line.split("\t")
                peak_mid_point_pos = int(fields[2])
                peak_value = float(fields[4])
                loc = fields[0] + fields[2] + fields[5]
                tx_start = int(fields[7])
                tx_strand = fields[11]
                cds_start = int(fields[12])
                cds_end = int(fields[13])
                exon_num = int(fields[15])
                exon_lengths = _parse_int_list(fields[16])
                exon_starts = _parse_int_list(fields[17])

                if cds_start != cds_end:
                    if peak_mid_point_pos <= cds_start:
                        _, bin_idx, region, _ = process_utr5_region(
                            peak_mid_point_pos, tx_start, tx_strand, cds_start,
                            exon_num, exon_lengths, exon_starts, calculator.bin_sum,
                        )
                    elif peak_mid_point_pos > cds_start and peak_mid_point_pos <= cds_end:
                        _, bin_idx, region, _ = process_cds_region(
                            peak_mid_point_pos, tx_start, tx_strand, cds_start, cds_end,
                            exon_num, exon_lengths, exon_starts, calculator.bin_sum,
                        )
                    else:
                        _, bin_idx, region, _ = process_utr3_region(
                            peak_mid_point_pos, tx_start, tx_strand, cds_end,
                            exon_num, exon_lengths, exon_starts, calculator.bin_sum,
                        )
                else:
                    _, bin_idx, region, _ = process_ncrna_region(
                        peak_mid_point_pos, tx_start, tx_strand,
                        exon_num, exon_lengths, exon_starts, calculator.bin_sum,
                    )

                if region and bin_idx > 0:
                    calculator.add_value(region, bin_idx, loc, peak_value)
                    if bin_nums and bin_idx in bin_nums and bin_out:
                        bin_out.write(f"{bin_idx}\t{region}\t{line}\n")
    finally:
        if bin_out:
            bin_out.close()


def process_utr5_region(peak_pos, tx_start, tx_strand, cds_start, exon_num, exon_lengths, exon_starts, bin_sum):
    pos = 0
    size_sum = 0
    flag_intron = False
    intron_length = 0
    for i in range(len(exon_lengths)):
        start = tx_start + exon_starts[i]
        end = start + exon_lengths[i]
        if exon_num > 1 and i < exon_num - 1:
            start_intron = end
            end_intron = tx_start + exon_starts[i + 1]
            if peak_pos <= end_intron and peak_pos > start_intron:
                pos = peak_pos - start_intron
                intron_length = end_intron - start_intron
                flag_intron = True
                break
        if peak_pos <= end and peak_pos > start:
            pos = peak_pos - start + size_sum
        if end <= cds_start:
            size_sum += exon_lengths[i]
        else:
            size_sum += cds_start - start
            break

    if tx_strand == "+":
        if flag_intron:
            return pos, math.ceil(pos / intron_length * bin_sum), "mRNAIntron", flag_intron
        return pos, math.ceil(pos / size_sum * bin_sum) if size_sum > 0 else 0, "UTR5", flag_intron
    if tx_strand == "-":
        if flag_intron:
            return pos, math.ceil((intron_length - pos + 1) / intron_length * bin_sum), "mRNAIntron", flag_intron
        return pos, math.ceil((size_sum - pos + 1) / size_sum * bin_sum) if size_sum > 0 else 0, "UTR3", flag_intron
    return 0, 0, None, False


def process_cds_region(peak_pos, tx_start, tx_strand, cds_start, cds_end, exon_num, exon_lengths, exon_starts, bin_sum):
    pos = 0
    size_sum = 0
    flag_intron = False
    intron_length = 0
    for i in range(len(exon_lengths)):
        start = tx_start + exon_starts[i]
        end = start + exon_lengths[i]
        if exon_num > 1 and i < exon_num - 1:
            start_intron = end
            end_intron = tx_start + exon_starts[i + 1]
            if peak_pos <= end_intron and peak_pos > start_intron:
                pos = peak_pos - start_intron
                intron_length = end_intron - start_intron
                flag_intron = True
                break
        if end > cds_start and start < cds_end:
            if start <= cds_start:
                if peak_pos <= end and peak_pos > start:
                    pos = peak_pos - cds_start + size_sum
                size_sum += exon_lengths[i] - (cds_start - start)
            elif end >= cds_end:
                if peak_pos <= end and peak_pos > start:
                    pos = peak_pos - start + size_sum
                size_sum += exon_lengths[i] - (end - cds_end)
            else:
                if peak_pos <= end and peak_pos > start:
                    pos = peak_pos - start + size_sum
                size_sum += exon_lengths[i]

    if tx_strand == "+":
        if flag_intron:
            return pos, math.ceil(pos / intron_length * bin_sum), "mRNAIntron", flag_intron
        return pos, math.ceil(pos / size_sum * bin_sum) if size_sum > 0 else 0, "CDS", flag_intron
    if tx_strand == "-":
        if flag_intron:
            return pos, math.ceil((intron_length - pos + 1) / intron_length * bin_sum), "mRNAIntron", flag_intron
        return pos, math.ceil((size_sum - pos + 1) / size_sum * bin_sum) if size_sum > 0 else 0, "CDS", flag_intron
    return 0, 0, None, False


def process_utr3_region(peak_pos, tx_start, tx_strand, cds_end, exon_num, exon_lengths, exon_starts, bin_sum):
    pos = 0
    size_sum = 0
    flag_intron = False
    intron_length = 0
    for i in range(len(exon_lengths) - 1, -1, -1):
        start = tx_start + exon_starts[i]
        end = start + exon_lengths[i]
        if exon_num > 1 and i > 0:
            start_intron = tx_start + exon_starts[i - 1] + exon_lengths[i - 1]
            end_intron = start
            if peak_pos <= end_intron and peak_pos > start_intron:
                pos = end_intron - peak_pos + 1
                intron_length = end_intron - start_intron
                flag_intron = True
                break
        if peak_pos <= end and peak_pos > start:
            pos = end - peak_pos + 1 + size_sum
        if start >= cds_end:
            size_sum += exon_lengths[i]
        else:
            size_sum += end - cds_end
            break

    if tx_strand == "+":
        if flag_intron:
            return pos, math.ceil((intron_length - pos + 1) / intron_length * bin_sum), "mRNAIntron", flag_intron
        return pos, math.ceil((size_sum - pos + 1) / size_sum * bin_sum) if size_sum > 0 else 0, "UTR3", flag_intron
    if tx_strand == "-":
        if flag_intron:
            return pos, math.ceil(pos / intron_length * bin_sum), "mRNAIntron", flag_intron
        return pos, math.ceil(pos / size_sum * bin_sum) if size_sum > 0 else 0, "UTR5", flag_intron
    return 0, 0, None, False


def process_ncrna_region(peak_pos, tx_start, tx_strand, exon_num, exon_lengths, exon_starts, bin_sum):
    pos = 0
    size_sum = 0
    flag_intron = False
    intron_length = 0
    exon_type = None
    for i in range(len(exon_lengths)):
        start = tx_start + exon_starts[i]
        end = start + exon_lengths[i]
        if exon_num > 1 and i < exon_num - 1:
            start_intron = end
            end_intron = tx_start + exon_starts[i + 1]
            if peak_pos <= end_intron and peak_pos > start_intron:
                pos = peak_pos - start_intron
                intron_length = end_intron - start_intron
                flag_intron = True
                break
        if peak_pos <= end and peak_pos > start:
            pos = peak_pos - start + size_sum
            if exon_num == 1:
                exon_type = "internal"
            elif exon_num == 2:
                exon_type = "first" if i == 0 else "last"
            elif i == 0:
                exon_type = "first"
            elif i == exon_num - 1:
                exon_type = "last"
            else:
                exon_type = "internal"
        size_sum += exon_lengths[i]

    if tx_strand == "+":
        if flag_intron:
            return pos, math.ceil(pos / intron_length * bin_sum), "ncRNAIntron", flag_intron
        bin_idx = math.ceil(pos / size_sum * bin_sum) if size_sum > 0 else 0
        if exon_type == "first":
            return pos, bin_idx, "exonFirst", flag_intron
        if exon_type == "last":
            return pos, bin_idx, "exonLast", flag_intron
        return pos, bin_idx, "exonInternal", flag_intron
    if tx_strand == "-":
        if flag_intron:
            return pos, math.ceil((intron_length - pos + 1) / intron_length * bin_sum), "ncRNAIntron", flag_intron
        bin_idx = math.ceil((size_sum - pos + 1) / size_sum * bin_sum) if size_sum > 0 else 0
        if exon_type == "first":
            return pos, bin_idx, "exonLast", flag_intron
        if exon_type == "last":
            return pos, bin_idx, "exonFirst", flag_intron
        return pos, bin_idx, "exonInternal", flag_intron
    return 0, 0, None, False


def write_output(output_file, calculator, weights, bin_sum):
    with open(output_file, "w") as f:
        for region in REGION_TYPES:
            for j in range(1, bin_sum + 1):
                value = calculator.region_bin[region].get(j, 0)
                if region == "UTR5":
                    weighted_pos = weights["UTR5"] * j
                elif region == "CDS":
                    weighted_pos = weights["CDS"] * j + weights["UTR5"] * bin_sum
                elif region == "UTR3":
                    weighted_pos = weights["UTR3"] * j + weights["UTR5"] * bin_sum + weights["CDS"] * bin_sum
                elif region == "exonFirst":
                    weighted_pos = weights["exonFirst"] * j
                elif region == "exonInternal":
                    weighted_pos = weights["exonInternal"] * j + weights["exonFirst"] * bin_sum
                elif region == "exonLast":
                    weighted_pos = weights["exonLast"] * j + weights["exonFirst"] * bin_sum + weights["exonInternal"] * bin_sum
                else:
                    weighted_pos = j
                f.write(f"{region}\t{weighted_pos}\t{value}\n")


def run_region2bin(
    input_bed,
    annotation_bed12,
    output_file,
    keep_tmp=False,
    strand=False,
    len_scale=False,
    loc=False,
    pct=False,
    rpm=False,
    bin_sum=100,
    bin_numbers=None,
    bin_output_file=None,
):
    if rpm and (loc or pct):
        raise ValueError("rpm is mutually exclusive with loc and pct")
    if bin_numbers and not bin_output_file:
        raise ValueError("bin_output_file is required when bin_numbers is set")

    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    weights = calculate_weights(annotation_bed12)[0] if len_scale else default_weights()
    tmp_file = run_intersect(input_bed, annotation_bed12, output_file, strand)
    peak_sum = count_peak_total(input_bed) if pct or rpm else 0
    pos_num = count_positions(tmp_file)
    calculator = RegionBinCalculator(
        bin_sum=bin_sum,
        pos_num=pos_num,
        peak_sum=peak_sum,
        rpm=rpm,
        loc_count=loc,
        loc_pct=pct,
    )
    bin_nums = set(int(x) for x in bin_numbers.split(",")) if isinstance(bin_numbers, str) else bin_numbers
    process_intersection(tmp_file, calculator, bin_nums, bin_output_file)
    write_output(output_file, calculator, weights, bin_sum)
    if not keep_tmp:
        os.remove(tmp_file)
    return output_file


def load_metagene_data(bin_file):
    return pd.read_csv(bin_file, sep="\t", header=None)


def compute_density(x_data, adjust=0.5):
    from scipy.stats import gaussian_kde

    if len(x_data) < 2:
        return np.array([]), np.array([])
    std = np.std(x_data, ddof=1)
    n = len(x_data)
    bw = std * (n ** (-1 / 5)) * adjust
    if bw <= 0:
        bw = 0.1
    kde = gaussian_kde(x_data, bw_method=bw / std if std > 0 else 0.1)
    x_min, x_max = x_data.min(), x_data.max()
    x_range = x_max - x_min
    x_eval = np.linspace(x_min - 0.1 * x_range, x_max + 0.1 * x_range, 512)
    return x_eval, kde(x_eval)


def expand_bin_data(bin_positions, bin_values):
    counts = np.round(bin_values).astype(int)
    counts = np.maximum(counts, 0)
    return np.repeat(bin_positions, counts)


def plot_mRNA_metagene(ax, metagene, labels, adjust=0.5, colors=None):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    if colors is None:
        colors = plt.cm.tab10.colors
    mRNA_data = metagene.iloc[0:300].copy()
    bin_positions = mRNA_data.iloc[:, 1].values
    utr5_end = metagene.iloc[99, 1]
    cds_end = metagene.iloc[199, 1]
    utr3_end = metagene.iloc[299, 1]
    y_max_all = 0
    density_data = []

    for i, label in enumerate(labels):
        col_idx = 2 + i
        if col_idx >= metagene.shape[1]:
            print(f"Warning: Column {col_idx} not found for label '{label}'", file=sys.stderr)
            continue
        expanded = expand_bin_data(bin_positions, mRNA_data.iloc[:, col_idx].values)
        if len(expanded) < 2:
            print(f"Warning: Not enough data for label '{label}'", file=sys.stderr)
            continue
        x_density, y_density = compute_density(expanded, adjust)
        y_max_all = max(y_max_all, y_density.max())
        density_data.append((x_density, y_density, label, colors[i % len(colors)]))

    for x_density, y_density, label, color in density_data:
        ax.plot(x_density, y_density, color=color, linewidth=1.5, label=label)
        ax.fill_between(x_density, y_density, alpha=0.2, color=color)

    ax.axvline(x=utr5_end, color="gray", linestyle="dotted", linewidth=1)
    ax.axvline(x=cds_end, color="gray", linestyle="dotted", linewidth=1)
    rect_y_min = -y_max_all / 30 if y_max_all else -0.01
    rect_y_max = -y_max_all / 60 if y_max_all else -0.005
    rect_height_small = rect_y_max - rect_y_min
    big_rect_y_min = rect_y_min - rect_height_small / 4
    big_rect_y_max = rect_y_max + rect_height_small / 4
    ax.add_patch(Rectangle((0, rect_y_min), utr5_end, rect_height_small, facecolor="black", edgecolor="black", alpha=0.99))
    ax.add_patch(Rectangle((utr5_end, big_rect_y_min), cds_end - utr5_end, big_rect_y_max - big_rect_y_min, facecolor="gray", edgecolor="black", alpha=0.2))
    ax.add_patch(Rectangle((cds_end, rect_y_min), utr3_end - cds_end, rect_height_small, facecolor="black", edgecolor="black", alpha=0.99))
    text_y = big_rect_y_min + big_rect_y_max * 3
    ax.text(utr5_end / 2, text_y, "5' UTR", ha="center", va="bottom", fontsize=8)
    ax.text((utr5_end + cds_end) / 2, text_y, "CDS", ha="center", va="bottom", fontsize=8)
    ax.text((cds_end + utr3_end) / 2, text_y, "3' UTR", ha="center", va="bottom", fontsize=8)
    ax.set_title("mRNA Metagene Profile", fontsize=10)
    ax.set_ylabel("Density", fontsize=9)
    ax.set_xlim(0, utr3_end)
    ax.set_ylim(big_rect_y_min * 1.5, y_max_all * 1.1 if y_max_all else 1)
    ax.set_xticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.8)


def plot_ncRNA_metagene(ax, metagene, labels, adjust=0.5, colors=None):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    if colors is None:
        colors = plt.cm.tab10.colors
    ncRNA_data = metagene.iloc[300:600].copy()
    bin_positions = ncRNA_data.iloc[:, 1].values
    exon_first_end = metagene.iloc[399, 1]
    exon_internal_end = metagene.iloc[499, 1]
    exon_last_end = metagene.iloc[599, 1]
    y_max_all = 0
    density_data = []

    for i, label in enumerate(labels):
        col_idx = 2 + i
        if col_idx >= metagene.shape[1]:
            continue
        expanded = expand_bin_data(bin_positions, ncRNA_data.iloc[:, col_idx].values)
        if len(expanded) < 2:
            continue
        x_density, y_density = compute_density(expanded, adjust)
        y_max_all = max(y_max_all, y_density.max())
        density_data.append((x_density, y_density, label, colors[i % len(colors)]))

    for x_density, y_density, label, color in density_data:
        ax.plot(x_density, y_density, color=color, linewidth=1.5, label=label)
        ax.fill_between(x_density, y_density, alpha=0.2, color=color)

    ax.axvline(x=exon_first_end, color="gray", linestyle="dotted", linewidth=1)
    ax.axvline(x=exon_internal_end, color="gray", linestyle="dotted", linewidth=1)
    rect_y_min = -y_max_all / 30 if y_max_all else -0.01
    rect_y_max = -y_max_all / 60 if y_max_all else -0.005
    rect_height_small = rect_y_max - rect_y_min
    big_rect_y_min = rect_y_min - rect_height_small / 4
    big_rect_y_max = rect_y_max + rect_height_small / 4
    ax.add_patch(Rectangle((0, rect_y_min), exon_first_end, rect_height_small, facecolor="black", edgecolor="black", alpha=0.99))
    ax.add_patch(Rectangle((exon_first_end, big_rect_y_min), exon_internal_end - exon_first_end, big_rect_y_max - big_rect_y_min, facecolor="gray", edgecolor="black", alpha=0.2))
    ax.add_patch(Rectangle((exon_internal_end, rect_y_min), exon_last_end - exon_internal_end, rect_height_small, facecolor="black", edgecolor="black", alpha=0.99))
    text_y = big_rect_y_min + big_rect_y_max * 3
    ax.text(exon_first_end / 2, text_y, "First exon", ha="center", va="bottom", fontsize=8)
    ax.text((exon_first_end + exon_internal_end) / 2, text_y, "Internal exon", ha="center", va="bottom", fontsize=8)
    ax.text((exon_internal_end + exon_last_end) / 2, text_y, "Last exon", ha="center", va="bottom", fontsize=8)
    ax.set_title("ncRNA Metagene Profile", fontsize=10)
    ax.set_ylabel("Density", fontsize=9)
    ax.set_xlim(0, exon_last_end)
    ax.set_ylim(big_rect_y_min * 1.5, y_max_all * 1.1 if y_max_all else 1)
    ax.set_xticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.8)


def plot_metagene(bin_file, output_dir, labels, adjust=0.5):
    import matplotlib.pyplot as plt

    if isinstance(labels, str):
        labels = labels.split(",")
    metagene = load_metagene_data(bin_file)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams["font.size"] = 10
    plt.rcParams["axes.grid"] = False
    fig, axes = plt.subplots(2, 1, figsize=(4, 6))
    plot_mRNA_metagene(axes[0], metagene, labels, adjust=adjust)
    plot_ncRNA_metagene(axes[1], metagene, labels, adjust=adjust)
    plt.tight_layout()

    prefix = "metagene" if len(labels) > 1 else labels[0]
    outputs = []
    for fmt in ["pdf", "png", "svg"]:
        output_file = output_dir / f"{prefix}.{fmt}"
        plt.savefig(output_file, format=fmt, dpi=300, bbox_inches="tight")
        outputs.append(str(output_file))
    plt.close()
    return outputs


__all__ = [
    "REGION_TYPES",
    "RegionBinCalculator",
    "calculate_weights",
    "run_region2bin",
    "plot_metagene",
    "process_utr5_region",
    "process_cds_region",
    "process_utr3_region",
    "process_ncrna_region",
    "write_output",
]
