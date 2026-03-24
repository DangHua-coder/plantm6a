#!/usr/bin/env python3
"""
保守m6A位点可视化模块（完整版）

绘制精美的Circos风格弦图展示物种间保守m6A分布

Features:
- 弧段严格按log(保守位点数)比例填满，无空白
- 弦颜色：两端物种颜色渐变
- 标准Circos弦图布局
"""
"""
plot_chord_v4.py
- 弧段严格按log(保守位点数)比例填满，无空白
- 弦颜色：两端物种颜色渐变
- 标准Circos弦图布局
"""

import argparse, os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from collections import defaultdict

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42
matplotlib.rcParams['font.family']  = 'Liberation Sans'
matplotlib.rcParams['font.sans-serif'] = ['Liberation Sans','DejaVu Sans']

# ── 配置 ──────────────────────────────────────────────────────────
SPECIES_ORDER = [
    "algae","marchantia",
    "arabidopsis","cucumber","tomato","potato","tobacco","soybean",
    "maize","rice",
    "fly","human",
]

LABELS = {
    "algae":"C. reinhardtii","arabidopsis":"A. thaliana",
    "cucumber":"C. sativus","fly":"D. melanogaster",
    "human":"H. sapiens","marchantia":"M. polymorpha",
    "maize":"Z. mays","potato":"S. tuberosum",
    "rice":"O. sativa","soybean":"G. max",
    "tobacco":"N. benthamiana","tomato":"S. lycopersicum",
}

SPECIES_COLORS = {
    "algae":       "#00897b",
    "marchantia":  "#43a047",
    "arabidopsis": "#e53935",
    "cucumber":    "#8e24aa",
    "tomato":      "#f06292",
    "potato":      "#ad1457",
    "tobacco":     "#6d4c41",
    "soybean":     "#fb8c00",
    "maize":       "#1e88e5",
    "rice":        "#039be5",
    "fly":         "#757575",
    "human":       "#3949ab",
}

R_OUTER  = 1.00
R_INNER  = 0.87
R_CHORD  = 0.84
R_LABEL  = 1.13
GAP_SP   = 0.025   # 物种间隔（弧度）
GAP_PAIR = 0.003   # 同物种内、不同弦之间的小间隔（弧度）


# ── 绘图工具 ──────────────────────────────────────────────────────
def arc_fill(ax, r_in, r_out, t1, t2, color, alpha=0.92):
    t = np.linspace(t1, t2, 120)
    xs = np.concatenate([r_out*np.cos(t), r_in*np.cos(t[::-1])])
    ys = np.concatenate([r_out*np.sin(t), r_in*np.sin(t[::-1])])
    ax.fill(xs, ys, color=color, alpha=alpha, lw=0, zorder=3)


def bezier_chord_gradient(ax, t1s, t1e, t2s, t2e, r, color_a, color_b, alpha):
    """分段渐变贝塞尔弦"""
    def bez(p0, p1, p2, p3, n=50):
        t = np.linspace(0, 1, n)[:, None]
        return (1-t)**3*p0 + 3*(1-t)**2*t*p1 + 3*(1-t)*t**2*p2 + t**3*p3

    ca = mcolors.to_rgba(color_a)
    cb = mcolors.to_rgba(color_b)
    n_seg = 8

    for i in range(n_seg):
        f0, f1 = i/n_seg, (i+1)/n_seg
        fmid = (f0+f1)/2
        ci = tuple(ca[c]*(1-fmid) + cb[c]*fmid for c in range(3))

        t1s_i = t1s + f0*(t1e-t1s)
        t1e_i = t1s + f1*(t1e-t1s)
        t2s_i = t2s + f0*(t2e-t2s)
        t2e_i = t2s + f1*(t2e-t2s)

        arc_a = np.column_stack([r*np.cos(np.linspace(t1s_i,t1e_i,8)),
                                  r*np.sin(np.linspace(t1s_i,t1e_i,8))])
        arc_b = np.column_stack([r*np.cos(np.linspace(t2s_i,t2e_i,8)),
                                  r*np.sin(np.linspace(t2s_i,t2e_i,8))])

        p0i = np.array([r*np.cos((t1s_i+t1e_i)/2), r*np.sin((t1s_i+t1e_i)/2)])
        p3i = np.array([r*np.cos((t2s_i+t2e_i)/2), r*np.sin((t2s_i+t2e_i)/2)])
        p1i, p2i = p0i*0.05, p3i*0.05

        fwd = bez(arc_a[-1], p1i, p2i, arc_b[0], n=25)
        bwd = bez(arc_b[-1], p2i, p1i, arc_a[0], n=25)

        xs = np.concatenate([arc_a[:,0], fwd[:,0], arc_b[:,0], bwd[:,0]])
        ys = np.concatenate([arc_a[:,1], fwd[:,1], arc_b[:,1], bwd[:,1]])
        ax.fill(xs, ys, color=ci, alpha=alpha, lw=0, zorder=2)


# ── 主函数 ────────────────────────────────────────────────────────
def plot_chord(df, output_prefix, min_sites=100, figsize=(9, 9)):
    species = SPECIES_ORDER

    # 1. 按物种对汇总总保守位点数
    pair_totals = defaultdict(int)
    for _, row in df.iterrows():
        a, b, n = row['spA'], row['spB'], int(row['count'])
        if a in species and b in species:
            pair_totals[(a, b)] += n

    pair_summary = {k: v for k, v in pair_totals.items() if v >= min_sites}
    if not pair_summary:
        print(f"无满足阈值 {min_sites} 的物种对")
        return

    # 2. 外圈弧长：按各物种参与的总保守位点数（线性）
    sp_total = defaultdict(int)
    for (a, b), n in pair_summary.items():
        sp_total[a] += n
        sp_total[b] += n

    grand = sum(sp_total[s] for s in species)
    total_angle = 2*np.pi - GAP_SP * len(species)

    arc = {}
    angle = np.pi/2
    for sp in species:
        span = (sp_total[sp] / grand) * total_angle
        arc[sp] = {
            'start': angle,
            'end':   angle + span,
            'mid':   angle + span/2,
            'span':  span,
        }
        angle += span + GAP_SP

    # 3. 为每个物种的弧段分配各弦的子段
    #    子段宽度按 log(n) 比例严格填满弧段（留 GAP_PAIR 小缝）
    #    排列顺序：按对方物种在 SPECIES_ORDER 中的顺序，保持对称美观
    chord_segs = {}   # (sp, partner) -> (seg_start, seg_end)

    for sp in species:
        # 找出该物种参与的所有对，按对方物种顺序排列
        partners = []
        for other in species:
            if other == sp:
                continue
            key = (sp, other) if (sp, other) in pair_summary else \
                  (other, sp) if (other, sp) in pair_summary else None
            if key:
                partners.append((other, pair_summary[key]))

        if not partners:
            continue

        # 计算 log 权重
        log_weights = {other: np.log1p(n) for other, n in partners}
        total_log   = sum(log_weights.values())

        # 可用弧宽（扣除间隔）
        n_gaps   = len(partners) - 1
        gap_total = GAP_PAIR * n_gaps
        usable   = arc[sp]['span'] - gap_total

        # 按权重分配，严格填满
        cursor = arc[sp]['start']
        for other, _ in partners:
            w    = log_weights[other]
            span = (w / total_log) * usable
            chord_segs[(sp, other)] = (cursor, cursor + span)
            cursor += span + GAP_PAIR

    # 4. 画布
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)

    # 5. 外圈弧
    for sp in species:
        a = arc[sp]
        arc_fill(ax, R_INNER, R_OUTER, a['start'], a['end'], SPECIES_COLORS[sp])
        for r_line in [R_INNER, R_OUTER]:
            t = np.linspace(a['start'], a['end'], 100)
            ax.plot(r_line*np.cos(t), r_line*np.sin(t),
                    color='white', lw=1.0, zorder=4)
        for t_end in [a['start'], a['end']]:
            ax.plot([R_INNER*np.cos(t_end), R_OUTER*np.cos(t_end)],
                    [R_INNER*np.sin(t_end), R_OUTER*np.sin(t_end)],
                    color='white', lw=1.0, zorder=4)

    # 6. 弦（按保守位点数从小到大画，大值在上层更清晰）
    sorted_pairs = sorted(pair_summary.items(), key=lambda x: x[1])

    for (a, b), n in sorted_pairs:
        # 获取两端弧段
        key_ab = (a, b) if (a, b) in chord_segs else None
        key_ba = (b, a) if (b, a) in chord_segs else None
        if key_ab is None or key_ba is None:
            continue

        a_s, a_e = chord_segs[(a, b)]
        b_s, b_e = chord_segs[(b, a)]

        ca = SPECIES_COLORS[a]
        cb = SPECIES_COLORS[b]

        # 透明度：log归一化，范围 0.35~0.80
        max_n = max(pair_summary.values())
        alpha = 0.35 + 0.45 * (np.log1p(n) / np.log1p(max_n))

        bezier_chord_gradient(ax, a_s, a_e, b_s, b_e,
                               R_CHORD, ca, cb, alpha)

    # 7. 标签
    for sp in species:
        mid = arc[sp]['mid']
        lx  = R_LABEL * np.cos(mid)
        ly  = R_LABEL * np.sin(mid)
        deg = np.degrees(mid)
        if -90 <= deg <= 90:
            ha, rot = 'left', deg
        else:
            ha, rot = 'right', deg + 180
        ax.text(lx, ly, LABELS[sp],
                ha=ha, va='center',
                rotation=rot, rotation_mode='anchor',
                fontsize=8.5, fontstyle='italic',
                color='#1a1a1a', zorder=5)

    # 8. 物种图例
    leg_patches = [
        mpatches.Patch(color=SPECIES_COLORS[sp], label=LABELS[sp])
        for sp in species
    ]
    leg = ax.legend(handles=leg_patches,
                    title='Species', title_fontsize=7.5,
                    fontsize=7,
                    loc='upper right',
                    bbox_to_anchor=(1.52, 1.02),
                    frameon=True, framealpha=0.95,
                    edgecolor='#cccccc',
                    handlelength=1.0, handleheight=0.9)
    leg.get_title().set_fontweight('bold')
    for t in leg.get_texts():
        t.set_fontstyle('italic')

    # 9. colorbar（弦透明度对应保守位点数）
    max_n = max(pair_summary.values())
    norm  = matplotlib.colors.LogNorm(vmin=max(min_sites, 1), vmax=max_n)
    sm    = plt.cm.ScalarMappable(cmap='Greys', norm=norm)
    sm.set_array([])
    cbar  = fig.colorbar(sm, ax=ax, shrink=0.22, aspect=14,
                         location='left', pad=0.12)
    cbar.set_label('Conserved m⁶A sites', fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    plt.tight_layout(pad=0.2)

    os.makedirs(os.path.dirname(output_prefix) or '.', exist_ok=True)
    for fmt in ['pdf', 'svg']:
        out = f"{output_prefix}.{fmt}"
        fig.savefig(out, format=fmt, dpi=300,
                    bbox_inches='tight', facecolor='white', edgecolor='none')
        print(f"✓ {out}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input',     default='summary/pair_region_breakdown.tsv')
    parser.add_argument('--output',    default='figures/conserved_m6A_chord')
    parser.add_argument('--min_sites', type=int, default=100)
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep='\t')
    print(f"读入 {len(df)} 行，共 {df['count'].sum()} 个位点")
    plot_chord(df, args.output, min_sites=args.min_sites)

# 导出的公共接口
__all__ = ['plot_chord', 'SPECIES_ORDER', 'LABELS', 'SPECIES_COLORS']

if __name__ == '__main__':
    main()
