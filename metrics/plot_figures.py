#!/usr/bin/env python3
import pandas as pd
import os
import sys

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except Exception as e:
    print("[error] matplotlib not available: ", e)
    print("Please install with: pip install matplotlib")
    sys.exit(1)

# Resolve paths robustly for both monorepo and standalone repo usage
AGG_PATH = sys.argv[1] if len(sys.argv) > 1 else 'metrics/results/aggregate_metrics.tsv'
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else 'figures'
if not os.path.exists(AGG_PATH):
    alt = 'AdaptiveInterpreter/metrics/results/aggregate_metrics.tsv'
    if os.path.exists(alt):
        AGG_PATH = alt

os.makedirs(OUT_DIR, exist_ok=True)

# Read aggregate metrics (expects a row with gene==ALL)
df = pd.read_csv(AGG_PATH, sep='\t')
if 'gene' not in df.columns:
    print('[error] aggregate file missing gene column:', df.columns.tolist())
    sys.exit(2)
row = df.loc[df['gene'] == 'ALL']
if row.empty:
    row = df.iloc[[0]]

row = row.iloc[0].to_dict()
N = float(row.get('N', 0)) or 0.0
agree = float(row.get('AGREE', 0))
better = float(row.get('BETTER_DATA', 0))
disagree = float(row.get('DISAGREE', 0))

if N <= 0:
    print('[error] N is zero in aggregate metrics; cannot plot')
    sys.exit(3)

# Compute percentages
p_ag = agree / N * 100.0
p_bd = better / N * 100.0
p_dg = disagree / N * 100.0

# Figure 2: Directional Agreement (DAL) stacked bar
fig, ax = plt.subplots(figsize=(7, 1.2))
left = 0.0
colors = {'AGREE': '#2ca02c', 'BETTER_DATA': '#1f77b4', 'DISAGREE': '#d62728'}

# Draw stacked horizontal bar
ax.barh([0], [p_ag], left=left, color=colors['AGREE'], label=f'AGREE {p_ag:.2f}% ({int(agree):,})')
left += p_ag
ax.barh([0], [p_bd], left=left, color=colors['BETTER_DATA'], label=f'BETTER_DATA {p_bd:.2f}% ({int(better):,})')
left += p_bd
ax.barh([0], [p_dg], left=left, color=colors['DISAGREE'], label=f'DISAGREE {p_dg:.2f}% ({int(disagree):,})')

ax.set_xlim(0, 100)
ax.set_yticks([])
ax.set_xlabel('Percentage of variants (N = {:,})'.format(int(N)))
ax.set_title('Figure 2. Directional Agreement (DAL): AGREE / BETTER_DATA / DISAGREE')

# Legend below
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.55), ncol=3, frameon=False)
plt.tight_layout()

svg_path = os.path.join(OUT_DIR, 'fig2_dal_stacked.svg')
png_path = os.path.join(OUT_DIR, 'fig2_dal_stacked.png')
plt.savefig(svg_path, bbox_inches='tight')
plt.savefig(png_path, dpi=300, bbox_inches='tight')
print('[write]', svg_path)
print('[write]', png_path)


# Figure 1: Confusion matrices (lenient and strict)
try:
    tp_l = int(row.get('lenient_TP', 0)); fp_l = int(row.get('lenient_FP', 0))
    tn_l = int(row.get('lenient_TN', 0)); fn_l = int(row.get('lenient_FN', 0))
    tp_s = int(row.get('strict_TP', 0)); fp_s = int(row.get('strict_FP', 0))
    tn_s = int(row.get('strict_TN', 0)); fn_s = int(row.get('strict_FN', 0))

    fig, axes = plt.subplots(1, 2, figsize=(8, 3))
    for ax, (tp, fp, fn, tn, name) in zip(
        axes,
        [
            (tp_l, fp_l, fn_l, tn_l, 'Lenient'),
            (tp_s, fp_s, fn_s, tn_s, 'Strict'),
        ],
    ):
        mat = [[tp, fp],[fn, tn]]
        im = ax.imshow(mat, cmap='Blues')
        for i in range(2):
            for j in range(2):
                ax.text(j, i, f"{mat[i][j]:,}", va='center', ha='center', color='black')
        ax.set_xticks([0,1]); ax.set_yticks([0,1])
        ax.set_xticklabels(['Pred +','Pred -'])
        ax.set_yticklabels(['True +','True -'])
        ax.set_title(f'Figure 1. {name} confusion matrix')
    plt.tight_layout()
    svg1 = os.path.join(OUT_DIR, 'fig1_confusion_matrices.svg')
    png1 = os.path.join(OUT_DIR, 'fig1_confusion_matrices.png')
    plt.savefig(svg1, bbox_inches='tight')
    plt.savefig(png1, dpi=300, bbox_inches='tight')
    print('[write]', svg1)
    print('[write]', png1)
except Exception as e:
    print('[warn] skipping confusion matrices:', e)

# Figure 3: Per-gene top-15 DISAGREE rate
PG_PATH = 'metrics/results/per_gene_metrics.tsv'
if not os.path.exists(PG_PATH):
    alt_pg = 'AdaptiveInterpreter/metrics/results/per_gene_metrics.tsv'
    if os.path.exists(alt_pg):
        PG_PATH = alt_pg
try:
    pg = pd.read_csv(PG_PATH, sep='\t')
    pg = pg[pg['gene'].astype(str) != 'ALL']
    pg = pg[pg['N'] > 0]
    # ensure numeric
    pg['DISAGREE'] = pd.to_numeric(pg['DISAGREE'], errors='coerce').fillna(0)
    pg['N'] = pd.to_numeric(pg['N'], errors='coerce').fillna(0)
    pg['dg_rate'] = pg['DISAGREE'] / pg['N'] * 100.0
    top = pg.sort_values('dg_rate', ascending=False).head(15)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.barh(top['gene'], top['dg_rate'], color='#d62728')
    ax.invert_yaxis()
    for i, (g, r, cnt) in enumerate(zip(top['gene'], top['dg_rate'], top['DISAGREE'])):
        ax.text(r + 0.5, i, f"{r:.2f}% ({int(cnt)})", va='center', fontsize=8)
    ax.set_xlabel('DISAGREE rate (%)')
    ax.set_title('Figure 3. Top-15 genes by DISAGREE rate (DAL)')
    plt.tight_layout()
    svg3 = os.path.join(OUT_DIR, 'fig3_per_gene_disagree_top15.svg')
    png3 = os.path.join(OUT_DIR, 'fig3_per_gene_disagree_top15.png')
    plt.savefig(svg3, bbox_inches='tight')
    plt.savefig(png3, dpi=300, bbox_inches='tight')
    print('[write]', svg3)
    print('[write]', png3)
except Exception as e:
    print('[warn] skipping per-gene DISAGREE plot:', e)

# Figure 4: DISAGREE ClinVar review quality distribution
try:
    pdg_path = 'metrics/results/per_gene_disagrees.tsv'
    if not os.path.exists(pdg_path):
        alt_pd = 'AdaptiveInterpreter/metrics/results/per_gene_disagrees.tsv'
        if os.path.exists(alt_pd):
            pdg_path = alt_pd
    if os.path.exists(pdg_path):
        dd = pd.read_csv(pdg_path, sep='\t')
        qcol = 'clinvar_quality_score' if 'clinvar_quality_score' in dd.columns else None
        if qcol is None:
            # try extract from review_flags if present
            if 'review_flags' in dd.columns:
                dd['__q__'] = dd['review_flags'].astype(str).str.extract(r'(\d+)')
                qcol = '__q__'
        if qcol:
            vals = pd.to_numeric(dd[qcol], errors='coerce')
            counts = vals.value_counts(dropna=False).sort_index()
            labels = [str(int(k)) if pd.notna(k) else 'NA' for k in counts.index.tolist()]
            fig, ax = plt.subplots(figsize=(4.5, 3))
            ax.bar(labels, counts.values, color='#9467bd')
            for i, v in enumerate(counts.values):
                ax.text(i, v + max(counts.values)*0.01 + 0.5, f"{int(v):,}", ha='center', fontsize=8)
            ax.set_xlabel('ClinVar review stars')
            ax.set_ylabel('Count of DISAGREE entries')
            ax.set_title('Figure 4. DISAGREE entries by ClinVar review status')
            plt.tight_layout()
            svg4 = os.path.join(OUT_DIR, 'fig4_disagree_review_quality.svg')
            png4 = os.path.join(OUT_DIR, 'fig4_disagree_review_quality.png')
            plt.savefig(svg4, bbox_inches='tight')
            plt.savefig(png4, dpi=300, bbox_inches='tight')
            print('[write]', svg4)
            print('[write]', png4)
        else:
            print('[warn] no quality column found in per_gene_disagrees.tsv')
    else:
        print('[warn] per_gene_disagrees.tsv not found; skipping Figure 4')
except Exception as e:
    print('[warn] skipping DISAGREE quality figure:', e)


# Figure 5: Top-15 genes by (AGREE + BETTER_DATA) rate (stacked)
try:
    pg = pd.read_csv(PG_PATH, sep='\t') if 'pg' not in locals() else pg
    pg = pg[pg['gene'].astype(str) != 'ALL']
    pg['AGREE'] = pd.to_numeric(pg['AGREE'], errors='coerce').fillna(0)
    pg['BETTER_DATA'] = pd.to_numeric(pg['BETTER_DATA'], errors='coerce').fillna(0)
    pg['N'] = pd.to_numeric(pg['N'], errors='coerce').fillna(0)
    pg = pg[pg['N'] > 0]
    pg['agree_rate'] = pg['AGREE'] / pg['N'] * 100.0
    pg['better_rate'] = pg['BETTER_DATA'] / pg['N'] * 100.0
    pg['combo_rate'] = pg['agree_rate'] + pg['better_rate']
    top = pg.sort_values('combo_rate', ascending=False).head(15)

    fig, ax = plt.subplots(figsize=(7, 5.0))
    y = range(len(top))
    ax.barh(top['gene'], top['agree_rate'], color='#2ca02c', label='AGREE %')
    ax.barh(top['gene'], top['better_rate'], left=top['agree_rate'], color='#1f77b4', label='BETTER_DATA %')
    ax.invert_yaxis()
    for i, (ga, bd) in enumerate(zip(top['agree_rate'], top['better_rate'])):
        total = ga + bd
        ax.text(total + 0.6, i, f"{total:.2f}%", va='center', fontsize=8)
    ax.set_xlabel('Rate (%)')
    ax.set_title('Figure 5. Top-15 genes by (AGREE + BETTER_DATA) rate (DAL)')
    ax.legend(loc='lower right', frameon=False)
    plt.tight_layout()
    svg5 = os.path.join(OUT_DIR, 'fig5_top15_agree_plus_better_dal.svg')
    png5 = os.path.join(OUT_DIR, 'fig5_top15_agree_plus_better_dal.png')
    plt.savefig(svg5, bbox_inches='tight')
    plt.savefig(png5, dpi=300, bbox_inches='tight')
    print('[write]', svg5)
    print('[write]', png5)
except Exception as e:
    print('[warn] skipping top-15 AGREE+BETTER figure:', e)
