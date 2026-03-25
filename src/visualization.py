import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd
import logging

logger = logging.getLogger(__name__)


def lollipop_plot(mut_df, output_path):
    """
    Lollipop plot — standard cancer genomics visualization.
    X axis = codon position along TP53 protein (1-393)
    Y axis = number of mutations at that position
    Color  = mutation effect type
    """
    color_map = {
        'missense': '#E24B4A',
        'nonsense': '#185FA5',
        'silent':   '#888780',
        'deletion': '#BA7517',
        'insertion': '#639922',
        'unknown':  '#CCCCCC'
    }

    # Count mutations per codon per effect
    counts = mut_df.groupby(['codon_number', 'effect']).size().reset_index(name='count')

    fig, ax = plt.subplots(figsize=(16, 5))

    # TP53 protein domains background
    domains = [
        (1,   67,  '#F0F4FF', 'Transactivation'),
        (68,  98,  '#FFF8E1', 'Proline-rich'),
        (102, 292, '#F1FFF1', 'DNA-binding'),
        (293, 325, '#FFF0F0', 'Linker'),
        (356, 393, '#F5F0FF', 'Tetramerization'),
    ]

    for start, end, color, label in domains:
        ax.axvspan(start, end, alpha=0.3, color=color, label=label)

    # Draw lollipops
    for _, row in counts.iterrows():
        color = color_map.get(row['effect'], '#CCCCCC')
        ax.vlines(row['codon_number'], 0, row['count'],
                  color=color, linewidth=1.5, alpha=0.8)
        ax.scatter(row['codon_number'], row['count'],
                   color=color, s=50, zorder=5)

    # Known hotspots
    hotspots = {175: 'R175H', 248: 'R248W/Q', 273: 'R273H', 249: 'R249S'}
    for codon, label in hotspots.items():
        ax.axvline(x=codon, color='red', linestyle='--', alpha=0.4, linewidth=1)
        ax.text(codon, ax.get_ylim()[1] * 0.9, label,
                fontsize=8, color='red', ha='center')

    ax.set_xlabel("Codon position", fontsize=12)
    ax.set_ylabel("Mutation count", fontsize=12)
    ax.set_title("TP53 Mutation Landscape", fontsize=14, fontweight='bold')
    ax.set_xlim(0, 394)

    # Legend for mutation effects
    legend_effects = [mpatches.Patch(color=c, label=e)
                      for e, c in color_map.items()]
    ax.legend(handles=legend_effects, loc='upper right', fontsize=9)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    logger.info(f"Lollipop plot saved to {output_path}")


def mutation_heatmap(mut_df, output_path):
    """
    Heatmap showing which samples have mutations at which positions.
    Rows = samples, Columns = codon positions with mutations
    """
    # Pivot table — samples vs codon positions
    pivot = mut_df.pivot_table(
        index='sample_id',
        columns='codon_number',
        values='effect',
        aggfunc='first'
    )

    # Convert to binary (1 = has mutation, 0 = no mutation)
    binary = pivot.notna().astype(int)

    plt.figure(figsize=(16, 8))
    sns.heatmap(
        binary,
        cmap='Blues',
        cbar_kws={'label': 'Mutation present'},
        linewidths=0.1,
        linecolor='lightgray'
    )
    plt.title("TP53 Mutation Heatmap — Samples vs Codon Positions", fontsize=13)
    plt.xlabel("Codon position", fontsize=11)
    plt.ylabel("Transcript ID", fontsize=11)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_path}")


def gc_comparison_plot(qc_df, output_path):
    """
    Compare GC content across all sequences.
    """
    fig, ax = plt.subplots(figsize=(14, 5))

    colors = ['steelblue' if p else 'coral' for p in qc_df['passed_QC']]
    ax.bar(qc_df['id'], qc_df['gc_pct'], color=colors, edgecolor='white')
    ax.axhline(y=qc_df['gc_pct'].mean(), color='black',
               linestyle='--', label=f"Mean: {qc_df['gc_pct'].mean():.1f}%")

    ax.set_title("GC Content per TP53 Transcript", fontsize=13)
    ax.set_xlabel("Transcript ID", fontsize=11)
    ax.set_ylabel("GC %", fontsize=11)
    ax.tick_params(axis='x', rotation=90)
    ax.legend()

    passed = mpatches.Patch(color='steelblue', label='Passed QC')
    failed = mpatches.Patch(color='coral', label='Failed QC')
    ax.legend(handles=[passed, failed])

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    print(f"Saved: {output_path}")