# %% [markdown]
# ## Notebook 06 — Visualization
# Publication quality figures for TP53 mutation analysis

# %%
import sys
sys.path.append("../src")

from visualization import lollipop_plot, mutation_heatmap, gc_comparison_plot
import pandas as pd

# %%
# Load results
mut_df = pd.read_csv("../results/mutations.csv")
qc_df = pd.read_csv("../results/qc_summary.csv")
print(f"Mutations loaded: {len(mut_df)}")
print(f"Sequences loaded: {len(qc_df)}")

# %%
# Lollipop plot
lollipop_plot(mut_df, "../results/figures/lollipop_tp53.png")

# %%
# Mutation heatmap
mutation_heatmap(mut_df, "../results/figures/mutation_heatmap.png")

# %%
# GC content comparison
gc_comparison_plot(qc_df, "../results/figures/gc_comparison.png")

print("All figures saved!")

# %%