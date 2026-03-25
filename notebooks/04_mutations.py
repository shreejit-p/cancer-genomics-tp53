# %% [markdown]
# ## Notebook 04 — Mutation Detection
# Compare all TP53 transcripts against reference NM_000546.6
# Detect and classify all mutations

# %%
import sys
sys.path.append("../src")

from mutations import detect_mutations
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# %%
# Detect all mutations
mut_df = detect_mutations("../data/aligned.fasta", reference_id="NM_000546.6")
print(f"Total mutations detected: {len(mut_df)}")
print(mut_df.head(10))

# %%
# Save mutations to CSV
mut_df.to_csv("../results/mutations.csv", index=False)
print("Saved: mutations.csv")

# %%
# Summary statistics
print("\nMutation types breakdown:")
print(mut_df['mutation_type'].value_counts())

print("\nMutation effects breakdown:")
print(mut_df['effect'].value_counts())

# %%
# Plot mutation type distribution
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

mut_df['mutation_type'].value_counts().plot.bar(
    ax=axes[0],
    color=['steelblue', 'seagreen', 'coral'],
    edgecolor='white'
)
axes[0].set_title("Mutation Types")
axes[0].set_xlabel("Type")
axes[0].set_ylabel("Count")
axes[0].tick_params(axis='x', rotation=0)

mut_df['effect'].value_counts().plot.bar(
    ax=axes[1],
    color=['coral', 'steelblue', 'seagreen', 'gray'],
    edgecolor='white'
)
axes[1].set_title("Mutation Effects")
axes[1].set_xlabel("Effect")
axes[1].set_ylabel("Count")
axes[1].tick_params(axis='x', rotation=0)

plt.tight_layout()
plt.savefig("../results/figures/mutation_distribution.png", dpi=150)
plt.show()
print("Saved: mutation_distribution.png")

# %%
# Mutations per sample
per_sample = mut_df.groupby('sample_id').size().reset_index(name='mutation_count')
per_sample = per_sample.sort_values('mutation_count', ascending=False)
print(per_sample)
per_sample.to_csv("../results/mutations_per_sample.csv", index=False)

# %%
# Plot mutations per sample
plt.figure(figsize=(14, 5))
plt.bar(per_sample['sample_id'], per_sample['mutation_count'],
        color='steelblue', edgecolor='white')
plt.title("Mutation Count per TP53 Transcript")
plt.xlabel("Transcript ID")
plt.ylabel("Number of Mutations")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig("../results/figures/mutations_per_sample.png", dpi=150)
plt.show()
print("Saved: mutations_per_sample.png")

# %%
