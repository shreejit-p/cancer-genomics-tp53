# %% [markdown]
# ## Notebook 02 — Data Cleaning & Quality Control
# We load raw sequences, filter out bad ones, and visualize basic stats.
# %%
# Setup path so python can find src/
import sys
sys.path.append("../src/")

from cleaning import load_and_filter
from cleaning import remove_duplicates
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt

# %%
# Run QC filtering on raw sequences
clean_records, qc_df = load_and_filter("../data/raw_sequences.fasta")

# remove duplicates
clean_records = remove_duplicates(clean_records)

# preview the stats table
print(qc_df)

# %%
# Save the cleaned sequences and QC report
SeqIO.write(clean_records, "../data/clean_sequences.fasta", "fasta")
qc_df.to_csv("../results/qc_summary.csv", index = False)
print(f"Saved {len(clean_records)} unique sequences")
print("Saved: qc_summary.csv")

# %%
# Bar charts — meaningful for small number of sequences
fig, axes = plt.subplots(2, 1, figsize=(14, 8))

# GC% per sequence
axes[0].bar(qc_df['id'], qc_df['gc_pct'], color='steelblue', edgecolor='white')
axes[0].set_title("GC Content per TP53 Transcript")
axes[0].set_ylabel("GC %")
axes[0].tick_params(axis='x', rotation=90)
axes[0].axhline(y=qc_df['gc_pct'].mean(), color='red', linestyle='--', label='mean')
axes[0].legend()

# Length per sequence
axes[1].bar(qc_df['id'], qc_df['length'], color='seagreen', edgecolor='white')
axes[1].set_title("Sequence Length per TP53 Transcript")
axes[1].set_ylabel("Length (bp)")
axes[1].tick_params(axis='x', rotation=90)
axes[1].axhline(y=qc_df['length'].mean(), color='red', linestyle='--', label='mean')
axes[1].legend()

plt.tight_layout()
plt.savefig("../results/figures/qc_distributions.png", dpi=150)
plt.show()
print("Saved: qc_distributions.png")
# %%

