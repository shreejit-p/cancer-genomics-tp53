# %% [markdown]
# ## Notebook 02 — Data Cleaning & Quality Control
# We load raw sequences, filter out bad ones, and visualize basic stats.
# %%
# Setup path so python can find src/
import sys
sys.path.append("../src/")

from cleaning import load_and_filter
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt

# %%
# Run QC filtering on raw sequences
clean_records, qc_df = load_and_filter("../data/raw_sequences.fasta")

# preview the stats table
print(qc_df)

# %%
# Save the cleaned sequences and QC report
SeqIO.write(clean_records, "../data/clean_sequences.fasta", "fasta")
qc_df.to_csv("../results/qc_summary.csv", index = False)
print("Saved: clean_sequences.fasta")
print("Saved: qc_summary.csv")

# %%
# Plot GC content and sequence lenghts
fig, axes = plt.subplots(1,2,figsize = (15,5))

qc_df['gc_pct'].plot.hist(ax = axes[0],bins = 20,color='steelblue', edgecolor='white')
axes[0].set_title("GC content Distribution")
axes[0].set_xlabel("GC %")

qc_df['length'].plot.hist(ax= axes[1], bins = 20, color='seagreen', edgecolor='white')
axes[1].set_title("Sequence Length Distribution")
axes[1].set_xlabel("Length (bp)")

plt.tight_layout()
plt.savefig("../results/figures/qc_distributions.png", dpi = 150)
plt.show()
print("Saved: qc_distributions.png")

# %%
