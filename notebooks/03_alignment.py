# %% [markdown]
# ## Notebook 03 — Multiple Sequence Alignment
# We align all 25 clean TP53 sequences using MAFFT
# This allows us to compare sequences position by position

# %%
#Setup
import sys
sys.path.append("../src")
from Bio import AlignIO
from alignment import run_mafft, alignment_summary
import pandas as pd

# %%
#Define file paths
input_fasta = "../data/clean_sequences.fasta"
output_fasta = "../data/aligned.fasta"

# %%
# Run MAFFT alignment
print("Running MAFFT alignment...")
alignment = run_mafft(input_fasta, output_fasta)
print("Alignment complete!")

# %%
# Print alignment summary
alignment_summary(alignment)

# %%
# Save alignment stats to CSV
stats = []
for rec in alignment:
    stats.append({
        "id": rec.id,
        "aligned_length": len(rec.seq),
        "gap_count": str(rec.seq).count('-'),
        "gap_pct": round(str(rec.seq).count('-') / len(rec.seq) * 100, 2)
    })

df = pd.DataFrame(stats)
print(df)
df.to_csv("../results/alignment_summary.csv", index=False)
print("Saved: alignment_summary.csv")

# %%


