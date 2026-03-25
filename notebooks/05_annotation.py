# %% [markdown]
# ## Notebook 05 — Variant Annotation
# Annotate detected mutations with ClinVar
# Export results in standard VCF format

# %%
import sys
sys.path.append("../src")

from annotation import annotate_with_clinvar, write_vcf
import pandas as pd
import matplotlib.pyplot as plt

# %%
# Load mutations from previous step
mut_df = pd.read_csv("../results/mutations.csv")
print(f"Loaded {len(mut_df)} mutations")
print(mut_df.head())

# %%
# Annotate with ClinVar
# Note: this takes a few minutes as we query NCBI API
print("Querying ClinVar... this may take a few minutes")
annotated_df = annotate_with_clinvar(mut_df)
print("Annotation complete!")
print(annotated_df.head(10))

# %%
# Save annotated results
annotated_df.to_csv("../results/mutations_annotated.csv", index=False)
print("Saved: mutations_annotated.csv")

# %%
# Summary
in_clinvar = annotated_df['in_clinvar'].sum()
total = len(annotated_df)
print(f"\nMutations found in ClinVar: {in_clinvar}/{total}")
print(f"Percentage: {round(in_clinvar/total*100, 2)}%")

# %%
# Export as VCF
write_vcf(mut_df, "../results/mutations.vcf")

# %%
# Plot ClinVar annotation results
labels = ['Found in ClinVar', 'Not in ClinVar']
values = [
    annotated_df['in_clinvar'].sum(),
    (~annotated_df['in_clinvar']).sum()
]

plt.figure(figsize=(6, 6))
plt.pie(values, labels=labels,
        colors=['steelblue', 'lightgray'],
        autopct='%1.1f%%',
        startangle=90)
plt.title("TP53 Mutations — ClinVar Coverage")
plt.tight_layout()
plt.savefig("../results/figures/clinvar_coverage.png", dpi=150)
plt.show()
print("Saved: clinvar_coverage.png")

# %%