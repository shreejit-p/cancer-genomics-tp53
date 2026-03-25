# %%
import pandas as pd

# Load all results
qc_df = pd.read_csv("/home/dragonfly/cancer-genomics-tp53/results/qc_summary.csv")
mut_df = pd.read_csv("/home/dragonfly/cancer-genomics-tp53/results/mutations.csv")
annotated_df = pd.read_csv("/home/dragonfly/cancer-genomics-tp53/results/mutations_annotated.csv")

# Summary
print("=== QC SUMMARY ===")
print(f"Total sequences: {len(qc_df)}")
print(f"Passed QC: {qc_df['passed_QC'].sum()}")
print(f"GC range: {qc_df['gc_pct'].min()} - {qc_df['gc_pct'].max()}")

print("\n=== MUTATION SUMMARY ===")
print(f"Total mutations: {len(mut_df)}")
print(mut_df['mutation_type'].value_counts().to_string())
print(mut_df['effect'].value_counts().to_string())

print("\n=== TOP MUTATED CODONS ===")
print(mut_df['codon_number'].value_counts().head(10).to_string())

print("\n=== CLINVAR SUMMARY ===")
print(f"Total annotated: {len(annotated_df)}")
print(f"Found in ClinVar: {annotated_df['in_clinvar'].sum()}")

print("\n=== MUTATIONS PER SAMPLE ===")
print(mut_df.groupby('sample_id').size().sort_values(ascending=False).to_string())

# %%

# %%
nonsense = mut_df[mut_df['effect'] == 'nonsense']
print(nonsense[['sample_id', 'position', 'codon_number', 'ref_base', 'alt_base']])

# %%