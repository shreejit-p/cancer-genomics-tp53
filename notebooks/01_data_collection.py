# %% [markdown]
# ## Notebook 01 — Data Collection
# Fetching TP53 mRNA sequences from NCBI using BioPython Entrez

# %%
# Import required libraries
from Bio import Entrez
from Bio import SeqIO

# %%
# Set email — NCBI requires this to track API usage
Entrez.email = "shreejit.data@gmail.com"

# %%
# Search NCBI for human TP53 mRNA RefSeq sequences
search_handle = Entrez.esearch(
    db="nucleotide",
    term="TP53[Gene] AND Homo sapiens[Organism] AND mRNA[Filter] AND RefSeq[Filter]",
    retmax=10
)
search_results = Entrez.read(search_handle)
search_handle.close()

id_list = search_results["IdList"]
print("Retrieved IDs:")
print(id_list)

# %%
# Fetch sequences in FASTA format
fetch_handle = Entrez.efetch(
    db="nucleotide",
    id=id_list,
    rettype="fasta",
    retmode="text"
)
fasta_data = fetch_handle.read()
fetch_handle.close()

# Preview first 500 characters
print(fasta_data[:500])

# %%
# Save raw sequences to data folder
output_path = "/home/dragonfly/cancer-genomics-tp53/data/raw_sequences.fasta"
with open(output_path, "w") as file:
    file.write(fasta_data)
print("Data saved to:", output_path)

# %%
# Parse and inspect sequences
sequences = list(SeqIO.parse("/home/dragonfly/cancer-genomics-tp53/data/raw_sequences.fasta", "fasta"))
print("Total Sequences:", len(sequences))

for seq in sequences:
    print(f"ID: {seq.id}  |  Length: {len(seq)}")