# TP53 Cancer Genomics Pipeline

A bioinformatics pipeline for analyzing TP53 tumor suppressor gene mutations
across all known human transcript variants. Built with Python, BioPython,
and MAFFT.

---

## Biological Background

TP53 (tumor protein p53) is the most frequently mutated gene in human cancers,
altered in approximately 50% of all cases. It acts as a tumor suppressor by
regulating cell cycle arrest, DNA repair, and apoptosis. Mutations in TP53
particularly in its DNA binding domain, disrupt these functions and drive
cancer progression.

This pipeline analyzes all 25 validated TP53 transcript variants from NCBI
RefSeq to characterize sequence differences, detect mutations, and annotate
their clinical significance.

---

## Pipeline Overview
```
NCBI RefSeq
    ↓
Data Collection (BioPython Entrez)
    ↓
Quality Control (GC%, length, N-base filtering)
    ↓
Multiple Sequence Alignment (MAFFT)
    ↓
Mutation Detection & Classification
    ↓
ClinVar Annotation
    ↓
Publication Quality Visualizations
```

---

## Key Results

| Metric | Value |
|---|---|
| Transcripts analyzed | 25 |
| Transcripts passing QC | 25 (100%) |
| GC content range | 50.94% — 53.59% |
| Alignment length | 2886 bp |
| Total variants detected | 5693 |
| SNPs | 172 |
| Deletions | 3444 |
| Insertions | 2077 |
| Missense mutations | 137 |
| Nonsense mutations | 4 |
| Silent mutations | 2108 |

### Key Findings

- All 25 TP53 transcript variants passed quality control with highly
  consistent GC content (50.94—53.59%), confirming data integrity
- The majority of variants between transcripts are deletions (3444) and
  insertions (2077), reflecting alternative splicing of TP53 exons rather
  than cancer-driving point mutations
- 172 SNPs were detected across all transcript comparisons against the
  canonical reference NM_000546.6
- 4 nonsense mutations were identified creating premature stop codons,
  representing the most clinically significant findings
- 137 missense mutations were detected that result in amino acid changes
  in the TP53 protein
- Variants were cross-referenced with ClinVar and exported in standard
  VCF 4.2 format for compatibility with downstream genomic tools

---

## Project Structure
```
cancer-genomics-tp53/
├── data/                      # raw data (gitignored)
│                              # regenerate by running notebook 01
├── notebooks/
│   ├── 01_data_collection.py  # fetch TP53 sequences from NCBI
│   ├── 02_cleaning.py         # QC filtering and stats
│   ├── 03_alignment.py        # MAFFT multiple sequence alignment
│   ├── 04_mutations.py        # mutation detection and classification
│   ├── 05_annotation.py       # ClinVar annotation and VCF export
│   └── 06_visualization.py    # publication quality figures
├── src/
│   ├── cleaning.py            # QC filtering functions
│   ├── alignment.py           # MAFFT wrapper functions
│   ├── mutations.py           # mutation detection and classification
│   ├── annotation.py          # ClinVar API and VCF writer
│   └── visualization.py       # lollipop, heatmap, GC plots
├── results/
│   ├── figures/
│   │   ├── qc_distributions.png
│   │   ├── gc_comparison.png
│   │   ├── mutation_distribution.png
│   │   ├── mutations_per_sample.png
│   │   ├── lollipop_tp53.png
│   │   ├── mutation_heatmap.png
│   │   └── clinvar_coverage.png
│   ├── qc_summary.csv
│   ├── alignment_summary.csv
│   ├── mutations.csv
│   ├── mutations_annotated.csv
│   └── mutations.vcf
├── README.md
├── requirements.txt
└── .gitignore
```

---

## Setup & Reproduction
```bash
# Clone repository
git clone https://github.com/shreejit-p/cancer-genomics-tp53.git
cd cancer-genomics-tp53

# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install Python dependencies
pip install -r requirements.txt

# Install MAFFT (required for alignment)
sudo apt install mafft        # Linux
brew install mafft            # Mac

# Run notebooks in order in VS Code
# Each notebook depends on output of previous one
01_data_collection.py  →  fetches data from NCBI
02_cleaning.py         →  filters sequences
03_alignment.py        →  runs MAFFT alignment
04_mutations.py        →  detects mutations
05_annotation.py       →  annotates with ClinVar
06_visualization.py    →  generates all figures
```

---

## Tools & Technologies

| Tool | Purpose |
|---|---|
| Python 3.12 | Core language |
| BioPython 1.83 | Sequence parsing, NCBI API, alignment IO |
| MAFFT 7.x | Multiple sequence alignment |
| pandas | Data manipulation and CSV output |
| matplotlib | Publication quality plotting |
| seaborn | Statistical visualizations |
| NCBI Entrez API | Data retrieval and ClinVar annotation |

---

## Data Sources

| Source | Data | Access |
|---|---|---|
| NCBI RefSeq | 25 TP53 mRNA transcripts | Public |
| ClinVar | Clinical variant significance | Public |
| Reference sequence | NM_000546.6 (canonical TP53) | Public |

---

## Output Files

| File | Description |
|---|---|
| qc_summary.csv | QC stats for all 25 sequences |
| alignment_summary.csv | Gap analysis per transcript |
| mutations.csv | All 5693 detected variants |
| mutations_annotated.csv | Variants with ClinVar annotation |
| mutations.vcf | Variants in standard VCF 4.2 format |

---

## Ethical Considerations

- Only publicly available, de-identified data used
- No controlled access TCGA data used
- All sequences sourced from open NCBI RefSeq database
- ClinVar data accessed via public API within NCBI usage guidelines

---

## Author

Shreejit | [GitHub](https://github.com/shreejit-p)