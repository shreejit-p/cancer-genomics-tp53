# TP53 Cancer Genomics Pipeline

A bioinformatics pipeline for analyzing TP53 tumor suppressor gene mutations
across all known human transcript variants. Built with Python, BioPython,
and MAFFT.

---

## Biological Background

TP53 (tumor protein p53) is the most frequently mutated gene in human cancers,
altered in approximately 50% of all cases. It acts as a tumor suppressor by
regulating cell cycle arrest, DNA repair, and apoptosis. Mutations in TP53,
particularly in its DNA binding domain, disrupt these functions and drive
cancer progression.

This pipeline analyzes all 25 validated TP53 transcript variants from NCBI
RefSeq to characterize sequence differences, detect mutations, and annotate
their clinical significance against the canonical reference NM_000546.6.

---

## Pipeline Overview
```
NCBI RefSeq (25 TP53 transcripts)
        ↓
Data Collection — BioPython Entrez
        ↓
Quality Control — GC%, length, N-base filtering
        ↓
Deduplication — remove identical sequences
        ↓
Multiple Sequence Alignment — MAFFT
        ↓
Mutation Detection & Classification
        ↓
Post-stop codon filtering
        ↓
ClinVar Annotation
        ↓
Publication Quality Visualizations
```

---

## Key Results

| Metric | Value |
|---|---|
| Transcripts fetched | 25 |
| Transcripts passing QC | 25 (100%) |
| GC content range | 50.94% — 53.59% |
| Alignment length | 2886 bp |
| Total variants detected | 5377 |
| SNPs | 138 |
| Deletions | 3428 |
| Insertions | 1811 |
| Missense mutations | 111 |
| Nonsense mutations | 2 |
| Silent mutations | 1836 |
| Variants found in ClinVar | 138/138 |

---

## Key Findings

**1. High sequence conservation across all transcripts**
All 25 TP53 transcripts passed QC with GC content ranging narrowly
from 50.94% to 53.59% — confirming these are high quality, biologically
consistent sequences from the same gene locus.

**2. Splicing differences dominate the variant landscape**
The majority of detected variants are deletions (3428) and insertions (1811),
reflecting alternative splicing of TP53 exons rather than point mutations.
This is expected as we are comparing transcript variants not tumor samples.

**3. Most divergent transcripts**
NM_001126116.2 and NM_001276698.3 showed the highest mutation counts
(642 each) against the canonical reference, indicating these represent
the most alternatively spliced TP53 variants.

**4. Missense mutations detected**
111 missense mutations were detected across all transcript comparisons,
representing amino acid changes relative to the canonical TP53 protein.

**5. Nonsense mutations at codon 32**
2 nonsense mutations were identified in transcripts NM_001407271.1 and
NM_001407270.1 at codon 32, causing premature protein truncation.
These likely reflect alternative transcription start sites rather than
pathogenic variants, as they occur in the divergent 5' regions of these
transcripts compared to NM_000546.6.

**6. ClinVar annotation**
All 138 detected SNPs returned ClinVar hits. This likely reflects the
broad nature of our gene-level ClinVar query rather than exact variant
matches — a known limitation of this annotation approach. Future work
should implement position-specific HGVS notation queries for precise
clinical significance matching.

---

## Visualizations

| Figure | Description |
|---|---|
| qc_distributions.png | GC% and length per transcript |
| gc_comparison.png | GC content with QC pass/fail coloring |
| mutation_distribution.png | Mutation types and effects breakdown |
| mutations_per_sample.png | Mutation count per transcript |
| lollipop_tp53.png | Mutation landscape across TP53 protein |
| mutation_heatmap.png | Heatmap of mutations across all transcripts |
| clinvar_coverage.png | ClinVar annotation coverage |

---

## Project Structure
```
cancer-genomics-tp53/
├── data/                      # raw data (gitignored)
│                              # regenerate by running notebook 01
├── notebooks/
│   ├── 01_data_collection.py  # fetch TP53 sequences from NCBI
│   ├── 02_cleaning.py         # QC filtering and deduplication
│   ├── 03_alignment.py        # MAFFT multiple sequence alignment
│   ├── 04_mutations.py        # mutation detection and classification
│   ├── 05_annotation.py       # ClinVar annotation and VCF export
│   └── 06_visualization.py    # publication quality figures
├── src/
│   ├── cleaning.py            # QC filtering and deduplication functions
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

## Limitations & Future Work

- ClinVar annotation uses gene-level queries — future work should
  implement exact HGVS notation for precise variant matching
- Analysis compares transcript variants not tumor vs normal samples —
  integrating TCGA controlled access data would enable true somatic
  mutation analysis
- Nonsense mutations at codon 32 require manual validation to confirm
  whether they represent true pathogenic variants or alternative
  transcription start sites
- Deduplication based on exact sequence identity — future work could
  use sequence similarity clustering (CD-HIT) for more robust deduplication

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

# Install MAFFT
sudo apt install mafft

# Run notebooks in order in VS Code
# Each notebook depends on output of previous one
01_data_collection.py  →  fetches data from NCBI
02_cleaning.py         →  QC filtering and deduplication
03_alignment.py        →  runs MAFFT alignment
04_mutations.py        →  detects and filters mutations
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

## Ethical Considerations

- Only publicly available de-identified data used
- No controlled access TCGA data used
- All sequences sourced from open NCBI RefSeq database
- ClinVar data accessed via public API within NCBI usage guidelines

---

## Author

Shreejit | [GitHub](https://github.com/shreejit-p)