# TP53 Cancer Genomics Pipeline

A bioinformatics pipeline for analyzing TP53 tumor suppressor gene mutations 
across all known human transcript variants. Built with Python, BioPython, 
and MAFFT.

## Project Overview

TP53 is mutated in approximately 50% of all human cancers. This pipeline:
- Fetches all validated TP53 transcripts from NCBI RefSeq
- Performs quality control filtering
- Aligns sequences using MAFFT (industry standard)
- Detects and classifies mutations against canonical reference NM_000546.6
- Annotates variants using ClinVar database
- Produces publication quality visualizations

## Results Summary

- 25 TP53 transcript variants analyzed
- Alignment length: 2886 bp
- Mutations detected and classified as silent, missense, nonsense
- Variants cross-referenced with ClinVar
- Key hotspots annotated: R175H, R248W/Q, R273H, R249S

## Project Structure
```
cancer-genomics-tp53/
├── data/              # raw data (gitignored - regenerate with notebook 01)
├── notebooks/         # step by step analysis
│   ├── 01_data_collection.py
│   ├── 02_cleaning.py
│   ├── 03_alignment.py
│   ├── 04_mutations.py
│   ├── 05_annotation.py
│   └── 06_visualization.py
├── src/               # reusable modules
│   ├── cleaning.py
│   ├── alignment.py
│   ├── mutations.py
│   ├── annotation.py
│   └── visualization.py
├── results/           # outputs
│   ├── figures/       # all plots
│   └── *.csv          # summary tables
└── requirements.txt
```

## Setup & Run
```bash
# Clone repository
git clone https://github.com/shreejit-p/cancer-genomics-tp53.git
cd cancer-genomics-tp53

# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Install MAFFT
sudo apt install mafft

# Run notebooks in order
# Open in VS Code and run cells top to bottom
```

## Tools & Technologies

| Tool | Version | Purpose |
|---|---|---|
| Python | 3.12 | Core language |
| BioPython | 1.83 | Sequence analysis |
| MAFFT | 7.x | Multiple sequence alignment |
| pandas | 2.x | Data manipulation |
| matplotlib | 3.x | Visualization |
| seaborn | 0.x | Statistical plots |

## Data Sources

- NCBI RefSeq — TP53 mRNA sequences
- ClinVar — Clinical variant significance
- Reference sequence: NM_000546.6 (canonical TP53)

## Ethical Considerations

- Only publicly available, de-identified data used
- No controlled access TCGA data used
- All data sourced from open NCBI databases

## Author

Shreejit | [GitHub](https://github.com/shreejit-p)
```

---

## `requirements.txt`
```
biopython==1.83
pandas==2.0.0
matplotlib==3.7.0
seaborn==0.12.0
numpy==1.24.0
requests==2.31.0