from Bio import SeqIO
from BioSeqUtils import gc_fraction
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def load_and_filter (fasta_path, min_len = 1000, max_n_pct = 0.01):
    """
    Load FASTA file and filter low quality sequences.
    Parameters:
        fasta_path = path to raw .fasta file
        min_len = minimum sequence length to keep (default = 1000)
        max_n_pct = maximum percentage of Ns to keep (default = 0.01)
    Returns:
        clean_records = list of sequences that passed QC
        df = summary dataframe with stats of all sequences
    """
    records = list(SeqIO.parse(fasta_path,"fasta"))
    logger.info(f"Loaded {len(records)}")

    clean =[]
    stats = []

    for rec in records:
        seq = str(rec.seq).upper()
        n_pct = seq.count('N')/len(seq)
        gc = round(gc_fraction(seq)*100,2)

        passed = len(seq)>= min_len and n_pct<= max_n_pct

        stats.append({
            "id": rec.id,
            "length": len(seq),
            "gc_pct": gc,
            "passed_QC": passed
        })

        if passed:
            clean.append(rec)

    df = pd.DataFrame(stats)
    logger.info(f"QC passed: {len(clean)}/{len(records)} sequences")
    return clean, df