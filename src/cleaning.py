from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
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

def remove_duplicates(records, keep_id="NM_000546.6"):
    """
    Remove duplicate sequences based on sequence content.
    Always keeps the reference sequence regardless of duplicates.
    
    Parameters:
        records : list of SeqRecord objects
        keep_id : ID of sequence to always keep (reference)
    """
    seen = set()
    unique = []
    duplicates = []

    # First pass — always add reference sequence
    for rec in records:
        if keep_id in rec.id:
            seen.add(str(rec.seq).upper())
            unique.append(rec)
            print(f"Reference protected: {rec.id}")
            break

    # Second pass — add remaining unique sequences
    for rec in records:
        if keep_id in rec.id:
            continue
        seq = str(rec.seq).upper()
        if seq not in seen:
            seen.add(seq)
            unique.append(rec)
        else:
            duplicates.append(rec.id)

    print(f"Removed {len(duplicates)} duplicate sequences")
    print(f"Duplicate IDs: {duplicates}")
    print(f"Remaining: {len(unique)} unique sequences")
    return unique