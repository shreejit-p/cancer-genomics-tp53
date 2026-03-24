from Bio import AlignIO
from Bio.Seq import Seq
import pandas as pd
import logging

logger = logging.getLogger(__name__)

def detect_mutations(alignment_path, reference_id = "NM_000546.6"):
    """
    Compare every sequence against reference and detect mutations.

    Parameters:
        alignment_path : path to aligned FASTA file
        reference_id   : ID of reference sequence (default NM_000546.6)

    Returns:
        DataFrame of all detected mutations with classification
    """
    alignment = AlignIO.read(alignment_path, "fasta")

    # Find reference sequence
    ref = None
    for rec in alignment:
        if reference_id in rec.id:
            ref = rec
            break

    if ref is None:
        raise ValueError(f"Reference sequence {reference_id} not found in alignment!")
    
    logger.info(f"Using reference: {ref.id}")

    ref_seq = str(ref.seq)
    mutations = []

    for record in alignment:
        if reference_id in record.id:
            continue

        seq = str(record.seq)
        for i, (ref_base, alt_base) in enumerate(zip(ref_seq, seq)):
            # Skip if same
            if ref_base == alt_base:
                continue

            # Skip if both are gaps
            if ref_base == "-" and alt_base == "-":
                continue

            mut_type = classify_mutation(ref_base, alt_base)
            effect = classify_effect(ref_seq, i, alt_base)


            mutations.append({
                "sample_id": record.id,
                "position": i + 1,
                "ref_base": ref_base,
                "alt_base": alt_base,
                "mutation_type": mut_type,
                "effect": effect,
                "codon_number": (i // 3) + 1
            })
    df = pd.DataFrame(mutations)
    logger.info(f"Total mutations detected: {len(df)}")
    return df

def classify_mutation(ref_base, alt_base):
    """Classify mutation as SNP, insertion or deletion."""
    if ref_base == '-':
        return "insertion"
    if alt_base == '-':
        return "deletion"
    return "SNP"

def classify_effect(ref_seq, pos, alt_base):
    """
    Translate codons and classify effect of mutation.
    Silent, Missense or Nonsense.
    """
    # Skip if alt_base is a gap
    if alt_base == '-':
        return "deletion"

    # Get codon start position
    codon_start = (pos // 3) * 3
    ref_codon = ref_seq[codon_start:codon_start + 3].replace('-', 'N')

    # Need complete codon
    if len(ref_codon) < 3:
        return "non-coding"

    # Build mutant codon
    mut_codon = list(ref_codon)
    mut_codon[pos % 3] = alt_base
    mut_codon = "".join(mut_codon)

    # Translate both codons
    try:
        aa_ref = str(Seq(ref_codon).translate())
        aa_mut = str(Seq(mut_codon).translate())
    except Exception:
        return "unknown"

    if aa_ref == aa_mut:
        return "silent"
    if aa_mut == '*':
        return "nonsense"
    return "missense"