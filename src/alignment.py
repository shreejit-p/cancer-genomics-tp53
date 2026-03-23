import subprocess
from Bio import AlignIO
import logging

logger = logging.getLogger(__name__)

def run_mafft(input_fasta, output_fasta):
    """
    Runs MAFFT alignment on input FASTA file.

    Parameters:
        input_fasta  : path to clean sequences FASTA
        output_fasta : path to save aligned sequences
    
    Returns:
        alignment object
    """
    cmd = ["mafft", "--auto", "--thread", "-1", input_fasta]
    logger.info(f"Running MAFFT on {input_fasta}")

    with open(output_fasta, "w") as out:
        subprocess.run(cmd, stdout = out, stderr=subprocess.DEVNULL, check= True)
    
    alignment = AlignIO.read(output_fasta, "fasta")
    logger.info(f"Alignment complete: {len(alignment)} sequences, {alignment.get_alignment_length()} columns")

    return alignment

def alignment_summary(alignment):
    """
    Print basic stats about the alignment.
    """
    print(f"Number of sequences  : {len(alignment)}")
    print(f"Alignment length     : {alignment.get_alignment_length()}")
    print(f"Sequence IDs         : ")
    for rec in alignment:
        print(f" {rec.id}")