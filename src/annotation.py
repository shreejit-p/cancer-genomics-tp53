# src/annotation.py

from Bio import Entrez
import pandas as pd
import logging
import time

logger = logging.getLogger(__name__)


def annotate_with_clinvar(mut_df, email="shreejit.data@gmail.com"):
    """
    Annotate mutations with ClinVar clinical significance.
    Uses NCBI Entrez API to look up each unique position.

    Parameters:
        mut_df : mutations DataFrame from detect_mutations()
        email  : your NCBI email

    Returns:
        annotated DataFrame with clinical significance column
    """
    Entrez.email = email

    # Only annotate SNPs — ClinVar mainly covers SNPs
    snp_df = mut_df[mut_df['mutation_type'] == 'SNP'].copy()
    logger.info(f"Annotating {len(snp_df)} SNPs with ClinVar")

    clinvar_results = []

    for _, row in snp_df.iterrows():
        # Build HGVS-like search query
        query = f"TP53[gene] AND {row['ref_base']}[ref] AND {row['alt_base']}[alt]"

        try:
            handle = Entrez.esearch(db="clinvar", term=query, retmax=1)
            result = Entrez.read(handle)
            handle.close()

            count = int(result['Count'])
            clinvar_results.append({
                "sample_id": row['sample_id'],
                "position": row['position'],
                "ref_base": row['ref_base'],
                "alt_base": row['alt_base'],
                "effect": row['effect'],
                "clinvar_hits": count,
                "in_clinvar": count > 0
            })

            # Be polite to NCBI — don't hammer their API
            time.sleep(0.34)

        except Exception as e:
            logger.warning(f"ClinVar lookup failed for position {row['position']}: {e}")
            clinvar_results.append({
                "sample_id": row['sample_id'],
                "position": row['position'],
                "ref_base": row['ref_base'],
                "alt_base": row['alt_base'],
                "effect": row['effect'],
                "clinvar_hits": 0,
                "in_clinvar": False
            })

    return pd.DataFrame(clinvar_results)


def write_vcf(mut_df, output_path):
    """
    Write mutations in standard VCF 4.2 format.
    VCF is the industry standard for variant data.

    Parameters:
        mut_df      : annotated mutations DataFrame
        output_path : where to save the VCF file
    """
    with open(output_path, 'w') as f:
        # VCF header
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=TP53_genomics_pipeline\n")
        f.write("##reference=NM_000546.6\n")
        f.write("##INFO=<ID=EFF,Number=1,Type=String,Description='Mutation effect'>\n")
        f.write("##INFO=<ID=TYPE,Number=1,Type=String,Description='Mutation type'>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for _, row in mut_df.iterrows():
            info = f"EFF={row['effect']};TYPE={row['mutation_type']}"
            f.write(f"17\t{row['position']}\t.\t{row['ref_base']}\t"
                   f"{row['alt_base']}\t.\tPASS\t{info}\n")

    logger.info(f"VCF saved to {output_path}")
    print(f"Saved: {output_path}")