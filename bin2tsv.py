import numpy as np
import pandas as pd
import os

def convert_bin_to_tsv(bin_filename, tsv_filename, num_rows, num_cols):
    """
    Reads a raw binary matrix (uint8) and converts it to a CSV.
    
    :param bin_filename: Path to your geno_pca.bin
    :param csv_filename: Path for the output CSV
    :param num_rows: Number of individuals
    :param num_cols: Number of SNPs
    """
    data = np.fromfile(bin_filename, dtype=np.uint8)
    
    if data.size != num_rows * num_cols:
        raise ValueError(f"Data size ({data.size}) does not match dimensions ({num_rows}x{num_cols})")

    matrix = data.reshape((num_rows, num_cols))

    geno = pd.DataFrame(matrix)

    geno.columns = [f'SNP_{i+1}' for i in range(num_cols)]

    geno.to_csv(tsv_filename, sep="\t", index=False)
    return geno

def merge(gen_party1, geno_party2):
    merged = pd.concat([gen_party1, geno_party2], axis=0)
    merged.reset_index(drop=True, inplace=True)

    merged_tsv = os.path.join(OUT_DIR, "geno_pca_merged.tsv")
    merged.to_csv(merged_tsv, sep="\t", index=False)


NUM_INDIVIDUALS = 1000 
NUM_SNPS = 10010

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(SCRIPT_DIR, "out")
os.makedirs(OUT_DIR, exist_ok=True)

bin_party1 = os.path.join(SCRIPT_DIR, "cache", "party1", "geno_pca1.bin")
bin_party2 = os.path.join(SCRIPT_DIR, "cache", "party2", "geno_pca2.bin")

geno_party1 = convert_bin_to_tsv(bin_party1, os.path.join(OUT_DIR, "geno_pca_party1.tsv"), NUM_INDIVIDUALS, NUM_SNPS)
geno_party2 = convert_bin_to_tsv(bin_party2, os.path.join(OUT_DIR, "geno_pca_party2.tsv"), NUM_INDIVIDUALS, NUM_SNPS)
merge(geno_party1, geno_party2)
