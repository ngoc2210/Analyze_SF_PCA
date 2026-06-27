import os
import numpy as np
import pandas as pd
import scipy.sparse.linalg as lsa

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(SCRIPT_DIR, "out")
os.makedirs(OUT_DIR, exist_ok=True)

N_PCS        = 5       # number of principal components
N_SNPS       = 10010

def _galinsky_scale(geno):
    """
    Scale a (n_samples x n_snps) genotype matrix using Galinsky et al.
    Returns (n_samples x n_snps) scaled matrix.
    """
    geno = geno.astype(np.float64)
    p = geno.mean(axis=0) / 2.0
    denom = np.sqrt(2.0 * p * (1.0 - p))
    monomorphic = denom <= 1e-10
    denom[monomorphic] = 1.0
    scaled = (geno - 2.0 * p) / denom
    scaled[:, monomorphic] = 0.0
    return scaled

def centralized_vertical_pca(data_file):
    if not os.path.isfile(data_file):
        print(f"{data_file} not found.")

    data_df = pd.read_csv(data_file, sep="\t", encoding="utf-8-sig")
    n_samples = len(data_df)

    snp_cols = data_df.columns[0: N_SNPS].tolist()
    geno     = data_df[snp_cols].values.astype(np.int8)

    scaled = _galinsky_scale(geno)

    k_eff = min(N_PCS, n_samples, N_SNPS)
    u, s, vt = lsa.svds(scaled, k=k_eff)

    idx = np.argsort(s)[::-1]
    s = s[idx]
    u = u[:, idx]
    vt = vt[idx]

    sample_eigenvec = u.T
    feature_eigenvec = vt

    np.savetxt(
        os.path.join(OUT_DIR, "sample_eigenvectors.tsv"),
        sample_eigenvec,
        delimiter="\t"
    )

    np.savetxt(
        os.path.join(OUT_DIR, "SNP_eigenvectors.tsv"),
        feature_eigenvec,
        delimiter="\t"
    )

    np.savetxt(
        os.path.join(OUT_DIR, "eigenvalues.tsv"),
        np.column_stack([
            s,
            s**2,
            (s**2) / np.sum(s**2) * 100
        ]),
        delimiter="\t"
    )

    print(f"{n_samples} individuals x {N_SNPS} SNPs  "
            f"| k={k_eff}  top-20 singular values: {np.round(s[:20], 4).tolist()}")

if __name__ == "__main__":
  data_file = os.path.join(OUT_DIR, "geno_pca_merged.tsv")
  centralized_vertical_pca(data_file)
