import os
import numpy as np
import pandas as pd
import scipy.sparse.linalg as lsa

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(SCRIPT_DIR, "out")
os.makedirs(OUT_DIR, exist_ok=True)

N_PCS        = 20       # number of principal components
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

    u  = np.flip(u,  axis=1)
    s  = np.flip(s)
    v  = np.flip(vt.T, axis=1)

    pc_labels = [f"PC{i+1}" for i in range(k_eff)]

    sample_eigenvec_df = pd.DataFrame(u, columns=pc_labels)
    sample_eigenvec_df.to_csv(os.path.join(OUT_DIR, "G_sample_eigenvectors.tsv"), sep="\t", index=False)

    feature_eigenvec_df = pd.DataFrame(v, index=snp_cols, columns=pc_labels)
    feature_eigenvec_df.index.name = "SNP"
    feature_eigenvec_df.to_csv(os.path.join(OUT_DIR, "H_SNP_eigenvectors.tsv"), sep="\t", index=False)

    ev_df = pd.DataFrame({
        "PC":                    pc_labels,
        "singular_value":        s,
        "eigenvalue_approx":     s ** 2,
        "variance_explained_percentage": (s ** 2) / (s ** 2).sum() * 100,
    })
    ev_df.to_csv(os.path.join(OUT_DIR, "eigenvalues.tsv"), sep="\t", index=False)

    print(f"{n_samples} individuals x {N_SNPS} SNPs  "
            f"| k={k_eff}  top-20 singular values: {np.round(s[:20], 4).tolist()}")

if __name__ == "__main__":
  data_file = os.path.join(OUT_DIR, "geno_pca_merged.tsv")
  centralized_vertical_pca(data_file)
