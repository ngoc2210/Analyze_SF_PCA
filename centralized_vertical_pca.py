import os
import numpy as np
import pandas as pd
import scipy.sparse.linalg as lsa

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(SCRIPT_DIR, "out")
os.makedirs(OUT_DIR, exist_ok=True)

N_PCS        = 5       # number of principal components
N_SNPS       = 10010
N_SAMPLES    = 2000
ALPHA        = 10

# def _galinsky_scale(geno):
#     """
#     Scale a (n_samples x n_snps) genotype matrix using Galinsky et al.
#     Returns (n_samples x n_snps) scaled matrix.
#     """
#     geno = geno.astype(np.float64)
#     p = geno.mean(axis=0) / 2.0
#     denom = np.sqrt(2.0 * p * (1.0 - p))
#     monomorphic = denom <= 1e-10
#     denom[monomorphic] = 1.0
#     scaled = (geno - 2.0 * p) / denom
#     scaled[:, monomorphic] = 0.0
#     return scaled

def centralized_vertical_pca(data_file):
    if not os.path.isfile(data_file):
        print(f"{data_file} not found.")

    data_df = pd.read_csv(data_file, sep="\t")
    A = data_df.to_numpy(dtype=float)
    A = A.T
    print(f"Loaded data matrix A: {A.shape[0]} features x {A.shape[1]} samples")

    std = A.std(axis=0)
    std = np.where(std == 0, 1.0, std)
    A_std = (A - A.mean(axis=0)) / std
    k_eff = min(N_SAMPLES, N_SNPS, N_PCS+ALPHA)
    u, s, vt = lsa.svds(A_std, k=k_eff)
    idx = np.argsort(s)[::-1]
    s = s[idx]
    u = u[:, idx]
    vt = vt[idx]


    print("\n--- ground-truth ---")
    print("Top eigenvalues (explained variance):", s**2)
    print("feature PCs shape:", u.shape)
    print("sample PCs shape:", vt.shape)

    np.savetxt(
        os.path.join(OUT_DIR, "sample_eigenvectors.tsv"),
        vt,
        delimiter="\t"
    )

    np.savetxt(
        os.path.join(OUT_DIR, "SNP_eigenvectors.tsv"),
        u,
        delimiter="\t"
    )

if __name__ == "__main__":
  data_file = os.path.join(OUT_DIR, "geno_pca_merged.tsv")
  centralized_vertical_pca(data_file)
