import os
import numpy as np
import pandas as pd
import scipy.sparse.linalg as lsa

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
POOL_DIR = os.path.join(SCRIPT_DIR, "data", "pool")

N_PCS        = 20       # number of principal components
BATCH_SIZES  = list(range(70, 170, 10))   # 70, 80, …, 160
SNP_START    = 3        # column index where SNP columns begin in pool.csv
N_SNPS       = 100

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

    data_df = pd.read_csv(data_file, encoding="utf-8-sig")
    n_samples = len(data_df)

    iid_col  = data_df["IID"].values
    snp_cols = data_df.columns[SNP_START: SNP_START + N_SNPS].tolist()
    geno     = data_df[snp_cols].values.astype(np.int8)

    scaled = _galinsky_scale(geno)

    k_eff = min(N_PCS, n_samples, N_SNPS)
    u, s, vt = lsa.svds(scaled, k=k_eff)

    u  = np.flip(u,  axis=1)
    s  = np.flip(s)
    v  = np.flip(vt.T, axis=1)

    pc_labels = [f"PC{i+1}" for i in range(k_eff)]

    sample_eigenvec_df = pd.DataFrame(u, columns=pc_labels)
    sample_eigenvec_df.insert(0, "IID", iid_col)
    sample_eigenvec_df.to_csv(os.path.join(POOL_DIR, "G_sample_eigenvectors.csv"), index=False)

    feature_eigenvec_df = pd.DataFrame(v, index=snp_cols, columns=pc_labels)
    feature_eigenvec_df.index.name = "SNP"
    feature_eigenvec_df.to_csv(os.path.join(POOL_DIR, "H_SNP_eigenvectors.csv"), index=False)

    ev_df = pd.DataFrame({
        "PC":                    pc_labels,
        "singular_value":        s,
        "eigenvalue_approx":     s ** 2,
        "variance_explained_percentage": (s ** 2) / (s ** 2).sum() * 100,
    })
    ev_df.to_csv(os.path.join(POOL_DIR, "eigenvalues.csv"), index=False)

    print(f"{n_samples} individuals x {N_SNPS} SNPs  "
            f"| k={k_eff}  top-20 singular values: {np.round(s[:20], 4).tolist()}")

    print(f"Results written to: {POOL_DIR}")


if __name__ == "__main__":
  data_file = os.path.join(POOL_DIR, "pool.csv")
  centralized_vertical_pca(data_file)
