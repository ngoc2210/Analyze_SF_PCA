"""
Run the cleartext Randomized PCA (RPCA) implementation (rpca.py) on the
pooled genomic dataset geno_pca_merged.tsv, and sanity-check the result
against scikit-learn's PCA -- the same kind of ground-truth comparison
SF-PCA itself uses in its accuracy evaluation (Section 7.5 of the paper).
"""

import numpy as np
import pandas as pd
import scipy.sparse.linalg as lsa
import os.path as path

from rpca import randomized_pca

# ---------------------------------------------------------------------
# Load data: 2000 samples (rows) x 10010 SNP features (columns), pooled
# from the 2 data providers (no further per-party split needed here).
# ---------------------------------------------------------------------
df = pd.read_csv("./out/geno_pca_merged.tsv", sep="\t")
A = df.to_numpy(dtype=float)
A = A.T
print(f"Loaded data matrix A: {A.shape[0]} features x {A.shape[1]} samples")

N_PCS = 5 
N_SNPS = 10010
N_SAMPLES = 2000
ALPHA = 10
P_ITERS = 20 

rpca_eigenvec, eigvals, A_proj = randomized_pca(A, N_PCS=N_PCS, alpha=ALPHA, p=P_ITERS, seed=0, standardize=True)

np.savetxt(
    path.join("out", "rpca_sample_eigenvectors.tsv"),
    rpca_eigenvec,
    delimiter="\t"
)

print("\n--- RPCA (this implementation, with standardization) ---")
print("Top eigenvalues (explained variance):", eigvals)
print("sample PCs (W) shape:", rpca_eigenvec.shape)
print("Projected data on sample eigenvectors shape:", A_proj.shape)

# std = A.std(axis=0)
# std = np.where(std == 0, 1.0, std)
# A_std = (A - A.mean(axis=0)) / std
# pca = PCA(n_components=N_PCS, svd_solver="full")
# true_A_proj = pca.fit_transform(A_std)
# true_eigenvec = pca.components_

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
print("sample PCs shape:", vt.shape)


print("\n--- Comparison with sklearn PCA (ground truth) ---")
for i in range(N_PCS):
    r = np.corrcoef(rpca_eigenvec[i, :], vt[i, :])[0, 1]
#     r = np.corrcoef(A_proj[:, i], true_A_proj[:, i])[0, 1]
    print(f"PC{i+1}: |Pearson r| = {abs(r):.4f}")
