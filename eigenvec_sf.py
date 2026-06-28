import math
import numpy as np
import pandas as pd

# --- Configuration (from configGlobal.toml) ---
NPC = 5
NUM_OVERSAMPLE = 10
KP = NPC + NUM_OVERSAMPLE   # 15
N_POWER_ITER = 20

# --- Load input data ---
X1 = pd.read_csv("./out/geno_pca_party1.tsv", sep="\t").to_numpy(dtype=float)  # n1 x nsnp
X2 = pd.read_csv("./out/geno_pca_party2.tsv", sep="\t").to_numpy(dtype=float)  # n2 x nsnp
n1, nsnp = X1.shape
n2 = X2.shape[0]
n = n1 + n2
print(f"X1 shape {X1.shape}")
print(f"X2 shape {X2.shape}")

# --- Compute mean and stdinv from combined data ---
pool_data = pd.read_csv("./out/geno_pca_merged.tsv", sep="\t").to_numpy(dtype=float)
mean = pool_data.mean(axis=0)                          # (nsnp,)
x2mean = (pool_data ** 2).mean(axis=0)
var = x2mean - mean ** 2
var = np.where(var < 1e-8, 1.0, var)                  # replace near-zero variance (mirrors SF-GWAS)
stdinv = 1.0 / np.sqrt(var)                           # 1/std as in SF-GWAS XStdInv

X = pool_data                                          # combined n x nsnp

# Scaling factors used in SF-GWAS power iteration
num_snp_sqrt_inv = 1.0 / math.sqrt(nsnp)
num_tot_ind_sqrt_inv = 1.0 / math.sqrt(n)

# Save mean and stdinv for comparison
np.savetxt("./out/Xmean.tsv", mean, delimiter="\t")
np.savetxt("./out/XStdInv.tsv", stdinv, delimiter="\t")


def matmult_lazynorm(Q, X, mean, stdinv, is_col):
    """
    Multiply Q by normalized X without materializing X_norm explicitly.
    Mirrors QXLazyNormStream (is_col=False) and QXtLazyNormStream (is_col=True) in SF-GWAS.

    is_col=True:  Q (kp x nind), X (nind x nsnp) -> (kp x nsnp)  [Q @ X_norm_cols]
    is_col=False: Q (kp x nsnp), X (nsnp x nind) -> (kp x nind)  [Q_scaled @ X_norm_rows]
    """
    if is_col:
        D = Q @ X
        correction = Q.sum(axis=1, keepdims=True) @ mean.reshape(1, -1)
        E = D - correction
        R = E * stdinv.reshape(1, -1)
        return R
    else:
        D = Q * stdinv.reshape(1, -1)
        E = D @ X
        correction = (D @ mean.reshape(-1, 1)) @ np.ones((1, X.shape[1]))
        R = E - correction
        return R


# --- Step 1: Load SF-GWAS sketch (aggregated, pre-normalization) from cache ---
Q_sketch = np.loadtxt("./cache/party1/Sketch.txt", delimiter=",")  # kp x nsnp

buckets1 = np.loadtxt("./cache/party1/SketchBucketId.txt").astype(int)
sgn1     = np.loadtxt("./cache/party1/SketchSign.txt")
buckets2 = np.loadtxt("./cache/party2/SketchBucketId.txt").astype(int)
sgn2     = np.loadtxt("./cache/party2/SketchSign.txt")
all_buckets = np.concatenate([buckets1, buckets2])
all_sgn     = np.concatenate([sgn1, sgn2])

bucket_count = np.bincount(all_buckets, minlength=KP)
pos_count    = np.bincount(all_buckets[all_sgn > 0].astype(int), minlength=KP)

# Normalize: (Q - meanWeight * mean) / sqrt(bucketCount) * stdinv  (SF-GWAS order)
for b in range(KP):
    mean_weight = float(2 * pos_count[b] - bucket_count[b])
    Q_sketch[b] = (Q_sketch[b] - mean_weight * mean) / math.sqrt(bucket_count[b]) * stdinv


# --- Step 2: Initial QXLazyNorm -> Q (kp x n) via QR ---
Q_init_x = matmult_lazynorm(Q_sketch, X.T, mean, stdinv, is_col=False) * num_snp_sqrt_inv
Q, _ = np.linalg.qr(Q_init_x.T)  # n x kp
Q = Q.T                            # kp x n


# --- Step 3: Power iterations ---
for it in range(N_POWER_ITER):
    # QXtLazyNorm: Q (kp x n) @ X_norm (n x nsnp) -> kp x nsnp, scale by 1/sqrt(n)
    Q_tmp = matmult_lazynorm(Q, X, mean, stdinv, is_col=True) * num_tot_ind_sqrt_inv

    # QXLazyNorm: Q_tmp (kp x nsnp) @ X_norm.T (nsnp x n) -> kp x n, scale by 1/sqrt(nsnp)
    Q_tmp = matmult_lazynorm(Q_tmp, X.T, mean, stdinv, is_col=False) * num_snp_sqrt_inv

    if it < N_POWER_ITER - 1:
        Q, _ = np.linalg.qr(Q_tmp.T)
        Q = Q.T
    else:
        Q = Q_tmp   # last iteration: skip QR (mirrors SF-GWAS)

# Q: kp x n


# --- Step 4: Gram matrix Z = Q @ Q' / n (kp x kp), aggregate across parties ---
Z = (Q @ Q.T) / n


# --- Step 5: Eigen decomposition, sort descending, take top NPC rows ---
eigenvalues, eigenvectors = np.linalg.eigh(Z)
idx = np.argsort(eigenvalues)[::-1]
V = eigenvectors[:, idx[:NPC]].T   # npc x kp


# --- Step 6: Per-party PC projections (npc x nind_local) ---
Qpc1 = V @ Q[:, :n1]   # npc x n1
Qpc2 = V @ Q[:, n1:]   # npc x n2

print(f"Qpc1 shape: {Qpc1.shape}")
print(f"Qpc2 shape: {Qpc2.shape}")

np.savetxt("./out/Qpc1.txt", Qpc1, delimiter=",")
np.savetxt("./out/Qpc2.txt", Qpc2, delimiter=",")

