import math
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from scipy.sparse.linalg import svds
from comparison import compute_all_metrics

# --- Configuration ---
NPC = 5                                      # fixed
OVERSAMPLE_VALUES = [10, 14, 18, 22, 26, 30]
POWER_ITER_VALUES = [20, 40, 60, 80, 100]

# --- Load input data (shared across all sweep runs) ---
X1 = pd.read_csv("./out/geno_pca_party1.tsv", sep="\t").to_numpy(dtype=float)
X2 = pd.read_csv("./out/geno_pca_party2.tsv", sep="\t").to_numpy(dtype=float)
nsample1, nsnp = X1.shape
nsample2 = X2.shape[0]
nsample = nsample1 + nsample2
print(f"X1 shape {X1.shape}")
print(f"X2 shape {X2.shape}")

X = pd.read_csv("./out/geno_pca_merged.tsv", sep="\t").to_numpy(dtype=float)
mean = X.mean(axis=0)
x2mean = (X ** 2).mean(axis=0)
var = x2mean - mean ** 2
var = np.where(var < 1e-8, 1.0, var)
stdinv = 1.0 / np.sqrt(var)

np.savetxt("./out/Xmean.tsv", mean, delimiter="\t")
np.savetxt("./out/XStdInv.tsv", stdinv, delimiter="\t")

num_snp_sqrt_inv = 1.0 / math.sqrt(nsnp)
num_tot_ind_sqrt_inv = 1.0 / math.sqrt(nsample)

# --- Cached SF-GWAS sketch was generated for KP = NPC + 10 = 15 ---
CACHED_KP = 15
buckets1_cached = np.loadtxt("./cache/party1/SketchBucketId.txt").astype(int)
sgn1_cached     = np.loadtxt("./cache/party1/SketchSign.txt")
buckets2_cached = np.loadtxt("./cache/party2/SketchBucketId.txt").astype(int)
sgn2_cached     = np.loadtxt("./cache/party2/SketchSign.txt")
all_buckets_cached = np.concatenate([buckets1_cached, buckets2_cached])
all_sgn_cached     = np.concatenate([sgn1_cached, sgn2_cached])
Q_sketch_cached_raw = np.loadtxt("./cache/party1/Sketch.txt", delimiter=",")


def matmult_lazynorm(Q, Xmat, mean, stdinv, is_col):
    """
    Multiply Q by normalized X without materializing X_norm explicitly.
    Mirrors QXLazyNormStream (is_col=False) and QXtLazyNormStream (is_col=True) in SF-GWAS.
    """
    if is_col:
        D = Q @ Xmat
        correction = Q.sum(axis=1, keepdims=True) @ mean.reshape(1, -1)
        E = D - correction
        R = E * stdinv.reshape(1, -1)
        return R
    else:
        D = Q * stdinv.reshape(1, -1)
        E = D @ Xmat
        correction = (D @ mean.reshape(-1, 1)) @ np.ones((1, Xmat.shape[1]))
        R = E - correction
        return R


def build_sketch(kp, rng):
    """
    Build the normalized random-projection sketch matrix Q_sketch (kp x nsnp).

    If kp matches the cached SF-GWAS sketch dimension (KP=15), reuse the
    cached cipher-domain sketch exactly. Otherwise, regenerate a fresh
    count-sketch-style random projection in plaintext, using the same
    bucket/sign/normalization scheme as SF-GWAS, since the cached cipher
    sketch was only ever generated for KP=15.
    """
    if kp == CACHED_KP:
        Q_sketch = Q_sketch_cached_raw.copy()
        all_buckets = all_buckets_cached
        all_sgn = all_sgn_cached
    else:
        n_total = nsample
        all_buckets = rng.integers(0, kp, size=n_total)
        all_sgn = rng.choice([-1.0, 1.0], size=n_total)
        sgn_col = all_sgn.reshape(-1, 1)
        Q_sketch = np.zeros((kp, nsnp))
        for b in range(kp):
            mask = all_buckets == b
            if mask.any():
                Q_sketch[b] = (sgn_col[mask] * X[mask]).sum(axis=0)

    bucket_count = np.bincount(all_buckets, minlength=kp)
    pos_count    = np.bincount(all_buckets[all_sgn > 0].astype(int), minlength=kp)

    for b in range(kp):
        bc = bucket_count[b] if bucket_count[b] > 0 else 1  # guard div-by-zero
        mean_weight = float(2 * pos_count[b] - bucket_count[b])
        Q_sketch[b] = (Q_sketch[b] - mean_weight * mean) / math.sqrt(bc) * stdinv

    return Q_sketch


def run_rpca(num_oversample, n_power_iter, rng):
    kp = NPC + num_oversample
    Q_sketch = build_sketch(kp, rng)

    Q_init_x = matmult_lazynorm(Q_sketch, X.T, mean, stdinv, is_col=False) * num_snp_sqrt_inv
    Q, _ = np.linalg.qr(Q_init_x.T)
    Q = Q.T

    for it in range(n_power_iter):
        Q_tmp = matmult_lazynorm(Q, X, mean, stdinv, is_col=True) * num_tot_ind_sqrt_inv
        Q_tmp = matmult_lazynorm(Q_tmp, X.T, mean, stdinv, is_col=False) * num_snp_sqrt_inv
        if it < n_power_iter - 1:
            Q, _ = np.linalg.qr(Q_tmp.T)
            Q = Q.T
        else:
            Q = Q_tmp  # last iteration: skip QR (mirrors SF-GWAS)

    Z = (Q @ Q.T) / nsample
    eigenvalues, eigenvectors = np.linalg.eigh(Z)
    idx = np.argsort(eigenvalues)[::-1]
    V = eigenvectors[:, idx[:NPC]].T

    Qpc1 = V @ Q[:, :nsample1]
    Qpc2 = V @ Q[:, nsample1:]
    return Qpc1, Qpc2


def unit_normalize_cols(M):
    norms = np.linalg.norm(M, axis=0, keepdims=True)
    norms = np.where(norms < 1e-8, 1.0, norms)
    return M / norms

def align_signs(ref, cmp):
    """Flip sign of each column in cmp to match the sign of correlation with ref."""
    aligned = cmp.copy()
    for i in range(ref.shape[1]):
        if np.corrcoef(ref[:, i], cmp[:, i])[0, 1] < 0:
            aligned[:, i] *= -1
    return aligned


# --- Ground truth ---
def compute_ground_truth(X_std, npc):
    U, S, Vt = svds(X_std, k=npc)              # ascending singular values
    print(f"U shape {U.shape}")
    print(f"S shape {S.shape}")
    print(f"Vt shape {Vt.shape}")
    order = np.argsort(S)[::-1]
    U = U[:, order]
    S = S[order]
    Vt = Vt[order, :]
    sample_score = U*S
    np.savetxt(f"./out/sample_pcs_gt.tsv", U, delimiter="\t")
    np.savetxt(f"./out/sample_score_gt.tsv", U*S, delimiter="\t")
    print("Top singular values:", S)
    return U, U*S

X_std = (X - mean) * stdinv
print(f"X_std shape {X_std.shape}")
sample_eigenvec_gt, sample_score_gt = compute_ground_truth(X_std, NPC)

# Qpc1_gt = sample_eigenvec_gt[:nsample1]
# Qpc2_gt = sample_eigenvec_gt[nsample1:]
Qpc1_gt =sample_score_gt[:nsample1]
Qpc2_gt =sample_score_gt[nsample1:]
print("compare projection on sample PCs")
Qpc1_gt_n = unit_normalize_cols(Qpc1_gt)
Qpc2_gt_n = unit_normalize_cols(Qpc2_gt)


#--- Sweep over NUM_OVERSAMPLE and N_POWER_ITER ---
rng = np.random.default_rng(seed=0)

for num_oversample in OVERSAMPLE_VALUES:
    kp = NPC + num_oversample
    n_power_iter = 20
    print(f"\n=== NPC={NPC}, NUM_OVERSAMPLE={num_oversample} (KP={kp}), "
            f"N_POWER_ITER={20} ===")

    Qpc1, Qpc2 = run_rpca(num_oversample, n_power_iter, rng)
    print(f"Qpc1 shape: {Qpc1.shape}, Qpc2 shape: {Qpc2.shape}")

    Qpc1_n = unit_normalize_cols(align_signs(Qpc1_gt_n, Qpc1.T))
    Qpc2_n = unit_normalize_cols(align_signs(Qpc2_gt_n, Qpc2.T))
    # Qpc1_n = unit_normalize_cols(Qpc1.T)
    # Qpc2_n = unit_normalize_cols(Qpc2.T)

    study_id_1 = f"qpc1_vs_gt_npc{NPC}_os{num_oversample}_kp{kp}_iter{n_power_iter}"
    study_id_2 = f"qpc2_vs_gt_npc{NPC}_os{num_oversample}_kp{kp}_iter{n_power_iter}"
    outdir_1 = f"out/metrics/qpc1_vs_gt/"
    outdir_2 = f"out/metrics/qpc2_vs_gt/"

    compute_all_metrics(Qpc1_gt_n, Qpc1_n, study_id=study_id_1, outfile=outdir_1)
    compute_all_metrics(Qpc2_gt_n, Qpc2_n, study_id=study_id_2, outfile=outdir_2)

for n_power_iter in POWER_ITER_VALUES:
    kp = CACHED_KP
    num_oversample = 10
    print(f"\n=== NPC={NPC}, NUM_OVERSAMPLE={num_oversample} (KP={kp}), "
            f"N_POWER_ITER={20} ===")

    Qpc1, Qpc2 = run_rpca(num_oversample, n_power_iter, rng)
    print(f"Qpc1 shape: {Qpc1.shape}, Qpc2 shape: {Qpc2.shape}")

    Qpc1_n = unit_normalize_cols(align_signs(Qpc1_gt_n, Qpc1.T))
    Qpc2_n = unit_normalize_cols(align_signs(Qpc2_gt_n, Qpc2.T))
    # Qpc1_n = unit_normalize_cols(Qpc1.T)
    # Qpc2_n = unit_normalize_cols(Qpc2.T)

    study_id_1 = f"qpc1_vs_gt_npc{NPC}_os{num_oversample}_kp{kp}_iter{n_power_iter}"
    study_id_2 = f"qpc2_vs_gt_npc{NPC}_os{num_oversample}_kp{kp}_iter{n_power_iter}"
    outdir_1 = f"out/metrics/qpc1_vs_gt/"
    outdir_2 = f"out/metrics/qpc2_vs_gt/"

    compute_all_metrics(Qpc1_gt_n, Qpc1_n, study_id=study_id_1, outfile=outdir_1)
    compute_all_metrics(Qpc2_gt_n, Qpc2_n, study_id=study_id_2, outfile=outdir_2)
