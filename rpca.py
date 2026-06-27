"""
Randomized PCA (RPCA), following the workflow described in Figure 1 / Section 3
of the SF-PCA paper (Froelicher, Cho et al., "Scalable and Privacy-Preserving
Federated Principal Component Analysis").

This script implements ONLY the algorithmic (cleartext) side of RPCA, i.e. it
ignores everything related to the "SF" (Secure & Federated) part of the paper:
  - no multiparty homomorphic encryption (MHE)
  - no secret-sharing / MPC
  - no splitting of the data across data providers (DPs)

The input matrix A is assumed to already be a single,
pooled, cleartext matrix -- exactly like geno_pca_merged.tsv, which already
contains the merged genotype data from the two parties.

The 8 steps below follow the paper's notation directly (Fig. 1, page 3):
  Step 1: Setup                  -- choose rho = N_PCS + alpha, build sketch Pi
  Step 2: Mean-Centering          -- O = column means of A
  Step 3: Random Projection       -- P = Pi x A
  Step 4: Power Iterations        -- P <- QR( (A^T A) x P ), repeated p times
  Step 5: Reduction                -- Z = P x A^T A x P^T   (small rho x rho matrix)
  Step 6: Eigendecomposition       -- Z = W L W^T
  Step 7: Reconstruction           -- W <- W x P (and re-orthonormalize), keep top N_PCS rows
  Step 8: Projection                -- A' = A x W^T (project samples onto the PCs)

Notes on faithfulness to the paper:
  - Random sketch matrix Pi: the paper uses a "count-sketch" matrix whose
    entries are drawn from {-1, 0, 1} (one nonzero +-1 entry per column).
    We implement this exact sketch (see `count_sketch_matrix`).
  - Mean-centering is treated explicitly via a separate mean vector `o`
    rather than by mutating A in place, mirroring the paper's "lazy
    mean-centering" idea (Step 2 in Fig. 1) even though, without
    encryption, there's no efficiency reason to delay it -- it's done here
    purely to mirror the paper's algorithmic structure.
  - STANDARDIZATION: SF-PCA's Fig. 1 Step 2 only mean-centers. However,
    for genotype/SNP data, this implementation additionally standardizes
    each feature column by its inverse standard deviation, following the
    companion SF-GWAS paper's Protocol 4 ("PCA-based GWAS"), lines 1-9:
    there, the parties jointly compute the column mean `m` and inverse
    standard deviation `s` of the genotype matrix, then mean-center AND
    multiply by `s` ("Step 9: Multiply each row of P with s for
    standardization") before using it in the random projection and power
    iterations. This is the standard genotype-PCA convention (cf. Price
    et al. 2006, "Principal components analysis corrects for
    stratification in genome-wide association studies") and matters here
    because rare SNPs have small variance and common SNPs have large
    variance -- without standardizing, PCA would be dominated by allele
    frequency/scale rather than genuine population structure. Pass
    `standardize=False` to recover the original SF-PCA Fig. 1 behavior
    (mean-centering only, no scaling).
  - QR re-orthogonalization (the "QR^T" box in Fig. 1) is applied to the
    rows of the sketch matrix P after each power iteration for numerical
    stability, exactly as described after Step 4 in the paper.
  - Step 5 forms the small (rho x rho) matrix Z = P (A^T A) P^T directly
    (i.e., it reuses the already mean-corrected covariance action).
  - Step 6 uses a standard symmetric eigendecomposition (the paper uses a
    QR-iteration-based eigensolver under encryption; in cleartext we can
    just call a numerically stable library routine, which solves the same
    mathematical problem on the same small matrix Z).
  - Step 7 reconstructs eigenvectors in the original feature space and
    keeps only the top `N_PCS` (=ψ) principal components, discarding the
    `alpha` (=α) oversampling components, as described in the text
    following Eigen in Step 6/7 of Fig. 1.
"""

import numpy as np



# random sketch matrix (count-sketch, entries in {-1, 0, 1})
def count_sketch_matrix(rho: int, n: int, rng: np.random.Generator) -> np.ndarray:
    """
    Build the public random sketch matrix Pi (rho x n) used in Step 1/3.

    Following the paper's count-sketch construction (Sec. 3, "We adopt the
    count-sketch approach for generating the latter, where the elements are
    drawn from {-1,0,1}"): each column has exactly one nonzero entry, placed
    at a uniformly random row, with a uniformly random sign (+1 or -1).
    """
    Pi = np.zeros((rho, n))
    rows = rng.integers(0, rho, size=n)
    signs = rng.choice([-1.0, 1.0], size=n)
    Pi[rows, np.arange(n)] = signs
    return Pi


def qrT(P: np.ndarray) -> np.ndarray:
    """
    The "QR^T" operation in Fig. 1: QR factorization applied to the ROWS of
    P (not the columns), used to re-orthonormalize the sketch after each
    power iteration. We obtain this by QR-factorizing P^T and transposing
    back, i.e. orthonormalizing the row space of P.
    """
    Q, _ = np.linalg.qr(P.T)   # Q: (cols x rank), orthonormal columns of P^T
    return Q.T                # rows of the returned matrix are orthonormal


def randomized_pca(
    A: np.ndarray,
    N_PCS: int,
    alpha: int = 4,
    p: int = 10,
    seed: int = 0,
    standardize: bool = True,
):
    """
    Cleartext Randomized PCA (RPCA), following Fig. 1 of the SF-PCA paper.

    Parameters
    ----------
    A     : (n_samples x m_features) data matrix (cleartext, pooled).
    N_PCS   : desired number of principal components (N_PCS in the paper).
    alpha : oversampling parameter (alpha in the paper); rho = N_PCS + alpha.
    p     : number of power iterations.
    seed  : RNG seed for the public random sketch matrix Pi.
    standardize : if True (default), also divide each mean-centered feature
                column by its standard deviation (per SF-GWAS Protocol 4,
                lines 1-9), as is standard for genotype/SNP-based PCA. If
                False, only mean-centers (original SF-PCA Fig. 1 behavior).

    Returns
    -------
    W        : (N_PCS x m) principal components (rows = PCs, in descending
                order of explained variance), as in the paper's output W.
    eigvals  : (N_PCS,) corresponding eigenvalues (explained variance), in
                descending order.
    A_proj   : (n x N_PCS) projected data, i.e. A' = A_std x W^T (Step 8), where
                A_std is A after mean-centering (and standardization, if
                `standardize=True`).
    """
    n, m = A.shape
    rho = N_PCS + alpha
    rng = np.random.default_rng(seed)

    Pi = count_sketch_matrix(rho, n, rng)          # (rho x n)

    # ---- Step 2: Mean-Centering (+ Standardization, per SF-GWAS Protocol 4) ----
    o = A.mean(axis=0)                                # (m,) column means ("m" in Protocol 4)
    if standardize:
        # Following SF-GWAS Protocol 4, lines 1-5: compute the column mean
        # `m` and the column variance `v` (= mean of squares - mean^2), then
        # the inverse standard deviation `s = 1/sqrt(v)` (MPC-SqrtInv there).
        var = A.var(axis=0)                            # (m,) per-column variance ("v - m^2")
        # Guard against (near-)monomorphic columns with ~zero variance,
        # which would otherwise blow up after dividing by std.
        std = np.sqrt(np.maximum(var, 1e-12))
        s = 1.0 / std                                   # (m,) inverse std ("s" in Protocol 4)
    else:
        s = np.ones(m)                                   # no-op scaling: recovers Fig. 1 behavior
    # (We keep A in "cleartext" form and apply -o and *s explicitly whenever
    #  needed, mirroring the paper's lazy mean-centering/standardization
    #  bookkeeping -- "on-the-fly standardization", Protocol 5 in SF-GWAS.)

    def standardize_matrix(M):
        """Compute ((M - 1*o^T) * s) for an (a x m) matrix M, i.e. mean-center
        and standardize its rows by the global per-feature mean/std."""
        return (M - o[np.newaxis, :]) * s[np.newaxis, :]

    # ---- Step 3: Random Projection ----
    # P = Pi x standardize(A)   -> (rho x n) x (n x m) = (rho x m)
    A_std = standardize_matrix(A)                     # (n x m), mean-centered & standardized once
    P = Pi @ A_std                                     # (rho x m)

    # ---- Step 4: Power Iterations ----
    # Repeatedly multiply with the covariance action A_std^T A_std to amplify
    # the spectral gap, re-orthonormalizing (QR^T) after each iteration.
    AtA = A_std.T @ A_std                              # (m x m) covariance (up to scale)
    for _ in range(p):
        P = P @ AtA                                   # (rho x m) x (m x m) = (rho x m)
        P = qrT(P)                                     # re-orthonormalize rows of P

    # ---- Step 5: Reduction ----
    # Z = P x (A^T A) x P^T  -> small (rho x rho) symmetric matrix.
    Z = P @ AtA @ P.T                                  # (rho x rho)
    Z = 0.5 * (Z + Z.T)                                 # enforce exact symmetry

    # ---- Step 6: Eigendecomposition ----
    # eigh returns eigenvalues in ascending order; flip to descending order.
    eigvals_all, U_all = np.linalg.eigh(Z)              # Z = U_all diag(eigvals_all) U_all^T
    order = np.argsort(eigvals_all)[::-1]
    eigvals_all = eigvals_all[order]
    U_all = U_all[:, order]                              # (rho x rho), columns = eigenvectors

    # ---- Step 7: Reconstruction ----
    # Reconstruct eigenvectors in the original feature space:
    #   W_full = U_all^T x P   -> (rho x rho) x (rho x m) = (rho x m)
    W_full = U_all.T @ P                                 # (rho x m)
    # Re-orthonormalize for numerical stability (final QR^T in Step 7).
    W_full = qrT(W_full)
    # Keep only the top N_PCS components (discard the alpha oversampling rows).
    W = W_full[:N_PCS, :]                                   # (N_PCS x m)
    eigvals = eigvals_all[:N_PCS]

    # ---- Step 8: Projection ----
    # Project the (mean-centered & standardized) samples onto the PCs:
    # A' = A_std x W^T
    A_proj = A_std @ W.T                                  # (n x N_PCS)

    return W, eigvals, A_proj
