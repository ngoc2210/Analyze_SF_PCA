import numpy as np
import scipy.linalg as la
import math as math
import os.path as path
import convenience as cv
import scipy.spatial.distance as d
import pandas as pd
import os



def angle(v1, v2):
    """
    Calculates the angle between to n-dimensional vectors
    and returns the result in degree. The angle returned
    can maximally be 90, otherwise the complement will be returned
    Args:
        v1: first vector
        v2: second vector

    Returns: angle in degree or NaN

    """
    dp = np.dot(v1, v2)
    norm_x = la.norm(v1)
    norm_y = la.norm(v2)
    co = np.clip(dp / (norm_x * norm_y), -1, 1)
    theta = np.arccos(co)
    a = math.degrees(theta)
    # angle can be Nan
    # print(angle)
    if math.isnan(a):
        return a
    # return the canonical angle
    if a > 90:
        return np.abs(a - 180)
    else:
        return a

def angle360(v1, v2):
    """
    Calculates the angle between to n-dimensional vectors
    and returns the result in degree.
    Args:
        v1: first vector
        v2: second vector

    Returns: angle in degree or NaN

    """
    dp = np.dot(v1, v2)
    norm_x = la.norm(v1)
    norm_y = la.norm(v2)
    co = np.clip(dp / (norm_x * norm_y), -1, 1)
    theta = np.arccos(co)
    a = math.degrees(theta)
    # angle can be Nan
    # print(angle)
    if math.isnan(a):
        return a
    # return the canonical angle
    return a


def euclidean_distance(V1, V2):
    res = []
    for i in range(5):
        res.append(d.euclidean(V1[:, i], V2[:, i]))
    return res

def mse(v1, v2):
    """
    Computes Mean Squared Error between two vectors
    """
    diff = v1 - v2
    return np.mean(diff ** 2)

def correlation(v1, v2):
    return abs(np.corrcoef(v1, v2)[0, 1])


def euclidean(v1, v2):
    return d.euclidean(v1, v2)


def compute_angles(canonical, split, reported=5):
    """
    Compute the angles of the vectors in a matrix with matching vectors
    in a second matrix.

    Args:
        canonical: The first matrix of eigenvectors
        split:  The second matrix of eigenvectors
        reported_angles: Limit for number of angles to compute

    Returns: A vector of angles

    """
    angles = list()
    for i in range(reported):
        a = angle(canonical[:, i], split[:, i])
        angles.append(np.around(a, 2))
    return angles

def compute_angles360(v1, v2, reported=5):
    """
    Compute the angles of the vectors in a matrix with matching vectors
    in a second matrix.

    Args:
        canonical: The first matrix of eigenvectors
        split:  The second matrix of eigenvectors
        reported_angles: Limit for number of angles to compute

    Returns: A vector of angles

    """
    angles = list()
    for i in range(reported):
        a = angle360(v1[:, i], v2[:, i])
        angles.append(np.around(a, 2))
    return angles

def compute_correlations(canonical, split, reported=5):
    """
        Compute the correlations of the vectors in a matrix with matching vectors
        in a second matrix.

        Args:
            canonical: The first matrix of eigenvectors
            split:  The second matrix of eigenvectors
            reported_angles: Limit for number of correlations to compute

        Returns: A vector of correlations

        """
    correlations = list()
    for i in range(reported):
        c = np.corrcoef(canonical[:, i], split[:, i])
        correlations.append(abs(c[0,1]))
    return correlations

def compute_mses(V1, V2, reported=5):
    """
    Computes MSE for matching column vectors in two matrices
    """
    res = []
    for i in range(reported):
        res.append(mse(V1[:, i], V2[:, i]))
    return res

def compute_save_angles360(W0, W1, study_id, filename, outfile, reported=5):
    angles = compute_angles360(W0, W1, reported=reported)
    with open(path.join(outfile, filename), 'a+') as handle:
        handle.write(cv.collapse_array_to_string(angles, study_id=study_id))

def compute_save_angles(W0, W1, study_id, filename, outfile, reported=5):
    angles = compute_angles(W0, W1, reported=reported)
    with open(path.join(outfile, filename), 'a+') as handle:
        handle.write(cv.collapse_array_to_string(angles, study_id=study_id))

def compute_save_mses(V1, V2, study_id, filename, outfile, reported=5):
    mses = compute_mses(V1, V2, reported=reported)
    with open(path.join(outfile, filename), 'a+') as handle:
        handle.write(cv.collapse_array_to_string(mses, study_id=study_id))

def compute_save_correlations(canonical, split, study_id, filename, outfile, reported=5):
    correlations = compute_correlations(canonical, split, reported=reported)
    with open(path.join(outfile, filename), 'a+') as handle:
        handle.write(cv.collapse_array_to_string(correlations, study_id=study_id))

def compute_save_euclidean_distance(V1, V2, study_id, filename, outfile):
    distances = euclidean_distance(V1, V2)
    with open(path.join(outfile, filename), 'a+') as handle:
        handle.write(cv.collapse_array_to_string(distances, study_id=study_id))

def compute_all_metrics(ref, cmp, study_id, outfile):
    os.makedirs(outfile, exist_ok=True)
    compute_save_angles(ref, cmp, study_id, "angles.txt", outfile)
    compute_save_angles360(ref, cmp, study_id, "angles360.txt", outfile)
    compute_save_correlations(ref, cmp, study_id, "correlations.txt", outfile)
    compute_save_euclidean_distance(ref, cmp, study_id, "euclidean.txt", outfile)
    compute_save_mses(ref, cmp, study_id, "mses.txt", outfile)

# ### reference A_proj on sample eigenvec
# A_proj_sample_ref = np.loadtxt("out/A_proj_on_sample_pcs_gt.tsv", delimiter="\t")
# A_proj_sample_ref = A_proj_sample_ref[:, :5]
# print(f"A_proj_sample_gt shape {A_proj_sample_ref.shape}")
# ### reference A_proj on feature eigenvec
# A_proj_feature_ref = np.loadtxt("out/A_proj_on_feature_pcs_gt.tsv", delimiter="\t")
# A_proj_feature_ref = A_proj_feature_ref[:, :5]
# print(f"A_proj_feature_gt shape {A_proj_feature_ref.shape}")

# ### reference sample eigenvec
# sample_eigenvec_ref = np.loadtxt("out/sample_eigenvectors.tsv", delimiter="\t")
# sample_eigenvec_ref = sample_eigenvec_ref[:5].T
# print(f"sample_eigenvec_gt shape {sample_eigenvec_ref.shape}")
# ### reference feature eigenvec
# feature_eigenvec_ref = np.loadtxt("out/feature_eigenvectors.tsv", delimiter="\t")
# feature_eigenvec_ref = feature_eigenvec_ref[:5].T
# print(f"feature_eigenvec_gt shape {feature_eigenvec_ref.shape}")

# ### rpca A_proj on sample eigenvec
# A_proj_sample_rpca = np.loadtxt("out/A_proj_on_sample_pcs_rpca.tsv", delimiter="\t")
# print(f"A_proj_sample_rpca shape {A_proj_sample_rpca.shape}")
# ### rpca A_proj on feature eigenvec
# A_proj_feature_rpca = np.loadtxt("out/A_proj_on_feature_pcs_rpca.tsv", delimiter="\t")
# print(f"A_proj_feature_rpca shape {A_proj_feature_rpca.shape}")

# ### rpca sample eigenvec
# sample_eigenvec_rpca = np.loadtxt("out/rpca_sample_eigenvectors.tsv", delimiter="\t")
# sample_eigenvec_rpca = sample_eigenvec_rpca.T
# print(f"sample_eigenvec_rpca shape {sample_eigenvec_rpca.shape}")
# ### rpca feature eigenvec
# feature_eigenvec_rpca = np.loadtxt("out/rpca_feature_eigenvectors.tsv", delimiter="\t")
# feature_eigenvec_rpca = feature_eigenvec_rpca.T
# print(f"feature_eigenvec_rpca shape {feature_eigenvec_rpca.shape}")

# compute_all_metrics(A_proj_sample_ref, A_proj_sample_rpca, study_id="comparison_A_proj_sample_ref_vs_rpca", outfile="out/metrics/A_proj_sample_ref_vs_rpca")
# compute_all_metrics(A_proj_feature_ref, A_proj_feature_rpca, study_id="comparison_A_proj_feature_ref_vs_rpca", outfile="out/metrics/A_proj_feature_ref_vs_rpca")

# compute_all_metrics(sample_eigenvec_ref, sample_eigenvec_rpca, study_id="comparison_sample_eigenvec_ref_vs_rpca", outfile="out/metrics/sample_eigenvec_ref_vs_rpca")
# compute_all_metrics(feature_eigenvec_ref, feature_eigenvec_rpca, study_id="comparison_feature_eigenvec_ref_vs_rpca", outfile="out/metrics/feature_eigenvec_ref_vs_rpca")

### Qpc with and without MHE, MPC
qpc1_crypto_sf = np.loadtxt("out/Qpc1.txt", delimiter=",")
qpc1_no_crypto_sf = np.loadtxt("cache/party1/Qpc1.txt", delimiter=",")
qpc2_crypto_sf = np.loadtxt("out/Qpc2.txt", delimiter=",")
qpc2_no_crypto_sf = np.loadtxt("cache/party2/Qpc2.txt", delimiter=",")
compute_all_metrics(qpc1_no_crypto_sf.T, qpc1_crypto_sf.T, study_id="comparison_qpc1_crypto_vs_no_crypto", outfile="out/metrics/qpc1_crypto_vs_no_crypto")
compute_all_metrics(qpc2_no_crypto_sf.T, qpc2_crypto_sf.T, study_id="comparison_qpc2_crypto_vs_no_crypto", outfile="out/metrics/qpc2_crypto_vs_no_crypto")


### mean and standard deviation
mean_sf = np.loadtxt("cache/party1/Xmean.txt", delimiter=",")[:10010]
mean = np.loadtxt("out/Xmean.tsv", delimiter="\t")
stdinv_sf = np.loadtxt("cache/party1/XStdInv.txt", delimiter=",")[:10010]
stdinv = np.loadtxt("out/XStdInv.tsv", delimiter="\t")

def compute_vector_metrics(v1, v2):
    return {
        "angle": angle(v1, v2),
        "angle360": angle360(v1, v2),
        "correlation": correlation(v1, v2),
        "mse": mse(v1, v2),
        "euclidean": euclidean(v1, v2),
    }

def compare_mean_stdinv(mean_a, mean_b, stdinv_a, stdinv_b):
    return {
        "mean": compute_vector_metrics(mean_a, mean_b),
        "stdinv": compute_vector_metrics(stdinv_a, stdinv_b),
    }

results = compare_mean_stdinv(mean, mean_sf, stdinv, stdinv_sf)

print(results["mean"])
print(results["stdinv"])

