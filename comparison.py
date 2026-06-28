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
        a = angle(canonical[i,:], split[i, :])
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
        c = np.corrcoef(canonical[i, :], split[i, :])
        correlations.append(abs(c[0,1]))
    return correlations


def compute_save_angles(W0, W1, study_id, filename, outfile, reported=5):
    angles = compute_angles(W0, W1, reported=reported)
    with open(path.join(outfile, filename), 'a+') as handle:
        handle.write(cv.collapse_array_to_string(angles, study_id=study_id))

def calculate_euclidean_distance(V1, V2):
    res = []
    for line in range(min(V1.shape[1], V2.shape[1])):
        res.append(d.euclidean(V1[:, line], V2[:, line]))
    return res

def subspace_reconstruction_error_element_wise(data, eigenvectors):
    '''
    returns the average elementwise subspace reconstruction error
    Args:
        data:
        eigenvectors:

    Returns:

    '''
    res = []
    # print(f"shape eigenvectors is {eigenvectors.shape}")
    # print(f"shape data is {data.shape}")
    for i in range(1, eigenvectors.shape[0]+1):
        G = eigenvectors[:i, :]

        proj = data @ G.T
        rec = proj @ G
        res.append(np.linalg.norm(data - rec) / (data.shape[1] * data.shape[0]))
    return res

def mev(u, truth):
    k = min(truth.shape[0], u.shape[0])  # number of eigenvectors in subspace
    m = np.dot(u.T, truth)
    total = 0
    for i in range(k):
        total = total + np.linalg.norm(m[:, i], 2)
    mev = total / k
    return mev

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
    for i in range(min(reported, min(v1.shape[0], v2.shape[0]))):
        a = angle360(v1[i, :], v2[i, :])
        angles.append(np.around(a, 2))
    return angles

def compute_save_angles360(W0, W1, study_id, filename, outfile, reported=5):
    angles = compute_angles360(W0, W1, reported=reported)
    with open(path.join(outfile, filename), 'a+') as handle:
        handle.write(cv.collapse_array_to_string(angles, study_id=study_id))

def mse(v1, v2):
    """
    Computes Mean Squared Error between two vectors
    """
    diff = v1 - v2
    return np.mean(diff ** 2)

def compute_mses(V1, V2, reported=5):
    """
    Computes MSE for matching column vectors in two matrices
    """
    res = []
    k = min(reported, min(V1.shape[0], V2.shape[0]))

    for i in range(k):
        res.append(mse(V1[i, :], V2[i, :]))

    return res

def compute_save_mses(V1, V2, study_id, filename, outfile, reported=5):
    mses = compute_mses(V1, V2, reported=reported)
    with open(path.join(outfile, filename), 'a+') as handle:
        handle.write(cv.collapse_array_to_string(mses, study_id=study_id))

def compute_save_correlations(canonical, split, study_id, filename, outfile, reported=5):
    correlations = compute_correlations(canonical, split, reported=reported)
    with open(path.join(outfile, filename), 'a+') as handle:
        handle.write(cv.collapse_array_to_string(correlations, study_id=study_id))

def compute_save_euclidean_distance(V1, V2, study_id, filename, outfile):
    distances = calculate_euclidean_distance(V1, V2)
    with open(path.join(outfile, filename), 'a+') as handle:
        handle.write(cv.collapse_array_to_string(distances, study_id=study_id))

def compute_save_subspace_reconstruction_error_element_wise(data, eigenvectors, study_id, filename, outfile):
    errors = subspace_reconstruction_error_element_wise(data, eigenvectors)
    with open(path.join(outfile, filename), 'a+') as handle:
        handle.write(cv.collapse_array_to_string(errors, study_id=study_id))

def compute_save_mev(u, truth, study_id, filename, outfile):
    value = mev(u, truth)
    with open(path.join(outfile, filename), 'a+') as handle:
        handle.write(cv.collapse_array_to_string([value], study_id=study_id))

def compute_all_metrics(data, eigenvectors_ref, eigenvectors_cmp, study_id, outfile):
    os.makedirs(outfile, exist_ok=True)
    compute_save_subspace_reconstruction_error_element_wise(data, eigenvectors_cmp, study_id, "rec_err_elem.txt", outfile)
    compute_save_angles(eigenvectors_ref, eigenvectors_cmp, study_id, "angles.txt", outfile)
    compute_save_angles360(eigenvectors_ref, eigenvectors_cmp, study_id, "angles360.txt", outfile)
    compute_save_correlations(eigenvectors_ref, eigenvectors_cmp, study_id, "correlations.txt", outfile)
    compute_save_euclidean_distance(eigenvectors_ref, eigenvectors_cmp, study_id, "euclidean.txt", outfile)
    compute_save_mses(eigenvectors_ref, eigenvectors_cmp, study_id, "mses.txt", outfile)
    compute_save_mev(eigenvectors_ref, eigenvectors_cmp, study_id, "mev.txt", outfile)
    compute_save_mev(eigenvectors_cmp, eigenvectors_ref, study_id, "mev_inverse.txt", outfile)

geno = pd.read_csv("out/geno_pca_merged.tsv", sep="\t")
data = geno.values.astype(np.float64)
data = data.T
print(f"data shape {data.shape}")

eigenvectors_ref = np.loadtxt("out/sample_eigenvectors.tsv", delimiter="\t")
eigenvectors_ref = eigenvectors_ref[:5]
print(f"eigenvectors_ref shape {eigenvectors_ref.shape}")

Qpc1 = np.loadtxt("cache/party1/Qpc1.txt", delimiter=",")
Qpc2 = np.loadtxt("cache/party2/Qpc2.txt", delimiter=",")
eigenvectors_sf = np.hstack([Qpc1, Qpc2])
print(f"eigenvectors_sf shape {eigenvectors_sf.shape}")

col_norms = np.linalg.norm(eigenvectors_sf, axis=0, keepdims=True)
col_norms = np.where(col_norms < 1e-10, 1.0, col_norms)  # guard against zero norm
eigenvectors_sf = eigenvectors_sf / col_norms

compute_all_metrics(data, eigenvectors_ref, eigenvectors_sf, study_id="comparison_ref_vs_sf", outfile="out/metrics/ref_vs_sf")

eigenvectors_rpca = np.loadtxt("out/rpca_sample_eigenvectors.tsv", delimiter="\t")
print(f"eigenvectors_rpca shape {eigenvectors_rpca.shape}")
study_id_rpca= "comparison_ref_vs_rpca"
outfile = "out/metrics/ref_vs_rpca"

compute_all_metrics(data, eigenvectors_ref, eigenvectors_rpca, study_id="comparison_ref_vs_rpca", outfile="out/metrics/ref_vs_rpca")

compute_all_metrics(data, eigenvectors_rpca, eigenvectors_sf, study_id="comparison_rpca_vs_sf", outfile="out/metrics/rpca_vs_sf")
