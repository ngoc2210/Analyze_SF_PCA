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
    for i in range(min(reported, min(canonical.shape[1], split.shape[1]))):
        a = angle(canonical[:, i], split[:, i])
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
    for i in range(min(reported, min(canonical.shape[1], split.shape[1]))):
        c = np.corrcoef(canonical[:, i], split[:, i])
        correlations.append(c[0,1])
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
    for i in range(1, eigenvectors.shape[1]+1):
        G = eigenvectors[:, :i]

        proj = G.T @ data
        rec = G @ proj
        res.append(np.linalg.norm(data - rec) / (data.shape[1] * data.shape[0]))
    return res

def subspace_reconstruction_error(data, eigenvectors):
    '''
    returns the average elementwise subspace reconstruction error
    Args:
        data:
        eigenvectors:

    Returns:

    '''
    res = []
    for i in range(1, eigenvectors.shape[1]+1):
        G = eigenvectors[:, :i]

        proj = G.T @ data
        rec = G @ proj

        res.append(np.round(np.linalg.norm(data - rec), 2))
    return res

def mev(u, truth):
    k = min(truth.shape[1], u.shape[1])  # number of eigenvectors in subspace
    m = np.dot(u.T, truth)
    total = 0
    for i in range(k):
        total = total + np.linalg.norm(m[:, i], 2)
    mev = total / k
    return mev

# def angle360(v1, v2):
#     """
#     Calculates the angle between to n-dimensional vectors
#     and returns the result in degree.
#     Args:
#         v1: first vector
#         v2: second vector

#     Returns: angle in degree or NaN

#     """
#     dp = np.dot(v1, v2)
#     norm_x = la.norm(v1)
#     norm_y = la.norm(v2)
#     co = np.clip(dp / (norm_x * norm_y), -1, 1)
#     theta = np.arccos(co)
#     a = math.degrees(theta)
#     # angle can be Nan
#     # print(angle)
#     if math.isnan(a):
#         return a
#     # return the canonical angle
#     return a

def mse(v1, v2):
    """
    Computes Mean Squared Error between two vectors
    """
    diff = v1 - v2
    return np.mean(diff ** 2)

def compute_mses(V1, V2, reported=20):
    """
    Computes MSE for matching column vectors in two matrices
    """
    res = []
    k = min(reported, min(V1.shape[1], V2.shape[1]))

    for i in range(k):
        res.append(mse(V1[:, i], V2[:, i]))

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


def compute_save_subspace_reconstruction_error(data, eigenvectors, study_id, filename, outfile):
    errors = subspace_reconstruction_error(data, eigenvectors)
    with open(path.join(outfile, filename), 'a+') as handle:
        handle.write(cv.collapse_array_to_string(errors, study_id=study_id))

def compute_save_mev(u, truth, study_id, filename, outfile):
    value = mev(u, truth)
    with open(path.join(outfile, filename), 'a+') as handle:
        handle.write(cv.collapse_array_to_string([value], study_id=study_id))

geno = pd.read_csv("out/geno_pca_merged.tsv", sep="\t")
data = geno.values.astype(np.float64)

G_ref = pd.read_csv("out/G_sample_eigenvectors.tsv", sep="\t")
eigenvectors_ref = G_ref.iloc[:,:5].values

Qpc1 = np.loadtxt("cache/party1/Qpc.txt", delimiter=",")
Qpc2 = np.loadtxt("cache/party2/Qpc.txt", delimiter=",")
Q = np.hstack([Qpc1, Qpc2])
eigenvectors_sf = Q.T

outfile = "out/metrics"
os.makedirs(outfile, exist_ok=True)

study_id = "comparison_ref_vs_sf"

compute_save_subspace_reconstruction_error_element_wise(
    data, eigenvectors_ref, study_id, "rec_err_elem_ref.txt", outfile
)

compute_save_subspace_reconstruction_error_element_wise(
    data, eigenvectors_sf, study_id, "rec_err_elem_sf.txt", outfile
)

compute_save_subspace_reconstruction_error(
    data, eigenvectors_ref, study_id, "rec_err_ref.txt", outfile
)

compute_save_subspace_reconstruction_error(
    data, eigenvectors_sf, study_id, "rec_err_sf.txt", outfile
)

# Angles
compute_save_angles(
    eigenvectors_ref, eigenvectors_sf,
    study_id, "angles_ref_vs_sf.txt", outfile
)

# Correlations
compute_save_correlations(
    eigenvectors_ref, eigenvectors_sf,
    study_id, "correlations_ref_vs_sf.txt", outfile
)

# Euclidean distances
compute_save_euclidean_distance(
    eigenvectors_ref, eigenvectors_sf,
    study_id, "euclidean_ref_vs_sf.txt", outfile
)

# MSE (column-wise)
compute_save_mses(
    eigenvectors_ref, eigenvectors_sf,
    study_id, "mses_ref_vs_sf.txt", outfile
)

# MEV (both directions, since it is asymmetric)
compute_save_mev(
    eigenvectors_ref, eigenvectors_sf,
    study_id, "mev_ref_vs_sf.txt", outfile
)

compute_save_mev(
    eigenvectors_sf, eigenvectors_ref,
    study_id, "mev_sf_vs_ref.txt", outfile
)

