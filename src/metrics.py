import numpy as np

def compute_zp_from_cov(cov_matrix):
    """
    Computes the pivot redshift from a 2x2 covariance matrix of (w0, wa).
    
    Parameters:
    - cov_matrix: 2D numpy array shaped (2, 2).
                  Assumes index 0 is w0 and index 1 is wa.
                  
    Returns:
    - zp: The pivot redshift.
    """
    # Extract the variance of wa and the covariance of w0, wa
    var_wa = cov_matrix[1, 1]
    cov_w0_wa = cov_matrix[0, 1]
    
    # Calculate pivot scale factor
    ap = 1.0 + (cov_w0_wa / var_wa)
    
    # Convert scale factor to redshift
    zp = (1.0 / ap) - 1.0
    
    return zp