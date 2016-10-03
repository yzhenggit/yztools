############################### 
# Gaussian Fitting to a line 
def gauss(x, A=1, mu=1, sigma=1):
    return np.real(A * np.exp(-(x-mu)**2/(2*sigma**2)))

def fit_direct(x, y, F=0, weighted=True, _weights=None):
    import numpy as np
    import scipy as sp
    import scipy.linalg as sl

    mask = (y > F)
    x = x[mask]
    y = y[mask]
    if _weights is None:
        _weights = y
    else:
        _weights = _weights[mask]
    # We do not want to risk working with negative values
    np.clip(y, 1e-10, np.inf, y)
    e = np.ones(len(x))
    if weighted:
        e = e * (_weights**2)
    v = (np.sum(np.vander(x, 5) * e[:, None], axis=0))[::-1]
    A = v[sl.hankel([0, 1, 2], [2, 3, 4])]
    ly = e * np.log(y)
    ls = np.sum(ly)
    x_ls = np.sum(ly*x)
    xx_ls = np.sum(ly*x**2)
    B = np.array([ls, x_ls, xx_ls])
    (a, b, c), res, rank, s=np.linalg.lstsq(A, B)
    A = np.exp(a-(b**2/(4*c)))
    mu = -b/(2*c)
    sigma = sp.sqrt(-1/(2*c))
    return A, mu, sigma

def fit_iterative(x, y, F=0, weighted=True, N=10):
    y_ = y
    for i in range(N):
        p = fit_direct(x, y, weighted=True, _weights=y_)
        A, mu, sigma = p
        y_ = gauss(x, A, mu, sigma)
    return np.real(A), np.real(mu), np.real(sigma)

