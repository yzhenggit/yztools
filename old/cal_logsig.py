##############################
# derive the sigma of the sum or mean of the log values through error propogation
def cal_logsig(val, sig, getsum=True):
    import numpy as np

    val = np.asarray(val)
    sig = np.asarray(sig)

    totval = (10**val).sum()
    totsig = np.sqrt(((10**val*sig)**2).sum()) / totval
    if getsum is True:
        return np.log10(totval), totsig
    else:
        return np.log10(totval/val.size), totsig

