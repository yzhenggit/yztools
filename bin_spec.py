###############################
# bin spectra 

def bin_spec(x1, y1, bins):
    import numpy as np

    b1 = np.mgrid[0:len(x1):bins]
    dg1 = np.digitize(np.mgrid[0:len(x1):1], b1)

    dg_x1 = np.array([np.mean(x1[dg1==j]) for j in np.arange(len(b1)+1)[1:]])
    dg_y1 = np.array([np.mean(y1[dg1==j]) for j in np.arange(len(b1)+1)[1:]])
    
    return dg_x1, dg_y1

def bin_spec_3Dcube(x1, y1, bins, ax=0):
    
    import numpy as np

    b1 = np.mgrid[0:len(x1):bins]
    dg1 = np.digitize(np.mgrid[0:len(x1):1], b1)

    dg_x1 = np.array([np.mean(x1[dg1==j], axis=ax) for j in np.arange(len(b1)+1)[1:]])
    dg_y1 = np.array([np.mean(y1[dg1==j], axis=ax) for j in np.arange(len(b1)+1)[1:]])

    return dg_x1, dg_y1
