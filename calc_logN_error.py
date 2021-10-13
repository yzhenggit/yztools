def calc_logN_error(logN, elogN): 
    import numpy as np 
    logN = np.asarray(logN)
    elogN = np.asarray(elogN)
    sum_logN = np.log10(np.sum(10**logN))
    
    top = np.sqrt(np.sum((elogN*10**logN*np.log(10))**2))
    bottom = np.sum(10**logN)*np.log(10)
    
    e_sum_logN = top/bottom 
    return sum_logN, e_sum_logN
