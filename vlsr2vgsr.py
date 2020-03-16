def vlsr2vgsr(vlsr, l_deg, b_deg):
    """
    How to use:
    vlsr2vgsr(vlsr, l_deg, b_deg)

    History:
    03/16/2020, YZ, UCB
    """

    import numpy as np
    l_rad = np.radians(l_deg)
    b_rad = np.radians(b_deg)
    vgsr = vlsr + 220*np.sin(l_rad)*np.cos(b_rad)

    print(">>> Input: vlsr=%.2f km/s, l=%.4f deg, b=%.4f deg"%(vlsr, l_deg, b_deg))
    print(">>> Output: vgsr=%.2f km/s"%(vgsr))

    return vgsr

if __name__ == '__main__':
    import sys
    import numpy as np
    vlsr = np.float(sys.argv[1])
    l_deg = np.float(sys.argv[2])
    b_deg = np.float(sys.argv[3])

    vlsr2vgsr(vlsr, l_deg, b_deg)
