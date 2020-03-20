# v{\rm GSR} - 62 \cos(l) \cos(b) + 40 \sin(l) \cos(b) - 35 \sin(b)
def vgsr2vlgsr(vgsr, l_deg, b_deg, do_print=False):
    """
    How to use:
    vgsr2vlgsr(vgsr, l_deg, b_deg)

    History:
    03/16/2020, YZ, UCB
    """

    import numpy as np
    l_rad = np.radians(l_deg)
    b_rad = np.radians(b_deg)
    vlgsr = vgsr - 62*np.cos(l_rad)*np.cos(b_rad) + \
            40*np.sin(l_rad)*np.cos(b_rad) - 35*np.sin(b_rad)
    if do_print == True:
        print(">>> Input: vgsr=%.2f km/s, l=%.4f deg, b=%.4f deg"%(vgsr, l_deg, b_deg))
        print(">>> Output: vlgsr=%.2f km/s"%(vlgsr))

    return vlgsr


if __name__ == '__main__':
    import sys
    import numpy as np

    vgsr = np.float(sys.argv[1])
    l_deg = np.float(sys.argv[2])
    b_deg = np.float(sys.argv[3])

    vgsr2vlgsr(vgsr,l_deg, b_deg, do_print=True)
