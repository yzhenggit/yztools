from __future__ import print_function

##############
def vhelio2vlsr_Westmeier(vel_init, l_deg, b_deg, reverse=False):
    '''
    - from http://www.atnf.csiro.au/people/Tobias.Westmeier/tools_hihelpers.php
    - l_deg:  should be in degree
    - b_deg: should be in degree
    - vel_init: velocity that need to be transformed
    -          vel_init = vhelio if reverse = False
    -          vel_init = vlsr   if reverse = True
    - Favor this one than the other one from Rosolowsky

    History:
    YZ. created a while ago, use to have ra, dec conversion as well
    03/16/2020, now only take l_deg, b_deg to avoid confusion

    How to use:
    vhelio2vlsr_Westmeier(vlsr, l_deg, b_deg, reverse=True)
    vhelio2vlsr_Westmeier(vhelio, l_deg, b_deg, reverse=False)
    '''
    import numpy as np

    l = np.radians(l_deg)
    b = np.radians(b_deg)
    # vlsr 00> vhelio
    if reverse == True:
        delv = -(9*np.cos(l)*np.cos(b)+12*np.sin(l)*np.cos(b)+7*np.sin(b))
        v_final = vel_init+delv
        print("Input: vlsr=%.2f km/s, l_deg=%.4f, b_deg=%.4f"%(vel_init, l_deg, b_deg))
        print("Output: vhelio=%.2f km/s"%(v_final))
    else:
        delv = +9*np.cos(l)*np.cos(b)+12*np.sin(l)*np.cos(b)+7*np.sin(b)
        v_final = vel_init+delv
        print("Input: vhelio=%.2f km/s, l_deg=%.4f, b_deg=%.4f"%(vel_init, l_deg, b_deg))
        print("Output: vlsr=%.2f km/s"%(v_final))

    # print 'Velocity correction at this (RA, DEC) is (km/s): ', delv
    return v_final

##########################
# Change from Vhelio to Vlsr, or the inverse.
def vhelio2vlsr_Rosolowsky(vel_init, obj_ra, obj_dec,
    			   LSRD=False, LSRK=True, reverse=False):

    '''
    - Modify from the IDL program helio2lsr.pro by Erick Rosolowsky
       https://people.ok.ubc.ca/erosolo/idl/lib/helio2lsr.pro
    - J2000 coordinates is required! ra, dec should be in radian.
    - vel_init: velocity that need to be transformed
                vel_init = vhelio if reverse = False
                vel_init = vlsr   if reverse = True
    - obj_ra, obj_dec should be in radian
    - LSRD: the dynamical LSR frame
        - consider rotation of the Sun in the galactic disk.
    - LSRK: the kinematic LSR frame, being used in almost all cases
        - the Sun has a random velocity relative to the average of stars
          in its neighbourhood of about 20 km/s in direction (ra, dec)
          = (18h03m50.24s, +30d00m16.8") (J2000)
        - see http://coursewiki.astro.cornell.edu/Astro4410/RadialVelocities
          for detailed explaination.

     J2000 coordinates is required

     value adopted from http://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html
     Table 5, as of May 2010.
    '''
    import numpy as np
    import sys
    if LSRK:
        solarmotion_ra   = (18+3/60.+50.24/3.6e3)*15   # in degree
        solarmotion_dec  = 30.+0/60.+16.8/3.6e3
        solarmotion_vmag = 20.0


    elif LSRD:  # dynamical LSR==True by default
        solarmotion_ra = (17+49/60.+58.66/3.6e3)*15
        solarmotion_dec = 28+7/60.+3.92/3.6e3
        solarmotion_vmag = 16.55

    else:
        print('Assign either LSRK or LSRD to be True')
        sys.exit(1)


    from astropy.coordinates import SkyCoord
    solar_cd = SkyCoord(solarmotion_ra, solarmotion_dec, unit='deg')
    obj_cd = SkyCoord(obj_ra, obj_dec, unit='deg')

    gc_dist = obj_cd.separation(solar_cd)
    gc_dist = gc_dist.radian

    # vlsr --> vhelio
    if reverse:
        delv = -solarmotion_vmag*np.cos(gc_dist)
    # vhelio --> vlsr, by default
    else:
        delv = solarmotion_vmag*np.cos(gc_dist)
    print('Velocity correction at this (RA, DEC) is (km/s)', delv)
    return vel_init+delv

if __name__ == '__main__':
    import sys
    import numpy as np

    vel_init = np.float(sys.argv[1])
    l_deg = np.float(sys.argv[2])
    b_deg = np.float(sys.argv[3])

    vhelio2vlsr_Westmeier(vel_init, l_deg, b_deg)
