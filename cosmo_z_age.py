import sys
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib as mpl
mpl.rcParams['font.family'] = 'stixgeneral'

from astropy.cosmology import Planck15 as cosmo
z = np.mgrid[0:6.01:0.01]
age = cosmo.age(z)

plt.figure(figsize=(5, 5))
plt.plot(z, age, color='k')
plt.grid('on', linestyle='--')
plt.xlabel('Redshift', fontsize=16)
plt.ylabel('Time since Big Bang (Gyrs)', fontsize=16)
plt.tight_layout()
plt.text(0.45, 13.2, "Planck15 Cosmology", fontsize=16)
figname = "/Users/Yong/Dropbox/databucket/redshift_age.pdf"
plt.savefig(figname)
plt.close()
print("Check ", figname)

if len(sys.argv) > 1:
    iz = np.float(sys.argv[1])
    iage = cosmo.age(iz)
    print("Age since BB is %.2f Gyrs at z=%.3f"%(iage.value, iz))
