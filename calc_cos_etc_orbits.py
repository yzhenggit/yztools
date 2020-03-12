import numpy as np
# only for cycle 25, 

etc_exptime = 5000 # seconds 
### etc_exptime is the exposure time estimated using ETC (not including buffer time)

first_exp_overhead = 14*60 #s# 6+(7or3)+5 mins for GS Acq + ACQ/Peaked + SciOverHead
subse_exp_overhead = 360  #s# 4+2 mins for GS reAcq + SciOverHead
change_config = 60 # Change of equipment (except to increment FT-POS; I don't get it)
cos_visorbit_time = 3240  # 54 min for any objects with -30<Dec<30. See table 6.1 in primier 

first_exp_time = cos_visorbit_time - first_exp_overhead
subse_exp_time = cos_visorbit_time - (subse_exp_overhead+change_config) ## assuming we change FP-POS??
orbits = np.ceil(1+(etc_exptime-first_exp_time)/subse_exp_time)

tot_obs_time = first_exp_time+(orbits-1)*subse_exp_time
tot_overhead = first_exp_overhead+(orbits-1)*subse_exp_overhead
print('%d orbits with %ds total exp. time and %ds total overhead time'%(orbits, tot_obs_time, tot_overhead))


