def read_johnson17_tb1(ion):
    from yztools.mstar2mhalo import mstar2mhalo
    from yztools.calc_r200 import calc_r200
    from astropy.table import Table
    import numpy as np

    if ion not in ['SiII', 'SiIII', 'SiIV', 'CIV', 'OVI']:
        print("Do not have this ion: {} in Johnson+2017, return []".format(ion))
        return []

    filename = '/Users/Yong/Dropbox/databucket/Johnson17_tb1.txt'
    johnson17 = Table.read(filename, format='ascii')

    ind = johnson17['flg_'+ion] != 'nan'
    johnson17_yz_ion = johnson17[ind].copy()
    # delete other ions
    for col_name in ['Wr_HI', 'eWr_HI', 'Wr_SiIII', 'flg_SiIII', 'eWr_SiIII',
                     'Wr_SiII', 'eWr_SiII', 'flg_SiII',
                     'Wr_SiIV', 'eWr_SiIV', 'flg_SiIV',
                     'Wr_CIV', 'eWr_CIV', 'flg_CIV',
                     'Wr_OVI', 'eWr_OVI', 'flg_OVI']:
        col_name_ion = col_name.split('_')[1]
        if ion != col_name_ion:
            johnson17_yz_ion.remove_column(col_name)

    johnson17_yz_ion['Wr_'+ion] = johnson17_yz_ion['Wr_'+ion]*1000 # from A to mA
    johnson17_yz_ion['eWr_'+ion] = johnson17_yz_ion['eWr_'+ion]*1000

    # for non detection, change from 2sigma to 3sigma
    j17_uplims = np.isnan(johnson17_yz_ion['eWr_'+ion])
    j17_det = np.logical_not(j17_uplims)
    johnson17_yz_ion['Wr_'+ion][j17_uplims] = johnson17_yz_ion['Wr_'+ion][j17_uplims]/2.*3
    johnson17_yz_ion['Wr_flg_'+ion] = [' '*8]*len(johnson17_yz_ion)
    johnson17_yz_ion['Wr_flg_'+ion][j17_uplims] = '<=(3sig)'
    johnson17_yz_ion['Wr_flg_'+ion][j17_det] = '='

    # process the rvir using my own code
    j17_mstar_Chabrier = 10**johnson17_yz_ion['logM*']

    # from Chabrier to Kroupa (because SMHM used Kroupa IMF),
    # conversion factor from Madau & Dickinson (2011)
    j17_mstar_yz = j17_mstar_Chabrier*0.66/0.61

    j17_mhalo_yz = mstar2mhalo(j17_mstar_yz)
    j17_r200c_yz, j17_r200m_yz = calc_r200(j17_mhalo_yz)
    j17_r200m_yz = j17_r200m_yz.value

    johnson17_yz_ion['logM*_yz'] = np.around(np.log10(j17_mstar_yz), decimals=2)
    johnson17_yz_ion['R200m_yz'] = np.around(j17_r200m_yz, decimals=1)

    johnson17_yz_ion.remove_columns(['logM*', 'Rh', 'flg_'+ion])
    johnson17_yz_ion.rename_column('d', 'impact_para_kpc')
    johnson17_yz_ion.rename_column('M_r', 'Mr(AB)')
    johnson17_yz_ion.rename_column('R.A.', 'RA(J2000)')
    johnson17_yz_ion.rename_column('Decl.', 'DEC(J2000)')
    johnson17_yz_ion.rename_column('ID', 'Galaxy')

    print('*'*30)
    print("Read table1 from Johnson+2017")
    print("including SiII1260, SiIII1206, SiIV1393, CIV1548, OVI1031, HI1215")
    print("these changes are made: only take {} values".format(ion))
    print("- updated M* from Chabrier to Kroupa")
    print("- Mhalo from M* based on Moster+2010")
    print("- Rvir based on r200, 200 times matter density")
    print("- Wr non-detections from 1sigma to 3 sigma")
    print("- No logN in tb1")

    return johnson17_yz_ion

def read_liang14_tb1_tb4(ion):
    """
    See their Table 4. They only have Wr values, no logN. For upper limits, they quote 2 sigma upper limits.
    They have SiII, SiIII, SiIV, CII, and CIV values. Their velocity integration windows are determined by eye.
    SiII's value is from 1260 instead of 1193; decide not to use SiII.
    Liang14's rvir value is calculated the same way as Johnson17 (mstar2mhalo from krastov,
    and rvir from Bryan and Norman). Even though they said the values are based on matter density,
    I believe the calcuation is meant for critical density.

    03/31/2021. YZ.

    """
    from yztools.mstar2mhalo import mstar2mhalo
    from yztools.calc_r200 import calc_r200
    from astropy.table import Table
    import numpy as np

    if ion not in ['SiII', 'SiIII', 'SiIV', 'CII', 'CIV']:
        print("Do not have this ion: {} in Liang+2014, return []".format(ion))
        return []

    # first read galaxy property
    liang14_tb1 = Table.read('/Users/Yong/Dropbox/databucket/Liang14_table1.dat', format='ascii')

    # ion measurements
    liang14_tb4 = Table.read('/Users/Yong/Dropbox/databucket/Liang14_table4.dat', format='ascii')

    l14_mstar_Chabrier = np.zeros(len(liang14_tb4))
    # l14_mhalo = np.zeros(len(liang14_tb4))
    # l14_rvir = np.zeros(len(liang14_tb4))

    for i in range(len(liang14_tb4)):
        igal = liang14_tb4['Galaxy'][i]
        ind = np.where(igal == liang14_tb1['Galaxy'])[0][0]
        l14_mstar_Chabrier[i] = 10**liang14_tb1['lgM*'][ind]
        #l14_mhalo[i] = 10**liang14_tb1['lgMh'][ind]
        #l14_rvir[i] = liang14_tb1['Rh'][ind]

    # from Chabrier to Kroupa, see Zheng+2020, IC1613
    l14_mstar_yz = l14_mstar_Chabrier * 0.66/0.61

    # from mstar to mhalo
    l14_mhalo_yz = mstar2mhalo(l14_mstar_yz)
    l14_r200c_yz, l14_r200m_yz = calc_r200(l14_mhalo_yz)
    l14_r200m_yz = l14_r200m_yz.value

    l14_impact_para = liang14_tb4['d_kpc']
    l14_impact_over_r200m = l14_impact_para/l14_r200m_yz

    ### write a new table
    liang14_yz = liang14_tb4.copy()
    liang14_yz['logM*_yz'] = np.around(np.log10(l14_mstar_yz), decimals=2)
    liang14_yz['R200m_yz'] = np.around(l14_r200m_yz, decimals=1)

    # only take targets with values
    has_wr = liang14_yz['Wr_'+ion] != '...'
    liang14_yz_ion = liang14_yz[has_wr].copy()

    # delete other ions
    for col_name in ['Wr_HI', 'eWr_HI', 'Wr_SiIII', 'eWr_SiIII', 'Wr_SiII', 'eWr_SiII',
                     'Wr_SiIV', 'eWr_SiIV', 'Wr_CII', 'eWr_CII', 'Wr_CIV', 'eWr_CIV']:
        col_name_ion = col_name.split('_')[1]
        if ion != col_name_ion:
            liang14_yz_ion.remove_column(col_name)

    # for non detection upper limits, it was quoted as 2 sigma, let's do 3 sigma instead
    ion_uplims = liang14_yz_ion['Wr_'+ion] == '<='
    ion_det = np.logical_not(ion_uplims)
    liang14_yz_ion['Wr_flg_'+ion] = [' '*8]*len(liang14_yz_ion)
    liang14_yz_ion['Wr_flg_'+ion][ion_uplims] = '<=(3sig)'
    liang14_yz_ion['Wr_flg_'+ion][ion_det] = '='

    liang14_yz_ion['Wr_'+ion][ion_uplims] = 'nan'
    liang14_yz_ion['Wr_'+ion] = np.asarray(liang14_yz_ion['Wr_'+ion], dtype='float')
    liang14_yz_ion['eWr_'+ion] = np.asarray(liang14_yz_ion['eWr_'+ion], dtype='float')

    liang14_yz_ion['Wr_'+ion][ion_uplims] = liang14_yz_ion['eWr_'+ion][ion_uplims]/2.*3
    liang14_yz_ion['eWr_'+ion][ion_uplims] = np.nan

    liang14_yz_ion.rename_column('z_gal', 'zgal')

    print('*'*30)
    print("Read from Liang+2014 tb1 & tb4")
    print("including SiII1260, SiIII1206, SiIV1393, CII1334, CIV1548, ")
    print("these changes are made: only take {} values".format(ion))
    print("- updated M* from Chabrier to Kroupa")
    print("- Mhalo from M* based on Moster+2010")
    print("- Rvir based on r200, 200 times matter density")
    print("- Wr non-detections from 1sigma to 3 sigma")
    print("- No logN")

    return liang14_yz_ion

def read_bordoloi14_tb1(ion):
    filename = '/Users/Yong/Dropbox/databucket/Bordoloi14_tb1.txt'
    from yztools.mstar2mhalo import mstar2mhalo
    from yztools.calc_r200 import calc_r200
    from astropy.table import Table
    import numpy as np

    if ion != 'CIV':
        print("Do not have this ion: {} in Bordoloi+2014, return []".format(ion))
        return []

    bordoloi14 = Table.read(filename, format='ascii')

    # update Wr' for upper limit,
    # multiple by 3 to get to 3 sig, to be
    # consistent with my measurements
    wr_nondet = np.isnan(bordoloi14['eWr'])
    wr_det = np.logical_not(wr_nondet)
    b14_Wr_CIV_yz = bordoloi14['Wr'].copy()  # in unit of mA
    b14_Wr_CIV_yz[wr_nondet] = 3*b14_Wr_CIV_yz[wr_nondet]

    # update logN, for upper limit, it was 2 sigma from the paper
    # now use 3 sigma for those non detections
    logN_nondet = bordoloi14['N_flg'] == '<'
    logN_det = np.logical_not(logN_nondet)
    b14_logN_yz = bordoloi14['logN'].copy()
    b14_logN_yz[logN_nondet] = np.around(np.log10((10**b14_logN_yz[logN_nondet])/2*3), decimals=2)

    # update stellar mass
    # from Salpeter to Kroupa (because SMHM used Kroupa IMF),
    # conversion factor of M_Kro/M_Sal = 0.66 from Madau & Dickinson (2011)
    b14_mstar_Salpeter = 10**bordoloi14['logM*']
    b14_mstar_yz = b14_mstar_Salpeter*0.66

    # update rvir
    b14_mhalo_yz = mstar2mhalo(b14_mstar_yz)
    b14_r200c, b14_r200m = calc_r200(b14_mhalo_yz)
    b14_r200c = b14_r200c.value
    b14_r200m_yz = b14_r200m.value

    b14_impact_para = bordoloi14['R'].copy()
    b14_impact_over_r200m = b14_impact_para/b14_r200m

    # write a new table with updated data
    bordoloi14_yz = bordoloi14.copy()
    bordoloi14_yz['logM*_yz'] = np.around(np.log10(b14_mstar_yz), decimals=2)
    bordoloi14_yz['R200m_yz'] = np.around(b14_r200m_yz, decimals=2)
    bordoloi14_yz['Wr_'+ion] = b14_Wr_CIV_yz # now non-detection is 3 sigma value
    bordoloi14_yz['Wr_flg_'+ion] = [' '*8]*len(bordoloi14_yz)
    bordoloi14_yz['Wr_flg_'+ion][wr_nondet] = '<=(3sig)'
    bordoloi14_yz['Wr_flg_'+ion][wr_det] = bordoloi14_yz['Wr_flg'][wr_det]

    bordoloi14_yz['logN_'+ion] = b14_logN_yz # now non detection is 3 sigma value
    bordoloi14_yz['N_flg_'+ion] = [' '*8]*len(bordoloi14_yz)
    bordoloi14_yz['N_flg_'+ion][logN_nondet] = '<=(3sig)'
    bordoloi14_yz['N_flg_'+ion][logN_det] = bordoloi14_yz['N_flg'][logN_det]
    bordoloi14_yz.remove_columns(['logM*', 'Rvir', 'Wr', '3sig_Wr', 'logN', 'N_flg', 'Wr_flg'])

    bordoloi14_yz.rename_column('RA', 'RA(J2000)')
    bordoloi14_yz.rename_column('DEC', 'DEC(J2000)')
    bordoloi14_yz.rename_column('z_sys', 'zgal')
    bordoloi14_yz.rename_column('R', 'impact_para_kpc')
    bordoloi14_yz.rename_column('elogN', 'elogN_'+ion)
    bordoloi14_yz.rename_column('eWr', 'eWr_'+ion)

    print('*'*30)
    print("Read table1 from Bordoloi+2014, only CIV, these changes are made: ")
    print("- updated M* from Salpeter to Kroupa")
    print("- Mhalo from M* based on Moster+2010")
    print("- Rvir based on r200, 200 times matter density")
    print("- Wr non-detections from 1sigma to 3 sigma")
    print("- logN non-detections from 2sigma to 3 sigma")
    return bordoloi14_yz
