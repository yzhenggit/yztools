from astropy.coordinates import SkyCoord
import astropy.units as u
from yztools.mstar2mhalo import mstar2mhalo
from yztools.calc_r200 import calc_r200
from astropy.table import Table, Column
import numpy as np
from glob import glob

def index_radius_limit(gal_qso_tb, rho_limit=100):
    """
    only choose sightlines within rvir
    """

    # only keep those within impact parameters
    #ind_a = gal_qso_tb['R200m_kpc'] >= rho_limit
    # ind_b = gal_qso_tb['impact_para_kpc'] <= rho_limit
    ind_c = gal_qso_tb['impact_para_kpc'] <= gal_qso_tb['R200m_kpc']

    ind_keep = ind_c

    # rvir or 100, which ever is smaller
    # ind_a = gal_qso_tb['R200m_kpc'] <= rho_limit
    # ind_b = gal_qso_tb['impact_para_kpc'] <= rho_limit
    # ind_c = gal_qso_tb['impact_para_kpc'] <= gal_qso_tb['R200m_kpc']

    # ind_within_rvir = ind_a & ind_c
    # ind_within_rho = (~ind_a) & ind_b
    # ind_keep = ind_within_rho | ind_within_rvir

    # print("index_radius_limit: {} pairs with impact parameters less than R200m or {} kpc, whichever is smaller".format(len(gal_qso_tb[ind_keep]), rho_limit))
    print("index_radius_limit: {} pairs with impact parameters less than R200m".format(len(gal_qso_tb[ind_keep])))
    return gal_qso_tb[ind_keep].copy()

def compile_literature(snr_limit=8, logmstar_limit=9,
                       read_all_ic1613=False, read_all_qu22=False,
                        print_info=True, skip_zheng_lgs3_det=True,
                        skip_zheng_ddo210_det=True):
    ############# process Bordoloi+2014 table #############
    from yztools.read_datasets.read_dwarfcgm_ions import read_bordoloi14_tb1
    tb_b14 = read_bordoloi14_tb1(snr_limit=snr_limit, logmstar_limit=logmstar_limit, print_info=print_info)

    tb_b14['References'] = 'Bordoloi+2014'
    # add other ions for later combining with other tables
    for ikey in ['Wr_flg_HI', 'Wr_HI', 'eWr_HI', 'logN_flg_HI', 'logN_HI', 'elogN_HI',
                 'Wr_flg_SiII','Wr_SiII', 'eWr_SiII', 'logN_flg_SiII', 'logN_SiII', 'elogN_SiII',
                 'Wr_flg_SiIII', 'Wr_SiIII', 'eWr_SiIII', 'logN_flg_SiIII','logN_SiIII', 'elogN_SiIII',
                 'Wr_flg_SiIV', 'Wr_SiIV', 'eWr_SiIV','logN_flg_SiIV', 'logN_SiIV', 'elogN_SiIV',
                 'Wr_flg_CII', 'Wr_CII','eWr_CII', 'logN_flg_CII', 'logN_CII', 'elogN_CII', 'gal_vlsr_km/s',
                 'gal_dmpc', 'logMHI', 'MHI_ref',
                 'Wr_flg_OVI','Wr_OVI', 'eWr_OVI', 'logN_flg_OVI', 'logN_OVI', 'elogN_OVI', 'qso_snr_stis']:
        if 'flg' in ikey or ikey == 'MHI_ref':
            tb_b14[ikey] = ' '*8
        else:
            tb_b14[ikey] = np.nan

    ############# process Liang+2014 table #############
    from yztools.read_datasets.read_dwarfcgm_ions import read_liang14_tb1_tb2_tb4
    tb_l14 = read_liang14_tb1_tb2_tb4(snr_limit=snr_limit, logmstar_limit=logmstar_limit, print_info=print_info)
    tb_l14['References'] = 'Liang&Chen2014'

    # add other ions for later combining with other tables
    for ikey in [ 'Wr_flg_OVI','Wr_OVI', 'eWr_OVI', 'logN_flg_OVI', 'logN_OVI', 'elogN_OVI',
                 'logSFR_flg', 'gal_vlsr_km/s', 'gal_dmpc', 'logMHI', 'MHI_ref']:
        if ikey in ['Wr_flg_OVI', 'logN_flg_OVI', 'MHI_ref']:
            tb_l14[ikey] = ' '*8
        elif ikey == 'logSFR_flg':
            tb_l14[ikey] = '='
        else:
            tb_l14[ikey] = np.nan

    ############# process Johnson+2017 table #############
    from yztools.read_datasets.read_dwarfcgm_ions import read_johnson17_tb1
    tb_j17 = read_johnson17_tb1(snr_limit=snr_limit, logmstar_limit=logmstar_limit, print_info=print_info)
    tb_j17['References'] = 'Johnson+2017'

    # add other ions for later combining with other tables
    for ikey in ['Wr_flg_CII', 'Wr_CII','eWr_CII', 'logN_flg_CII', 'logN_CII', 'elogN_CII',
                 'logSFR_flg', 'logSFR(Msun/yr)', 'SFR_ref', 'logMHI', 'MHI_ref',
                 'qso_snr_stis', 'gal_dmpc', 'gal_vlsr_km/s']:
        if 'flg' in ikey or 'ref' in ikey:
            tb_j17[ikey] = ' '*8
        else:
            tb_j17[ikey] = np.nan

    ########## clean up columns ##########
    save_keys = ['gal_name', 'gal_ra_deg', 'gal_dec_deg', 'gal_z',
                 'logM*', 'M*_ref', 'R200m_kpc' , 'logMHI', 'MHI_ref',
                 'qso_hsla_name', 'qso_ra_deg', 'qso_dec_deg', 'qso_hsla_snr','impact_para_kpc',
                 'Wr_flg_HI', 'Wr_HI', 'eWr_HI', 'logN_flg_HI', 'logN_HI', 'elogN_HI',
                 'Wr_flg_SiII', 'Wr_SiII', 'eWr_SiII', 'logN_flg_SiII', 'logN_SiII', 'elogN_SiII',
                 'Wr_flg_SiIII', 'Wr_SiIII', 'eWr_SiIII', 'logN_flg_SiIII','logN_SiIII', 'elogN_SiIII',
                 'Wr_flg_SiIV', 'Wr_SiIV', 'eWr_SiIV','logN_flg_SiIV', 'logN_SiIV', 'elogN_SiIV',
                 'Wr_flg_CII', 'Wr_CII','eWr_CII', 'logN_flg_CII', 'logN_CII', 'elogN_CII',
                 'Wr_flg_CIV', 'Wr_CIV','eWr_CIV', 'logN_flg_CIV', 'logN_CIV', 'elogN_CIV',
                 'Wr_flg_OVI','Wr_OVI', 'eWr_OVI', 'logN_flg_OVI', 'logN_OVI', 'elogN_OVI',
                'References', 'logSFR_flg', 'logSFR(Msun/yr)', 'SFR_ref',
                'qso_snr_stis', 'gal_dmpc', 'gal_vlsr_km/s']

    for key in tb_b14.colnames:
        if key not in save_keys:
            #print("Removing {} from Bordoloi+2014".format(key))
            tb_b14.remove_column(key)

    for key in tb_l14.colnames:
        if key not in save_keys:
            #print("Removing {} from Liang+2014".format(key))
            tb_l14.remove_column(key)

    for key in tb_j17.colnames:
        if key not in save_keys:
            #print("Removing {} from Johnson+2017".format(key))
            tb_j17.remove_column(key)

    ###### read in Zheng+2019a, WLM ####
    zheng19_wlm = read_zheng19_wlm(print_info=print_info)

    ##### read in Zheng+2020b, IC1613 ###
    if read_all_ic1613 == True:
        zheng20_ic1613 = read_zheng20_ic1613(print_info=print_info)
    else:
        zheng20_ic1613 = read_zheng20_ic1613(snr_limit=snr_limit, print_info=print_info)

    #### read in Qu & Bregman 2022 ####
    if read_all_qu22 == True:
        qu22_tb = read_qu22_NGC3109_SexA_SexB(print_info=print_info)
    else:
        qu22_tb = read_qu22_NGC3109_SexA_SexB(snr_limit=snr_limit, print_info=print_info)

    #### read archival work ###
    zheng_newdw_tb = read_zheng_newdw_tb(snr_limit=snr_limit,
                                              logmstar_limit=logmstar_limit,
                                              print_info=print_info,
                                              skip_zheng_lgs3_det=skip_zheng_lgs3_det,
                                              skip_zheng_ddo210_det=skip_zheng_ddo210_det)

    ###### combine all values together ######
    from astropy.table import vstack
    tb_main = vstack([zheng_newdw_tb, zheng19_wlm, zheng20_ic1613, qu22_tb, tb_l14, tb_b14,  tb_j17])
    tb_main['gal_ra_deg'] = np.around(tb_main['gal_ra_deg'], decimals=4)
    tb_main['gal_dec_deg'] = np.around(tb_main['gal_dec_deg'], decimals=4)
    tb_main['qso_ra_deg'] = np.around(tb_main['qso_ra_deg'], decimals=4)
    tb_main['qso_dec_deg'] = np.around(tb_main['qso_dec_deg'], decimals=4)
    tb_main['qso_hsla_snr'] = np.around(tb_main['qso_hsla_snr'], decimals=1)

    ind_sort = np.argsort(tb_main['logM*'])
    tb_main = tb_main[ind_sort].copy()

    ### add new order ####
    ########## clean up columns ##########
    new_order = ['References', 'gal_name', 'gal_ra_deg', 'gal_dec_deg', 'gal_vlsr_km/s', 'gal_dmpc', 'gal_z',
                 'logM*', 'M*_ref', 'logMHI', 'MHI_ref', 'R200m_kpc', 'logSFR_flg', 'logSFR(Msun/yr)', 'SFR_ref',
                 'qso_hsla_name', 'qso_ra_deg', 'qso_dec_deg', 'qso_hsla_snr', 'qso_snr_stis', 'impact_para_kpc',
                 'Wr_flg_HI', 'Wr_HI', 'eWr_HI', 'logN_flg_HI', 'logN_HI', 'elogN_HI',
                 'Wr_flg_SiII', 'Wr_SiII', 'eWr_SiII', 'logN_flg_SiII', 'logN_SiII', 'elogN_SiII',
                 'Wr_flg_SiIII', 'Wr_SiIII', 'eWr_SiIII', 'logN_flg_SiIII','logN_SiIII', 'elogN_SiIII',
                 'Wr_flg_SiIV', 'Wr_SiIV', 'eWr_SiIV','logN_flg_SiIV', 'logN_SiIV', 'elogN_SiIV',
                 'Wr_flg_CII', 'Wr_CII','eWr_CII', 'logN_flg_CII', 'logN_CII', 'elogN_CII',
                 'Wr_flg_CIV', 'Wr_CIV','eWr_CIV', 'logN_flg_CIV', 'logN_CIV', 'elogN_CIV',
                 'Wr_flg_OVI','Wr_OVI', 'eWr_OVI', 'logN_flg_OVI', 'logN_OVI', 'elogN_OVI']
    tb_main_neworder = tb_main[new_order].copy()

    print("compile_literature: {} pairs meeting input criteria".format(len(tb_main_neworder)))
    return tb_main_neworder

def read_johnson17_tb1(snr_limit=-999, logmstar_limit=999, print_info=True):
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from yztools.mstar2mhalo import mstar2mhalo
    from yztools.calc_r200 import calc_r200
    from astropy.table import Table
    import numpy as np


    filename = '/Users/Yong/Dropbox/databucket/Johnson17_tb1.txt'
    johnson17 = Table.read(filename, format='ascii')
    gal_coord = SkyCoord(ra=johnson17['R.A.'], dec=johnson17['Decl.'],
                     frame='icrs', unit=(u.hour, u.deg))

    # first find qso sightlines
    # cross match johnson17_tb1 with HSLA and see what's the nearest qso
    hsla_tb = Table.read('/Users/Yong/Dropbox/databucket/hsla_all_targets_COS_sample.fits',
                         format='fits')
    hsla_coord = SkyCoord(ra=hsla_tb['RA'], dec=hsla_tb['DEC'], unit=(u.deg, u.deg), frame='icrs')
    qso_name = []
    qso_ra = []
    qso_dec = []
    qso_snr = []
    for i in range(len(johnson17)):
        igal_coord = gal_coord[i]
        sep_arcsec = igal_coord.separation(hsla_coord).arcsec
        ind_min = np.argmin(sep_arcsec)
        # print(johnson17_tb1['qso_name'][i], hsla_tb['Target Name'][ind_min])
        # checked with Johnson+2015, this is right
        qso_name.append(hsla_tb['Target Name'][ind_min])
        qso_ra.append(np.around(hsla_tb['RA'][ind_min], decimals=4))
        qso_dec.append(np.around(hsla_tb['DEC'][ind_min], decimals=4))
        qso_snr.append(np.around(hsla_tb['median S/N'][ind_min], decimals=1))

    johnson17['qso_hsla_name'] = np.asarray(qso_name)
    johnson17['qso_ra_deg'] = np.asarray(qso_ra)
    johnson17['qso_dec_deg'] = np.asarray(qso_dec)
    johnson17['qso_hsla_snr'] = np.asarray(qso_snr)

    for ion in ['HI', 'SiIII', 'SiII', 'SiIV', 'CIV', 'OVI']:
        #print(ion)
        johnson17['Wr_'+ion] = (johnson17['Wr_'+ion])*1000 # from A to mA
        johnson17['eWr_'+ion] = (johnson17['eWr_'+ion])*1000

        # for non detection, change from 2sigma to 3sigma
        j17_uplims = np.isnan(johnson17['eWr_'+ion])
        j17_hasvalues = np.isfinite(johnson17['Wr_'+ion])
        j17_det = np.logical_not(j17_uplims)
        johnson17['Wr_'+ion][j17_uplims] = johnson17['Wr_'+ion][j17_uplims]/2.*3
        johnson17['Wr_flg_'+ion] = [' '*8]*len(johnson17)
        johnson17['Wr_flg_'+ion][j17_uplims & j17_hasvalues] = '<=(3sig)'
        johnson17['Wr_flg_'+ion][j17_det] = '='
        if 'flg_'+ion in johnson17.colnames:
            johnson17.remove_columns('flg_'+ion)

    # process the rvir using my own code
    j17_mstar_Chabrier = 10**johnson17['logM*'].copy()

    # from Chabrier to Kroupa (because SMHM used Kroupa IMF),
    # conversion factor from Madau & Dickinson (2011)
    j17_mstar_kroupa = j17_mstar_Chabrier*0.66/0.61
    log_j17_mstar_kroupa = np.log10(j17_mstar_kroupa)

    log_j17_mhalo_yz = mstar2mhalo(log_j17_mstar_kroupa, ref='B13-GK14')
    j17_mhalo_yz = 10**log_j17_mhalo_yz
    j17_r200c_yz, j17_r200m_yz = calc_r200(j17_mhalo_yz)
    j17_r200m_yz = j17_r200m_yz.value

    johnson17['logM*'] = np.around(log_j17_mstar_kroupa, decimals=2)
    johnson17['M*_ref'] = 'Johnson+2017(Kroupa IMF)'
    johnson17['R200m_kpc'] = np.around(j17_r200m_yz, decimals=1)

    gal_coord = SkyCoord(ra=johnson17['R.A.'], dec=johnson17['Decl.'],
                         frame='icrs', unit=(u.hour, u.deg))
    johnson17['gal_ra_deg'] = np.around(gal_coord.icrs.ra.deg, decimals=4)
    johnson17['gal_dec_deg'] = np.around(gal_coord.icrs.dec.deg, decimals=4)

    johnson17.remove_column('Rh')
    johnson17.rename_column('d', 'impact_para_kpc')
    johnson17.rename_column('M_r', 'Mr(AB)')
    #johnson17.remove_column('R.A.'')
    #johnson17.remove_column('Decl.')
    johnson17.rename_column('ID', 'gal_name')
    johnson17.rename_column('zgal', 'gal_z')

    # convert HI Wr to logN assuming optically thin
    from yztools.llist import load_atomic_data
    for ion, ion_line in zip(['HI', 'SiII', 'SiIII', 'SiIV', 'CIV', 'OVI'],
                             ['HI1215', 'SiII1260', 'SiIII1206', 'SiIV1393', 'CIV1548', 'OVI1031']):
        line_info = load_atomic_data.morton03(ion)
        line_ind = np.where(np.asarray(line_info['line_name']) == ion_line)[0][0]
        # print(line_info['line_name'][line_ind])
        line_f = line_info['f'][line_ind]
        line_lambda = line_info['lambda'][line_ind]
        ion_N = 1.13e17*johnson17['Wr_'+ion]/line_f/line_lambda**2
        ion_eN = 1.13e17*johnson17['eWr_'+ion]/line_f/line_lambda**2
        ion_logN = np.log10(ion_N)
        ion_elogN = ion_eN/ion_N/np.log(10)
        johnson17['logN_flg_'+ion] = johnson17['Wr_flg_'+ion].copy()
        johnson17['logN_'+ion] = np.around(ion_logN, decimals=2)
        johnson17['elogN_'+ion] = np.around(ion_elogN, decimals=2)
    #for ion in ['HI', 'SiII', 'SiIII', 'SiIV', 'CIV', 'OVI']:
    #    johnson17['logN_'+ion] = np.zeros(len(johnson17))+np.nan
    #    johnson17['elogN_'+ion] = np.zeros(len(johnson17))+np.nan
    #    johnson17['logN_flg_'+ion] = johnson17['Wr_flg_'+ion]

    # for D1 and D2, add values from table 2
    # for errors, we take the mean of the two
    ind_d1 = np.where(johnson17['gal_name'] == 'D1')[0][0]
    johnson17[ind_d1]['logN_HI'] = 15.7
    johnson17[ind_d1]['elogN_HI'] = 0.4
    johnson17[ind_d1]['logN_SiIII'] = 13.14
    johnson17[ind_d1]['elogN_SiIII'] = 0.03
    johnson17[ind_d1]['logN_CIV'] = 13.73
    johnson17[ind_d1]['elogN_CIV'] = 0.04
    johnson17[ind_d1]['logN_OVI'] = 14.10
    johnson17[ind_d1]['elogN_OVI'] = 0.03

    ind_d2 = np.where(johnson17['gal_name'] == 'D2')[0][0]
    johnson17[ind_d2]['logN_HI'] = 15.06
    johnson17[ind_d2]['elogN_HI'] = 0.02
    johnson17[ind_d2]['logN_SiIII'] = 12.48
    johnson17[ind_d2]['elogN_SiIII'] = 0.05
    johnson17[ind_d2]['logN_OVI'] = 14.17
    johnson17[ind_d2]['elogN_OVI'] = 0.02
    johnson17.show_in_notebook()

    new_order = ['gal_name', 'gal_ra_deg', 'gal_dec_deg','gal_z', 'Mr(AB)',
                 'impact_para_kpc', 'logM*', 'M*_ref', 'R200m_kpc',
                 'qso_hsla_name', 'qso_ra_deg', 'qso_dec_deg', 'qso_hsla_snr',
                 'Wr_flg_HI', 'Wr_HI', 'eWr_HI',
                 'logN_flg_HI', 'logN_HI', 'elogN_HI',
                 'Wr_flg_SiII','Wr_SiII', 'eWr_SiII',
                 'logN_flg_SiII', 'logN_SiII', 'elogN_SiII',
                 'Wr_flg_SiIII', 'Wr_SiIII', 'eWr_SiIII',
                 'logN_flg_SiIII', 'logN_SiIII', 'elogN_SiIII',
                 'Wr_flg_SiIV', 'Wr_SiIV', 'eWr_SiIV',
                 'logN_flg_SiIV', 'logN_SiIV', 'elogN_SiIV',
                 'Wr_flg_CIV', 'Wr_CIV', 'eWr_CIV',
                 'logN_flg_CIV', 'logN_CIV', 'elogN_CIV',
                 'Wr_flg_OVI', 'Wr_OVI', 'eWr_OVI',
                 'logN_flg_OVI', 'logN_OVI', 'elogN_OVI']
    johnson17_neworder =johnson17[new_order].copy()

    # only return table with certain snr and mstar
    snr_cut = johnson17_neworder['qso_hsla_snr'] >= snr_limit
    mstar_cut = johnson17_neworder['logM*'] <= logmstar_limit
    johnson17_final = johnson17_neworder[snr_cut & mstar_cut].copy()

    if print_info == True:
        print('*'*30)
        print("Read table1 from Johnson+2017, qso SNR from HSLA")
        print("including SiII1260, SiIII1206, SiIV1393, CIV1548, OVI1031, HI1215")
        print("these changes are made:")
        print("- updated M* from Chabrier to Kroupa")
        print("- Mhalo from M* based on Moster+2010")
        print("- Rvir based on r200, 200 times matter density")
        print("- Wr non-detections from 1sigma to 3 sigma")
        print("- logN from tb2 when available, otherwise converted from Wr assuming optically thin.")
        print("- only those with SNR>={}, logMstar<={}".format(snr_limit, logmstar_limit))

    return johnson17_final

def read_liang14_tb1_tb2_tb4(snr_limit=-999, logmstar_limit=999, print_info=True):
    """
    See their Table 4. They only have Wr values, no logN. For upper limits, they quote 2 sigma upper limits.
    They have SiII, SiIII, SiIV, CII, and CIV values. Their velocity integration windows are determined by eye.
    SiII's value is from 1260 instead of 1193; decide not to use SiII.
    Liang14's rvir value is calculated the same way as Johnson17 (mstar2mhalo from krastov,
    and rvir from Bryan and Norman). Even though they said the values are based on matter density,
    I believe the calcuation is meant for critical density.

    03/31/2021. YZ.

    """
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from yztools.mstar2mhalo import mstar2mhalo
    from yztools.calc_r200 import calc_r200
    from astropy.table import Table
    import numpy as np

    # first read galaxy property
    liang14_tb1 = Table.read('/Users/Yong/Dropbox/databucket/Liang14_table1.dat', format='ascii')

    # read QSO data
    liang14_tb2 = Table.read('/Users/Yong/Dropbox/databucket/Liang14_table2.dat', format='ascii')

    # ion measurements; final table is based on this one
    liang14_tb4 = Table.read('/Users/Yong/Dropbox/databucket/Liang14_table4.dat', format='ascii')
    liang14_yz = liang14_tb4.copy()

    ########## first based on tb4, find same galaxy in tb1 ###########
    ind_in_tb1 = np.zeros(len(liang14_tb4), dtype=int)
    for i in range(len(liang14_tb4)):
        igal = liang14_tb4['Galaxy'][i]
        ind_tb1 = np.where(igal == liang14_tb1['Galaxy'])[0][0]
        ind_in_tb1[i] = ind_tb1

    liang14_yz['gal_ra_deg'] = liang14_tb1['RA_deg'][ind_in_tb1]
    liang14_yz['gal_dec_deg'] = liang14_tb1['Dec_deg'][ind_in_tb1]
    # liang14_yz['gal_z'] = liang14_tb1['z_spec'][ind_in_tb1]
    liang14_yz['logSFR(Msun/yr)'] = liang14_tb1['lgSFR'][ind_in_tb1]
    sfr_nan = liang14_yz['logSFR(Msun/yr)'] == -9.99
    liang14_yz['logSFR(Msun/yr)'][sfr_nan] = np.nan
    liang14_yz['lgM*'] = liang14_tb1['lgM*'][ind_in_tb1]
    liang14_yz['M*_flg'] = liang14_tb1['flag'][ind_in_tb1]
    liang14_yz.rename_column('z_gal', 'gal_z')
    liang14_yz.rename_column('Galaxy', 'gal_name')
    liang14_yz.rename_column('d_kpc', 'impact_para_kpc')
    liang14_yz.rename_column('theta_arcsec', 'impact_para_arcsec')

    ############## now add qso information from tb2 ########
    ind_in_tb2 = np.zeros(len(liang14_yz), dtype=int)
    all_qso_ra = liang14_tb2['RA(deg)']
    all_qso_dec = liang14_tb2['Dec(deg)']
    all_qso_coord = SkyCoord(ra=all_qso_ra, dec=all_qso_dec,
                         unit=(u.deg, u.deg), frame='icrs')

    for i in range(len(liang14_yz)):
        i_gal_theta = liang14_yz['impact_para_arcsec'][i]
        i_gal_ra = liang14_yz['gal_ra_deg'][i]
        i_gal_dec = liang14_yz['gal_dec_deg'][i]
        i_gal_coord = SkyCoord(ra=i_gal_ra, dec=i_gal_dec,
                                unit=(u.deg, u.deg), frame='icrs')


        sep_arcsec = i_gal_coord.separation(all_qso_coord).arcsec
        ind_min = np.argmin(sep_arcsec)
        if abs(i_gal_theta - sep_arcsec[ind_min]) > 1:
            print(liang14_yz['gal_name'][i])
        else:
            ind_in_tb2[i] = ind_min

    liang14_yz['qso_hsla_name'] = liang14_tb2['QSO'][ind_in_tb2]
    liang14_yz['qso_ra_deg'] = liang14_tb2['RA(deg)'][ind_in_tb2]
    liang14_yz['qso_dec_deg'] = liang14_tb2['Dec(deg)'][ind_in_tb2]
    liang14_yz['qso_z'] = liang14_tb2['z_QSO'][ind_in_tb2]
    liang14_yz['qso_snr_g130m'] = np.asarray(liang14_tb2['SN_G130M'][ind_in_tb2], dtype=float)
    liang14_yz['qso_snr_g160m'] = np.asarray(liang14_tb2['SN_G160M'][ind_in_tb2], dtype=float)
    liang14_yz['qso_snr_stis'] = np.asarray(liang14_tb2['SN_STIS'][ind_in_tb2], dtype=float)

    for key in ['qso_snr_g130m', 'qso_snr_g160m', 'qso_snr_stis']:
        ind_nan = liang14_yz[key] == -999
        liang14_yz[key][ind_nan] = np.nan

    # change mass from chabrier to Kroupa, see Zheng+2020, IC1613
    l14_mstar_Chabrier = 10**liang14_yz['lgM*']
    l14_mstar_kroupa = l14_mstar_Chabrier * 0.66/0.61
    log_l14_mstar_kroupa = np.log10(l14_mstar_kroupa)
    liang14_yz.remove_columns(['lgM*', 'z_Lya'])

    # from mstar to mhalo
    log_l14_mhalo_yz = mstar2mhalo(log_l14_mstar_kroupa, ref='B13-GK14')
    l14_mhalo_yz = 10**log_l14_mhalo_yz
    l14_r200c_yz, l14_r200m_yz = calc_r200(l14_mhalo_yz)
    l14_r200m_yz = l14_r200m_yz.value

    ### write a new table
    liang14_yz['logM*'] = np.around(log_l14_mstar_kroupa, decimals=2)
    liang14_yz['M*_ref'] = 'Liang+2014(Kroupa IMF)'
    liang14_yz['R200m_kpc'] = np.around(l14_r200m_yz, decimals=1)

    # for non detection upper limits, it was quoted as 2 sigma, let's do 3 sigma instead
    for ion in ['HI', 'SiII', 'SiIII', 'SiIV', 'CII', 'CIV']:
        # ion = 'SiII'
        ion_novals = liang14_yz['Wr_'+ion] == '...'
        ion_uplims = liang14_yz['Wr_'+ion] == '<='
        ion_det = np.logical_not(np.logical_or(ion_uplims, ion_novals))

        liang14_yz['Wr_flg_'+ion] = [' '*8]*len(liang14_yz)
        liang14_yz['Wr_flg_'+ion][ion_uplims] = '<=(3sig)'
        liang14_yz['Wr_flg_'+ion][ion_det] = '='

        liang14_yz['Wr_'+ion][ion_novals] = np.nan
        liang14_yz['Wr_'+ion][ion_uplims] = np.nan
        liang14_yz['Wr_'+ion] = np.asarray(liang14_yz['Wr_'+ion],  dtype='float')

        liang14_yz['eWr_'+ion][ion_novals] = np.nan
        liang14_yz['eWr_'+ion] = np.asarray(liang14_yz['eWr_'+ion],  dtype='float')

        liang14_yz['Wr_'+ion][ion_uplims] = liang14_yz['eWr_'+ion][ion_uplims]/2.*3
        liang14_yz['eWr_'+ion][ion_uplims] = np.nan

    # cross match Liang14 QSOs with HSLA and find the HSLA SNR
    l14_qso_coord = SkyCoord(ra=liang14_yz['qso_ra_deg'], dec=liang14_yz['qso_dec_deg'],
                             frame='icrs', unit=(u.deg, u.deg))

    # this is all COS targets (including QSO, gal, stars, etc.) from HSLA 2nd
    hsla_tb = Table.read('/Users/Yong/Dropbox/databucket//hsla_all_targets_COS_sample.fits',
                        format='fits')
    hsla_coord = SkyCoord(ra=hsla_tb['RA'], dec=hsla_tb['DEC'],
    frame='icrs', unit=(u.deg, u.deg))

    l14_qso_snr = np.zeros(len(liang14_yz))+np.nan

    for i in range(len(liang14_yz)):
        i_qso_coord = l14_qso_coord[i]
        i_qso_name = liang14_yz['qso_hsla_name'][i]

        sep_arcsec = i_qso_coord.separation(hsla_coord).arcsec
        ind_min = np.argmin(sep_arcsec)
        hsla_name = hsla_tb['Target Name'][ind_min]
        # checked
        if (i_qso_name == hsla_name) or (sep_arcsec[ind_min] < 5):
            l14_qso_snr[i] = hsla_tb['median S/N'][ind_min]
            liang14_yz['qso_hsla_snr'] = np.around(l14_qso_snr, decimals=1)

    # convert HI Wr to logN assuming optically thin
    from yztools.llist import load_atomic_data
    for ion, ion_line in zip(['HI', 'SiII', 'SiIII', 'SiIV', 'CII', 'CIV'],
                             ['HI1215', 'SiII1260', 'SiIII1206', 'SiIV1393', 'CII1334', 'CIV1548']):

        line_info = load_atomic_data.morton03(ion)
        line_ind = np.where(np.asarray(line_info['line_name']) == ion_line)[0][0]

        # print(line_info['line_name'][line_ind])
        line_f = line_info['f'][line_ind]
        line_lambda = line_info['lambda'][line_ind]
        ion_N = 1.13e17*liang14_yz['Wr_'+ion]/line_f/line_lambda**2
        ion_eN = 1.13e17*liang14_yz['eWr_'+ion]/line_f/line_lambda**2
        ion_logN = np.log10(ion_N)
        ion_elogN = ion_eN/ion_N/np.log(10)

        liang14_yz['logN_flg_'+ion] = liang14_yz['Wr_flg_'+ion].copy()
        liang14_yz['logN_'+ion] = np.around(ion_logN, decimals=2)
        liang14_yz['elogN_'+ion] = np.around(ion_elogN, decimals=2)

    liang14_yz['SFR_ref'] = '?'
    new_order = ['gal_name', 'gal_ra_deg', 'gal_dec_deg', 'gal_z',
                 'M*_flg', 'logM*', 'M*_ref', 'logSFR(Msun/yr)', 'SFR_ref', 'R200m_kpc',
                 'impact_para_arcsec', 'impact_para_kpc',
                 'qso_hsla_name', 'qso_ra_deg', 'qso_dec_deg', 'qso_z',
                 'qso_snr_g130m', 'qso_snr_g160m','qso_hsla_snr', 'qso_snr_stis',
                 'Wr_flg_HI', 'Wr_HI', 'eWr_HI',
                 'logN_flg_HI', 'logN_HI', 'elogN_HI',
                 'Wr_flg_SiII','Wr_SiII', 'eWr_SiII',
                 'logN_flg_SiII', 'logN_SiII', 'elogN_SiII',
                 'Wr_flg_SiIII', 'Wr_SiIII', 'eWr_SiIII',
                 'logN_flg_SiIII', 'logN_SiIII', 'elogN_SiIII',
                 'Wr_flg_SiIV', 'Wr_SiIV', 'eWr_SiIV',
                 'logN_flg_SiIV', 'logN_SiIV', 'elogN_SiIV',
                 'Wr_flg_CII', 'Wr_CII', 'eWr_CII',
                 'logN_flg_CII', 'logN_CII', 'elogN_CII',
                 'Wr_flg_CIV', 'Wr_CIV', 'eWr_CIV',
                 'logN_flg_CIV', 'logN_CIV', 'elogN_CIV']

    liang14_yz_neworder = liang14_yz[new_order].copy()

    # only those within certains snr and mstar limit
    #g130m_cut = liang14_yz_neworder['qso_snr_g130m'] >= snr_limit
    #g160m_cut = liang14_yz_neworder['qso_snr_g160m'] >= snr_limit
    cos_cut = liang14_yz_neworder['qso_hsla_snr'] >= snr_limit
    stis_cut = liang14_yz_neworder['qso_snr_stis'] >= snr_limit
    spec_cut = cos_cut | stis_cut
    mstar_cut = liang14_yz_neworder['logM*'] <= logmstar_limit
    liang14_yz_final = liang14_yz_neworder[spec_cut & mstar_cut].copy()

    if print_info == True:
        print('*'*30)
        print("Read from Liang+2014 tb1 & tb4")
        print("including HI1215, SiII1260, SiIII1206, SiIV1393, CII1334, CIV1548, ")
        print("these changes are made: only take {} values".format(ion))
        print("- updated M* from Chabrier to Kroupa")
        print("- Mhalo from M* based on Moster+2010")
        print("- Rvir based on r200, 200 times matter density")
        print("- Wr non-detections from 1sigma to 3 sigma")
        print("- No logN values given in paper, convert values from Wr assuming optically thin.")
        print("- only those with SNR>={}, logMstar<={}".format(snr_limit, logmstar_limit))

    return liang14_yz_final

def read_bordoloi14_tb1(snr_limit=-999, logmstar_limit=999, print_info=True):
    filename = '/Users/Yong/Dropbox/databucket/Bordoloi14_tb1.txt'
    from yztools.mstar2mhalo import mstar2mhalo
    from yztools.calc_r200 import calc_r200
    from astropy.table import Table
    import numpy as np


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
    logN_nondet = bordoloi14['logN_flg'] == '<'
    logN_det = np.logical_not(logN_nondet)
    b14_logN_yz = bordoloi14['logN'].copy()
    b14_logN_yz[logN_nondet] = np.around(np.log10((10**b14_logN_yz[logN_nondet])/2*3), decimals=2)

    # update stellar mass
    # from Salpeter to Kroupa (because SMHM used Kroupa IMF),
    # conversion factor of M_Kro/M_Sal = 0.66 from Madau & Dickinson (2011)
    b14_mstar_Salpeter = 10**bordoloi14['logM*'].copy()
    b14_mstar_kroupa = b14_mstar_Salpeter*0.66
    log_b14_mstar_kroupa = np.log10(b14_mstar_kroupa)

    # update star formation rate
    b14_sfr_flg = [' ']*len(bordoloi14)
    b14_log_ssfr = np.zeros(len(bordoloi14))+np.nan
    for i in range(len(bordoloi14)):
        i_log_ssfr = bordoloi14['logsSFR'][i]
        if '<' in i_log_ssfr:
            b14_sfr_flg[i] = '<'
            b14_log_ssfr[i] = float(i_log_ssfr.split('<')[1])
        else:
            b14_sfr_flg[i] = '='
            b14_log_ssfr[i] = float(i_log_ssfr)

    b14_sfr = (10**b14_log_ssfr) * b14_mstar_Salpeter
    b14_logsfr = np.around(np.log10(b14_sfr), decimals=2)

    # update rvir
    log_b14_mhalo_yz = mstar2mhalo(log_b14_mstar_kroupa, ref='B13-GK14')
    b14_mhalo_yz = 10**log_b14_mhalo_yz
    b14_r200c, b14_r200m = calc_r200(b14_mhalo_yz)
    b14_r200c = b14_r200c.value
    b14_r200m_yz = b14_r200m.value

    b14_impact_para = bordoloi14['R'].copy()
    b14_impact_over_r200m = b14_impact_para/b14_r200m

    # write a new table with updated data
    ion = 'CIV'
    bordoloi14_yz = bordoloi14.copy()
    bordoloi14_yz['logSFR_flg'] = b14_sfr_flg
    bordoloi14_yz['logSFR(Msun/yr)'] = b14_logsfr
    bordoloi14_yz['logM*'] = np.around(log_b14_mstar_kroupa, decimals=2)
    bordoloi14_yz['M*_ref'] = 'Bordoloi14(Kroupa IMF)'
    bordoloi14_yz['R200m_kpc'] = np.around(b14_r200m_yz, decimals=2)
    bordoloi14_yz['Wr_'+ion] = b14_Wr_CIV_yz # now non-detection is 3 sigma value
    bordoloi14_yz['Wr_flg_'+ion] = [' '*8]*len(bordoloi14_yz)
    bordoloi14_yz['Wr_flg_'+ion][wr_nondet] = '<=(3sig)'
    bordoloi14_yz['Wr_flg_'+ion][wr_det] = bordoloi14_yz['Wr_flg'][wr_det]

    bordoloi14_yz['logN_'+ion] = b14_logN_yz # now non detection is 3 sigma value
    bordoloi14_yz['logN_flg_'+ion] = [' '*8]*len(bordoloi14_yz)
    bordoloi14_yz['logN_flg_'+ion][logN_nondet] = '<=(3sig)'
    bordoloi14_yz['logN_flg_'+ion][logN_det] = bordoloi14_yz['logN_flg'][logN_det]
    bordoloi14_yz.remove_columns(['Rvir', 'Wr', '3sig_Wr', 'logsSFR',
                                  'logN', 'logN_flg', 'Wr_flg'])
    bordoloi14_yz.rename_column('z_sys', 'gal_z')
    bordoloi14_yz.rename_column('R', 'impact_para_kpc')
    bordoloi14_yz.rename_column('elogN', 'elogN_'+ion)
    bordoloi14_yz.rename_column('eWr', 'eWr_'+ion)
    bordoloi14_yz.rename_column('QSO-Name', 'qso_name')
    bordoloi14_yz.rename_column('Galaxy', 'gal_name')

    # change ra/dec from hour/deg to deg deg
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    gal_coord = SkyCoord(ra=bordoloi14_yz['RA'], dec=bordoloi14_yz['DEC'],
                         unit=(u.hour, u.deg),  frame='icrs')
    bordoloi14_yz['gal_ra_deg'] = np.around(gal_coord.icrs.ra.deg, decimals=4)
    bordoloi14_yz['gal_dec_deg'] = np.around(gal_coord.icrs.dec.deg, decimals=4)

    # add SNR for QSO from HSLA,tested
    # this is all COS targets (including QSO, gal, stars, etc.) from HSLA 2nd
    hsla_tb = Table.read('/Users/Yong/Dropbox/databucket//hsla_all_targets_COS_sample.fits',
                         format='fits')
    import re
    ind_in_hsla = np.zeros(len(bordoloi14_yz), dtype=int)

    for i in range(len(bordoloi14_yz)):
        qso_name = bordoloi14_yz['qso_name'][i]
        qso_ss = re.split('(\d+)', qso_name)
        qso_s1 = qso_ss[1]
        qso_s2 = qso_ss[3]
        for h in range(len(hsla_tb)):
            hsla_qso_name = hsla_tb['Target Name'][h]
            if (qso_s1 in hsla_qso_name) and (qso_s2 in hsla_qso_name):
                ind_in_hsla[i] = h
                #print(qso_name, hsla_qso_name)
    # print(ind_in_hsla[ind_in_hsla == 0])
    bordoloi14_yz['qso_hsla_name'] = hsla_tb['Target Name'][ind_in_hsla]
    bordoloi14_yz['qso_ra_deg'] = hsla_tb['RA'][ind_in_hsla]
    bordoloi14_yz['qso_dec_deg'] = hsla_tb['DEC'][ind_in_hsla]
    bordoloi14_yz['qso_hsla_snr'] = np.around(hsla_tb['median S/N'][ind_in_hsla], decimals=1)

    bordoloi14_yz['SFR_ref'] = '?'

    # new order
    neworder = ['gal_name', 'RA', 'DEC', 'gal_ra_deg', 'gal_dec_deg', 'gal_z','phi',
                'logSFR_flg', 'logSFR(Msun/yr)', 'SFR_ref', 'logM*', 'M*_ref', 'R200m_kpc',
                'L/L*', 'impact_para_kpc',  'qso_name',
                'qso_hsla_name', 'qso_ra_deg','qso_dec_deg','qso_hsla_snr',
                'logN_flg_CIV','logN_CIV', 'elogN_CIV',
                'Wr_flg_CIV','Wr_CIV','eWr_CIV']
    bordoloi14_yz_neworder = bordoloi14_yz[neworder].copy()

    # only return table with certain snr and mstar
    snr_cut = bordoloi14_yz_neworder['qso_hsla_snr'] >= snr_limit
    mstar_cut = bordoloi14_yz_neworder['logM*'] <= logmstar_limit
    bordoloi14_yz_final = bordoloi14_yz_neworder[snr_cut & mstar_cut].copy()

    if print_info == True:
        print('*'*30)
        print("Read table1 from Bordoloi+2014, only CIV, these changes are made: ")
        print("- updated M* from Salpeter to Kroupa")
        print("- Mhalo from M* based on Moster+2010")
        print("- Rvir based on r200, 200 times matter density")
        print("- Wr non-detections from 1sigma to 3 sigma")
        print("- logN non-detections from 2sigma to 3 sigma")
        print("- added SNR from HSLA 2nd data release.")
        print("- only those with SNR>={}, logMstar<={}".format(snr_limit, logmstar_limit))
    return bordoloi14_yz_final

def read_zheng19_wlm(print_info=True):
    from astropy.table import Table, Column
    wlm_tb = Table()
    gal_col = Column(['Zheng+2019b'], name='References')
    wlm_tb.add_column(gal_col, index=0)
    wlm_tb['gal_name'] = 'WLM'
    wlm_tb['gal_ra_deg'] = 0.4921
    wlm_tb['gal_dec_deg'] = -15.4611
    wlm_tb['gal_vlsr_km/s'] = -124.82
    wlm_tb['gal_dmpc'] = 9.8
    wlm_tb['gal_z'] = np.nan
    wlm_tb['logM*'] = 7.42
    wlm_tb['M*_ref'] = 'Cook+2014'
    wlm_tb['logMHI'] = 7.79
    wlm_tb['MHI_ref'] = 'Putman+2021'

    wlm_tb['R200m_kpc'] = 98.0 # checked, using the new B13-GK14 relation
    wlm_tb['logSFR_flg'] = '=' # check
    wlm_tb['logSFR(Msun/yr)'] = -2.43 # GALEX, FUV, Lee+2011
    wlm_tb['SFR_ref'] = 'FUV,Lee11'

    wlm_tb['qso_hsla_name'] = 'PHL2525'
    wlm_tb['qso_ra_deg'] = 0.1018
    wlm_tb['qso_dec_deg'] = -12.7633
    wlm_tb['qso_hsla_snr'] = 18.0 # HSLA
    wlm_tb['qso_snr_stis'] = np.nan
    wlm_tb['impact_para_kpc'] = 46.6

    wlm_tb['Wr_flg_HI'] = ' '*8
    wlm_tb['Wr_HI'] = np.nan
    wlm_tb['eWr_HI'] = np.nan
    wlm_tb['logN_flg_HI'] = '<=(3sig)'
    wlm_tb['logN_HI'] = 18.02
    wlm_tb['elogN_HI'] = np.nan

    # logN values based on the paper,
    # Wr values are recalculated for the normalized spectra
    # over [-60, 5] of the wlm_vhel
    # see WLM_PHL2525_aod_bin3.pdf
    # SiII 1260
    wlm_tb['Wr_flg_SiII'] = '='
    wlm_tb['Wr_SiII'] = 118.7
    wlm_tb['eWr_SiII'] = 5
    wlm_tb['logN_flg_SiII'] = '='
    wlm_tb['logN_SiII'] = 12.97
    wlm_tb['elogN_SiII'] = 0.07

    # SiIII 1206
    wlm_tb['Wr_flg_SiIII'] = '='
    wlm_tb['Wr_SiIII'] = 233.1
    wlm_tb['eWr_SiIII'] = 4.9
    wlm_tb['logN_flg_SiIII'] = '='
    wlm_tb['logN_SiIII'] = 13.30
    wlm_tb['elogN_SiIII'] =0.09

    # SiIV 1393
    wlm_tb['Wr_flg_SiIV'] = '='
    wlm_tb['Wr_SiIV'] = 72.6
    wlm_tb['eWr_SiIV'] = 8.7
    wlm_tb['logN_flg_SiIV'] = '='
    wlm_tb['logN_SiIV'] = 12.95
    wlm_tb['elogN_SiIV'] = 0.10

    # CII 1334
    wlm_tb['Wr_flg_CII'] = '='
    wlm_tb['Wr_CII'] = 169.0
    wlm_tb['eWr_CII'] = 5.7
    wlm_tb['logN_flg_CII'] = '='
    wlm_tb['logN_CII'] = 13.93
    wlm_tb['elogN_CII'] = 0.04

    # CIV 1548
    wlm_tb['Wr_flg_CIV'] = '='
    wlm_tb['Wr_CIV'] = 151.4
    wlm_tb['eWr_CIV'] = 10.8
    wlm_tb['logN_flg_CIV'] = '='
    wlm_tb['logN_CIV'] = 13.67
    wlm_tb['elogN_CIV'] = 0.09

    wlm_tb['Wr_flg_OVI'] = ' '*8
    wlm_tb['Wr_OVI'] = np.nan
    wlm_tb['eWr_OVI'] = np.nan
    wlm_tb['logN_flg_OVI'] = '='
    wlm_tb['logN_OVI'] = np.nan
    wlm_tb['elogN_OVI'] = np.nan

    if print_info == True:
        print('*'*30)
        print("Read table1 from Zheng+2019b, qso SNR from HSLA")
        print("including HI21cm, SiII, SiIII, SiIV, CIV, CII")
        print('logN from table1. No Wr given in paper. ')
        print("Calculated Wr directly from spectra instead.")

    return wlm_tb

###### now reading in IC1613 #####
def read_zheng20_ic1613(snr_limit=-999, print_info=True):
    from astropy.table import Table, Column
    import numpy as np
    ic1613_tb = Table()
    gal_col = Column(['Zheng+2020b']*6, name='References')
    ic1613_tb.add_column(gal_col, index=0)

    ic1613_tb['gal_name'] = ['IC1613']*6
    ic1613_tb['gal_ra_deg'] = [16.1992]*6
    ic1613_tb['gal_dec_deg'] = [2.1333]*6
    ic1613_tb['gal_vlsr_km/s'] = [-236.39]*6
    ic1613_tb['gal_dmpc'] = [0.76]*6
    ic1613_tb['gal_z'] = [np.nan]*6

    ic1613_tb['logM*'] = [np.around(np.log10(1e8), decimals=2)]*6 # McConnachie 2012, checked
    ic1613_tb['M*_ref'] = ['McConnachie2012']*6
    ic1613_tb['logMHI'] = 7.81 # Putman 21
    ic1613_tb['MHI_ref'] = ['Putman+2021']*6
    ic1613_tb['R200m_kpc'] = [123.6]*6 # using the new b13-GK14 relation
    ic1613_tb['logSFR_flg'] = ['=']*6 # checked
    ic1613_tb['logSFR(Msun/yr)'] = -2.23
    ic1613_tb['SFR_ref'] = 'FUV,Lee11'


    ic1613_tb['qso_hsla_name'] = ['LBQS-0100+0205', 'LBQS-0101+0009', '2MASX J01022632-0039045',
                                  'PG 0044+030', 'HB89-0107-025-NED05', 'LBQS-0107-0235']
    ic1613_tb['qso_ra_deg'] = [15.8041, 15.9281, 15.6097,
                               11.7746, 17.5677, 17.5547]
    ic1613_tb['qso_dec_deg'] = [2.3528, 0.4270, -0.6513,
                                3.3319, -2.3142, -2.3314]
    ic1613_tb['qso_hsla_snr'] = [7.6, 7.5, 8.0,  # checked using HSLA similar algorithm
                                 7.6, 11.4, 10.8] # last 3 checked with HSLA
    ic1613_tb['qso_snr_stis'] = [np.nan]*6
    ic1613_tb['impact_para_kpc'] = [6.3, 22.8, 37.6,
                                    60.7, 61.2, 61.4]


    ic1613_tb['Wr_flg_HI'] = [' '*8]*6
    ic1613_tb['Wr_HI'] = [np.nan]*6
    ic1613_tb['eWr_HI'] = [np.nan]*6
    ic1613_tb['logN_flg_HI'] = [' '*8]*6
    ic1613_tb['logN_HI'] = [np.nan]*6
    ic1613_tb['elogN_HI'] = [np.nan]*6

    # if there are multiple absorbers, we show the combined values here and errors with proporgation
    ic1613_tb['Wr_flg_SiII'] = ['<=(3sig)', '=', '=', '=', '<=(3sig)', '='] # note that Wr is in 1193
    #ic1613_tb['Wr_SiII'] = [49.8, 77.8, 183.3, 102.5, 38.7, 111.0]
    #ic1613_tb['eWr_SiII'] = [np.nan, 16.1, 16.9, 25.2, np.nan, 17.7]
    ic1613_tb['logN_flg_SiII'] = ['<=(3sig)', '=', '=', '=', '<=(3sig)', '=']
    ic1613_tb['logN_SiII'] = [13.03, 13.19, 13.53, 13.31, 13.04, 13.30]
    ic1613_tb['elogN_SiII'] = [np.nan, 0.06, 0.03, 0.1, np.nan, 0.07]

    from yztools.llist import load_atomic_data
    ion_line = 'SiII1260'
    line_info = load_atomic_data.morton03('SiII')
    line_ind = np.where(np.asarray(line_info['line_name']) == ion_line)[0][0]
    # print(line_info['line_name'][line_ind])
    line_f = line_info['f'][line_ind]
    line_lambda = line_info['lambda'][line_ind]
    N = 10**ic1613_tb['logN_SiII']
    eN = ic1613_tb['elogN_SiII'] * N * np.log(10)
    # based on Savage+1996, eq 3
    line_Wr = N * line_f * (line_lambda**2) / 1.13e17
    line_eWr = eN * line_f * (line_lambda**2) / 1.13e17
    ic1613_tb['Wr_SiII'] = np.around(line_Wr, decimals=2)
    ic1613_tb['eWr_SiII'] = np.around(line_eWr, decimals=2)

    ic1613_tb['Wr_flg_SiIII'] = ['=', '=', '=', '=', '=', '=']
    ic1613_tb['Wr_SiIII'] = [126.2, 271.3, 299.4, 318.9, 100.7, 105.1]
    ic1613_tb['eWr_SiIII'] = [16.6, 23.1, 17.5, 27.6, 14.0, 15.4]
    ic1613_tb['logN_flg_SiIII'] = ['=', '=', '=', '=', '=', '=']
    ic1613_tb['logN_SiIII'] = [12.96, 13.30, 13.38, 13.39, 12.76, 12.77]
    ic1613_tb['elogN_SiIII'] = [0.06, 0.05, 0.06, 0.07, 0.06, 0.06]

    ic1613_tb['Wr_flg_SiIV'] = ['=', '<=(3sig)', '=', '<=(3sig)', '<=(3sig)', '<=(3sig)']
    ic1613_tb['Wr_SiIV'] = [65.6, 35.1, 84.0, 68.7, 34.2, 42.6]
    ic1613_tb['eWr_SiIV'] = [13.2, np.nan, 15.1, np.nan, np.nan, np.nan]
    ic1613_tb['logN_flg_SiIV'] = ['=', '<=(3sig)', '=', '<=(3sig)', '<=(3sig)', '<=(3sig)']
    ic1613_tb['logN_SiIV'] = [13.00, 12.79, 13.02, 13, 13.02, 12.80]
    ic1613_tb['elogN_SiIV'] = [0.07, np.nan, 0.07, np.nan, np.nan, np.nan]

    ic1613_tb['Wr_flg_CII'] = ['<=(3sig)', '=', '=', '=', '<=(3sig)', '='] # all good
    ic1613_tb['Wr_CII'] = [48.9, 189.9, 390.7, 116.7, 29.7, 71.3]
    ic1613_tb['eWr_CII'] = [np.nan, 14.8, 13.7, 15.8, np.nan, 14.0]
    ic1613_tb['logN_flg_CII'] = ['<=(3sig)', '=', '=', '=', '<=(3sig)', '=']
    ic1613_tb['logN_CII'] = [13.74, 14.21, 14.52, 13.89, 13.38, 13.60]
    ic1613_tb['elogN_CII'] = [np.nan, 0.05, 0.04, 0.09, np.nan, 0.08]

    ic1613_tb['Wr_flg_CIV'] = ['=',  '=', '=', ' '*8,  '=', '=']
    ic1613_tb['Wr_CIV'] = [  117.3, 144.4, 205.7, np.nan, 121.5, 221.5]
    ic1613_tb['eWr_CIV'] = [  27.5,  23.2, 22.4, np.nan, 15.9, 10.7]
    ic1613_tb['logN_flg_CIV'] = ['=',  '=', '=', ' '*8,  '=', '=']
    ic1613_tb['logN_CIV'] = [ 13.57, 13.64, 13.78, np.nan, 13.56, 13.91]
    ic1613_tb['elogN_CIV'] = [ 0.09,  0.07, 0.05, np.nan, 0.06, 0.04]

    ic1613_tb['Wr_flg_OVI'] = [' '*8]*6
    ic1613_tb['Wr_OVI'] = [np.nan]*6
    ic1613_tb['eWr_OVI'] = [np.nan]*6
    ic1613_tb['logN_flg_OVI'] = [' '*8]*6
    ic1613_tb['logN_OVI'] = [np.nan]*6
    ic1613_tb['elogN_OVI'] = [np.nan]*6

    if print_info == True:
        print('*'*30)
        print("Read table1 from Zheng+2020b, qso SNR from HSLA when available")
        print("including SiII, SiIII, SiIV, CIV, CII")
        print("- only those with SNR>={}".format(snr_limit))
        print('- Wr of SiII converted from logN from SiII 1193 to SiII 1260.')

    return ic1613_tb[ic1613_tb['qso_hsla_snr']>=snr_limit]

def read_qu22_NGC3109_SexA_SexB(snr_limit=-999, print_info=True):
    ### now Qu & Bregman ####
    from astropy.table import Table, Column
    import numpy as np
    ddir = '/Users/Yong/Dropbox/GitRepo/yztools/read_datasets/cgm_tables'
    qu22_ew_tb = Table.read(ddir+'/Qu22_table_EW_vpm50_EW3sig.dat', format='ascii')
    qu22_logN_tb = Table.read(ddir+'/Qu22_table_N_vpm50_N3sig.dat', format='ascii')

    qu22_tb = Table()
    gal_col = Column(['Qu&Bregman2022']*6, name='References')
    qu22_tb.add_column(gal_col, index=0)
    # qu22_tb['gal_name'] = ['NGC3109', 'NGC3109', 'NGC3109',
    #                        'Sextans A', 'Sextans A', 'Sextans B']
    qu22_tb['gal_name'] = qu22_ew_tb['Galaxy']
    qu22_tb['gal_ra_deg'] = [150.78, 150.78, 150.78, 152.7533, 152.7533, 150.0004]
    qu22_tb['gal_dec_deg'] = [-26.16, -26.16, -26.16, -4.6928, -4.6928, 5.3322]
    qu22_tb['gal_vlsr_km/s'] = [393.7, 393.7, 393.7, 317.3, 317.3, 294.0]
    qu22_tb['gal_dmpc'] = [1.34, 1.34, 1.34, 1.45, 1.45, 1.43]
    qu22_tb['gal_z'] = [np.nan]*6

    qu22_tb['logM*'] = [8.28, 8.28, 8.28, 7.35, 7.35, 7.56] # cook+2014
    qu22_tb['M*_ref'] = ['Cook+2014']*6
    qu22_tb['logMHI'] = [8.65, 8.65,8.65, 7.89, 7.89, 7.71] # Putman21
    qu22_tb['MHI_ref'] = ['Putman+2021']*6

    from yztools.calc_r200 import calc_r200
    from yztools.mstar2mhalo import mstar2mhalo
    log_mhalo = mstar2mhalo(qu22_tb['logM*'], ref='B13-GK14')
    mhalo = 10**log_mhalo
    r200m = calc_r200(mhalo)[1].value
    qu22_tb['R200m_kpc'] = np.around(r200m, decimals=1)
    qu22_tb['logSFR_flg'] = ['=']*6
    qu22_tb['logSFR(Msun/yr)'] = [-1.71, -1.71, -1.71, -2.12, -2.12, -2.58]
    qu22_tb['SFR_ref'] = 'FUV,Lee11'

    qu22_tb['qso_hsla_name'] = ['CTS M00.02', 'ESO 499-G 041', '1RXS J1015-2748',
                                'MARK 1253', 'PG 1011-040', 'PG 1001+054']
    # CTS M00.02, hsla_name = '2MASX-J10053271-2417161'
    # 'ESO 499-G 041', ESO-499-G-041
    # 1RXS J1015-2748, 2MASS-J10155924-2748289
    # MARK 1253, MRK-1253
    # PG 1011-040, PG1011-040
    #'PG 1001+054', 2XMM-J100420.0+051300

    qu22_tb['qso_ra_deg'] = [151.386, 151.481, 153.997, 154.887, 153.586, 151.084]
    qu22_tb['qso_dec_deg'] = [-24.2878, -23.0568, -27.808, -3.3375, -4.31124, 5.2168]
    qu22_tb['qso_hsla_snr'] = [6.7, 4.3, 5.3, 8.6, 29.7, 15.6]
    qu22_tb['qso_snr_stis'] = [np.nan]*6
    qu22_tb['impact_para_kpc'] = [44, 73, 75, 64, 21, 28]

    # column densities are 3 sigma upper limits over +/- 50 km/s in their table
    # EW are 1 sigma upper limits in their table
    ion = 'HI'
    qu22_tb['Wr_flg_'+ion] = [' '*8]*6
    qu22_tb['Wr_'+ion] = [np.nan]*6
    qu22_tb['eWr_'+ion] = [np.nan]*6
    qu22_tb['logN_flg_'+ion] = [' '*8]*6
    qu22_tb['logN_'+ion] = [np.nan]*6
    qu22_tb['elogN_'+ion] = [np.nan]*6

    ion = 'SiII'
    qu22_tb['Wr_flg_'+ion] = ['<=(3sig)', '<=(3sig)', '<=(3sig)', '<=(3sig)', '<=(3sig)', '<=(3sig)']
    qu22_tb['Wr_'+ion] = np.around(qu22_ew_tb['SiII 1260']*3, decimals=1)
    qu22_tb['eWr_'+ion] = [np.nan]*6
    qu22_tb['logN_flg_'+ion] = ['<=(3sig)', '<=(3sig)', '<=(3sig)', '<=(3sig)', '<=(3sig)', '<=(3sig)']
    qu22_tb['logN_'+ion] = qu22_logN_tb['SiII 1260']
    qu22_tb['elogN_'+ion] = [np.nan]*6

    ion = 'SiIII'
    qu22_tb['Wr_flg_'+ion] = ['<=(3sig)', '<=(3sig)', '<=(3sig)', '<=(3sig)', '<=(3sig)', '<=(3sig)']
    qu22_tb['Wr_'+ion] = np.around(qu22_ew_tb['SiIII 1206']*3, decimals=1)
    qu22_tb['eWr_'+ion] = [np.nan]*6
    qu22_tb['logN_flg_'+ion] = ['<=(3sig)', '<=(3sig)', '<=(3sig)', '<=(3sig)', '<=(3sig)', '<=(3sig)']
    qu22_tb['logN_'+ion] = qu22_logN_tb['SiIII 1206']
    qu22_tb['elogN_'+ion] = [np.nan]*6

    ion = 'SiIV'
    qu22_tb['Wr_flg_'+ion] = ['<=(3sig)', '<=(3sig)', '<=(3sig)', '<=(3sig)', '<=(3sig)', '<=(3sig)']
    qu22_tb['Wr_'+ion] = np.around(qu22_ew_tb['SiIV 1393']*3, decimals=1)
    qu22_tb['eWr_'+ion] = [np.nan]*6
    qu22_tb['logN_flg_'+ion] = ['<=(3sig)', '<=(3sig)', '<=(3sig)', '<=(3sig)', '<=(3sig)', '<=(3sig)']
    qu22_tb['logN_'+ion] = qu22_logN_tb['SiIV 1393']
    qu22_tb['elogN_'+ion] = [np.nan]*6

    ion = 'CII'
    qu22_tb['Wr_flg_'+ion] = ['<=(3sig)', '<=(3sig)', '<=(3sig)', ' '*8, ' '*8,  ' '*8]
    qu22_tb['Wr_'+ion] = np.around(qu22_ew_tb['CII 1334']*3, decimals=1)
    qu22_tb['eWr_'+ion] = [np.nan]*6
    qu22_tb['logN_flg_'+ion] = ['<=(3sig)', '<=(3sig)', '<=(3sig)', ' '*8, ' '*8,  ' '*8]
    qu22_tb['logN_'+ion] = qu22_logN_tb['CII 1334']
    qu22_tb['elogN_'+ion] = [np.nan]*6

    ion = 'CIV'
    qu22_tb['Wr_flg_'+ion] = ['<=(3sig)', '<=(3sig)', '<=(3sig)', '<=(3sig)', '=', '<=(3sig)']
    qu22_tb['Wr_'+ion] = np.around(qu22_ew_tb['CIV 1548']*3, decimals=1)
    qu22_tb['eWr_'+ion] = [np.nan]*6
    qu22_tb['logN_flg_'+ion] = ['<=(3sig)', '<=(3sig)', '<=(3sig)', '<=(3sig)', '=', '<=(3sig)']
    qu22_tb['logN_'+ion] = qu22_logN_tb['CIV 1548']
    qu22_tb['elogN_'+ion] = [np.nan]*6

    qu22_tb['logN_'+ion][4] = 13.04
    qu22_tb['elogN_'+ion][4] = 0.08

    from yztools.llist import load_atomic_data
    ion_line = 'CIV1548'
    line_info = load_atomic_data.morton03('CIV')
    line_ind = np.where(np.asarray(line_info['line_name']) == ion_line)[0][0]
    # print(line_info['line_name'][line_ind])
    line_f = line_info['f'][line_ind]
    line_lambda = line_info['lambda'][line_ind]
    N = 10**qu22_tb['logN_'+ion][4]
    eN = qu22_tb['elogN_'+ion][4] * N * np.log(10)
    # based on Savage+1996, eq 3
    qu22_tb['Wr_'+ion][4] = np.around(N * line_f * (line_lambda**2) / 1.13e17, decimals=1)
    qu22_tb['eWr_'+ion][4] = np.around(eN * line_f * (line_lambda**2) / 1.13e17, decimals=1)

    ion = 'OVI'
    qu22_tb['Wr_flg_'+ion] = [' '*8]*6
    qu22_tb['Wr_'+ion] = [np.nan]*6
    qu22_tb['eWr_'+ion] = [np.nan]*6
    qu22_tb['logN_flg_'+ion] = [' '*8]*6
    qu22_tb['logN_'+ion] = [np.nan]*6
    qu22_tb['elogN_'+ion] = [np.nan]*6

    if print_info == True:
        print('*'*30)
        print("Read table1 from Qu&Bregman 2022, qso SNR from HSLA when available")
        print("including SiII1260, SiIII1206, SiIV1393, CIV1548, CII1334")
        print("EW and logN non-detections updated by Qu for +/-50 km/s, and 3sigma values")
        print("EW for CIV detection was converted from logN values for CIV1548")
        print("- only those with SNR>={}".format(snr_limit))
        print('Still need to double check M* and logN')

    return qu22_tb[qu22_tb['qso_hsla_snr']>=snr_limit]

def read_zheng_newdw_tb(snr_limit=-999, logmstar_limit=999, print_info=True,
                           skip_zheng_lgs3_det=True, skip_zheng_ddo210_det=True):
    from astropy.table import Table, Column

    """
    read AOD data from recent dwarf archival analysis

    skip_zheng_lgs3_det: skip two QSOs "PG0052+251" and "RXJ00537+2232", there are
             detections along these two sightlines, but not sure if they are related
             to LGS3 given that there are some loose HI clouds nearby. skip in this
             compilation.
    skip_zheng_ddo210_det: skip MARK509 near DDO190. Very high SNR sightline, but not sure
             if the little dip is actually detection, voigt profile fitting isn't great. skip
             this one in compiliation.

    History:
    April 26, 2022. YZ. UCB.
    """

    gal_qso_tb = Table.read('/Users/Yong/Dropbox/DwarfsCGM_8Mpc/tables/tb_final_gal_qso_2022-03-31.csv',
                            format='ascii')
    if skip_zheng_lgs3_det==True:
        ind_lgs3 = gal_qso_tb['gal_name'] == 'LGS3'
        ind_lgs3_qso1 = gal_qso_tb['qso_name'] == 'PG0052+251'
        ind_lgs3_qso2 = gal_qso_tb['qso_name'] == 'RXJ00537+2232'
        ind_lgs3_skip = (ind_lgs3 & ind_lgs3_qso1) | ((ind_lgs3 & ind_lgs3_qso2))
        keep = ~ind_lgs3_skip
        gal_qso_tb = gal_qso_tb[keep].copy()

    if skip_zheng_ddo210_det==True:
        ind_ddo = gal_qso_tb['gal_name'] == 'DDO210'
        ind_ddo_qso1 = gal_qso_tb['qso_name'] == 'MARK509'
        ind_ddo_skip = ind_ddo & ind_ddo_qso1 # | ((ind_lgs3 & ind_lgs3_qso2))
        keep = ~ind_ddo_skip
        gal_qso_tb = gal_qso_tb[keep].copy()

    zheng_newdw_tb = gal_qso_tb.copy()
    zheng_newdw_tb.remove_columns(['filters', 'proposal_pi', 'proposal_id'])
    zheng_newdw_tb.rename_column('qso_name', 'qso_hsla_name')
    zheng_newdw_tb.rename_column('qso_snr', 'qso_hsla_snr')
    zheng_newdw_tb.rename_column('gal_vlsr_kms', 'gal_vlsr_km/s')
    zheng_newdw_tb.rename_column('gal_dkpc', 'gal_dmpc')
    zheng_newdw_tb['gal_dmpc'] = np.around(zheng_newdw_tb['gal_dmpc']/1e3, decimals=2)

    zheng_newdw_tb['qso_snr_stis'] = np.nan
    zheng_newdw_tb['gal_z'] = np.nan

    ref_col = Column(['ThisWork']*len(gal_qso_tb), name='References')
    zheng_newdw_tb.add_column(ref_col, index=0)
    # zheng_newdw_tb.show_in_notebook()

    # pull galaxy information K13
    #ddir = '/Users/Yong/Dropbox/databucket/11HUGS_LVL'
    #k13_tb = Table.read(ddir+'/11HUGS_K13online_Lee11_K08_combined_final.csv', format='ascii')

    #gal_logMstar = np.zeros(len(zheng_newdw_tb))+np.nan
    #gal_logMHI = np.zeros(len(zheng_newdw_tb))+np.nan
    #gal_logsfr = np.zeros(len(zheng_newdw_tb))+np.nan
    #gal_logsfr_ref = [' ']*len(zheng_newdw_tb)
    #for i in range(len(zheng_newdw_tb)):
    #    gal_name = zheng_newdw_tb['gal_name'][i]
    #    ind = np.where(k13_tb['Name'] == gal_name)[0][0]
    #    gal_logMstar[i] = k13_tb['logM*_Ks(K13)'][ind]
    #    gal_logMHI[i] = k13_tb['logMHI(K13)'][ind]

    #    logsfr_fuv = k13_tb['logSFR_Msun/yr (FUV)'][ind]
    #    logsfr_ha = k13_tb['logSFR_Msun/yr (LHa)'][ind]
    #    if np.isfinite(logsfr_fuv):
    #        gal_logsfr[i] = logsfr_fuv
    #        gal_logsfr_ref[i] = 'FUV'
    #    else:
    #        gal_logsfr[i] = logsfr_ha
    #        gal_logsfr_ref[i] = 'Halpha'

    # pull galaxy information, see
    # https://docs.google.com/spreadsheets/d/1qGMlRgKu3ot7U7PRQqoHU0y9Ib28MRwb8k4j5Xkw234/edit#gid=1017337321

    gal_info = Table.read('/Users/Yong/Dropbox/DwarfsCGM_8Mpc/tables/tb_final_gal_property_2022-05-12.csv',
                          format='ascii')
    gal_logMstar = np.zeros(len(zheng_newdw_tb))+np.nan
    gal_logMstar_ref = [' ']*len(zheng_newdw_tb)

    gal_logMHI = np.zeros(len(zheng_newdw_tb))+np.nan
    gal_logMHI_ref = [' ']*len(zheng_newdw_tb)

    gal_logsfr = np.zeros(len(zheng_newdw_tb))+np.nan
    gal_logsfr_ref = [' ']*len(zheng_newdw_tb)

    for i in range(len(zheng_newdw_tb)):
        gal_name = zheng_newdw_tb['gal_name'][i]
        ind = np.where(gal_info['gal_name'] == gal_name)[0][0]
        gal_logMstar[i] = gal_info['logM*_Msun'][ind]
        gal_logMstar_ref[i] = gal_info['ref_M*'][ind]

        gal_logMHI[i] = gal_info['logMHI'][ind]
        gal_logMHI_ref[i] = gal_info['ref_MHI'][ind]

        gal_logsfr[i] = gal_info['logSFR'][ind]
        gal_logsfr_ref[i] = gal_info['ref_SFR'][ind]

    zheng_newdw_tb['logM*'] = gal_logMstar
    zheng_newdw_tb['M*_ref'] =  gal_logMstar_ref
    zheng_newdw_tb['logMHI'] = gal_logMHI
    zheng_newdw_tb['MHI_ref'] =  gal_logMHI_ref
    zheng_newdw_tb['logSFR_flg'] = '='
    zheng_newdw_tb['logSFR(Msun/yr)'] = gal_logsfr
    zheng_newdw_tb['SFR_ref'] = gal_logsfr_ref

    # calculate halo mass
    from yztools.mstar2mhalo import mstar2mhalo
    log_gal_mhalo = mstar2mhalo(zheng_newdw_tb['logM*'], ref='B13-GK14')
    gal_mhalo = 10**log_gal_mhalo

    # Impact parameters
    from yztools.calc_r200 import calc_r200
    gal_r200m = calc_r200(gal_mhalo)[1].value
    zheng_newdw_tb['R200m_kpc'] = np.around(gal_r200m, decimals=1)

    ion_lines = ['SiII1260', 'SiIII1206', 'SiIV1393', 'CII1334', 'CIV1548']
    from astropy.table import vstack
    for i in range(len(zheng_newdw_tb)):
        gal_name = zheng_newdw_tb['gal_name'][i]
        if gal_name == '[KK2000]03':
            gal_name = 'KK2000-03'
        qso_name = zheng_newdw_tb['qso_hsla_name'][i]
        this_aod_tb = read_gal_aod_info(gal_name, qso_name, ion_lines=ion_lines)
        if i == 0:
            ion_aod_tb = this_aod_tb
        else:
            ion_aod_tb = vstack([ion_aod_tb, this_aod_tb])

    for ion in ['HI', 'SiII', 'SiIII', 'SiIV', 'CII', 'CIV', 'OVI']:
        zheng_newdw_tb['Wr_flg_'+ion] = ion_aod_tb['Wr_flg_'+ion]
        zheng_newdw_tb['Wr_'+ion] = ion_aod_tb['Wr_'+ion]
        zheng_newdw_tb['eWr_'+ion] = ion_aod_tb['eWr_'+ion]
        zheng_newdw_tb['logN_flg_'+ion] = ion_aod_tb['logN_flg_'+ion]
        zheng_newdw_tb['logN_'+ion] = ion_aod_tb['logN_'+ion]
        zheng_newdw_tb['elogN_'+ion] = ion_aod_tb['elogN_'+ion]

    snr_cut = zheng_newdw_tb['qso_hsla_snr'] >= snr_limit
    mstar_cut = zheng_newdw_tb['logM*'] <= logmstar_limit
    return zheng_newdw_tb[snr_cut & mstar_cut]

def read_gal_aod_info(gal_name, qso_name,
                      ion_lines=['SiII1260', 'SiIII1206', 'SiIV1393', 'CII1334', 'CIV1548']):
    import re
    if '{}_{}'.format(gal_name, qso_name) in ['DDO210_MARK509', 'DDO187_NGC-5548']:
        nbin = 'bin1'
    else:
        nbin = 'bin3'
    aod_file = glob('/Users/Yong/Dropbox/DwarfsCGM_8Mpc/step1_aod_results/{}_{}_aod_{}.txt'.format(gal_name, qso_name, nbin))
    line_aod_tb = Table()
    for ion in ['HI', 'OVI']:
        line_aod_tb.add_column(Column([' '*8], name='Wr_flg_'+ion))
        line_aod_tb.add_column(Column([' '*8], name='logN_flg_'+ion))
        line_aod_tb.add_column(Column([np.nan], name='Wr_'+ion))
        line_aod_tb.add_column(Column([np.nan], name='eWr_'+ion))
        line_aod_tb.add_column(Column([np.nan], name='logN_'+ion))
        line_aod_tb.add_column(Column([np.nan], name='elogN_'+ion))

    if len(aod_file) == 0:
        print('Can\'t find {}_{}_aod_bin3.txt'.format(gal_name, qso_name))
        for iline in ion_lines:
            ion = re.split('(\d+)', iline)[0]
            line_aod_tb.add_column(Column([' '*8], name='Wr_flg_'+ion))
            line_aod_tb.add_column(Column([' '*8], name='logN_flg_'+ion))
            line_aod_tb.add_column(Column([np.nan], name='Wr_'+ion))
            line_aod_tb.add_column(Column([np.nan], name='eWr_'+ion))
            line_aod_tb.add_column(Column([np.nan], name='logN_'+ion))
            line_aod_tb.add_column(Column([np.nan], name='elogN_'+ion))
    else:
        aod_file = aod_file[0]
        aod_tb = Table.read(aod_file, format='ascii')
        for iline in ion_lines:
            ion = re.split('(\d+)', iline)[0]
            ind_line = np.where(aod_tb['line'] == iline)[0][0]

            if aod_tb['flag'].mask[ind_line] == True:
                line_aod_tb['Wr_flg_'+ion] = ' '*8
            else:
                line_aod_tb['Wr_flg_'+ion] = aod_tb['flag'][ind_line]

            line_aod_tb['logN_flg_'+ion] = line_aod_tb['Wr_flg_'+ion]
            if line_aod_tb['Wr_flg_'+ion] == '<=(3sig)':
                line_aod_tb['Wr_'+ion] = np.around(aod_tb['Wr(mA)(3sig)'][ind_line], decimals=1)
                line_aod_tb['eWr_'+ion] = np.nan
                line_aod_tb['logN_'+ion] = np.around(aod_tb['logN(3sig)'][ind_line], decimals=2)
                line_aod_tb['elogN_'+ion] = np.nan
            else:
                line_aod_tb['Wr_'+ion] = np.around(aod_tb['Wr(mA)'][ind_line], decimals=1)
                line_aod_tb['eWr_'+ion] = np.around(aod_tb['eWr(mA)'][ind_line], decimals=1)
                line_aod_tb['logN_'+ion] = np.around(aod_tb['logN'][ind_line], decimals=2)
                line_aod_tb['elogN_'+ion] = np.around(aod_tb['elogN'][ind_line], decimals=2)
    return line_aod_tb
