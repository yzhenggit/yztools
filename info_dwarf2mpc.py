def info_dwarf2mpc(gal_name, yes_print=True):
    import numpy as np
    from astropy.table import Table
    table_name = '/Users/Yong/Dropbox/databucket/putman20_dwarfs_2Mpc_yzedits.csv'
    dw_2mpc = Table.read(table_name, format='csv')
    ind = np.where(dw_2mpc['GalaxyName'] == gal_name)[0][0]
    print("Info from Putman20, mostly based on McConnachie12")
    if yes_print == True:
        for icol in dw_2mpc.colnames:
            print('%20s'%icol+'%20s'%dw_2mpc[icol][ind])
    return dw_2mpc[ind]

if __name__ == "__main__":
    import sys
    gal_name = sys.argv[1]
    res = info_dwarf2mpc(gal_name, yes_print=True)
