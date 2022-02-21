from __future__ import print_function
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

def discrete_map(cmap, cmin, cmax, cnum):
    '''
    Self defined discrete map. 
    
    An example shown below:
    ========================
    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.set_axis_bgcolor('black')

    cmap, val2color, norm, bounds = discrete_map(plt.cm.bwr, 0.15, 0.4, 20)
    im = ax.imshow(m0, origin='lower', interpolation='nearest', 
                   cmap=cmap, norm=norm)
    # val2color.to_rgba() is used for scattered plot. 
    
    cax = fig.add_axes([0.83, 0.1, 0.03, 0.8])
    fig.colorbar(im, cax=cax, spacing='proportional', format='%.2f', 
                 ticks=bounds[::2], boundaries=bounds)
    fig.savefig('new.pdf')
    ========================
    Original code source: 
    http://stackoverflow.com/questions/14777066/matplotlib-discrete-colorbar
 
    Modified by Y. Zheng. 09, 2017. Columbia. 
    '''

    # cmap can be the name, or the actual map, e.g., plt.cm.bwr
    if type(cmap) == str: cmap = plt.get_cmap(cmap)

    # extract all colors from the map
    cmaplist = [cmap(i) for i in range(cmap.N)]

    # create the new map
    new_cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

    # define the bins and normalize
    bounds = np.linspace(cmin, cmax, cnum)
    norm = mpl.colors.BoundaryNorm(bounds, new_cmap.N)
    val2color = mpl.pyplot.cm.ScalarMappable(norm=norm, cmap=new_cmap)

    return new_cmap, val2color, norm, bounds

# customize the matplotlib colormap 
def yz_colors():

    # plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    # plt.rc('text', usetex=True)
    # setup the color for plotting 
    '''
    Some colors that Yong likes to use in her work. 
    To be updated. 
    '''   
 
    cs = [plt.cm.gist_rainbow(.75), plt.cm.gist_rainbow(.6), plt.cm.brg(0.85), 
          plt.cm.autumn(0.45), plt.cm.gist_rainbow(.05), plt.cm.gist_rainbow(.95), 
          plt.cm.gnuplot(.3), 'k']
    
    for ics in range(len(cs)):
        plt.plot([0, 1], [ics, ics], color=cs[ics], linewidth=4)
        plt.text(0.1, ics+0.2, '%d'%(ics))
    plt.ylim(-1, len(cs))
    
    return cs

def yz_normcmap(cmap, vmin, vmax):
    '''
    To customize the color map based on the dynamic range of the variables of interest.
    - Input: colormap name, e.g., plt.cm.Reds
    -        vmin, vmax: the minimum/maximum of the variable of interest
    - Output: 
             nmap: normalized colormap, used the same way as cmap
             val2color: a function that returns the color of a specific val; can be used 
                        to plot individual points. 
                        
    - Example: 
            > cmap = plt.cm.Reds
            > zmin, zmax = 0, 100
            > nmap, val2color = yz_normcmap(cmap, zmin, zmax)
            > plt.scatter(x, y, color=val2color.to_rgba(some_z_value)) # z must be zmin<z<zmax

    - History: updated as of 2016.10.03. Yong Zheng @ Columbia Astro.
    '''
    # to call: cmap = yz_normcmap(cmap, vmin, vmax)
    #          cmap.to_rgba(the_value_you_want_to_color)
    nmap = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    val2color = mpl.pyplot.cm.ScalarMappable(norm=nmap, cmap=cmap)

    return nmap, val2color
