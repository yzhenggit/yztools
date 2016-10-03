from __future__ import print_function
import matplotlib as mpl
import matplotlib.pyplot as plt

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
