import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['figure.figsize'] = (7,5)
plt.rcParams['font.size'] = 14
plt.rcParams['lines.linewidth'] = 2.
plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"

plt.style.use('seaborn-v0_8-colorblind')

# a function for formatting the axes' ticks:
def format_axes_ticks(fig):
    for i, ax in enumerate(fig.axes):
        ax.tick_params(labelsize=14)

cycle_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

def heat_grid ( xgrid, ygrid, x, y ) : 
    return np.ma.log10(np.histogram2d(x,y,(xgrid,ygrid))[0])

def contour_grid (x, y, **kwargs) :
    Z, xgrid, ygrid = np.histogram2d(x,y,**kwargs)
    Z = np.ma.log10(Z)
    xmesh = 0.5 * (xgrid[1:]+xgrid[:-1])
    ymesh = 0.5 * (ygrid[1:]+ygrid[:-1])
    X, Y = np.meshgrid(xmesh,ymesh); X, Y = X.T, Y.T
    return X, Y, Z

def plot_triangle (x, y, 
                   cgrid_kw = {}, contour_kw = {}, 
                   x1D_kw = {}, y1D_kw = {},
                   hist_kw = {},
                   filled = True, fig = None, axes = None) :
    """ Function from plotting the gridded scatter
    
    Parameters
    ----------
    x :
    y :
    cgrid_kw : dict
    contour_kw : dict
    x1D_kw : dict, 
    y1D_kw : dict,
    hist_kw : dict
    filled : bool
    fig : matplotlib.Figure
    axes : matplotlib.Axes
    
    Returns
    -------
    axes : matplotlib.Axes
    (X,Y,Z) : tuple
    (xhist, xbins) : tuple
    (yhist, ybins) : tuple
    """
    if fig is None and axes is None :
        fig, axes = plt.subplots(2,2, figsize=(7,7), 
                                 # sharey = True,
                                 # sharex = True,
                                 gridspec_kw = {'hspace':0.,
                                                'wspace':0.,
                                                'height_ratios':(0.5,1),
                                                'width_ratios':(1,0.5)})
    if axes is None :
        axes = fig.axes
        
    # Name axes    
    (xHist, off),(cont, yHist) = axes
    
    # Turn off upper-left panel
    off.axis('off')
    # Turn-off upper-right x-ticks
    xHist.set_xticks([])
    # Turn-off lower-left y-ticks
    yHist.set_yticks([])
        
    # Make contour-plot
    X, Y, Z = contour_grid(x,y, **cgrid_kw)
    if filled :
        cont.contourf(X,Y,Z,**contour_kw)
    else :
        c = cont.contour(X,Y,Z,**contour_kw)
        cont.clabel(c, inline=True, fontsize=8)
        
    # Plot histograms
    xhist, xbins = np.histogram(x, **x1D_kw)
    xhist = np.append(xhist, xhist[-1])
    xHist.fill_between(xbins, np.ma.log10(xhist), #np.zeros_like(xbins), 
                       step='mid', **hist_kw )
    xHist.set_ylabel('$\\log($counts$)$')
    
    yhist, ybins = np.histogram(y, **y1D_kw)
    yhist = np.append(yhist, yhist[-1])
    yHist.fill_betweenx(ybins, np.ma.log10(yhist), #np.zeros_like(xbins), 
                        step='mid', **hist_kw)
    yHist.set_xlabel('$\\log($counts$)$')
    
    return axes, (X, Y, Z), (xhist, xbins), (yhist, ybins)

def plot_comparison ( x, data, model, dist, axs, labels ) :
    axs[0].plot(np.ma.log10(x), np.ma.log10(data), 'o', label=labels[0])
    axs[0].plot(np.ma.log10(x), np.ma.log10(model), label=labels[1])
    axs[1].plot(np.ma.log10(x), 
                [dist([pd],[pm]) for pd,pm in zip(data,model)],
                'o', label=labels[2])
    return axs
