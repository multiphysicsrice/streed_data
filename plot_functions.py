# Libraries
import numpy as np
import matplotlib.pyplot as plt

# Plot setup function
def plot_setup(xlabel, ylabel, 
               title='None', legend=False, twin_axis=False, y2label='None',
               xscale='linear',yscale='linear',y2scale='linear',
               xticks=[], yticks=[], y2ticks=[], l_maj = 6, l_min = 4,
               xlim=[], ylim=[], y2lim=[],
               xclr = 'black', yclr = 'black', y2clr='black',
               fig_size=(5,4), font = 'Arial', fontsize = 10):

    plt.rcParams['font.family'] = font

    fig = plt.figure(figsize=fig_size);
    ax = fig.add_axes((0.3,0.3,0.5,0.5))
    ax.set_xlabel(xlabel,  color=xclr, fontsize=fontsize)
    ax.set_ylabel(ylabel,  color=yclr, fontsize=fontsize)
    ax.tick_params(axis='x', labelcolor=xclr, labelsize=fontsize)
    ax.tick_params(axis='y', labelcolor=yclr, labelsize=fontsize)
    ax.tick_params(which='major',length=l_maj,direction='in')
    ax.tick_params(which='minor',length=l_min,direction='in')

    if title != 'None':
        ax.set_title(title,fontsize=fontsize)

    if legend is True:
        ax.legend(loc='lower left', fontsize=fontsize-6)
        plt.rcParams["legend.fancybox"] = False
    
    if xticks != []:
        ax.set_xticks(xticks)

    if yticks != []:
        ax.set_yticks(yticks)

    if xlim != []:
        ax.set_xlim(xlim)
    
    if ylim != []:
        ax.set_ylim(ylim)

    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    
    if twin_axis == True:
        ax2 = ax.twinx()
        ax2.set_ylabel(y2label,  color=y2clr, fontsize=fontsize)
        ax2.tick_params(axis='y', labelcolor=y2clr, labelsize=fontsize)
        ax2.tick_params(which='major',length=l_maj, direction='in')
        ax2.tick_params(which='minor',length=l_min, direction='in')
        
        if y2ticks != []:
            ax2.set_yticks(y2ticks)

        if y2lim != []:
            ax2.set_ylim(y2lim)
        
        ax2.set_yscale(y2scale)
        return ax, ax2
    else:
        return ax

# Curve plotting function
def plot(axis, xvar, yvar, linestyle='-', color='black', label='_nolabel_',linewidth=1, markersize=5):
    axis.plot(xvar,yvar,linestyle,color=color,label=label,linewidth=linewidth,markersize=markersize,zorder=3)


