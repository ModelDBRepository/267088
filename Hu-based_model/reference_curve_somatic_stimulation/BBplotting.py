def BarlowRangeVarPlotter(rangeVars=['v'], vertical_stack=False, show=True, sections=False, ylims=[None, None], save_PNG=False):
    from neuron import h
    import numpy as np
    from matplotlib import pyplot as plt


    fig, (ax) = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=(10, 9))
    for rangeVar in rangeVars:
        yData_array = np.array([])
        xData_array = np.array([])
        x_prev = 0.0
        if sections==False:
            sections = h.allsec()
        for sec in sections:
        # for sec in [dend, soma, hill, ais:
            yvec = h.Vector()
            xvec = h.Vector()
            h.RangeVarPlot(rangeVar, sec(0), sec(1)).to_vector(yvec, xvec)
            ydata = np.array(yvec)
            xdata = x_prev + np.array(xvec)
            label = rangeVar+' '+sec.name()
            if vertical_stack == True:
                plt.plot(xdata/sec.L, ydata, label=label, linewidth=0, markersize=5, marker='o', alpha=0.6)
            else:
                x_prev = xdata[-1]
                plt.plot(xdata, ydata, label=label, linewidth=1, markersize=5, marker='o', alpha=0.6)
            
            yData_array = np.hstack((yData_array, ydata))
            xData_array = np.hstack((xData_array, xdata))
            # ax.set_ylabel(rangeVar, fontsize=40)
            ax.set_xlabel('position', fontsize=40)

    ax.set_title(r'$\Delta t = $'+str(h.t)+r'[ms]')
    ax.set_ylim(ylims)
    ax.legend()
    if show == True:
        if save_PNG==True:
            from time import time
            import os
            figurename='figure'+str(int(time()))+'.png'
            plt.savefig(figurename, dpi=600)
            os.system("open " + figurename)
            # os.open(figurename)
        else:
            plt.show()

def plot_CLSparams(sections):
    from BBplotting import BarlowRangeVarPlotter
    from mech_settings import mech_name
    params=[]
    params.append('gl'+'_'+mech_name)
    params.append('gnabar'+'_'+mech_name)
    params.append('gnabar2'+'_'+mech_name)
    params.append('gkbar'+'_'+mech_name)
    params.append('INaKmax'+'_'+mech_name)
    params.append('gnal'+'_'+mech_name)
    params.append('gkl'+'_'+mech_name)
    print(params)
    BarlowRangeVarPlotter(params, sections=sections, ylims=[0.0, None], save_PNG=False)

