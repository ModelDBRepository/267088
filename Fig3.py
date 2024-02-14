"""
This script generates backproagation threshold and AP threshold figures from the main text.
(See end of script.)
"""


import os, sys
filepath = os.getcwd()+'/'
import itertools
import numpy as np
from sorcery import print_args
import sqlite3
import pandas as pd
pd.set_option('display.max_colwidth', None)
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
xlim = ylim = (None, None) #<---default (don't touch)
label_size = 20
legend_size = 17
line_width = 3
# matplotlib.rcParams['figure.titlesize'] = label_size
# matplotlib.rcParams['axes.titlesize'] = label_size
matplotlib.rcParams['xtick.labelsize'] = label_size
matplotlib.rcParams['ytick.labelsize'] = label_size
matplotlib.rcParams['axes.labelsize'] = label_size*1.3
matplotlib.rcParams['legend.fontsize'] = legend_size 
matplotlib.rcParams['legend.title_fontsize'] = legend_size 
matplotlib.rcParams['lines.linewidth'] = line_width

############################################################################
###matplotlib LaTeX settings################################################
rc('text', usetex=True) ; latex_preamble = r''
latex_preamble += r'\usepackage{stix}'
latex_preamble += r'\usepackage{amsmath}'
rc('text.latex', preamble=latex_preamble) ; 
rc('font', **{
    'family': 'STIXGeneral',
    'serif': ['Times'],
    'sans-serif': ['Helvetica']
})
# rc('font', **{
#     'family': 'sans-serif',
#     'sans-serif': ['Helvetica']  # You can replace 'Helvetica' with the specific sans-serif font you want to use
# })
############################################################################
matplotlib_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
# matplotlib_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
plotly_colors = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
plotly_symbols = ['o', 's', 'D', 'P', 'X', '^', 'v']




markers_dict = {}
markers_dict['0.2'] = 'p'
markers_dict['0.3'] = 'v'
markers_dict['0.4'] = 'P'
markers_dict['0.5'] = 'o'
markers_dict['0.6'] = 'X'
markers_dict['0.7'] = 'd'
markers_dict['0.8'] = '*'

### Get the default color cycle
default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
### Create colors_dict
colors_dict = {}
color_iterator = itertools.cycle(default_colors)
for key in reversed(markers_dict.keys()):
    if key == '0.5':
        colors_dict[key] = 'black'
    else:
        colors_dict[key] = next(color_iterator)


################################################################################################
###Build the figure#############################################################################
################################################################################################
# xlabel=r'\textrm{Na\textsubscript{V} Separation} $x$'
xlabel=r'Na\textsubscript{V} Separation $x$'
legend_title = r'Crossover Location: $\kappa$'



fig = plt.figure(figsize=(9,6))
ax = fig.add_subplot(1, 1, 1)
ax.tick_params(left=True, right=True)

def plotThresh(originalDF, propString, marker='o', identical_markers=False, alpha=1.0, legend=True, title=False, flip_legend_order=False, show=True, fix_crossOverPosition=False, scan_vals=None, filename=None, xlim=None, ylim=None):
    if fix_crossOverPosition:
        scan_vals = [0.5]
    elif scan_vals==None: 
        scan_vals = np.sort(np.unique(originalDF['CrossOverPosition'][originalDF['CrossOverPosition']>0.001].to_numpy()))
    scan_vals = list(scan_vals)
    scan_vals.sort()
    # scan_vals = reversed(scan_vals)


    ##############################################################################################################
    ###In some later simulations, we only simulated the x=0 setup once, since kappa has no effect in that case. 
    ###This saves simulation time. 
    ###The curves other than kappa=0.5 reuse this point. 
    ###In this ugly loop, we find that point.
    ###Note, some past simulations repeated the simulations for this identical point, but this loop had to be added
    ###after the code was changed to only simulate the x=0 threshold once. 
    yvalue_at__kappa_equals_zero_point_five__xequalszero_value = None
    for scan_value in scan_vals:
        DF = originalDF[(originalDF['CrossOverPosition']==scan_value) & (originalDF['Qthresh_'+propString]!=np.inf)]
        y = DF['Qthresh_'+propString].to_numpy()
        x = DF['deviation'].to_numpy()

        sort = np.argsort(x)
        sort = np.flip(sort)
        x = x[sort]
        y = y[sort]
        data = np.vstack((x, y))
        keep_looking_for_kappa0point5 = True
        threshold_point_at__x_equals_zero = None
        if len(data[0]):
            kappa = np.unique(DF['CrossOverPosition'].to_numpy())[0]
            if keep_looking_for_kappa0point5 == True:
                if kappa == 0.5:
                    keep_looking_for_kappa0point5 = False
                x0_index = np.argwhere(x==0).flatten()
                ###For axonal stimulation, the threshold may be undefined at x=0.0 in the Hu-based model
                ###so we exclude that case from this procedure. 
                if x0_index.size>0:
                    y_at_x0 = y[x0_index]
                    # print(y)
                    yvalue_at__kappa_equals_zero_point_five__xequalszero_value = y_at_x0[0]
                    threshold_point_at__x_equals_zero = [0.0, yvalue_at__kappa_equals_zero_point_five__xequalszero_value]
    print_args(threshold_point_at__x_equals_zero)
    ##############################################################################################################

    legend_varname_flag = False
    legend_varname = r'$ \kappa = '
    for scan_value in scan_vals:
        if legend_varname_flag == False:
            legend_varname = r'$'
        DF = originalDF[(originalDF['CrossOverPosition']==scan_value) & (originalDF['Qthresh_'+propString]!=np.inf)]

        y = DF['Qthresh_'+propString].to_numpy()
        x = DF['deviation'].to_numpy()

        sort = np.argsort(x)
        sort = np.flip(sort)
        x=x[sort]
        y = y[sort]
        data = np.vstack((x, y))
        print_args(data)
        
        print_args(data)
        # exit()
        print_args(scan_value)
        print_args(y)
        print_args(x)

        print_args(data[1])
        print_args(data[0])
        print_args(DF['CrossOverPosition'].to_numpy())
        print_args(np.unique(DF['CrossOverPosition'].to_numpy()))

        if type(xlim) is list or type(xlim) is tuple:
            plt.xlim(xlim) 
            xmin = xlim[0]
            select = np.argwhere(x>=xmin).flatten()
            print_args(select)
            x = x[select]
            y = y[select]

        print_args(x, y)

        if type(ylim) is list or type(ylim) is tuple:
            ymin, ymax = ylim
            # if ymin is not None:
            #     plt.ylim(ymin=ymin)
            if ymax is not None:
                select = np.argwhere(y<=ymax).flatten()
                x = x[select]
                y = y[select]
            print_args(x, y)

        if len(data[0]):
            kappa = np.unique(DF['CrossOverPosition'].to_numpy())[0]
            if identical_markers:
                marker=marker
            else:
                marker = markers_dict[str(kappa)]
                color = colors_dict[str(kappa)]
            if kappa == 0.5:
                """
                plot the data at kappa=0.5 in black
                """
                ax.plot(x, y, marker=marker, color='black', markersize=10, linewidth=1.0,  label=legend_varname+str(scan_value)+r'$', alpha=alpha)
            else:
                """Add the threshold point at x=0 if found...
                  This point is collected (as described above) from the kappa=0.5 data
                  and added to the other data sets, since more recent simulations
                  do not repeat the x=0.0 simulations for the other curves as it is
                  useless repetition: kappa has no effect at x=0."""  
                if threshold_point_at__x_equals_zero and 0 not in list(x):
                    x = np.array(list(x)+[threshold_point_at__x_equals_zero[0]])
                    y = np.array(list(y)+[threshold_point_at__x_equals_zero[1]])
                """
                plot the data where kappa is NOT 0.5 using a variety of coloured symbols. 
                """
                ax.plot(x, y, marker=marker, color=color, markersize=10, linewidth=1.0,  label=legend_varname+str(scan_value)+r'$', alpha=alpha)
            legend_varname_flag = False
            
        # plt.plot(data[0], data[1], marker='*',markersize=10, linewidth=1,  label=r'gammaAIS = '+str(scan_value))
        if propString=="backProp":
            threshtext='BP'
        else:
            threshtext='FP'
        plt.ylabel(r'$I_{'+threshtext+r'}'+r'\ \ [nA]$', fontsize=27)
        # plt.xlabel(r'$𝒳𝔁𝑥$', fontsize=70)
        # plt.xlabel(r'$𝑥$', fontsize=70)
        plt.xlabel(xlabel)


       


        if title:
            # plt.title('BackProp threshold: stimulating at the soma', fontsize=20)
            plt.title('stimulating at '+str(np.unique(originalDF['stimLoc'].to_numpy())), fontsize=20)
        plt.tight_layout()
    

    if type(ylim) is list or type(ylim) is tuple:
        ymin, ymax = ylim
        if ymin is not None:
            plt.ylim(ymin=ymin)
        if ymax is not None:
            plt.ylim(ymax=ymax) 


    if legend:
        if flip_legend_order:
            # Get the handles and labels
            handles, labels = ax.get_legend_handles_labels()
            # Reverse the order
            handles, labels = handles[::-1], labels[::-1]
            # Create the legend with the reversed handles and labels
            ax.legend(handles, labels, title=legend_title)
        else:
            ax.legend(title=legend_title)

    if filename:
        plt.savefig(filepath+filename+'.pdf')
    if show:
        plt.show()


###Figure 2 data
backprop_soma_apicalTipsBPcriterionMinimumVtipMinus63mV_newtune = pd.read_csv('backprop_soma_apicalTipsBPcriterionMinimumVtipMinus63mV_newtune.csv')
###Figure 3 data
backprop_soma_apicalTipsBPcriterionMinimumVtipMinus63mV_newtune_BACKWARDais = pd.read_csv('backprop_soma_apicalTipsBPcriterionMinimumVtipMinus63mV_newtune_BACKWARDais.csv')
###Figure 4 data
backprop_ais_apicalTipsBPcriterionMinimumVtipMinus63mV_newtune = pd.read_csv('backprop_ais_apicalTipsBPcriterionMinimumVtipMinus63mV_newtune.csv')
###Figure 5 data
forwardProp_AIS_crossOverPosition_newtune = pd.read_csv('forwardProp_AIS_crossOverPosition_newtune.csv')



###Plot Figure 2
# plotThresh(backprop_soma_apicalTipsBPcriterionMinimumVtipMinus63mV_newtune,'backProp', filename='backprop_soma_apicalTipsBPcriterionMinimumVtipMinus63mV_newtune', flip_legend_order=True, show=True, legend=True)

###Plot Figure 3
plotThresh(backprop_soma_apicalTipsBPcriterionMinimumVtipMinus63mV_newtune_BACKWARDais,'backProp', filename='backprop_soma_apicalTipsBPcriterionMinimumVtipMinus63mV_newtune_BACKWARDais', flip_legend_order=False, show=True, legend=True)

###Plot Figure 4
# plotThresh(backprop_ais_apicalTipsBPcriterionMinimumVtipMinus63mV_newtune,'backProp', filename='backprop_ais_apicalTipsBPcriterionMinimumVtipMinus63mV_newtune', flip_legend_order=False, show=True, legend=True, xlim=(0.37, None), ylim=(None, 10.0), scan_vals=[0.4, 0.5, 0.6, 0.7, 0.8])

###Plot Figure 5
# plotThresh(forwardProp_AIS_crossOverPosition_newtune,'forwardProp', filename='forwardProp_AIS_crossOverPosition_newtune', flip_legend_order=True, show=True, legend=True, ylim=(None, 0.541))



