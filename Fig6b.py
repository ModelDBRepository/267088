from Delta_Vrs_data import axonal_stimulation, somatic_stimulation
"""Select whether you wish to plot axonal stimulation or somatic stimulation data:
"""
data = axonal_stimulation
# data = somatic_stimulation

"""
Or import some data from a dictionary that was created using database_to_dictionary.py:
"""
# from data_dict_2024_02_09__23_22_26 import data


from sorcery import print_args
import numpy as np
import sqlite3
import pandas as pd
pd.set_option('display.max_colwidth', None)
from pylatex import NoEscape
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
xlim = ylim = (None, None); show_legend=True #<---default (don't touch)
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
rc('font', **{'family': 'STIXGeneral', 'serif': ['Times']})
############################################################################


# plotly_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
plotly_colors = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
plotly_symbols = ['o', 's', 'D', 'P', 'X', '^', 'v', '*']


fig = plt.figure(figsize=(6.8,7))
ax = fig.add_subplot(1, 1, 1)
ax.tick_params(left=True, right=True)
################################################################################################
###Build the figure#############################################################################
################################################################################################ = 'Qthresh'  ###'vpeak_SDtips', 'Qthresh'
ylabel = r'${I_{BP} \ \ [nA]}$'
xlabel = r'${ \Delta V_{RS}} \ \ [mV]$'


### Customize figure settings for axonal versus somatic stimulation data
if data == axonal_stimulation:
    xlim = (-2.0, 2.0)
    ylim = (-0.2, 10.0)
    show_legend=False
    filename=r'VRSaxonal_BackpropThresh'
    pass

elif data == somatic_stimulation:
    xlim = (-5.2, 5.2)
    ylim = (None, None)
    show_legend=True
    filename=r'VRSsomatic_BackpropThresh'
    pass

else:
    filename=r'Vrs_new'



###Settings specific to "Reference Curve" (the curve where all gating properties are right-shifted)
ref_curve_width=5.5/2
ref_curve_markerSize=27.0/2

###Width of other curves in the plot
standard_curve_width=0.63*ref_curve_width
standard_curve_markerSize=0.63*ref_curve_markerSize


###Count the plots
index=0
###Load, sort and plot the data
for j in range(len(data['stimLoc'])):

    stimLoc = data['stimLoc'][j]
    CrossOverPosition = data['CrossOverPosition'][j]
    activate_DeltaVrs_minf2__inAIS = data['activate_DeltaVrs_minf2__inAIS'][j] 
    activate_DeltaVrs_mtau2__inAIS = data['activate_DeltaVrs_mtau2__inAIS'][j] 
    activate_DeltaVrs_hinf2__inAIS = data['activate_DeltaVrs_hinf2__inAIS'][j] 
    activate_DeltaVrs_htau2__inAIS = data['activate_DeltaVrs_htau2__inAIS'][j] 
    print('-----------------------------------')
    print_args(CrossOverPosition)
    print_args(activate_DeltaVrs_minf2__inAIS)
    print_args(activate_DeltaVrs_mtau2__inAIS)
    print_args(activate_DeltaVrs_hinf2__inAIS)
    print_args(activate_DeltaVrs_htau2__inAIS)
    print('-----------------------------------')




    ###Get the data for this curve
    ydata = np.array(data['Qthresh'][j], dtype=float)
    vRS = np.array(data['DeltaVrs'][j], dtype=float)


    Nav_text = r'\text{Na\textsubscript{V}1.2}'
    m_inf = False
    m_tau = False
    h_inf = False
    h_tau = False

    title = 'Stimulation Site: '+str(stimLoc)
    addOn_text = r' \huge${'
    width=standard_curve_width
    marker_size=standard_curve_markerSize
    check_ref_curve=0
    if activate_DeltaVrs_minf2__inAIS:
        addOn_text+=  r' m^{^{'+Nav_text+r'}}_{\infty},'#\cdot'
        check_ref_curve+=1
        m_inf = True
    if activate_DeltaVrs_mtau2__inAIS:
        addOn_text+= r' \tau^{^{'+Nav_text+r'}}_{m},'#\cdot'
        check_ref_curve+=1
        m_tau = True
    if activate_DeltaVrs_hinf2__inAIS:
        addOn_text+= r' h^{^{'+Nav_text+r'}}_{\infty},'#\cdot'
        check_ref_curve+=1
        h_inf = True
    if activate_DeltaVrs_htau2__inAIS:
        addOn_text+= r' \tau^{^{'+Nav_text+r'}}_{h},'#\cdot'
        check_ref_curve+=1
        h_tau = True
    if addOn_text[-1] == r',':
        addOn_text = addOn_text[:-1]
    # if addOn_text[-5:] == r'\cdot':
    #     addOn_text = addOn_text[:-5]
    addOn_text+= r'\qquad }$'

    if check_ref_curve==4:
        width=ref_curve_width
        marker_size = ref_curve_markerSize
        addOn_text = 'Reference Curve'

    ###Select range of x-values to be plottted
    if (xlim[0] is not None) and (xlim[1] is not None):
        data_range = np.where( ((vRS >= xlim[0]) & (vRS <= xlim[1])))
        vRS_plot = vRS[data_range]
        ydata_plot = ydata[data_range]
    else:
        vRS_plot = vRS
        ydata_plot = ydata

    ###Plot the data
    plt.plot(vRS_plot, ydata_plot, color=plotly_colors[index], marker=plotly_symbols[index], markersize=marker_size, linewidth=width, label=addOn_text)#label=r'Separation: $x$')#, label=NoEscape(addOn_text))
    index+=1

plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.xlim(xlim)
plt.ylim(ylim)
if show_legend:
    plt.legend()
plt.savefig(filename+'.pdf')
plt.show()





