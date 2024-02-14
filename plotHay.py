"""This program is imported by Fig7a.py, Fig7b.py, Fig8a.py, and Fig8b.py, wherein it is used to produce
the plots from the manuscript having the same names as those programs."""

import pprint
from datetime import datetime
import matplotlib
import matplotlib.pyplot as plt
import tikzplotlib
from matplotlib import rc
import numpy as np
import itertools

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
from matplotlib.ticker import FuncFormatter
# Define a custom function for formatting
def format_ticks(x, _):
    return f"{x:g}"
     # return f"{x:.1f}"

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


somatic_thresh_data05 = {'electrode_location': ['L5PCtemplate[0].soma[0](0.5)',
'L5PCtemplate[0].soma[0](0.5)',
'L5PCtemplate[0].soma[0](0.5)'],
'separation': [0.01, 0.5, 1.0],
'crossover': [0.5, 0.5, 0.5],
'Ithresh': [4.160763959350585, 4.39141350830078, 4.727652821655273],
'vpeak': [40.617710197552775, 40.505886505892306, 40.26155485962343]}

axonal_thresh_data05 = {'electrode_location': ['L5PCtemplate[0].axon[1](0.983871)',
'L5PCtemplate[0].axon[1](0.983871)',
'L5PCtemplate[0].axon[1](0.983871)'],
'separation': [0.01, 0.5, 1.0],
'crossover': [0.5, 0.5, 0.5],
'Ithresh': [0.6352254846191403, 0.6031213134765622, 0.5802333587646481],
'vpeak': [42.88338802147963, 42.888373722049415, 42.892018771266955]}

class DisplayHayData():
    """docstring for DislayHayData"""
    def __init__(self, data_dict=None, xAxisVariableName=None, filename=None, legend_location=None, linewidth=0.0):
        self.legend_location = legend_location
        self.symbol_scaleNavTwo = r'$\alpha\textsubscript{\textsubscript{Na\textsubscript{V}1.2}}$'
        if filename is None:
            self.filename = ''
        else:
            self.filename = filename
        if data_dict is not None:
            self.data_dict = data_dict
        else:
            self.data_dict = somatic_thresh_data05

        if 'scaleNavTwo' in self.data_dict.keys():
            self.scaleNavTwo = np.array(self.data_dict['scaleNavTwo'])
        if 'DeltaVrs' in self.data_dict.keys():
            self.DeltaVrs = np.array(self.data_dict['DeltaVrs'])
        self.separation = np.array(self.data_dict['separation'])
        self.crossover = np.array(self.data_dict['crossover'])
        self.Ithresh = np.array(self.data_dict['Ithresh'])
        self.electrode_location = self.data_dict['electrode_location'][0]
        if 'axon' in self.electrode_location:
            self.electrode_location = 'axon'
            self.filename+='HayAxonal'
        elif 'soma' in self.electrode_location:
            self.electrode_location = 'soma'
            self.filename += 'HaySomatic'

        if xAxisVariableName is not None:
            self.xAxisVariableName = xAxisVariableName
            self.xAxisVariable = np.array(self.data_dict[self.xAxisVariableName])
            if self.xAxisVariableName == 'DeltaVrs':
                self.xAxisLabel = r'$\Delta V_{RS}$'
            elif self.xAxisVariableName == 'scaleNavTwo':
                self.xAxisLabel = r'Na\textsubscript{V}1.2 density scaling factor:\ '+self.symbol_scaleNavTwo
        else:
            self.xAxisVariableName = 'separation'
            self.xAxisVariable = np.array(self.data_dict[self.xAxisVariableName])
            # self.xAxisLabel = r'$𝒙$'
            self.xAxisLabel = r'Na\textsubscript{V} Separation $x$'
            # self.xAxisLabel = r'\NaV Separation\ $x$'
        self.filename += '__'+self.xAxisVariableName

        self.alpha=0.7
        self.linewidth = linewidth
        if self.linewidth==0.0: 
            self.linestyle=''
        else:
            self.linestyle ='--'
        self.marker = r'o'
        self.markersize = 10.0
        self.color = None
        self.previous_color = None
        self.plotted_kappa_valsandcolors_dict = {}

    def replace_addplot_options(self, text):
        old_str = r'\addplot [line width=0pt,'
        new_str = r'\addplot [only marks, line width=0pt,'
        return text.replace(old_str, new_str)

    def plot_thresh(self, default_color=None, default_linewidth=None, title=False, select_kappa_vals=None, legend_location=None, xlim=(None, None), ylim=(None, None), vline=None, use_ratio_of_total_AIS_NaV_conductances_from_each_subtype=False):
        if default_color is not None:
            self.color = default_color
        if legend_location is not None:
            self.legend_location = legend_location

        ###Create and Close Figure 1 to clear any pre-existing figures: 
        self.fig1, self.ax1 = plt.subplots(figsize=(9,6))
        plt.close(self.fig1)

        ### Create Figure 1 and Axes 1
        self.fig1, self.ax1 = plt.subplots(figsize=(9, 6))
        # self.fig1 = plt.figure(figsize=(9,6))
        # self.ax1 = self.fig1.add_subplot(1, 1, 1)
        self.ax1.tick_params(left=True, right=True)
        # Apply the custom format function to the x and y axes
        self.formatter = FuncFormatter(format_ticks)
        self.ax1.xaxis.set_major_formatter(self.formatter)
        self.ax1.yaxis.set_major_formatter(self.formatter)
        

        for i, threshold in enumerate(self.Ithresh):
            if threshold is not None:
                if select_kappa_vals is None or self.crossover[i] in select_kappa_vals:
                    self.color = colors_dict[str(self.crossover[i])]
                    if str(self.crossover[i]) in self.plotted_kappa_valsandcolors_dict.keys():
                        self.color = self.plotted_kappa_valsandcolors_dict[str(self.crossover[i])]

                    if use_ratio_of_total_AIS_NaV_conductances_from_each_subtype==True:
                        self.marker = f"${self.xAxisVariable[i]}$"
                        self.line, = self.ax1.plot(self.lookup_conductances(self.xAxisVariable[i], self.crossover[i]), threshold, linewidth=self.linewidth, marker=self.marker, markersize=self.markersize, color=self.color, alpha=self.alpha)
                    else:
                        self.line, = self.ax1.plot(self.xAxisVariable[i], threshold, linewidth=self.linewidth, marker=self.marker, markersize=self.markersize, color=self.color, alpha=self.alpha)
                    

                    if str(self.crossover[i]) not in self.plotted_kappa_valsandcolors_dict.keys():
                        self.previous_color = self.ax1.lines[-1].get_color()
                        self.plotted_kappa_valsandcolors_dict[str(self.crossover[i])] = self.previous_color
                        self.line.set_label(r'$\kappa = '+f"{self.crossover[i]:g}"+r'$')
                    # if self.linewidth !=0:
                    #     self.line.set_linewidth(self.linewidth)
                    #     self.line.set_linestyle(self.linestyle)
                


                    # self.ax1.lines[-1].set_color(self.previous_color)
            
        print(self.plotted_kappa_valsandcolors_dict.items())
        # indicate stimulation site in title
        if title==True:
            self.ax1.set_title(r'BackProp threshold as a function of $Na_{V}$ distribution | Electrode: '+self.electrode_location)
        # axis labels
        self.ax1.set_ylabel(r'$I_{BP}\ [nA]$', fontsize=27)
        # self.ax1.set_ylabel(r'$\Qbp\ [\nano\Ampere]$', fontsize=27)

        if use_ratio_of_total_AIS_NaV_conductances_from_each_subtype==True:
            # self.ax1.set_xlabel(r'$\bar{G}_{\text{Na\textsubscript{V}1.2}}\fracslash \bar{G}_{\text{Na\textsubscript{V}1.6}}$', fontsize=23)
            # self.ax1.set_xlabel(r'$\bar{G}_{\text{Na\textsubscript{V}1.2}}\fracslash \bar{G}_{\text{Na\textsubscript{V}, total}}$', fontsize=23)
            # self.ax1.set_xlabel(r'$\bar{G}_{\text{Na\textsubscript{V}1.6}}\fracslash \bar{G}_{\text{Na\textsubscript{V}, total}}$', fontsize=23)
            self.ax1.set_xlabel(r'$\bar{G}_{\text{Na\textsubscript{V}1.2}} \ \fracslash\ \left(\bar{G}_{\text{Na\textsubscript{V}1.2}} + \bar{G}_{\text{Na\textsubscript{V}1.6}}\right)$', fontsize=17)
        else:
            self.ax1.set_xlabel(self.xAxisLabel, fontsize=27)

        ###Plot a vertical line at x=vline
        if vline is not None:
            plt.axvline(vline, linestyle=':', color='black', alpha=1.0, label=self.symbol_scaleNavTwo+f"$ = {vline}"+r'$')

        ### Display the legend with (reversed order)
        handles, labels = self.ax1.get_legend_handles_labels()
        if self.legend_location is None:
            self.ax1.legend(reversed(handles), reversed(labels))
        else:
            self.ax1.legend(reversed(handles), reversed(labels), loc=legend_location)

        
        self.ax1.set_xlim(xlim)
        self.ax1.set_ylim(ylim)

        # Adjust layout
        if self.xAxisVariableName=='scaleNavTwo':
            self.fig1.subplots_adjust(top=0.98, bottom=0.15)
        else:
            self.fig1.subplots_adjust(top=0.98, bottom=0.13)
        ###Save figure

        # Get the current figure size
        # fig = plt.gcf()
        fig_width, fig_height = self.fig1.get_size_inches()

        # # Get the TikZ code for the current figure
        # tikz_code = tikzplotlib.get_tikz_code(
        #     axis_width=f"{fig_width*(7.0/9.0)}in",
        #     axis_height=f"{fig_height*(7.0/9.0)}in",
        # )
        # if self.linewidth==0:
        #     tikz_code = self.replace_addplot_options(tikz_code) #<---removes lines from legend handles in the TiKz plots

        # # Save the TikZ code to a .tex file
        # with open(self.filename+".tex", "w") as f:
        #     f.write(tikz_code)
        self.fig1.savefig(self.filename+'.pdf')
        
        # Show Figure 1
        # self.fig1.show()
        self.fig1.show()
        # plt.show(block=False)
        
        # ###Close Figure 1
        # plt.waitforbuttonpress()
        # plt.close(self.fig1)


    def lookup_conductances(self, x, kappa):
        """This code uses a lookup table (called conductances_dict) to 
            map "deviation" and "CrossOverPosition" ("x" and "kappa", respectively) settings
            to the ratio ("ratioNaVTwoPerNaVSix"), of total NaVTwo conductance in the AIS to 
            the corresponding total NaVSix conductance added up segment by segment.
        """ 
            
        conductances_dict = {
        'x=0.0|𝜅=0.2|GbarNaVSix': 72.11545740387473,
        'x=0.0|𝜅=0.2|GbarNaVTwo': 72.11545740387473,
        'x=0.0|𝜅=0.2|ratioNaVSixPerNaVTwo': 1.0,
        'x=0.0|𝜅=0.2|ratioNaVSixPerTotalNaV': 0.5,
        'x=0.0|𝜅=0.2|ratioNaVTwoPerNaVSix': 1.0,
        'x=0.0|𝜅=0.2|ratioNaVTwoPerTotalNaV': 0.5,
        'x=0.0|𝜅=0.3|GbarNaVSix': 72.11545740387473,
        'x=0.0|𝜅=0.3|GbarNaVTwo': 72.11545740387473,
        'x=0.0|𝜅=0.3|ratioNaVSixPerNaVTwo': 1.0,
        'x=0.0|𝜅=0.3|ratioNaVSixPerTotalNaV': 0.5,
        'x=0.0|𝜅=0.3|ratioNaVTwoPerNaVSix': 1.0,
        'x=0.0|𝜅=0.3|ratioNaVTwoPerTotalNaV': 0.5,
        'x=0.0|𝜅=0.4|GbarNaVSix': 72.11545740387473,
        'x=0.0|𝜅=0.4|GbarNaVTwo': 72.11545740387473,
        'x=0.0|𝜅=0.4|ratioNaVSixPerNaVTwo': 1.0,
        'x=0.0|𝜅=0.4|ratioNaVSixPerTotalNaV': 0.5,
        'x=0.0|𝜅=0.4|ratioNaVTwoPerNaVSix': 1.0,
        'x=0.0|𝜅=0.4|ratioNaVTwoPerTotalNaV': 0.5,
        'x=0.0|𝜅=0.5|GbarNaVSix': 72.11545740387473,
        'x=0.0|𝜅=0.5|GbarNaVTwo': 72.11545740387473,
        'x=0.0|𝜅=0.5|ratioNaVSixPerNaVTwo': 1.0,
        'x=0.0|𝜅=0.5|ratioNaVSixPerTotalNaV': 0.5,
        'x=0.0|𝜅=0.5|ratioNaVTwoPerNaVSix': 1.0,
        'x=0.0|𝜅=0.5|ratioNaVTwoPerTotalNaV': 0.5,
        'x=0.0|𝜅=0.6|GbarNaVSix': 72.11545740387473,
        'x=0.0|𝜅=0.6|GbarNaVTwo': 72.11545740387473,
        'x=0.0|𝜅=0.6|ratioNaVSixPerNaVTwo': 1.0,
        'x=0.0|𝜅=0.6|ratioNaVSixPerTotalNaV': 0.5,
        'x=0.0|𝜅=0.6|ratioNaVTwoPerNaVSix': 1.0,
        'x=0.0|𝜅=0.6|ratioNaVTwoPerTotalNaV': 0.5,
        'x=0.0|𝜅=0.7|GbarNaVSix': 72.11545740387473,
        'x=0.0|𝜅=0.7|GbarNaVTwo': 72.11545740387473,
        'x=0.0|𝜅=0.7|ratioNaVSixPerNaVTwo': 1.0,
        'x=0.0|𝜅=0.7|ratioNaVSixPerTotalNaV': 0.5,
        'x=0.0|𝜅=0.7|ratioNaVTwoPerNaVSix': 1.0,
        'x=0.0|𝜅=0.7|ratioNaVTwoPerTotalNaV': 0.5,
        'x=0.0|𝜅=0.8|GbarNaVSix': 72.11545740387473,
        'x=0.0|𝜅=0.8|GbarNaVTwo': 72.11545740387473,
        'x=0.0|𝜅=0.8|ratioNaVSixPerNaVTwo': 1.0,
        'x=0.0|𝜅=0.8|ratioNaVSixPerTotalNaV': 0.5,
        'x=0.0|𝜅=0.8|ratioNaVTwoPerNaVSix': 1.0,
        'x=0.0|𝜅=0.8|ratioNaVTwoPerTotalNaV': 0.5,
        'x=0.1|𝜅=0.2|GbarNaVSix': 76.42930037171149,
        'x=0.1|𝜅=0.2|GbarNaVTwo': 67.80161443603785,
        'x=0.1|𝜅=0.2|ratioNaVSixPerNaVTwo': 1.1272489749313088,
        'x=0.1|𝜅=0.2|ratioNaVSixPerTotalNaV': 0.5299092810552224,
        'x=0.1|𝜅=0.2|ratioNaVTwoPerNaVSix': 0.8871154662712708,
        'x=0.1|𝜅=0.2|ratioNaVTwoPerTotalNaV': 0.47009071894477755,
        'x=0.1|𝜅=0.3|GbarNaVSix': 74.99829155717508,
        'x=0.1|𝜅=0.3|GbarNaVTwo': 69.23262325057425,
        'x=0.1|𝜅=0.3|ratioNaVSixPerNaVTwo': 1.083279645287065,
        'x=0.1|𝜅=0.3|ratioNaVSixPerTotalNaV': 0.5199876299553604,
        'x=0.1|𝜅=0.3|ratioNaVTwoPerNaVSix': 0.9231226713717157,
        'x=0.1|𝜅=0.3|ratioNaVTwoPerTotalNaV': 0.4800123700446395,
        'x=0.1|𝜅=0.4|GbarNaVSix': 73.55752920255324,
        'x=0.1|𝜅=0.4|GbarNaVTwo': 70.67338560519605,
        'x=0.1|𝜅=0.4|ratioNaVSixPerNaVTwo': 1.0408094726559294,
        'x=0.1|𝜅=0.4|ratioNaVSixPerTotalNaV': 0.509998354379162,
        'x=0.1|𝜅=0.4|ratioNaVTwoPerNaVSix': 0.9607906406234064,
        'x=0.1|𝜅=0.4|ratioNaVTwoPerTotalNaV': 0.49000164562083803,
        'x=0.1|𝜅=0.5|GbarNaVSix': 72.11545754101068,
        'x=0.1|𝜅=0.5|GbarNaVTwo': 72.1154572667387,
        'x=0.1|𝜅=0.5|ratioNaVSixPerNaVTwo': 1.0000000038032344,
        'x=0.1|𝜅=0.5|ratioNaVSixPerTotalNaV': 0.5000000009508085,
        'x=0.1|𝜅=0.5|ratioNaVTwoPerNaVSix': 0.9999999961967657,
        'x=0.1|𝜅=0.5|ratioNaVTwoPerTotalNaV': 0.49999999904919135,
        'x=0.1|𝜅=0.6|GbarNaVSix': 70.67338663646102,
        'x=0.1|𝜅=0.6|GbarNaVTwo': 73.55752817128831,
        'x=0.1|𝜅=0.6|ratioNaVSixPerNaVTwo': 0.9607906681133821,
        'x=0.1|𝜅=0.6|ratioNaVSixPerTotalNaV': 0.49000165277093444,
        'x=0.1|𝜅=0.6|ratioNaVTwoPerNaVSix': 1.0408094428764694,
        'x=0.1|𝜅=0.6|ratioNaVTwoPerTotalNaV': 0.5099983472290655,
        'x=0.1|𝜅=0.7|GbarNaVSix': 69.2326307033145,
        'x=0.1|𝜅=0.7|GbarNaVTwo': 74.99828410443483,
        'x=0.1|𝜅=0.7|ratioNaVSixPerNaVTwo': 0.9231228624765377,
        'x=0.1|𝜅=0.7|ratioNaVSixPerTotalNaV': 0.4800124217169197,
        'x=0.1|𝜅=0.7|ratioNaVTwoPerNaVSix': 1.0832794210265984,
        'x=0.1|𝜅=0.7|ratioNaVTwoPerTotalNaV': 0.5199875782830802,
        'x=0.1|𝜅=0.8|GbarNaVSix': 67.80166774928324,
        'x=0.1|𝜅=0.8|GbarNaVTwo': 76.42924705846617,
        'x=0.1|𝜅=0.8|ratioNaVSixPerNaVTwo': 0.8871167826292586,
        'x=0.1|𝜅=0.8|ratioNaVSixPerTotalNaV': 0.47009108858290566,
        'x=0.1|𝜅=0.8|ratioNaVTwoPerNaVSix': 1.1272473022505283,
        'x=0.1|𝜅=0.8|ratioNaVTwoPerTotalNaV': 0.5299089114170945,
        'x=0.2|𝜅=0.2|GbarNaVSix': 80.7431433395483,
        'x=0.2|𝜅=0.2|GbarNaVTwo': 63.48777146820104,
        'x=0.2|𝜅=0.2|ratioNaVSixPerNaVTwo': 1.271790479840514,
        'x=0.2|𝜅=0.2|ratioNaVSixPerTotalNaV': 0.5598185621104448,
        'x=0.2|𝜅=0.2|ratioNaVTwoPerNaVSix': 0.7862930379266581,
        'x=0.2|𝜅=0.2|ratioNaVTwoPerTotalNaV': 0.44018143788955516,
        'x=0.2|𝜅=0.3|GbarNaVSix': 77.88112571047546,
        'x=0.2|𝜅=0.3|GbarNaVTwo': 66.34978909727387,
        'x=0.2|𝜅=0.3|ratioNaVSixPerNaVTwo': 1.173796130629681,
        'x=0.2|𝜅=0.3|ratioNaVSixPerTotalNaV': 0.5399752599107206,
        'x=0.2|𝜅=0.3|ratioNaVTwoPerNaVSix': 0.8519366982949175,
        'x=0.2|𝜅=0.3|ratioNaVTwoPerTotalNaV': 0.46002474008927924,
        'x=0.2|𝜅=0.4|GbarNaVSix': 74.99960100123192,
        'x=0.2|𝜅=0.4|GbarNaVTwo': 69.2313138065174,
        'x=0.2|𝜅=0.4|ratioNaVSixPerNaVTwo': 1.0833190485281747,
        'x=0.2|𝜅=0.4|ratioNaVSixPerTotalNaV': 0.5199967087583244,
        'x=0.2|𝜅=0.4|ratioNaVTwoPerNaVSix': 0.9230890949057213,
        'x=0.2|𝜅=0.4|ratioNaVTwoPerTotalNaV': 0.48000329124167557,
        'x=0.2|𝜅=0.5|GbarNaVSix': 72.1154576781466,
        'x=0.2|𝜅=0.5|GbarNaVTwo': 72.11545712960273,
        'x=0.2|𝜅=0.5|ratioNaVSixPerNaVTwo': 1.0000000076064672,
        'x=0.2|𝜅=0.5|ratioNaVSixPerTotalNaV': 0.5000000019016168,
        'x=0.2|𝜅=0.5|ratioNaVTwoPerNaVSix': 0.9999999923935328,
        'x=0.2|𝜅=0.5|ratioNaVTwoPerTotalNaV': 0.49999999809838314,
        'x=0.2|𝜅=0.6|GbarNaVSix': 69.23131586904739,
        'x=0.2|𝜅=0.6|GbarNaVTwo': 74.99959893870191,
        'x=0.2|𝜅=0.6|ratioNaVSixPerNaVTwo': 0.9230891477917235,
        'x=0.2|𝜅=0.6|ratioNaVSixPerTotalNaV': 0.48000330554186926,
        'x=0.2|𝜅=0.6|ratioNaVTwoPerNaVSix': 1.0833189864622155,
        'x=0.2|𝜅=0.6|ratioNaVTwoPerTotalNaV': 0.5199966944581309,
        'x=0.2|𝜅=0.7|GbarNaVSix': 66.34980400275434,
        'x=0.2|𝜅=0.7|GbarNaVTwo': 77.881110804995,
        'x=0.2|𝜅=0.7|ratioNaVSixPerNaVTwo': 0.8519370527326751,
        'x=0.2|𝜅=0.7|ratioNaVSixPerTotalNaV': 0.46002484343383954,
        'x=0.2|𝜅=0.7|ratioNaVTwoPerNaVSix': 1.1737956422864786,
        'x=0.2|𝜅=0.7|ratioNaVTwoPerTotalNaV': 0.5399751565661605,
        'x=0.2|𝜅=0.8|GbarNaVSix': 63.4878780946918,
        'x=0.2|𝜅=0.8|GbarNaVTwo': 80.74303671305756,
        'x=0.2|𝜅=0.8|ratioNaVSixPerNaVTwo': 0.7862953968441058,
        'x=0.2|𝜅=0.8|ratioNaVSixPerTotalNaV': 0.44018217716581154,
        'x=0.2|𝜅=0.8|ratioNaVTwoPerNaVSix': 1.2717866644185176,
        'x=0.2|𝜅=0.8|ratioNaVTwoPerTotalNaV': 0.5598178228341885,
        'x=0.3|𝜅=0.2|GbarNaVSix': 85.05698630738512,
        'x=0.3|𝜅=0.2|GbarNaVTwo': 59.17392850036422,
        'x=0.3|𝜅=0.2|ratioNaVSixPerNaVTwo': 1.4374064467743017,
        'x=0.3|𝜅=0.2|ratioNaVSixPerTotalNaV': 0.5897278431656673,
        'x=0.3|𝜅=0.2|ratioNaVTwoPerNaVSix': 0.6956974502543175,
        'x=0.3|𝜅=0.2|ratioNaVTwoPerTotalNaV': 0.4102721568343327,
        'x=0.3|𝜅=0.3|GbarNaVSix': 80.76395986377587,
        'x=0.3|𝜅=0.3|GbarNaVTwo': 63.466954943973484,
        'x=0.3|𝜅=0.3|ratioNaVSixPerNaVTwo': 1.2725356043167915,
        'x=0.3|𝜅=0.3|ratioNaVSixPerTotalNaV': 0.5599628898660811,
        'x=0.3|𝜅=0.3|ratioNaVTwoPerNaVSix': 0.7858326294429205,
        'x=0.3|𝜅=0.3|ratioNaVTwoPerTotalNaV': 0.440037110133919,
        'x=0.3|𝜅=0.4|GbarNaVSix': 76.44167279991062,
        'x=0.3|𝜅=0.4|GbarNaVTwo': 67.78924200783875,
        'x=0.3|𝜅=0.4|ratioNaVSixPerNaVTwo': 1.127637225845826,
        'x=0.3|𝜅=0.4|ratioNaVSixPerTotalNaV': 0.5299950631374869,
        'x=0.3|𝜅=0.4|ratioNaVTwoPerNaVSix': 0.886810028154146,
        'x=0.3|𝜅=0.4|ratioNaVTwoPerTotalNaV': 0.4700049368625131,
        'x=0.3|𝜅=0.5|GbarNaVSix': 72.11545781528258,
        'x=0.3|𝜅=0.5|GbarNaVTwo': 72.11545699246676,
        'x=0.3|𝜅=0.5|ratioNaVSixPerNaVTwo': 1.0000000114097012,
        'x=0.3|𝜅=0.5|ratioNaVSixPerTotalNaV': 0.5000000028524253,
        'x=0.3|𝜅=0.5|ratioNaVTwoPerNaVSix': 0.9999999885902988,
        'x=0.3|𝜅=0.5|ratioNaVTwoPerTotalNaV': 0.4999999971475747,
        'x=0.3|𝜅=0.6|GbarNaVSix': 67.78924510163374,
        'x=0.3|𝜅=0.6|GbarNaVTwo': 76.44166970611555,
        'x=0.3|𝜅=0.6|ratioNaVSixPerNaVTwo': 0.8868101045183006,
        'x=0.3|𝜅=0.6|ratioNaVSixPerTotalNaV': 0.47000495831280364,
        'x=0.3|𝜅=0.6|ratioNaVTwoPerNaVSix': 1.1276371287437936,
        'x=0.3|𝜅=0.6|ratioNaVTwoPerTotalNaV': 0.5299950416871964,
        'x=0.3|𝜅=0.7|GbarNaVSix': 63.466977302194174,
        'x=0.3|𝜅=0.7|GbarNaVTwo': 80.76393750555515,
        'x=0.3|𝜅=0.7|ratioNaVSixPerNaVTwo': 0.7858331238224827,
        'x=0.3|𝜅=0.7|ratioNaVSixPerTotalNaV': 0.44003726515075936,
        'x=0.3|𝜅=0.7|ratioNaVTwoPerNaVSix': 1.2725348037453013,
        'x=0.3|𝜅=0.7|ratioNaVTwoPerTotalNaV': 0.5599627348492406,
        'x=0.3|𝜅=0.8|GbarNaVSix': 59.17408844010038,
        'x=0.3|𝜅=0.8|GbarNaVTwo': 85.05682636764901,
        'x=0.3|𝜅=0.8|ratioNaVSixPerNaVTwo': 0.6957006388214713,
        'x=0.3|𝜅=0.8|ratioNaVSixPerTotalNaV': 0.4102732657487173,
        'x=0.3|𝜅=0.8|ratioNaVTwoPerNaVSix': 1.4373998587869878,
        'x=0.3|𝜅=0.8|ratioNaVTwoPerTotalNaV': 0.5897267342512826,
        'x=0.4|𝜅=0.2|GbarNaVSix': 89.37082927522198,
        'x=0.4|𝜅=0.2|GbarNaVTwo': 54.86008553252737,
        'x=0.4|𝜅=0.2|ratioNaVSixPerNaVTwo': 1.6290683546643883,
        'x=0.4|𝜅=0.2|ratioNaVSixPerTotalNaV': 0.61963712422089,
        'x=0.4|𝜅=0.2|ratioNaVTwoPerNaVSix': 0.613847784309833,
        'x=0.4|𝜅=0.2|ratioNaVTwoPerTotalNaV': 0.3803628757791101,
        'x=0.4|𝜅=0.3|GbarNaVSix': 83.64679401707627,
        'x=0.4|𝜅=0.3|GbarNaVTwo': 60.58412079067312,
        'x=0.4|𝜅=0.3|ratioNaVSixPerNaVTwo': 1.380671914115714,
        'x=0.4|𝜅=0.3|ratioNaVSixPerTotalNaV': 0.5799505198214413,
        'x=0.4|𝜅=0.3|ratioNaVTwoPerNaVSix': 0.7242850309158894,
        'x=0.4|𝜅=0.3|ratioNaVTwoPerTotalNaV': 0.4200494801785587,
        'x=0.4|𝜅=0.4|GbarNaVSix': 77.88374459858925,
        'x=0.4|𝜅=0.4|GbarNaVTwo': 66.34717020916013,
        'x=0.4|𝜅=0.4|ratioNaVSixPerNaVTwo': 1.1738819357790236,
        'x=0.4|𝜅=0.4|ratioNaVSixPerTotalNaV': 0.5399934175166491,
        'x=0.4|𝜅=0.4|ratioNaVTwoPerNaVSix': 0.8518744258010664,
        'x=0.4|𝜅=0.4|ratioNaVTwoPerTotalNaV': 0.46000658248335097,
        'x=0.4|𝜅=0.5|GbarNaVSix': 72.11545795241855,
        'x=0.4|𝜅=0.5|GbarNaVTwo': 72.11545685533083,
        'x=0.4|𝜅=0.5|ratioNaVSixPerNaVTwo': 1.0000000152129345,
        'x=0.4|𝜅=0.5|ratioNaVSixPerTotalNaV': 0.5000000038032335,
        'x=0.4|𝜅=0.5|ratioNaVTwoPerNaVSix': 0.9999999847870658,
        'x=0.4|𝜅=0.5|ratioNaVTwoPerTotalNaV': 0.4999999961967664,
        'x=0.4|𝜅=0.6|GbarNaVSix': 66.34717433422014,
        'x=0.4|𝜅=0.6|GbarNaVTwo': 77.88374047352924,
        'x=0.4|𝜅=0.6|ratioNaVSixPerNaVTwo': 0.8518745238843518,
        'x=0.4|𝜅=0.6|ratioNaVSixPerTotalNaV': 0.4600066110837381,
        'x=0.4|𝜅=0.6|ratioNaVTwoPerNaVSix': 1.1738818006203897,
        'x=0.4|𝜅=0.6|ratioNaVTwoPerTotalNaV': 0.5399933889162618,
        'x=0.4|𝜅=0.7|GbarNaVSix': 60.584150601634015,
        'x=0.4|𝜅=0.7|GbarNaVTwo': 83.64676420611532,
        'x=0.4|𝜅=0.7|ratioNaVSixPerNaVTwo': 0.7242856454356996,
        'x=0.4|𝜅=0.7|ratioNaVSixPerTotalNaV': 0.42004968686767913,
        'x=0.4|𝜅=0.7|ratioNaVTwoPerNaVSix': 1.3806707426852871,
        'x=0.4|𝜅=0.7|ratioNaVTwoPerTotalNaV': 0.5799503131323208,
        'x=0.4|𝜅=0.8|GbarNaVSix': 54.86029878550896,
        'x=0.4|𝜅=0.8|GbarNaVTwo': 89.37061602224041,
        'x=0.4|𝜅=0.8|ratioNaVSixPerNaVTwo': 0.6138516352159489,
        'x=0.4|𝜅=0.8|ratioNaVSixPerTotalNaV': 0.3803643543316233,
        'x=0.4|𝜅=0.8|ratioNaVTwoPerNaVSix': 1.6290581349485316,
        'x=0.4|𝜅=0.8|ratioNaVTwoPerTotalNaV': 0.6196356456683767,
        'x=0.5|𝜅=0.2|GbarNaVSix': 93.6846722430587,
        'x=0.5|𝜅=0.2|GbarNaVTwo': 50.54624256469062,
        'x=0.5|𝜅=0.2|ratioNaVSixPerNaVTwo': 1.8534448356504085,
        'x=0.5|𝜅=0.2|ratioNaVSixPerTotalNaV': 0.6495464052761118,
        'x=0.5|𝜅=0.2|ratioNaVTwoPerNaVSix': 0.539535885161147,
        'x=0.5|𝜅=0.2|ratioNaVTwoPerTotalNaV': 0.3504535947238882,
        'x=0.5|𝜅=0.3|GbarNaVSix': 86.52962817037667,
        'x=0.5|𝜅=0.3|GbarNaVTwo': 57.701286637372654,
        'x=0.5|𝜅=0.3|ratioNaVSixPerNaVTwo': 1.4996134958684288,
        'x=0.5|𝜅=0.3|ratioNaVSixPerTotalNaV': 0.5999381497768019,
        'x=0.5|𝜅=0.3|ratioNaVTwoPerNaVSix': 0.66683849055446,
        'x=0.5|𝜅=0.3|ratioNaVTwoPerTotalNaV': 0.4000618502231981,
        'x=0.5|𝜅=0.4|GbarNaVSix': 79.3258163972679,
        'x=0.5|𝜅=0.4|GbarNaVTwo': 64.9050984104815,
        'x=0.5|𝜅=0.4|ratioNaVSixPerNaVTwo': 1.2221815903518853,
        'x=0.5|𝜅=0.4|ratioNaVSixPerTotalNaV': 0.5499917718958113,
        'x=0.5|𝜅=0.4|ratioNaVTwoPerNaVSix': 0.818209018933172,
        'x=0.5|𝜅=0.4|ratioNaVTwoPerTotalNaV': 0.4500082281041887,
        'x=0.5|𝜅=0.5|GbarNaVSix': 72.11545808955454,
        'x=0.5|𝜅=0.5|GbarNaVTwo': 72.11545671819492,
        'x=0.5|𝜅=0.5|ratioNaVSixPerNaVTwo': 1.0000000190161678,
        'x=0.5|𝜅=0.5|ratioNaVSixPerTotalNaV': 0.5000000047540418,
        'x=0.5|𝜅=0.5|ratioNaVTwoPerNaVSix': 0.9999999809838327,
        'x=0.5|𝜅=0.5|ratioNaVTwoPerTotalNaV': 0.4999999952459581,
        'x=0.5|𝜅=0.6|GbarNaVSix': 64.90510356680645,
        'x=0.5|𝜅=0.6|GbarNaVTwo': 79.32581124094294,
        'x=0.5|𝜅=0.6|ratioNaVSixPerNaVTwo': 0.8182091371201327,
        'x=0.5|𝜅=0.6|ratioNaVSixPerTotalNaV': 0.45000826385467224,
        'x=0.5|𝜅=0.6|ratioNaVTwoPerNaVSix': 1.2221814138127574,
        'x=0.5|𝜅=0.6|ratioNaVTwoPerTotalNaV': 0.5499917361453277,
        'x=0.5|𝜅=0.7|GbarNaVSix': 57.701323901073806,
        'x=0.5|𝜅=0.7|GbarNaVTwo': 86.52959090667545,
        'x=0.5|𝜅=0.7|ratioNaVSixPerNaVTwo': 0.6668392083733098,
        'x=0.5|𝜅=0.7|ratioNaVSixPerTotalNaV': 0.40006210858459884,
        'x=0.5|𝜅=0.7|ratioNaVTwoPerNaVSix': 1.499611881609967,
        'x=0.5|𝜅=0.7|ratioNaVTwoPerTotalNaV': 0.5999378914154012,
        'x=0.5|𝜅=0.8|GbarNaVSix': 50.546509130917535,
        'x=0.5|𝜅=0.8|GbarNaVTwo': 93.68440567683173,
        'x=0.5|𝜅=0.8|ratioNaVSixPerNaVTwo': 0.5395402657009943,
        'x=0.5|𝜅=0.8|ratioNaVSixPerTotalNaV': 0.35045544291452946,
        'x=0.5|𝜅=0.8|ratioNaVTwoPerNaVSix': 1.853429787489088,
        'x=0.5|𝜅=0.8|ratioNaVTwoPerTotalNaV': 0.6495445570854705,
        'x=0.6|𝜅=0.2|GbarNaVSix': 97.99851521089559,
        'x=0.6|𝜅=0.2|GbarNaVTwo': 46.23239959685379,
        'x=0.6|𝜅=0.2|ratioNaVSixPerNaVTwo': 2.1196934631436393,
        'x=0.6|𝜅=0.2|ratioNaVSixPerTotalNaV': 0.6794556863313346,
        'x=0.6|𝜅=0.2|ratioNaVTwoPerNaVSix': 0.47176632724853385,
        'x=0.6|𝜅=0.2|ratioNaVTwoPerTotalNaV': 0.32054431366866554,
        'x=0.6|𝜅=0.3|GbarNaVSix': 89.41246232367705,
        'x=0.6|𝜅=0.3|GbarNaVTwo': 54.818452484072324,
        'x=0.6|𝜅=0.3|ratioNaVSixPerNaVTwo': 1.6310650569651912,
        'x=0.6|𝜅=0.3|ratioNaVSixPerTotalNaV': 0.6199257797321619,
        'x=0.6|𝜅=0.3|ratioNaVTwoPerNaVSix': 0.6130963297445843,
        'x=0.6|𝜅=0.3|ratioNaVTwoPerTotalNaV': 0.3800742202678381,
        'x=0.6|𝜅=0.4|GbarNaVSix': 80.76788819594653,
        'x=0.6|𝜅=0.4|GbarNaVTwo': 63.4630266118028,
        'x=0.6|𝜅=0.4|ratioNaVSixPerNaVTwo': 1.2726762732259191,
        'x=0.6|𝜅=0.4|ratioNaVSixPerTotalNaV': 0.5599901262749737,
        'x=0.6|𝜅=0.4|ratioNaVTwoPerNaVSix': 0.785745771361988,
        'x=0.6|𝜅=0.4|ratioNaVTwoPerTotalNaV': 0.4400098737250262,
        'x=0.6|𝜅=0.5|GbarNaVSix': 72.11545822669046,
        'x=0.6|𝜅=0.5|GbarNaVTwo': 72.11545658105892,
        'x=0.6|𝜅=0.5|ratioNaVSixPerNaVTwo': 1.0000000228194013,
        'x=0.6|𝜅=0.5|ratioNaVSixPerTotalNaV': 0.5000000057048503,
        'x=0.6|𝜅=0.5|ratioNaVTwoPerNaVSix': 0.9999999771805992,
        'x=0.6|𝜅=0.5|ratioNaVTwoPerTotalNaV': 0.4999999942951497,
        'x=0.6|𝜅=0.6|GbarNaVSix': 63.463032799392835,
        'x=0.6|𝜅=0.6|GbarNaVTwo': 80.76788200835651,
        'x=0.6|𝜅=0.6|ratioNaVSixPerNaVTwo': 0.7857459081671443,
        'x=0.6|𝜅=0.6|ratioNaVSixPerTotalNaV': 0.4400099166256071,
        'x=0.6|𝜅=0.6|ratioNaVTwoPerNaVSix': 1.2726760516419762,
        'x=0.6|𝜅=0.6|ratioNaVTwoPerTotalNaV': 0.5599900833743928,
        'x=0.6|𝜅=0.7|GbarNaVSix': 54.81849720051365,
        'x=0.6|𝜅=0.7|GbarNaVTwo': 89.41241760723568,
        'x=0.6|𝜅=0.7|ratioNaVSixPerNaVTwo': 0.6130971364773552,
        'x=0.6|𝜅=0.7|ratioNaVSixPerTotalNaV': 0.3800745303015184,
        'x=0.6|𝜅=0.7|ratioNaVTwoPerNaVSix': 1.631062910757756,
        'x=0.6|𝜅=0.7|ratioNaVTwoPerTotalNaV': 0.6199254696984815,
        'x=0.6|𝜅=0.8|GbarNaVSix': 46.232719476326075,
        'x=0.6|𝜅=0.8|GbarNaVTwo': 97.99819533142326,
        'x=0.6|𝜅=0.8|ratioNaVSixPerNaVTwo': 0.4717711312945116,
        'x=0.6|𝜅=0.8|ratioNaVSixPerTotalNaV': 0.3205465314974349,
        'x=0.6|𝜅=0.8|ratioNaVTwoPerNaVSix': 2.1196718783026425,
        'x=0.6|𝜅=0.8|ratioNaVTwoPerTotalNaV': 0.679453468502565,
        'x=0.7|𝜅=0.2|GbarNaVSix': 102.31235817873241,
        'x=0.7|𝜅=0.2|GbarNaVTwo': 41.91855662901696,
        'x=0.7|𝜅=0.2|ratioNaVSixPerNaVTwo': 2.440741437836376,
        'x=0.7|𝜅=0.2|ratioNaVSixPerTotalNaV': 0.709364967386557,
        'x=0.7|𝜅=0.2|ratioNaVTwoPerNaVSix': 0.4097115673532637,
        'x=0.7|𝜅=0.2|ratioNaVTwoPerTotalNaV': 0.29063503261344303,
        'x=0.7|𝜅=0.3|GbarNaVSix': 92.29529647697746,
        'x=0.7|𝜅=0.3|GbarNaVTwo': 51.9356183307719,
        'x=0.7|𝜅=0.3|ratioNaVSixPerNaVTwo': 1.777109803317628,
        'x=0.7|𝜅=0.3|ratioNaVSixPerTotalNaV': 0.6399134096875224,
        'x=0.7|𝜅=0.3|ratioNaVTwoPerNaVSix': 0.5627114307360934,
        'x=0.7|𝜅=0.3|ratioNaVTwoPerTotalNaV': 0.36008659031247764,
        'x=0.7|𝜅=0.4|GbarNaVSix': 82.20995999462514,
        'x=0.7|𝜅=0.4|GbarNaVTwo': 62.02095481312422,
        'x=0.7|𝜅=0.4|ratioNaVSixPerNaVTwo': 1.3255190966074701,
        'x=0.7|𝜅=0.4|ratioNaVSixPerTotalNaV': 0.5699884806541358,
        'x=0.7|𝜅=0.4|ratioNaVTwoPerNaVSix': 0.7544214206791869,
        'x=0.7|𝜅=0.4|ratioNaVTwoPerTotalNaV': 0.43001151934586435,
        'x=0.7|𝜅=0.5|GbarNaVSix': 72.11545836382642,
        'x=0.7|𝜅=0.5|GbarNaVTwo': 72.11545644392287,
        'x=0.7|𝜅=0.5|ratioNaVSixPerNaVTwo': 1.000000026622636,
        'x=0.7|𝜅=0.5|ratioNaVSixPerTotalNaV': 0.5000000066556589,
        'x=0.7|𝜅=0.5|ratioNaVTwoPerNaVSix': 0.9999999733773646,
        'x=0.7|𝜅=0.5|ratioNaVTwoPerTotalNaV': 0.4999999933443411,
        'x=0.7|𝜅=0.6|GbarNaVSix': 62.02096203197923,
        'x=0.7|𝜅=0.6|GbarNaVTwo': 82.2099527757702,
        'x=0.7|𝜅=0.6|ratioNaVSixPerNaVTwo': 0.7544215747349111,
        'x=0.7|𝜅=0.6|ratioNaVSixPerTotalNaV': 0.43001156939654167,
        'x=0.7|𝜅=0.6|ratioNaVTwoPerNaVSix': 1.3255188259314832,
        'x=0.7|𝜅=0.6|ratioNaVTwoPerTotalNaV': 0.5699884306034584,
        'x=0.7|𝜅=0.7|GbarNaVSix': 51.93567049995351,
        'x=0.7|𝜅=0.7|GbarNaVTwo': 92.2952443077959,
        'x=0.7|𝜅=0.7|ratioNaVSixPerNaVTwo': 0.5627123140467885,
        'x=0.7|𝜅=0.7|ratioNaVSixPerTotalNaV': 0.3600869520184382,
        'x=0.7|𝜅=0.7|ratioNaVTwoPerNaVSix': 1.7771070137214944,
        'x=0.7|𝜅=0.7|ratioNaVTwoPerTotalNaV': 0.6399130479815618,
        'x=0.7|𝜅=0.8|GbarNaVSix': 41.91892982173466,
        'x=0.7|𝜅=0.8|GbarNaVTwo': 102.31198498601478,
        'x=0.7|𝜅=0.8|ratioNaVSixPerNaVTwo': 0.40971670941058014,
        'x=0.7|𝜅=0.8|ratioNaVSixPerTotalNaV': 0.29063762008034066,
        'x=0.7|𝜅=0.8|ratioNaVTwoPerNaVSix': 2.440710805860477,
        'x=0.7|𝜅=0.8|ratioNaVTwoPerTotalNaV': 0.7093623799196594,
        'x=0.8|𝜅=0.2|GbarNaVSix': 106.62620114656924,
        'x=0.8|𝜅=0.2|GbarNaVTwo': 37.60471366118014,
        'x=0.8|𝜅=0.2|ratioNaVSixPerNaVTwo': 2.835447760811989,
        'x=0.8|𝜅=0.2|ratioNaVSixPerTotalNaV': 0.7392742484417794,
        'x=0.8|𝜅=0.2|ratioNaVTwoPerNaVSix': 0.3526779839927749,
        'x=0.8|𝜅=0.2|ratioNaVTwoPerTotalNaV': 0.26072575155822053,
        'x=0.8|𝜅=0.3|GbarNaVSix': 95.17813063027799,
        'x=0.8|𝜅=0.3|GbarNaVTwo': 49.05278417747145,
        'x=0.8|𝜅=0.3|ratioNaVSixPerNaVTwo': 1.9403206612274335,
        'x=0.8|𝜅=0.3|ratioNaVSixPerTotalNaV': 0.6599010396428833,
        'x=0.8|𝜅=0.3|ratioNaVTwoPerNaVSix': 0.515378730939971,
        'x=0.8|𝜅=0.3|ratioNaVTwoPerTotalNaV': 0.3400989603571167,
        'x=0.8|𝜅=0.4|GbarNaVSix': 83.6520317933038,
        'x=0.8|𝜅=0.4|GbarNaVTwo': 60.57888301444558,
        'x=0.8|𝜅=0.4|ratioNaVSixPerNaVTwo': 1.3808777519611286,
        'x=0.8|𝜅=0.4|ratioNaVSixPerTotalNaV': 0.579986835033298,
        'x=0.8|𝜅=0.4|ratioNaVTwoPerNaVSix': 0.724177066782194,
        'x=0.8|𝜅=0.4|ratioNaVTwoPerTotalNaV': 0.42001316496670194,
        'x=0.8|𝜅=0.5|GbarNaVSix': 72.11545850096239,
        'x=0.8|𝜅=0.5|GbarNaVTwo': 72.11545630678692,
        'x=0.8|𝜅=0.5|ratioNaVSixPerNaVTwo': 1.0000000304258696,
        'x=0.8|𝜅=0.5|ratioNaVSixPerTotalNaV': 0.5000000076064672,
        'x=0.8|𝜅=0.5|ratioNaVTwoPerNaVSix': 0.9999999695741314,
        'x=0.8|𝜅=0.5|ratioNaVTwoPerTotalNaV': 0.4999999923935327,
        'x=0.8|𝜅=0.6|GbarNaVSix': 60.57889126456557,
        'x=0.8|𝜅=0.6|GbarNaVTwo': 83.65202354318373,
        'x=0.8|𝜅=0.6|ratioNaVSixPerNaVTwo': 0.724177236827904,
        'x=0.8|𝜅=0.6|ratioNaVSixPerTotalNaV': 0.42001322216747644,
        'x=0.8|𝜅=0.6|ratioNaVTwoPerNaVSix': 1.3808774277140714,
        'x=0.8|𝜅=0.6|ratioNaVTwoPerTotalNaV': 0.5799867778325236,
        'x=0.8|𝜅=0.7|GbarNaVSix': 49.05284379939331,
        'x=0.8|𝜅=0.7|GbarNaVTwo': 95.178071008356,
        'x=0.8|𝜅=0.7|ratioNaVSixPerNaVTwo': 0.5153796802110729,
        'x=0.8|𝜅=0.7|ratioNaVSixPerTotalNaV': 0.340099373735358,
        'x=0.8|𝜅=0.7|ratioNaVTwoPerNaVSix': 1.9403170873761486,
        'x=0.8|𝜅=0.7|ratioNaVTwoPerTotalNaV': 0.659900626264642,
        'x=0.8|𝜅=0.8|GbarNaVSix': 37.60514016714322,
        'x=0.8|𝜅=0.8|GbarNaVTwo': 106.62577464060608,
        'x=0.8|𝜅=0.8|ratioNaVSixPerNaVTwo': 0.352683394741051,
        'x=0.8|𝜅=0.8|ratioNaVSixPerTotalNaV': 0.2607287086632467,
        'x=0.8|𝜅=0.8|ratioNaVTwoPerNaVSix': 2.835404260340142,
        'x=0.8|𝜅=0.8|ratioNaVTwoPerTotalNaV': 0.7392712913367534,
        'x=0.9|𝜅=0.2|GbarNaVSix': 110.94004411440605,
        'x=0.9|𝜅=0.2|GbarNaVTwo': 33.29087069334332,
        'x=0.9|𝜅=0.2|ratioNaVSixPerNaVTwo': 3.332446457658708,
        'x=0.9|𝜅=0.2|ratioNaVSixPerTotalNaV': 0.7691835294970019,
        'x=0.9|𝜅=0.2|ratioNaVTwoPerNaVSix': 0.3000798400531766,
        'x=0.9|𝜅=0.2|ratioNaVTwoPerTotalNaV': 0.23081647050299814,
        'x=0.9|𝜅=0.3|GbarNaVSix': 98.06096478357826,
        'x=0.9|𝜅=0.3|GbarNaVTwo': 46.169950024171094,
        'x=0.9|𝜅=0.3|ratioNaVSixPerNaVTwo': 2.123913167162645,
        'x=0.9|𝜅=0.3|ratioNaVSixPerTotalNaV': 0.6798886695982432,
        'x=0.9|𝜅=0.3|ratioNaVTwoPerNaVSix': 0.47082904115892327,
        'x=0.9|𝜅=0.3|ratioNaVTwoPerTotalNaV': 0.3201113304017568,
        'x=0.9|𝜅=0.4|GbarNaVSix': 85.09410359198246,
        'x=0.9|𝜅=0.4|GbarNaVTwo': 59.136811215766926,
        'x=0.9|𝜅=0.4|ratioNaVSixPerNaVTwo': 1.438936287611241,
        'x=0.9|𝜅=0.4|ratioNaVSixPerTotalNaV': 0.5899851894124604,
        'x=0.9|𝜅=0.4|ratioNaVTwoPerNaVSix': 0.6949578022388238,
        'x=0.9|𝜅=0.4|ratioNaVTwoPerTotalNaV': 0.4100148105875396,
        'x=0.9|𝜅=0.5|GbarNaVSix': 72.11545863809835,
        'x=0.9|𝜅=0.5|GbarNaVTwo': 72.11545616965104,
        'x=0.9|𝜅=0.5|ratioNaVSixPerNaVTwo': 1.000000034229102,
        'x=0.9|𝜅=0.5|ratioNaVSixPerTotalNaV': 0.5000000085572753,
        'x=0.9|𝜅=0.5|ratioNaVTwoPerNaVSix': 0.9999999657708991,
        'x=0.9|𝜅=0.5|ratioNaVTwoPerTotalNaV': 0.4999999914427246,
        'x=0.9|𝜅=0.6|GbarNaVSix': 59.13682049715195,
        'x=0.9|𝜅=0.6|GbarNaVTwo': 85.09409431059728,
        'x=0.9|𝜅=0.6|ratioNaVSixPerNaVTwo': 0.6949579871113016,
        'x=0.9|𝜅=0.6|ratioNaVSixPerTotalNaV': 0.41001487493841127,
        'x=0.9|𝜅=0.6|ratioNaVTwoPerNaVSix': 1.4389359048259187,
        'x=0.9|𝜅=0.6|ratioNaVTwoPerTotalNaV': 0.5899851250615887,
        'x=0.9|𝜅=0.7|GbarNaVSix': 46.17001709883314,
        'x=0.9|𝜅=0.7|GbarNaVTwo': 98.0608977089162,
        'x=0.9|𝜅=0.7|ratioNaVSixPerNaVTwo': 0.47083004722110683,
        'x=0.9|𝜅=0.7|ratioNaVSixPerTotalNaV': 0.32011179545227764,
        'x=0.9|𝜅=0.7|ratioNaVTwoPerNaVSix': 2.1239086288186475,
        'x=0.9|𝜅=0.7|ratioNaVTwoPerTotalNaV': 0.6798882045477224,
        'x=0.9|𝜅=0.8|GbarNaVSix': 33.29135051255177,
        'x=0.9|𝜅=0.8|GbarNaVTwo': 110.93956429519753,
        'x=0.9|𝜅=0.8|ratioNaVSixPerNaVTwo': 0.30008546296402683,
        'x=0.9|𝜅=0.8|ratioNaVSixPerTotalNaV': 0.23081979724615242,
        'x=0.9|𝜅=0.8|ratioNaVTwoPerNaVSix': 3.332384015282594,
        'x=0.9|𝜅=0.8|ratioNaVTwoPerTotalNaV': 0.7691802027538477,
        'x=1.0|𝜅=0.2|GbarNaVSix': 115.25388708224288,
        'x=1.0|𝜅=0.2|GbarNaVTwo': 28.977027725506503,
        'x=1.0|𝜅=0.2|ratioNaVSixPerNaVTwo': 3.977422673368005,
        'x=1.0|𝜅=0.2|ratioNaVSixPerTotalNaV': 0.7990928105522243,
        'x=1.0|𝜅=0.2|ratioNaVTwoPerNaVSix': 0.25141909274460367,
        'x=1.0|𝜅=0.2|ratioNaVTwoPerTotalNaV': 0.20090718944777566,
        'x=1.0|𝜅=0.3|GbarNaVSix': 100.94379893687871,
        'x=1.0|𝜅=0.3|GbarNaVTwo': 43.28711587087063,
        'x=1.0|𝜅=0.3|ratioNaVSixPerNaVTwo': 2.3319594504286947,
        'x=1.0|𝜅=0.3|ratioNaVSixPerTotalNaV': 0.699876299553604,
        'x=1.0|𝜅=0.3|ratioNaVTwoPerNaVSix': 0.42882392308158074,
        'x=1.0|𝜅=0.3|ratioNaVTwoPerTotalNaV': 0.300123700446396,
        'x=1.0|𝜅=0.4|GbarNaVSix': 86.53617539066106,
        'x=1.0|𝜅=0.4|GbarNaVTwo': 57.69473941708828,
        'x=1.0|𝜅=0.4|ratioNaVSixPerNaVTwo': 1.4998971529288232,
        'x=1.0|𝜅=0.4|ratioNaVSixPerTotalNaV': 0.5999835437916226,
        'x=1.0|𝜅=0.4|ratioNaVTwoPerNaVSix': 0.6667123796103736,
        'x=1.0|𝜅=0.4|ratioNaVTwoPerTotalNaV': 0.4000164562083774,
        'x=1.0|𝜅=0.5|GbarNaVSix': 72.11545877523432,
        'x=1.0|𝜅=0.5|GbarNaVTwo': 72.11545603251504,
        'x=1.0|𝜅=0.5|ratioNaVSixPerNaVTwo': 1.0000000380323362,
        'x=1.0|𝜅=0.5|ratioNaVSixPerTotalNaV': 0.500000009508084,
        'x=1.0|𝜅=0.5|ratioNaVTwoPerNaVSix': 0.9999999619676652,
        'x=1.0|𝜅=0.5|ratioNaVTwoPerTotalNaV': 0.4999999904919162,
        'x=1.0|𝜅=0.6|GbarNaVSix': 57.6947497297383,
        'x=1.0|𝜅=0.6|GbarNaVTwo': 86.5361650780111,
        'x=1.0|𝜅=0.6|ratioNaVSixPerNaVTwo': 0.6667125782350919,
        'x=1.0|𝜅=0.6|ratioNaVSixPerTotalNaV': 0.4000165277093452,
        'x=1.0|𝜅=0.6|ratioNaVTwoPerNaVSix': 1.4998967060846218,
        'x=1.0|𝜅=0.6|ratioNaVTwoPerTotalNaV': 0.5999834722906547,
        'x=1.0|𝜅=0.7|GbarNaVSix': 43.287190398272955,
        'x=1.0|𝜅=0.7|GbarNaVTwo': 100.94372440947646,
        'x=1.0|𝜅=0.7|ratioNaVSixPerNaVTwo': 0.4288249779914918,
        'x=1.0|𝜅=0.7|ratioNaVSixPerTotalNaV': 0.3001242171691972,
        'x=1.0|𝜅=0.7|ratioNaVTwoPerNaVSix': 2.3319537138059174,
        'x=1.0|𝜅=0.7|ratioNaVTwoPerTotalNaV': 0.699875782830803,
        'x=1.0|𝜅=0.8|GbarNaVSix': 28.977560857960345,
        'x=1.0|𝜅=0.8|GbarNaVTwo': 115.25335394978902,
        'x=1.0|𝜅=0.8|ratioNaVSixPerNaVTwo': 0.25142488148834813,
        'x=1.0|𝜅=0.8|ratioNaVSixPerTotalNaV': 0.20091088582905814,
        'x=1.0|𝜅=0.8|ratioNaVTwoPerNaVSix': 3.977331098180684,
        'x=1.0|𝜅=0.8|ratioNaVTwoPerTotalNaV': 0.7990891141709418
        }


        ratio = conductances_dict[f"x={x}|𝜅={kappa}|ratioNaVTwoPerTotalNaV"]
        # ratio = conductances_dict[f"x={x}|𝜅={kappa}|ratioNaVSixPerTotalNaV"]
        # ratio = conductances_dict[f"x={x}|𝜅={kappa}|ratioNaVTwoPerNaVSix"]
        # ratio = conductances_dict[f"x={x}|𝜅={kappa}|ratioNaVSixPerNaVTwo"]
        # ratio = conductances_dict[f"x={x}|𝜅={kappa}|GbarNaVTwo"]
        # ratio = conductances_dict[f"x={x}|𝜅={kappa}|GbarNaVSix"]
        
        # print_args(ratio)
        return ratio











"""Demo"""
if __name__=='__main__':
    import Hay_data as hd

    ###____________________________________________________________________________________
    """Plotting effects of Nav distribution on backpropagation threshold""" 

    ###Fig 7a
    pSomatic = DisplayHayData(hd.somaticStim_negative70mV_1msPulse_new, filename='Fig7A')
    pSomatic.plot_thresh(select_kappa_vals=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
    # pSomatic.plot_thresh(select_kappa_vals=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], use_ratio_of_total_AIS_NaV_conductances_from_each_subtype=True)
    
    ###Fig 7b
    pAxonal = DisplayHayData(hd.axonalStim_negative70mV_1msPulse_new, filename='Fig7B')
    pAxonal.plot_thresh(select_kappa_vals=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
    # pAxonal.plot_thresh(select_kappa_vals=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], use_ratio_of_total_AIS_NaV_conductances_from_each_subtype=True)


    ###____________________________________________________________________________________
    """ploting 'scaleNavTwo' results"""

    ###Fig 8a
    pSomatic = DisplayHayData(hd.somaticStim_scaleNavTwo__sweepKappa, xAxisVariableName='scaleNavTwo', filename='Fig8A')
    pSomatic.plot_thresh(legend_location='lower right')

    ###Fig 8b
    data = hd.axonalStim_scaleNavTwo__sweepKappa
    data_lowend = hd.lowend_axonalStim_scaleNaVTwo
    for key, value in data.items():
        data[key].extend(data_lowend[key])
    pAxonal = DisplayHayData(data, xAxisVariableName='scaleNavTwo', filename='Fig8B')
    pAxonal.plot_thresh(legend_location='lower right', vline=0.32, ylim=(0.57, 0.606), xlim=(-0.03, None))




