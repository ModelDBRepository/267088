"""Credit to Robert McDougal for providing the h.SaveState() example code upon which the initialize() (below) is based,
and for his assistance with numerous other implementation hurdles."""

import copy
import math
from sorcery import print_args
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from neuron import h, gui
from neuron.units import ms, mV
import pprint
from datetime import datetime
#====================== General files and tools =====================
h.load_file("nrngui.hoc")
#====================== cvode =======================================
h.dt = 1.0e-01
cvode = h.CVode()
cvode.active(1)
cvode.atol(1.0e-06)


state_filename='warmed_up.dat'
v_init = -80.53 * mV
t_warmup = 3000.0 * ms
t_stop = t_warmup + 1000.0 * ms


###define a python function to restore from saved state. 
def restore_state(state_filename=state_filename):
    s = h.SaveState()
    f = h.File(state_filename) #<---Requires that this script be run at least once, otherwise this file will not exist yet.
    s.fread(f)
    f.close()
    s.restore()

###create a function to setup and optionally warmup and save the simulation state
def initialize(warm_up_and_save=False, state_filename=state_filename, load_state=True):
    global fih
###save the warmed-up simulation (calcium equilibrated)
    if warm_up_and_save is True:
        h.finitialize(v_init)
        h.continuerun(t_warmup)
        # plt.plot(cell1.t, cell1.vsoma, label="original", linestyle='--', linewidth=3.0, color='black')
        s = h.SaveState()
        s.save()
        f = h.File(state_filename)
        s.fwrite(f)
        f.close()
    if load_state is True:
        cvode.active(1)
        fih = h.FInitializeHandler(restore_state)


class Pyramidal:
    cells = [] # Class variable to store all neurons (instances) created by Pyramidal
    def __init__(self, gid):
        self.gid = gid
        Pyramidal.cells.append(self) 
        self.tropical_colors = [
        '#FF4500',  # Coral
        '#FFD700',  # Gold
        '#32CD32',  # LimeGreen
        '#20B2AA',  # LightSeaGreen
        '#FF69B4',  # HotPink
        '#FF6347',  # Tomato
        '#FFDAB9',  # PeachPuff
        '#98FB98',  # PaleGreen
        '#F08080',  # LightCoral
        '#00CED1',  # DarkTurquoise
        '#ADFF2F',  # GreenYellow
        '#87CEEB',  # SkyBlue
        '#FFA07A',  # LightSalmon
        '#B0E0E6',  # PowderBlue
        '#FFDEAD',  # NavajoWhite
        '#FFB6C1',  # LightPink
        '#FFA500',  # Orange
        '#BDB76B',  # DarkKhaki
        '#8A2BE2',  # BlueViolet
        '#66CDAA',  # MediumAquamarine
        '#FA8072',  # Salmon
        '#7FFFD4',  # Aquamarine
        '#F0E68C',  # Khaki
        '#F5DEB3',  # Wheat
        '#EE82EE',  # Violet
        '#FAFAD2',  # LightGoldenrodYellow
        '#D8BFD8',  # Thistle
        '#DDA0DD',  # Plum
        '#CD853F',  # Peru
        '#FFC0CB',  # Pink
        '#DA70D6',  # Orchid
        '#EEE8AA',  # PaleGoldenrod
        '#40E0D0',  # Turquoise
        '#FFFAF0',  # FloralWhite
        '#AFEEEE',  # PaleTurquoise
        '#F5F5DC',  # Beige
        '#FFEFDB',  # PapayaWhip
        '#FFEBCD',  # BlanchedAlmond
        '#FFE4C4',  # Bisque
        '#FFB5C5',  # HotPinkLight
        '#F4A460',  # SandyBrown
        '#D2B48C',  # Tan
        '#C0C0C0',  # Silver
        '#A52A2A',  # Brown
        '#9ACD32',  # YellowGreen
        '#FAEBD7',  # AntiqueWhite
        '#8B4513',  # SaddleBrown
        '#6B8E23',  # OliveDrab
        '#808080',  # Gray
        '#BC8F8F'   # RosyBrown
        ]

        self.v_init = v_init
        self.t_warmup = t_warmup
        self.t_stop = t_stop
        h.load_file("import3d.hoc")
        # h.load_file("models/L5PCbiophys1.hoc")
        h.load_file("models/L5PCbiophys4_mixedAISNavs.hoc")
        h.load_file("models/L5PCtemplate.hoc")
        self.cell = h.L5PCtemplate("morphologies/cell1.asc")
        self.soma = self.cell.soma[0]
        # self.all = self.soma.wholetree()

        ### Identify apical segments for observation
         # self.apic = self.cell.apic[37]
        self.apic = self.cell.apic[36]
        self.apic_BAC = self.cell.apic[37]
        self.apic_BAC_seg = self.cell.apic[37](0.6)

        ### Identify axon and AIS
        self.axon0 = self.cell.axon[0]
        self.axon1 = self.cell.axon[1]
        self.apic0 = self.cell.apic[0]
    
        self.axon0.nseg = 31
        self.axon1.nseg = 31
        self.axon0.L = 30.0 #um
        self.axon1.L = 30.0 #um 

        self.ais_segments = [seg for seg in self.axon0] + [seg for seg in self.axon1]
        self.sigma = 10.0 #<---called "flatness" in Thresher.py, this parameter multiplies the argument of the tanh() function in the Nav profiles.
        self.aisNavSeparation = 0.0
        self.aisNavCrossover = 0.5
        self.aisDeltaRightShift = 0.0
        self.get_aisNavProperties()

        ###Passive cable: AIS in original Hay model is a dead end: there was no axon attached to it.
        ###Adding a passive cable is the least intrusive means of adding an axon,
        ###i.e. providing the proper boundary conditions at the end of the AIS, without interfering with the model tuning. 
        self.passive = h.Section('passive')
        self.passive.Ra = 100.0 # Ohm cm
        self.passive.cm = 1.0
        self.passive.diam = 1.0 #um
        self.passive.L = 400.0 #um
        self.passive.nseg = 2*math.ceil(0.5*self.passive.L) + 1
        self.passive.connect(self.ais_segments[-1].sec(1))
        h.pas.insert(self.passive)
        for seg in self.passive:
            seg.pas.e = v_init
            seg.pas.g = 0.01

        self.area_soma = self.get_area_of_section(self.soma)
        self.area_ais = self.get_area_of_section(self.ais_segments)
        self.gNaVTotal_soma = self.get_gNaVTotal(self.soma(0.5))
        self.gNaVTotal_ais = self.get_gNaVTotal(self.ais_segments[1])




        ### initialize current clamp with zero amplitude
        self.ic = h.IClamp(self.soma(0.5)) #h.IClamp(axon1(0.5)) #h.IClamp(axon1(1.0))
        self.ic.amp = 0.0 #nA
        self.amplitudes = [self.ic.amp]
        self.vpeak_data = {}
        self.tpeak_data = {}

        ### Insert mechanism to record vpeak in every segment of the model.
        self.all = self.soma.wholetree()
        for sec in self.all:
            h.peak.insert(sec)

        ### default parameters for generating threshold curves
        self.initialize_threshold_dict()
        self.IMAX = 20.0 #nA  #<---for threshold finder
        self.vSpikeThreshold = -70.0 #10.0<---somatic BAP criterion  #mV
        self.spike_recording_site = self.apic_BAC_seg #self.soma(0.5) #self.apic(0.99)
        self.spikeFound = False
        self.epsilon = 1.0e-4
        self.i_spikeFound = None
        self.vpeak = None
        self.separation_vals = [0.0, 1.0]
        self.crossover_vals = [0.5]
        self.scaleNavTwo_vals = [1.0]
        self.scaleNavTwo = 1.0
        self.deltaRightShift_vals = [0.0]
        self.DeltaVrs = 0.0

        
        self.data_dict = {}
        self.setup_recording_vectors()


        h.define_shape()

    def __str__(self):
        return f"Pyramidal[{self.gid}]"

    def get_area_of_section(self, section=None):
        if section is not None:
            area_section = 0.0
            for seg in section:
                area_section+=seg.area()
            return area_section




    def stim(self, location=False, amplitude=False, duration=False, delay=False):
        ###move current clamp to "location"
        if location is not False:
            self.ic.loc(location)

        if amplitude is not False:
            self.ic.amp = amplitude
        else:
            self.ic.amp = 0.0 #nA

        if duration is not False:
            self.ic.dur = duration #ms

        # else:
        #     self.ic.dur = 1.0 * ms

        if delay is not False:
            self.ic.delay = delay #ms
        else:
            self.ic.delay = t_warmup + 100.0*ms

    def run_amplitudes(self, amplitudes=None):
        if amplitudes is not None:
            self.amplitudes = amplitudes

        import matplotlib.pyplot as plt
        # Create Figure 1 and Axes 1
        self.fig1, self.ax1 = plt.subplots()
        
        for amplitude in self.amplitudes:
            self.stim(amplitude=amplitude)
            h.finitialize(v_init)
            self.reset_peak()
            h.continuerun(t_stop)
            # plt.plot(self.t, self.caisoma, label=f"gnabar={self.soma(0.5).NaTs2_t.gNaTs2_tbar}")
            self.ax1.plot(self.t, self.vsoma, label=f"{self.ic.amp:g}")
            # plt.plot(self.t, self.currentvec)

        # zoom the x-axis to view the spikes
        self.ax1.set_xlim([3098.0, 3110.0])

        # create legend
        self.ax1.legend().set_title('1ms pulse amplitude [nA]')
        # indicate stimulation site in title
        self.ax1.set_title(f"Nav Separation = {self.aisNavSeparation:g}"+f" | Crossover = {self.aisNavCrossover:g}"+r' | Stimulating at '+str(self.ic.get_segment().sec.name()).replace('L5PCtemplate', '').replace('[0]', '').replace('[1]', '').replace('.', ''))
        # axis labels
        self.ax1.set_ylabel(r'$V_{soma} \  \ \ [mV]$', fontsize=20)
        self.ax1.set_xlabel(r'$t \ \ \ [ms]$', fontsize=20)

        # Show Figure 1
        self.fig1.show()

    def spike(self, plot=False):       
        self.stim(amplitude=self.ic.amp)
        h.finitialize(v_init)
        # self.initialize()
        self.reset_peak()
        h.continuerun(t_stop)


        self.t_postWarmUp = np.array(self.t)[np.where(np.array(self.t)>self.t_warmup)]
        self.spikeSignal = np.array(self.vSpikeRecording)[np.where(np.array(self.t)>self.t_warmup)]
        self.vMax = np.max(self.spikeSignal)
        # print_args(self.aisNavSeparation, self.aisNavCrossover)
        print(f"x={self.aisNavSeparation:g}, 𝜿={self.aisNavCrossover:g}, I={self.ic.amp:g}")
        # print_args(self.ic.amp)
        print_args(self.vMax)
        # print_args(self.apic_BAC_seg.peak.vpeak)
        if self.vMax >= self.vSpikeThreshold:
            self.spikeFound = True
            self.i_spikeFound = self.ic.amp
            self.vpeak = self.vMax
        else:
            self.spikeFound = False
            

        if plot==True:
            import matplotlib.pyplot as plt
            # Create Figure 1 and Axes 1
            self.fig1, self.ax1 = plt.subplots()
            # plt.plot(self.t, self.caisoma, label=f"gnabar={self.soma(0.5).NaTs2_t.gNaTs2_tbar}")
            # self.ax1.plot(self.t, self.vsoma, label=f"{self.ic.amp:g}")
            # plt.plot(self.t, self.currentvec)
            self.ax1.plot(self.t_postWarmUp, self.spikeSignal, label=f"{self.ic.amp:g}")

            # zoom the x-axis to view the spikes
            self.ax1.set_xlim([3098.0, 3110.0])

            # create legend
            self.ax1.legend().set_title('1ms pulse amplitude [nA]')
            # indicate stimulation site in title
            self.ax1.set_title(f"Nav Separation = {self.aisNavSeparation:g}"+f" | Crossover = {self.aisNavCrossover:g}"+r' | Stimulating at '+str(self.ic.get_segment().sec.name()).replace('L5PCtemplate', '').replace('[0]', '').replace('[1]', '').replace('.', ''))
            # axis labels
            self.ax1.set_ylabel(r'$V_{soma} \  \ \ [mV]$', fontsize=20)
            self.ax1.set_xlabel(r'$t \ \ \ [ms]$', fontsize=20)

        
            # Show Figure 1
            self.fig1.show()
        # else:
        #     plt.close(fig1)

        return self.spikeFound



    def ithresh(self, IMAX=None, epsilon=None):
        """Theshold finder"""
        self.i_spikeFound = None
        self.spikeFound = False
        if epsilon is not None:
            self.epsilon = epsilon # Tolerance for convergence
        if IMAX is not None:
            self.IMAX = IMAX
        delta = (self.IMAX - self.epsilon)/5.1  # Initial step size

        self.ic.amp = self.IMAX
        while True:
            if self.spike():
                self.ic.amp -= delta
                continue
            else:
                self.ic.amp += delta

            if self.ic.amp <= 0.0 or self.ic.amp >= self.IMAX:
                print(f"Error, no threshold found. Pulse amplitude was set to {self.ic.amp:g} upon quitting")
                break

            if abs(delta) < self.epsilon:
                break

            delta = delta / 2.0  # Update step size (bisection method)


    def threshCurve(self, separation_vals=None, crossover_vals=None, deltaRightShift_vals=None, scaleNavTwo_vals=None):
        if separation_vals is not None:
            self.separation_vals = separation_vals
        if crossover_vals is not None:
            self.crossover_vals = crossover_vals
        if deltaRightShift_vals is not None:
            self.deltaRightShift_vals = deltaRightShift_vals
        if scaleNavTwo_vals is not None:
            self.scaleNavTwo_vals = scaleNavTwo_vals
        self.initialize_threshold_dict()
        for separation in self.separation_vals:
            for crossover in self.crossover_vals:
                for DeltaVrs in self.deltaRightShift_vals:
                    for scaleNavTwo in self.scaleNavTwo_vals:
                        self.NaV2scaleNavTwo(scaleNavTwo)
                        self.NaV2DeltaVrs(DeltaVrs)
                        self.set_aisNavProfile(separation=separation, crossover=crossover)
                        self.ithresh()
                        self.update_threshold_dict()

    def initialize_threshold_dict(self):
        self.threshold_dict = {}
        self.threshold_dict['spike_recording_site'] = []
        self.threshold_dict['electrode_location'] = []
        self.threshold_dict['separation'] = []
        self.threshold_dict['crossover'] = []
        self.threshold_dict['Ithresh'] = []
        self.threshold_dict['vpeak'] = []
        self.threshold_dict['DeltaVrs'] = []
        self.threshold_dict['scaleNavTwo'] = []

    def update_threshold_dict(self):
        self.threshold_dict['spike_recording_site'] += [str(self.spike_recording_site)]
        self.threshold_dict['electrode_location'] += [str(self.ic.get_segment())]
        self.threshold_dict['separation'] += [self.aisNavSeparation]
        self.threshold_dict['crossover'] += [self.aisNavCrossover]
        self.threshold_dict['Ithresh'] += [self.i_spikeFound]
        self.threshold_dict['vpeak'] += [self.vpeak]
        self.threshold_dict['DeltaVrs'] += [self.DeltaVrs]
        self.threshold_dict['scaleNavTwo'] += [self.scaleNavTwo]


    def plot_threshold_dict(self, threshold_dict=None):
        if threshold_dict is not None:
            data = self.threshold_dict
        else:
            data = threshold_dict
        import matplotlib.pyplot as plt
        # Create Figure 1 and Axes 1
        self.fig1, self.ax1 = plt.subplots()
        
        self.ax1.plot(data['separation'], data['Ithresh'], label='')

        # indicate stimulation site in title
        self.ax1.set_title(r'BackProp theshold as a function of $Na_{V}$ distribution')
        # axis labels
        self.ax1.set_ylabel(r'$\kappa$', fontsize=20)
        self.ax1.set_xlabel(r'$x$', fontsize=20)

        # Show Figure 1
        self.fig1.show()




    def create_recording(self, vecname, vecvar, vecsource):
        """method for recording data and automatically adding the new vector to data_dict. Adds the data vector called self.str(vecname) to the class"""
        setattr(self, vecname, h.Vector().record(getattr(vecsource, '_ref_'+vecvar)))
        datavec = getattr(self, vecname)
        self.data_dict[vecname] = datavec
        # return datavec
    def setup_recording_vectors(self):
        # self.t = h.Vector().record(h._ref_t)
        # self.vsoma = h.Vector().record(self.soma(0.5)._ref_v)

        self.create_recording(vecname='t', vecvar='t', vecsource=h)
        self.create_recording(vecname='currentvec', vecvar = 'i', vecsource = self.ic)
        self.create_recording(vecname='vsoma', vecvar='v', vecsource=self.soma(0.5)) #<--- creates the vector self.vsoma = h.Vector().record(self.soma(0.5)._ref_v) and adds it to self.data_dict
        self.create_recording(vecname='caisoma', vecvar='cai', vecsource=self.soma(0.5))
        self.create_recording(vecname='vaxon0', vecvar='v', vecsource=self.axon0(0.5))
        self.create_recording(vecname='vaxon1', vecvar='v', vecsource=self.axon1(0.5))
        self.create_recording(vecname='vapic', vecvar='v', vecsource=self.apic(0.99))
        self.create_recording(vecname='vSpikeRecording', vecvar='v', vecsource=self.spike_recording_site)

    def get_gNaVTotal(self, seg):
        gNaV2 = getattr(getattr(seg, 'NaTa2_t', None), 'gNaTa2_tbar', 0)
        gNaV6 = getattr(getattr(seg, 'NaTa6_t', None), 'gNaTa6_tbar', 0)
        gNaVSomatic_variant = getattr(getattr(seg, 'NaTs2_t', None), 'gNaTs2_tbar', 0)
        gNaVaxonal_variant = getattr(getattr(seg, 'NaTa_t', None), 'gNaTa_tbar', 0)
        return gNaV2 + gNaV6 + gNaVSomatic_variant + gNaVaxonal_variant



    def get_aisNavProperties(self):
        self.aisLength = h.distance(self.ais_segments[0].sec(0), self.ais_segments[-1].sec(1))
        self.ais_positions = []
        self.ais_normalized_positions = []
        self.Nav2profile = []
        self.Nav6profile = []
        self.aisDeltaRightShift_profile = []
        for seg in self.ais_segments:
            position_of_seg = h.distance(self.ais_segments[0].sec(0), seg)
            self.ais_positions += [position_of_seg]
            self.ais_normalized_positions += [position_of_seg/self.aisLength]
            self.Nav2profile += [seg.NaTa2_t.gNaTa2_tbar]
            self.Nav6profile += [seg.NaTa6_t.gNaTa6_tbar]
            if seg.NaTa2_t.vRS03 == seg.NaTa6_t.vRS03:
                self.aisDeltaRightShift_profile += [seg.NaTa2_t.vRS03]
            else:
                print(f"Error: aisDeltaRightShift values do NOT match for each NaV subtype present in AIS: seg.NaTa2_t.vRS03={seg.NaTa2_t.vRS03:g}, seg.NaTa6_t.vRS03={seg.NaTa6_t.vRS03:g}")
                exit()

        self.ais_gNavTotal = self.Nav2profile[0] + self.Nav6profile[0]

    ###Define a function to alter Nav2 right-shift
    def NaV2DeltaVrs(self, DeltaVrs=False):
        if DeltaVrs is not False:
            self.DeltaVrs = DeltaVrs
        for seg in self.ais_segments:
            seg.NaTa2_t.DeltaVrs = self.DeltaVrs

    def NaV2scaleNavTwo(self, scaleNavTwo=False):
        if scaleNavTwo is not False:
            self.scaleNavTwo = scaleNavTwo
        for seg in self.ais_segments:
            seg.NaTa2_t.scaleNavTwo = self.scaleNavTwo

    ###Define separate functions to create the density profile of each NaV subtype
    def Nav2Profile_func(self, segment):
            import numpy as np
            s = h.distance(self.ais_segments[0].sec(0), segment)/self.aisLength
            return self.ais_gNavTotal*(0.5)*(1.0 - self.aisNavSeparation*np.tanh(self.sigma*(s - self.aisNavCrossover)))
    def Nav6Profile_func(self, segment):
            import numpy as np
            s = h.distance(self.ais_segments[0].sec(0), segment)/self.aisLength
            return self.ais_gNavTotal*(0.5)*(1.0 + self.aisNavSeparation*np.tanh(self.sigma*(s - self.aisNavCrossover)))

    def set_aisNavProfile(self, separation=False, crossover=False, deltaRightShift=False):
        self.get_aisNavProperties()
        if separation is not False:
            self.aisNavSeparation = separation
        if crossover is not False:
            self.aisNavCrossover = crossover
        if deltaRightShift is not False:
            self.aisDeltaRightShift = deltaRightShift


        for seg in self.ais_segments:
            seg.NaTa2_t.gNaTa2_tbar = self.Nav2Profile_func(seg) #self.ais_gNavTotal
            seg.NaTa6_t.gNaTa6_tbar = self.Nav6Profile_func(seg) #self.ais_gNavTotal
            seg.NaTa2_t.vRS03 = self.aisDeltaRightShift
            seg.NaTa6_t.vRS03 = self.aisDeltaRightShift

        self.get_aisNavProperties()



    def plot_aisNavProfiles(self):
        import matplotlib.pyplot as plt
        self.get_aisNavProperties()


        # Create Figure 1 and Axes 1
        self.fig1, self.ax1 = plt.subplots()

        # Add data to Figure 1, Axes 1
        self.ax1.plot(self.ais_positions, np.array(self.Nav2profile)+np.array(self.Nav6profile), linewidth=1.0, color='black', label=r'$\bar{g}_{Na_{V}, Total} = \bar{g}_{Na_{V}1.2}+\bar{g}_{Na_{V}1.6}$')
        self.ax1.plot(self.ais_positions, self.Nav2profile, linewidth=3.0, label=r'$\bar{g}_{Na_{V}1.2}$')
        self.ax1.plot(self.ais_positions, self.Nav6profile, linewidth=2.0, linestyle='--', label=r'$\bar{g}_{Na_{V}1.6}$')

        # create legend
        self.ax1.legend()
        # create title
        self.ax1.set_title(f"Nav Separation = {self.aisNavSeparation:g}"+f" | Crossover = {self.aisNavCrossover:g}")
        #set y-limit for NaV density axis
        self.ax1.set_ylim([0.0, self.ais_gNavTotal*1.05])
        # Show Figure 1
        self.fig1.show()

        # Close Figure 1
        # plt.close(fig1)

    def get_vpeak(self):
        data = {}
        for sec in self.all:
        # for sec in [self.soma, self.axon0, self.axon1]:
            multiplier = 1.0
            if 'dend' in str(sec) or 'apic' in str(sec):
                multiplier = -1.0
            # if 'apic[36]' in str(sec):
            data[str(sec)] = [[multiplier*h.distance(self.soma(0.5), seg) for seg in sec], [seg.peak.vpeak for seg in sec]]
        self.vpeak_data = data

    def get_tpeak(self):
        data = {}
        for sec in self.all:
        # for sec in [self.soma, self.axon0, self.axon1]:
            multiplier = 1.0
            if 'dend' in str(sec) or 'apic' in str(sec):
                multiplier = -1.0
            # if 'apic[36]' in str(sec):
            data[str(sec)] = [[multiplier*h.distance(self.soma(0.5), seg) for seg in sec], [seg.peak.tpeak for seg in sec]]
        self.tpeak_data = data


    def plot_peak_everywhere(self, amplitude=None, separation=False, crossover=False, deltaRightShift=False, plot_tpeak=False, filename=None, legend1_location=None, legend2_location=None):
        import matplotlib
        import matplotlib.pyplot as plt
        from matplotlib import rc
        import tikzplotlib 
        ### global formatting
        label_size = 20
        legend_size = 17
        line_width = 3
        matplotlib.rcParams['xtick.labelsize'] = label_size
        matplotlib.rcParams['ytick.labelsize'] = label_size
        matplotlib.rcParams['axes.labelsize'] = label_size*1.3
        matplotlib.rcParams['legend.fontsize'] = legend_size 
        matplotlib.rcParams['legend.title_fontsize'] = legend_size 
        matplotlib.rcParams['lines.linewidth'] = line_width
        ###Create and Close Figure 1 and Figure 2, to clear any pre-existing figures: 
        self.fig1, self.ax1 = plt.subplots(figsize=(9,6))
        plt.close(self.fig1)
        if plot_tpeak==True:
            self.fig2, self.ax2 = plt.subplots(figsize=(9,6))
            plt.close(self.fig2)

        ### Apply the custom format function to the x and y axes
        from matplotlib.ticker import FuncFormatter
        ### Define a custom function for formatting
        def format_ticks(x, _):
            return f"{x:g}"
             # return f"{x:.1f}"
        self.formatter = FuncFormatter(format_ticks)
        self.fig1, self.ax1 = plt.subplots(figsize=(9,6))
        self.ax1.tick_params(left=True, right=True)
        self.ax1.xaxis.set_major_formatter(self.formatter)
        self.ax1.yaxis.set_major_formatter(self.formatter)
        if plot_tpeak==True:
            self.fig2, self.ax2 = plt.subplots(figsize=(9,6))
            self.ax2.tick_params(left=True, right=True)
            self.ax2.xaxis.set_major_formatter(self.formatter)
            self.ax2.yaxis.set_major_formatter(self.formatter)

        if separation or crossover or deltaRightShift:
            self.set_aisNavProfile(separation=separation, crossover=crossover, deltaRightShift=deltaRightShift)
        if amplitude is None:
            amplitude=8.2
        self.stim(amplitude=amplitude)
        h.finitialize(v_init)
        self.reset_peak()
        h.continuerun(t_stop)
        if filename is None:
            self.figname1 = 'vpeak_'+str(self.ic.amp)+'_'+str(self.ic.dur)+'_'+str(self.aisNavSeparation)+'_'+str(self.aisNavCrossover)+'_'+str(self.ic.get_segment())
            self.figname2 = 'tpeak_'+str(self.ic.amp)+'_'+str(self.ic.dur)+'_'+str(self.aisNavSeparation)+'_'+str(self.aisNavCrossover)+'_'+str(self.ic.get_segment())
        else:
            self.figname1 = filename+'_vpeak'
            self.figname2 = filename+'_tpeak'

        self.get_vpeak()
        print_args(self.cell.apic[37].nseg)
        self.ax1.plot(-1.0*h.distance(self.soma(0.5), self.cell.apic[37](0.6)), self.cell.apic[37](0.6).peak.vpeak, linewidth=0.0, linestyle='', marker='d', markersize=10.0, color='black', alpha=0.7, label='record')#str(self.cell.apic[37])[-8:])
        self.ax1.plot(-1.0*h.distance(self.soma(0.5), self.apic(1.0)), self.apic(0.99).peak.vpeak, linewidth=0.0, linestyle='', marker='p', markersize=10.0, color='black', alpha=0.7, label='bifurcation') #, label=str(self.apic)[-8:])
        for sec, segdata in self.vpeak_data.items():
            # if 'apic[36]' in str(sec): 
            #     self.ax1.plot(segdata[0], segdata[1], linewidth=0.0, linestyle='', marker='s', alpha=0.7, label=str(sec))
            # elif 'apic' in str(sec) and max(segdata[1])>23.0 and max(segdata[0])<-500.0:
            #     self.ax1.plot(segdata[0], segdata[1], linewidth=0.0, linestyle='', marker='p', alpha=0.7, label=str(sec))
            # elif 'apic' in str(sec) or 'dend' in str(sec):
            #     self.ax1.plot(segdata[0], segdata[1], linewidth=0.0, linestyle='', marker='o', color='blue', alpha=0.3)
            # else:
            #     self.ax1.plot(segdata[0], segdata[1], linewidth=0.0, linestyle='', marker='o', alpha=0.3, label=str(sec))
            if 'axon' in str(sec):
                self.ax1.plot(segdata[0], segdata[1], linewidth=0.0, linestyle='', marker='o', alpha=0.3, label='AIS_'+str(sec)[-2:-1])
            else:    
                self.ax1.plot(segdata[0], segdata[1], linewidth=0.0, linestyle='', marker='o', alpha=0.3)

        

        # create legend
        if legend1_location is not None:
            self.ax1.legend(loc=legend1_location)
        else:
            self.ax1.legend(loc='upper left')
        self.ax1.set_ylabel(r'$V_{peak}\ [mV]$', fontsize=27)
        self.ax1.set_xlabel(r'Distance to Soma $s\ [\mu m]$', fontsize=27)
        # Show Figure 1
        self.fig1.subplots_adjust(top=0.98, bottom=0.13)
        self.ax1.set_ylim(-86.0, 54.0)
        self.fig1.savefig(self.figname1+'.pdf')
        self.fig1.show()

        if plot_tpeak==True:
            self.get_tpeak()
            self.ax2.plot(-1.0*h.distance(self.soma(0.5), self.cell.apic[37](0.6)), self.cell.apic[37](0.6).peak.tpeak-self.ic.delay, linewidth=0.0, linestyle='', marker='d', markersize=10.0, color='black', alpha=0.7, label=str(self.cell.apic[37])[-8:])
            self.ax2.plot(-1.0*h.distance(self.soma(0.5), self.apic(1.0)), self.apic(0.99).peak.tpeak-self.ic.delay, linewidth=0.0, linestyle='', marker='p', markersize=10.0, color='black', alpha=0.7, label=str(self.apic)[-8:])
            # self.ax2.legend(loc='upper right')
            for sec, segdata in self.tpeak_data.items():
                # self.ax2.plot(segdata[0], np.array(segdata[1])-self.ic.delay, linewidth=0.0, linestyle='', marker='o', alpha=0.3, label=str(sec))
                if 'axon' in str(sec):
                    self.ax2.plot(segdata[0], np.array(segdata[1])-self.ic.delay, linewidth=0.0, linestyle='', marker='o', alpha=0.3, label='AIS_'+str(sec)[-2:-1])
                else:    
                    self.ax2.plot(segdata[0], np.array(segdata[1])-self.ic.delay, linewidth=0.0, linestyle='', marker='o', alpha=0.3)
            if legend2_location is not None:
                self.ax2.legend(loc=legend2_location)
            else:    
                self.ax2.legend(loc='upper right')
            self.ax2.set_ylabel(r'$t_{peak}\ [ms]$', fontsize=27)
            self.ax2.set_xlabel(r'Distance to Soma $s\ [\mu m]$', fontsize=27)
            # Show Figure 2
            self.fig2.subplots_adjust(top=0.98, bottom=0.13)
            self.fig2.savefig(self.figname2+'.pdf')
            self.fig2.show()


        # Close Figure 1
        # plt.close(fig1)

    def reset_peak(self):
        for sec in self.all:
            for seg in sec:
                seg.peak.vpeak = -1000.0
                seg.peak.tpeak = -1000.0


    def get_neuron_coordinates_and_colors(self):
        x, y, z, diameters, colors = [], [], [], [], []
        color_idx = 0

        for sec in self.all:
            for i in range(int(sec.n3d()) - 1):  # -1 to avoid an out-of-range error
                x1, y1, z1 = sec.x3d(i), sec.y3d(i), sec.z3d(i)
                x2, y2, z2 = sec.x3d(i + 1), sec.y3d(i + 1), sec.z3d(i + 1)
                diameter = (sec.diam3d(i) + sec.diam3d(i + 1)) / 2.0  # average diameter for this segment

                x.append(np.array([x1, x2]))
                y.append(np.array([y1, y2]))
                z.append(np.array([z1, z2]))
                diameters.append(diameter)
                colors.append(self.tropical_colors[color_idx])
            
            color_idx = (color_idx + 1) % len(self.tropical_colors)
        
        return x, y, z, diameters, colors

   


    def show3Dcell(self, colors=None, figsize=(10, 10), xlim=None, ylim=None, zlim=None):
        from mpl_toolkits.mplot3d.art3d import Line3DCollection
        x, y, z, diameters, colors = self.get_neuron_coordinates_and_colors()

        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')

        # Remove grid and axes
        ax.grid(False)
        ax.axis('off')
        ax.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))

        # Draw the neuron segments as pipes
        for xi, yi, zi, di, color in zip(x, y, z, diameters, colors):
            ax.add_collection3d(Line3DCollection([list(zip(xi, yi, zi))], linewidths=di, colors=color))

        # Set axis limits if provided
        if xlim:
            ax.set_xlim(xlim)
        else:
            ax.set_xlim(min([xi.min() for xi in x]), max([xi.max() for xi in x]))
        if ylim:
            ax.set_ylim(ylim)
        else:
            ax.set_ylim(min([yi.min() for yi in y]), max([yi.max() for yi in y]))
        if zlim:
            ax.set_zlim(zlim)
        else:
            ax.set_zlim(min([zi.min() for zi in z]), max([zi.max() for zi in z]))
        
        # Adjust aspect ratio and margins
        ax.set_box_aspect([ub - lb for lb, ub in (ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d())])
        plt.tight_layout()


        plt.show()




def export_dict_to_python_script(dict_to_export):
    # Use pprint to get a nicely formatted string of the dictionary
    formatted_dict_str = pprint.pformat(dict_to_export, width=80, sort_dicts=True)
    # Generate a file name with the current date and time
    datetime_now = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    file_name = f'data_dict_{datetime_now}.py'

    # Write to a file, including the variable assignment
    with open(file_name, 'w') as file:
        file.write(f"data = {formatted_dict_str}\n")
    print(f'Dictionary has been saved to {file_name}')

    


###*******************************************************************************
###*** Demo **********************************************************************
###*******************************************************************************
if __name__ == '__main__':
    ###
    ###Create cells
    ###
    cell1 = Pyramidal(1)
    cell2 = Pyramidal(2)
    ###
    ### Initialize/resume the simulation
    ###
    initialize(warm_up_and_save=True) 
    # initialize(warm_up_and_save=False)  
    for key, item in cell1.data_dict.items():
        print(key, item)



    ###Display cell1 in 3D matplotlib
    # cell1.show3Dcell()
    cell1.set_aisNavProfile(separation=0.0, crossover=0.5)
    cell1.plot_aisNavProfiles()
    # exit()


    # # cell1.stim(location=cell1.soma(0.5), duration=1.0)
    # cell1.stim(location=cell1.ais_segments[-1], duration=1.0)
    # cell1.threshCurve(separation_vals=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], crossover_vals=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    # thresh_data = copy.deepcopy(cell1.threshold_dict)


  
    




