"""
This script computes the threshold curves in Fig 2.
If you get an error saying 'ValueError: argument not a density mechanism name.', 
try running the script again.  If that doesn't work, you need to compile the MOD-file called CLSTmainen.mod manually.
"""

###Quick test:
# x_list = [0.0, 1.0]
# kappa_list = [0.2, 0.8]

###Full simulation
x_list = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
kappa_list = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

import os, sys
from mech_settings import mech_name
mech_filename = mech_name+'.mod'
os.system('nrnivmodl '+mech_filename)
print("--------NMODL mechanism compiled----------------------------------------"+'\n\n')

from run_settings import L_hill_extend
from run_settings import runs_dict, charges, timed_run, use_CVODE, N_ranvier, chloride, Hu_somaDend, Naked_axon_postAIS, cell_connected__FLAG, taper_hillock__FLAG, taper_dendrite__FLAG, FlipNavs_LeftRight, vThreshold, temperature, vLeftShift, AC, timeLS, CVODE_pre_equilm_dur, equilm_dur, resume_dur, total_dur, equilibration_timestep, stimulation_timestep, stimON, stimDur, stimOFF, fixed_step_dur
from mech_settings import mech_name, rescale_ranvier, rescale_soma, connect_cell, LeftShift_Sections, RightShift_Nav12, initialize_stimON, rescaleCLS_ΑΙSandNakedAxon, DeMyelinate, plot_CLS_vRS
from mech_settings import gamma as gamma0 
from measurement import return_spikeTimes, exit_times, cell_length, classify_grid_point, section_report, section_min_tpeak, ZeroCrossingTimes
from time import time
import datetime

cwd = os.getcwd()
sys.path.insert(1, cwd+'/morphology')   #<---insert at 1, 0 is the script path (or '' in REPL)
data_folder = cwd+'/data'

from neuron.units import ms, mV
from neuron import h, rxd#, gui
# rxd.nthread(4)

h.celsius = temperature
h.dt = equilibration_timestep



h.load_file('stdrun.hoc')
# cvode = h.CVode()
h.load_file('stdgui.hoc')
cvode = h.cvode
if use_CVODE == True:
    h.cvode_active(0)
    # cvode.use_long_double(1)
    cvode.atol(1.0e-6)



import numpy as np
from matplotlib import pyplot as plt
from BBplotting import BarlowRangeVarPlotter, plot_CLSparams

import mech_settings

if Hu_somaDend==True:
    from Hu_somaDend import makeSomatodendritic
else:
    from dendrite import makeDendrite
    from soma import makeSoma

from hillock import makeHillock
from ais import makeAIS, ais_tanh_distribution, diminish_AIS_Nav12, plot_ais_channels, plot_NaV_channels
from axon import makeAxon
if Naked_axon_postAIS==True:
    from naked import makeNaked
from tapering import conical_taper, cosine_taper
from setup_RxD import setup_RxD_active_outside2, setup_RxD_NaKCl, RxD_setConcentrations
from stimulate import squareCurrentPulse
from copy import deepcopy
import pandas as pd
from sorcery import dict_of, print_args, unpack_keys, unpack_attrs, maybe #<--- columns = dict_of(n_jobs, users, queues, priorities) gives a dictionary
from snoop import pp, spy #<--- import snoop ; snoop.install() #<---Then snoop, pp, and spy will be available in every file without needing to import them.


def count_segments():
    i=0
    for sec in h.allsec():
        for seg in sec:
            i+=1
    print('Total Number of Segments = ', i)

def apply_d_lambda_rule():
    ###Apply d_lambda rule:
    count_segments()
    h.xopen("fixnseg.hoc")
    soma(0.5).area()
    h.geom_nseg()
    count_segments()
    exit()



def getname(x, Vars=vars()):
    ### Returns variable name as string, eg: getname(a) == 'a'
    for k in Vars:
        if type(x) == type(Vars[k]):
            if x is Vars[k]:
                return k
    return None



def listify_dict_scalars(dict):
    for key in dict.keys():
        dict[key] = [dict[key]]

def append_dicts(d0_input, d1_input):
        from copy import deepcopy
        d0 = deepcopy(d0_input) ; l0 = len(d0)
        d1 = deepcopy(d1_input) ; l1 = len(d1)
        d0.update(d1)
        # print_args(l0, l1, len(d0))
        if len(d0) != (l0+l1):
            print('ERROR: DUPLICATE VARIABLE NAMES IN data_dicts (overwritten by .update()), CHECK VARIABLE NAMES IN measurement.py (Because of how measurement.py is used here, functions that return dictionaries must not use the same keys.)')
            exit()
        return deepcopy(d0)



#___________________________________________________________________________________________
###Build the cell
naked_list = []
connect_these = []
if Hu_somaDend == True:
    dend, soma, sd_list = makeSomatodendritic(h)
    # soma.nseg *= 3
    connect_these+=[soma]
    naked_list+=sd_list    
else:
    dend = makeDendrite(h)
    connect_these+=[dend]
    naked_list += [dend] 
    if taper_dendrite__FLAG == True:
        dend_taper = makeDendrite(h, name='dend_taper', L=11.0) 
        soma = makeSoma(h)
        cosine_taper(left_section=dend, tapering_section=dend_taper, right_section=soma)
        connect_these+=[dend_taper]
        naked_list += [dend_taper]
    else:
        soma = makeSoma(h)
    
    connect_these+=[soma]
    naked_list += [soma] 

hill = makeHillock(h) 
ais = makeAIS(h)
connect_these+=[hill, ais]
naked_list += [hill, ais]
if Naked_axon_postAIS==True:
    naked = makeNaked(h) 
    naked_list.append(naked)
    connect_these.append(naked)
axon_list, ranvier_list, myelin_list = makeAxon(h, N_ranvier=N_ranvier) 
naked_list += ranvier_list
connect_these += axon_list


if cell_connected__FLAG==True:
    full_cell = connect_cell(h, connect_these)
h.define_shape()

if taper_hillock__FLAG == True:
    # conical_taper(left_section=soma, tapering_section=hill, right_section=ais)
    cosine_taper(left_section=soma, tapering_section=hill, right_section=ais)

    # print_args(soma(1.0).diam)
    # print_args(hill(0.0).diam)

if L_hill_extend > 0.0:
    diameters = []
    positions = []
    L_hill_0 = hill.L
    L_hill_1 = L_hill_0+L_hill_extend
    initial_hill_seg_length = L_hill_0/hill.nseg
    New_nseg = int(hill.nseg*(L_hill_1/L_hill_0))
    for seg in hill:
        positions += [seg.x*L_hill_0]
        diameters += [seg.diam]
        
    
    print_args(positions)
    print_args(diameters)
    hill.diam = diameters[-1]
    hill.nseg = New_nseg
    hill.L = L_hill_1
    for i in range(len(positions)):
        for seg in reversed(hill):
            if seg.x*hill.L <= reversed(positions[i]):
                seg.diam = reversed(diameters[i])


        
    print_args(positions)
    print_args(diameters)

    exit()

def change_Ra():
    from run_settings import Ra, Ra_axon
    for sec in h.allsec():
        sec.Ra = Ra
    for sec in axon_list:
        sec.Ra = Ra_axon
change_Ra()

def change_cm():
    from run_settings import cm, cm_myelin
    for sec in h.allsec():
        for seg in sec.allseg():
            seg.cm = cm
    for sec in myelin_list:
        for seg in sec.allseg():
            seg.cm = cm_myelin
change_cm()





from run_settings import stimLoc, Hu_somaDend
if stimLoc=='ais':
    stimLoc=ais(0.9)
elif stimLoc=='soma':
    stimLoc=soma(0.3)

else:
    print('ERROR: Need to select stimLoc in run_settings.py')
    exit()



from run_settings import stimON, stimDur, charge
initialize_stimON(h) 
amplitude = charge/stimDur
ic, stimLoc = squareCurrentPulse(stimLoc=stimLoc, delay=stimON, amplitude=amplitude, duration=stimDur)



#___________________________________________________________________________________________
###RxD
rxd_sections = h.allsec() 
if chloride==True:
    ki, ko, nai, nao, cli, clo, cyt, ecs, na, k, cl = setup_RxD_NaKCl(h, rxd, cvode, sections=rxd_sections, diffusion_on=True, d_Na=1.79, d_K=2.62, d_Cl=2.72) #<---contains h.finitialize() [DO NOT CALL h.finitalize() AFTER THIS POINT! WRECKS THE CONCENTRATION INTIALIZATION!]
else:
    ki, ko, nai, nao, cyt, ecs, na, k = setup_RxD_active_outside2(h, rxd, cvode, sections=rxd_sections, diffusion_on=True, d_Na=0.6, d_K=1.0) #<---contains h.finitialize() [DO NOT CALL h.finitalize() AFTER THIS POINT! WRECKS THE CONCENTRATION INTIALIZATION!]





#___________________________________________________________________________________________
###Setup Recording Vectors)
from recording_vectors_setup import seg_timeseries, setup_recording_vectors, seg_mech_timeseries
select_node =  ranvier_list[-1] #ranvier_list[N_ranvier-8]
data_dict = setup_recording_vectors(h, [soma(0.5), ais(0.5), select_node(0.5)])







class PyData:
    fig = plt.figure(figsize=(8, 6))
    p = None
    segments = []
    data_objects = []  #<---These lists will be shared by all instances of the Class "PyData". Updating in one updates them all.
    PyData_dict = {}
    import matplotlib.pyplot as plt
    def __init__(self, segment):
        """If True, infinite values for vpeak and t_vpeak et cetera are returned. This is meant to prevent unreliable measurements from being recorded."""
        self.return_null_data = False

        self.sensitivity=1.0*mV
        self.postStimCutOff=100.0*ms
        self.dvdt_thresh = 30.0*mV

        self.seg = segment
        self.segname = str(self.seg)
        self.t_vec = h.Vector().record(h._ref_t)
        self.v_vec = h.Vector().record(self.seg._ref_v)
        self.data_objects.append(self)


        self.sec = self.seg.sec
        self.all_xyz = []
        self.xyz = []
        self.seg_dist = h.distance(self.sec(0), self.seg)
        dist = 0.0
        dist_low = 0.0
        dist_high = 0.0
        for i in np.arange(self.sec.n3d()):
            xyz_i = np.array([self.sec.x3d(i), self.sec.y3d(i), self.sec.z3d(i)])
            self.all_xyz.append(xyz_i)
            if self.sec.arc3d(i) < dist:
                print('Error: segment position calculator assumes sec.n3d(i) starts at sec(0)')
                exit()
            dist = self.sec.arc3d(i)
            if dist < self.seg_dist:
                dist_low = dist
                xyz_low = xyz_i
            elif dist == self.seg_dist:
                self.xyz = xyz_i
            elif dist_high == 0.0:
                dist_high = dist
                xyz_high = xyz_i
        if len(self.xyz)==0:
            self.xyz = (xyz_low + xyz_high)/2.0
        self.BlenderString = 'bpy.ops.surface.primitive_nurbs_surface_sphere_add(radius=7, enter_editmode=False, location=('+str(self.xyz[0])+','+str(self.xyz[1])+','+str(self.xyz[2])+'))' +'\n'+'bpy.context.object.name = \"'+self.segname+'\"'



        if segment==stimLoc:
            print('WARNING: YOU ARE RECORDING FROM THE EXACT SAME LOCATION (same segment) WHERE CURRENT IS INJECTED. THIS MAY CAUSE ERRORS IN SPIKE/PEAK ANALYSIS.')
        elif segment.sec == stimLoc.sec:
            print('WARNING: YOU ARE RECORDING FROM THE SAME SECTION WHERE CURRENT IS INJECTED.')
            print('THIS MAY CAUSE ERRORS IN SPIKE/PEAK ANALYSIS: ensure that the recording SEGMENT (subcompartment of this section) is not right next to the stimulation segment.')
    


    def add_segment(self, segment):
        self.segments.append(segment)



    def process_postStim(self, plotPeak=False):
        self.plotPeak = plotPeak

        self.t_full = self.t_vec.as_numpy()
        self.v_full = self.v_vec.as_numpy()

        self.t = self.t_vec.as_numpy()
        self.v = self.v_vec.as_numpy()
        select = np.where((self.t>=stimON) & (self.t<=stimOFF+self.postStimCutOff))
        self.t = self.t[select]
        self.v = self.v[select]

        self.vrest = self.v[0]
        self.vmin = min(self.v)
        self.vmax = max(self.v)


        self.vpeak_save = max(self.v)
        ###First approximation
        self.vpeak = max(self.v)
        self.t_vpeak = self.t[np.where(self.v==self.vpeak)][0]


        self.depolarization=(self.vpeak-self.vrest)
        self.t_depolarization=self.t[np.where(self.v==self.vpeak)][0]
        print(self.segname+': depolarization='+str(self.depolarization))#+', '+'t_depolarization='+str(self.t_depolarization))
        if self.depolarization>=self.sensitivity:
            """First and second derivatives of voltage."""
            # self.dv = np.gradient(self.v)
            # self.d2v = np.gradient(self.dv)
            # self.dvdt = self.dv
            # self.d2vdt2 = self.d2v
            
            self.dv = np.gradient(self.v)
            # self.dv = self.v[1:] - self.v[:-1]
            # self.d2vdt2 = np.gradient(self.dvdt, self.t)


            if self.t_vpeak<stimOFF:
                print('ERROR, PULSE DURATION IS TOO LONG (stimOFF>=t_vpeak). Analysis assumes vpeak occurs after the pulse. Null data will be returned.')
                self.return_null_data=True
            elif np.absolute(self.t_vpeak-stimOFF)<=0.2*ms:
                reselect = np.where((self.t>=stimOFF+0.5*ms))
                if (max(self.v[reselect]) > min(self.v[np.where(self.t <= stimOFF+4.0*ms)]) ):
                    ###First approximation
                    if max(self.v[reselect]) != self.v[reselect][0]:
                        self.vpeak = max(self.v[reselect])
                        self.t_vpeak = self.t[np.where(self.v==self.vpeak)][0]
                    if np.absolute(self.t_vpeak-stimOFF)<=0.2*ms:
                            print('vpeak too close to stimOFF: t_vpeak will not be interpolated. Null data will be returned.')
                            self.return_null_data=True

                else:
                    print('vpeak too close to stimOFF: t_vpeak will not be interpolated. Null data will be returned.')
                    self.return_null_data=True

            else:
                from scipy.signal import argrelmax, argrelmin

                ###interpolate t_vpeak using dv:
                from scipy import interpolate
                f = interpolate.interp1d(self.t, self.dv)
                tnew = np.linspace(self.t_vpeak-0.1, self.t_vpeak+0.1, 299)
                dvnew = f(tnew)
                f = interpolate.interp1d(dvnew, tnew)
                dvnew = 0.0
                self.tnew = float(f(dvnew))
                if np.absolute(self.tnew-self.t_vpeak)>=1.0:
                    print('ERROR: dv interpolation of peak has failed.')
                    self.return_null_data=True
            
                else:
                    self.t_vpeak=self.tnew

            ###Check for multiple spikes
            select_postPeak = np.where(self.t>self.t_vpeak)
            self.t_postPeak = self.t[select_postPeak]
            self.v_postPeak = self.v[select_postPeak]
            self.dt_postPeak = np.gradient(self.t_postPeak)
            self.dv_postPeak = np.gradient(self.v_postPeak)
            if max(self.dv_postPeak - self.dvdt_thresh*self.dt_postPeak)>0.0:
                print('ERROR: Multiple spikes.')
                print('SPIKE ANALYSIS NOT YET EQUIPPED FOR MULTIPLE SPIKES. Null data will be returned.')
                self.return_null_data=True
        else:
            print('Post-stimulation voltage change (depolarization= '+str(self.depolarization)+'mV) is too small. (Sensitivity set to '+str(self.sensitivity)+'mV.)')
            print('THUS no voltage peak or spike data will be saved at this segment (i.e. this location).')
            self.return_null_data=True

        if self.plotPeak==True:
            self.plot()
        

        if self.return_null_data==True:
            self.t_vpeak=np.inf
            self.vpeak=np.inf
        # print(self.segname)
        # print_args(self.vpeak, self.t_vpeak)
        self.PyData_dict[self.segname+'vpeak'] = self.vpeak
        self.PyData_dict[self.segname+'t_vpeak'] = self.t_vpeak
        self.PyData_dict[self.segname+'depolarization'] = self.depolarization
        self.PyData_dict[self.segname+'t_depolarization'] = self.t_depolarization


    def plot(self, linestyle='-'):
        # import matplotlib.pyplot as plt

        self.p = self.plt.plot(self.t_full, self.v_full, label=self.segname, linestyle=linestyle)
        # self.p = self.plt.plot(self.t, self.v, label=self.segname, linestyle=linestyle)
        color = self.p[0].get_color()
        if self.return_null_data==False:
            # self.plt.plot(self.t_vpeak, self.vpeak, marker=r'$'+str(self.seg.x)+r'$', color=color, markersize=25.0, alpha=0.7)
            self.plt.plot(self.t_vpeak, self.vpeak, marker='*', color=color, markersize=17.0, alpha=0.7)
        else:
            print('peak discarded for '+self.segname)
        self.plt.legend()

    @classmethod
    def process_all(cls, plotPeak=False):
        for obj in cls.data_objects:
            obj.process_postStim(plotPeak)


    @classmethod
    def print_BlenderScript(cls):
        print('Blender script: ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■')
        print('import bpy'+'\n')
        for obj in cls.data_objects:
            print(obj.BlenderString)
        print('■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■')


class ConductionVelocity(PyData):
    def __init__(self, seg0, seg1):
        self.seg0 = seg0
        self.seg1 = seg1

        self.data0 = PyData(seg0)
        self.data1 = PyData(seg1)


    def get_speed(self):
        self.data0.process_postStim()
        self.data1.process_postStim()

        self.length = h.distance(self.seg0, self.seg1)/1000.0   #<---Converts distance to mm
        self.delta_t = np.absolute(self.data1.t_depolarization - self.data0.t_depolarization)
        self.speed = self.length/self.delta_t


        print('delta_t (ms) = ', self.delta_t)
        print('length (mm) = ', self.length)
        print('speed (m/s) = ', self.speed)



somaData = PyData(soma(0.5))
PyData.print_BlenderScript()

#__________________________________________________________________________________________________________________________________________________________
#■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
#■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
#■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
#■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
#■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
#■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
#■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
class StartUp():
    """docstring for ClassName"""
            
    from run_settings import stimON, stimDur, amplitude, charge, FlipNavs_LeftRight, Naked_axon_postAIS, temperature, vLeftShift, AC
    print_args(stimON)
    
    global ic, stimLoc, RunUp_duration
    global cvode

        

    RunUp_duration = 30.0*ms
    ic, stimLoc = squareCurrentPulse(stimLoc=stimLoc, delay=RunUp_duration+10.0*ms, amplitude=0.0, duration=stimDur)

    def __init__(self, gamma=None, deviation=None, CrossOverPosition=None, fraction=None, vRS0=None, vRS03 = None):

        if gamma is None:
            from run_settings import gamma as rs_gamma
            self.gamma=rs_gamma
        else:
            self.gamma=gamma


        if deviation is None:
            from run_settings import deviation as rs_deviation
            self.deviation=rs_deviation
        else:
            self.deviation=deviation


        if CrossOverPosition is None:
            from run_settings import CrossOverPosition as rs_CrossOverPosition
            self.CrossOverPosition=rs_CrossOverPosition
        else:
            self.CrossOverPosition=CrossOverPosition


        if fraction is None:
            from run_settings import fraction as rs_fraction
            self.fraction=rs_fraction
        else:
            self.fraction=fraction


        if vRS0 is None:
            from run_settings import vRS0 as rs_vRS0
            self.vRS0=rs_vRS0
        else:
            self.vRS0=vRS0

        if vRS03 is None:
            from run_settings import vRS03 as rs_vRS03
            self.vRS03=rs_vRS03
        else:
            self.vRS03=vRS03


        self.now = datetime.datetime.now()
        self.date_str = self.now.strftime('%Y-%m-%d')
        self.time_str = self.now.strftime('%H-%M-%S')
        


        from mech_settings import rescale_soma
        self.rescale_soma = rescale_soma
        self.rescale_ais = None 
        self.gammaAIS = None


        for sec in h.allsec():
            seg_len = sec.L/sec.nseg
            for seg in sec:
                setattr(seg, 'seg_len_'+mech_name, seg_len)


        self.reset_channels_and_ions()
        self.selectiveRightShift_NaV12(section=ais, vRS03=self.vRS03)

        
        from time import time ; runStart = time()
        if use_CVODE == True:
            timed_run(h, duration=1.0*ms)
            h.cvode_active(1)
        else:
            timed_run(h, duration=1.0*ms)
            h.cvode_active(1)
            timed_run(h, duration=2.0*ms)
            h.cvode_active(0)
        if chloride==True:
            RxD_setConcentrations(ions=[nai, nao, ki, ko, cli, clo], cvode=cvode)
        else:
            RxD_setConcentrations(ions=[nai, nao, ki, ko], cvode=cvode)

        timed_run(h, duration=RunUp_duration)
        print_args(ais.L)
        print_args(ais.diam)

        
        self.dt_animate = 1.0e-3*ms
        self.chloride = chloride
        self.propFunc = None
        self.Always_reset_channels_and_ions = True


        self.h = h
        self.sd_list = sd_list
        self.myelin_list = myelin_list
        self.ranvier_list = ranvier_list
        self.soma = soma
        self.hill = hill
        self.ais = ais
        self.naked = naked
        self.myelin = myelin_list[0]
        self.ranvier = ranvier_list[0]
        self.axon_list = axon_list
        self.send_dict = {}
        self.Ra = self.soma.Ra
        self.T = h.celsius
        self.Cm_myelin = self.myelin.cm
        self.Cm_soma = self.soma.cm
        

        self.vThreshold = -63.0 #-63.2 #-63.0 #-60.0
        self.vThreshold_data = []
        self.vpeak_SDtips = None
        self.stimLoc_data = []

        self.vpeak_data = []
        # self.gamma_data = []
        self.gammaAIS_data = []
        self.deviation_data = []
        self.fraction_data = []
        self.CrossOverPosition_data = []
        self.vRS03_data = []
        self.DeltaVrs_data = []
        

        self.delta_min = 1.0e-4
        self.Qthresh = None
        self.Qthresh_data = []
        self.iThresh_backprop = None
        self.Qthresh_backProp = None
        self.Qthresh_backProp_data = []
        self.iThresh_forwardProp = None
        self.Qthresh_forwardProp = None
        self.Qthresh_forwardProp_data = []
        self.amplitude = 5.0
        self.maximumAmplitude=10.0
        self.backProp_criterion = self.min_vpeak_apicalSDtips #self.max_vpeak_apicalSDtips #self.min_vpeak_apicalSDtips #self.mean_vpeak_SDtips #self.get_mainApicaltip_vpeak #self.get_mainApicalMiddle_vpeak 
        print('**************************************')
        print('**************************************')
        print('**************************************')
        print_args(self.backProp_criterion.__name__)
        print('**************************************')
        print('**************************************')
        print('**************************************')
        self.ic = ic
        self.stimLoc = self.ic.get_segment()
        self.delay = ic.delay
        self.amp = self.ic.amp
        self.concentrationReinit=True
        self.preStim_delay = 1000.0
        self.postStim_run = 1000.0
        
        self.RunUp_duration = RunUp_duration
        

        

        self.pandas_data = {}
        self.pandas_df = pd.DataFrame.from_dict(self.pandas_data)
        self.df_filename = 'pandas_data'+'_Lais='+str(self.ais.L)+'_'+self.date_str+self.time_str+'.csv'
        self.pandas_df.to_csv(self.df_filename)

    # def getConcentrations(self, sec):
        self.soma_nai = self.soma(0.5).na_ion.nai
        self.soma_nao = self.soma(0.5).na_ion.nao
        self.soma_ki = self.soma(0.5).k_ion.ki
        self.soma_ko = self.soma(0.5).k_ion.ko
        self.soma_cli = self.soma(0.5).cl_ion.cli
        self.soma_clo = self.soma(0.5).cl_ion.clo

        self.hill_nai = self.hill(0.5).na_ion.nai
        self.hill_nao = self.hill(0.5).na_ion.nao
        self.hill_ki = self.hill(0.5).k_ion.ki
        self.hill_ko = self.hill(0.5).k_ion.ko
        self.hill_cli = self.hill(0.5).cl_ion.cli
        self.hill_clo = self.hill(0.5).cl_ion.clo

        self.ais_nai = self.ais(0.5).na_ion.nai
        self.ais_nao = self.ais(0.5).na_ion.nao
        self.ais_ki = self.ais(0.5).k_ion.ki
        self.ais_ko = self.ais(0.5).k_ion.ko
        self.ais_cli = self.ais(0.5).cl_ion.cli
        self.ais_clo = self.ais(0.5).cl_ion.clo

        self.naked_nai = self.naked(0.5).na_ion.nai
        self.naked_nao = self.naked(0.5).na_ion.nao
        self.naked_ki = self.naked(0.5).k_ion.ki
        self.naked_ko = self.naked(0.5).k_ion.ko
        self.naked_cli = self.naked(0.5).cl_ion.cli
        self.naked_clo = self.naked(0.5).cl_ion.clo

        self.myelin_nai = self.myelin(0.5).na_ion.nai
        self.myelin_nao = self.myelin(0.5).na_ion.nao
        self.myelin_ki = self.myelin(0.5).k_ion.ki
        self.myelin_ko = self.myelin(0.5).k_ion.ko
        self.myelin_cli = self.myelin(0.5).cl_ion.cli
        self.myelin_clo = self.myelin(0.5).cl_ion.clo

        self.ranvier_nai = self.ranvier(0.5).na_ion.nai
        self.ranvier_nao = self.ranvier(0.5).na_ion.nao
        self.ranvier_ki = self.ranvier(0.5).k_ion.ki
        self.ranvier_ko = self.ranvier(0.5).k_ion.ko
        self.ranvier_cli = self.ranvier(0.5).cl_ion.cli
        self.ranvier_clo = self.ranvier(0.5).cl_ion.clo






    def reset_channels_and_ions(self, gamma=None, gammaAIS=None, deviation=None, CrossOverPosition=None, fraction=None, vRS0=None):
        if gamma is not None:
            self.gamma = gamma
        if gammaAIS is not None:
            self.gammaAIS = gammaAIS
        if deviation is not None:
            self.deviation = deviation
        if CrossOverPosition is not None:
            self.CrossOverPosition = CrossOverPosition
        if fraction is not None:
            self.fraction = fraction
        if vRS0 is not None:
            self.vRS0 = vRS0

        RightShift_Nav12(h, mech_name=mech_name, vRS0=self.vRS0)
        sections = [ranvier_list[-5]] #sd_list + [ranvier_list[-5]] #sd_list + [ais] + ranvier_list #[ranvier_list[-3]] #list(h.allsec())#[select_node] #ranvier_list #naked_list #[select_node]
        LeftShift_Sections(sections=sections, mech_name=mech_name, vLeftShift=vLeftShift, AC=AC, timeLS=timeLS)
        vLS = seg_mech_timeseries(h, recording_segment=sections[0](0.5), varname='vLS')
        ### ais channel distributions: MUST FOLLOW h.finitialize()
        if Naked_axon_postAIS==True:
            gnabarTotal = rescaleCLS_ΑΙSandNakedAxon([hill, ais, naked], self.gamma)
        else:
            gnabarTotal = rescaleCLS_ΑΙSandNakedAxon([hill, ais], self.gamma)
        gnabarTotal = ais_tanh_distribution(ais, gnabarTotal=gnabarTotal, flatness=10.0, deviation=self.deviation, FlipNavs_LeftRight=FlipNavs_LeftRight, do_nothing=False, CrossOverPosition=self.CrossOverPosition)
        if self.fraction!=1.0:
            gnabar2 = diminish_AIS_Nav12(ais, fraction=self.fraction)



        from mech_settings import linear_Nav_matcher, all_one_Nav
        linear_Nav_matcher(hill, ais, gamma=self.gamma)

        



        if (self.gammaAIS != None) and (self.gammaAIS != self.gamma):
            # print('rescaling AIS NaVs only by gammaAIS=', self.gammaAIS, ', whereas gamma=', self.gamma)
            self.rescale_AIS_voltage_gated_ONLY(self.gammaAIS)



    def selectiveRightShift_NaV12(self, section, vRS03=None, add_MODfilename=True):
        """Adjusts the local right-shift in section, selectively to any combination of minf2, mtau2, hinf2, htau2, depending on the disableRightShift_minf2 (etc.) flags from run_settings.py"""
        if vRS03 is not None:
            self.vRS03 = vRS03

        """Only has an effect if the flags from run_settings.py are set to True"""
        from run_settings import disableRightShift_minf2__inAIS, disableRightShift_mtau2__inAIS, disableRightShift_hinf2__inAIS, disableRightShift_htau2__inAIS
        from mech_settings import mech_name
        import numpy as np
        from neuron import h

        for seg in section.allseg():
            if disableRightShift_minf2__inAIS:
                setattr(seg, 'disableRightShift_minf2'+'_'+mech_name, 1)
            if disableRightShift_mtau2__inAIS:
                setattr(seg, 'disableRightShift_mtau2'+'_'+mech_name, 1)
            if disableRightShift_hinf2__inAIS:
                setattr(seg, 'disableRightShift_hinf2'+'_'+mech_name, 1)
            if disableRightShift_htau2__inAIS:
                setattr(seg, 'disableRightShift_htau2'+'_'+mech_name, 1)
        

        MODfilename = r''
        if add_MODfilename:
            MODfilename = '_'+'CLSTmainen'
        for seg in section.allseg():
                setattr(seg, 'vRS03_'+mech_name, self.vRS03)
        print('⁍vRS03 = ', self.vRS03)

    def rescale_AIS_voltage_gated_ONLY(self, gammaAIS=None, add_MODfilename=True):
        MODfilename = r''
        if add_MODfilename:
            MODfilename = '_'+'CLSTmainen'

        if gammaAIS is None:
            self.gammaAIS = self.gamma     
        else:
            self.gammaAIS = gammaAIS

            self.rescale_ais = self.gammaAIS/self.gamma
            self.rescale_naked = (27.0/34.0)*self.rescale_ais
            for seg in self.ais:
                gna16 = self.rescale_ais*getattr(seg, 'gnabar'+MODfilename)
                gna12 =  self.rescale_ais*getattr(seg, 'gnabar2'+MODfilename)
                setattr(seg, 'gnabar'+MODfilename, gna16)
                setattr(seg, 'gnabar2'+MODfilename, gna12)

                gk = self.rescale_ais*getattr(seg, 'gkbar'+MODfilename)
                setattr(seg, 'gkbar'+MODfilename, gk)

            for seg in self.naked:
                gna16 = self.rescale_naked*getattr(seg, 'gnabar'+MODfilename)
                gna12 =  self.rescale_naked*getattr(seg, 'gnabar2'+MODfilename)
                setattr(seg, 'gnabar'+MODfilename, gna16)
                setattr(seg, 'gnabar2'+MODfilename, gna12)

                gk = self.rescale_naked*getattr(seg, 'gkbar'+MODfilename)
                setattr(seg, 'gkbar'+MODfilename, gk)

            gna16_left = getattr(soma(1.0), 'gnabar'+MODfilename)
            gna16_right = getattr(ais(0.0), 'gnabar'+MODfilename)
            gna16_vals = np.linspace(gna16_left, gna16_right, self.hill.nseg)

            gna12_left = getattr(soma(1.0), 'gnabar2'+MODfilename)
            gna12_right = getattr(ais(0.0), 'gnabar2'+MODfilename)
            gna12_vals = np.linspace(gna12_left, gna12_right, self.hill.nseg)

            gk_left = getattr(soma(1.0), 'gkbar'+MODfilename)
            gk_right = getattr(ais(0.0), 'gkbar'+MODfilename)
            gk_vals = np.linspace(gk_left, gk_right, self.hill.nseg)
            i=0
            for seg in self.hill:
                setattr(seg, 'gnabar'+MODfilename, gna16_vals[i])
                setattr(seg, 'gnabar2'+MODfilename, gna12_vals[i])
                
                setattr(seg, 'gkbar'+MODfilename, gk_vals[i])
                i+=1

        gNavTotal_AIS = getattr(self.ais(0.7), 'gnabar'+MODfilename) + getattr(self.ais(0.7), 'gnabar2'+MODfilename)
        print_args(gNavTotal_AIS)



    def get_segVar(self, varname, segment, add_MODfilename=True):
        MODfilename = r''
        if add_MODfilename:
            MODfilename = '_'+'CLSTmainen'
        return getattr(segment, varname+MODfilename) 

        

    def get_rangeVar(self, varname, sections, home_segment=soma(0), add_MODfilename=True):
        MODfilename = r''
        if add_MODfilename:
            MODfilename = '_'+'CLSTmainen'
        data = []
        position = []
        for sec in sections:
            # home_segment = sec.subtree()[-1](1.0)
            for seg in sec.allseg():
                # position.append(h.distance(sections[0](0), seg))
                position.append(self.h.distance(home_segment, seg))
                data.append(getattr(seg, varname+MODfilename))
        return [position, data]


    def get_endVar(self, varname, sections, home_segment, add_MODfilename=True):
        MODfilename = r''
        if add_MODfilename:
            MODfilename = '_'+'CLSTmainen'
        data = []
        position = []
        for sec in sections:
            if sec.subtree() == [sec]:
                seg = sec(1.0)
                position.append(self.h.distance(home_segment, seg))
                data.append(getattr(seg, varname+MODfilename))
        return [position, data]    


    def get_tip_segments(self, sections=None, add_MODfilename=True):
        if sections is None:
            sections=self.sd_list
        MODfilename = r''
        if add_MODfilename:
            MODfilename = '_'+'CLSTmainen'

        tip_segments = []
        for sec in sections:
            if sec.subtree() == [sec]:
                seg = sec(1.0)
                tip_segments.append(seg)
        return tip_segments 


    def reset_segVar(self, varname, segments, reset_value=-1.0e6, add_MODfilename=True):
        MODfilename = r''
        if add_MODfilename:
            MODfilename = '_'+'CLSTmainen'
        for seg in segments:
            setattr(seg, varname+MODfilename, reset_value)
        return 0

    def reset_tip_vpeak(self, varname='vpeak'):
        self.reset_segVar('vpeak', self.get_tip_segments(), reset_value=-1000.0)
        self.reset_segVar('tpeak', self.get_tip_segments(), reset_value=-1.0)
        self.reset_segVar('tup', self.get_tip_segments(), reset_value=-1.0)
        self.reset_segVar('tdown', self.get_tip_segments(), reset_value=-1.0)
        self.reset_segVar('tspike', self.get_tip_segments(), reset_value=-1.0)
        self.reset_segVar('spike', self.get_tip_segments(), reset_value=0.0)
        # self.reset_segVar('stimON', self.get_tip_segments(), reset_value=self.ic.delay)


    def reset_vpeak(self, segments=None, varname='vpeak'):
        if segments is None:
            segments = []
            for sec in h.allsec():
                segments+=[seg for seg in sec.allseg()]
        self.reset_segVar('vpeak', segments, reset_value=-1000.0)
        self.reset_segVar('tpeak', segments, reset_value=-1.0)
        self.reset_segVar('tup', segments, reset_value=-1.0)
        self.reset_segVar('tdown', segments, reset_value=-1.0)
        self.reset_segVar('tspike', segments, reset_value=-1.0)
        self.reset_segVar('spike', segments, reset_value=0.0)
        # self.reset_segVar('stimON', segments, reset_value=self.ic.delay)


    def get_SDtips_vpeak(self):
        data = self.get_endVar('vpeak', self.sd_list, home_segment=self.soma(0.5))
        return data

    def get_apical_SDtips_vpeak(self):
        data = self.get_endVar('vpeak', self.sd_list, home_segment=self.soma(0.5))
        vpeak_array = np.array(data[1])
        position_array = np.array(data[0])
        select_elements = np.where(position_array>=1000.0)
        return [position_array[select_elements], vpeak_array[select_elements]]

    def get_mainApicaltip_vpeak(self):
        data = self.get_rangeVar('vpeak', [self.h.dend11[22]], home_segment=self.soma(0.5))
        return data[1][-1]

    def get_mainApicalMiddle_vpeak(self):
        return self.get_segVar('vpeak', self.h.dend11[22](0.5))



    def mean_vpeak(self, sections=None):
        if sections is None:
            sections = self.h.allsec()
        data = self.get_rangeVar('vpeak', sections, home_segment=self.soma(0.5))
        return np.mean(data[1])

    def mean_vpeak_SDtips(self):
        data = self.get_endVar('vpeak', self.sd_list, home_segment=self.soma(0.5))
        return np.mean(data[1])

    def mean_vpeak_apicalSDtips(self):
        data = self.get_apical_SDtips_vpeak()
        return np.mean(data[1])

    def min_vpeak_apicalSDtips(self):
        data = self.get_apical_SDtips_vpeak()
        min_vpeak = np.min(data[1])
        print('min_vpeak_apicalSDtips = ', min_vpeak)
        return min_vpeak

    def max_vpeak_apicalSDtips(self):
        data = self.get_apical_SDtips_vpeak()
        max_vpeak = np.max(data[1])
        print_args(max_vpeak)
        return max_vpeak



    def stimulate(self):
        self.reset_tip_vpeak()
        self.ic.amp = self.amplitude
        self.ic.delay = self.h.t+self.delay #self.RunUp_duration*2.0
        self.h.cvode.re_init()
        self.h.continuerun(self.h.t+self.RunUp_duration*3.0)
        data_vpeak_SDtips = self.get_SDtips_vpeak()
        print(np.mean(data_vpeak_SDtips[1]))
        print(self.h.t)
        return data_vpeak_SDtips

    def concReinit(self):
        if self.chloride==True:
            RxD_setConcentrations(ions=[nai, nao, ki, ko, cli, clo], cvode=cvode)
        else:
            RxD_setConcentrations(ions=[nai, nao, ki, ko], cvode=cvode)

    def runStim(self):
        self.ic.amp = self.amplitude
        self.ic.delay = self.h.t+self.delay
        if self.concentrationReinit: self.concReinit()
        self.h.cvode.re_init()
        print_args(self.df_filename)
        print(' = pd.read_csv(\''+self.df_filename+'\')')
        print('+= [\'sqlite3'+self.df_filename+'\']')
        print_args(self.pandas_df)
        print('____________________________________________________')
        print_args(self.h.t)
        print_args(self.ic.delay)
        print_args(self.ic.amp)
        self.h.continuerun(self.ic.delay+self.ic.dur+self.postStim_run)

    def backPropCheck(self, amplitude=None, preStim_delay=None, postStim_run=None, vThreshold=None):
        if self.Always_reset_channels_and_ions:
            self.reset_channels_and_ions()
        if preStim_delay:
            self.preStim_delay = preStim_delay
        if postStim_run:
            self.postStim_run = postStim_run
            
        self.delay=self.preStim_delay
        if amplitude is None:
            self.amplitude=self.maximumAmplitude
        else:
            self.amplitude=amplitude
        if vThreshold is not None:
            self.vThreshold = vThreshold
        self.reset_vpeak()
        self.runStim()
        number = self.backProp_criterion()
        if number>=self.vThreshold:
            backProp=True
            self.vpeak_SDtips = self.backProp_criterion()
        else:
            backProp=False
        print_args(number)
        print_args(backProp)
        print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')
        return backProp




    def forwardPropCheck(self, amplitude=None, sections=None, preStim_delay=None, postStim_run=None, vThreshold=None):
        if self.Always_reset_channels_and_ions:
            self.reset_channels_and_ions()
        if preStim_delay:
            self.preStim_delay = preStim_delay
        if postStim_run:
            self.postStim_run = postStim_run
            
        self.delay=self.preStim_delay
        if amplitude is None:
            self.amplitude=self.maximumAmplitude
        else:
            self.amplitude=amplitude
        if vThreshold is not None:
            self.vThreshold = vThreshold

        if sections is None:
            sections=[self.ranvier_list[-3]]
        segments = []
        for sec in sections:
            segments+=[seg for seg in sec.allseg()]
        
        # print_args(segments)
        self.reset_vpeak()
        self.runStim()
        number = self.mean_vpeak(sections=sections)
        if number>=self.vThreshold:
            forwardProp=True
            self.vpeak_SDtips = self.mean_vpeak_SDtips()
        else:
            forwardProp=False
        print_args(number)
        print_args(forwardProp)
        print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')
        return forwardProp

    def propCheck(self, amplitude=None, backProp_sections=None, forwardProp_sections=None, preStim_delay=None, postStim_run=None, vThreshold=None):
        if self.Always_reset_channels_and_ions:
            self.reset_channels_and_ions()
        if preStim_delay:
            self.preStim_delay = preStim_delay
        if postStim_run:
            self.postStim_run = postStim_run
            
        self.delay=self.preStim_delay
        if amplitude is None:
            self.amplitude=self.maximumAmplitude
        else:
            self.amplitude=amplitude
        if vThreshold is not None:
            self.vThreshold = vThreshold

        if backProp_sections is None:
            segments_BP = self.get_tip_segments()
            backProp_sections = []
            for seg in segments_BP:
                sec = seg.sec
                if sec not in backProp_sections:
                    backProp_sections+=[sec]
        else:
            segments_BP = []
            for sec in backProp_sections:
                segments_BP+=[seg for seg in sec.allseg()]
        

        if forwardProp_sections is None:
            forwardProp_sections=[self.ranvier_list[-3]]
        segments_FP = []
        for sec in forwardProp_sections:
            segments_FP+=[seg for seg in sec.allseg()]



        self.reset_vpeak()
        self.runStim()

        backProp=False
        forwardProp=False
        number_BP = self.mean_vpeak(sections=backProp_sections)
        number_FP = self.mean_vpeak(sections=forwardProp_sections)
        if number_BP>=self.vThreshold:
            backProp=True
            self.vpeak_SDtips = self.mean_vpeak_SDtips()
        if number_FP>=self.vThreshold:
            forwardProp=True
            self.vpeak_SDtips = self.mean_vpeak_SDtips()

        print_args(number_BP)
        print_args(backProp)
        print_args(number_FP)
        print_args(forwardProp)
        print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')
        return [backProp, forwardProp]


    def print_thresh_settings(self, show_it_all=False):
        print_args(self.gamma)
        print_args(self.deviation)
        print_args(self.CrossOverPosition)
        print_args(self.stimLoc)
        print_args(self.preStim_delay)
        print_args(self.postStim_run)
        if show_it_all == True:
            import AllParams as ap
            ap.show_it_all(cell)


    def backProp_thresh(self, startingAmplitude=None, delta_min=None):
        self.Always_reset_channels_and_ions = False
        self.reset_channels_and_ions()
        self.print_thresh_settings()
        if startingAmplitude is None:
            startingAmplitude = self.maximumAmplitude
        self.amplitude = startingAmplitude
        if delta_min is None:
            delta_min = self.delta_min

        backProp = self.backPropCheck(amplitude=self.amplitude)
        prev_amplitude = self.ic.amp
        if not backProp:
            print('No backProp at startingAmplitude = ', startingAmplitude)
            BP_amplitude = np.inf
        else:
            while backProp:
                BP_amplitude = self.ic.amp
                self.amplitude=prev_amplitude*0.5
                backProp = self.backPropCheck(amplitude=self.amplitude)
                prev_amplitude = self.ic.amp
            delta = 0.5*(BP_amplitude - prev_amplitude)
            keepHoning=True
            while keepHoning is True:
                while not backProp:
                    prev_amplitude = self.ic.amp
                    backProp = self.backPropCheck(amplitude=prev_amplitude+delta)
                BP_amplitude = self.ic.amp
                delta = 0.5*(BP_amplitude - prev_amplitude)
                if delta < delta_min:
                    keepHoning = False
                    break

                while backProp:
                    BP_amplitude = self.ic.amp
                    backProp = self.backPropCheck(amplitude=BP_amplitude-delta)
                    prev_amplitude = self.ic.amp
                delta = 0.5*(BP_amplitude - prev_amplitude)
                if delta < delta_min:
                    keepHoning = False
                    break
            
        self.Always_reset_channels_and_ions = True
        self.iThresh_backprop = BP_amplitude
        self.Qthresh_backProp = self.iThresh_backprop*self.ic.dur
        self.print_thresh_settings()
        print_args(self.Qthresh_backProp)
        return self.Qthresh_backProp


    def forwardProp_thresh(self, startingAmplitude=None, delta_min=None):
        self.Always_reset_channels_and_ions = False
        self.propFunc = self.forwardPropCheck 
        self.reset_channels_and_ions()
        self.print_thresh_settings()
        if startingAmplitude is None:
            startingAmplitude = self.maximumAmplitude
        self.amplitude = startingAmplitude
        if delta_min is None:
            delta_min = self.delta_min

        propagation = self.propFunc(amplitude=self.amplitude)
        prev_amplitude = self.ic.amp
        if not propagation:
            print('No propagation at startingAmplitude = ', startingAmplitude)
            thresh_amplitude = np.inf
        else:
            while propagation:
                thresh_amplitude = self.ic.amp
                self.amplitude=prev_amplitude*0.5
                propagation = self.propFunc(amplitude=self.amplitude)
                prev_amplitude = self.ic.amp
            delta = 0.5*(thresh_amplitude - prev_amplitude)
            keepHoning=True
            while keepHoning is True:
                while not propagation:
                    prev_amplitude = self.ic.amp
                    propagation = self.propFunc(amplitude=prev_amplitude+delta)
                thresh_amplitude = self.ic.amp
                delta = 0.5*(thresh_amplitude - prev_amplitude)
                if delta < delta_min:
                    keepHoning = False
                    break

                while propagation:
                    thresh_amplitude = self.ic.amp
                    propagation = self.propFunc(amplitude=thresh_amplitude-delta)
                    prev_amplitude = self.ic.amp
                delta = 0.5*(thresh_amplitude - prev_amplitude)
                if delta < delta_min:
                    keepHoning = False
                    break
            
        self.Always_reset_channels_and_ions = True
        self.iThresh_forwardProp = thresh_amplitude
        self.Qthresh_forwardProp = self.iThresh_forwardProp*self.ic.dur
        self.print_thresh_settings()
        print_args(self.Qthresh_forwardProp)
        return self.Qthresh_forwardProp

    def axialCurrent_data(self, sections=None):
        if sections is None:
            sections=[self.ais]


        segments = []
        vectors = []
        axial_resistance = []
        for sec in sections:
            for seg in sec:
                segments+=seg
                vectors+=h.Vector().record(seg._ref_v)
                axial_resistance.append(seg.ri())
        self.axialCurrent_data = dict(segments=segments, vectors=vectors, axial_resistance=axial_resistance)
        return self.axialCurrent_data




    def sweeper(self, gamma_vals=None, deviation_vals=None, CrossOverPosition_vals=None, vRS03_vals=None):
        if gamma_vals is None:
            gamma_vals = [20.0]

        if deviation_vals is None:
            deviation_vals = [1.0, 0.0] #[1.0, 0.9, 0.8, 0.7, 0.5, 0.0] #[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] #[1.0, 0.5, 0.1] #[1.0, 0.0, 0.5, 0.2, 0.1]

        if CrossOverPosition_vals is None:
            CrossOverPosition_vals = [self.CrossOverPosition]

        if vRS03_vals is None:
            vRS03_vals = [self.vRS03]


        thresHolding_function = self.backProp_thresh #self.backProp_thresh  #self.forwardProp_thresh
        for vRS03 in vRS03_vals:
            for gamma in gamma_vals:
                for CrossOverPosition in CrossOverPosition_vals:
                    for deviation in deviation_vals:
                        self.reset_channels_and_ions(gamma=self.gamma, gammaAIS=gamma, deviation=deviation, CrossOverPosition=CrossOverPosition)
                        self.selectiveRightShift_NaV12(section=ais, vRS03=vRS03)
                        self.Qthresh = thresHolding_function(startingAmplitude=self.maximumAmplitude)

                        self.save_to_database()
                        self.update_pandas_df()
                        print_args(self.stimLoc_data)
                        # print_args(self.gamma_data)
                        print_args(self.gammaAIS_data)
                        print_args(self.deviation_data)
                        print_args(self.CrossOverPosition_data)
                        print_args(self.Qthresh_data)
                        print_args(self.Qthresh_backProp_data)
                        print_args(self.Qthresh_forwardProp_data)
                        print_args(self.vpeak_data)
                        print_args(self.vRS03_data)
                        print_args(self.DeltaVrs_data)


    def spaceArray(self, sec, varname):
        a_t = []
        return 0 

    def ignitionData(self, dt_animate=None):
        if dt_animate is not None:
            self.dt_animate = dt_animate

        self.ic.amp = self.amplitude
        self.ic.delay = self.h.t+self.delay
        if self.concentrationReinit:
            self.concReinit()
            self.h.cvode.re_init()
        print_args(self.df_filename)
        print_args(self.pandas_df)
        print('____________________________________________________')
        print_args(self.h.t)
        print_args(self.ic.delay)
        print_args(self.ic.amp)
        self.h.continuerun(self.ic.delay)




        return 0




    def update_data_lists(self):
        self.Qthresh_data.append(self.Qthresh)
        self.Qthresh_backProp_data.append(self.Qthresh_backProp)
        self.Qthresh_forwardProp_data.append(self.Qthresh_forwardProp)
        self.vpeak_data.append(self.vpeak_SDtips)
        self.vThreshold_data.append(self.vThreshold)
        self.stimLoc_data.append(self.stimLoc)
        # self.gamma_data.append(self.gamma)
        self.gammaAIS_data.append(self.gammaAIS)
        self.deviation_data.append(self.deviation)
        self.CrossOverPosition_data.append(self.CrossOverPosition)
        
        self.vRS03_data.append(self.vRS03)
        self.DeltaVrs_data.append(self.vRS03 - 13.0)
        

    def update_pandas_df(self, print_the_df=True):
        self.update_data_lists()

        self.pandas_data['stimLoc'] = self.stimLoc_data
        self.pandas_data['Qthresh'] = self.Qthresh_data
        self.pandas_data['Qthresh_backProp'] = self.Qthresh_backProp_data
        self.pandas_data['Qthresh_forwardProp'] = self.Qthresh_forwardProp_data
        self.pandas_data['vpeak'] = self.vpeak_data
        self.pandas_data['vThreshold'] = self.vThreshold_data
        # self.pandas_data['gamma'] = self.gamma_data
        self.pandas_data['gammaAIS'] = self.gammaAIS_data
        self.pandas_data['deviation'] = self.deviation_data
        self.pandas_data['CrossOverPosition'] = self.CrossOverPosition_data
        self.pandas_data['vRS03'] = self.vRS03_data
        self.pandas_data['DeltaVrs'] = self.DeltaVrs_data
        # self.pandas_data['fraction'] = self.fraction_data
        self.pandas_df = pd.DataFrame.from_dict(self.pandas_data)
        self.pandas_df.to_csv(self.df_filename)
        print_args(self.df_filename)
        if print_the_df is True:
            print_args(self.pandas_df)

    def save_to_database(self):
        import pandas as pd
        import sqlite3
        import time

        import run_settings

        start = time.perf_counter()
        data = pd.DataFrame(
            {
            'simulation_time': [self.h.t],
            'stimLoc' : str(self.stimLoc),
            'vThreshold' : self.vThreshold,
            'temperature': self.h.celsius,
            'ais.L' : self.ais.L,
            'vRS0' : self.vRS0,
            'disableRightShift_minf2__inAIS' : run_settings.disableRightShift_minf2__inAIS,
            'disableRightShift_mtau2__inAIS' : run_settings.disableRightShift_mtau2__inAIS,
            'disableRightShift_hinf2__inAIS' : run_settings.disableRightShift_hinf2__inAIS,
            'disableRightShift_htau2__inAIS' : run_settings.disableRightShift_htau2__inAIS,
            'vRS03' : self.vRS03,
            'DeltaVrs' : self.vRS03 - 13.0,
            'deviation' : self.deviation,
            'CrossOverPosition' : self.CrossOverPosition,
            'Qthresh' : self.Qthresh,
            'Qthresh_backProp' : self.Qthresh_backProp,
            'Qthresh_forwardProp' : self.Qthresh_forwardProp,
            'vpeak_SDtips' : self.vpeak_SDtips,
            }
        )
        with sqlite3.connect('sqlite3'+self.df_filename[:-4]+'.sqlite') as conn: 
            data.to_sql("data", conn, if_exists="append", index=False)




if __name__=='__main__':
    cell = StartUp(gamma=50.0, deviation=1.0, CrossOverPosition=0.5)
    cell.sweeper(gamma_vals=[20.0], deviation_vals=x_list, CrossOverPosition_vals=kappa_list)
    




































