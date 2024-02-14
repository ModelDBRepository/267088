"""This script is included because it is used by too many of the scripts involved in the Hu-based model. 
It contains a lot of unnecessary pieces that have been left in to avoid breaking the model unexpectedly. 
Aside from the settings outlined below, everything else should be left alone. 
Many of the other settings are overwritten at other points in the simulator, but some are necessary during initialization. 

This code was written before I had much familiarity with object-oriented programming. 
The code for the Hay-based model is more recent and much more developer friendly.

One can still use this program to control the following settings:

For generating "backward AIS" results:
-FlipNavs_LeftRight

The following flags should be set to "False", unless running simulations related to Fig 6:
-disableRightShift_minf2__inAIS
-disableRightShift_mtau2__inAIS
-disableRightShift_hinf2__inAIS
-disableRightShift_htau2__inAIS

The following variable can be set to "ais" or "soma". 
It sets the location of the stimulating electrode:
-stimLoc

"""

from neuron.units import ms, mV, mM
import numpy as np
from settings_checker import LS_AC__value_checker
from sorcery import print_args

#______________________________________________________________________________
use_CVODE = True

cell_connected__FLAG = True
taper_hillock__FLAG = True
taper_dendrite__FLAG = True
FlipNavs_LeftRight = False #<---'backward' AIS Nav distributions (preserves crossover location)


disableRightShift_minf2__inAIS = True
disableRightShift_mtau2__inAIS = True

disableRightShift_hinf2__inAIS = True
disableRightShift_htau2__inAIS = True




chloride=True
Hu_somaDend = True
N_ranvier = 15
temperature = 37.0
vThreshold = 0.0*mV
gamma = 30.0 #50.0 #<---Kole, Hu, gamma~30,40,50
gl = 0.0001 #0.000036
gl_rescaling=True


#_____________________________________________________________________________

###Left-Shift (Nav1.2 and Nav1.6)
vLeftShift = 0.0*mV #20.0 #15.0 #10.0 #6.0 #8.0 #1.0 #5.5
AC = 1.0

vRS0 = 13.0*mV #<---Right-Shift (Nav1.2)
vRS03 = 0.0*mV
deviation = 1.0
CrossOverPosition = 0.5
fraction = 1.0

L_ais = 25.0 #50.0 #100.0 
L_hill_extend = int(0) #int(10.0) #<---Increases ais distance
diam_ais = 1.22 #1.22 #1.6


stimLoc = 'soma' #'ais' #'soma'
stimDur = 1.0*ms  


amplitude = 1.0 
charge=amplitude*stimDur


Ra = 150.0 
Ra_axon = Ra 

cm = 0.75
cm_myelin = 0.02

Vrest = -70.0*mV
vShiftNa = 0.0 
vShift_h = 0.0 
vShiftK = 0.0 
Ecl = Vrest
Ena = 60.0*mV 
Ek = -90.0*mV 

gk_leak = 0.0001



def conc_io(zx, Ex, xi, temperature=temperature):
    from neuron.units import mV
    elementary_charge = 1.60217662e-19
    qx = zx*elementary_charge
    Tkelvin = temperature+273.15
    Kboltzmann = 1.38064852e-23
    KT = Kboltzmann*Tkelvin
    KT_over_q = (KT/qx)*1000.0*mV #<---multiply by 1000 for conversion to millivolts
    # print_args(KT_over_q)
    xo = xi*np.exp(Ex/KT_over_q)
    return xo

conc_nai = 15.0*mM  #(5-15 mM)
conc_nao = conc_io(+1.0, Ena, conc_nai)
conc_ki = 140.0*mM #140.0
conc_ko = conc_io(+1.0, Ek, conc_ki)

if chloride==True:
    conc_cli =  11.0*mM #11.65*mM #11.5*mM #9.9*mM
    conc_clo = conc_io(-1.0, Ecl, conc_cli)




Kmko = 3.5
Kmnai = 10.0
ImaxPump =  0.7*9.09e-2 
gk_leak = ImaxPump/( ((Vrest-Ek)/2.0)*((1.0 + Kmko/conc_ko)**2.0)*((1.0 + Kmnai/conc_nai)**3.0) )


gna_leak = (3.0/2.0)*gk_leak*((Vrest - Ek)/(Ena - Vrest))


def classic_settings():
    global gamma
    global chloride
    global Vrest
    global vShiftNa
    global Ena
    global Ek
    global ImaxPump
    global gna_leak
    global gk_leak
    global conc_nai
    global conc_nao
    global conc_ki
    global conc_ko
    global conc_cli
    global conc_clo

    gamma = 100.0

    Vrest = -59.9
    vShiftNa = 0.0*mV
    vShiftK = 0.0*mV
    vShift_h = 0.0*mV
    Ena = 51.5
    Ek = -81.3
    ImaxPump = 9.09e-2
    gk_leak = 0.0001
    gna_leak = 0.00025

    conc_nai = 20.0*mM 
    conc_nao = 154.0*mM  
    conc_ki = 150.0*mM 
    conc_ko = 6.0*mM  
    if chloride==True:
        conc_cli =  11.65*mM 
        conc_clo = 125.0*mM 
# classic_settings()






LS_AC__value_checker(vLeftShift, AC)


equilm_dur = 7.0*ms #1500000.0*ms
resume_dur = 10.0*ms #5.0*ms
equilibration_timestep = 0.01*ms
stimulation_timestep = 0.001*ms 
CVODE_pre_equilm_dur = 3.0*ms




if vLeftShift!=0.0*mV:
    timeLS = 50.0*ms #200.0*ms
else:
    timeLS = 0.0*ms
equilm_dur+=CVODE_pre_equilm_dur+timeLS
stimON = equilm_dur+1.0*ms
total_dur = equilm_dur+int(stimDur)+resume_dur
fixed_step_dur = stimON+11.0*ms
# fixed_step_dur = stimON+int(stimDur)+10.0*ms
stimOFF=stimON+stimDur

#‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾


#__________________________________________________________________________________________________________________________________________________________
#■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
#■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
#■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
###Grid building (copy_then_import_then_run.py)
N_charges = 3
N_equilibrations = 50 #<---Number of values for the interesting (and expensive) parameter. Eg. deviation, AISNav12inhibition(fractions), gamma. It's what we're SWEEPING.

charges = np.linspace(0.001, 5.0, N_charges)
def create_runs_dict():
    indexes = np.arange(N_equilibrations)
    runs_dict = {} #<---Parameters that are being swept (system is equilibrated each time these are changed)
    runs_dict['indexes'] = indexes 
    static = np.ones(len(indexes))

    gammas = np.linspace(0.001, 60.0, len(indexes)) #np.linspace(0.001, 60.0, len(indexes)) #gamma*static
    deviations = 1.0*static #np.linspace(0.001, 1.0, len(indexes)) #1.0*static
    fractions = 1.0*static #np.linspace(1.0, 0.001, len(indexes)) #1.0*static
    vRS0s = vRS0*static #vRS0*static #np.linspace(0.0, 30.0, len(indexes))*mV
    vLeftShifts = vLeftShift*static #np.linspace(0.0, 20.0, len(indexes))*mV #vLeftShift*static
    ACs = AC*static #1.0*static

    for index in indexes:
        temporary_dict = {}
        temporary_dict['gamma'] = gammas[index]
        temporary_dict['deviation'] = deviations[index] 
        temporary_dict['fraction'] = fractions[index] 
        temporary_dict['vRS0'] = vRS0s[index]
        temporary_dict['vLeftShift'] = vLeftShifts[index]
        temporary_dict['AC'] = ACs[index]
        # temporary_dict[''] = 
        runs_dict[str(index)] = temporary_dict
    return runs_dict, indexes

Naked_axon_postAIS = True #<---Deprecated: Leave as "True" since the simulator expects all Sections to be present. Sections can be shortened. 

def create3D_runs_dict():
    global charges
    charges = np.array([0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 20.0, 50.0])

    # charges = np.array([2.0, 4.0, 5.0, 7.0, 10.0, 20.0]) #<---stimLoc = h.dend11[1](0.9)
    # charges = np.array([0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0]) #<---stimLoc = naked(0.5)

    N_gammas = 20
    N_vRS0s = 20
    N_deviations = 5


    indexes = np.arange(N_gammas*N_vRS0s*N_deviations)
    gamma_array = np.linspace(0.001, 200, N_gammas)
    deviations_array = np.linspace(0.0, 1.0, N_deviations)
    vRS0s_array = np.linspace(0.0, 20.0, N_vRS0s)



    gammas = []
    deviations = []
    vRS0s = []
    for vRS0 in vRS0s_array:
        for gamma in gamma_array:
            for deviation in deviations_array:
                vRS0s.append(vRS0)
                gammas.append(gamma)
                deviations.append(deviation)



                
    runs_dict = {} #<---Parameters that are being swept (system is equilibrated each time these are changed)
    runs_dict['indexes'] = indexes 
    static = np.ones(len(indexes))

    vLeftShifts = vLeftShift*static #np.linspace(0.0, 20.0, len(indexes))*mV #vLeftShift*static
    ACs = AC*static #1.0*static
    fractions = 1.0*static #np.linspace(1.0, 0.001, len(indexes)) #1.0*static

    for index in indexes:
        temporary_dict = {}
        temporary_dict['gamma'] = gammas[index]
        temporary_dict['deviation'] = deviations[index] 
        temporary_dict['fraction'] = fractions[index] 
        temporary_dict['vRS0'] = vRS0s[index]
        temporary_dict['vLeftShift'] = vLeftShifts[index]
        temporary_dict['AC'] = ACs[index]
        # temporary_dict[''] = 
        runs_dict[str(index)] = temporary_dict

    return runs_dict, indexes

# runs_dict, indexes = create_runs_dict()
runs_dict, indexes = create3D_runs_dict()



def timed_run(h, duration):
    from time import time
    start = time()

    h.continuerun(duration)
    # h.run(duration)

    end = time()
    RunTime = end - start
    print('RunTime = ', RunTime, ' seconds')
#■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
#■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
#■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
#‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾


if __name__ == '__main__':
    import os, sys
    from mech_settings import mech_name
    mech_filename = mech_name+'.mod'
    os.system('nrnivmodl '+mech_filename)
    print("--------NMODL mechanism compiled----------------------------------------"+'\n\n')
    
    # from plot_mhn import plot_mhn ; plot_mhn()

    from overhead_GridBuilder import single_run, single_blender, blender_run, show_equilibration, find_thresh
    single_run()

    


