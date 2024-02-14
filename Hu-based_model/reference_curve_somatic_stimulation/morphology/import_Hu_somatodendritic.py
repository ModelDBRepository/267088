import os, sys
import numpy as np
from sorcery import dict_of, print_args
import matplotlib.pyplot as plt
from neuron import h, gui, crxd as rxd                                                                                    
from neuron.units import ms, mV
h.load_file('P_Soma_Dendrites.hoc')
soma = h.soma


def tuneup_somatodend():
    ###build a sectionlist for soma and dendrites
    somatodendritic = h.SectionList()
    for sec in h.allsec():
        if sec.L/sec.nseg>40:
            sec.nseg = int(sec.L/40.0) + 1
            ###Keep the number of segments odd in every section
            if sec.nseg % 2 == 0:
                sec.nseg+=1
        ###make sure no segments exceed 40 uM length. Note, soma.nseg remains 10.
        somatodendritic.append(sec)

    ###build a sectionlist for dendrites only
    dendritic = h.SectionList()
    for sec in somatodendritic:
        dendritic.append(sec)
    dendritic.remove(soma) #<---remove soma for pure dendritic sectionlist
    h.distance(0, soma(0))  #<--- set the point where axon seated on soma as the origin
    return somatodendritic, dendritic
somatodendritic, dendritic = tuneup_somatodend()

def add_spines():
    ###Translated from P_Spines.hoc
    spine_dens = 1    #<--- just using a simple spine density model due to lack of data on some neuron types.
    spine_area = 0.83 #<--- um^2  -- K Harris
    for sec in dendritic:
        a = 0
        for seg in sec:
            a+= seg.area()
        F = (sec.L*spine_area*spine_dens + a)/a

        sec.L = sec.L*(F**(2.0/3.0))
        for seg in sec:
            seg.diam = seg.diam*(F**(1.0/3.0))
add_spines()

sd_list = list(somatodendritic)
d_list = list(dendritic)
  

hill = h.Section(name='hill')
hill.L = 10.0
hill.diam = 10.0
hill.connect(soma(1), 0)


# h.dend11[20].diam = 40.0


h.define_shape()
from blenderneuron.quick import bn # import the BlenderNEURON library
assert bn.is_blender_ready()       # check if communication with Blender can be established
bn.to_blender()                    # send section morphology data to Blender

