def makeRanvier(h, j):
    # from neuron import h
    import numpy as np
    from run_settings import diam_ais
    sec = h.Section(name='ranvier'+str(j))

    Area = 6.0; VolumeIn = 3.0 #; VolumeOut = 3.0
    length = (Area**2.0)/(4.0*np.pi*VolumeIn)
    radius = 2.0*VolumeIn/Area ; diameter = 2.0*radius

    sec.L = 1.0 #length #1.0  #<---Kole        #length
    sec.diam = diam_ais #1.1 #diameter #1.1 #<---Kole        #diameter
    sec.Ra = 150.0
    sec.cm = 1.0
    sec.nseg = 11


    ###__________________________________________________________________________
    from mech_settings import mech_name_ranvier as mech_name
    from mech_settings import rescale_ranvier as rescale
    from mech_settings import insertCLS
    insertCLS(sec, mech_name, rescale)

    return sec

