def makeNaked(h):
    # from neuron import h
    from run_settings import diam_ais
    sec = h.Section(name='naked')
    sec.L = 400.0 #0.1 #400.0
    sec.diam = diam_ais #1.02 #1.6  
    sec.Ra = 150
    sec.cm = 1.0 
    sec.nseg = 45 #3 #45
    

    ###__________________________________________________________________________
    from mech_settings import mech_name_naked as mech_name
    from mech_settings import rescale_naked as rescale
    from mech_settings import insertCLS
    insertCLS(sec, mech_name, rescale)

    return sec