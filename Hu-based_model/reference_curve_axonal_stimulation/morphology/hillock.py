def makeHillock(h):
    # from neuron import h
    from run_settings import diam_ais
    sec = h.Section(name='hill')
    
    sec.L = 10.0
    sec.diam = diam_ais
    sec.Ra = 150.0
    sec.cm = 1.0 
    sec.nseg = 51 


    ###__________________________________________________________________________
    from mech_settings import mech_name_hill as mech_name
    from mech_settings import rescale_hill as rescale
    from mech_settings import insertCLS
    insertCLS(sec, mech_name, rescale)

    return sec




