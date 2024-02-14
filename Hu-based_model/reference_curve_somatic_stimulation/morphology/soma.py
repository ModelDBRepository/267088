def makeSoma(h):
    # from neuron import h
    sec = h.Section(name='soma')

    sec.L = 20 
    sec.diam = 20
    sec.Ra = 150
    sec.cm = 1.0 
    sec.nseg = 11


    ###__________________________________________________________________________
    from mech_settings import mech_name_soma as mech_name
    from mech_settings import rescale_soma as rescale
    from mech_settings import insertCLS
    insertCLS(sec, mech_name, rescale)

    return sec

