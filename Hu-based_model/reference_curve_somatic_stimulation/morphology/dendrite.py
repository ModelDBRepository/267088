def makeDendrite(h, name='dend', L=100.0):
    # from neuron import h
    sec = h.Section(name=name)

    sec.L = L # microns
    sec.diam = 1 # microns
    sec.Ra = 150
    sec.cm = 1.0
    sec.nseg = 11



    ###__________________________________________________________________________
    from mech_settings import mech_name_dend as mech_name
    from mech_settings import rescale_dend as rescale
    from mech_settings import insertCLS
    insertCLS(sec, mech_name, rescale)

    return sec

