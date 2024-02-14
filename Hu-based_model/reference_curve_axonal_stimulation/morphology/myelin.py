def makeMyelin(h, j):
    # from neuron import h
    from run_settings import diam_ais
    sec = h.Section(name='myelin'+str(j))

    sec.L = 100.0 #60.0#<---Kole
    sec.diam = diam_ais #1.6 #<---Kole
    sec.Ra = 150.0
    sec.cm = 0.02 #<---Kole, Hu
    sec.nseg = 21


    ###__________________________________________________________________________
    from mech_settings import mech_name_myelin as mech_name
    from mech_settings import rescale_myelin as rescale
    from mech_settings import insertCLS
    insertCLS(sec, mech_name, rescale)

    return sec

