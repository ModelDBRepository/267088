def makeSomatodendritic(h):
    from sorcery import dict_of, print_args  
    if __name__=='__main__': 
        h.load_file('P_Soma_Dendrites.hoc')                    
    else:                                                           
        h.load_file('morphology/'+'P_Soma_Dendrites.hoc')
    soma = h.soma
    dend = h.dend11[0]  #<--- dend11: primary apical dendrite


    def tuneup_somatodend():
        ###build a sectionlist for soma and dendrites
        somatodendritic = h.SectionList()
        for sec in h.allsec():
            ###make sure no segments exceed 10 uM length. Note, soma.nseg remains 10.
            if sec.L/sec.nseg>10:
                sec.nseg = int(sec.L/10.0) + 1
                ###Keep the number of segments odd in every section
            if sec.nseg % 2 == 0:
                sec.nseg+=1
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

    for sec in somatodendritic:
        sec.Ra = 150
        sec.cm = 1.0
    


    ###__________________________________________________________________________
    from mech_settings import mech_name_soma as mech_name
    from mech_settings import rescale_soma, rescale_dend
    from mech_settings import insertCLS
    insertCLS(soma, mech_name, rescale_soma)
    for sec in dendritic:
        insertCLS(sec, mech_name, rescale_dend)
    return dend, soma, sd_list

