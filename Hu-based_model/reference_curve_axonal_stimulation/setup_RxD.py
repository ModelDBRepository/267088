def RxD_setConcentrations(ions, cvode):
    import numpy as np
    from run_settings import chloride
    from run_settings import conc_nai, conc_nao, conc_ki, conc_ko
    if chloride==True:
        from run_settings import conc_cli, conc_clo
        concentrations=[conc_nai, conc_nao, conc_ki, conc_ko, conc_cli, conc_clo]
    else:
        concentrations=[conc_nai, conc_nao, conc_ki, conc_ko]
    if len(ions) != len(concentrations):
        print("ERROR: ions list does not match concentrations list")
        exit()
    else:
        N_concentrations = len(concentrations)
        for j in np.arange(N_concentrations):
            ions[j].nodes.concentration = concentrations[j]
        # print("Concentrations have been reinitialized")

    from run_settings import use_CVODE
    if use_CVODE == True:
        cvode.re_init()




def setup_RxD_NaKCl(h, rxd, cvode, sections, diffusion_on, d_Na, d_K, d_Cl):
    import numpy as np
    from run_settings import Vrest
    if diffusion_on==False:
        d_Na = 0.0
        d_K = 0.0
        d_Cl = 0.0
    from measurement import cell_length
    ecs_volume_multiplier = 100.0 #0.1 #100.0
    ###Where?
    cyt = rxd.Region(sections, name='cyt', nrn_region='i')
    # ecs = rxd.Region(sections, name='ecs', nrn_region='o')
    ecs = rxd.Region(sections, name='ecs', nrn_region='o', geometry=rxd.Shell(lo=1.0, hi=np.sqrt(ecs_volume_multiplier+1.0)))

    # ecs = rxd.Extracellular(-100, -100, -100, cell_length(h)+100, 100, 100, dx=20)

    ###Who?
    na = rxd.Species([cyt, ecs], name='na', d=d_Na, charge=1, atolscale=1.0)
    k = rxd.Species([cyt, ecs], name='k', d=d_K, charge=1, atolscale=1.0)
    cl = rxd.Species([cyt, ecs], name='cl', d=d_Cl, charge=-1, atolscale=1.0)


    # naecs = rxd.Species(ecs, name='naecs', d=0.6, charge=1.0, initial=154.0)
    # kecs = rxd.Species(ecs, name='kecs', d=1, charge=1.0, initial=6.0)

    ki, ko, nai, nao, cli, clo = k[cyt], k[ecs], na[cyt], na[ecs], cl[cyt], cl[ecs]

    h.finitialize(Vrest)
    RxD_setConcentrations(ions=[nai, nao, ki, ko, cli, clo], cvode=cvode)

    # nai.nodes.concentration = 10.0 #20.0
    # nao.nodes.concentration = 145.0 #154.0

    # ki.nodes.concentration = 140.0 #150.0
    # ko.nodes.concentration = 5.0 #6.0

    # cli.nodes.concentration = 4.0 #4.0
    # clo.nodes.concentration = 110.0 #110.0

    # from run_settings import use_CVODE
    # if use_CVODE == True:
    #     cvode.re_init()

    return ki, ko, nai, nao, cli, clo, cyt, ecs, na, k, cl


    


def setup_RxD_active_outside2(h, rxd, cvode, sections, diffusion_on, d_Na, d_K):
    if diffusion_on==False:
        d_Na = 0.0
        d_K = 0.0
    from measurement import cell_length
    from run_settings import Vrest

    ###Where?
    cyt = rxd.Region(sections, name='cyt', nrn_region='i')
    ecs = rxd.Region(sections, name='ecs', nrn_region='o')
    # ecs = rxd.Region(sections, name='ecs', nrn_region='o', geometry=rxd.Shell(lo=1.0, hi=2.0))
    # ecs = rxd.Extracellular(-100, -100, -100, cell_length(h)+100, 100, 100, dx=20)

    ###Who?
    na = rxd.Species([cyt, ecs], name='na', d=d_Na, charge=1, atolscale=1.0)
    k = rxd.Species([cyt, ecs], name='k', d=d_K, charge=1, atolscale=1.0)


    # naecs = rxd.Species(ecs, name='naecs', d=0.6, charge=1.0, initial=154.0)
    # kecs = rxd.Species(ecs, name='kecs', d=1, charge=1.0, initial=6.0)

    ki, ko, nai, nao = k[cyt], k[ecs], na[cyt], na[ecs]

    h.finitialize(Vrest)
    RxD_setConcentrations(ions=[nai, nao, ki, ko], cvode=cvode)
    # for node in nai.nodes:
    #     node.concentration = 20.0
    # for node in nao.nodes:
    #     node.concentration = 154.0

    # for node in ki.nodes:
    #     node.concentration = 150.0
    # for node in ko.nodes:
    #     node.concentration = 6.0

    # from run_settings import use_CVODE
    # if use_CVODE == True:
    #     cvode.re_init()


    return ki, ko, nai, nao, cyt, ecs, na, k







def setup_RxD_old(h, rxd, diffusion_on, d_Na, d_K):
    if diffusion_on==False:
        d_Na = 0.0
        d_K = 0.0

    ###Where?
    length = 0
    for sec in h.allsec():
        length+=sec.L
    print(length)

    cyt = rxd.Region(h.allsec(), name='cyt', nrn_region='i')
    ecs = rxd.Extracellular(-100, -100, -100, (length+100), 100, 100, dx=33)

    ###Who?
    na = rxd.Species(cyt, name='na', d=d_Na, charge=1.0, initial=20.0, atolscale=1e-6)
    k = rxd.Species(cyt, name='k', d=d_K, charge=1.0, initial=150.0, atolscale=1e-6)

    # extracellular parameters provide a constant concentration for the Nernst potential and reactions.
    naecs = rxd.Parameter(ecs, name='na', charge=1.0, value=154.0)
    kecs = rxd.Parameter(ecs, name='k', charge=1.0, value=6.0)


    ki, ko, nai, nao = k[cyt], kecs[ecs], na[cyt], naecs[ecs]
    return ki, ko, nai, nao, cyt, ecs, na, k, naecs, kecs
