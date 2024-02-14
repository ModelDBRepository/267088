
###Reference values (Cooling Reverses Pathological[Chaos paper], PAB,...)
gk_node = 0.036   
gk_0 = gk_node


from run_settings import gl, gna_leak, gk_leak, ImaxPump, gamma, Hu_somaDend, chloride

if chloride==True:
    mech_name='CLSTmainen'
else:
    mech_name = 'III_RxD'

mech_name_dend = mech_name
mech_name_soma = mech_name
mech_name_hill = mech_name
mech_name_ais = mech_name
mech_name_myelin = mech_name
mech_name_ranvier = mech_name
mech_name_naked = mech_name


rescale_ranvier = 6.69 
rescale_pump_relative = 1.0
rescale_myelin = 0.02*rescale_ranvier 


if Hu_somaDend==True:
    rescale_soma = 1.25*0.025*rescale_ranvier 
else:
    rescale_soma = 0.025*rescale_ranvier #0.02*rescale_ranvier#<---Kole #0.025*rescale_ranvier #<---Hu
###Decrease dendritic excitability:
rescale_dend = 0.1*rescale_soma # 0.025*rescale_ranvier
###Increase somatic excitability: 
rescale_soma *= 3.0 #<--- If you do this, make sure to divide gamma and gammaAIS by three in Thresher.py
rescale_hill = rescale_soma


rescale_ais = gamma*rescale_soma



rescale_naked = (27.0/34.0)*rescale_ais #<---Hu


def plot_CLS_vRS(h):
    from BBplotting import BarlowRangeVarPlotter
    from mech_settings import mech_name
    BarlowRangeVarPlotter(['vRS'+'_'+mech_name, 'vLS'+'_'+mech_name, 'vLS2'+'_'+mech_name], sections=list(h.allsec()), ylims=[None, None], save_PNG=False)



def initialize_stimON(h, mech_name=mech_name):
    from run_settings import stimON, vThreshold
    for sec in h.allsec():
        for seg in sec.allseg():
            setattr(seg, 'spike'+'_'+mech_name, 0.0)
            setattr(seg, 'vThreshold'+'_'+mech_name, vThreshold)
            setattr(seg, 'stimON'+'_'+mech_name, stimON)
            setattr(seg, 'tpeak'+'_'+mech_name, -1.0)
            setattr(seg, 'tup'+'_'+mech_name, -1.0)
            setattr(seg, 'tdown'+'_'+mech_name, -1.0)
            setattr(seg, 'tspike'+'_'+mech_name, -1.0)





def insertCLS(sec, mech_name, rescale):
    from neuron.units import mV, ms
    from sorcery import print_args
    sec.insert(mech_name)

    from run_settings import gl, gl_rescaling, vShiftNa, vShiftK, vShift_h
    if gl_rescaling==True:
        gl*=rescale

    if any(secName in sec.name() for secName in ('soma', 'dend', 'hill')):
        gna16 = 0.0*rescale
        gna12 = 0.12*rescale
        
    else:
        gna16 = 0.12*rescale
        gna12 = 0.0*rescale



    
    gk = gk_0*rescale
    rescale_pump = rescale_pump_relative*rescale
    INaKmax = ImaxPump*rescale_pump
    gnal = gna_leak*rescale_pump
    gkl = gk_leak*rescale_pump


    for seg in sec.allseg():
        setattr(seg, 'gl'+'_'+mech_name, gl)
        setattr(seg, 'gnabar'+'_'+mech_name, gna16)
        setattr(seg, 'gnabar2'+'_'+mech_name, gna12)
        setattr(seg, 'gkbar'+'_'+mech_name, gk)
        setattr(seg, 'INaKmax'+'_'+mech_name, INaKmax)
        setattr(seg, 'gnal'+'_'+mech_name, gnal)
        setattr(seg, 'gkl'+'_'+mech_name, gkl)

        setattr(seg, 'vShiftNa'+'_'+mech_name, vShiftNa)
        setattr(seg, 'vShiftK'+'_'+mech_name, vShiftK)
        setattr(seg, 'vShift_h'+'_'+mech_name, vShift_h)



def print_params(seg):
    import numpy as np
    from sorcery import print_args
    gl = getattr(seg, 'gl'+'_'+mech_name)
    gnabar = getattr(seg, 'gnabar'+'_'+mech_name)
    gnabar2 = getattr(seg, 'gnabar2'+'_'+mech_name)
    gkbar = getattr(seg, 'gkbar'+'_'+mech_name)
    INaKmax = getattr(seg, 'INaKmax'+'_'+mech_name)
    gnal = getattr(seg, 'gnal'+'_'+mech_name)
    gkl = getattr(seg, 'gkl'+'_'+mech_name)


    gna = getattr(seg, 'gna'+'_'+mech_name)
    gna2 = getattr(seg, 'gna2'+'_'+mech_name)
    gk = getattr(seg, 'gk'+'_'+mech_name)



    r_mPassive = 1.0/(gnal + gkl + gl)
    l_mPassive = np.sqrt((seg.diam*r_mPassive)/(4.0*seg.sec.Ra))
    tau_mPassive = r_mPassive*seg.cm
    g_mPassive = 1.0/r_mPassive

    r_m = 1.0/(gna+gna2+gnal + gk+gkl + gl)
    l_m = np.sqrt((seg.diam*r_m)/(4.0*seg.sec.Ra))
    tau_m = r_m*seg.cm
    g_m = 1.0/r_m

    print('_______________________________________________________')
    print('Parameters for '+str(seg)+ ' are: ') 


    
    conversion_factor = 10.0**(4.0)   #<---S/cm^2 to S/m^2
    mainen_conversion_factor = (1.0/3.21)*conversion_factor

    print( 'gl =', gl*conversion_factor, 'S/m²')
    print( 'gnabar =', gnabar*mainen_conversion_factor, 'S/m²')
    print( 'gnabar2 =', gnabar2*mainen_conversion_factor, 'S/m²')
    print( 'gkbar =', gkbar*mainen_conversion_factor, 'S/m²')
    print( 'INaKmax =', INaKmax, 'mA/cm²')
    print( 'gnal =', gnal*conversion_factor, 'S/m²')
    print( 'gkl =', gkl*conversion_factor, 'S/m²')

    # print('g_m = ', g_m*conversion_factor, 'S/m^2')
    print('g_mPassive = ', g_mPassive*conversion_factor, 'S/m²')
    # print('r_m = ', r_m)
    print('r_mPassive = ', r_mPassive, '𝛺⋅cm²')
    # print('l_m = ', l_m)
    print('l_mPassive = ', l_mPassive)
    # print('tau_m = ', tau_m)
    print('tau_mPassive = ', tau_mPassive)
    # print('ri() ÷ (segment length) = ', seg.ri()/(seg.sec.L/seg.sec.nseg))
    






def rescaleCLS_ΑΙSandNakedAxon(sections, gamma, mech_name=mech_name_ais):
    from run_settings import Naked_axon_postAIS
    rescale_ais = gamma*rescale_soma
    rescale_naked = (27.0/34.0)*rescale_ais #<---Hu
    sec_names = []


    def rescale_section(sec, rescale):
        from run_settings import gl, gl_rescaling
        if gl_rescaling==True:
            gl*=rescale

        gna16 = 0.12*rescale
        gna12 = 0.0*rescale
        gk = gk_0*rescale

        rescale_pump = rescale_pump_relative*rescale
        INaKmax = ImaxPump*rescale_pump
        gnal = gna_leak*rescale_pump
        gkl = gk_leak*rescale_pump
        for seg in sec.allseg():
            setattr(seg, 'gl'+'_'+mech_name, gl)
            setattr(seg, 'gnabar'+'_'+mech_name, gna16)
            setattr(seg, 'gnabar2'+'_'+mech_name, gna12)
            setattr(seg, 'gkbar'+'_'+mech_name, gk)
            setattr(seg, 'INaKmax'+'_'+mech_name, INaKmax)
            setattr(seg, 'gnal'+'_'+mech_name, gnal)
            setattr(seg, 'gkl'+'_'+mech_name, gkl)
        gnabarTotal = gna12 + gna16
        return gnabarTotal




    if not any('hill' in sec.name() for sec in sections):
        print('ERROR: rescaleCLS_ΑΙSandNakedAxon() is designed only for hillock AIS and Bare axon. First argument should be either: [hill, ais, naked] or [hill, ais]')
        exit()
    for sec in sections:
        if sec.name() == 'hill':
            hill_sec = sec
            sec_names.append(sec.name())

        elif sec.name() == 'ais':
            ais_sec = sec
            sec_names.append(sec.name())
            
            
        elif sec.name() == 'naked' and Naked_axon_postAIS==True:
            naked_sec = sec
            sec_names.append(sec.name())

        else:
            print('ERROR: rescaleCLS_ΑΙSandNakedAxon() is designed only for hillock AIS and Bare axon. First argument should be either: [hill, ais, naked] or [hill, ais]')
            exit()

        
        

    gnabarTotal_AIS = rescale_section(ais_sec, rescale_ais)
    gnabarTotal_naked = rescale_section(naked_sec, rescale_naked)


    def channelMatch_hillock_to_AIS(sec):
        from mech_settings import mech_name_hill as mech_name
        import numpy as np
        rescaleVals = np.linspace(rescale_hill, rescale_ais, len(list(sec.allseg())), endpoint=True)
        rescaleVals_NaV16 = np.linspace(0.0, rescale_ais, len(list(sec.allseg())), endpoint=True)

        i_seg = 0
        proximal_Nav16_fraction = getattr(ais_sec(0), 'gnabar'+'_'+mech_name)/gnabarTotal_AIS
        proximal_Nav12_fraction = getattr(ais_sec(0), 'gnabar2'+'_'+mech_name)/gnabarTotal_AIS
        for seg in sec.allseg():
            rescale = rescaleVals[i_seg]
            rescale_NaV16 = rescaleVals_NaV16[i_seg]
            i_seg+=1

            from run_settings import gl, gl_rescaling
            if gl_rescaling==True:
                gl*=rescale




            gna16 = getattr(ais_sec(0), 'gnabar'+'_'+mech_name)*rescale_NaV16/rescale_ais
            gna12 =  getattr(ais_sec(0), 'gnabar2'+'_'+mech_name)*rescale/rescale_ais
            gk = gk_0*rescale
            rescale_pump = rescale_pump_relative*rescale
            INaKmax = ImaxPump*rescale_pump
            gnal = gna_leak*rescale_pump
            gkl = gk_leak*rescale_pump
            setattr(seg, 'gl'+'_'+mech_name, gl)
            setattr(seg, 'gnabar'+'_'+mech_name, gna16)
            setattr(seg, 'gnabar2'+'_'+mech_name, gna12)
            setattr(seg, 'gkbar'+'_'+mech_name, gk)
            setattr(seg, 'INaKmax'+'_'+mech_name, INaKmax)
            setattr(seg, 'gnal'+'_'+mech_name, gnal)
            setattr(seg, 'gkl'+'_'+mech_name, gkl)

    channelMatch_hillock_to_AIS(sec=hill_sec)


    if Naked_axon_postAIS==True:
        naked_str = '; and gnabarTotal_naked = '+str(gnabarTotal_naked)+', which is (27/34)*g_AIS as in Hu2009.' 
    else:
        naked_str=''
    # print('⁍The CLS parameters of the AIS have been rescaled with gamma =', gamma, ' and gnabarTotal =', gnabarTotal_AIS, naked_str)
    # print('⁍sections rescaled = ', sec_names)
    return gnabarTotal_AIS


def linear_Nav_matcher2(soma_sec, hill_sec, ais_sec, gamma=gamma):
        from mech_settings import mech_name_hill as mech_name
        import numpy as np
        rescale_ais = gamma*rescale_soma

        sec = hill_sec
        rescaleVals = np.linspace(rescale_hill, rescale_ais, len(list(sec.allseg())), endpoint=True)
        print(rescaleVals/rescale_ais)
        print(rescale_hill)
        print(rescale_soma)
        # exit()
        i_seg = 0
        gnabarTotal_AIS = getattr(ais_sec(0), 'gnabar'+'_'+mech_name) + getattr(ais_sec(0), 'gnabar2'+'_'+mech_name)
        proximal_Nav16_fraction = getattr(ais_sec(0), 'gnabar'+'_'+mech_name)/gnabarTotal_AIS
        proximal_Nav12_fraction = getattr(ais_sec(0), 'gnabar2'+'_'+mech_name)/gnabarTotal_AIS
        for seg in sec.allseg():
            rescale = rescaleVals[i_seg]
            i_seg+=1

            from run_settings import gl, gl_rescaling
            if gl_rescaling==True:
                gl*=rescale





            gna16 = getattr(ais_sec(0), 'gnabar'+'_'+mech_name)*rescale/rescale_ais
            gna12 =  getattr(ais_sec(0), 'gnabar2'+'_'+mech_name)*rescale/rescale_ais
            gk = gk_0*rescale
            rescale_pump = rescale_pump_relative*rescale
            INaKmax = ImaxPump*rescale_pump
            gnal = gna_leak*rescale_pump
            gkl = gk_leak*rescale_pump
            setattr(seg, 'gl'+'_'+mech_name, gl)
            setattr(seg, 'gnabar'+'_'+mech_name, gna16)
            setattr(seg, 'gnabar2'+'_'+mech_name, gna12)
            setattr(seg, 'gkbar'+'_'+mech_name, gk)
            setattr(seg, 'INaKmax'+'_'+mech_name, INaKmax)
            setattr(seg, 'gnal'+'_'+mech_name, gnal)
            setattr(seg, 'gkl'+'_'+mech_name, gkl)

def linear_Nav_matcher(hill_sec, ais_sec, gamma=gamma):
        from mech_settings import mech_name_hill as mech_name
        import numpy as np
        rescale_ais = gamma*rescale_soma

        sec = hill_sec
        rescaleVals = np.linspace(rescale_hill, rescale_ais, len(list(sec.allseg())), endpoint=True)
        rescaleVals_NaV16 = np.linspace(0.0, rescale_ais, len(list(sec.allseg())), endpoint=True)

        i_seg = 0
        gnabarTotal_AIS = getattr(ais_sec(0), 'gnabar'+'_'+mech_name) + getattr(ais_sec(0), 'gnabar2'+'_'+mech_name)
        proximal_Nav16_fraction = getattr(ais_sec(0), 'gnabar'+'_'+mech_name)/gnabarTotal_AIS
        proximal_Nav12_fraction = getattr(ais_sec(0), 'gnabar2'+'_'+mech_name)/gnabarTotal_AIS
        for seg in sec.allseg():
            rescale = rescaleVals[i_seg]
            rescale_NaV16 = rescaleVals_NaV16[i_seg]
            i_seg+=1

            from run_settings import gl, gl_rescaling
            if gl_rescaling==True:
                gl*=rescale





            gna16 = getattr(ais_sec(0), 'gnabar'+'_'+mech_name)*rescale_NaV16/rescale_ais
            gna12 =  getattr(ais_sec(0), 'gnabar2'+'_'+mech_name)*rescale/rescale_ais
            gk = gk_0*rescale
            rescale_pump = rescale_pump_relative*rescale
            INaKmax = ImaxPump*rescale_pump
            gnal = gna_leak*rescale_pump
            gkl = gk_leak*rescale_pump
            setattr(seg, 'gl'+'_'+mech_name, gl)
            setattr(seg, 'gnabar'+'_'+mech_name, gna16)
            setattr(seg, 'gnabar2'+'_'+mech_name, gna12)
            setattr(seg, 'gkbar'+'_'+mech_name, gk)
            setattr(seg, 'INaKmax'+'_'+mech_name, INaKmax)
            setattr(seg, 'gnal'+'_'+mech_name, gnal)
            setattr(seg, 'gkl'+'_'+mech_name, gkl)

def flat_Nav_matcher(hill_sec, ais_sec, gamma=gamma):
        from mech_settings import mech_name_hill as mech_name
        import numpy as np
        rescale_ais = gamma*rescale_soma

        sec = hill_sec
        rescaleVals = 0.8*rescale_ais*np.ones(len(list(sec.allseg())))
        print(rescaleVals/rescale_ais)
        print(rescale_hill)
        print(rescale_soma)
        # exit()
        i_seg = 0
        gnabarTotal_AIS = getattr(ais_sec(0), 'gnabar'+'_'+mech_name) + getattr(ais_sec(0), 'gnabar2'+'_'+mech_name)
        proximal_Nav16_fraction = getattr(ais_sec(0), 'gnabar'+'_'+mech_name)/gnabarTotal_AIS
        proximal_Nav12_fraction = getattr(ais_sec(0), 'gnabar2'+'_'+mech_name)/gnabarTotal_AIS
        for seg in sec.allseg():
            rescale = rescaleVals[i_seg]
            i_seg+=1

            from run_settings import gl, gl_rescaling
            if gl_rescaling==True:
                gl*=rescale





            gna16 = getattr(ais_sec(0), 'gnabar'+'_'+mech_name)*rescale/rescale_ais
            gna12 =  getattr(ais_sec(0), 'gnabar2'+'_'+mech_name)*rescale/rescale_ais
            gk = gk_0*rescale
            rescale_pump = rescale_pump_relative*rescale
            INaKmax = ImaxPump*rescale_pump
            gnal = gna_leak*rescale_pump
            gkl = gk_leak*rescale_pump
            setattr(seg, 'gl'+'_'+mech_name, gl)
            setattr(seg, 'gnabar'+'_'+mech_name, gna16)
            setattr(seg, 'gnabar2'+'_'+mech_name, gna12)
            setattr(seg, 'gkbar'+'_'+mech_name, gk)
            setattr(seg, 'INaKmax'+'_'+mech_name, INaKmax)
            setattr(seg, 'gnal'+'_'+mech_name, gnal)
            setattr(seg, 'gkl'+'_'+mech_name, gkl)



def all_one_Nav():
    from mech_settings import mech_name_hill as mech_name
    import numpy as np
    from neuron import h

    for sec in h.allsec():
        for seg in sec.allseg():

            gnabarTotal = getattr(seg, 'gnabar'+'_'+mech_name) + getattr(seg, 'gnabar2'+'_'+mech_name)

            gna16 = gnabarTotal
            gna12 = 0.0

            setattr(seg, 'gnabar'+'_'+mech_name, gna16)
            setattr(seg, 'gnabar2'+'_'+mech_name, gna12)


   

def switch_NaVs(sections):
    import numpy as np
    from neuron import h

    for sec in sections:
        for seg in sec.allseg():

            gna16 = getattr(seg, 'gnabar'+'_'+mech_name)
            gna12 = getattr(seg, 'gnabar2'+'_'+mech_name)
            gnabarTotal = gna16+gna12

            new_gna16 = gna12
            new_gna12 = gna16

            setattr(seg, 'gnabar'+'_'+mech_name, new_gna16)
            setattr(seg, 'gnabar2'+'_'+mech_name, new_gna12)


def connect_cell(h, connect_these):
    import numpy as np
    # connect_these = list(h.allsec())
    for i in np.arange(1, len(connect_these)):
        connect_these[i].connect(connect_these[i-1](1), 0)

    return connect_these


def LeftShift_a_section(sec, mech_name, vLeftShift, AC=1.0, timeLS=200.0):
    for seg in sec.allseg():
        setattr(seg, 'timeLS_'+mech_name, timeLS)

        setattr(seg, 'vLeftShift_'+mech_name, vLeftShift)
        setattr(seg, 'AC_'+mech_name, AC)

        setattr(seg, 'vLeftShift2_'+mech_name, vLeftShift)
        setattr(seg, 'AC2_'+mech_name, AC)

def LeftShift_Sections(sections, mech_name, vLeftShift, AC=1.0, timeLS=200.0):
    from neuron import h
    for sec in sections:
        for seg in sec.allseg():
            setattr(seg, 'timeLS_'+mech_name, timeLS)
            
            setattr(seg, 'vLeftShift_'+mech_name, vLeftShift)
            setattr(seg, 'AC_'+mech_name, AC)

            setattr(seg, 'vLeftShift2_'+mech_name, vLeftShift)
            setattr(seg, 'AC2_'+mech_name, AC)

    print('⁍Temperature = '+str(h.celsius)+'℃')
    print('⁍vLeftShift = ', vLeftShift)
    print('⁍AC = ', AC)
    print('⁍timeLS = ', timeLS)
    print('⁍Left-shifted sections: ', sections)


def RightShift_Nav12(h, mech_name, vRS0):
    for sec in h.allsec():
        for seg in sec.allseg():
            setattr(seg, 'vRS0_'+mech_name, vRS0)
    print('⁍vRS0 = ', vRS0)

def DeMyelinate(myelin_list):
    for internode in myelin_list[-9:-3]:
        internode.cm = 1.0



