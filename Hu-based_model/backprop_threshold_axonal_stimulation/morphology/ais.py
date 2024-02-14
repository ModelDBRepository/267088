def makeAIS(h): 
    from run_settings import diam_ais, L_ais
    sec = h.Section(name='ais')

    sec.L = L_ais #100.0
    sec.diam = diam_ais #1.6 
    sec.Ra = 150
    sec.cm = 1.0 
    sec.nseg = 221


    ###__________________________________________________________________________
    from mech_settings import mech_name_ais as mech_name
    from mech_settings import rescale_ais as rescale
    from mech_settings import insertCLS
    insertCLS(sec, mech_name, rescale)


    return sec



def plot_NaV_channels(sections):
    from BBplotting import BarlowRangeVarPlotter
    from mech_settings import mech_name_ais
    BarlowRangeVarPlotter(['gnabar_'+mech_name_ais, 'gnabar2_'+mech_name_ais], sections=sections, ylims=[0.0, None], save_PNG=False)

def plot_ais_channels(ais):
    from BBplotting import BarlowRangeVarPlotter
    from mech_settings import mech_name_ais
    BarlowRangeVarPlotter(['gnabar_'+mech_name_ais, 'gnabar2_'+mech_name_ais], sections=[ais], ylims=[0.0, None], save_PNG=False)


def diminish_AIS_Nav12(sec, fraction):
    from mech_settings import mech_name_ais as mech_name
    for seg in sec.allseg():
        gnabar2 = getattr(seg, 'gnabar2_'+mech_name)
        gnabar2 *= fraction
        setattr(seg, 'gnabar2_'+mech_name, gnabar2)
    return gnabar2



def ais_tanh_distribution(sec, gnabarTotal=0.12, flatness=5.0, deviation=1.0, FlipNavs_LeftRight=False, do_nothing=False, CrossOverPosition=0.5):
    ### MUST FOLLOW h.finitialize()
    from mech_settings import mech_name_ais as mech_name
    import numpy as np
    if do_nothing==False:
        print('⁍Implementing tanh-distribution in ais')
        if np.absolute(deviation)>1.0:
            print('ERROR: deviation must not exceed 1.0 !!!!')
            exit()
        elif deviation <0.0:
            print('WARNING: negative deviation value: flips ais Nav channel profiles left-to-right!!! (If this is desired, use FlipNavs_LeftRight=True to do it correctly.)')
            exit()
        if FlipNavs_LeftRight==True:
                print('⁍Flipping the AIS Nav profiles left-to-right: Nav16 now peaks on the LEFT, and Nav12 peaks on the RIGHT.')
                CrossOverPosition = 1.0 - CrossOverPosition
        i=0
        for seg in sec.allseg():
            if gnabarTotal <= 0.0:
                print('ERROR, gnabarTotal='+str(gnabarTotal))
                exit()

            x = seg.x
            if FlipNavs_LeftRight==True:
                x = (1.0 - x)

            
            gnabar2 = 0.5*gnabarTotal*(   1.0 - deviation*np.tanh(flatness*(x-CrossOverPosition))   )
            gnabar = gnabarTotal - gnabar2 #gnabarTotal*((1.0 + np.tanh(flatness*(x-0.5)))/2.0)
            setattr(seg, 'gnabar_'+mech_name, gnabar)
            setattr(seg, 'gnabar2_'+mech_name, gnabar2)
            
            i+=1
        
    # print('sec(0.3).gnabar + sec(0.3).gnabar2 = ', getattr(sec(0.3), 'gnabar_'+mech_name)+getattr(sec(0.3), 'gnabar2_'+mech_name))
    return getattr(sec(0.3), 'gnabar_'+mech_name)+getattr(sec(0.3), 'gnabar2_'+mech_name)









