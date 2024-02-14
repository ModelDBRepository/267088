def cell_length(h):
    print('WARNING: cell_length() ONLY WORKS FOR TOY MODEL (RADIALLY SYMMETRIC MORPHOLOGY)')
    length = 0.0
    for sec in h.allsec():
        # print(sec.L)
        length+=sec.L
    # print('cell length = ', length)
    return length 
    


def verify_stim_delay_was_initialized_everywhere(h, mech_name):
    from run_settings import equilm_dur, resume_dur, total_dur, stimON
    for sec in h.allsec():
        for seg in sec.allseg():
            seg_stimON = getattr(seg, 'stimON_'+mech_name)
            if seg_stimON ==-1.0 or seg_stimON!=stimON:
                print('ERROR: stimulation delay (stimON) not properly passed to NMODL mechanism: spike detection will not be reliable')
                exit()



def exit_times(h, dend, soma, hill, ais, aisNeighbour):
    from sorcery import dict_of
    texit_selector = section_exitTimes
    texit_dend0, texit_dend1 = texit_selector(h, dend)
    texit_soma0, texit_soma1 = texit_selector(h, soma)
    texit_hill0, texit_hill1 = texit_selector(h, hill)
    texit_ais0, texit_ais1 = texit_selector(h, ais)
    texit_aisNeighbour0, texit_aisNeighbour1 = texit_selector(h, aisNeighbour)

    return dict_of(texit_dend0, texit_dend1, texit_soma0, texit_soma1, texit_hill0, texit_hill1, texit_ais0, texit_ais1, texit_aisNeighbour0, texit_aisNeighbour1)


def ZeroCrossingTimes(h, dend, soma, hill, ais, endNode):
    from sorcery import dict_of
    from numpy import inf
    length = cell_length(h)
    tspike_selector = section_min_tZeroCrossing

    tspike_dend, xspike_dend = tspike_selector(h, dend)
    tspike_soma, xspike_soma = tspike_selector(h, soma)
    tspike_hill, xspike_hill = tspike_selector(h, hill)
    tspike_ais, xspike_ais = tspike_selector(h, ais)
    tspike_endNode, xspike_endNode = tspike_selector(h, endNode)

    return dict_of(tspike_dend, xspike_dend, tspike_soma, xspike_soma, tspike_hill, xspike_hill, tspike_ais, xspike_ais, tspike_endNode, xspike_endNode)




def section_min_tZeroCrossing(h, section):
    from numpy import inf
    from run_settings import equilm_dur, resume_dur, total_dur, stimON
    from mech_settings import mech_name
    verify_stim_delay_was_initialized_everywhere(h, mech_name)

    tspikeMin = inf
    x_tspikeMin = inf
    spike_detected=False

    spikeTimes = [inf]
    spikeSegPositions = [inf]
    spikeSegnames = ['placeholder']
    for seg in section.allseg():
        spike = getattr(seg, 'spike_'+mech_name)
        if spike==1.0:
            spike_detected=True
            tspike = getattr(seg, 'tspike_'+mech_name)
            x = seg.x
            if tspike < min(spikeTimes):
                tspikeMin = tspike
                x_tspikeMin = x

            spikeTimes.append(tspike)
            spikeSegPositions.append(x)
            spikeSegnames.append(str(seg))
            # print('tspike = ', tspike)
            

    return tspikeMin, x_tspikeMin


def classify_grid_point(h, dend, soma, hill, ais, endNode):
    from run_settings import equilm_dur, total_dur
    from sorcery import dict_of
    from numpy import inf
    length = cell_length(h)
    tpeak_selector = section_min_tpeak
    tpeak_dend = tpeak_selector(h, dend)
    tpeak_soma = tpeak_selector(h, soma)
    tpeak_hill = tpeak_selector(h, hill)
    tpeak_ais = tpeak_selector(h, ais)
    tpeak_endNode = tpeak_selector(h, endNode)

    
    action_potential=False
    delta_t = inf
    propagation_type='undetermined'
    if equilm_dur<tpeak_endNode<total_dur:
        action_potential=True

    if (equilm_dur<tpeak_soma<total_dur) and (equilm_dur<tpeak_ais<total_dur):
        delta_t = (tpeak_soma - tpeak_ais)
        if delta_t>0.0:
            propagation_type = 'back-prop'
        elif delta_t==0.0:
            propagation_type = 'simultaneous-soma-ais'
        elif delta_t<0.0:
            propagation_type = 'Z-prop'


    print('______________classify_grid_point_______________')
    print('action_potential = ', action_potential)
    print('delta_t = ', delta_t)
    print('propagation_type = ', propagation_type)
    return dict_of(action_potential, delta_t, propagation_type, tpeak_dend, tpeak_soma, tpeak_hill, tpeak_ais, tpeak_endNode)




def section_exitTimes(h, section):
    from numpy import inf
    ###records dictionary of post-stimulation spike data. keys are the section names.
    from run_settings import equilm_dur, resume_dur, total_dur, stimON
    from mech_settings import mech_name
    verify_stim_delay_was_initialized_everywhere(h, mech_name)

    exit_times = []
    for x in [0, 1]:
        spike_detected=False
        seg = section(x)
        spike = getattr(seg, 'spike_'+mech_name)
        if spike==1.0:
            spike_detected=True
            # tspike = getattr(seg, 'tspike_'+mech_name)
            tpeak = getattr(seg, 'tpeak_'+mech_name)
            
            tspike = tpeak
        elif spike_detected==False:
            tspike = inf
        exit_times.append(tspike)
    return exit_times[0], exit_times[1]

def section_min_tpeak(h, section):
    from numpy import inf
    from run_settings import equilm_dur, resume_dur, total_dur, stimON
    from mech_settings import mech_name
    verify_stim_delay_was_initialized_everywhere(h, mech_name)
    spike_detected=False
    spikeTimes = []
    spikeSegnames = []
    for seg in section.allseg():
        spike = getattr(seg, 'spike_'+mech_name)
        if spike==1.0:
            spike_detected=True
            tpeak = getattr(seg, 'tpeak_'+mech_name)
            spikeTimes.append(tpeak)
            spikeSegnames.append(str(seg))
            # print('tpeak = ', tpeak)

    if spike_detected:
        min_tpeak = min(spikeTimes)
    else:
        min_tpeak = inf
    return min_tpeak
    

def section_report(h, section):
    from run_settings import equilm_dur, resume_dur, total_dur, stimON
    from mech_settings import mech_name


    verify_stim_delay_was_initialized_everywhere(h, mech_name)


    spike_detected=False
    spikeTimes = []
    spikeSegnames = []
    for seg in section.allseg():
        spike = getattr(seg, 'spike_'+mech_name)
        tspike = getattr(seg, 'tspike_'+mech_name)
        vpeak = getattr(seg, 'vpeak_'+mech_name)
        tpeak = getattr(seg, 'tpeak_'+mech_name)
        tup = getattr(seg, 'tup_'+mech_name)
        tdown = getattr(seg, 'tdown_'+mech_name)
        # print('spike = ', spike)
        # print('tspike = ', tspike)
        # print('vpeak = ', vpeak)
        print('tpeak = ', tpeak)
        # print('tup = ', tup)
        # print('tdown = ', tdown)



def section_spikeTimes(h, section):
    ###records dictionary of post-stimulation spike data. keys are the section names.
    from numpy import inf
    from run_settings import equilm_dur, resume_dur, total_dur, stimON
    from mech_settings import mech_name

    verify_stim_delay_was_initialized_everywhere(h, mech_name)

    spikeData_dict = {}
    min_tspike_dict = {}

    spike_detected=False
    spikeTimes = []
    spikeSegnames = []
    for seg in section.allseg():
        spike = getattr(seg, 'spike_'+mech_name)
        if spike==1.0:
            spike_detected=True
            # tspike = getattr(seg, 'tspike_'+mech_name)
            tpeak = getattr(seg, 'tpeak_'+mech_name)
            tspike = tpeak
            spikeTimes.append(tspike)
            spikeSegnames.append(str(seg))

    if spike_detected:
        min_tspike = min(spikeTimes)
    else:
        min_tspike = inf
    spikeData_dict[str(section)] = dict(min_tspike=min_tspike, spikeTimes=spikeTimes, spikeSegnames=spikeSegnames)
    min_tspike_dict[str(section)] = min_tspike
    # print('spikeData_dict = ', spikeData_dict)
    return spikeData_dict, min_tspike_dict




def return_spikeTimes(h, sections):
    ###records dictionary of post-stimulation spike data. keys are the section names.
    from numpy import inf
    from run_settings import equilm_dur, resume_dur, total_dur, stimON
    from mech_settings import mech_name

    verify_stim_delay_was_initialized_everywhere(h, mech_name)

    spikeData_dict = {}
    min_tspike_dict = {}
    for sec in sections:
        spike_detected=False
        spikeTimes = []
        spikeSegnames = []
        for seg in sec.allseg():
            spike = getattr(seg, 'spike_'+mech_name)
            if spike==1.0:
                spike_detected=True
                # tspike = getattr(seg, 'tspike_'+mech_name)
                tpeak = getattr(seg, 'tpeak_'+mech_name)
                tspike = tpeak
                spikeTimes.append(tspike)
                spikeSegnames.append(str(seg))

        if spike_detected:
            min_tspike = min(spikeTimes)            
        else:
            min_tspike = inf
        spikeData_dict[str(sec)] = dict(min_tspike=min_tspike, spikeTimes=spikeTimes, spikeSegnames=spikeSegnames)
        min_tspike_dict[str(sec)] = min_tspike
    # print('spikeData_dict = ', spikeData_dict)
    return spikeData_dict, min_tspike_dict
