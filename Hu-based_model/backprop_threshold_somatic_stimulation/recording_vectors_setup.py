def setup_recording_vectors(h, recording_segments):
    from mech_settings import mech_name
    from run_settings import chloride
    data_dict = {}

    varnames= ['v', 'ena', 'ek', 'nai', 'nao', 'ki', 'ko']
    if chloride==True:
        varnames= ['v', 'ena', 'ek', 'ecl', 'nai', 'nao', 'ki', 'ko', 'cli', 'clo']
        

    mech_varnames = ['m', 'h', 'n', 'mLS', 'hLS', 'nLS', 'm2', 'h2', 'mLS2', 'hLS2', 'I_Nav', 'I_Nav2']
    for i in range(len(mech_varnames)):
        mech_varnames[i] += '_'+mech_name
    varnames+=mech_varnames

    for seg in recording_segments:
        data_dict[str(seg)] = seg_timeseries(h, seg, variable_names=varnames)
    return data_dict


def seg_timeseries(h, recording_segment, variable_names=['v', 'ena', 'ek', 'nai', 'nao', 'ki', 'ko']):
    vectors = {}
    t = h.Vector().record(h._ref_t)
    vectors['t'] = t
    for var in variable_names:
        # vectors[var+'_'+str(recording_segment)] = h.Vector().record(getattr(recording_segment, '_ref_'+var))
        vectors[var] = h.Vector().record(getattr(recording_segment, '_ref_'+var))
    return vectors


def seg_mech_timeseries(h, recording_segment, varname='vLS'):
    ###Grabs timeseries for NMODL mechanism variable 'varname' at this segment
    from mech_settings import mech_name
    vector = h.Vector().record(getattr(recording_segment, '_ref_'+varname+'_'+mech_name))
    #getattr(seg, varname+'_'+mech_name)       
    return vector


