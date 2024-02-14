def squareCurrentPulse(stimLoc, delay=10, amplitude=10, duration=5):
    from neuron import h 
    stim = h.IClamp(stimLoc)
    stim.delay = delay
    stim.amp = amplitude
    stim.dur = duration
    return stim, stimLoc