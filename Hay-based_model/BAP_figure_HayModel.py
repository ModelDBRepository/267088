from HayClass import Pyramidal, initialize


cell1 = Pyramidal(1)
cell2 = Pyramidal(2)
initialize(warm_up_and_save=True) 


### Axonal theshold
cell1.stim(location=cell1.ais_segments[-1], duration=1.0) #<---1ms axonal pulse
cell1.plot_peak_everywhere(separation=1.0, crossover=0.5, amplitude=0.5802063404158035, plot_tpeak=False, filename='BAP_axonalHay_suprathresh') #<---1ms axonal pulse
cell1.plot_peak_everywhere(separation=1.0, crossover=0.5, amplitude=0.58, plot_tpeak=False, filename='BAP_axonalHay_subthresh') #<---1ms axonal pulse


### Somatic threshold
cell2.stim(location=cell2.soma(0.5), duration=1.0) #<---1ms somatic pulse
cell2.plot_peak_everywhere(separation=1.0, crossover=0.5, amplitude=4.7277, plot_tpeak=False, filename='BAP_somaticHay_suprathresh') #<---1ms axonal pulse
cell2.plot_peak_everywhere(separation=1.0, crossover=0.5, amplitude=4.7276, plot_tpeak=False, filename='BAP_somaticHay_subthresh') #<---1ms axonal pulse

# ###Vrest
# cell2.plot_peak_everywhere(separation=1.0, crossover=0.5, amplitude=0.0, plot_tpeak=False, filename='BAP_somaticHay_vrest') #<---1ms axonal pulse
