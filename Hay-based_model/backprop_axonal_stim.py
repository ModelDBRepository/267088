from HayClass import Pyramidal, initialize, export_dict_to_python_script
from sorcery import print_args

cell1 = Pyramidal(1)
initialize(warm_up_and_save=True) 
cell1.stim(location=cell1.ais_segments[-1], duration=1.0) #<---1ms axonal pulses

###Quick test run:
# cell1.threshCurve(separation_vals=[1.0], crossover_vals=[0.8])
# cell1.threshCurve(separation_vals=[0.0, 1.0], crossover_vals=[0.2, 0.5, 0.8])

###Full simulation:
cell1.threshCurve(separation_vals=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], crossover_vals=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
thresh_data_axonal_stim = cell1.threshold_dict

print_args(thresh_data_axonal_stim)

export_dict_to_python_script(thresh_data_axonal_stim)