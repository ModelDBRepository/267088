from plotHay import DisplayHayData
import Hay_data as hd

###____________________________________________________________________________________
"""Plotting effects of Nav distribution on backpropagation threshold""" 

###Fig 7b
pAxonal = DisplayHayData(hd.axonalStim_negative70mV_1msPulse_new, filename='Fig7B_')
pAxonal.plot_thresh(select_kappa_vals=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
# pAxonal.plot_thresh(select_kappa_vals=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], use_ratio_of_total_AIS_NaV_conductances_from_each_subtype=True)



