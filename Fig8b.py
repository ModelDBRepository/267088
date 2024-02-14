from plotHay import DisplayHayData
import Hay_data as hd

###____________________________________________________________________________________
"""ploting 'scaleNavTwo' results"""

###Fig 8b
data = hd.axonalStim_scaleNavTwo__sweepKappa
data_lowend = hd.lowend_axonalStim_scaleNaVTwo
for key, value in data.items():
    data[key].extend(data_lowend[key])
pAxonal = DisplayHayData(data, xAxisVariableName='scaleNavTwo', filename='Fig8B_')
pAxonal.plot_thresh(legend_location='lower right', vline=0.32, ylim=(0.57, 0.606), xlim=(-0.03, None))



