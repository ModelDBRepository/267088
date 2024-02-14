from plotHay import DisplayHayData
import Hay_data as hd

###____________________________________________________________________________________
"""ploting 'scaleNavTwo' results"""

###Fig 8a
pSomatic = DisplayHayData(hd.somaticStim_scaleNavTwo__sweepKappa, xAxisVariableName='scaleNavTwo', filename='Fig8A_')
pSomatic.plot_thresh(legend_location='lower right')




