# Impact on backpropagation of the spatial heterogeneity of sodium channel kinetics in the axon initial segment

## The manuscript associated with this code is in production at PLOS Computational Biology. 


## Click the "Files" tab to see the source code.
* Python scripts whose names match those of figures from the paper regenerate said figures from saved data. 
Simply enter, e.g., `python Fig2.py` in the terminal. Please use a recent version of Python 3.


* If there are any issues using the code please contact BBarl039@uOttawa.ca. 


## This ModelDB page includes two models, in separate folders:

1. The __Hay-based__ model, which was built using the complete biophysics and morphology from [Hay et al. (2011)](https://modeldb.science/139653), with an added section of passive cable attached at the end of the AIS, and 

2. The __Hu-based__ model, which was built using the Na<sub>V</sub>1.2, Na<sub>V</sub>1.6, and K<sub>V</sub> channel models, and pyramidal cell morphology from [Hu et al. (2009)](https://modeldb.science/123897). 

* The simulations compute backpropagation and forwardpropagation (AP) thresholds with a range of AIS Na<sub>V</sub> distributions, from scratch, albeit slowly in the case of the Hu-based model (hours) and more rapidly in the Hay-based model.


### In the Hay-based model folder, the user can run the following simulations:

* `python backprop_somatic_stim.py`

    This script generates the data for Fig 7A.


* `python backprop_axonal_stim.py`

    This script generates the data for Fig 7B.


* `python BAP_figure_HayModel.py`

    This script generates space plots of backpropagation with axonal and somatic stimulation, just above and just below threshold. It should create four PDFs in total.
 


### The Hu-based model folder has four subfolders in which the user can run simulations:
(Inside a given folder, enter the corresponding command into the terminal.)

* *__backprop_threshold_somatic_stimulation__*:

    `python backprop_somatic_stim.py`

* *__backprop_threshold_axonal_stimulation__*:

    `python backprop_axonal_stim.py`

* *__reference_curve_somatic_stimulation__*:

    `python reference_curve_somatic.py`

* *__reference_curve_axonal_stimulation__*:

    `python reference_curve_axonal.py`


* Running the Hu-based model code will produce CSV files containing the threshold and various model settings, such as the NaV distribution parameters 'x' and '𝜅' (called 'deviation' and 'CrossOverPosition', respectively, in the CSV files). 
To plot the data in files with names beginning "`sqlite3`...", they must first be processed into Python dictionaries using the script `database_to_dictionary.py`.  



