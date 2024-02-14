"""
Export simulation data contained in SQL files to human-readable (and plotable) data_dictionary.py files. 
The dictionaries in these python files can then be imported as "data" in, e.g., Fig6a.py.
That is, if the python file that is created by this program is data_dict_2024_02_09__23_22_26.py, 
one can write: 
from data_dict_2024_02_09__23_22_26 import data
at the top of Fig6a.py
"""
NavTwo_rightshift=13.0 #mV

data_output = {}
data_output['stimLoc'] = []
data_output['AIS_NaV_subtype_separation'] = []
data_output['CrossOverPosition'] = []
data_output['activate_DeltaVrs_minf2__inAIS'] = []
data_output['activate_DeltaVrs_mtau2__inAIS'] = []
data_output['activate_DeltaVrs_hinf2__inAIS'] = []
data_output['activate_DeltaVrs_htau2__inAIS'] = []
data_output['DeltaVrs'] = []
data_output['Qthresh'] = []

import numpy as np
from sorcery import print_args
import sqlite3
import pandas as pd
import pprint
from datetime import datetime

DB_FILENAME_list = []
# DB_FILENAME_list += ['sqlite3pandas_data_EXAMPLE_FILENAME_1.csv']
# DB_FILENAME_list += ['sqlite3pandas_data_EXAMPLE_FILENAME_2.csv']
# DB_FILENAME_list += ['sqlite3pandas_data_EXAMPLE_FILENAME_3.csv']  #(and so on...)
DB_FILENAME_list += ['Hu-based_model/reference_curve_somatic_stimulation/sqlite3pandas_data_Lais=25.0_2024-02-0918-39-17.sqlite']


for DB_FILENAME in DB_FILENAME_list:
    with sqlite3.connect(DB_FILENAME) as conn:
        data = pd.read_sql("SELECT * FROM data", conn)

    # data.to_csv(".csv", index=False)
    # data.to_csv(DB_FILENAME+"_data.csv", index=False)

    stimLoc = data['stimLoc'][0]
    deviation = data['deviation'][0]
    vThreshold = data['vThreshold'][0]
    CrossOverPosition = data['CrossOverPosition'][0]
    disableRightShift_minf2__inAIS = data['disableRightShift_minf2__inAIS'][0] 
    disableRightShift_mtau2__inAIS = data['disableRightShift_mtau2__inAIS'][0] 

    disableRightShift_hinf2__inAIS = data['disableRightShift_hinf2__inAIS'][0] 
    disableRightShift_htau2__inAIS = data['disableRightShift_htau2__inAIS'][0] 
    print('-----------------------------------')
    print_args(vThreshold)
    print_args(CrossOverPosition)
    print_args(disableRightShift_minf2__inAIS)
    print_args(disableRightShift_mtau2__inAIS)
    print_args(disableRightShift_hinf2__inAIS)
    print_args(disableRightShift_htau2__inAIS)

    data_output['stimLoc'] += [stimLoc]
    data_output['AIS_NaV_subtype_separation'] += [deviation]
    data_output['CrossOverPosition'] += [CrossOverPosition]

    ###Change "disableRightShift..." to a more meaningful name:
    data_output['activate_DeltaVrs_minf2__inAIS'] += [disableRightShift_minf2__inAIS]
    data_output['activate_DeltaVrs_mtau2__inAIS'] += [disableRightShift_mtau2__inAIS]
    data_output['activate_DeltaVrs_hinf2__inAIS'] += [disableRightShift_hinf2__inAIS]
    data_output['activate_DeltaVrs_htau2__inAIS'] += [disableRightShift_htau2__inAIS]
    print('-----------------------------------')


    ###Get the data for this curve
    data = data[(data['Qthresh']!= np.inf)]
    stimLoc = np.array(data['stimLoc'])[0]
    Qthresh = np.array(data['Qthresh'], dtype=float)
    Vpeak = np.array(data['vpeak_SDtips'], dtype=float)
    vRS = np.array(data['vRS03'] - NavTwo_rightshift, dtype=float)


    data_output['DeltaVrs'] += [list(vRS)]
    data_output['Qthresh'] += [list(Qthresh)] 


print_args(data_output)


def export_dict_to_python_script(dict_to_export):
    # Use pprint to get a nicely formatted string of the dictionary
    formatted_dict_str = pprint.pformat(dict_to_export, width=80, sort_dicts=True)
    # Generate a file name with the current date and time
    datetime_now = datetime.now().strftime("%Y_%m_%d__%H_%M_%S")
    file_name = f'data_dict_{datetime_now}.py'

    # Write to a file, including the variable assignment
    with open(file_name, 'w') as file:
        file.write(f"data = {formatted_dict_str}\n")
    print(f'Dictionary has been saved to {file_name}')
export_dict_to_python_script(data_output)