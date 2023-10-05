#%%

import os
import sys
sys.path.append('/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/ieegpy/ieeg')
sys.path.append('tools')

# sets path to one directory up from where code is
path = "/".join(os.path.abspath(os.getcwd()).split('/')[:-1])

import json
import numpy as np
from get_iEEG_data import get_iEEG_data
from plot_iEEG_data import plot_iEEG_data

# %%
with open("../credentials.json") as f:
    credentials = json.load(f)
    username = credentials['username']
    password = credentials['password']

iEEG_filename = "HUP166_phaseII" #test on HUP166 first 
start_time_usec = 402580 * 1e6
stop_time_usec = 402800 * 1e6
electrodes = ["LE10","LE11","LH01","LH02","LH03","LH04"]

# %%
data, fs = get_iEEG_data(username, password, iEEG_filename, start_time_usec, stop_time_usec, select_electrodes=electrodes)

# %% Plot the data
t_sec = np.linspace(start_time_usec, stop_time_usec, num=data.shape[0]) / 1e6
fig, ax = plot_iEEG_data(data, t_sec)
fig.set_size_inches(18.5, 10.5)
ax.set_title(iEEG_filename)
fig.show()

# %%
