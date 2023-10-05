#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 15:35:28 2021

@author: ninaghosn
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 14:48:03 2021

implement change point on spike data during AED taper - start with whole emu stay

@author: ninaghosn
"""
import matplotlib.pyplot as plt
import ruptures as rpt
import scipy.io as sio
from os.path import dirname, join as pjoin
from pathlib import Path
import numpy as np

# get path of data
my_path = Path("/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA")
data_dir = pjoin(dirname(sio.__file__), my_path)
mat_fname = pjoin(data_dir, 'all_spike_rate.mat')

# load the spike rate data
mat =sio.loadmat(mat_fname)
spike_rates = mat['all_spike_rate']

# do the change point detection on all the patients in the list
for i in range(0,spike_rates.size-1):    # load signal
    signal = spike_rates[0][i];
    signal= signal.transpose();
    
    # determine the optimal number of breakpoints - try var(segments)v #segments
   
    pen = 20; #penalty value? dont know the units 
    n_bkps = 2  # number of breakpoints 

    
    # detection 

    algo = rpt.KernelCPD(kernel="linear").fit(signal)
    result = algo.predict(n_bkps=n_bkps, pen=pen)
    
    
    # display
    rpt.display(signal, result)
    plt.show()
    
    
    
    def test_bkpts(signal, max_bkpts):
    """Gets the tuning curve of variance of data segments vs the number of segments 
    ----------
    Signal : full signal, unsegmenred (1D)
    max_bkpts : highest number of breakpoints to test
    -------
    returns:
    bkpt_curve : 2 x max_bpts array containing the tuning curve. should be the input of get_nBkpts()
    """
    
    
    
    return bkpt_curve
    # other code fragments
    
    # n=len(signal)
    # sigma = np.std(signal)
    # epsilon=.00005 * n * sigma ** 2
