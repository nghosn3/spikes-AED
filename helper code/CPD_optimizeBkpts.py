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
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
  
def test_bkpts(signal, max_bkpts,plot=True):
    """Gets the tuning curve of variance of data segments vs the number of segments 
    ----------
    Signal : full signal, unsegmenred (1D)
    max_bkpts : highest number of breakpoints to test
    -------
    returns:
    bkpt_curve : 2 x max_bpts array containing the tuning curve. should be the input of get_nBkpts()
    """
    bkpt_curve = np.zeros((2,max_bkpts-1))
    for i in range(1,max_bkpts-1):
            
        algo = rpt.KernelCPD(kernel="linear").fit(signal)
        result = algo.predict(n_bkps=i)
        result.insert(0,0)
        result.insert(-1,len(signal)-1)
        # compute the sum of variances between all segments 
        var = 0;
        for ii in range(0,len(result)-1):
            seg = signal[result[ii]:result[ii+1]]
            var =var + np.var(seg)
            
        bkpt_curve[0,i] = i;
        bkpt_curve[1,i] = var/i;
    
    if plot:
        plt.plot(bkpt_curve[0,:], bkpt_curve[1,:])
        plt.show()
    return bkpt_curve

def get_nBkpts(bkpts_curve):
    """Gets the tuning curve of variance of data segments vs the number of segments 
    ----------
     bkpt_curve : 2 x max_bpts array containing the tuning curve. 
    ----------
    returns: 
        n_bkpts: the 'optimal' number of breakpoints at the elbow point of the tuning curve 
    """
    bkpts_curve = bkpts_curve[:,1:]
    x = bkpts_curve[0,:]
    y =  bkpts_curve[1,:]
    errors = []
    for i in range(2,len(x)-2):
        y1 = y[0:i].reshape(-1, 1);
        y2 = y[i+1:-1].reshape(-1, 1);

        x1 = x[0:i].reshape(-1, 1);
        x2 = x[i+1:-1].reshape(-1, 1);  
        
        #do linear regression for two segments of curve
        reg1 = LinearRegression().fit(x1,y1)
        reg2 = LinearRegression().fit(x2,y2)
        
        y1_pred = reg1.predict(x1);
        y2_pred = reg2.predict(x2);
        
        #find the mean squared error of the line for each segment and add
        res1 =mean_squared_error(y1, y1_pred);
        res2 =mean_squared_error(y2, y2_pred) ;
        
        errors.append(res1 +res2);
    
    # find the index (i) with the smallest error, that should be the n_bkpts
    n_bkpts = np.argmin(errors) + 2

    
    return n_bkpts


#################################################################################

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
    max_bkpts =20;
    bkpt_curve = test_bkpts(signal,max_bkpts,plot=True)
    n_bkpts = get_nBkpts(bkpt_curve)  # number of breakpoints 
    
    # detection 
    algo = rpt.KernelCPD(kernel="linear").fit(signal)
    result = algo.predict(n_bkps=n_bkpts)
    
    # display
    rpt.display(signal, result)
    plt.show()
    tle = "number of break points = {n_bkpts}".format(n_bkpts = n_bkpts)
    plt.title(tle)
    
    
  
