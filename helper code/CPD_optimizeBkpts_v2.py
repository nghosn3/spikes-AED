#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 16:53:12 2021

@author: ninaghosn
"""

import matplotlib.pyplot as plt
import ruptures as rpt
import scipy.io as sio
from scipy.io import savemat
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
    bkpt_curve = []
    for i in range(1,max_bkpts):
            
        algo = rpt.KernelCPD(kernel="linear").fit(signal)
        result = algo.predict(n_bkps=i)
        
        result.insert(0,0)
        result.insert(-1,len(signal)-1)
        
        # compute the sum of variances between all segments 
        sse = 0;
        for ii in range(0,len(result)-1):
            seg = signal[result[ii]:result[ii+1]]
            mu = np.mean(seg)
            sse_seg = np.sum((seg-mu)**2)
            sse = sse+sse_seg
            #var =var + np.var(seg)
            
        
        #bkpt_curve.append(var/i);
        bkpt_curve.append(sse);
    
    if plot:
        plt.figure()
        plt.plot(range(1,max_bkpts), bkpt_curve)
    return bkpt_curve

def get_nBkpts(bkpts_curve):
    """Gets the tuning curve of variance of data segments vs the number of segments 
    ----------
     bkpt_curve : 2 x max_bpts array containing the tuning curve. 
    ----------
    returns: 
        n_bkpts: the 'optimal' number of breakpoints at the elbow point of the tuning curve 
    """
    x = np.arange(len(bkpts_curve))
    y =  np.array(bkpts_curve)
    errors = []
    for i in range(2,len(x)-2):
        y1 = y[0:i].reshape(-1, 1);
        y2 = y[i:-1].reshape(-1, 1);

        x1 = x[0:i].reshape(-1, 1);
        x2 = x[i:-1].reshape(-1, 1);  
        
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

# get path of the cohort data
my_path = Path("/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA")
data_dir = pjoin(dirname(sio.__file__), my_path)
mat_fname = pjoin(data_dir, 'ptIDs_10_13.mat')

#load the ptIDS
mat =sio.loadmat(mat_fname)
pt_ids_raw = mat['ptIDs']
pt_ids = [str(element) for element in pt_ids_raw[0]]

# get path of data
my_path = Path("/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/DATA")
data_dir = pjoin(dirname(sio.__file__), my_path)
#mat_fname = pjoin(data_dir, 'all_spike_rate_ds.mat') 

mat_fname = pjoin(data_dir, 'all_spike_rate_new.mat') 

# load the spike rate data
mat =sio.loadmat(mat_fname)
spike_rates = mat['all_spike_rate']



for i in range(0,spike_rates.size):    # load signal
    signal = spike_rates[0][i];
    signal= signal.transpose();
    
    if signal.size>0:
    # determine the optimal number of breakpoints - try var(segments)v #segments
        max_bkpts =10;
        bkpt_curve = test_bkpts(signal,max_bkpts,plot=True)
        n_bkpts = get_nBkpts(bkpt_curve)  # number of breakpoints 
    
        # detection - use binary segmentation
        algo = rpt.KernelCPD(kernel="linear").fit(signal)
        #algo = rpt.Binseg(model='linear').fit(signal)
        result = algo.predict(n_bkps=n_bkpts)
    
        # display
        rpt.display(signal, result)
        
        # save the breakpoint indices in a matfile for plotting
        result_array = np.array(result)
        adict = {}
        adict['bkpts'] = result_array
        filename = r'/Volumes/borel.seas.upenn.edu/public/USERS/nghosn3/Pioneer/spikes-AED/results/change_point/downsamp_spikes_new/' + pt_ids[i] + "_bkpts" '.mat'
        savemat(filename, adict)   
       
        tle =  "{ID} number of break points = {n_bkpts}".format(n_bkpts = n_bkpts,ID = pt_ids[i])
        plt.title(tle)
        plt.xticks(ticks = [0, 100, 200, 300, 400],labels=["0", "50","100","150","200"])
   

    
    
    
  
