#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 20:25:19 2020

MultiModes is a code to extract all significant pulse modes 
from one or more signals. For each of them, it calculates 
the LombScargle periodogram and performs a non linear 
simultaneous fit, using a multisine function, 
to a set number of frequencies. As a stop criterion, 
you can choose between the FAP or the SNR criterion, 
depending on the type of analysis you want to perform. 

@author: David Pamos Ortega

"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.timeseries import LombScargle
from astropy.io import fits
from lmfit import Parameters, minimize
from scipy import stats
import os
import argparse
import glob
from timeit import default_timer as timer
# Initial parameters


osratio = 5 # Oversampling ratio, by default 
max_freq = 100 # Maximum frequency in the periodogram, by default
sim_fit_n = 20 # Maximum number of frequencies of the simultaneous fit, by default
stop = 'SNR'
min_snr = 4
max_fap = 0.01
tail_per = 80

# Reading the initial file with the values of the parameters, if it exists 

dash = 100*'-'
print(dash)
print('Running MultiModes'.center(110))
print(dash)

filename = 'ini.txt'
if os.path.isfile(filename):
    with open(filename, 'r') as file:
        for line in file:
            line.replace("\n", "")
            if line.startswith("sim_fit_n"):
                sim_fit_n = int(line.split(' ')[1])
            if line.startswith("osratio"):
                osratio = float(line.split(' ')[1])
            if line.startswith("maxfreq"):
                max_freq = float(line.split(' ')[1])
            if line.startswith("maxfap"):
                max_fap = float(line.split(' ')[1])  
            if line.startswith("tail_per"):
                tail_per = float(line.split(' ')[1])
            if line.startswith("minsnr"):
                min_snr = float(line.split(' ')[1])
            if line.startswith("stop"):
                stop = line.split(' ')[1]
else:
    print('Not ini.txt. Values of the parameters by default:')     




print('Number of frequencies of the simultaneous fit: ' + str(sim_fit_n))
print('Samples per peak: ' + str(osratio))
print('Maximum frequency: ' + str(max_freq))
print('Stop Criterion: ' + stop)


# Creating the necessary lists

all_best_freqs = []
all_max_amps = []
all_phs = []
all_sigma_amps = []
all_sigma_freqs = []
all_sigma_phs = []
params = Parameters()
all_faps = []
S_N = []
    
# Defining all the necessary functions

def sinusoid(t, A, f, ph): 
    ''' Sine function to fit a frequency of a light curve'''
    sin = A*np.sin(2*np.math.pi*(f*t + ph))
    return sin


def periodogram(time, flux): 
    '''Fast Lomb Scargle to calculate the periodogram'''
    ls = LombScargle(time, flux, normalization = 'psd', center_data = True, fit_mean = False)
    frequency, power = ls.autopower(method = 'fast', maximum_frequency = max_freq, samples_per_peak = osratio)
    best_frequency = frequency[np.argmax(power)]
    amps = 2*np.sqrt(power/len(time))
    noise = []
    for (f, a) in zip(frequency, amps):
        if f > tail_per:
            noise += [a]
    noise = np.mean(noise)
    return ls, frequency, amps, best_frequency, 2*np.sqrt(power.max()/len(time)), power, noise

def fit(t, params):
    '''Multi sine fit function with all the parameters of frequencies, amplitudes, phases and mean floats'''
    y = 0
    pars = params.valuesdict()
    amps_dict = {k:pars[k] for k in pars if 'a' in k}
    freqs_dict = {k:pars[k] for k in pars if 'b' in k}
    phs_dict = {k:pars[k] for k in pars if 'c' in k}
    freqs = list([pars[k] for k in freqs_dict])
    amps = list([pars[k] for k in amps_dict])
    phs = list([pars[k] for k in phs_dict])
    for (a, f, p) in zip(amps, freqs, phs):
        y += sinusoid(t, a, f, p) 
    return y, amps, freqs, phs

def residual(params, t, flux): 
    '''Residual between the model and data'''
    return flux - fit(t, params)[0]


def lightcurve(file):
    '''Reading the .fits file to extract all data'''
    hdul = fits.open(file)
    data = hdul[1].data
    time = np.array(data['col1'])
    time = time - time[0]
    T = time[-1]
    N = len(time)
    r = 1/T
    fluxes = np.array(data['col2'])
    mean_flux = np.mean(fluxes)
    fluxes = (fluxes-mean_flux)/mean_flux*1000
    return time, fluxes, T, N, r

def snr_or_fap(stop):
    '''Choosing between the SNR or the FAP stop criterion. 
       SNR by default'''
    if stop == 'SNR':
        return min_snr, S_N
    elif stop == 'FAP':
        return max_fap,  all_faps

def comb_freqs(pd):
    ''' Function to detect the combination frequencies once all of them
        having been extracted'''
    freqs = pd['Freqs']
    amps = pd['Amps']
    f1 = freqs[0]
    f2 = freqs[1]
    f_lim = 100.
    
    n1_list = np.arange(-2, 2).tolist()
    n2_list =np.arange(-2, 2).tolist()
    children1_2 = []
    
    n1_n2 = []
    
    
    for i in n1_list:
        for j in n2_list:
            candidate = i*f1 + j*f2
            if 0 < candidate < f_lim and candidate != f1 and candidate !=f2:
                children1_2.append(round(i*f1 + j*f2, 2))
                n1_n2.append((i,j,0))
         

    
    comb_freqs = []
    
    comb_amps = []
    
    n1_n2_values = []
    
    error = rayleigh
    
    for (ch, v) in zip(children1_2, n1_n2):
        for (f,a) in zip(freqs,amps):
            if abs(ch-f)<error:
                comb_freqs.append(f)
                comb_amps.append(a)
                n1_n2_values.append(v)
    
    return comb_freqs, comb_amps, n1_n2_values

    
# Creating the directory with the fits files as argument to the command line    
parser = argparse.ArgumentParser(description='Open directory with .fits files')
parser.add_argument('--d', required = 'True', type = str, help ='Select a .fits files directory to be opened')
args = parser.parse_args()
d = args.d
# Crating the list with all the files and filenames to do massive analysis

pth = './'+d+'/'
fits_files = [f for f in glob.glob(pth+'*.fits')]
fits_names = [os.path.basename(f) for f in glob.iglob(pth+'*.fits')]
fits_names = [os.path.splitext(f)[0] for f in fits_names]
fits_names = sorted(fits_names)
fits_files = sorted(fits_files)


# Starting with the massive analysis


columns = ['Number',
           'f (c/d)',
           'A (mmag)',
           stop]

# Setting stop criterion

parade = snr_or_fap(stop)[0]
snr_or_faps = snr_or_fap(stop)[1]

for (f, nm) in zip(fits_files, fits_names):
    start = timer()
    n = 1
    num = 1
    newpath = './results/' + nm + '/'
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    data = lightcurve(f) # Extracting data of every light curve
    time = data[0]
    lc = data[1] # Light curve
    T = data[2] # Total time of sampling
    N = data[3] # Number of points of the light curve
    rayleigh = data[4] # Rayleigh resolution of the sampling
    sigma_lc = stats.sem(list(lc)) # 1-sigma error of the magnitude
    sigma_amp = np.sqrt(2/N)*sigma_lc # 1-sigma error of the amplitude
    lc_df = pd.DataFrame({'Time (d)':time, 'Amplitude (mmag)':list(lc)})
    lc_df.to_csv(newpath+'lc.dat', sep = ' ', index=False, header = None)
    lc_df.plot(kind='scatter', x='Time (d)', y = 'Amplitude (mmag)', color='blue', s = 1, title=nm)
    plt.savefig(newpath+'LC.png')
    plt.show()
    lc0 = lc
    ls0 = periodogram(time, lc0)
    periodograms = [ls0,] # Calculating the initial periodogram
    n_per = [0,]
    per = pd.DataFrame({'Frequency (c/d)':ls0[1], 'Amplitude (mmag)':ls0[2]})
    per.to_csv(newpath+'pg.dat', sep=' ', index=False, header = None)
    per.plot(kind = 'line', x='Frequency (c/d)', y='Amplitude (mmag)', title = nm, legend = False)
    plt.xlabel('Frequency (c/d)')
    plt.ylabel('Amplitude (mmag)')
    plt.savefig(newpath+'LS.png')
    plt.show()
    # Initialization of the lists to save the extracted frequencies, amplitudes and phases
    all_best_freqs = []
    all_max_amps = []
    all_phs = []
    snr_or_faps = []
    all_sigma_amps = []
    all_sigma_freqs = []
    all_sigma_phs = []
    params = Parameters()
    print(dash)
    print('{:<15s}{:>18s}{:>25s}{:>32}' .format(columns[0],
                                                columns[1],
                                                columns[2],
                                                columns[3]))
    print(dash)
    while n <= sim_fit_n:
        ls = periodogram(time, lc)
        freq = ls[3]
        amp = ls[4]
        sigma_freq = np.sqrt(6/N)/(np.math.pi*T)*sigma_lc/amp
        ph = 0.5
        snr = amp/ls[6]
        var = len(time)*(ls[6])**2/4
        fap = 1 - (1 - np.exp(-ls[5].max()/var))**(N/2)
        if parade == min_snr:
            if snr > parade: # Stop criterion
                snr_or_faps.append(snr)
                all_sigma_amps.append(sigma_amp)
                new_guesses = [amp, freq, ph]
                params.add('p_'+str(n)+'a', value = new_guesses[0], min=0, max=2*amp)
                params.add('p_'+str(n)+'b', value = new_guesses[1], min=freq-rayleigh/2, max=freq+rayleigh/2)
                params.add('p_'+str(n)+'c', value = new_guesses[2], min=0, max=1)
                #best_freqs = fit(time, params)[2]
                max_amps = fit(time, params)[1]
                res = minimize(residual, params, args=(time, lc0), method = 'least_squares')
                lc = res.residual
                params = res.params
                sigma_freqs = [np.sqrt(6/N)/(np.math.pi*T)*sigma_lc/np.abs(a) for a in max_amps]
                sigma_freq = np.mean(sigma_freqs)
                sigma_phs = [sigma_amp/np.abs(a) for a in max_amps]
                sigma_ph = np.mean(sigma_phs)
                data = [num,freq, amp, snr]
                print('{:<15d}{:>18f}{:>25f}{:>32}'.format(data[0],
                                                           data[1],
                                                           data[2],
                                                           data[3]))
                n += 1
                num+=1
                if n > sim_fit_n:
                    all_best_freqs += fit(time, params)[2]
                    all_max_amps += fit(time, params)[1]
                    all_phs += fit(time, params)[3]
                    all_sigma_freqs += sigma_freqs
                    all_sigma_phs += sigma_phs
                    lc0 = lc
                    ls0 = periodogram(time, lc0)
                    periodograms.append(ls0)
                    n_per.append(sim_fit_n+n_per[-1])
                    n = 1                
                    params = Parameters()        
            else:
                if n!= 1:
                    all_best_freqs += fit(time, params)[2]
                    all_max_amps += fit(time, params)[1]
                    all_phs += fit(time, params)[3]
                    all_sigma_freqs += sigma_freqs
                    all_sigma_phs += sigma_phs
                    
                break
            
        elif parade == max_fap:
            if fap < parade: # Stop criterion
                snr_or_faps.append(fap)
                all_sigma_amps.append(sigma_amp)
                new_guesses = [amp, freq, ph]
                params.add('p_'+str(n)+'a', value = new_guesses[0], min=0, max=2*amp)
                params.add('p_'+str(n)+'b', value = new_guesses[1], min=freq-rayleigh/2, max=freq+rayleigh/2)
                params.add('p_'+str(n)+'c', value = new_guesses[2], min=0, max=1)
                #best_freqs = fit(time, params)[2]
                max_amps = fit(time, params)[1]
                res = minimize(residual, params, args=(time, lc0), method = 'least_squares')
                lc = res.residual
                params = res.params
                sigma_freqs = [np.sqrt(6/N)/(np.math.pi*T)*sigma_lc/np.abs(a) for a in max_amps]
                sigma_freq = np.mean(sigma_freqs)
                sigma_phs = [sigma_amp/np.abs(a) for a in max_amps]
                sigma_ph = np.mean(sigma_phs)
                data = [num,freq, amp, fap]
                print('{:<15d}{:>18f}{:>25f}{:>32}'.format(data[0],
                                                           data[1],
                                                           data[2],
                                                           data[3]))
                n += 1
                num+=1
                if n > sim_fit_n:
                    all_best_freqs += fit(time, params)[2]
                    all_max_amps += fit(time, params)[1]
                    all_phs += fit(time, params)[3]
                    all_sigma_freqs += sigma_freqs
                    all_sigma_phs += sigma_phs
                    lc0 = lc
                    ls0 = periodogram(time, lc0)
                    periodograms.append(ls0)
                    n_per.append(sim_fit_n+n_per[-1])
                    n = 1                
                    params = Parameters()        
            else:
                if n!= 1:
                    all_best_freqs += fit(time, params)[2]
                    all_max_amps += fit(time, params)[1]
                    all_phs += fit(time, params)[3]
                    all_sigma_freqs += sigma_freqs
                    all_sigma_phs += sigma_phs
                break
        
        
    l = len(all_best_freqs)        
    
    
    prew_df = pd.DataFrame({'Freqs': all_best_freqs, 'Amps': all_max_amps, 
                            'Phases':all_phs, 'Amplitude 1-sigma error (mmag)': all_sigma_amps, 
                            'Frequency 1-sigma error (c/d)': all_sigma_freqs, 
                            'Phase 1-sigma error (c/d)': all_sigma_phs, 
                            'SNR/FAP':snr_or_faps})
    
    
    res = pd.DataFrame({'Frequencies': ls[1], 'Amplitudes': ls[2]})
    res.to_csv(newpath+'res.dat', sep=' ', index=False, header = None)
    for (p, n) in zip(periodograms, n_per):
        per = pd.DataFrame({'Frequency (c/d)':p[1], 'Amplitude (mmag)':p[2]})
        per.plot(kind = 'line', x='Frequency (c/d)', y='Amplitude (mmag)', title = 'Periodogram after subtracting ' + str(n) + ' frecuencies', legend=False)
        plt.xlabel('Frequency (c/d)')
        plt.ylabel('Amplitude (mmag)')
        plt.savefig(newpath+'LS_' +str(n) + '.png')
        
    
    # Filtering the frequencies that are closer than the Rayleigh resolution

    added_freqs = []
    added_amps = []
    copied_freqs = all_best_freqs.copy()
    copied_amps = all_max_amps.copy()
    filtered_freqs = all_best_freqs.copy()
    
    for (freq, amp) in zip(all_best_freqs, all_max_amps):
        
        if freq in filtered_freqs:
            copied_freqs.remove(freq)
            copied_amps.remove(amp)
            added_freqs.append(freq)
            added_amps.append(amp)
            for (f,a)  in zip(copied_freqs, copied_amps):
                if added_freqs[-1]-rayleigh < f < added_freqs[-1]+rayleigh:
                    if added_amps[-1] > a:
                        if f in filtered_freqs:
                            filtered_freqs.remove(f)
                            prew_df = prew_df[prew_df.Freqs != f]
                    else:
                        if added_freqs[-1] in filtered_freqs:
                            filtered_freqs.remove(added_freqs[-1])
                            prew_df = prew_df[prew_df.Freqs != added_freqs[-1]]    
        else:
            copied_freqs.remove(freq)
            copied_amps.remove(amp)
        

    # Filtering the linear combinations of frequencies
    try:
        combineds = comb_freqs(prew_df)[0]
    
        for f in combineds:
            prew_df = prew_df[prew_df.Freqs != f]
    except KeyError:
        pass
    prew_df.to_csv(newpath+'best_modes.dat', sep = ' ', index = False, header = None)
    #time_df.to_csv(newpath+'velocity.dat', sep = ' ', index = False, header = None)
    end = timer()

    print(end-start)

print(dash)
print('Finished'.center(110))
print(dash)
end = timer()

print(end-start)
