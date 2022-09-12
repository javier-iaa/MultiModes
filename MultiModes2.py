#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 20:25:19 2020

MultiModes is a code to extract the most significant frequencies 
from a sample of variable stars. For each of them, it calculates 
the LombScargle periodogram and performs a nonlinear 
simultaneous fit, using a multisine function, 
to a set number of frequencies. As a stop criterion, 
you can choose between the FAP or the SNR criterion, 
depending on the type of analysis you want to perform. 

@author: David Pamos Ortega (UGR)
Modified by Javier Pascual Granado (IAA-CSIC) and Cristian Rodrigo on Jul 2022.
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

# Default initial parameters

osratio = 5 # Oversampling ratio 
max_freq = 100 # Max frequency in the periodogram
sim_fit_n = 20 # Max number of frequencies of the simultaneous fit
stop = 'SNR' # Stop criterion
min_snr = 4.0 # Min value of Signal-to-Noise Ratio (SNR) for detection of a frequency
max_fap = 0.01 # Max value of False Alarm Probability (FAP)
#tail_per = 80 # Min frequency of the tail of the periodogram used for noise estimation
timecol = 1 # column order for time and fluxes 
fluxcol = 2
save_plot_per = 0 # save plots of periodogram every xx iterations
#save_plot_resps = 0 # save plot of the power spectrum of the residual after last iteration

# Initializing the necessary lists

all_best_freqs = []  # extracted frequencies
all_max_amps = []    # extracted amplitudes
all_phs = []         # extracted phases
all_sigma_amps = []  # min errors in the amplitudes
all_sigma_freqs = [] # min errors in the frequencies
all_sigma_phs = []   # min errors in the phases
params = Parameters()
all_faps = []  # FAP values for each extracted frequency
S_N = []       # SNR values for each extracted frequency
all_rms = []   # RMS of the residuals for each step of the pre-whitening cascade
    
# Function definitions

def sinusoid(t, A, f, ph): 
    ''' Sine function to fit a frequency of a light curve'''
    return A*np.sin(2*np.math.pi*(f*t + ph))

def noise_estimate(time, flux, n=20): 
    '''Estimation of noise using the mean of near-Nyquist frequency bins'''
    ls = LombScargle(time, flux, normalization = 'psd', center_data = True, fit_mean = False)
    frequency, power = ls.autopower()
    l = len(time)
    amps = 2*np.sqrt(power/l)
    if l<=n:
        n=l-1

    return np.mean( amps[(n-1):-1] ) # last n bins excluding nyquist

def periodogram(time, flux):
    '''Fast Lomb Scargle to calculate the periodogram (Press & Ribicky 1989)'''
    ls = LombScargle(time, flux, normalization = 'psd', center_data = True, fit_mean = False)
    frequency, power = ls.autopower(method = 'fast', maximum_frequency = max_freq, samples_per_peak = osratio)
    indmpow = np.argmax(power)
    best_frequency = frequency[indmpow]
    amps = 2*np.sqrt(power/len(time))
    ampmax = amps[indmpow]
    #noise = []
    #for (f, a) in zip(frequency, amps):
    #    if f > tail_per:
    #        noise += [a]
    #noise = np.mean(noise)
    return ls, frequency, amps, best_frequency, ampmax, power#, noise

def fit(t, params):
    '''Multi-sine fit function with all the parameters of frequencies, amplitudes, phases'''
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
    return flux + fit(t, params)[0]
#    return fit(t, params)[0] - flux

def lightcurve(file):
    '''Reading the file to extract all data'''
    if isfits:
        hdul = fits.open(file)
        data = hdul[1].data
        data = pd.DataFrame(np.array(data))
        time = np.array( data.iloc[:,timecol-1] )     # extract times
        fluxes = np.array( data.iloc[:,fluxcol-1] )   # extract fluxes
    elif isascii:
        data = np.loadtxt(file, skiprows=1)
        time = data[:,timecol-1]     # extract times
        fluxes = data[:,fluxcol-1]  # extract fluxes

    time = time - time[0]
    T = time[-1] - time[0]
    N = len(time)
    r = 1/T   # Rayleigh frequency resolution

#   mean_flux = np.mean(fluxes)
#   fluxes = (fluxes-mean_flux)/mean_flux*1000 # convert fluxes to mmag

    return time, fluxes, T, N, r

def snr_or_fap(par):
    '''Choosing between the SNR or the FAP stop criterion. 
       SNR by default'''
    if 'SNR' in par:
        return min_snr, S_N
    elif 'FAP' in par:
        return max_fap,  all_faps

'''
def comb_freqs(pd):
    # Function to detect the combination frequencies once all of them having been extracted
    
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
'''

# Here begins the main part of the code

dash = 100*'-'
print(dash)
print('Running MultiModes'.center(110))
print(dash)
    
# Creating the directory with the fits files as argument to the command line    
parser = argparse.ArgumentParser(description='Open the directory with the files to be processed')
parser.add_argument('--p', type = str, help = 'Select a parameters file')
command_group = parser.add_mutually_exclusive_group()
command_group.add_argument('--d', type = str, help ='Select a directory containing a list of files')
command_group.add_argument('--file', type = str, help ='Select a single file to be opened')
args = parser.parse_args()

# Reading the initial file with the values of the parameters, if it exists 
if args.p:
    paramname = args.p
else:
    paramname = 'ini.txt'

if os.path.isfile(paramname):
    with open(paramname, 'r') as file:
        for line in file:
            line.replace("\n", "")
            if line.startswith("sim_fit_n"):
                sim_fit_n = int(line.split(' ')[1])
            if line.startswith("osratio"):
                osratio = float(line.split(' ')[1])
            if line.startswith("maxfreq"):
                max_freq = float(line.split(' ')[1])
            if line.startswith("max_fap"):
                max_fap = float(line.split(' ')[1])  
            #if line.startswith("tail_per"):
            #    tail_per = float(line.split(' ')[1])
            if line.startswith("min_snr"):
                min_snr = float(line.split(' ')[1])
            if line.startswith("stop"):
                stop = str(line.split(' ')[1])
            if line.startswith("timecol"):
                timecol = int(line.split(' ')[1])
            if line.startswith("fluxcol"):
                fluxcol = int(line.split(' ')[1])
            if line.startswith("save_plot_per"):
            	save_plot_per = int(line.split(' ')[1])
            if line.startswith("save_plot_resps"):
            	save_plot_resps = int(line.split(' ')[1])
else:
    print('Not ini.txt. Default values will be used:')     


print('Number of frequencies of the simultaneous fit: ' + str(sim_fit_n))
print('Samples per peak: ' + str(osratio))
print('Maximum frequency: ' + str(max_freq))
if 'SNR' in stop:
    print('Stop Criterion: SNR > ' + str(min_snr))
elif 'FAP' in stop:
    print('Stop Criterion: FAP < ' + str(max_fap))


# Creating the list with all the files and filenames
isfits = False
isascii = False
if args.d:
    pth = './'+args.d+'/'
    fits_files = [f for f in glob.glob(pth+'*.fits')]
    ascii_files = [f for f in glob.glob(pth+'*.dat')]
    if len(ascii_files)==0:
        if len(fits_files)==0:
            print("No input files found.")

        fname = [os.path.splitext(f)[0] for f in fits_files]
        fits_files = sorted(fits_files)
        filepath = fits_files
        isfits = True
    else:
        fname = [os.path.splitext(f)[0] for f in ascii_files]
        ascii_files = sorted(ascii_files)
        filepath = ascii_files
        isascii = True
    
    fname = sorted(fname)

else:
    fname = list( os.path.splitext(args.file) )
    if fname[1]=='.fits':
        fits_files = [ './'+args.file ]
        filepath = fits_files
        isfits = True
    elif fname[1]=='.dat' or fname[1]=='.txt':
        ascii_files = [ './'+args.file ]
        filepath = ascii_files
        isascii = True

    fname = [ fname[0] ]

# Starting with the iterative analysis
columns = ['Number',
           'f',
           'A',
           stop]

# Setting stop criterion
parade = snr_or_fap(stop)[0]
snr_or_faps = snr_or_fap(stop)[1]

for (f, nm) in zip(filepath, fname):
    start = timer() # Counting executing time for every analysed light curve
    n = 1
    num = 1
    if os.path.isfile(paramname):
        newpath = './results/' + nm + '_' + paramname[:-4] + '/'
    else:
        newpath = './results/' + nm + '/'
        
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    data = lightcurve(f) # Extracting data of every light curve
    time = data[0]     # Time vector
    lc = data[1]       # Light curve
    T = data[2]        # Total time span
    N = data[3]        # Number of points
    rayleigh = data[4] # Rayleigh resolution

    sigma_lc = stats.sem(list(lc)) # 1-sigma error of fluxes
    sigma_amp = np.sqrt(2/N)*sigma_lc # 1-sigma error of the amplitude
    
    # Save lightcurve as ASCII file (for FITS) and a scatter plot
    """
    lc_df = pd.DataFrame({'Time':time, 'Flux':list(lc)})
    lc_df.to_csv(newpath+'lc.dat', sep = ' ', index=False, header = None)
    lc_df.plot(kind='scatter', x='Time', y = 'Flux', color='blue', s = 1, title=nm)
    plt.savefig(newpath+'LC.png')
    plt.close()
    """
    
    # Calculating the initial periodogram
    lc0 = lc
    ls0 = periodogram(time, lc0)
    periodograms = [ls0,]
    n_per = [0,]
    per = pd.DataFrame({'Frequency':ls0[1], 'Amplitude':ls0[2]})
    
    # Noise level estimation
    no = noise_estimate(time, lc)

    # Save the initial periodogram as ASCII file and line plot
    """
    per.to_csv(newpath+'pg.dat', sep=' ', index=False, header = None)
    per.plot(kind = 'line', x='Frequency', y='Amplitude', title = nm, legend = False)
    plt.xlabel('Frequency')
    plt.ylabel('Amplitude')
    plt.savefig(newpath+'LS.png')
    plt.close()
    """

    # Initialization of the lists to save the extracted frequencies, amplitudes and phases
    all_best_freqs = []
    all_max_amps = []
    all_phs = []
    snr_or_faps = []
    all_sigma_amps = []
    all_sigma_freqs = []
    all_sigma_phs = []
    all_rms = []
    params = Parameters()
    
    # Print header table
    print(dash)
    print('{:<15s}{:>18s}{:>25s}{:>32}' .format(columns[0],
                                                columns[1],
                                                columns[2],
                                                columns[3]))
    print(dash)
    
    while n <= sim_fit_n:
        ls = periodogram(time, lc)
        freq = ls[3]  # freq at max power
        amp = ls[4]   # amp at max power
        rms = np.sqrt( sum(lc**2)/N )
        sigma_freq = np.sqrt(6/N)/(np.math.pi*T)*sigma_lc/amp  # freq error
        ph = 0.5
        snr = amp/no
        var = (N/4)*no**2
        fap = 1 - (1 - np.exp(-ls[5].max()/var))**(N/2)
        
        if parade == min_snr:
            if snr > parade: # Stop criterion
                snr_or_faps.append(snr)
                all_rms.append(rms)
                all_sigma_amps.append(sigma_amp)
                new_guesses = [amp, freq, ph]
                params.add('p_'+str(n)+'a', value = new_guesses[0]) # changed: params.add('p_'+str(n)+'a', value = new_guesses[0], min=0, max=2*amp)
                params.add('p_'+str(n)+'b', value = new_guesses[1], min=freq-rayleigh/2, max=freq+rayleigh/2)
                params.add('p_'+str(n)+'c', value = new_guesses[2]) # changed: params.add('p_'+str(n)+'c', value = new_guesses[2], min=0, max=1)
                #best_freqs = fit(time, params)[2]
                max_amps = fit(time, params)[1]
                res = minimize(residual, params, args=(time, lc0), method = 'leastsq') # method changed from 'least_squares' to 'leastsq'
                lc = res.residual
                params = res.params
                
                # Error estimation
                sigma_freqs = [np.sqrt(6/N)/(np.math.pi*T)*sigma_lc/np.abs(a) for a in max_amps]
                sigma_freq = np.mean(sigma_freqs)
                sigma_phs = [sigma_amp/np.abs(a) for a in max_amps]
                sigma_ph = np.mean(sigma_phs)
                
                # Print table of extracted frequencies from the LS periodogram
                data = [num, freq, amp, snr]
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
                all_rms.append(rms)
                all_sigma_amps.append(sigma_amp)
                new_guesses = [amp, freq, ph]
                params.add('p_'+str(n)+'a', value = new_guesses[0])
                params.add('p_'+str(n)+'b', value = new_guesses[1])
                params.add('p_'+str(n)+'c', value = new_guesses[2])
               # best_freqs = fit(time, params)[2] # comentar
                max_amps = fit(time, params)[1]
                res = minimize(residual, params, args=(time, lc0), method = 'leastsq') # method changed from 'least_squares' to 'leastsq'
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
                            'SNR/FAP':snr_or_faps,
                            'rms':all_rms}
                           )
    
    # residual spectrum after last iteration
    """
    if save_plot_resps !=0:
    	resps = pd.DataFrame({'Frequencies': ls[1], 'Amplitudes': ls[2]})
    	resps.to_csv(newpath+'res_ps.dat', sep=' ', index=False, header = None)
    """
        
    # residual lightcurve after last iteration
    reslc = pd.DataFrame({'Time': time, 'Residuals': lc})
    reslc.to_csv(newpath+'res_lc.dat', sep=' ', index=False, header = None)

    # save plots
    if save_plot_per !=0:
        for (p, n) in zip(periodograms, n_per):
    	    if np.mod(n, save_plot_per)==0:
                per = pd.DataFrame({'Frequency':p[1], 'Amplitude':p[2]})
                per.plot(kind = 'line', x='Frequency', y='Amplitude', title = 'Periodogram after subtracting ' + str(n) + ' frequencies', legend=False)
                plt.xlabel('Frequency')
                plt.ylabel('Amplitude')
                plt.savefig(newpath+'LS_' +str(n) + '.png')
                plt.close()
        
    
    # Filtering the frequencies that are closer than the Rayleigh resolution
    """
    added_freqs = []
    added_amps = []
    copied_freqs = all_best_freqs.copy()
    copied_amps = all_max_amps.copy()
    filtered_freqs = all_best_freqs.copy()
    
    for (freq, amp) in zip(all_best_freqs, all_max_amps):
        
        if freq in filtered_freqs:
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
    """    
    prew_df.to_csv(newpath+'best_modes.dat', sep = ' ', index = False, header = None)
    #time_df.to_csv(newpath+'velocity.dat', sep = ' ', index = False, header = None)
    end = timer()

    print('Executing Time: ' + str(end-start))

print(dash)
print('Finished'.center(110))
print(dash)
