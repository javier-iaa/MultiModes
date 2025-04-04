#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Clean close frequencies from Multimodes output files best_modes.dat
Multimodes2 already implement this but this program can be used to filter files 
where the option was not activated at running time.

Filtering the frequencies that are closer than the Rayleigh resolution
It compares amplitudes and removes less significant frequencies from the 
list of filtered frequencies.

The program is called with an argument that is used as Rayleigh frequency.
Then it imports the best_modes.dat file and outputs best_modes_clean.dat

Version: 0.1
12:00 mar 27, 2025
"""

import sys
import numpy as np
import pandas as pd

# inputs: filepath, rayleigh
filepath = sys.argv[1]
R = sys.argv[2]     # rayleigh R
R = np.float64(R)
data = np.loadtxt(filepath+'best_modes.dat')
all_best_freqs = list(data[:,0])    # extract freqs
all_max_amps = list(data[:,1])  # extract amps

# The rest of frequency parameters
all_phs = data[:,2]
all_sigma_amps = data[:,3]
all_sigma_freqs = data[:,4]
all_sigma_phs = data[:,5]
snr_or_faps = data[:,6]
all_rms = data[:,7]

prew_df = pd.DataFrame(
                {'Freqs': all_best_freqs,
                'Amps': all_max_amps,
                'Phases': all_phs,
                'Amplitude 1-sigma error': all_sigma_amps,
                'Frequency 1-sigma error': all_sigma_freqs,
                'Phase 1-sigma error': all_sigma_phs,
                'SNR/FAP': snr_or_faps,
                'rms': all_rms}
                )

added_freqs = []
added_amps = []
copied_freqs = all_best_freqs.copy()
copied_amps = all_max_amps.copy()
filtered_freqs = set(all_best_freqs)  # faster membership testing

for (freq, amp) in zip(all_best_freqs, all_max_amps):
    if freq in filtered_freqs:
        # Membership checks
        if freq in copied_freqs:
            copied_freqs.remove(freq)
        if amp in copied_amps:
            copied_amps.remove(amp)

        added_freqs.append(freq)
        added_amps.append(amp)

        for (f, a) in zip(copied_freqs, copied_amps):
            if added_freqs and (added_freqs[-1] - R < f < added_freqs[-1] + R):
                if added_amps[-1] > a:
                    if f in filtered_freqs:
                        filtered_freqs.remove(f)
                        prew_df = prew_df[prew_df.Freqs != f]
                else:  # If the current amplitude is not greater, remove the current frequency
                    if freq in filtered_freqs:
                        filtered_freqs.remove(freq)
                        prew_df = prew_df[prew_df.Freqs != freq]
    else:
        # If the frequency is not in filtered_freqs, remove it from copied lists
        if freq in copied_freqs:
            copied_freqs.remove(freq)
        if amp in copied_amps:
            copied_amps.remove(amp)

# Write the frequency list
prew_df.to_csv(filepath+'best_modes_clean.dat', sep=' ', index=False, header=None)