## What is it?
MultiModes2 is a Python code to extract the most significant frequencies of a set of lightcurve (or other time series) files contained in a folder.
The original Multimodes code is authored by David Pamos Ortega (University of Granada).  

Modifications by Javier Pascual Granado (IAA-CSIC) are aimed to make the code more efficient, accesible and extensible.

Other contributors: Antonio García Hernández (UGR), Sebastiano de Franciscis (IAA-CSIC) and Cristian Rodrigo.

## Description
MultiModes2 takes as input a directory with light curves, corrected from 'outliers' and 'nan' values, in fits or ascii format and the initial parameters written in a text file named ini.txt

With every light curve, the code calculates the frequencies spectrum, or periodogram, with the Fast Lomb Scargle algorithm (Press & Ribicky 1989). It extracts the higher amplitude peak and evaluates if it is real signal or due to noise, either by the False Alarm Probability or by the Signal to Noise criterion, it is a decision of the user at the time of choosing the initial parameters. By default it is chosen to use as  stop criterion that S/N is greater than 4, (Breger 1993). Then, Multimodes fits frequency, amplitude and phase through non-linear optimization, using a multisine function. This function is redefined with the new calculated parameters at each iteration. It does a simultaneous fit of a number of peaks (20 by default). Then, they are subtracted from the original signal and goes back to the beginning of the loop with the residual, repeating the same process, until the stop criterion is reached. 
 
Multimodes make use of [astropy](https://www.astropy.org) for the calculation of the periodograms and [lmfit](https://lmfit.github.io/lmfit-py/) for the non-linear and simultaneous fitting of the extracted signals, using the non-linear least squares method for python.

## Citations
If MultiModes2 is used, please, cite the author this way: Pamos Ortega, D. et al. 2022 (https://doi.org/10.1093/mnras/stac864)

## Requirements
Python >=3.8 with the following modules installed:
- numpy >=1.19.2
- matplotlib >=3.3.2
- pandas >=1.1.2
- astropy >=4.0.2
- lmfit >=1.0.2
- scipy >=1.5.2

## Input
- Directory with light curves in fits or ASCII format (.dat or .txt), corrected from 'outliers' and 'nan' values. Note that the first row (header) in ASCII files is dropped by default.
- ini.txt with optional parameters

## Optional parameters:
- sim_fit_n: Number of simultaneous peaks to be fitted before extracting to the original light curve for obtaining the residual: 20 by default
- max_freq: Maximum value of the analysed frequencies domain: 100 c/d by default (delta Scuti stars)
- os_ratio: oversampling factor, 5 by default
- stop: Stop criterion, FAP or SNR, SNR by default
- min_snr: Minimum signal to noise ratio, 4 by default (Breger 1993)
- max_fap: Maximum value of the False Alarm Probability, 0.01 by default (Balona et al. 2014)
- timecol: column for time 
- fluxcol: column for fluxes
- save_data_res: save data of residual every xx iterations
- save_freq_list: parameter allows to save intermediate results files each xx iterations. It is recommended that this number coincides with save_data_res being a multiple of sim_fit_n also.
- save_plot_resps: write residual spectrum after last iteration (binary flag)
- max_iter: maximum number of iterations
- header_lines: skip xx header lines
  
## Output
- Directory 'results', containing subdirectories corresponding to every analysed light curve. Each subdirectory contains:
  - file best_modes.dat with 8 columns containing the values of the most significant frequencies, amplitudes, phases, the corresponding errors, SNR/FAP and the rms
  - res_lc.dat (optional), with the final residual after extracting all the most significant frequencies
  - Other intermediate residual files if save_data_res and/or save_freq_list parameter are activated
  - res_ps.dat (optional), with the periodogram of the final residual after extracting all the most significant frequencies.
  
Screen output shows the parameters for the peak at maximum amplitude of the Lomb-Scargle periodogram for each iteration.

# How to run it
- Copy MultiModes2.py and ini.txt to your working directory.
- Put also the directory or file with the light curves to be analysed inside your working directory.
- Enter your working directory and type the command: 
`python MultiModes2.py --d <lightcurves_directory>`
or alternatively
`python MultiModes2.py --file <lightcurve_file>`

If you want to use other parameter file than ini.txt you can call the program using the flag --p <ini_file>. For example, 
`python MultiModes2.py --file lightcurve.dat --p ini2.txt`
In this example a suffix "ini2" will be added to the output folder such as `\results\lightcurve_ini2` thus allowing to run the program with different configurations.

Since v0.1.1 MultiModes2 can be also imported from a Jupyter Notebook. For example:

```
>>> import MultiModes2 as mm2
>>> args = mm2.arguments(['--file', 'example.txt'])
>>> mm2.multimodes(args)
```

The output, though, will not appear in the notebook yet.

