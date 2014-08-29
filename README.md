TRIUMF-griffin_data_sorting
===========================

#####Author: Alison Cheeseman, 2014

C program to sort high rated data from GRIFFIN. Accepts multiple csv files containing long signal traces as input. The program can
be updated to run on different lenghs of files by changing the constants defined at the top of main.c. Waveform properties and
signal processing parameters can also be updated.

The code implements three main algorithms (triggering, baseline restoration, and the moving window deconvolution method) to
detect pulses in the input signal and evaluate their pulse height. The resulting pulse heights are histogrammed to produce a spectrum
from the input data. The energy resolution (FWHM) of a specific peak in the spectrum is calculated by fitting a Gaussian to the
peak data (FWHM = 2.35 * sigma). Depending on the input data, constants also defined in main.c can be changed to alter which peak
the program should fit.

This implementation of the the code differentiates between waveforms that have noticeable effects of pileup (either multiple
triggers in the 6us processing window or a significant negative slope to the baseline) and only evaluates the pulse height of
waveforms that are not deemed to be affected by pileup. All other waveforms are rejected.
