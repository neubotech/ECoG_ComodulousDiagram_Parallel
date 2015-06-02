#%config InlineBackend.figure_format = 'retina'
#%pylab inline
import math
import numpy as np
from scipy.signal import butter, lfilter, hilbert, filtfilt

import time 

def AmpFiltering(Signal, FilterWeights):
	transformed = np.zeros([FilterWeights.shape[0], len(Signal)]);
	# tmp = 
	for i in range(FilterWeights.shape[0]):
		# transformed[i, :] = filtfilt(FilterWeights[i, 0][0, :], [1], Signal);
		tmp = filtfilt(FilterWeights[i, 0][0, :], [1], Signal);
		transformed[i, :] = abs(hilbert(tmp));
	return transformed;

def PhaseFiltering(Signal, FilterWeights):
	transformed = np.zeros([FilterWeights.shape[0], len(Signal)]);

	for i in range(FilterWeights.shape[0]):
		# transformed[i, :] = filtfilt(FilterWeights[i, 0][0, :], [1], Signal);
		tmp = filtfilt(FilterWeights[i, 0][0, :], [1], Signal);
		transformed[i, :] = np.angle(hilbert(tmp));
	return transformed;

def ModIndex_v7(Phase, Amps, nbin):
	# t=time.clock();
	Q_Phase = np.ceil(Phase/(2*np.pi)*nbin + np.ceil(nbin/2)) - 1;

	idx = Q_Phase.argsort();
	Amps1 = Amps[:, idx];
	Q_Phase1 = Q_Phase[idx];
	idx1 = np.concatenate((np.concatenate(([0], np.diff(Q_Phase1).ravel().nonzero()[0]+1)),[ Phase.shape[0]]))	
	if (len(idx1) != nbin + 1):
		sys.exit('need longer recordings');
	else:
		MeanAmp = np.zeros([Amps.shape[0], nbin]);
		for j in range(nbin):
			MeanAmp[:, j] = np.mean(Amps1[:, idx1[j]:idx1[j+1]], axis = 1);

	p_MeanAmp = np.divide(MeanAmp.transpose(),MeanAmp.sum(axis = 1));
	MI = (np.log(nbin) - (-((p_MeanAmp)*np.log((p_MeanAmp))).sum(axis = 0)))/np.log(nbin);
	# print time.clock() - t;
	return MI, MeanAmp

def ModIndex_v6(Phase, Amps, nbin):
	# t=time.clock();
	Q_Phase = np.ceil(Phase/(2*np.pi)*nbin + np.ceil(nbin/2)) - 1;
	MeanAmp = np.zeros([Amps.shape[0], nbin]);
	count = np.zeros([1, nbin], dtype=int);

	for i in range(Phase.shape[0] - 1):
		MeanAmp[:, Q_Phase[i]] += Amps[:, i];
		count[0, Q_Phase[i]] += 1;

	MeanAmp /= count;
	p_MeanAmp = np.divide(MeanAmp.transpose(),MeanAmp.sum(axis = 1));
	MI = (np.log(nbin) - (-((p_MeanAmp)*np.log((p_MeanAmp))).sum(axis = 0)))/np.log(nbin);
	# print time.clock() - t;
	return MI, MeanAmp

def MIComodulogram(Signal, AmpFreqVec, PhaseFreqVec, AmpWs, PhaseWs):
	## Filtering with coef from matlab
	t=time.clock();
	AmpFreqTransformed = np.zeros([AmpFreqVec.shape[1], len(D)]);
	PhaseFreqTransformed = np.zeros([PhaseFreqVec.shape[1], len(D)]);

	AmpFreqTransformed = AmpFiltering(D, AmpWs);
	PhaseFreqTransformed = PhaseFiltering(D, PhaseWs);
	print time.clock() - t;

	## Comodulation loop
	t=time.clock();
	nbin = 18;

	Comodulogram = np.zeros([PhaseFreqVec.shape[1], AmpFreqVec.shape[1]], dtype = np.float32);
	for i in range(PhaseWs.shape[0] - 1):
		MI, MeanAmp = ModIndex_v7(PhaseFreqTransformed[i, :], AmpFreqTransformed, nbin);
		Comodulogram[i, :] = MI;
	print time.clock() - t;
	return Comodulogram.transpose()


if __name__ == "__main__":
	import scipy as sp
	import pandas as pd
	import matplotlib.pyplot as plt

	import os
	import sys
	import scipy.io

	## load data
	filename = 'data_ON_OFF_DBS_rest_23pts.mat';
	data = sp.io.loadmat(filename);
	print(data.keys())
	DBS_ON = data['DBS_ON'];
	D = DBS_ON[0][:];

	## set sample rate (not used currently)
	srate = 1024;

	## load filter coef from mat file
	filter = sp.io.loadmat('filter_coef_linear.mat');
	AmpFreqVec = filter['AmpFreqVector'];
	AmpWs = filter['AmpWs'];
	PhaseFreqVec = filter['PhaseFreqVector'];
	PhaseWs = filter['PhaseWs'];

	## Comodulogram
	# t=time.clock();
	Comodulogram = MIComodulogram(D, AmpFreqVec, PhaseFreqVec, AmpWs, PhaseWs);
	# print time.clock() - t;

	# plt.figure(1)
	# # im = plt.imshow(Comodulogram, aspect = 'auto', origin='lower', interpolation='spline16',extent=[np.min(PhaseFreqVec), np.max(PhaseFreqVec), np.min(AmpFreqVec), np.max(AmpFreqVec)]);
	# im = plt.imshow(Comodulogram, aspect = 'auto', origin='lower', interpolation='none',extent=[np.min(PhaseFreqVec), np.max(PhaseFreqVec), np.min(AmpFreqVec), np.max(AmpFreqVec)]);
	# # im = plt.pcolormesh(PhaseFreqVec.transpose(), AmpFreqVec.transpose(), Comodulogram)
	# plt.xlabel('Phase Frequencies (Hz)')
	# plt.ylabel('Amplitude Frequencies (Hz)')
	# plt.colorbar(im)


	# savefig('test.png')



    
