
# Pytest
# import pytest
# from xscale.filtering import linearfilters
# from xscale.signal.generator import example_xyt
# Pandas
import pandas as pd
# Xarray
import xarray as xr
# Xscale
import  xscale
# Dask
import dask.array as da
# Numpy
import numpy as np
# Scipy
import scipy.signal as sig
# Matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec

import sys
import time

def get_xarray_netcdf(filename,chunks=None):
 	return xr.open_dataset(filename, chunks=chunks)


def myplot(w):
	"""
	Plot the weights distribution of the window and the associated
	spectrum (work only for 1D and 2D windows).
	"""
	if w.ndim == 1:

		w.plot()

	elif w.ndim == 2:

		w.plot()

	elif w.ndim == 3:

		# Compute 2D spectral response
		nx = w.n[0]
		ny = w.n[1]
		spectrum = (np.fft.fft2(w.coefficients[0,:,:].squeeze(), [1024, 1024]) /
		            (np.size(w.coefficients[0,:,:].squeeze()) / 2.0))
		response = np.abs(np.fft.fftshift(spectrum / abs(spectrum).max()))
		fx = np.fft.fftshift(np.fft.fftfreq(1024, w.dx[w.dims[1]]))
		fy = np.fft.fftshift(np.fft.fftfreq(1024, w.dx[w.dims[0]]))
		# gs = gridspec.GridSpec(1, 3, width_ratios=[2, 1, 2, 1], height_ratios=[1, 2])
		gs = gridspec.GridSpec(1, 3)
		plt.figure(figsize=(11., 4.))

		# Weight disribution along x
		ax_nx = plt.subplot(gs[0])
		ax_nx.plot(np.arange(0, nx), w.coefficients[0,:,:].squeeze()[:, ny-1])
		ax_nx.set_xlim((0, nx))

		# Weight disribution along y
		ax_nx = plt.subplot(gs[1])
		ax_nx.plot(w.coefficients[0,:,:].squeeze()[nx-1, :], np.arange(0, ny))
		ax_nx.set_ylim((0, ny))

		# Frequency response for fy = 0
		ax_fx = plt.subplot(gs[2])
		# ax_fx.semilogx(fx, np.mean(response[:,:].squeeze(),axis=1), lw=1.5)
		ax_fx.semilogx(fx, response[:,512].squeeze(), lw=1.5)
		# ax_fx.set_ylim((-120, 0))
		ax_fx.set_ylabel("Normalized magnitude [dB]")
		ax_fx.set_xlabel("Frequency [cycles per sample]")
		ax_fx.grid(True)

	
		plt.tight_layout()
			
	else:
		raise ValueError("This number of dimension is not supported by the plot function")



if __name__ == "__main__":

	start_time = time.time()

	case=1

	if case==0:

		datadir = "./"
		filein = datadir+"psi2d.nc"
		cur_time = time.time()
		ds = get_xarray_netcdf(filein,chunks={'x':100, 'y':100})
		print ds
		print ds.psi
		w = ds.psi.window
		w.set(n={'x':50, 'y':50}, window={'x': 'hanning', 'y': 'hanning'})
		w.plot()	
		# plt.show()
		# sys.exit()
		res = w.convolve(compute=True)
		# myds = ds
		# myds['psi']= res
		# myds.to_netcdf(infile+'filtered')plt.figure()

		plt.pcolormesh(ds.psi)
		plt.colorbar()
		plt.figure()
		plt.pcolormesh(res)
		plt.colorbar()
		plt.show()

	elif case==1:

		datadir = "./"
		filein = datadir+"expsi.nc"
		cur_time = time.time()
		ds = get_xarray_netcdf(filein,chunks={'x':100, 'y':100, 'zt':4, 'zw':4})
		# print 'Elapsed time for read ', str(time.time() - cur_time)
		# print ds.psi
		w = ds['psi'].window
		# w.set(n={'x':50, 'y':50, 'zt':2}, window={'x': 'boxcar', 'y': 'boxcar', 'z': 'boxcar'})
		w.set(n={'x':50, 'y':50',zt':2}, window={'x': 'hanning', 'y': 'hanning', 'z':None}) #, cutoff={'x':200, 'y':200, 'z':200})
		# w.set(n={'x':50, 'y':50, 'zt':2}, window={'x': 'hanning', 'y': 'hanning', 'z': 'boxcar'}) #, cutoff={'x':200, 'y':200, 'z':200})
		# w.set(n={'x':5, 'y':5, 'zt':2}, window={'x': 'hanning', 'y': 'hanning', 'z': 'boxcar'}) #, cutoff={'x':200, 'y':200, 'z':200})
		# w.set(n={'x':50, 'y':50, 'zt':2}, window={'x': 'cosine', 'y': 'cosine', 'z': 'boxcar'})) #, cutoff={'x':200, 'y':200, 'z':200})
		print w
		myplot(w)
		res = w.convolve(compute=True)
		# myds = ds
		# myds['psi']= res
		# myds.to_netcdf(infile+'filtered')

		print 'Elapsed time for all ', str(time.time() - start_time)

		plt.figure()
		plt.pcolormesh(ds.psi[:,:,-1])
		plt.colorbar()
		plt.figure()
		plt.pcolormesh(res[:,:,-1])
		plt.colorbar()
		plt.show()


	else:

		datadir = "./"
		filein = datadir+"nemo_psi.nc"
		cur_time = time.time()
		ds = get_xarray_netcdf(filein,chunks={'x':1032, 'y':756, 'zt':4, 'zw':4})
		print 'Elapsed time for read ', str(time.time() - cur_time)
		print ds.psi
		# w = ds.psi.window
		w = ds['psi'].window
		# w.set(n={'x':129, 'y':126, 'zt':2}, window={'x': 'hanning', 'y': 'hanning'}) # hyper long
		w.set(n={'x':43, 'y':42, 'zt':2}, window={'x': 'hanning', 'y': 'hanning'})
		# w.set(n={'x':129, 'y':126}, window='hanning', cutoff={'x':20, 'y':20})
		# w.set(n={'x':129, 'y':126, 'zt':1},dim=['x', 'y'], window='hanning', cutoff=None, chunks={'x': 1032, 'y': 756, 'zt':1})
		print w
		print dir(w)
		# w.plot()
		# w.set(n={'x':25, 'y':25, 'zt':25}, window='hanning')
		res = w.convolve(compute=True)
		# myds = ds
		# myds['psi']= res
		# myds.to_netcdf(infile+'filtered')
		print ds

		print 'Elapsed time for all ', str(time.time() - start_time)

		plt.figure()
		plt.pcolormesh(ds.psi[:,:,-1])
		plt.colorbar()
		plt.figure()
		plt.pcolormesh(res[:,:,-1])
		plt.colorbar()
		plt.show()
