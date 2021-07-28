from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt #basic plotting
from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel
import csv

def readfits(fname):
    sp = fits.open(fname)
    header = sp[0].header
    wcs = WCS(header)
    index = np.arange(header['NAXIS1'])
    wavelength = wcs.wcs_pix2world(index[:,np.newaxis], 0)
    wavelength = wavelength.flatten()
    wavelength = wavelength*u.AA
    flux = sp[0].data*u.flx
    sp.close()
    return wavelength,flux

def telluric(fwhm,vel,wobs,scale):
    data=np.array(list(csv.reader(open('transdata_0.5_1_mic.txt','r'),delimiter=' ')))
    wavenum=[]
    fa=[]
    for row in data:
        wavenum.append(float(row[0]))
        if row[1] == '':
            fa.append(float(row[2]))
        else:
            fa.append(float(row[1]))
    w=[]
    fl=[]
    for i,val in enumerate(fa):
        if val > 0:
            fl.append(val)
            w.append(wavenum[i])

    fa=np.array(fl)
    f=fa[::-1]

    vacw=1e8/np.array(w)
    airwav=vacw/(1 + 2.735182e-4+ (11.46377774e0/vacw)**2+ (128.9214492e0/vacw)**4) / 1e4
    wavea=airwav*10000
    waveang=wavea[::-1]
    f1=np.interp(wobs.value,waveang,f)


    gauss_kernel = Gaussian1DKernel(fwhm)  #here the 2 is stdev
    smoothed_data_gauss = convolve(f1, gauss_kernel)

    return smoothed_data_gauss*scale+1-scale


w,f=readfits('../2011-cropn.fits')
w=w-.28*u.AA
plt.plot(w,f,label='data')

ftell=telluric(2,0,w,0.7)
plt.plot(w,ftell)
plt.plot(w,f/(ftell))
plt.show()
