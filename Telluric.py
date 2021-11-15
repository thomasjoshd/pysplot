from PyQt5 import QtWidgets

import numpy as np

from astropy.wcs import WCS
from astropy.io import fits
import astropy.units as u
from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel

import csv
# import matplotlib.pyplot as plt #basic plotting


class Telluric(QtWidgets.QMainWindow):
    def __init__(self,parent):
        super(Telluric, self).__init__(parent)


    def telluric(self,fwhm,vel,scale):
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
        f1=np.interp(self.parent().wavelength.value,waveang,f)


        gauss_kernel = Gaussian1DKernel(fwhm)
        smoothed_data_gauss = convolve(f1, gauss_kernel)

        return smoothed_data_gauss*scale+1-scale

#what it should do is create the telluric profile as a new spectrum in the database
#Telluric-num as one generate multiple.
#The vel,scale, and fwhm need to be adjustable parameters, perhaps requireing its
#own subwindow like headerwin, or arithwin.
