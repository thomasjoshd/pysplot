import os
import sys
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMessageBox,QInputDialog, QLineEdit
from PyQt5.QtCore import Qt, QSize

from astropy.wcs import WCS
from astropy.io import fits
from astropy.table import Table
from astropy.convolution import convolve, Box1DKernel, Gaussian1DKernel


import astropy.units as u
from astropy.constants import c

import numpy as np
from specutils import Spectrum1D
import csv

class Spectra(QtWidgets.QMainWindow):
    def __init__(self,parent):
        super(Spectra, self).__init__(parent)

    def open(self,drop=False):
        """Open up spectrum or lists of spectra"""
        # badfile=False
        filez=False
        if drop == False:
            filez=self.openFileNamesDialog()
        else:
            filez=drop
        types=['.fit','.fits','.FIT','.FITS','.txt','.TXT','.csv','.CSV','.dat','.DAT','.s']
        if filez != False:
            lst=list(filez)
            for item in lst:
                if '.list' in item or '.lst' in item :
                    self.parent().listname=item
                    self.read_list()
                    for listitem in self.parent().listedfiles:
                        self.parent().database[listitem]={}
                else:
                    for ft in types:
                        if ft in item:
                            self.parent().database[item]={}
                        else:
                            pass
            self.parent().region_clear()
            self.loadSpectra()
            self.parent().stackrebuild()
            try:
                self.parent().plotSpectra(spec=self.parent().fname)
            except:
                pass
            self.parent().stackpane()
            self.parent().singleplottoggle()
            self.parent().reset()


    def openFileNamesDialog(self):
        files, _ = QtWidgets.QFileDialog.getOpenFileNames(self,"Open Spectra or list of spectra")
        return files

    def saveFileDialog(self):
        fileName, _ = QtWidgets.QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()")
        return fileName

    def loadSpectra(self):
        for item in self.parent().database:
            self.parent().fname=item
            try:
                self.parent().database[self.parent().fname]['wavelength']
            except:
                if '.fit' in item or '.FIT' in item:
                    self.read_fits()
                elif '.txt' in item or '.TXT' in item or '.dat' in item or '.DAT' in item:
                    self.read_txt()
                elif '.csv' in item or '.CSV' in item:
                    self.read_txt()
                elif '.s' in item:
                    self.read_espadons()
                else:
                    self.parent().message.append("File extension unexpected for spectra: fit,fits,txt,dat are currently recognized.")
                    self.parent().outputupdate()
                    self.parent().stack.remove(self.parent().fname)

    #
    # def read_fits(self):
    #     """Reads a fits file into the dictionary of stored spectra."""
    #     #need to add a way to read multispec fits files
    #     #Based on Read a UVES spectrum from the ESO pipeline
    #     try:
    #         multiextension=False
    #         try: #assume the first extnsion is the one we want
    #             tab=Table.read(self.parent().fname, hdu=1).columns
    #             # print(tab)
    #             multiextension=True
    #         except:
    #             pass
    #
    #         if multiextension == True:
    #             print("multiextension")
    #             self.read_tabledata()
    #         elif multiextension == False:
    #             # print("Assuming Spectrum to be opened is a standard 1-D fits file.")
    #             # print('1d')
    #             self.read_1Dfits()
    #         else:
    #             print("Unexpected problem opening Fits File Spectra.read_fits.")
    #             self.parent().stack.remove(self.parent().fname)
    #     except:
    #         print('problem in Spectra.read_fits')

    def read_fits(self):
        """Reads a fits file into the dictionary of stored spectra."""
        #need to add a way to read multispec fits files
        #Based on Read a UVES spectrum from the ESO pipeline
        try:
            self.read_1Dfits()

        except:
            print('problem in Spectra.read_fits')


    def read_1Dfits(self):
        print(self.parent().fname)
        sp = fits.open(self.parent().fname)
        sp.info()
        i=0
        self.setheader(sp[i].header)
        try:
            try:
                # print(self.parent().fname)
                s=Spectrum1D.read(self.parent().fname)
                # self.parent().wavelength = s.spectral_axis #specdata['loglam'] * u.AA
                # self.parent().flux = s.flux #specdata['flux'] *u.flx #* 10**-17 * u.Unit('erg cm-2 s-1 AA-1')
                self.writespec(s.spectral_axis,s.flux)
            except:
                print('Exception at Spectra.read_1Dfits, specutils not able to open spectrum.')
                try:
                    print('moving on to old method')
                    self.old()
                except:
                    print('still an issue')



            # self.parent().database[self.parent().fname]['header']=s.meta
        except:
            print('specutils unable to read spectrum')
            self.parent().poplast=True

    def old(self):
        sp = fits.open(self.parent().fname)
        sp.info()
        # layers=np.arange(len(sp))
        layers=[0]

        #try looping over each "frame" of the spectra cube.
        for i in layers:
            print(i)
            try:
                header=sp[i].header
                self.setheader(header)
                # for i,val in enumerate(sp):
                # header = sp[i].header
                if header['CUNIT1'] == '0.1 nm':
                    header['CUNIT1']='Angstrom'
                if header['NAXIS1'] > 1 or header['NAXIS'] == 1:
                    wcs = WCS(header)
                    index = np.arange(header['NAXIS1'])
                    wavelength = wcs.wcs_pix2world(index[:,np.newaxis], 0)
                    # wavelength = wcs.wcs_pix2world(wcs, sp[i].data,mode='wcs')
                    try:
                        if header['CUNIT1'] == 'Angstrom':
                            wavelength = wavelength.flatten()/1e-10
                    except:
                        pass

                    wavelength = wavelength.flatten()*u.AA
                    # print(wavelength)
                    #need to be checking wavelength uinits and type here.
                    try:
                        if 'Wavelength units=angstroms' in header['WAT1_001']:
                            self.doppler(wavelength,header)
                    except:
                        pass
                    flux = sp[i].data*u.flx
                    self.writespec(wavelength,flux)
                elif header['NAXIS'] > 1:
                    print('here')
                    wcs = WCS(header)
                    index = np.arange(header['NAXIS1'])
                    wavelength= index*u.pixel
                    flux = sp[i].data[0].flatten()*u.flx
                    self.writespec(wavelength,flux)
                    # self.parent().database[self.parent().fname]['header']=header
                else:
                    self.parent().poplast=True
            # sp.close()
            except:
                print("excepted in Spectra.old %s"%str(i))
            # # sp = fits.open(self.parent().fname)
            # i=1
            # self.setheader(sp[i].header)
            # # print('we have excepted')
            # # for i,val in enumerate(sp):
            # # header = sp[i].header
            # if header['NAXIS'] == 1:
            #     wcs = WCS(header)
            #     index = np.arange(header['NAXIS1'])
            #     wavelength = wcs.wcs_pix2world(index[:,np.newaxis], 0)
            #     # wavelength = wcs.wcs_pix2world(wcs, sp[i].data,mode='wcs')
            #     wavelength = wavelength.flatten()
            #     #need to be checking wavelength uinits and type here.
            #     if 'Wavelength units=angstroms' in header['WAT1_001']:
            #         wavelength = wavelength*u.AA
            #         self.doppler(wavelength,header)
            #     # wavelength = wavelength*u.AA
            #
            #     # self.doppler(wavelength,header)
            #
            #     flux = sp[i].data*u.flx
            #     self.writespec(wavelength,flux)
            #     # self.parent().database[self.parent().fname]['header']=header
            #     # sp.close()
            # # elif header['NAXIS'] > 1:
            # #     wcs = WCS(header)
            # #     index = np.arange(header['NAXIS1'])
            # #     wavelength= index*u.pixel
            # #     flux = sp[i].data[0].flatten()*u.flx
            # #     self.writespec(wavelength,flux)
            # #     self.parent().database[self.parent().fname]['header']=header
            #     # sp.close()
            # else:
            #     # sp.close()
            #     self.parent().poplast=True
        sp.close()
    #
    def doppler(self):
        try:
            wavelength=self.parent().wavelength/(1.0+self.parent().header['BSS_VHEL']/2.997925e05)
            self.parent().database[self.parent().fname]['heliocentric']=self.parent().header['BSS_VHEL']
            self.parent().database[self.parent().fname]['wavelength']=wavelength
            # BESS keyword is plus not minus.
        except:
            pass
        try:
            wavelength=self.parent().wavelength/(1.0-self.parent().header['VHELIO']/2.997925e05)
            self.parent().database[self.parent().fname]['heliocentric']=self.parent().header['VHELIO']
            self.parent().database[self.parent().fname]['wavelength']=wavelength
        except:
            pass

    def read_tabledata(self):
        try:
            extension=1
            sp = fits.open(self.parent().fname)
            header = sp[extension].header
            sp.close()
            tab=Table.read(self.parent().fname, hdu=extension).columns
            try:
                wavelength= tab['wavelength']
            except:
                pass
            try:
                wavelength=10**tab['loglam']
            except:
                pass
            self.writespec(wavelength*u.Angstrom,tab['flux'].quantity*u.flx)
            self.parent().database[self.parent().fname]['header']=header
        except:
            print('Problem with Read Tables')


    def read_txt(self):
        try:
            f1=open(self.parent().fname,'r')

            data=[]
            for line in f1:
                if ',' in line:
                    data.append(line.split(","))
                else:
                    data.append(line.split())
            trash=[]
            for i,row in enumerate(data):
                if len(trash)>0 and i < trash[-1]+1:
                    pass
                else:
                    try:
                        w=float(row[0])
                    except:
                        trash.append(i)

            w=[]
            fx=[]
            for i,row in enumerate(data):
                if len(trash) > 0 and i < trash[-1]+1:
                    pass
                elif len(trash) > 0 and i > trash[-1]+1:
                    w.append(float(row[0]))
                    try:
                        fx.append(float(row[1]))
                    except:
                        try:
                            fx.append(float(row[2]))
                        except:
                            pass
                if len(trash) == 0:
                    w.append(float(row[0]))
                    try:
                        fx.append(float(row[1]))
                    except:
                        try:
                            fx.append(float(row[2]))
                        except:
                            pass

            wave=np.array(w)
            f=np.array(fx)
            f=f*u.flx
            self.writespec(wave*u.AA,f)
            hdu=fits.PrimaryHDU()
            self.parent().database[self.parent().fname]['header']=hdu.header
        except:
            print('Problem opening text spectra')

    def read_espadons(self):
        try:
            f1=open(self.parent().fname,'r')
            data=np.array(list(csv.reader(f1,delimiter=' ')),dtype='object')
            wave=[]
            f=[]
            err=[]
            for i,row in enumerate(data):
                if i < 2:
                    pass
                else:
                    wave.append(float(row[2]))
                    f.append(float(row[3]))
                    err.append(float(row[5]))
            W=np.array(wave)
            idx=np.argsort(W)
            F=np.array(f)

            self.writespec(W[idx]*10*u.AA,F[idx]*u.flx)
            f1.close()
        except:
            print('problem with opening espadons data')

    def read_list(self):
        try:
            path=os.path.dirname(self.parent().listname)
            f1=open(self.parent().listname,'r')
            names=list(csv.reader(f1))
            self.parent().listedfiles=[]
            for item in names:
                if '\\' in item[0] or '/' in item[0]:
                    self.parent().listedfiles.append(item[0])
                else:
                    self.parent().listedfiles.append(os.path.abspath(os.path.join(path,item[0])))
            f1.close()
        except:
            print('problem reading list')

    def writespec(self,w,f):
        try:
            self.parent().database[self.parent().fname]['wavelength']=w
            self.parent().database[self.parent().fname]['wavelength_orig']=w
            self.parent().database[self.parent().fname]['flux']=f
            self.parent().database[self.parent().fname]['flux_orig']=f
        except:
            print('Exception occured in Spectra.writespec')
    def setheader(self,header):
        try:
            self.parent().database[self.parent().fname]['header']=header
        except:
            print('Exception in Spectra.setheader')


    def updatespectrum(self):
        self.parent().database[self.parent().fname]['wavelength']=self.parent().wavelength
        self.parent().database[self.parent().fname]['flux']=self.parent().flux
        try:
            self.parent().database[self.parent().fname]['header']=self.parent().header
        except:
            hdu = fits.PrimaryHDU(self.parent().flux.value)
            self.parent().header = hdu.header
            self.parent().header['ORIGIN'] = 'PySplot Version: %s'%VERSION

        self.parent().header['CRVAL1'] = self.parent().wavelength[0].value
        self.parent().header['CD1_1'] = self.parent().wavelength[1].value-self.parent().wavelength[0].value
        self.parent().header['CDELT1']= self.parent().wavelength[1].value-self.parent().wavelength[0].value
        self.parent().header['CUNIT1']='%s'%self.parent().wavelength.unit
        self.parent().header['CTYPE1']='WAVELENG'
        self.parent().header['CRPIX1'] = 1.
        self.parent().header['WAT0_001']='system=equispec'
        self.parent().header['WAT1_001']='wtype=linear label=waveleng Wavelength unit=angstroms'



    def saveFits(self,extend=False):
        """Save a new fits file with header."""
        try:
            self.updatespectrum()
            hdu = fits.PrimaryHDU(data=self.parent().flux.value,header=self.parent().header)

            path=os.path.dirname(self.parent().fname)
            basename=os.path.basename(self.parent().fname)
            path_wo_ext=os.path.splitext(self.parent().fname)[0]
            if extend == False:
                savename, _ = QtWidgets.QFileDialog.getSaveFileName(self,"Save fits spectra",os.path.join(path_wo_ext,basename,"%s.fits"%(path_wo_ext)))
                # hdu.writeto(savename,overwrite=True)
            else:
                savename=path_wo_ext+extend+".fits"
            hdu.writeto(savename,overwrite=True,output_verify='fix')
            self.parent().message.append("Saved to: %s"%savename)
            self.parent().outputupdate()
        except:
            self.parent().message.append("Nothing to save.")
            self.parent().outputupdate()


    def save1DText(self,extend=False,header=True):
        """Save a headerless text spectrum."""
        try:
            path=os.path.dirname(self.parent().fname)
            basename=os.path.basename(self.parent().fname)
            path_wo_ext=os.path.splitext(self.parent().fname)[0]
            if extend == False:
                savename, _ = QtWidgets.QFileDialog.getSaveFileName(self,"Save Text Spectra",os.path.join(path_wo_ext,basename,"%s.txt"%(path_wo_ext)))
            else:
                savename=path_wo_ext+extend
            # print(savename)
            dataout=open(savename,'w')
            if header==True:
                dataout.write('Wavelength %s, Flux %s\n'%(self.parent().wavelength[0].unit,self.parent().flux[0].unit))
            #add a header line

            for i,val in enumerate(self.parent().flux):
                dataout.write('%s,%s\n'%(self.parent().wavelength[i].value,self.parent().flux[i].value))
            dataout.close()
            self.parent().message.append("Saved to: %s"%savename)
            self.parent().utputupdate()
        except:
            self.parent().message.append("Nothing to save.")
            self.parent().outputupdate()


    def restore(self):
        """Resets the last/current spectrum to the original state."""
        try:
            self.parent().database[self.parent().fname]['wavelength']=self.parent().database[self.parent().fname]['wavelength_orig']
            self.parent().database[self.parent().fname]['flux']=self.parent().database[self.parent().fname]['flux_orig']
            self.parent().wavelength=self.parent().database[self.parent().fname]['wavelength_orig']
            self.parent().flux=self.parent().database[self.parent().fname]['flux_orig']
            # self.plotRegions()
            if self.parent().script == False:
                self.parent().message.append("Flux restored to the original.")
                self.parent().outputupdate()
            else:
                pass
            self.updatespectrum()
            self.parent().plotSpectra(spec=self.parent().fname) #included a canvas redraw
            self.parent().reset()
        except:
            print('Exception occred in restore')

    def smooth(self,func="boxcar"):
        """boxcar smooth"""
        try:
            self.parent().grabSpectra(self.parent().fname)
            if func == "boxcar":
                if self.parent().script == False:
                    self.parent().boxwidth, okPressed = QInputDialog.getInt(self, "Integer","Smoothing Box Width:", int(self.parent().boxwidth), 1, 100, 1)
                    if okPressed:
                        self.parent().flux=convolve(self.parent().flux,Box1DKernel(int(self.parent().boxwidth)))*u.flx
                else:
                    if self.parent().boxwidth == False:
                        self.parent().boxwidth, okPressed = QInputDialog.getInt(self, "Integer","Smoothing Box Width:", int(self.parent().boxwidth), 1, 100, 1)
                        if okPressed:
                            self.parent().flux=convolve(self.parent().flux,Box1DKernel(int(self.parent().boxwidth)))*u.flx
                    else:
                        self.parent().flux=convolve(self.parent().flux,Box1DKernel(int(self.parent().boxwidth)))*u.flx

            elif func == "gaussian":
                if self.parent().script == False:
                    self.parent().boxwidth, okPressed = QInputDialog.getDouble(self, "Float","Standard deviation in units of the wavelenth: ", float(self.parent().boxwidth), 0, 10000, 30)
                    if okPressed:
                        self.parent().flux=convolve(self.parent().flux,Gaussian1DKernel(float(self.parent().boxwidth)))*u.flx
                else:
                    if self.parent().boxwidth == False:
                        self.parent().boxwidth, okPressed = QInputDialog.getDouble(self, "Float","Standard deviation in units of the wavelenth: ", float(self.parent().boxwidth), 0, 10000, 30)
                        if okPressed:
                            self.parent().flux=convolve(self.parent().flux,Gaussian1DKernel(float(self.parent().boxwidth)))*u.flx
                    else:
                        self.parent().flux=convolve(self.parent().flux,Gaussian1DKernel(float(self.parent().boxwidth)))*u.flx



            self.parent().getlims()
            self.updatespectrum()
            self.parent().plotSpectra(spec=self.parent().fname)
        except:
            print('Exception Occured in Spectra.smooth')
