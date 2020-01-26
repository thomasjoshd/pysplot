"""Author: Dr. Joshua Thomas
thomas.joshd@gmail.com
This program was designed to emulate some basic IRAF splot functions.
Tested on Python 3.6.5 and 3.6.7, Linux Mint 19.1 and Windows 10.
Uses Astropy library, and some parts are directly modified from the UVES tutorial.
Tested using astropy-3.2.3 numpy-1.17.4"""

version="0.6.3"

import sys
if sys.version_info < (3, 5):
    raise "must use python 3.6 or greater"
else:
    pass
import platform
OperatingSys=platform.system() #File dropdown not supported on Mac

from functools import partial
import numpy as np #arrays and math
import csv
import os
import tkinter as tk
#clean up the messy calls!!
import matplotlib.pyplot as plt #basic plotting

try:
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
except:
    pass
try:
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
except:
    pass

from astropy.wcs import WCS
from astropy.io import fits

import astropy.units as u
from astropy.constants import c

from astropy.modeling import models, fitting
from astropy.convolution import convolve, Box1DKernel

from scipy.interpolate import interp1d

import datetime
import webbrowser

UPDATED='{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())



plot_params = {'axes.linewidth': 1,
               'xtick.labelsize': 'medium',
               'ytick.labelsize': 'medium',
               'xtick.major.pad': 8,
               'xtick.major.size': 12,
               'xtick.minor.size': 6,
               'ytick.major.size': 12,
               'ytick.minor.size': 6,
               'xtick.direction': 'in',
               'ytick.direction': 'in',
               'patch.linewidth': 0,
               'font.size': 12}#,
               #'axes.grid':True}
plt.rcParams.update(plot_params)


def find_nearest_index(A,value):
    return (abs(A-value)).argmin()

def fill_nan(A):
    '''
    interpolate to fill nan values
    '''
    inds = np.arange(A.shape[0])
    good = np.where(np.isfinite(A))
    if len(good[0]) == 0:
        return np.nan_to_num(A)
    f = interp1d(inds[good], A[good], bounds_error=False)
    B = np.where(np.isfinite(A), A, f(inds))
    return B


class App:
    def __init__(self,master):
        self.master=master
        master.title("PySplot, Version %s"%str(version))

        self.promptframe=tk.Frame()
        self.promptframe.pack(side="top")
        l1=tk.Label(self.promptframe, text="Prompt:").pack( side = "left")
        self.output=tk.Entry(self.promptframe,width=100)
        self.output.pack(side = "left")

        self.inputframe=tk.Frame()
        self.inputframe.pack(side="top")
        l2=tk.Label(self.inputframe, text="Input:  ").pack( side = "left")
        self.w1v=tk.StringVar()
        self.w1 = tk.Entry(self.inputframe,text=self.w1v,width=100)
        self.w1.pack(side = "left")
        self.w1v.set(1)	#default value of so the program doesn't barf if accidentally used.

        #--------------------------------------------------
        #file menu
        menu = tk.Menu(self.master)
        master.config(menu=menu)

        filemenu = tk.Menu(menu)
        menu.add_cascade(label="File", menu=filemenu)
        filemenu.add_command(label="Open 1-D Fits (o)",command=self.openSpectra)
        filemenu.add_command(label="Restore Original",command=self.restore)
        filemenu.add_command(label="Image Header (h)",command=self.imhead)
        filemenu.add_separator()
        filemenu.add_command(label="Save New Fits",command=self.save_fits)
        filemenu.add_separator()
        filemenu.add_command(label="Save Headerless Text Spectra",command=self.save1DText)
        filemenu.add_separator()
        filemenu.add_command(label="Close All Spectra",command=self.stackreset)
        filemenu.add_command(label="Exit", command=self._quit)

        viewmenu = tk.Menu(menu)
        menu.add_cascade(label="View", menu=viewmenu)
        viewmenu.add_command(label="Reset View (r)",command=self.reset)
        viewmenu.add_command(label="Grid Toggle (|)",command=self.gridtoggle)
        viewmenu.add_command(label="Over Plot ([)",command=self.overplottoggle)
        viewmenu.add_command(label="Stack Plot (])",command=self.stackplottoggle)
        viewmenu.add_command(label="Single Plot (\\)",command=self.singleplottoggle)
        viewmenu.add_command(label="Dynamical--WIP (|)",command=self.dynamical)


        modmenu = tk.Menu(menu)
        menu.add_cascade(label="Modify", menu=modmenu)
        modmenu.add_command(label="Set Continuum (s)", command=self.continuum)
        modmenu.add_command(label="Normalize (t)", command=self.normalize)
        modmenu.add_command(label="Reset Normalization Parameters (q)", command=self.norm_clear)
        modmenu.add_command(label="Save Norm Parameters",command=self.SaveNorm)
        modmenu.add_command(label="Load Norm Parameters",command=self.LoadNorm)
        modmenu.add_separator()
        modmenu.add_command(label="Boxcar Smooth (b)", command=self.smooth)
        modmenu.add_command(label="Crop Spectra (c)", command=self.scopy)
        modmenu.add_separator()
        modmenu.add_command(label="Convert Wavelength <-> Velocity (u)", command=self.velocity)

        regionmenu = tk.Menu(menu)
        menu.add_cascade(label="Region", menu=regionmenu)
        regionmenu.add_command(label="Define Region (x)", command=self.regionload)
        regionmenu.add_command(label="Clear Region (r)", command=self.reset)
        regionmenu.add_command(label="Equivalent Width (e)", command=self.eqw)
        regionmenu.add_command(label="Gaussian (g)", command=partial(self.fit,func="gauss"))
        regionmenu.add_command(label="Voigt (v)", command=partial(self.fit,func="voigt"))
        regionmenu.add_command(label="Lorentzian (l)", command=partial(self.fit,func="lorentz"))
        regionmenu.add_command(label="Crop Spectra (c)", command=self.scopy)
        regionmenu.add_command(label="Save EQW/Fit Region", command=self.SaveRegion)
        regionmenu.add_command(label="Load EQW/Fit Region", command=self.LoadRegion)

        regionmenu.add_separator()
        regionmenu.add_command(label="Bisect Feature (w)",command=self.BisectLine)
        regionmenu.add_command(label="Save Bisection Regions", command=self.SaveBisect)
        regionmenu.add_command(label="Load Bisection Regions", command=self.LoadBisect)

        stackmenu = tk.Menu(menu)
        menu.add_cascade(label="Stack", menu=stackmenu)
        stackmenu.add_command(label="Stack Plot Mode (])",command=self.stackplottoggle)
        stackmenu.add_command(label="Show Stack Pane ({)",command=self.stackpane)
        stackmenu.add_command(label="Hide Stack Pane (})",command=self.hidepane)
        stackmenu.add_command(label="Remove Selected From Stack",command=self.removefromstack)
        stackmenu.add_command(label="Detailed Database Entry",command=self.currentspectra)
        stackmenu.add_command(label="Print Stack List",command=self.stackprint)
        stackmenu.add_command(label="Save Stack List",command=self.savestack)
        stackmenu.add_command(label="Stack Clear",command=self.stackreset)
        stackmenu.add_separator()
        # stackmenu.add_command(label="Load Norm Parameters",command=self.LoadNorm)
        stackmenu.add_command(label="Normalize",command=partial(self.stacker,func="norm"))
        stackmenu.add_separator()
        # stackmenu.add_command(label="Load Region",command=self.LoadRegion)
        stackmenu.add_command(label="Equivalent Width",command=partial(self.stacker,func="eqw"))
        stackmenu.add_command(label="Gaussian",command=partial(self.stacker,func="gauss"))
        stackmenu.add_command(label="Voigt",command=partial(self.stacker,func="voigt"))
        stackmenu.add_command(label="Lorentzian",command=partial(self.stacker,func="lorentz"))
        stackmenu.add_command(label="Crop",command=partial(self.stacker,func="scopy"))
        stackmenu.add_separator()
        # stackmenu.add_command(label="Load Bisection Regions", command=self.LoadBisect)
        stackmenu.add_command(label="Bisect",command=partial(self.stacker,func="bisect"))


        helpmenu = tk.Menu(menu)
        menu.add_cascade(label="Help", menu=helpmenu)
        helpmenu.add_command(label="About", command=self.About)

        #--------------------------------------------------
        #some initial parameters
        self.output.insert(tk.END,"Get started by opening a 1-D Spectrum or List of 1-D Spectra. File>open (or press o)")
        self.generate_plot()
        self.gridvalue=True
        self.overplot=False
        self.stackplot=False
        self.stack=[]
        self.stackint=0
        self.listedfiles=[] #used for lists files
        self.database={}
        self.loadedregions=False
        self.loadednorm=False
        self.loadedbisect=False
        self.pane=False



        #keyboard shortcuts (listed alphabetically)
        self.master.bind('b',self.smooth)
        self.master.bind('c',self.scopy)

        self.master.bind('e', self.eqw)

        self.master.bind('g', partial(self.fit,"gauss"))

        self.master.bind('h', self.imhead)



        self.master.bind('l', partial(self.fit,"lorentz"))


        self.master.bind('o', self.openSpectra)

        self.master.bind('q', self.norm_clear)
        self.master.bind('r', self.reset)
        self.master.bind('s', self.continuum)
        self.master.bind('t', self.normalize)
        self.master.bind('u', self.velocity)
        self.master.bind('v', partial(self.fit,"voigt"))
        self.master.bind('w', self.BisectLine)
        self.master.bind('x', self.regionload)
        self.master.bind('X', self.region_clear)



        self.master.bind('|', self.gridtoggle)
        self.master.bind('[', self.overplottoggle)
        self.master.bind(']', self.stackplottoggle)
        self.master.bind('\\', self.singleplottoggle)
        self.master.bind('|', self.dynamical)
        self.master.bind('{', self.stackpane)
        self.master.bind('}', self.hidepane)
        self.master.bind('<space>', self.coord)

        self.master.bind("<MouseWheel>",self.MouseWheelHandler)
        self.master.bind("<Button-4>",self.MouseWheelHandler)
        self.master.bind("<Button-5>",self.MouseWheelHandler)

    def MouseWheelHandler(self,event=None):
        print("Nyan")

    def captainslog(self):
        """Save the output from measuremetns to a csv"""
        path=os.path.dirname(self.fname)
        time='{0:%Y%m%d.%H%M%S}'.format(datetime.datetime.now())
        basename=os.path.basename("pysplot-%s.log"%time)
        savename=os.path.join(path,basename)
        print("Log being written to: %s"%savename)
        self.starlog=open(savename,'w')
        self.starlog.write('PySplot Log %s'%time)
        self.starlog.write('\n')

    def endoflog(self):
        try:
            if self.starlog.closed == False:
                self.starlog.close()
                del self.starlog
                print('del log mem')
        except:
            pass

    def checklog(self):
        try:
            self.starlog
        except:
            self.captainslog()

    def openSpectra(self,event=None):
        """Open up spectrum or lists of spectra"""
        if self.overplot == False and self.stackplot == False:
            if 'Darwin' not in OperatingSys:
                filez=tk.filedialog.askopenfilenames(title='Choose a single spectrum',filetypes=(("Fits Files", "*.fit*"),
                                                            ("Fits Files", "*.FIT* "),
                                                            ("Spectra List", "*.list"),
                                                            ("Spectra List", "*.lst"),
                                                            ("Text Files", "*.txt*"),
                                                            ("All files", "*.*") )) #file dialog
            else:
                filez=tk.filedialog.askopenfilenames(title='Choose a single spectrum')
            if len(filez) > 0:
                lst=list(filez)
                for item in lst:
                    if '.list' in item or '.lst' in item :
                        self.stackplot=True
                        self.stackplottoggle()
                        self.listname=item
                        self.read_list()
                        for listitem in self.listedfiles:
                            self.database[listitem]={}
                    else:
                        self.database[item]={}

                self.norm_clear()
                self.loadSpectra()
                self.plotSpectra()
                self.plotRegions()
                self.stack=[]
                for row in self.database:
                    self.stack.append(row)
                if self.pane == True:
                    self.stackpane()
        else:
            if 'Darwin' not in OperatingSys:
                filez = tk.filedialog.askopenfilenames(title='Choose a list of spectra',filetypes=(("Spectra List", "*.list"),
                                                            ("Spectra List", "*.lst"),
                                                            ("Fits Files", "*.fit*"),
                                                            ("Fits Files", "*.FIT* "),
                                                            ("Text Files", "*.txt*"),
                                                            ("All files", "*.*") )) #file dialog
            else:
                filez=tk.filedialog.askopenfilenames(title='Choose a single spectrum')
            if len(filez) > 0:
                lst=list(filez)
                for item in lst:
                    if '.list' in item or '.lst' in item :
                        self.listname=item
                        self.read_list()
                        for listitem in self.listedfiles:
                            self.database[listitem]={}
                    else:
                        self.database[item]={}
                self.stackplot=True
                self.loadSpectra()
                self.stack=[]
                for row in self.database:
                    self.stack.append(row)
                self.stackplottoggle()



    def loadSpectra(self):
        for item in self.database:
            self.fname=item
            if '.fit' in item or '.FIT' in item:
                self.read_fits()
            elif '.txt' in item or '.TXT' in item:
                self.read_txt()


    def plotSpectra(self,spec=None):
        if spec == None:
            for i,f in enumerate(self.database):
                self.fname=f
                self.stackint=i
                self.wavelength=self.database[self.fname]['wavelength']
                self.flux=self.database[self.fname]['flux']
                self.flux_orig=self.database[self.fname]['flux_orig']
                try:
                    self.header=self.database[self.fname]['header']
                except:
                    pass
                self.splot()
        elif spec != None:
            self.fname=spec
            self.wavelength=self.database[self.fname]['wavelength']
            self.flux=self.database[self.fname]['flux']
            self.flux_orig=self.database[self.fname]['flux_orig']
            try:
                self.header=self.database[self.fname]['header']
            except:
                pass
            self.splot()

    def read_fits(self):
        """Reads a fits file into the dictionary of stored spectra."""
        #need to add a way to read multispec fits files
        #Based on Read a UVES spectrum from the ESO pipeline
        self.sp = fits.open(self.fname)
        header = self.sp[0].header

        if header['NAXIS'] == 1:
            wcs = WCS(header)
            index = np.arange(header['NAXIS1'])
            wavelength = wcs.wcs_pix2world(index[:,np.newaxis], 0)
            wavelength = wavelength.flatten()
            wavelength = wavelength*u.AA
            try:
                wavelength=wavelength/(1.0-header['VHELIO']/2.997925e05)
                self.database[self.fname]['heliocentric']=header['VHELIO']
            except:
                pass
            flux = self.sp[0].data*u.flx
            self.database[self.fname]['wavelength']=wavelength
            self.database[self.fname]['flux']=flux
            self.database[self.fname]['flux_orig']=flux
            self.database[self.fname]['header']=header
            self.sp.close()
        elif header['NAXIS'] > 1:
            wcs = WCS(header)
            index = np.arange(header['NAXIS1'])
            wavelength= index*u.pixel
            flux = self.sp[0].data[0].flatten()*u.flx
            self.database[self.fname]['wavelength']=wavelength
            self.database[self.fname]['flux']=flux
            self.database[self.fname]['flux_orig']=flux
            self.database[self.fname]['header']=header
            self.sp.close()
        else:
            tk.messagebox.showerror(title="Dimension Error",message="PySplot was only designed to work with 1D extracted spectra.")




    def read_txt(self):
            f1=open(self.fname,'r')
            data=np.array(list(csv.reader(f1,delimiter=' ')))
            wavelength=data[:,0].astype(float)
            try:
                flux=data[:,2].astype(float)
            except:
            #spec,=self.ax.plot(self.wavelength,self.flux)
            # self.ax.set_ylim([max(0,min(self.flux)),min(100,max(self.flux))])
                flux=data[:,1].astype(float)
            flux=flux*u.flx
            self.database[self.fname]['wavelength']=wavelength*u.AA
            self.database[self.fname]['flux']=flux
            self.database[self.fname]['flux_orig']=flux
            f1.close()


    def read_list(self):
        path=os.path.dirname(self.listname)
        f1=open(self.listname,'r')
        names=list(csv.reader(f1))
        self.listedfiles=[]
        for item in names:
            if '\\' in item[0] or '/' in item[0]:
                self.listedfiles.append(item[0])
            else:
                self.listedfiles.append(os.path.abspath(os.path.join(path,item[0])))
        f1.close()

    def generate_plot(self,exp=True):
        try:
            self.figframe.destroy()
        except:
            pass
        if exp == True:
            self.figframe=tk.Frame()
            self.figframe.pack(side="left", fill="both",expand=1)
        elif exp == False:
            self.figframe=tk.Frame()
            self.figframe.pack(side="left")
        else:
            print("Halp!")
        self.fig=plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.ax.tick_params(right= True,top= True,which="both")
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.figframe)
        self.canvas.get_tk_widget().pack(side="left", fill="both")
        try:
            self.toolbar = NavigationToolbar2TkAgg( self.canvas, self.figframe )
        except:
            self.toolbar = NavigationToolbar2Tk( self.canvas, self.figframe )

        self.ax.set_title("Single Spectra Mode",fontsize=12)
        self.toolbar.update()
        self.canvas.draw()
        self.canvas._tkcanvas.pack(side="top", fill="both", expand=1)

    def splot(self,xlim=None,ylim=None):
        #Normal plot mode
        if self.overplot == False and self.stackplot == False:
            self.ax.clear()
            self.ax.set_title("%s"%(self.fname),fontsize=10)
            spec,=self.ax.step(self.wavelength,self.flux)
            self.ax.set_ylabel("Flux")

        #Overplot Mode
        elif self.overplot == True and self.stackplot == False:
            self.output.delete(0,tk.END)
            self.ax.set_title("Overplot Mode, Display For Qualitative Comparison Only",fontsize=10)
            spec,=self.ax.step(self.wavelength,self.flux)
            self.ax.set_ylabel("Flux")

        elif self.overplot == False and self.stackplot == True:
            spec,=self.ax.plot(self.wavelength,self.flux+self.stackint*u.flx)
            self.ax.set_title("Stack Plot Mode, Display For Qualitative Comparison Only",fontsize=10)
            self.ax.set_ylabel("Stack Plot Flux")


        #Conflict
        elif self.overplot == True and self.stackplot == True:
            self.output.delete(0,tk.END)
            self.overplot=False
            self.stackplot=False
            tk.messagebox.showerror(title="Display Conflict",message="Overplot and stackplot can't be used at the same time, program reverting to single spectra mode.")
            self.restore()

        self.fig.suptitle("PySplot - Date: "+'{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now(),fontsize=10))
        if xlim != None:
            self.ax.set_xlim(xlim)
            self.ax.set_ylim(ylim)

        plt.grid(self.gridvalue)
        self.xaxislabel()
        if xlim == None:
            self.toolbar.update()
        self.canvas.draw()

    def xaxislabel(self):
        if self.wavelength.unit == 'km / s':
            self.ax.set_xlabel("Velocity (%s)"%self.wavelength.unit)
        else:
            if 'heliocentric' in self.database[self.fname]:
                    self.ax.set_xlabel("Heliocentric Wavelength (%s)"%self.wavelength.unit)
            else:
                self.ax.set_xlabel("Wavelength (%s)"%self.wavelength.unit)



    def zoomout(self,event=None):
        """Not accessible from menu.  Was designed to replot a region, but the matplotlib toolbar is better."""
        self.ax.set_xlim([min(self.wavelength.value),max(self.wavelength.value)])
        self.ax.set_ylim([min(self.flux),max(self.flux)])
        self.canvas.draw()


    def gridtoggle(self,event=None):
        """Toggles the plot grid on and off."""
        if self.gridvalue == False:
            self.gridvalue=True
        elif self.gridvalue == True:
            self.gridvalue=False
        else:
            pass
        plt.grid(self.gridvalue)
        self.canvas.draw()

    def overplottoggle(self,event=None):
        """Switches to overplot mode and replots the stack."""
        self.overplot=True
        self.stackplot=False
        self.ax.clear()
        self.ax.set_title("Overplot Plot Mode, Display For Qualitative Comparison Only",fontsize=12)
        self.plotSpectra()
        self.stackpane()
        self.plotRegions()


    def stackplottoggle(self,event=None):
        """Switches to stackplot mode and replots the stack."""
        self.output.delete(0,tk.END)
        self.stackplot=True
        self.overplot=False
        # self.stack=list(dict.fromkeys(self.stack)) #removes any duplicate file added to the list.
        self.ax.clear()
        self.ax.set_title("Stack Plot Mode, Display For Qualitative Comparison Only",fontsize=12)
        self.plotSpectra()
        self.stackpane()
        self.plotRegions()


    def singleplottoggle(self,event=None):
        """Single spectrum display mode, replots the last active spectrum."""
        self.output.delete(0,tk.END)
        self.stackplot=False
        self.overplot=False
        self.ax.set_title("Single Spectra Mode",fontsize=12)
        try:
            self.plotSpectra(spec=self.fname)
        except:
            pass
        try:
            self.plotSpectra(spec=self.stack[0])
        except:
            pass
        self.plotRegions()
        self.stackpane()

    def dynamical(self,event=None):
        self.output.delete(0,tk.END)

        t=tk.Toplevel(self.master,height=600,width=600)
        t.wm_title("Dynamical Spectrum")
        dynamicalframe=tk.Frame(t)
        dynamicalframe.pack(side="top", fill="both",expand=1)
        dynamical=plt.figure()
        ax = dynamical.add_subplot(111)
        ax.tick_params(right= True,top= True,which="both")
        canvas = FigureCanvasTkAgg(dynamical, master=dynamicalframe)
        canvas.get_tk_widget().pack(side="top", fill="both")
        try:
            self.toolbar = NavigationToolbar2TkAgg(canvas, dynamicalframe )
        except:
            self.toolbar = NavigationToolbar2Tk(canvas, dynamicalframe )

        jd=[]
        for i,row in enumerate(self.database):
            self.header=self.database[row]['header']
            try:
                jd.append(self.header['JD'])
            except:
                jd.append(i)
        wavelengthbase=self.database[self.stack[0]]['wavelength']
        x=len(wavelengthbase)

        dynamical.suptitle("PySplot - Date: "+'{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now(),fontsize=10))
        ax.set_ylabel("JD-%f3"%float(min(jd)))
        # self.xaxislabel()
        tt=np.arange(0, max(jd)-min(jd),.5) #phase steps

        data=np.empty([len(tt),x],dtype=float)
        data.fill(np.nan)
        jd=np.array(jd)-min(jd)
        for i,row in enumerate(self.database):
            flux=self.database[row]['flux']
            wavelength=self.database[row]['wavelength']
            temp_spec=np.empty(x,dtype=float)
            for j,val in enumerate(wavelength):
              temp_spec[j]=np.interp(val,x,flux)
            try:
                data[find_nearest_index(tt,self.jd[i])-1,:]=temp_spec
            except:
                pass
            data[find_nearest_index(tt,self.jd[i]),:]=temp_spec
            try:
                data[find_nearest_index(tt,self.jd[i])+1,:]=temp_spec
            except:
                pass
        # data needs to be resampled so all is on same sampling.
        # interpolate, but only along the time axis
        # for i,junk in enumerate(data[0,:]):
        #     data[:,i]=fill_nan(data[:,i])

        cmap=plt.get_cmap('Spectral')
        # cmap.set_under(color='white')
        # cbaxes = self.fig.add_axes([.88, .32, 0.03, .58])
        # ax.colorbar(cax = cbaxes).ax.tick_params(axis='y', direction='out')  #not sure where defined.

        ax.imshow(data,cmap=cmap, origin='lower',aspect='auto',alpha=1)#, extent=(self.wavelength[1]/u.Angstrom,self.wavelength[-1]/u.Angstrom,tt[0],tt[-1])),interpolation='nearest'


        self.toolbar.update()
        canvas.draw()


    def stackreset(self):
        """Reset the stacking parameters"""
        self.output.delete(0,tk.END)
        self.output.insert(tk.END,"Get started by opening a 1-D Spectrum or List of 1-D Spectra. File>open (or press o)")
        self.stack=[]
        self.database.clear()
        self.stackint=0
        del self.wavelength
        del self.flux
        del self.flux_orig
        del self.header
        self.loadednorm=False
        self.loadedbisect=False
        self.loadedregions=False
        self.ax.clear()
        self.canvas.draw()
        self.endoflog()
        if self.pane == True:
            self.stackpane()

    def currentspectra(self):
        print(self.database[self.fname].keys())

    def stackprint(self):
        """Prints stack to terminal."""
        print("Stack:")
        for f in self.database:
            print(f)
        print("Number in Stack: %s"%(str(len(self.database))))
        print("End Stack.")

    def savestack(self,name="stack.list",stack=None):
        """Save the stack for later re-loading."""
        if stack == None:
            stack=self.stack
        path=os.path.dirname(self.fname)
        basename=os.path.basename(name)
        savename=tk.filedialog.asksaveasfilename(initialdir=path,initialfile=basename, defaultextension=".list")
        dataout=open(savename,'w')
        for row in stack:
            dataout.write('%s\n'%(row))
        dataout.close()
        print("Stack saved to:  %s"%(savename))


    def smooth(self,event=None):
        self.output.delete(0,tk.END)
        self.output.insert(tk.END,"Enter an integer for the boxcar smoothing:")
        self.boxwidth=3
        self.flux=convolve(self.flux,Box1DKernel(self.boxwidth))
        self.output.delete(0,tk.END)
        self.splot()

    def measuremode(self,event=None):
        self.stackplot=False
        self.overplot=False

    def coord(self,event=None):
        """uses one click to print mouse position in the text box"""
        self.output.delete(0,tk.END)
        self.output.insert(tk.END,"Left Click (x) on a point to display the coordinates.")
        clicks= plt.ginput(1)
        clicks= np.array(clicks)
        x=clicks[0][0]
        y=clicks[0][1]
        print("Mouse Position: Wavelength = "+"{0.value:0.03f} {0.unit:FITS}".format(x)+", Flux = "+"{0.value:0.03f} {0.unit:FITS}".format(y))
        self.output.delete(0,tk.END)
        self.output.insert(tk.END,"Mouse Position: Wavelength = %s, Flux = %s" %(x,y))

    def velocity(self,event=None):
        """#uses one click to convert wavelength to velocity and vice versa"""
        self.measuremode()
        if self.wavelength.unit == u.AA:
            self.output.delete(0,tk.END)
            self.output.insert(tk.END,"Left Click (x) on where you want zero to be.")
            clicks= plt.ginput(1)
            clicks= np.array(clicks)
            self.wavecenter=clicks[0][0]*u.AA
            self.wavelength_backup=self.wavelength
            self.wavelength=(self.wavelength-self.wavecenter)/self.wavecenter*c.to('km/s')
            self.splot()
        elif self.wavelength.unit == 'km / s':
            self.wavelength=self.wavelength_backup
            self.splot()


    def restore(self,event=None):
        """Resets the last/current spectrum to the original state."""
        try:
            self.flux=self.flux_orig
            self.splot()
            self.norm_clear()
        except:
            pass

    def reset(self, event=None):
        """Just replots and clears the norm parameters.  A handy way to clear the graph of clutter."""
        if self.stackplot == False:
            if len(self.stack) > 0:
                self.splot()
            self.norm_clear()
        elif self.stackplot == True:
            self.stackplottoggle()
        self.region_clear()
        self.loadedbisect=False


    def pltregion(self,x,y,sym='-',c='black'):
        """takes two clicks, and plots a line between them"""
        self.ax.plot(x,y,sym,color=c)
        self.canvas.draw()

    def region(self,message=None):
        """uses two clicks to define a region for fitting or measuring."""
        if message == None:
            regionmessage="Left Click (x) on a point to the left and to the right of what you wish to select. Right Click (backspace) to remove a point."
        else:
            regionmessage=message
        self.output.delete(0,tk.END)
        self.output.insert(tk.END,regionmessage)
        clicks= plt.ginput(2)
        clicks= np.array(clicks)
        x=[]
        y=[]
        for i,j in enumerate(clicks):
            x.append(clicks[i][0])
            y.append(clicks[i][1])
        self.pltregion(x,y,sym='s',c='black')
        self.output.delete(0,tk.END)
        return x,y

    def points2fit(self):
        """collects many points. for fitting"""
        self.output.delete(0,tk.END)
        self.output.insert(tk.END,"Left Click (x) to add points along the continuum.  Right click (backspace) removes a point, Middle Click (enter) ends selection.")

        clicks= plt.ginput(0)
        clicks= np.array(clicks)
        x=[]
        y=[]
        for i,j in enumerate(clicks):
            x.append(clicks[i][0])
            y.append(clicks[i][1])
        self.pltregion(x,y,sym='s',c='red')
        self.output.delete(0,tk.END)
        return x,y

    def BisectLine(self,event=None):
        """Designed for measuring the center of very broad Wolf-Rayet star emission lines, may work for other situtaitons."""
        self.measuremode()
        print("enter bisection mode.")
        if self.loadedbisect == False:
            lx,ly=self.region(message="Select two points along the left edge of the feature.") #left edges of the feature to bisect
            rx,ry=self.region(message="Select two points along the right edge of the feature.") #right edges of the feature to bisect

            for xi in lx:
                self.x_norm.append(xi)
            for yi in ly:
                self.y_norm.append(yi)
            for xi in rx:
                self.x_norm.append(xi)
            for yi in ry:
                self.y_norm.append(yi)
        elif self.loadedbisect == True:
            print("Loaded bisection parameters.")
            lx=np.array([self.x_norm[0],self.x_norm[1]])
            rx=np.array([self.x_norm[2],self.x_norm[3]])
            ly=np.array([self.y_norm[0],self.y_norm[1]])
            ry=np.array([self.y_norm[2],self.y_norm[3]])

        rxg,ryg=self.chop(self.wavelength,self.flux,rx[0],rx[1])
        lxg,lyg=self.chop(self.wavelength,self.flux,lx[0],lx[1])

        center=(np.average(lxg)+np.average(rxg))/2.
        stderror=np.sqrt( (np.std(lxg)/np.sqrt(len(lxg)))**2 + (np.std(rxg)/np.sqrt(len(rxg)))**2 )
        # print("Bisected Center: %s , standard error: %s" %(center,stderror  ))

        outstring="Bisected Click Center: "+"{0.value:0.03f} {0.unit:FITS}".format(center)+\
                   ", Stnd Error: "+"{0.value:0.03f} {0.unit:FITS}".format(stderror)
        self.output.delete(0,tk.END)
        self.output.insert(tk.END,outstring)
        print(outstring)
        self.checklog()
        self.starlog.write(outstring)
        self.starlog.write('\n')
        self.ax.vlines(center.value,min(self.flux.value),max(self.flux.value))
        self.canvas.draw()



    def regionload(self,message=None):
        """Loads saved regions and gets the cloest values in the data"""
        if self.loadedregions == False :
            x,y=self.region(message) #click points to slice data
            self.saveregions_x=x
            self.saveregions_y=y
            self.loadedregions=True
        elif self.loadedregions == True:
            x,y=(self.saveregions_x,self.saveregions_y)
            # print("using loaded regions")
            # self.output.delete(0,tk.END)
            # self.output.insert(tk.END,"Using a loaded region. Use View>reset(r) to remove selected region.")
            # self.pltregion(x,y,sym='s',c='black')
            self.plotRegions()
        else:
            print("Unexpected Error in regionload()")
        xg,yg=self.chop(self.wavelength,self.flux,x[0],x[1])
        return xg,yg


    def eqw(self,event=None):
        """Measure equivalent width between two points IRAF style"""
        self.measuremode()
        xg,yg=self.regionload()
        continuum=(yg[0]+yg[-1])/2 #sets the continuum to the average of the left and right click.
        dwidth=[]
        for i,f in enumerate(yg):
          if i ==0:
            pass
          else:
            dwidth.append((1-f/continuum)*abs(xg[i]-xg[i-1]))
        width=sum(dwidth)
        bisect=(xg[0].value+xg[-1].value)/2.*xg[0].unit
        self.ax.vlines(bisect.value,min(yg.value),max(yg.value))
        self.canvas.draw()
        self.output.delete(0,tk.END)
        t=self.filedate()
        outstring=t+"Equivalent Width, "+"{0.value:0.03f}, {0.unit:FITS}".format(width)+\
                   ", Bisected Click Center, "+"{0.value:0.03f}, {0.unit:FITS}".format(bisect)
        self.output.insert(tk.END,outstring)
        print(outstring)
        self.checklog()
        self.starlog.write(outstring)
        self.starlog.write('\n')

    def filedate(self):
        basename=os.path.basename(self.fname)
        d=self.database[self.fname]
        try:
            jd=str(d['header']['JD'])
        except:
            jd=''
        try:
            hjd=str(d['header']['HJD'])
        except:
            hjd=''
        return "%s,JD, %s, HJD, %s,"%(basename,jd,hjd)



    def find_nearest_index(self,array,value):
        """Finds the closest data points to the input array."""
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx

    def chop(self,array1,array2,start1,stop1):
        """Cuts out a region of the spectrum for fitting or other functions that should only be run on a section of the spectrum."""
        small=min([start1,stop1])
        large=max([start1,stop1])
        startidx=self.find_nearest_index(array1,small)
        stopidx=self.find_nearest_index(array1,large)
        c1=array1[startidx:stopidx]
        c2=array2[startidx:stopidx]
        return c1,c2

    def scopy(self,event=None,script=None):
        """Copy out a section of a spectrum to a new spectrum"""
        self.measuremode()
        xg,yg=self.regionload()
        self.wavelength=xg
        self.flux=yg
        self.header['CRVAL1']=xg[0].value
        # self.database[self.fname]['crop-header']=self.header
        # self.database[self.fname]['crop-flux']=self.flux
        # self.database[self.fname]['crop-wave']=self.wavelength
        self.splot()
        if script == None:
            self.norm_clear()
            self.loadedregions=False
        else:
            self.plotRegions()

    def fit(self,func="gauss",event=None):
        """wrapper function for fitting line profiles"""
        # Fit the data using a Gaussian, Voigt, or Lorentzian profile
        self.measuremode()
        xg,yg=self.regionload() #get regions
        xgf=np.linspace(xg[0],xg[-1],len(xg)*3)
        ygf=np.interp(xgf,xg,yg)*u.flx

        # mid=len(yg)//2 #get guess for line peak by taking center of clicks
        xa=np.array([xgf[0].value,xgf[-1].value])
        ya=np.array([ygf[0].value,ygf[-1].value])
        linecoeff = np.polyfit(xa,ya,1) #does a linear polynomail fit to the data.
        ygfn=ygf/np.polyval(linecoeff,xgf.value)#*u.flx #normalized
        ygn=yg/np.polyval(linecoeff,xg/xg[0])*u.flx #xg/xg[0] to remove the unit
        invert=False
        reflevel=yg[0]
        if yg[len(yg)//2] < reflevel:
            ygf=reflevel-ygf
            invert=True
        # self.ax.clear()

        if func=="gauss":
          g_init = models.Gaussian1D(amplitude=np.max(yg), mean=np.mean(xg), stddev=xg[2]-xg[0])
          fit_g = fitting.LevMarLSQFitter()
          if invert == False:
            g = fit_g(g_init, xg, yg)
            plt.plot(xgf, g(xgf), label='%s'%(func))#*np.polyval(linecoeff,xgf*xg[0]/xg[0].value)
            amp=g.amplitude
            # self.ax.plot(xgf,ygfn)
            # self.ax.plot(xgf,g(xgf))
          elif invert == True:
            g = fit_g(g_init, xg, reflevel-yg)
            plt.plot(xgf, (reflevel-g(xgf)), label='%s'%(func))
            amp=reflevel-g.amplitude
            # self.ax.plot(xgf,-ygfn)
            # self.ax.plot(xgf,-g(xgf))
          t=self.filedate()
          outstring=t+"Gaussian Center, "+"{0.value:0.03f}, {0.unit:FITS}".format(g.mean)+", FWHM, "+ \
                           "{0.value:0.03f}, {0.unit:FITS}".format(g.fwhm)+", Amplitude, "+"{0.value:0.03f}, {0.unit:FITS}".format(amp)
          print(outstring)
          self.checklog()
          self.starlog.write(outstring)
          self.starlog.write('\n')
          self.output.delete(0,tk.END)
          self.output.insert(tk.END, outstring)


        elif func=="voigt":
          g_init = models.Voigt1D(x_0=np.mean(xg),amplitude_L=np.max(yg)-reflevel , fwhm_L=xg[1]-xg[0], fwhm_G=xg[1]-xg[0])
          fit_g = fitting.LevMarLSQFitter()
          if invert == False:
            g = fit_g(g_init, xg, yg)
            plt.plot(xgf, g(xgf))#*np.polyval(linecoeff,xgf), label='%s'%(func))
            amp=g.amplitude_L
            # self.ax.plot(xg,ygn)
            # self.ax.plot(xgf,g(xgf))
          elif invert == True:
            g = fit_g(g_init, xg, reflevel-yg)
            plt.plot(xgf, (reflevel-g(xgf)))#*np.polyval(linecoeff,xgf), label='%s'%(func))
            # self.ax.plot(xg,-ygn)
            # self.ax.plot(xgf,-g(xgf))
            amp=reflevel-g.amplitude_L
          t=self.filedate()
          outstring=t+"Voigt Center, "+"{0.value:0.03f}, {0.unit:FITS}".format(g.x_0)+", Lorentzian_FWHM, "+\
                  "{0.value:0.03f}, {0.unit:FITS}".format(g.fwhm_L)+\
                  ", Gaussian_FWHM, "+"{0.value:0.03f}, {0.unit:FITS}".format(g.fwhm_G)+\
                  ", Amplitude, "+"{0.value:0.03f}, {0.unit:FITS}".format(amp)
          print(outstring)
          self.checklog()
          self.starlog.write(outstring)
          self.starlog.write('\n')
          self.output.delete(0,tk.END)
          self.output.insert(tk.END,outstring)

        elif func=="lorentz":
          g_init = models.Lorentz1D(x_0=np.mean(xg),amplitude=np.max(yg)-reflevel, fwhm=xg[1]-xg[0])
          fit_g = fitting.LevMarLSQFitter()
          if invert == False:
            g = fit_g(g_init, xg, yg)
            plt.plot(xgf, g(xgf))#*np.polyval(linecoeff,xgf), label='%s'%(func))
            # self.ax.plot(xg,ygn)
            # self.ax.plot(xgf,g(xgf))
            amp=g.amplitude
          elif invert == True:
            g = fit_g(g_init, xg, reflevel-yg)
            plt.plot(xgf, (reflevel-g(xgf)))#*np.polyval(linecoeff,xgf), label='%s'%(func))
            # self.ax.plot(xg,-ygn)
            # self.ax.plot(xgf,-g(xgf))
            amp=reflevel-g.amplitude

          t=self.filedate()
          outstring=t+"Lorentz Center, "+"{0.value:0.03f}, {0.unit:FITS}".format(g.x_0)+", FWHM, "+\
                  "{0.value:0.03f}, {0.unit:FITS}".format(g.fwhm)+\
                  ", Amplitude, "+"{0.value:0.03f}, {0.unit:FITS}".format(amp)
          print(outstring)
          self.checklog()
          self.starlog.write(outstring)
          self.starlog.write('\n')
          self.output.delete(0,tk.END)
          self.output.insert(tk.END,outstring)


        else:
          print("Error with Fitting Function Selection")

        self.canvas.draw()

    def SaveRegion(self):
        """Save the regions used for equivalent width measurements and for fitting line profiles."""
        path=os.path.dirname(self.fname)
        basename=os.path.basename("region.par")
        savename=tk.filedialog.asksaveasfilename(initialdir=path,initialfile=basename, defaultextension=".par")
        dataout=open(savename,'w')
        dataout.write('%s\n'%(self.fname))
        for row in self.saveregions_x:
            dataout.write('%s '%(row))
        dataout.write('\n')
        for row in self.saveregions_y:
            dataout.write('%s '%(row))
        dataout.write('\n')
        dataout.close()
        print("Region saved to:  %s"%(savename))
        self.norm_clear()

    def plotRegions(self):
        message="Using a loaded region, when finished use View>Reset or press r."
        if self.loadednorm == True:
            self.ax.plot(self.x_norm,self.y_norm,'s',color='black')
            self.output.delete(0,tk.END)
            self.output.insert(tk.END,message)
        if self.loadedbisect == True:
            self.ax.plot(self.x_norm,self.y_norm,'s',color='black')
            self.output.delete(0,tk.END)
            self.output.insert(tk.END,message)
        if self.loadedregions == True:
            self.ax.plot(self.saveregions_x,self.saveregions_y,'s',color='black')
            self.output.delete(0,tk.END)
            self.output.insert(tk.END,message)
        self.canvas.draw()


    def LoadRegion(self):
        """Load the regions used for equivalent width measurements and for fitting line profiles."""
        self.norm_clear()
        if 'Darwin' not in OperatingSys:
            file=tk.filedialog.askopenfilename(title='Choose a region parameter file (.par)',filetypes=(("Parameter", "*.par"),
                                                        ("All files", "*.*") ))
        else:
            file=tk.filedialog.askopenfilename(title='Choose a region parameter file (.par)')
        dataout=open(file)
        for i,line in enumerate(dataout):
            if i == 0 :
                pass
            elif i == 1 :
                self.saveregions_x=np.array(line.split(),dtype=float)
                # print(line.split()[0])
            elif i == 2 :
                self.saveregions_y=np.array(line.split(),dtype=float)
            else:
                pass
        dataout.close()
        self.loadedregions=True
        self.plotRegions()

    def SaveBisect(self):
        """Save the regions for bisection"""
        path=os.path.dirname(self.fname)
        basename=os.path.basename("bisect.par")
        savename=tk.filedialog.asksaveasfilename(initialdir=path,initialfile=basename, defaultextension=".par")
        dataout=open(savename,'w')
        dataout.write('%s\n'%(self.fname))
        for row in self.x_norm:
            dataout.write('%s '%(row))
        dataout.write('\n')
        for row in self.y_norm:
            dataout.write('%s '%(row))
        dataout.write('\n')
        dataout.close()
        print("Bisection Regions Saved to:  %s"%(savename))
        self.norm_clear()

    def LoadBisect(self):
        """Load the regions bisection."""
        self.norm_clear()
        if 'Darwin' not in OperatingSys:
            file=tk.filedialog.askopenfilename(title='Choose a bisection parameter file (.par)',filetypes=(("Parameter", "*.par"),
                                                        ("All files", "*.*") ))
        else:
            file=tk.filedialog.askopenfilename(title='Choose a bisection parameter file (.par)')
        dataout=open(file)

        for i,line in enumerate(dataout):
            if i == 0 :
                pass
            elif i == 1 :
                self.x_norm=np.array(line.split(),dtype=float)
                # print(line.split()[0])
            elif i == 2 :
                self.y_norm=np.array(line.split(),dtype=float)
            else:
                pass
        dataout.close()
        self.loadedbisect=True
        self.plotRegions()


    def norm_clear(self,event=None):
        """Clear and reset normalization parameters and region selections."""
        self.x_norm=[]
        self.y_norm=[]
        self.goodfit=False
        self.order=1
        self.output.delete(0,tk.END)

    def region_clear(self,event=None):
        self.output.delete(0,tk.END)
        self.saveregions_x=[0,0]
        self.saveregions_y=[0,0]
        self.loadedregions=False
        xlim=self.ax.get_xlim()
        ylim=self.ax.get_ylim()

        self.splot(xlim=xlim,ylim=ylim)


    def SaveNorm(self):
        """Save the regions and powerlaw for the normalization."""
        path=os.path.dirname(self.fname)
        basename=os.path.basename("norm.par")
        savename=tk.filedialog.asksaveasfilename(initialdir=path,initialfile=basename, defaultextension=".par")
        dataout=open(savename,'w')
        dataout.write('%s\n'%(self.fname))
        dataout.write('%s\n'%(self.order))
        for row in self.x_norm:
            dataout.write('%s '%(row))
        dataout.write('\n')
        for row in self.y_norm:
            dataout.write('%s '%(row))
        dataout.write('\n')
        dataout.close()
        print("Normalization Parameters Saved to:  %s"%(savename))
        self.norm_clear()

    def LoadNorm(self):
        """Load the regions and powerlaw for the normalization."""
        self.norm_clear()
        if 'Darwin' not in OperatingSys:
            file=tk.filedialog.askopenfilename(title='Choose a normalization parameter file (.par)',filetypes=(("Parameter", "*.par"),
                                                        ("All files", "*.*") ))
        else:
            file=tk.filedialog.askopenfilename(title='Choose a normalization parameter file (.par)')
        dataout=open(file)

        for i,line in enumerate(dataout):
            if i == 0 :
                pass
            elif i == 1 :
                self.order=int(line)
            elif i == 2 :
                self.x_norm=np.array(line.split(),dtype=float)
                # print(line.split()[0])
            elif i == 3 :
                self.y_norm=np.array(line.split(),dtype=float)
            else:
                pass
        dataout.close()
        self.loadednorm=True
        self.plotRegions()


    def continuum(self,event=None):
        """Create array of click to use as the continuum."""
        self.measuremode()
        x,y=self.region() #click points to slice data
        self.pltregion(x,y,sym='-s',c='black')
        xg,yg=self.chop(self.wavelength,self.flux,x[0],x[1])
        for xi in xg:
            self.x_norm.append(xi.value)
        for yi in yg:
            self.y_norm.append(yi.value)
        self.loadedregions=True

    def normalize(self,event=None,script=None):
        """Continuum normalize by using selected points as continuum."""
        self.measuremode()
        xn=np.array(self.x_norm)
        yn=np.array(self.y_norm)
        if len(self.x_norm) < 2:
            self.output.delete(0,tk.END)
            self.output.insert(tk.END,"Please select region(s) for the continuum before useing normalize to fit.  Use Set Continuum or the s-key to set ranges.")
        else:
            self.output.delete(0,tk.END)
            self.output.insert(tk.END,"The integer below will be used for the order of the polynomial fit.")
            self.order=int(self.w1.get())
            while self.goodfit == False:
                linecoeff = np.polyfit(xn,yn,self.order)
                nflux=self.flux/np.polyval(linecoeff,self.wavelength.value)
                self.normtest,=self.ax.plot(self.wavelength,nflux)
                self.ax.set_ylim([max(0,min(nflux.value)),min(100,max(nflux.value))])
                self.canvas.draw()
                if script == None:
                    answer=tk.messagebox.askyesno("Question","Proceed with the fit?")
                    if answer == True:
                       self.goodfit = True
                       self.ax.lines.remove(self.normtest)
                       self.output.delete(0,tk.END)
                       self.flux=nflux
                       self.splot()

                    else:
                       self.ax.lines.remove(self.normtest)
                       self.output.delete(0,tk.END)
                       self.output.insert(tk.END,"Change the integer below and re-run the normalize command.")
                       goodfit = False
                       break
                else:
                    self.goodfit = True
                    self.ax.lines.remove(self.normtest)
                    self.output.delete(0,tk.END)
                    self.flux=nflux
                    self.splot()

    def abort_stack(self):
        self.output.delete(0,tk.END)
        self.output.insert(tk.END,"Please select region.")
        print("No region loaded.")

    def stacker(self,func=None):
        self.stackforsaving=[]
        """A helper function to automate tasks for many spectra."""
        if func == None:
            print("No function selected for stacker.")
        else:
            for f in self.database:
                self.fname=f
                self.plotSpectra(spec=self.fname)
                self.measuremode()
                if func == "norm":
                    if self.loadednorm == True:
                        self.normalize(script='Yes')
                        self.goodfit = False
                        self.save_fits(extend='-norm.fits')
                    else:
                        self.output.delete(0,tk.END)
                        self.output.insert(tk.END,"Please select continuum (s) and set the fit order for the continuum before useing normalize to fit.")
                        print("No normalization parameters loaded.")

                elif func == "scopy":
                    if self.loadedregions == True:
                        self.scopy(script='Yes')
                        self.save_fits(extend='-crop.fits')
                    else:
                        self.abort_stack()
                elif func == "eqw":
                    if self.loadedregions == True:
                        self.eqw()
                    else:
                        self.abort_stack()
                elif func == "gauss":
                    if self.loadedregions == True:
                        self.fit(func="gauss")
                    else:
                        self.abort_stack()
                elif func == "voigt":
                    if self.loadedregions == True:
                        self.fit(func="voigt")
                    else:
                        self.abort_stack()
                elif func == "lorentz":
                    if self.loadedregions == True:
                        self.fit(func="lorentz")
                    else:
                        self.abort_stack()
                elif func == "bisect":
                    if self.loadedbisect == True:
                        self.BisectLine()
                    else:
                        self.abort_stack()
                else:
                    print("Stacker cannot handle a function.")
            if func == "norm":
                self.savestack(name="norm.list",stack=self.stackforsaving)
            elif func == "scopy":
                self.savestack(name="crop.list",stack=self.stackforsaving)

    def stackwindowplot(self,spec):
        self.singleplottoggle()
        self.plotSpectra(spec=self.fname)
        self.plotRegions()

    def stackwindowopenspectra(self):
        self.openSpectra()

    def stackwindowclearstack(self):
        self.stackreset()

    def removefromstack(self):
        if self.stackplot==True:
            pass
        else:
            self.stack.remove(self.fname)
            self.database.pop(self.fname)
            self.ax.clear()
            self.canvas.draw()
            self.stackpane()
            self.output.delete(0,tk.END)
            self.output.insert(tk.END,"You may replot the stack (]), overplot ([), or choose a new spectrum from the stack.")

    def hidepane(self,event=None):
        self.pane=False
        self.canvasframe.destroy()
        self.stackcanvas.destroy()
        self.buttonframe.destroy()

        if self.stackplot == True or self.overplot == True:
            self.plotSpectra()
        else:
            self.plotSpectra(spec=self.fname)


    def stackpane(self,event=None):
        self.pane=True
        try:
            self.canvasframe.destroy()
        except:
            pass
        try:
            self.stackcanvas.destroy()
        except:
            pass
        try:
            self.vsb.destroy()
        except:
            pass
        try:
            self.buttonframe.destroy()
        except:
            pass

        self.buttonframe=tk.Frame()
        self.buttonframe.pack(side="top")
        tk.Button(self.buttonframe,text="Remove Selected", command=self.removefromstack).pack(side="left",padx=2)
        tk.Button(self.buttonframe,text="Hide Pane (})", command=self.hidepane).pack(side="left",padx=2)
        self.canvasframe=tk.Frame()
        self.canvasframe.pack(side = "top", fill="both",expand="yes")
        self.stackcanvas = tk.Canvas(self.canvasframe,background="#ffffff")
        self.stackcanvas.pack(side="top", fill="both",expand="yes")
        self.stackframe=tk.Frame(self.stackcanvas)
        self.stackframe.pack(side = "left", fill="both",expand="yes")

        self.stackcanvas.create_window((0,0), window=self.stackframe, anchor="nw",tags="self.stackframe")

        self.vsb = tk.Scrollbar(self.stackcanvas, orient="vertical", command=self.stackcanvas.yview)
        self.hsb = tk.Scrollbar(self.stackcanvas, orient="horizontal", command=self.stackcanvas.xview)
        self.stackcanvas.configure(yscrollcommand=self.vsb.set,xscrollcommand=self.hsb.set)
        self.vsb.pack(side="right", fill="y")
        self.hsb.pack(side="bottom", fill="x")

        self.canvasframe.bind("<Configure>", self.onFrameConfigure)


        specnumber=list(np.arange(len(self.database))+1)
        specnumber.reverse()
        i=0
        for i,s in enumerate(self.stack[::-1]): #the syntax [::-1] reverses the list without modifying it so that  the button list is the same vertical order as the stack plotted spectra.
            t="%s"%(specnumber[i])
            tk.Button(self.stackframe,text="%s: %s"%(t,os.path.basename(s)),command=partial(self.stackwindowplot,s)).pack(side="top")



    def onFrameConfigure(self, event):
        '''Reset the scroll region to encompass the inner frame'''
        self.stackcanvas.configure(scrollregion=self.stackcanvas.bbox("all"))


    def save_fits(self,extend=None):
        """Save a new fits file with header."""
        try:
            hdu = fits.PrimaryHDU(self.flux.value,self.header)
            path=os.path.dirname(self.fname)
            basename=os.path.basename(self.fname)
            path_wo_ext=os.path.splitext(self.fname)[0]
            if extend == None:
                savename=tk.filedialog.asksaveasfilename(initialdir='./',initialfile=basename, defaultextension=".fits")
                hdu.writeto(savename)
            else:
                savename=os.path.join(path_wo_ext,basename,extend)
                # savename=path_wo_ext+extend #probably need to use an os.join here
                self.stackforsaving.append(savename)
                hdu.writeto(savename,overwrite=True)
                print("Saved to: ", savename)
        except:
            pass

    def save1DText(self):
        """Save a headerless text spectrum."""
        try:
            savename=tk.filedialog.asksaveasfilename(initialfile=self.fname,defaultextension=".txt")
            dataout=open(savename,'w')
            for i,val in enumerate(self.flux):
                dataout.write('%s %s\n'%(self.wavelength[i].value,self.flux[i]))
            w.destroy()
            dataout.close()
        except:
            pass

    def _quit(self):
        try:
            self.endoflog()
        except:
            pass
        self.master.destroy()  # this is necessary on Windows to prevent
                        # Fatal Python Error: PyEval_RestoreThread: NULL tstate
        self.master.quit()     # stops mainloop


    def imhead(self,event=None):
        """Display the image header as loaded from the file."""
        try:
            self.headerwin=tk.Toplevel(self.master,height=600,width=600)
            self.headerwin.wm_title("FITS Header:  %s"%(self.fname))

            keywordframe=tk.Frame(self.headerwin)
            keywordframe.pack(side="top")
            l1=tk.Label(keywordframe, text="Header Keyword:").pack( side = "left")
            self.entryvar1=tk.StringVar()
            self.header_keyword = tk.Entry(keywordframe,text=self.entryvar1,width=60)
            self.header_keyword.pack(side = "left")
            self.entryvar1.set('')
            b = tk.Button(keywordframe,text="Set Value", command=self.header_set)
            b.pack(side = "left")

            valueframe=tk.Frame(self.headerwin)
            valueframe.pack(side="top")
            l2=tk.Label(valueframe, text="Value:  ").pack( side = "left")
            self.entryvar2=tk.StringVar()
            self.header_value = tk.Entry(valueframe,text=self.entryvar2,width=60)
            self.header_value.pack(side = "left")
            self.entryvar2.set('')

            b2 = tk.Button(valueframe,text="Close Window", command=lambda: self.destroychild(self.headerwin))
            b2.pack(side = "left")

            self.header_list()

            self.header_disp.bind("<Double-Button-1>",self.grab_keyword)

        except:
            self.destroychild(self.headerwin)

    def grab_keyword(self,event=None):
        self.entryvar1.set(self.header_disp.get(tk.ACTIVE).split()[0])

    def header_list(self,event=None):
        self.listframe=tk.Frame(self.headerwin)
        self.listframe.pack(side="top")#,fill="both",exapand=1)
        s=tk.Scrollbar(self.listframe)
        s.pack(side="right",fill=tk.Y)#,expand=1)
        self.header_disp=tk.Listbox(self.listframe,yscrollcommand=s.set,height=40,width=80)
        for keys in self.header.tostring(sep=',').split(','):
            self.header_disp.insert(tk.END,"%s "%(keys))

        s.config(command=self.header_disp.yview)
        self.header_disp.pack(side=tk.LEFT,fill=tk.BOTH,expand=1)


    def header_set(self,event=None):
        """Write the new keyword/value to the header stored in the dictionary."""
        self.header[self.header_keyword.get()]=self.header_value.get()
        self.entryvar1.set('')
        self.entryvar2.set('')
        if self.listframe is not None:
            self.listframe.destroy()
            self.header_list()


    def destroychild(self,w):
        w.destroy()

    def About(self):
        t=tk.Toplevel(self.master,height=600,width=600)
        t.wm_title("About")
        tk.Label(t,text="Author: Dr. Joshua Thomas \n Clarkson University \n jthomas@clarkson.edu \n thomas.joshd@gmail.com").pack()
        tk.Label(t,text="This program was designed to emulate some basic IRAF splot functions.").pack()
        tk.Label(t,text="Astropy, scipy, and numpy libraries are used.").pack()
        tk.Label(t,text="Version %s"%version).pack()
        tk.Label(t,text="Last Updated %s"%UPDATED).pack()
        link2=tk.Label(t,text="Check for updates at:  https://github.com/thomasjoshd/pysplot", fg="blue", cursor="hand2")
        link2.pack()
        link2.bind("<Button-1>", lambda e: self.callback("https://github.com/thomasjoshd/pysplot"))
        b = tk.Button(t,text="Close", command=lambda: self.destroychild(t))
        b.pack()

    def callback(self,url):
        webbrowser.open_new(url)

#----------------------------------------------------------------------------------
#Begin GUI



root = tk.Tk() #main GUI window
program=App(root)
root.protocol("WM_DELETE_WINDOW", program._quit)
root.mainloop() #lets the GUI run

#need a way to stort the spectra in a stack by date, will require a way to examine the header and store to a list.
#this will be nice for stack plots, but will also allow dynamical spectra.

#implement mouse scroll wheel.

#right now heliocentric correction is done when the file is read, which doesn't allow for updating if the keyword is changed in the new header editing mode.
