"""Author: Dr. Joshua Thomas
thomas.joshd@gmail.com
This program was designed to emulate some basic IRAF splot functions.
Tested on Python 3.6.5 and 3.6.7, Linux Mint 19.1 and Windows 10.
Uses Astropy library, and some parts are directly modified from the UVES tutorial."""
UPDATED="16-APR-2019"
version="0.3-alpha"
from functools import partial
import numpy as np #arrays and math
import csv
import os
import tkinter as tk
import tkinter.messagebox
from tkinter.filedialog import askopenfilename,askopenfilenames,asksaveasfilename
from functools import partial

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
from astropy.modeling.functional_models import Voigt1D,Lorentz1D

import datetime
#import zoom

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

class App:
    def __init__(self,master):
        self.master=master
        master.title("PySplot, Version %s"%str(version))

        tk.Label(self.master, text="Prompt:").pack( side = tk.TOP)

        self.output=tk.Entry(self.master,width=100)
        self.output.pack(side = tk.TOP)

        tk.Label(self.master, text="Input:").pack( side = tk.TOP)

        self.w1v=tk.StringVar()
        self.w1 = tk.Entry(self.master,text=self.w1v,width=10)
        self.w1.pack(side = tk.TOP)
        self.w1v.set(1)	#default value of so the program doesn't barf if accidentally pushed.

        #--------------------------------------------------
        #file menu
        menu = tk.Menu(self.master)
        master.config(menu=menu)

        filemenu = tk.Menu(menu)
        menu.add_cascade(label="File", menu=filemenu)
        filemenu.add_command(label="Open 1-D Fits (o)",command=self.openSpectra)
        filemenu.add_command(label="Restore Original (r)",command=self.restore)
        filemenu.add_command(label="Image Header (h)",command=self.imhead)
        filemenu.add_separator()
        filemenu.add_command(label="Save New Fits",command=self.save_fits)
        filemenu.add_separator()
        filemenu.add_command(label="Save Headerless Text Spectra",command=self.save1DText)
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=self._quit)

        viewmenu = tk.Menu(menu)
        menu.add_cascade(label="View", menu=viewmenu)
        viewmenu.add_command(label="Reset View",command=self.splot)
        viewmenu.add_command(label="Grid Toggle (|)",command=self.gridtoggle)
        viewmenu.add_command(label="Over Plot ([)",command=self.overplottoggle)
        viewmenu.add_command(label="Stack Plot (])",command=self.stackplottoggle)
        viewmenu.add_command(label="Zoom to Fit (z)",command=self.zoomout)

        modmenu = tk.Menu(menu)
        menu.add_cascade(label="Modify", menu=modmenu)
        modmenu.add_command(label="Set Continuum (s)", command=self.continuum)
        modmenu.add_command(label="Normalize (t)", command=self.normalize)
        modmenu.add_command(label="Reset Normalization parameters (q)", command=self.norm_clear)
        modmenu.add_separator()
        modmenu.add_command(label="Cut Spectra (c)", command=self.scopy)
##        modmenu.add_command(label="Trim Spectra (type limits)", command=self.scopy_manual)

        fitmenu = tk.Menu(menu)
        menu.add_cascade(label="Measure", menu=fitmenu)
        fitmenu.add_command(label="Equivalent Width (e)", command=self.eqw)
        fitmenu.add_separator()
        fitmenu.add_command(label="Gaussian (g)", command=partial(self.fit,func="gauss"))
        fitmenu.add_command(label="Voigt (v)", command=partial(self.fit,func="voigt"))
        fitmenu.add_command(label="Lorentzian (l)", command=partial(self.fit,func="lorentz"))
        fitmenu.add_separator()
        fitmenu.add_command(label="Convert Wavelength <-> Velocity (u)", command=self.velocity)



        helpmenu = tk.Menu(menu)
        menu.add_cascade(label="Help", menu=helpmenu)
        helpmenu.add_command(label="About", command=self.About)

        #--------------------------------------------------
        self.output.insert(tk.END,"Get started by opening a 1-D Spectrum or List of 1-D Spectra.")
        self.generate_plot()
        self.gridvalue=True
        self.overplot=False
        self.stackplot=False
        self.stackint=0

##        self.ZP=zoom.ZoomPan()
        self.master.bind('e', self.eqw)
        self.master.bind('o', self.openSpectra)
        self.master.bind('h', self.imhead)
        self.master.bind('r', self.restore)
        self.master.bind('z', self.zoomout)
        self.master.bind('s', self.continuum)
        self.master.bind('t', self.normalize)
        self.master.bind('q', self.norm_clear)
        self.master.bind('g', partial(self.fit,"gauss"))
        self.master.bind('v', partial(self.fit,"voigt"))
        self.master.bind('l', partial(self.fit,"lorentz"))
        self.master.bind('|', self.gridtoggle)
        self.master.bind('[', self.overplottoggle)
        self.master.bind(']', self.stackplottoggle)
        self.master.bind('<space>', self.coord)
        self.master.bind('u', self.velocity)
        self.master.bind('c',self.scopy)

    def openSpectra(self,event=None):
        filez = askopenfilenames(title='Choose a list of spectra',filetypes=(("Fits Files", "*.fit*"),
                                                            ("Fits Files", "*.FIT* "),
                                                            ("Text Files", "*.txt*"),
                                                            ("Spectra List", "*.list"),
                                                            ("Spectra List", "*.lst"),
                                                            ("All files", "*.*") )) #file dialog
        lst = list(filez)
##        if len(lst) > 1:
##            self.overplot=True
##            self.stackplot=False
        for item in lst:
            self.fname=item
            print(self.fname)
            if '.fit' in item or '.FIT' in item:
                self.read_fits()
                self.splot()
            elif '.txt' in item or '.TXT' in item:
                self.read_txt()
                self.splot()
            elif '.list' in item or '.lst' in item :
                self.overplot=True
                self.stackplot=False
                self.listname=item
                self.read_list()
        self.norm_clear()



    def read_fits(self):
        #Based on Read a UVES spectrum from the ESO pipeline
        self.sp = fits.open(self.fname)
        header = self.sp[0].header

        wcs = WCS(header)
        #make index array
        index = np.arange(header['NAXIS1'])

        wavelength = wcs.wcs_pix2world(index[:,np.newaxis], 0)
        wavelength = wavelength.flatten()
        wavelength = wavelength*u.AA
        flux = self.sp[0].data
        self.wavelength=wavelength
        self.flux=flux
        self.header=header
        self.flux_orig=flux


    def read_txt(self):
            data=np.array(list(csv.reader(open(self.fname,'r'),delimiter=' ')))
            wavelength=data[:,0].astype(float)
            self.wavelength=wavelength*u.AA
            try:
                flux=data[:,2].astype(float)
            except:
                flux=data[:,1].astype(float)
            self.flux=flux
            self.flux_orig=flux


    def read_list(self):
        path=os.path.dirname(self.listname)
        names=list(csv.reader(open(self.listname,'r')))
        for item in names:
            if '\\' in item[0] or '/' in item[0]:
                self.fname=item[0]
            else:
                self.fname=os.path.abspath(os.path.join(path,item[0]))
            if '.txt' in self.fname or '.TXT' in self.fname:
                self.read_txt()
            elif '.fit' in self.fname or '.FIT' in self.fname:
                self.read_fits()

    def generate_plot(self):
        self.fig=plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.ax.tick_params(right= True,top= True,which='both')
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.master)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        try:
            self.toolbar = NavigationToolbar2TkAgg( self.canvas, self.master )
        except:
            self.toolbar = NavigationToolbar2Tk( self.canvas, self.master )

        self.toolbar.update()
        self.canvas.draw()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    def splot(self):
        if self.overplot == False and self.stackplot == False:
            self.ax.clear()
            self.ax.set_title("%s"%(self.fname),fontsize=12)
            spec,=self.ax.step(self.wavelength,self.flux)
            #spec,=self.ax.plot(self.wavelength,self.flux)
            self.ax.set_ylim([max(0,min(self.flux)),min(100,max(self.flux))])

        elif self.overplot == True and self.stackplot == False:
            self.output.delete(0,tk.END)
            self.ax.set_title("Overplot Mode, For Qualitative Comparison Only",fontsize=12)
            spec,=self.ax.step(self.wavelength,self.flux)
##            self.ax.set_ylim([max(0,min(self.flux)),min(100,max(self.flux))])

        elif self.overplot == False and self.stackplot == True:
            self.output.delete(0,tk.END)
            self.output.insert(tk.END,"The number entered below will change the stack spacing.")
##            self.ax.set_ylim([max(0,min(self.flux)),min(100,max(self.flux))])
            spec,=self.ax.plot(self.wavelength,self.flux+self.stackint)
            self.stackint=self.stackint+float(self.w1.get())
            self.ax.set_title("Stack Plot Mode, For Qualitative Comparison Only",fontsize=12)

        elif self.overplot == True and self.stackplot == True:
            self.output.delete(0,tk.END)
            self.overplot=False
            self.stackplot=False
            tkinter.messagebox.showerror(title="Display Conflict",message="Overplot and stackplot can't be used at the same time, program reverting to single spectra mode.")
            self.restore()

        self.ax.set_xlabel("%s"%self.wavelength.unit)
        self.ax.set_ylabel("Flux")
        self.fig.suptitle("PySplot - Date: "+'{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now(),fontsize=12))
        plt.grid(self.gridvalue)
        self.toolbar.update()
        self.canvas.draw()

    def zoomout(self,event=None):
        self.ax.set_xlim([min(self.wavelength.value),max(self.wavelength.value)])
        self.ax.set_ylim([min(self.flux),max(self.flux)])
        self.canvas.draw()


    def gridtoggle(self,event=None):
        if self.gridvalue == False:
            self.gridvalue=True
        elif self.gridvalue == True:
            self.gridvalue=False
        else:
            pass
        plt.grid(self.gridvalue)
        self.canvas.draw()

    def overplottoggle(self,event=None):
        if self.overplot == False:
            self.overplot=True
        elif self.overplot == True:
            self.overplot = False
        else:
            pass
        self.splot()

    def stackplottoggle(self,event=None):
        if self.stackplot == False:
            self.stackplot=True
        elif self.stackplot == True:
            self.stackplot = False
        else:
            pass
        self.splot()

    def measuremode(self,event=None):
        self.stackplot=False
        self.overplot=False
        #self.splot()

    def coord(self,event=None):
        #uses one click to print mouse position in the text box
        self.output.delete(0,tk.END)
        self.output.insert(tk.END,"Left Click (x) on a point to the left and to the right of what you wish to measure. Right Click (backspace) to remove a point.")
        clicks= plt.ginput(1)
        clicks= np.array(clicks)
        x=clicks[0][0]
        y=clicks[0][1]
        print("Mouse Position: Wavelength = "+"{0.value:0.03f} {0.unit:FITS}".format(x)+", Flux = "+"{0.value:0.03f} {0.unit:FITS}".format(y))
        self.output.delete(0,tk.END)
        self.output.insert(tk.END,"Mouse Position: Wavelength = %s, Flux = %s" %(x,y))

    def velocity(self,event=None):
        self.measuremode()
        #uses one click to convert wavelength to velocity and vice versa
        self.output.delete(0,tk.END)
        self.output.insert(tk.END,"Left Click (x) on a point to the left and to the right of what you wish to measure. Right Click (backspace) to remove a point.")
        clicks= plt.ginput(1)
        clicks= np.array(clicks)
        x=clicks[0][0]*u.AA
        y=clicks[0][1]
        self.wavelength=(self.wavelength-x)/x*c.to('km/s')
        print(x,c.to('km/s'))
        self.splot()


    def restore(self,event=None):
        self.flux=self.flux_orig
        self.splot()

    def pltregion(self,x,y,sym='-',c='black'):
        #takes two clicks
        self.ax.plot(x,y,sym,color=c)
        self.canvas.draw()

    def region(self):
        #uses two clicks to define a region for fitting or measuring.
        self.output.delete(0,tk.END)
        self.output.insert(tk.END,"Left Click (x) on a point to the left and to the right of what you wish to select. Right Click (backspace) to remove a point.")
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
        #collects many points.
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

    def eqw(self,event=None):
        # Measure equivalent width between two points IRAF style
        self.measuremode()
        x,y=self.region() #click points to slice data
        xg,yg=self.chop(self.wavelength,self.flux,x[0],x[1])

        continuum=(yg[0]+yg[-1])/2
        dwidth=[]
        for i,f in enumerate(yg):
          if i ==0:
            pass
          else:
            dwidth.append((1-f/continuum)*abs(xg[i]-xg[i-1]))
        width=sum(dwidth)
        bisect=(xg[0].value+xg[-1].value)/2.*xg[0].unit
        self.ax.vlines(bisect.value,min(yg),max(yg))
        self.canvas.draw()
        self.output.delete(0,tk.END)
        outstring="Equivalent Width = "+"{0.value:0.03f} {0.unit:FITS}".format(width)+\
                   ", Bisected Click Center = "+"{0.value:0.03f} {0.unit:FITS}".format(bisect)
        self.output.insert(tk.END,outstring)
        print(outstring)


    def find_nearest_index(self,array,value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx

    def chop(self,array1,array2,start1,stop1):
        small=min([start1,stop1])
        large=max([start1,stop1])
        startidx=self.find_nearest_index(array1,small)
        stopidx=self.find_nearest_index(array1,large)
        c1=array1[startidx:stopidx]
        c2=array2[startidx:stopidx]
        return c1,c2

    def scopy(self,event=None):
        self.measuremode()
        x,y=self.region() #click points to slice data
        xg,yg=self.chop(self.wavelength,self.flux,x[0],x[1])
        self.wavelength=xg
        self.flux=yg
        self.header['CRVAL1']=x[0]
        self.splot()

    def scopy_manual(self):
        #not working yet, the code doesn' wait for the dialog box entry.
        #turned off in the menu options so not accessible to the user.
        self.EntryDialog(message="Enter the left and right wavelengths to trim spectra, separated by a space.")
        print(self.response)
##        xg,yg=self.chop(self.wavelength,self.flux,x1,x2)
##        self.wavelength=xg
##        self.flux=yg
##        self.splot()

    def fit(self,func="gauss",event=None):
        # Fit the data using a Gaussian, Voigt, or Lorentzian profile
        self.measuremode()
        x,y=self.region() #click points to slice data
        xg,yg=self.chop(self.wavelength,self.flux,x[0],x[1])
        xgf=np.linspace(xg[0],xg[-1],len(xg)*3)
        ygf=np.interp(xgf,xg,yg)

        mid=len(yg)//2
        linecoeff = np.polyfit(xgf,ygf,1)
        ygf=ygf/np.polyval(linecoeff,xgf)
        invert=False
        reflevel=yg[0]
        if yg[mid] < reflevel:
            ygf=reflevel-ygf
            invert=True

        if func=="gauss":
          g_init = models.Gaussian1D(amplitude=np.max(yg)-reflevel, mean=np.mean(xg), stddev=xg[1]-xg[0])
          fit_g = fitting.LevMarLSQFitter()
          if invert == False:
            g = fit_g(g_init, xgf, ygf)
            plt.plot(xgf, g(xgf)*np.polyval(linecoeff,xgf), label='%s'%(func))
            amp=g.amplitude
          elif invert == True:
            g = fit_g(g_init, xgf, ygf)
            plt.plot(xgf, (reflevel-g(xgf))*np.polyval(linecoeff,xgf), label='%s'%(func))
            amp=-g.amplitude
          outstring="Gaussian Center: "+"{0.value:0.03f} {0.unit:FITS}".format(g.mean)+", FWHM: "+ \
                           "{0.value:0.03f} {0.unit:FITS}".format(g.fwhm)+", Amplitude: "+"{0.value:0.03f} {0.unit:FITS}".format(amp)
          print(outstring)
          self.output.delete(0,tk.END)
          self.output.insert(tk.END, outstring)
        elif func=="voigt":
          g_init = Voigt1D(x_0=np.mean(xg),amplitude_L=np.max(yg)-reflevel , fwhm_L=xg[1]-xg[0], fwhm_G=xg[1]-xg[0])
          fit_g = fitting.LevMarLSQFitter()
          if invert == False:
            g = fit_g(g_init, xgf, ygf)
            plt.plot(xgf, g(xgf)*np.polyval(linecoeff,xgf), label='%s'%(func))
            amp=g.amplitude_L
          elif invert == True:
            g = fit_g(g_init, xgf, ygf)
            plt.plot(xgf, (reflevel-g(xgf))*np.polyval(linecoeff,xgf), label='%s'%(func))
            amp=-g.amplitude_L
          outstring=  "Voigt Center: "+"{0.value:0.03f} {0.unit:FITS}".format(g.x_0)+", Lorentzian_FWHM: "+\
                  "{0.value:0.03f} {0.unit:FITS}".format(g.fwhm_L)+\
                  ", Gaussian_FWHM: "+"{0.value:0.03f} {0.unit:FITS}".format(g.fwhm_G)+\
                  ", Amplitude: "+"{0.value:0.03f} {0.unit:FITS}".format(amp)
          print(outstring)
          self.output.delete(0,tk.END)
          self.output.insert(tk.END,outstring)
        elif func=="lorentz":
          g_init = Lorentz1D(x_0=np.mean(xg),amplitude=np.max(yg)-reflevel, fwhm=xg[1]-xg[0])
          fit_g = fitting.LevMarLSQFitter()
          if invert == False:
            g = fit_g(g_init, xgf, ygf)
            plt.plot(xgf, g(xgf)*np.polyval(linecoeff,xgf), label='%s'%(func))
            amp=g.amplitude
          elif invert == True:
            g = fit_g(g_init, xgf, ygf)
            plt.plot(xgf, (reflevel-g(xgf))*np.polyval(linecoeff,xgf), label='%s'%(func))
            amp=-g.amplitude
          outstring= "Lorentz Center: "+"{0.value:0.03f} {0.unit:FITS}".format(g.x_0)+", FWHM: "+\
                  "{0.value:0.03f} {0.unit:FITS}".format(g.fwhm)+\
                  ", Amplitude: "+"{0.value:0.03f} {0.unit:FITS}".format(amp)
          print(outstring)
          self.output.delete(0,tk.END)
          self.output.insert(tk.END,outstring)


        else:
          print("Error with Fitting Function Selection")

        self.canvas.draw()


    def norm_clear(self,event=None):
        self.x_norm=[]
        self.y_norm=[]
        self.goodfit=False
        self.output.delete(0,tk.END)
        self.splot()

    def continuum(self,event=None):
        #Create array of click to use as the continuum.
        self.measuremode()
        x,y=self.region() #click points to slice data
        self.pltregion(x,y,sym='-s',c='black')
        xg,yg=self.chop(self.wavelength,self.flux,x[0],x[1])
        for xi in xg:
            self.x_norm.append(xi.value)
        for yi in yg:
            self.y_norm.append(yi)

    def normalize(self,event=None):
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
                self.ax.set_ylim([max(0,min(nflux)),min(100,max(nflux))])
                self.canvas.draw()
                answer=tkinter.messagebox.askyesno("Question","Proceed with the fit?")
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

    def save_fits(self):
        self.sp[0].header=self.header
        self.sp[0].data=self.flux
        savename=asksaveasfilename(defaultextension=".fits")
        self.sp.writeto(savename)

    def save1DText(self):
        savename=asksaveasfilename(defaultextension=".txt")
        dataout=open(savename,'w')
        for i,val in enumerate(self.flux):
            dataout.write('%s %s\n'%(self.wavelength[i].value,self.flux[i]))
        dataout.close()

    def _quit(self):
        self.master.destroy()  # this is necessary on Windows to prevent
                        # Fatal Python Error: PyEval_RestoreThread: NULL tstate
        self.master.quit()     # stops mainloop


    def imhead(self,event=None):
        t=tk.Toplevel(self.master,height=600,width=600)
        t.wm_title("FITS Header")
        s=tk.Scrollbar(t)
        s.pack(side=tk.RIGHT,fill=tk.Y)
        hlist=tk.Listbox(t,yscrollcommand=s.set,height=20,width=80)

        for keys in self.header.tostring(sep=',').split(','):
            hlist.insert(tk.END,"%s "%(keys))
        hlist.pack(side=tk.LEFT,fill=tk.BOTH)
        s.config(command=hlist.yview)
        b = tk.Button(t,text="Close", command=lambda: self.destroychild(t))
        b.pack()

    def EntryDialog(self,message="Text Input"):
        t=tk.Toplevel(self.master,height=200,width=200)
        t.wm_title("Text Input Needed")
        tk.Label(t,text=message).pack()
        wv=tk.StringVar()
        self.w = tk.Entry(t,text=wv,width=20)
        self.w.pack()
        b = tk.Button(t,text="Set Value", command=self.fetch)
        b.pack()
        b2 = tk.Button(t,text="Close", command=lambda: self.destroychild(t))
        b2.pack()

    def fetch(self):
        self.response=self.w.get()

    def destroychild(self,w):
        w.destroy()

    def About(self):
        t=tk.Toplevel(self.master,height=600,width=600)
        t.wm_title("About")
        tk.Label(t,text="Author: Dr. Joshua Thomas \n thomas.joshd@gmail.com").pack()
        tk.Label(t,text="This program was designed to emulate some basic IRAF splot functions.").pack()
        tk.Label(t,text="Uses Astropy library, and some parts are directly modified from the UVES tutorial.").pack()
        tk.Label(t,text="Version %s"%version).pack()
        tk.Label(t,text="Last Updated %s"%UPDATED).pack()
        b = tk.Button(t,text="Close", command=lambda: self.destroychild(t))
        b.pack()

#----------------------------------------------------------------------------------
#Begin GUI
root = tk.Tk() #main GUI window
program=App(root)
root.mainloop() #lets the GUI run
#root.destroy()

#Need to have measureing/fit/norm functions check to see if overplot or stackplot is engaged. Or it
#needs to replot the last spectrum and switch out of overplot/stack
#test some measurments of eqw, and gauss with iraf on same spectra.
#doppler correct spectra
#right now can't handle text spectra with headers
