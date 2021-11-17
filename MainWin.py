import os
import sys


from functools import partial

from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMessageBox,QInputDialog, QLineEdit
from PyQt5.QtCore import Qt, QSize

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.image as img
import matplotlib.cm as cm
import matplotlib.patches as patches
import matplotlib.colors as colors
from matplotlib import animation

import numpy as np
import csv


from astropy.wcs import WCS
from astropy.io import fits
from astropy.table import Table

import astropy.units as u
from astropy.constants import c

from astropy.modeling import models, fitting
from astropy.convolution import convolve, Box1DKernel, Gaussian1DKernel


import datetime

from ArithWin import ArithWin
from HeaderWin import HeaderWin
from version import VERSION,UPDATED
from plotdefaults import plotdefaults
from helperfunctions import find_nearest_index,julian
from Menu import Menu
from LogFiles import LogFiles
from Spectra import Spectra
from about import about
# from MainWidget import MainWidget

plotdefaults()

class MainWin(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        super(MainWin, self).__init__(parent)
        self.setAcceptDrops(True)
        self.setWindowTitle("PySplot, Version %s"%str(VERSION))
        self.initparams()
        self.message.append("Welcome to PySplot! Get started by opening a 1-D Spectrum or List of 1-D Spectra File>open or by drag and drop!")
        self.create_mainwidget()
        # MainWidget(self).create_mainwidget()
        self.resize(1200,800)
        self.log=LogFiles(self)
        self.Spectra=Spectra(self)
        Menu(self).menu()
        self.cmd()


    def cmd(self):
        try:
            self.Spectra.open(drop=sys.argv[1:])
        except:
            print("Command Line Argument Not Understood.")


    def initparams(self):
        self.listedfiles=[] #used for lists files
        self.slist=[]
        self.message=[]
        self.database={}
        self.loadedregions=False
        self.loadednorm=False
        self.loadedbisect=False
        self.loadedalign=False
        self.gridvalue=True
        self.titlevalue=''
        self.titletog=True
        self.legend=False
        self.overplot=False
        self.stackplot=False
        self.stack=[]
        self.stackint=0
        self.delimdef=False
        self.delimdef_question=False
        self.height=1
        self.order=1
        self.orderset=False
        self.ref_wave=1
        self.script=False
        self.stacknum=0
        self.limreset()
        self.poplast=False
        self.plotstyle='line'#'step','point'
        self.clickx=False
        self.clicky=False
        self.x=[]
        self.y=[]
        self.firstclick=False
        self.heightcheck=False
        self.boxwidth=False
        self.suffix=False
        self.color=True

    def dragEnterEvent(self, event):
        # print('drag-enter')
        if event.mimeData().hasUrls():
            # print('has urls')
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event):
        files = []
        for url in event.mimeData().urls():
            files.append(url.toLocalFile())
        # spec=['.fit','.fits','.FIT','.FITS','.txt','.TXT','.csv','.CSV','.dat','.DAT','.s']
        self.Spectra.open(drop=files)
        try:
            if '.par' in files[0]:
                self.paropen(files[0])
        except:
            pass

    def reset(self):
        """Just replots and clears the norm parameters.  A handy way to clear the graph of clutter."""
        self.ax.clear()
        self.limreset()
        self.replot()
        # self.region_clear()

    def region_clear(self):
        """clears region selection and associated overplots."""
        try:
            del self.x_norm
        except:
            pass
        try:
            del self.y_norm
        except:
            pass
        self.goodfit=False
        self.loadednorm=False
        self.orderset=False
        try:
            del self.saveregions_x
        except:
            pass
        try:
            del self.saveregions_y
        except:
            pass
        self.loadedregions=False
        self.loadedbisect=False
        self.loadedalign=False
        self.clickx=False
        self.clicky=False
        self.height=0
        self.x=[]
        self.y=[]
        self.firstclick=False
        self.heightcheck=False
        self.boxwidth=False
        self.replot()

    def replot(self):
        try:
            self.getlims()
            if self.overplot == False and self.stackplot == False:
                self.plotSpectra(spec=self.fname)
            else:
                # self.ax.clear()
                try:
                    self.plotSpectra(spec=self.slist)
                except:
                    self.plotSpectra()
        except:
            pass


    def create_mainwidget(self):
        self.create_pysplot()
        self.create_sidepane()
        self.create_output()
        self.mainwidget=QtWidgets.QWidget()
        self.setCentralWidget(self.mainwidget)

        # self.vbox = QtWidgets.QVBoxLayout()
        # self.vbox.addWidget(self.outputbox)
        # self.vbox.addWidget(self.pysplotbox)
        # self.hbox = QtWidgets.QHBoxLayout()
        # self.hbox.addLayout(self.vbox)
        # self.hbox.addWidget(self.sidepanebox)
        self.grid_layout=QtWidgets.QGridLayout()
        self.grid_layout.addWidget(self.outputbox,0,0)
        self.grid_layout.addWidget(self.pysplotbox,1,0)
        self.grid_layout.addWidget(self.sidepanebox,0,1,2,1)

        # self.grid_layout.setColumnStretch(1,3)
        # self.grid_layout.setRowStretch(1,4)

        self.grid_layout.setColumnStretch(0,4)
        self.grid_layout.setColumnStretch(1,1)
        self.grid_layout.setRowStretch(0,1)
        self.grid_layout.setRowStretch(1,7)


        self.grid_layout.setSpacing(0)
        self.grid_layout.setHorizontalSpacing(0)
        self.grid_layout.setVerticalSpacing(0)

        # self.grid_layout.setColumnStretch(0, 100)
        # self.grid_layout.setColumnStretch(1, 1)
        # self.grid_layout.setRowStretch(0, 1)
        # self.grid_layout.setRowStretch(1, 100)

        self.mainwidget.setLayout(self.grid_layout)




    def stackpane(self):
        self.sidepanebox.close()
        self.create_sidepane()
        self.grid_layout.addWidget(self.sidepanebox,0,1,2,1)
        self.mainwidget.setLayout(self.grid_layout)

    def outputupdate(self):
        try:
            self.outputbox.close()
        except:
            pass
        self.create_output()
        self.grid_layout.addWidget(self.outputbox,0,0)
        self.mainwidget.setLayout(self.grid_layout)


    def create_pysplot(self):
        self.pysplotbox = QtWidgets.QFrame()
        self.fig=Figure()
        self.ax = self.fig.add_subplot(111)
        self.ax.tick_params(right= True,top= True,which="both")
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar( self.canvas, self)

        layout = QtWidgets.QVBoxLayout()
        # layout.addWidget(button)
        layout.addWidget(self.canvas)
        layout.addWidget(self.toolbar)

        # self.ax.set_title("Single Spectra Mode",fontsize=12)
        # self.logo()
        self.frontispiece()

        self.toolbar.update()
        self.canvas.draw()
        # self.pysplotbox=QtWidgets.QWidget()
        self.pysplotbox.setLayout(layout)
        # self.setCentralWidget(widget)
    #
    # def logo(self):
          # '''Display logo on load.  Breaks the display.'''
    #     self.ax.imshow(img.imread('icon.png'))


    def frontispiece(self):
        self.a=self.fig.text(0.2,0.8,'^Message box^')
        self.b=self.fig.text(0.6,0.8,'Stack of spectra open-->')
        self.c=self.fig.text(0.2,0.2,'v Interactive plot controls v')
        self.d=self.fig.text(0.35,0.5,'All measruments saved to CSV file')

        base=.65
        shift=-.3
        lg=40
        sm=30
        self.e=self.fig.text(.2,base,"P",color='#7f4fc9',fontsize=lg)
        self.f=self.fig.text(.3,base,"y",color='#3e49bb',fontsize=sm)
        self.g=self.fig.text(.4,base+shift,"S",color='#526eff',fontsize=lg)
        self.h=self.fig.text(.5,base+shift,"p",color='#32c12c',fontsize=sm)
        self.i=self.fig.text(.6,base+shift,"l",color='#ffef00',fontsize=sm)
        self.j=self.fig.text(.7,base+shift,"o",color='#ff9a00',fontsize=sm)
        self.k=self.fig.text(.8,base+shift,"t",color='#d40c00',fontsize=sm)

    def frontispiece_clear(self):
        try:
            self.a.remove()
            self.b.remove()
            self.c.remove()
            self.d.remove()
            self.e.remove()
            self.f.remove()
            self.g.remove()
            self.h.remove()
            self.i.remove()
            self.j.remove()
            self.k.remove()
        except:
            pass

    def create_sidepane(self):
        self.sidepanebox=QtWidgets.QFrame()
        layout=QtWidgets.QVBoxLayout()
        scroll = QtWidgets.QScrollArea()             # Scroll Area which contains the widgets, set as the centralWidget
        widget = QtWidgets.QWidget()                 # Widget that contains the collection of Vertical Box
        self.vbox = QtWidgets.QVBoxLayout()           # The Vertical Box that contains the Horizontal Boxes of  labels and buttons
        # button1=QtWidgets.QPushButton("Delete Selected")
        # button1.clicked.connect(self.removefromstack)
        self.buildstack()
        # layout.addWidget(button1)
        widget.setLayout(self.vbox)

        #Scroll Area Properties
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        # scroll.setWidgetResizable(False)
        scroll.setWidget(widget)
        layout.addWidget(scroll)
        self.sidepanebox.setLayout(layout)

    def buildstack(self):
        try:
            # specnumber=list(np.arange(len(self.database))+1)
            # specnumber.reverse()
            # i=0
            specforward=list(self.database)
            speclist=specforward.copy()
            speclist.reverse()

            for f in speclist: #the syntax [::-1] reverses the list without modifying it so that  the button list is the same vertical order as the stack plotted spectra.
                b=QtWidgets.QPushButton("%s: %s"%(self.database[f]['stacknumber'],os.path.basename(f)))
                # b=QtQidgets.QPushButton("<font color=%s> %s: %s</font>"%(colors.to_hex(self.database[self.fname]['plotcolor']),t,os.path.basename(s)))
                # b.setText("<font color=%s> %s: %s</font>"%(colors.to_hex(self.database[self.fname]['plotcolor']),t,os.path.basename(s)))

                b.setStyleSheet("Text-align:left")
                # b.setStyleSheet("QPushButton::pressed"
                #                 "{"
                #                 "background-color : red;"
                #                 "}")
                b.clicked.connect(partial(self.stackwindowplot,f))
                # b=QtWidgets.QLabel()
                # b.setBackground(colors.to_hex(self.database[self.fname]['plotcolor']))#QColor(self.database[self.fname]['plotcolor']))
                # b.setColor(colors.to_hex(self.database[self.fname]['plotcolor']))
                # b.setStyleSheet("background-color: white")
                # b.setTextInteractionFlags(Qt.TextSelectableByMouse)
                self.vbox.addWidget(b)
        except:
            print('Exception occured in buildstack')

    def create_output(self):
        self.outputbox=QtWidgets.QWidget()
        layout=QtWidgets.QVBoxLayout()
        scroll = QtWidgets.QScrollArea()             # Scroll Area which contains the widgets, set as the centralWidget
        widget = QtWidgets.QWidget()                 # Widget that contains the collection of Vertical Box
        vbox = QtWidgets.QVBoxLayout()           # The Vertical Box that contains the Horizontal Boxes of  labels and buttons
        for m in self.message[::-1]: #the minus 1 puts the list inverted because the scrollbox always starts at the top....
            b=QtWidgets.QLabel(m)
            # b.setStyleSheet("background-color: white")
            b.setTextInteractionFlags(Qt.TextSelectableByMouse)
            vbox.addWidget(b)
        widget.setLayout(vbox)

        #Scroll Area Properties
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)

        scroll.setWidgetResizable(True)
        scroll.setWidget(widget)
        layout.addWidget(scroll)
        self.outputbox.setLayout(layout)

    def imhead(self):
        try:
            self.winpopup = HeaderWin(self)
            self.winpopup.show()
        except:
            self.message.append("No header to display.")
            self.outputupdate()

    def imarith(self):
        """Spectrum arithmatic window"""
        try:
            self.winpopup = ArithWin(self)
            self.winpopup.show()
        except:
            print("Unable to load Arithmetic Window.")

    def close_out(self):
        """what to do when closed."""
        self.log.endoflog()
        QtWidgets.QApplication.instance().quit()

    def About(self):
        about()


    def measuremode(self):
        """Puts us into single spectrum mode."""
        self.log.checklog()
        self.stackplot=False
        self.overplot=False

    def plotting_guts(self,sym='-',basename=''):
        if self.plotstyle == 'step':
            if self.stackplot == True:
                self.spec,=self.ax.step(self.wavelength,self.flux.value+(self.database[self.fname]['stacknumber']-1),sym,color=self.database[self.fname]['plotcolor'],label="%s: %s"%(self.database[self.fname]['stacknumber'],basename))
            else:
                if self.color == True:
                    self.spec,=self.ax.step(self.wavelength,self.flux,color=self.database[self.fname]['plotcolor'],\
                    label="%s: %s"%(self.database[self.fname]['stacknumber'],basename))
                else:
                    spec,=self.ax.step(self.wavelength,self.flux,color="black",\
                    label="%s: %s"%(self.database[self.fname]['stacknumber'],basename))
        else:
            if self.stackplot == True:
                self.spec,=self.ax.plot(self.wavelength,self.flux.value+(self.database[self.fname]['stacknumber']-1),sym,color=self.database[self.fname]['plotcolor'],label="%s: %s"%(self.database[self.fname]['stacknumber'],basename))
            else:
                if self.color == True:
                    self.spec,=self.ax.plot(self.wavelength,self.flux,sym,color=self.database[self.fname]['plotcolor'],\
                    label="%s: %s"%(self.database[self.fname]['stacknumber'],basename))
                else:
                    self.spec,=self.ax.plot(self.wavelength,self.flux,sym,color="black",\
                    label="%s: %s"%(self.database[self.fname]['stacknumber'],basename))




    def splot(self):
        """Sets up the plot"""
        # print(self.plotstyle)
        basename=os.path.basename(self.fname)
        #Normal plot mode
        if self.plotstyle=='point':
            sym='.-'
        elif self.plotstyle=='line':
            sym='-'
        else:
            sym='-'
        if self.overplot == False and self.stackplot == False:
            self.ax.clear()
            self.titlevalue="%s:%s"%(self.database[self.fname]['stacknumber'],basename)
            self.plotting_guts(sym,basename)
            self.ax.set_ylabel("Flux")

        #Overplot Mode
        elif self.overplot == True and self.stackplot == False:
            self.titlevalue="Overplot Mode, Display For Qualitative Comparison Only"
            self.plotting_guts(sym,basename)
            self.ax.set_ylabel("Flux")

        #Stackplot Mode
        elif self.overplot == False and self.stackplot == True:
            # if self.plotstyle == 'step':
            #     spec,=self.ax.step(self.wavelength,self.flux+(self.database[self.fname]['stacknumber']-1)*u.flx,\
            #     color=self.database[self.fname]['plotcolor'],label="%s: %s"%(self.database[self.fname]['stacknumber'],basename))
            # else:
            self.titlevalue="Stack Plot Mode, Display For Qualitative Comparison Only"
            self.plotting_guts(sym,basename)
            # spec,=self.ax.plot(self.wavelength,self.flux+(self.database[self.fname]['stacknumber']-1)*u.flx,sym,color=self.database[self.fname]['plotcolor'],label="%s: %s"%(self.database[self.fname]['stacknumber'],basename))
            self.ax.set_ylabel("Stack Plot Flux")

        if self.titletog == True:
            self.fig.suptitle("PySplot - Date: "+'{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now(),fontsize=10))
            self.ax.set_title(self.titlevalue)#,fontsize=10)
        else:
            pass


    def xaxislabel(self):
        """Determines the x-axis label."""
        if self.wavelength.unit == 'km / s':
            self.ax.set_xlabel("Velocity (%s)"%self.wavelength.unit)
        else:
            if 'heliocentric' in self.database[self.fname]:
                    self.ax.set_xlabel("Heliocentric Wavelength (%s)"%self.wavelength.unit)
            else:
                self.ax.set_xlabel("Wavelength (%s)"%self.wavelength.unit)

    def animated(self):
        anim = animation.FuncAnimation(self.fig, self.plotSpectra(), init_func=self.plotSpectra(),frames=100,interval=20, blit=True)



    def grabSpectra(self,spec):
        """Grab spectra from database, separate from plotting"""
        self.fname=spec
        self.wavelength=self.database[self.fname]['wavelength']
        self.flux=self.database[self.fname]['flux']
        self.flux_orig=self.database[self.fname]['flux_orig']
        try:
            self.header=self.database[self.fname]['header']
        except:
            print('Exception occured getting header')



    def plotSpectra(self,spec=None,savename=False):
        """Plot a spectrum or spectra."""
        self.frontispiece_clear()
        try:
            if spec == None :
                specforward=list(self.database)
                speclist=specforward.copy()
                speclist.reverse() #done so everything is in the same vertical order
                for f in speclist:
                    self.grabSpectra(f)
                    self.splot()

                    if self.legend == True:
                        self.ax.legend()
            elif spec != None : #plot singular spectrum
                if type(spec) is str:
                    self.grabSpectra(spec)
                    self.splot()
                elif type(spec) is list:
                    for val in spec:
                        self.grabSpectra(val)
                        self.splot()
                else:
                    pass
            self.plotRegions()


            if self.xlim[0] != 0:
                #setting previously set limits
                self.ax.set_xlim(self.xlim)
                self.ax.set_ylim(self.ylim)

            self.ax.grid(self.gridvalue)
            self.xaxislabel()
            if savename != False:
                self.canvas.print_figure(savename)
            if self.xlim[0] == 0:
                # valmax=max(self.flux)
                # if valmax > 2*average(self.flux):
                #     valmax=average(self.flux)
                # self.ax.set_ylim((min(min(self.flux),0),1000))
                self.toolbar.update()
            if self.legend == True:
                self.ax.legend()

            self.canvas.draw()
            self.stackpane()
        except:
            print('Exception occured in plotSpectra')

    def gridtoggle(self):
        """Toggles the plot grid on and off."""
        try:
            if self.gridvalue == False:
                self.gridvalue=True
            elif self.gridvalue == True:
                self.gridvalue=False
            else:
                pass
            self.ax.grid(self.gridvalue)
            self.canvas.draw()
        except:
            print('Exception occured in gridtoggle')

    def titletoggle(self):
        """Toggles the plot title on and off."""
        try:
            if self.titletog == False:
                self.ax.set_title(self.titlevalue)
                self.fig.suptitle("PySplot - Date: "+'{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now(),fontsize=10))
                self.titletog=True
            elif self.titletog == True:
                self.ax.set_title('')
                self.fig.suptitle('')
                self.titletog=False
            self.canvas.draw()
        except:
            print('Exception occured in titletoggle')

    def legendtoggle(self):
        """Toggles the plot legend on and off."""
        try:
            if self.legend == False:
                self.legend=True
                self.ax.legend()
                self.canvas.draw()
            elif self.legend == True:
                self.legend=False
                self.getlims()
                self.ax.clear()
                if self.overplot == False and self.stackplot == False:
                    self.plotSpectra(spec=self.fname)
                else:
                    self.plotSpectra()
            else:
                pass
        except:
            print('Exception occured in legendtoggle')

    def colortoggle(self):
        """Toggles the plot colors on and off."""
        try:
            if self.color == False:
                self.color=True
                self.plotSpectra()
            elif self.color == True:
                self.color=False
                self.getlims()
                self.ax.clear()
                if self.overplot == False and self.stackplot == False:
                    self.plotSpectra(spec=self.fname)
                else:
                    self.plotSpectra()
            else:
                pass
        except:
            print('Exception occured in colortoggle')



    def set_style_line(self):
        """Toggles the plot style."""
        try:
            self.ax.clear()
            self.plotstyle='line'
            self.plotSpectra(spec=self.fname)
        except:
            print("Failed to change line style.")

    def set_style_point(self):
        """Toggles the plot style."""
        try:
            self.ax.clear()
            self.plotstyle='point'
            self.plotSpectra(spec=self.fname)
        except:
            print("Failed to change line style.")

    def set_style_step(self):
        """Toggles the plot style."""
        try:
            self.ax.clear()
            self.plotstyle='step'
            self.plotSpectra(spec=self.fname)
        except:
            print("Failed to change line style.")



    def singleplottoggle(self):
        """Single spectrum display mode, replots the last active spectrum."""
        self.getlims()
        self.ax.clear()
        self.overplot=False
        self.stackplot=False
        try:
            self.plotSpectra(spec=self.fname)
        except:
            pass

    def overplottoggle(self):
        """Switches to overplot mode and replots the stack."""
        try:
            self.overplot=True
            self.stackplot=False
            self.getlims()
            self.ax.clear()
            self.plotSpectra()
        except:
            pass


    def stackplottoggle(self):
        """Switches to stackplot mode and replots the stack."""
        try:
            self.stackplot=True
            self.overplot=False
            self.reset()
            self.plotSpectra()
        except:
            pass

    def limreset(self):
        self.xlim=(0,1)
        self.ylim=(0,1)


    def getlims(self):
        """gets current display limits"""
        try:
            self.xlim=self.ax.get_xlim()
            self.ylim=self.ax.get_ylim()
        except:
            print('Exception at getlims')


    def stackreset(self):
        """Reset the stacking parameters"""
        self.initparams()
        try:
            del self.wavelength
        except:
            pass
        try:
            del self.flux
        except:
            pass
        try:
            del self.flux_orig
        except:
            pass
        try:
            del self.header
        except:
            pass
        self.region_clear()
        self.ax.clear()
        self.log.endoflog()
        self.stackpane()
        self.message.append("All spectra closed.")
        self.outputupdate()
        self.canvas.draw()

    def currentspectra(self):
        """prints current spectrum"""
        try:
            self.message.append("%s"%self.database[self.fname].keys())
            self.outputupdate()
        except:
            print('Exception at currentspectra')

    def rname(self):
        """renames the current spectrum"""
        try:
            Newname, okPressed = QInputDialog.getText(self, "Set new name for spectrum","New name: ", QLineEdit.Normal, "%s"%str(self.fname))
            if okPressed:
                self.database[Newname]=self.database[self.fname]
                del self.database[self.fname]
                self.fname=Newname
                self.Spectra.updatespectrum()
                self.stackrebuild()
                self.ax.clear()
                self.plotSpectra(spec=self.fname)
            else:
                pass
        except:
            print('Excpetion at rname')

    def stackprint(self):
        """Prints stack info."""
        try:
            self.message.append("Stack:")
            for f in self.database:
                self.message.append("%s"%f)
            self.message.append("Number in Stack: %s"%(str(len(self.database))))
            self.message.append("End Stack.")
            self.outputupdate()
        except:
            print('Exception at stackprint')


    def savestack(self,name="stack.list",stack=False):
        """Save the stack for later re-loading."""
        try:
            if stack == False:
                stack=self.stack
            path=os.path.dirname(self.fname)
            # basename=os.path.basename(self.fname)
            path_wo_ext=os.path.splitext(self.fname)[0]
            initialname=os.path.join(path_wo_ext,"stack.list")
            savename, _ = QtWidgets.QFileDialog.getSaveFileName(self,"Save Stack List",initialname)
            dataout=open(savename,'w')
            for row in stack:
                dataout.write('%s\n'%(row))
            dataout.close()
            self.message.append("Stack saved to:  %s"%(savename))
            self.outputupdate()
        except:
            self.message.append("Stack is empty.")
            self.outputupdate()


    def click_1_text(self):
        try:
            self.message.append("Click on a point to display the coordinates.")
            self.outputupdate()
            self.region_clear()
        except:
            print('Exception occred in click_1_text')

    def coord(self):
        try:
            self.measuremode()
            self.click_1_text()
            self.click=self.fig.canvas.mpl_connect('button_press_event', self.mouseclick_coord)
        except:
            print('Exception occred in coord')

    def mouseclick_coord(self,event):
        try:
            self.clickx,self.clicky=event.xdata,event.ydata
            self.coord_plot()
            self.fig.canvas.mpl_disconnect(self.click)
        except:
            print('Exception occred in mouseclick_coord')

    def coord_plot(self):
        """uses one click to print mouse position in the text box, and log.  Only
        cares about the x-coordinate"""
        try:
            idx=find_nearest_index(self.wavelength,self.clickx)
            xg=self.wavelength[idx]
            yg=self.flux[idx]
            self.ax.plot(xg,yg,'ks')
            t=self.filedate()
            #this doesn't work correctly.  Check that its using data or mouse position.
            self.message.append(t+"Mouse Position Wavelength,"+"{0.value:0.03f}, {0.unit:FITS}".format(xg)+\
            ", Flux, "+"{0.value:0.03f}, {0.unit:FITS}".format(yg))
            self.outputupdate()
            # print(self.message)
            self.log.checklog()
            self.log.write(self.message[-1])
            #self.log.write('\n')
            self.canvas.draw()
        except:
            print('Exception occred in coord_plot')

    def velocity(self):
        try:
            self.message.append("Click on where you want zero to be.")
            self.outputupdate()
            self.click=self.fig.canvas.mpl_connect('button_press_event', self.mouseclick_vel)
        except:
            print('Exception occred in  velocity')

    def mouseclick_vel(self,event):
        try:
            self.clickx,self.clicky=event.xdata,event.ydata
            self.vel_plot()
            self.fig.canvas.mpl_disconnect(self.click)
        except:
            print('Exception occred in mouseclick_vel')

    def vel_plot(self):
        """uses one click to convert wavelength to velocity and vice versa"""
        try:
            # self.measuremode()
            if self.wavelength.unit == u.AA:
                # self.message.append("Left Click on where you want zero to be.")
                # self.outputupdate()
                # clicks= plt.ginput(1)
                # clicks= np.array(clicks)
                # print(self.clickx)
                self.wavecenter=self.clickx*u.AA
                self.wavelength_backup=self.wavelength
                self.wavelength=(self.wavelength-self.wavecenter)/self.wavecenter*c.to('km/s')
                # self.splot()
                self.Spectra.updatespectrum()
                self.plotSpectra()
                self.reset()
                # self.plotRegions()
                # self.canvas.draw()
                # self.stackpane()
            elif self.wavelength.unit == 'km / s':
                self.wavelength=self.wavelength_backup
                self.Spectra.updatespectrum()
                self.plotSpectra()
                self.reset()
                # self.plotRegions()
                # self.canvas.draw()
                # self.stackpane()
        except:
            print('Exception occred in vel_plot')





    def pltregion(self,x,y,sym='-',c='black'):
        """takes two clicks, and plots a line between them"""
        try:
            self.ax.plot(x,y,sym,color=c)
            self.canvas.draw()
        except:
            print('Exception occred in pltregion')


    def click_height(self):
        try:
            self.measuremode()
            self.message.append("Click to set vertical height.")
            self.outputupdate()
            self.click=self.fig.canvas.mpl_connect('button_press_event', self.mouseclick_height)
        except:
            print('Exception occred in click_height')


    def bisectheight(self):
        try:
            stop=False
            # self.height, okPressed = QInputDialog.getDouble(self, "Float","Height To Measure Bisection: (Note, can use shortcut '1' (see region menu) to select a height before running task.", self.height, -10000, 10000, 30)
            height, okPressed = QInputDialog.getText(self, "Set Bisect Height(s)","Height(s) To Measure Bisection comma separated.", QLineEdit.Normal, "%s"%str(self.height))
            self.height=[]#float(height)
            hs=list(height.split(","))
            for h in hs:
                self.height.append(float(h))
            if okPressed:
                pass
            else:
                stop=True #if cancelled don't try to run the rest of the bisect.
            return stop
        except:
            print('Exception occred in bisectheight')
            return True

    def BisectLine(self):
        """Designed for measuring the center of very broad Wolf-Rayet star emission lines, may work for other situtaitons."""
        self.measuremode()
        self.regionload()
        xg,yg=self.chopclick()
        reflevel,invert=self.invert_check(xg,yg)
        stop = False
        try:
            while stop == False:
                if self.loadedbisect == False and self.heightcheck == False:
                    if self.script == False:
                        stop = self.bisectheight()
                        # print('no script')
                    elif self.script == True and self.stacknum == 0:
                        stop = self.bisectheight()
                        # print('script')
                elif self.heightcheck == True:
                    self.heightcheck = False

                #split spectrum at peak
                if invert == False:
                    split=find_nearest_index(yg.value,np.max(yg.value))
                elif invert == True:
                    split=find_nearest_index(yg.value,np.min(yg.value))
                else:
                    print('Invert Problem in Bisect')
                    stop=True
                yg_left=yg[0:split+1]
                yg_right=yg[split::]
                xg_left=xg[0:split+1]
                xg_right=xg[split::]
                c1=0
                c2=0
                # print(self.height,type(self.height))
                for h in self.height:
                    if invert == False:
                        for i,val in enumerate(yg_left):
                            if yg_left[i].value < h and yg_left[i+1].value > h:
                                idx=i
                        x1=xg_left[idx:idx+2].value
                        y1=yg_left[idx:idx+2].value

                        for i,val in enumerate(yg_right):
                            if yg_right[i].value > h and yg_right[i+1].value < h:
                                idx=i
                        x2=xg_right[idx:idx+2].value
                        y2=yg_right[idx:idx+2].value

                    elif invert == True:
                        for i,val in enumerate(yg_left):
                            if yg_left[i].value > h and yg_left[i+1].value < h:
                                idx=i
                        x1=xg_left[idx:idx+2].value
                        y1=yg_left[idx:idx+2].value

                        # linecoeff = np.polyfit(x,y,1) #does a linear polynomail fit to the data.
                        # c1=(h-linecoeff[1])/linecoeff[0] #x=slope*y+intercept, swap x,y becauase here we want y to be the indep. variable.
                        # # self.ax.plot(x,y,'s')
                        # # self.ax.plot(c1,h,'^')

                        for i,val in enumerate(yg_right):
                            if yg_right[i].value < h and yg_right[i+1].value > h:
                                idx=i
                        x2=xg_right[idx:idx+2].value
                        y2=yg_right[idx:idx+2].value
                        # linecoeff = np.polyfit(x,y,1) #does a linear polynomail fit to the data.
                        # c2=(h-linecoeff[1])/linecoeff[0]  #x=slope*y+intercept, swap x,y becauase here we want y to be the indep. variable.
                        # self.ax.plot(x,y,'s')
                        # self.ax.plot(c2,h,'^')
                    else:
                        print('Invert Problem in Bisects')

                    linecoeff = np.polyfit(x1,y1,1) #does a linear polynomail fit to the data.
                    c1=(h-linecoeff[1])/linecoeff[0] #x=slope*y+intercept, swap x,y becauase here we want y to be the indep. variable.

                    linecoeff = np.polyfit(x2,y2,1) #does a linear polynomail fit to the data.
                    c2=(h-linecoeff[1])/linecoeff[0]  #x=slope*y+intercept, swap x,y becauase here we want y to be the indep. variable.

                    center=(c1+c2)/2.
                    if self.script == False:
                        # self.ax.plot(xg_right,yg_right,'bo')
                        # self.ax.plot(xg_left,yg_left,'r^')
                        self.ax.plot(x1,y1,'s')
                        self.ax.plot(c1,h,'^')
                        self.ax.plot(x2,y2,'s')
                        self.ax.plot(c2,h,'^')
                        self.ax.hlines(h,min(self.wavelength.value),max(self.wavelength.value)) #plots the cross-cut height
                        self.ax.vlines(center,min(self.flux.value),max(self.flux.value)) #plot line at bisect location

                    t=self.filedate()
                    self.message.append(t+"Bisected Center, %s"%(float(center))+", {0.unit:FITS}".format(self.wavelength)+\
                               ", Height, %s"%(h)+", {0.unit:FITS}".format(self.flux))

                    if self.script == False:
                        self.outputupdate()
                        self.canvas.draw()
                    self.log.checklog()
                    self.log.write(self.message[-1])
                    stop=True
        except:
            print('Exception occured in BisectLine')

    def mouseclick_height(self,event):
        try:
            self.height=event.ydata
            self.message.append("Height Set: %s"%(self.height))
            self.heightcheck=True
            self.outputupdate()
            self.fig.canvas.mpl_disconnect(self.click)
            self.firstclick=False
        except:
            print('Exception occred in mouseclick_height')



    def mouseclick_region(self,event):
        try:
            if self.firstclick is False:
                self.x.append(event.xdata)
                self.y.append(event.ydata)
                self.firstclick=True
            else:
                self.x.append(event.xdata)
                self.y.append(event.ydata)
                self.fig.canvas.mpl_disconnect(self.click)
                self.plotRegions()
                self.loadedregions=True
                self.firstclick=False
        except:
            print('Exception occred in mouseclick_region')


    def region(self):
        """uses two clicks to define a region for fitting or measuring."""
        try:
            #self.measuremode()
            # if message == ' ':
            regionmessage="Define region first: Click on left and right edges of your region.  Then run your function."
            # else:
            #     regionmessage=message
            self.message.append(regionmessage)
            self.outputupdate()
            self.click=self.fig.canvas.mpl_connect('button_press_event', self.mouseclick_region)
        except:
            print('Exception occred in region')

    def chopclick(self):
        return self.chop(self.wavelength,self.flux,self.x[0],self.x[1])

    def regionload(self):
        """Loads saved regions and gets the cloest values in the data"""

        if self.loadedregions == False :
            self.region() #click points to slice data
            self.saveregions_x=self.x
            self.saveregions_y=self.y
            self.loadedregions=True
        else:
            self.x,self.y=(self.saveregions_x,self.saveregions_y)
        # return self.chopclick()
        # except:
        #     # print('Exception occred in regionload')
        #     return None,None

    def signal2noise(self):
        """Measure the Signal To Noise in a region."""
        try:
            self.measuremode()
            # xg,yg=self.regionload()
            self.regionload()
            xg,yg=self.chopclick()
            mu=np.average(yg)
            sigma=np.std(yg)
            snr=mu/sigma
            t=self.filedate()
            self.message.append(t+"Signal to Noise Ratio, "+"{0.value:0.03f}".format(snr)+\
            ",Average, {0.value:0.03f}".format(mu)+",STD, {0.value:0.03f}".format(sigma))

            self.outputupdate()
            self.log.checklog()
            self.log.write(self.message[-1])
            #self.log.write('\n')
        except:
            print('Exception occred in signal2noise')

    def eqw(self):
        """Measure equivalent width between two points IRAF style"""
        try:
            # self.log.checklog()
            self.measuremode()
            # xg,yg=self.regionload()
            self.regionload()
            xg,yg=self.chopclick()

            continuum=(yg[0]+yg[-1])/2. #linearly normalize the feature from the horizontal click locations.
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
            t=self.filedate()
            self.message.append(t+"Equivalent Width, "+"{0.value:0.03f}, {0.unit:FITS}".format(width)+\
                       ", Bisected Click Center, "+"{0.value:0.03f}, {0.unit:FITS}".format(bisect))
            self.outputupdate()
            self.log.checklog()
            self.log.write(self.message[-1])
            # self.log.write('\n')
        except:
            print('Exception Occured in EQW')

    def filedate(self):
        basename=os.path.basename(self.fname)
        d=self.database[self.fname]
        test=False
        try:
            jd=str(d['header']['JD'])
            # print(jd)
        except:
            jd=''
            test=True
        try:
            obsdate=str(d['header']['DATE-OBS'])
            # print(obsdate)
        except:
            obsdate=''

        try:
            if test == True and obsdate != '':
                y=int(obsdate[0:4])
                m=int(obsdate[5:7])
                d=int(obsdate[8:10])
                h=int(obsdate[11:13])
                mm=int(obsdate[14:16])
                s=int(obsdate[17:19])
                date=(y,m,d,h,mm,s)

                jd=str(julian(date))
        except:
            pass

        try:
            hjd=str(d['header']['HJD'])
            # print(hjd)
        except:
            hjd=''
        return "%s,JD, %s, HJD, %s, OBS-DATE, %s,"%(basename,jd,hjd,obsdate)


    def chop(self,array1,array2,start1,stop1):
        """Cuts out a region of the spectrum for fitting or other functions
        that should only be run on a section of the spectrum."""
        try:
            small=min([start1,stop1])
            large=max([start1,stop1])
            startidx=find_nearest_index(array1,small)
            stopidx=find_nearest_index(array1,large)
            c1=array1[startidx:stopidx]
            c2=array2[startidx:stopidx]
            return c1,c2
        except:
            print('Exception occred in chop')
            return None,None

    def clickalign(self):
        try:
            # self.measuremode()
            self.message.append("Click on a point to select reference piont.")
            self.outputupdate()
            self.click=self.fig.canvas.mpl_connect('button_press_event', self.mouseclick_align)
        except:
            print('Exception occured in clickalign')

    def mouseclick_align(self,event):
        try:
            self.clickx,self.clicky=event.xdata,event.ydata
            self.align_1()
            self.fig.canvas.mpl_disconnect(self.click)
        except:
            print('Exception occured in mouseclick_align')

    def align_1(self):
        """Align a spectrum to a feature.
        Aligns off one click instead of gaussian, then alignes to reference wavelength."""
        # ref_wave=5889.950
        try:
            if self.loadedalign == False:
                if self.script == False:
                    check=self.getalign()
                elif self.script == True and self.stacknum == 0:
                    check=self.getalign()
            # print(lines)
            #shift the spectrum and plot
            if check == True:
                shift=self.ref_wave-self.clickx
                self.wavelength=self.wavelength+shift*self.wavelength[0]/self.wavelength[0].value
                self.Spectra.updatespectrum()
                if self.script == False:
                    self.plotSpectra()
                    self.reset()
        except:
            print('Exception occured in align_1')

    def getalign(self):
        try:
            self.ref_wave, okPressed = QInputDialog.getDouble(self, "Float","Wavelength to align feature to: ", self.ref_wave, 0, 100000, 4)
            if okPressed:
                return True
            else:
                self.fig.canvas.mpl_disconnect(self.click)
                return False
        except:
            print('Exception occured in getalign')

    def align(self):
        """Align a spectrum to a feature. Based off sodium_shifted_norm_v3.py
        Fits a gaussian to a feature, then alignes to reference wavelength."""
        try:
            #self.measuremode()
            self.regionload()
            xg,yg=self.chopclick()
            if self.loadedalign == False:
                if self.script == False:
                    self.getalign()
                elif self.script == True and self.stacknum == 0:
                    self.getalign()

            min_peakind=find_nearest_index(xg.value,min(xg.value))
            lines=[]
            g_init = models.Gaussian1D(amplitude=1./min(yg), mean=xg[min_peakind], stddev=(xg[1]-xg[0])*4) #fits the peak

            # print(1/min(yg),xg[min_peakind])
            fit_g = fitting.LevMarLSQFitter()
            # fit_g = fitting.SLSQPLSQFitter()
            # fit_g = fitting.SimplexLSQFitter()
            # fit_g = fitting.LinearLSQFitter()

            x=xg
            y=1/yg-1/yg[0]

            g = fit_g(g_init, x, y)
            if self.script==False:
                self.ax.plot(x,y,drawstyle='steps-mid')
                self.ax.plot(x, g(x), label='Gaussian')
            lines.append(g.mean.value)

            # print(lines)
            #shift the spectrum and plot
            shift=self.ref_wave-lines[0] #sodium D1
            self.wavelength=self.wavelength+shift*self.wavelength[0]/self.wavelength[0].value
            self.Spectra.updatespectrum()
            if self.script == False:
                self.plotSpectra()
                self.reset()
        except:
            print('Exception occured in align')

    def linrescale(self):
        """Strech the wavelength scale with a lienear function"""
        try:
            self.grabSpectra(self.fname)
            stop=False
            while stop == False:
                slope, okPressed = QInputDialog.getDouble(self, "Float","Slope of linear correction function: ", 0, -1000, 1000, 7)
                if okPressed:
                    pass
                else:
                    stop = True
                    break
                intercept, okPressed = QInputDialog.getDouble(self, "Float","Intercept linear correction function: ", 0, -1000, 1000, 7)
                if okPressed:
                    pass
                else:
                    stop = True
                    break
                newwave=self.wavelength*slope+intercept*self.wavelength.unit
                self.wavelength=newwave
                self.Spectra.updatespectrum()
                self.plotSpectra(spec=self.fname)
                stop=True
        except:
            print('Exception occured in linrescale')

    def scopy(self,script=False):
        """Copy out a section of a spectrum to a new spectrum"""
        try:
            #self.measuremode()
            self.regionload()
            xg,yg=self.chopclick()
            self.wavelength=xg
            self.flux=yg
            # self.xlim=(0,1)
            # self.ylim=(0,1)
            self.Spectra.updatespectrum()
            self.plotSpectra(spec=self.fname)
        except:
            print('Exception occured in scopy (crop)')

    def invert_check(self,xg,yg):
        try:
            # self.regionload()
            # xg,yg=self.chopclick()
            xgf=np.linspace(xg[0],xg[-1],len(xg)*3)
            ygf=np.interp(xgf,xg,yg)#*u.flx
            # mid=len(yg)//2 #get guess for line peak by taking center of clicks
            xa=np.array([xgf[0].value,xgf[-1].value])
            try:
                ya=np.array([ygf[0].value,ygf[-1].value])
            except:
                ya=np.array([ygf[0],ygf[-1]])
            linecoeff = np.polyfit(xa,ya,1) #does a linear polynomail fit to the data.
            ygfn=ygf/np.polyval(linecoeff,xgf.value)#*u.flx #normalized
            ygn=yg/np.polyval(linecoeff,xg/xg[0])*u.flx #xg/xg[0] to remove the unit
            invert=False
            reflevel=yg[0]
            if yg[int(len(yg)//2)] < reflevel:
                ygf=reflevel.value-ygf.value
                invert=True

            return (reflevel,invert)
        except:
            print('Exception occured in invert_check')
            return (None,None)


    def fit(self,func="gauss"):
        """wrapper function for fitting line profiles"""
        # Fit the data using a Gaussian, Voigt, or Lorentzian profile
        try:
            self.measuremode()
            self.regionload()
            xg,yg=self.chopclick()
            reflevel,invert=self.invert_check(xg,yg)

            xgf=np.linspace(xg[0],xg[-1],len(xg)*3)
            ygf=np.interp(xgf,xg,yg)#*u.flx

            if func=="gauss":
              g_init = models.Gaussian1D(amplitude=np.max(yg), mean=np.mean(xg), stddev=xg[2]-xg[0])
              fit_g = fitting.LevMarLSQFitter()
              if invert == False:
                g = fit_g(g_init, xg, yg)
                if self.script == False:
                    self.ax.plot(xgf, g(xgf), label='%s'%(func))#*np.polyval(linecoeff,xgf*xg[0]/xg[0].value)
                amp=g.amplitude
                # self.ax.plot(xgf,ygfn)
                # self.ax.plot(xgf,g(xgf))
              elif invert == True:
                g = fit_g(g_init, xg, reflevel-yg)
                if self.script == False:
                    self.ax.plot(xgf, (reflevel-g(xgf)), label='%s'%(func))
                amp=reflevel-g.amplitude
                # self.ax.plot(xgf,-ygfn)
                # self.ax.plot(xgf,-g(xgf))
              t=self.filedate()
              self.message.append(t+"Gaussian Center, "+"{0.value:0.03f}, {0.unit:FITS}".format(g.mean)+", FWHM, "+ \
                               "{0.value:0.03f}, {0.unit:FITS}".format(g.fwhm)+", Amplitude, "+"{0.value:0.03f}, {0.unit:FITS}".format(amp))
              # print(self.message)
              self.log.checklog()
              self.log.write(self.message[-1])
              #self.log.write('\n')
              #
              self.outputupdate()


            elif func=="voigt":
              g_init = models.Voigt1D(x_0=np.mean(xg),amplitude_L=np.max(yg)-reflevel , fwhm_L=xg[1]-xg[0], fwhm_G=xg[1]-xg[0])
              fit_g = fitting.LevMarLSQFitter()
              if invert == False:
                g = fit_g(g_init, xg, yg)
                if self.script == False:
                    self.ax.plot(xgf, g(xgf))#*np.polyval(linecoeff,xgf), label='%s'%(func))
                amp=g.amplitude_L
                # self.ax.plot(xg,ygn)
                # self.ax.plot(xgf,g(xgf))
              elif invert == True:
                g = fit_g(g_init, xg, reflevel-yg)
                if self.script == False:
                    self.ax.plot(xgf, (reflevel-g(xgf)))#*np.polyval(linecoeff,xgf), label='%s'%(func))
                # self.ax.plot(xg,-ygn)
                # self.ax.plot(xgf,-g(xgf))
                amp=reflevel-g.amplitude_L
              t=self.filedate()
              self.message.append(t+"Voigt Center, "+"{0.value:0.03f}, {0.unit:FITS}".format(g.x_0)+", Lorentzian_FWHM, "+\
                      "{0.value:0.03f}, {0.unit:FITS}".format(g.fwhm_L)+\
                      ", Gaussian_FWHM, "+"{0.value:0.03f}, {0.unit:FITS}".format(g.fwhm_G)+\
                      ", Amplitude, "+"{0.value:0.03f}, {0.unit:FITS}".format(amp))
              # print(self.message)
              self.log.checklog()
              self.log.write(self.message[-1])
              #self.log.write('\n')
              self.outputupdate()

            elif func=="lorentz":
              g_init = models.Lorentz1D(x_0=np.mean(xg),amplitude=np.max(yg)-reflevel, fwhm=xg[1]-xg[0])
              fit_g = fitting.LevMarLSQFitter()
              if invert == False:
                g = fit_g(g_init, xg, yg)
                if self.script == False:
                    self.ax.plot(xgf, g(xgf))#*np.polyval(linecoeff,xgf), label='%s'%(func))
                # self.ax.plot(xg,ygn)
                # self.ax.plot(xgf,g(xgf))
                amp=g.amplitude
              elif invert == True:
                g = fit_g(g_init, xg, reflevel-yg)
                if self.script == False:
                    self.ax.plot(xgf, (reflevel-g(xgf)))#*np.polyval(linecoeff,xgf), label='%s'%(func))
                # self.ax.plot(xg,-ygn)
                # self.ax.plot(xgf,-g(xgf))
                amp=reflevel-g.amplitude

              t=self.filedate()
              self.message.append(t+"Lorentz Center, "+"{0.value:0.03f}, {0.unit:FITS}".format(g.x_0)+", FWHM, "+\
                      "{0.value:0.03f}, {0.unit:FITS}".format(g.fwhm)+\
                      ", Amplitude, "+"{0.value:0.03f}, {0.unit:FITS}".format(amp))
              # print(self.message)
              self.log.checklog()
              self.log.write(self.message[-1])
              #self.log.write('\n')
              #
              self.outputupdate()

            elif func=="moffat":
              g_init = models.Moffat1D(x_0=np.mean(xg),amplitude=np.max(yg)-reflevel, gamma=xg[1]-xg[0],alpha=1)
              fit_g = fitting.LevMarLSQFitter()
              if invert == False:
                g = fit_g(g_init, xg, yg)
                if self.script == False:
                    self.ax.plot(xgf, g(xgf))#*np.polyval(linecoeff,xgf), label='%s'%(func))
                # self.ax.plot(xg,ygn)
                # self.ax.plot(xgf,g(xgf))
                amp=g.amplitude
              elif invert == True:
                g = fit_g(g_init, xg, reflevel-yg)
                if self.script == False:
                    self.ax.plot(xgf, (reflevel-g(xgf)))#*np.polyval(linecoeff,xgf), label='%s'%(func))
                # self.ax.plot(xg,-ygn)
                # self.ax.plot(xgf,-g(xgf))
                amp=reflevel-g.amplitude
              # print('amp',amp)
              # print('fwhm',g.fwhm)
              # print('x0',g.x_0)
              t=self.filedate()
              self.message.append(t+"Moffat Center, "+"{0.value:0.03f}, {0.unit:FITS}".format(g.x_0)+", FWHM, "+\
                      "{0.value:0.03f}, {0.unit:FITS}".format(g.fwhm)+\
                      ", Amplitude, "+"{0.value:0.03f}, {0.unit:FITS}".format(amp))
              # print(self.message)
              self.log.checklog()
              self.log.write(self.message[-1])
              #self.log.write('\n')
              #
              self.outputupdate()

            else:
              print("Error with Fitting Function Selection")

            self.canvas.draw()
        except:
            print('Exception occured in fit')

    def plotRegions(self):
        try:
            rect=patches.Rectangle((min(self.saveregions_x),min(self.flux.value)),max(self.saveregions_x)-min(self.saveregions_x),max(self.flux.value)-min(self.flux.value),edgecolor='b',facecolor='b',alpha=0.15)
            self.ax.add_patch(rect)
            self.canvas.draw()
        except:
            pass
            #print('Exception occured in plotRegions')


    def paropen(self,file):
        # self.region_clear()
        dataout=open(file)
        # func=False
        for i,line in enumerate(dataout):
            if i == 0 :
                line=line.lower()
                line=line.strip()
                dataout.close()
                # print(line)
                if 'bisect' == line:
                    self.LoadBisect(file=file)
                    break
                elif 'norm' == line:
                    self.LoadNorm(file=file)
                    break
                elif 'region' == line:
                    self.LoadRegion(file=file)
                    break
                elif 'align' == line:
                    func='align'
                    break
                else:
                    break
                    # pass
            else:
                break

    def SaveRegion(self):
        """Save the regions used for equivalent width measurements and for fitting line profiles."""
        try:
            path=os.path.dirname(self.fname)
            path_wo_ext=os.path.splitext(self.fname)[0]
            initialname=os.path.join(path_wo_ext,"region.par")
            savename, _ = QtWidgets.QFileDialog.getSaveFileName(self,"Save Normalization Parameters",initialname)
            dataout=open(savename,'w')
            # dataout.write('%s\n'%(self.fname))
            dataout.write('region\n')
            for row in self.saveregions_x:
                dataout.write('%s '%(row))
            dataout.write('\n')
            for row in self.saveregions_y:
                dataout.write('%s '%(row))
            dataout.write('\n')
            dataout.close()
            self.message.append("Region saved to:  %s"%(savename))
            self.outputupdate()
            self.region_clear()
        except:
            self.message.append("Nothing to save.")
            self.outputupdate()


    def LoadRegion(self,file=False):
        """Load the regions used for equivalent width measurements and for fitting line profiles."""
        self.region_clear()
        try:
            if file==False:
                filename, _= QtWidgets.QFileDialog.getOpenFileNames(self,"Open Region File.")
                dataout=open(filename[0])
            else:
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
            # print(self.saveregions_x,self.saveregions_y)
            dataout.close()
            self.loadedregions=True
            self.message.append("Regions loaded from:  %s"%(file))
            self.outputupdate()
            self.plotRegions()
        except:
            self.message.append("Not a properly formatted input file.")
            self.outputupdate()

    def SaveBisect(self):
        """Save the regions for bisection"""
        try:
            path=os.path.dirname(self.fname)
            path_wo_ext=os.path.splitext(self.fname)[0]
            initialname=os.path.join(path_wo_ext,"bisect.par")
            savename, _ = QtWidgets.QFileDialog.getSaveFileName(self,"Save Bisection Parameters",initialname)
            dataout=open(savename,'w')
            # dataout.write('%s\n'%(self.fname))
            dataout.write('bisect\n')
            for row in self.saveregions_x:
                dataout.write('%s '%(row))
            dataout.write('\n')
            # for row in self.saveregions_y:
            #     dataout.write('%s '%(row))
            # dataout.write('\n')
            dataout.write(str(self.height)[1:-1])
            dataout.write('\n')
            dataout.close()
            self.message.append("Bisection Regions Saved to:  %s"%(savename))
            self.outputupdate()
            self.region_clear()
        except:
            self.message.append("Nothing to save.")
            self.outputupdate()


    def LoadBisect(self,file=False):
        """Load the regions bisection."""
        try:
            self.region_clear()
            if file ==  False:
                filename, _ = QtWidgets.QFileDialog.getOpenFileNames(self,"Open Bisect File.")
                dataout=open(filename[0])
            else:
                dataout=open(file)

            for i,line in enumerate(dataout):
                if i == 0 :
                    pass
                if i == 1 :
                    self.saveregions_x=np.array(line.split(),dtype=float)
                    self.saveregions_y=np.array([0,0],dtype=float)
                    # print(line.split()[0])
                # elif i == 2 :
                #     self.saveregions_y=np.array(line.split(),dtype=float)
                elif i == 2 :
                    self.height=[]#float(height)
                    try:
                        hs=list(line.split(" "))
                    except:
                        pass
                    try:
                        hs=list(line.split(","))
                    except:
                        pass
                    for h in hs:
                        self.height.append(float(h))
                    # print(self.height)
                else:
                    pass
            dataout.close()
            self.loadedbisect=True
            self.loadedregions=True
            self.plotRegions()
            self.message.append("Bisection Paremeters Loaded from  %s"%(file))
            self.outputupdate()
        except:
            self.message.append("Not a properly formatted input file.")
            self.outputupdate()

    def SaveAlign(self):
        """Save the regions for alignment"""
        try:
            path=os.path.dirname(self.fname)
            path_wo_ext=os.path.splitext(self.fname)[0]
            initialname=os.path.join(path_wo_ext,"align.par")
            savename, _ = QtWidgets.QFileDialog.getSaveFileName(self,"Save Alignment Parameters",initialname)
            dataout=open(savename,'w')
            # dataout.write('%s\n'%(self.fname))
            dataout.write('align\n')
            for row in self.x_norm:
                dataout.write('%s '%(row))
            dataout.write('\n')
            for row in self.y_norm:
                dataout.write('%s '%(row))
            dataout.write('\n')
            dataout.write(str(self.ref_wave))
            dataout.write('\n')
            dataout.close()
            self.message.append("Alignment Paremeters Saved to:  %s"%(savename))
            self.outputupdate()
            self.region_clear()
        except:
            self.message.append("Nothing to save.")
            self.outputupdate()

    def LoadAlign(self,file=False):
        """Load the regions alignment."""
        try:
            self.region_clear()
            if file == False:
                filename, _ = QtWidgets.QFileDialog.getOpenFileNames(self,"Open Align File.")
                dataout=open(filename[0])
            else:
                dataout=open(file)

            for i,line in enumerate(dataout):
                if i == 0 :
                    pass
                elif i == 1:
                    self.x_norm=np.array(line.split(),dtype=float)
                    # print(line.split()[0])
                elif i == 2 :
                    self.y_norm=np.array(line.split(),dtype=float)
                elif i == 3 :
                    self.ref_wave=float(line.split()[0])
                else:
                    pass
            dataout.close()
            self.loadedalign=True
            self.loadedregions=True
            self.plotRegions()
            self.message.append("Alignment Paremeters Loaded From:  %s, \n Wavelength: %s"%(file,self.ref_wave))
            self.outputupdate()
        except:
            self.message.append("Not a properly formatted input file.")
            self.outputupdate()

    def SaveNorm(self):
        """Save the regions and powerlaw for the normalization."""
        try:
            path=os.path.dirname(self.fname)
            path_wo_ext=os.path.splitext(self.fname)[0]
            initialname=os.path.join(path_wo_ext,"norm.par")
            savename, _ = QtWidgets.QFileDialog.getSaveFileName(self,"Save Normalization Parameters",initialname)
            dataout=open(savename,'w')
            # dataout.write('%s\n'%(self.fname))
            dataout.write('norm\n')
            dataout.write('%s\n'%(self.order))
            for row in self.x_norm:
                dataout.write('%s '%(row))
            dataout.write('\n')
            # for row in self.y_norm:
            #     dataout.write('%s '%(row))
            # dataout.write('\n')
            dataout.close()
            self.message.append("Normalization Parameters Saved to:  %s"%(savename))
            self.outputupdate()
            self.region_clear()
        except:
            self.message.append("Nothing to save.")
            self.outputupdate()

    def LoadNorm(self,file=False):
        """Load the regions and powerlaw for the normalization."""
        self.message.append("Function Not working at this time.")
        self.outputupdate()
        try:
            self.region_clear()
            if file == False:
                filename, _ = QtWidgets.QFileDialog.getOpenFileNames(self,"Open Normalization File.")
                dataout=open(filename[0])
            else:
                dataout=open(file)

            for i,line in enumerate(dataout):
                if i == 0 :
                    pass
                if i == 1:
                    self.order=int(line)
                elif i >= 2 :
                    self.x_norm=np.array(line.split(),dtype=float)
                    # print(line.split()[0])
                # elif i == 3 :
                #     self.y_norm=np.array(line.split(),dtype=float)
                else:
                    pass
            # print(self.x_norm)
            dataout.close()
            self.loadednorm=True
            # self.plotRegions()
            self.plotcontinuum()

            self.message.append("Normalization Parameters Loaded From:  %s, \n Polynomial Order: %s"%(file,self.order))
            self.outputupdate()
        except:
            self.message.append("Not a properly formatted input file.")
            self.outputupdate()

    def mouseclick_cont(self,event):
        try:
            self.x.append(event.xdata)
            self.y.append(event.ydata)
        except:
            print('Exception occured in mouseclick_cont')

    def continuum(self,message=None):
        """uses two clicks to define a region for fitting or measuring."""
        try: #if it doesn't exist, create it.
            self.x_norm
        except:
            self.x_norm=[]
            self.y_norm=[]
        try:
            if (len(self.x) % 2) != 0: #if not even, pop last value.  must have matched pairs of points.
                self.x.pop()
                self.y.pop()
        except:
            print('Exception occured in continuum')

        try:
            if self.firstclick == False:
                if message == None or message == False:
                    regionmessage="Define continuum region: Click on left and right edges of your regions. Press c again to stop."
                else:
                    regionmessage=message
                self.message.append(regionmessage)
                self.outputupdate()
                self.click=self.fig.canvas.mpl_connect('button_press_event', self.mouseclick_cont)
                self.firstclick=True
            else:
                self.fig.canvas.mpl_disconnect(self.click)
                self.firstclick=False
                self.continuum_region()
        except:
            print('Exception occured in continuum')

    def continuum_region(self):
        for i,value in enumerate(self.x):
            if (i % 2) == 0: #every other pair reprsents a region of continuum.
                try:
                    xg,yg=self.chop(self.wavelength,self.flux,self.x[i],self.x[i+1])
                except:
                    pass
            for xi in xg:
                self.x_norm.append(xi.value)
            for yi in yg:
                self.y_norm.append(yi.value)
        self.plotcontinuum()

    def plotcontinuum(self):
        self.ax.plot(self.x_norm,self.y_norm,'ks')
        self.canvas.draw()

    def getCorrespondingFlux(self,array):
        out=[]
        for value in array:
            idx=find_nearest_index(self.wavelength,value)
            out.append(self.flux[idx].value)
        return out

    def normalize(self):
        """Continuum normalize by using selected points as continuum."""
        try:
            # self.measuremode()
            self.goodfit=False

            try:
                xn=np.array(self.x_norm)
                self.y_norm=self.getCorrespondingFlux(xn)
                yn=np.array(self.y_norm)
                self.canvas.draw()
                self.loadednorm=True
            except:
                self.message.append("Please select continuum region (c) first.")
                self.outputupdate()
                stop=True
                self.loadednorm=False

            while self.goodfit == False:
                if self.loadednorm == True:
                    if self.orderset == False:
                        i, okPressed = QInputDialog.getInt(self, "Integer","Polynomial Order:", self.order, 1, 20, 1)
                        if okPressed:
                            self.order=i
                            self.orderset=True
                        else:
                            # self.goodfit=True
                            self.orderset=False
                            break

                    linecoeff = np.polyfit(xn,yn,self.order)
                    nflux=self.flux/np.polyval(linecoeff,self.wavelength.value)

                    if self.script == False:
                        # self.reset()
                        self.normtest,=self.ax.plot(self.wavelength,nflux)
                        self.ax.set_ylim([max(0,min(nflux.value)),min(100,max(nflux.value))]) #need to change this to be over the selected region.
                        self.canvas.draw()
                        answer=QMessageBox.question(self, 'Question', "Proceed with the fit?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
                        if answer == QMessageBox.Yes:
                           self.ax.lines.remove(self.normtest)
                           self.flux=nflux
                           self.splot()
                           self.ax.grid(self.gridvalue)
                           self.canvas.draw()
                           # self.stackpane()
                           # self.plotRegions()
                           self.Spectra.updatespectrum()
                           self.loadednorm=True
                           self.orderset=True
                           # self.reset()
                           self.goodfit = True
                           # break
                        else:
                           self.ax.lines.remove(self.normtest)
                           # self.reset()
                           self.canvas.draw()
                           self.orderset=False
                           self.goodfit = True
                           # break

                    elif self.script == True and self.orderset == True:
                        self.flux=nflux
                        self.Spectra.updatespectrum()
                        self.goodfit = True
                        # break
                    else:
                        pass

                else:
                    self.goodfit=True
        except:
            print('Exception occured in normalize')

    def saveImage(self):
        try:
            # self.stackforsaving=[]
            if self.suffix == False:
                self.suffix, okPressed = QInputDialog.getText(self, "Save Image","File format, any matplotlib accepted filetype.", QLineEdit.Normal, "")
            if '.' not in self.suffix:
                self.suffix='.'+self.suffix
            # for f in self.database:
            self.overplot=False
            self.stackplot=False
            # self.ax.clear()
            self.plotSpectra(spec=self.fname,savename=self.fname[:-5]+self.suffix)

        except:
            print('Exception occured in savepng')

    def savefits(self):
        try:
            # self.stackforsaving=[]
            if self.suffix == False:
                self.suffix, okPressed = QInputDialog.getText(self, "Save Stack Spectra","File suffix.", QLineEdit.Normal, "")

            self.overplot=False
            self.stackplot=False
            self.plotSpectra(spec=self.fname)
            self.Spectra.saveFits(extend=self.suffix+'.fits')

        except:
            print('Exception occured in MainWin.savefits')


    def savetext(self):
        try:
            # self.stackforsaving=[]
            if self.suffix == False:
                self.suffix, okPressed = QInputDialog.getText(self, "Save Stack Spectra","File suffix.", QLineEdit.Normal, "")
            self.overplot=False
            self.stackplot=False
            self.plotSpectra(spec=self.fname)
            self.Spectra.save1DText(extend=self.suffix+'.txt')

        except:
            print('Exception occured in MainWin.savetext')

#

    def stackadd(self,stack=False,outname='AddSpec'):
        '''Add spectra'''
        while outname in self.database:
            outname=outname+'1'
        try:
            wlist=[]
            flist=[]
            if stack == False:
                stack=self.database
                print('oops')

            for i,f in enumerate(stack):
                self.fname=f
                self.grabSpectra(self.fname)
                wlist.append(self.wavelength.value)
            wflat=[y for x in wlist for y in x]
            winterp=np.arange(min(wflat),max(wflat),abs(wflat[1]-wflat[0]))
            addflux=[]
            for i,f in enumerate(stack):
                self.fname=f
                self.grabSpectra(self.fname)
                newflux=np.interp(winterp,self.wavelength.value,self.flux.value)
                if i == 0:
                    addflux=newflux
                else:
                    addflux=addflux+newflux

            self.database[outname]={}
            self.database[outname]['wavelength']=winterp*self.wavelength.unit
            self.database[outname]['wavelength_orig']=winterp*self.wavelength.unit
            self.database[outname]['flux']=addflux*u.flx
            self.database[outname]['flux_orig']=addflux*u.flx
            self.database[outname]['header']=self.header
            self.Spectra.updatespectrum()
            self.stackrebuild()
            self.ax.clear()
            self.plotSpectra(spec=outname)
            self.reset()
            self.singleplottoggle()
            # self.plotRegions()
        except:
            pass

    def stackadd_const(self,stack=False,const=0,outname='AddSpec'):
        while outname in self.database:
            outname=outname+'1'
        try:
            if stack == False:
                stack=self.database
            self.fname=stack[0]
            self.grabSpectra(self.fname)
            newflux=self.flux.value+const
            self.database[outname]={}
            self.database[outname]['wavelength']=self.wavelength
            self.database[outname]['wavelength_orig']=self.wavelength
            self.database[outname]['flux']=newflux*u.flx
            self.database[outname]['flux_orig']=newflux*u.flx
            self.Spectra.updatespectrum()
            self.stackrebuild()
            self.ax.clear()
            self.plotSpectra(spec=outname)
            self.reset()
            self.singleplottoggle()
        except:
            print('Exception occured in stackadd_const')

    def stacksubtract(self,stack=False,outname='SubSpec'):
        '''Subtract spectra, subtract second spectrum from first spectrum.'''
        while outname in self.database:
            outname=outname+'1'
        try:
            if stack == False:
                stack=self.database
                print('oops')
            wlist=[]
            flist=[]
            for i,f in enumerate(stack):
                self.fname=f
                self.grabSpectra(self.fname)
                wlist.append(self.wavelength.value)
            wflat=[y for x in wlist for y in x]
            winterp=np.arange(min(wflat),max(wflat),abs(wflat[1]-wflat[0]))
            addflux=[]
            for i,f in enumerate(stack):
                self.fname=f
                self.grabSpectra(self.fname)
                newflux=np.interp(winterp,self.wavelength.value,self.flux.value)
                if i == 0:
                    addflux=newflux
                else:
                    addflux=addflux-newflux

            self.database[outname]={}
            self.database[outname]['wavelength']=winterp*self.wavelength.unit
            self.database[outname]['wavelength_orig']=winterp*self.wavelength.unit
            self.database[outname]['flux']=addflux*u.flx
            self.database[outname]['flux_orig']=addflux*u.flx
            self.database[outname]['header']=self.header
            self.Spectra.updatespectrum()
            self.stackrebuild()
            self.ax.clear()
            self.plotSpectra(spec=outname)
            self.reset()
            self.singleplottoggle()
        except:
            print('Exception occured in stacksubtract')

    def stackcombine(self,stack=False,func='ave',outname='Combined'):
        '''Average spectra'''
        while outname in self.database:
            outname=outname+'1'
        print('stack')
        try:
            wlist=[]
            flist=[]
            if stack == False:
                stack=self.database
            for i,f in enumerate(stack):
                self.fname=f
                self.grabSpectra(self.fname)
                wlist.append(self.wavelength.value)
            wflat=[y for x in wlist for y in x]
            winterp=np.arange(min(wflat),max(wflat),abs(wflat[1]-wflat[0]))
            addflux=[]
            for i,f in enumerate(stack):
                self.fname=f
                self.grabSpectra(self.fname)
                newflux=np.interp(winterp,self.wavelength.value,self.flux.value)
                addflux.append(newflux)
            if func == 'ave':
                addflux=np.mean(np.array(addflux),axis=0)
                # self.database[outname]={}
                # self.database[outname]['wavelength']=winterp*self.wavelength.unit
                # self.database[outname]['wavelength_orig']=winterp*self.wavelength.unit
                # self.database[outname]['flux']=addflux*u.flx
                # self.database[outname]['flux_orig']=addflux*u.flx
            elif func == 'median':
                addflux=np.median(np.array(addflux),axis=0)

            self.database[outname]={}
            self.database[outname]['wavelength']=winterp*self.wavelength.unit
            self.database[outname]['wavelength_orig']=winterp*self.wavelength.unit
            self.database[outname]['flux']=addflux*u.flx
            self.database[outname]['flux_orig']=addflux*u.flx
            self.database[outname]['header']=self.header

            self.Spectra.updatespectrum()
            self.stackrebuild()
            self.ax.clear()
            self.plotSpectra(spec=outname)
            self.reset()
            self.singleplottoggle()
            # self.plotRegions()
        except:
            print('Exception occured in stackcombine')


    def divide_const(self,spec=False,divisor=0,outname='DividedSpec'):
        """divde a spectrum by a constant"""
        while outname in self.database:
            outname=outname+'1'
        try:
            if spec != False:
                self.fname=spec
            else:
                spec=self.fname
            self.grabSpectra(self.fname)
            if divisor == 0:
                divisor, okPressed = QInputDialog.getDouble(self, "Divisor","divide by", 1, -100000, 100000, 3)
            if divisor != 0:
                newflux=self.flux.value/divisor

                self.database[outname]={}
                self.database[outname]['wavelength']=self.wavelength
                self.database[outname]['wavelength_orig']=self.wavelength
                self.database[outname]['flux']=newflux*u.flx
                self.database[outname]['flux_orig']=newflux*u.flx
                self.Spectra.updatespectrum()
                self.stackrebuild()
                self.ax.clear()
                self.plotSpectra(spec=outname)
                self.reset()
                self.singleplottoggle()
        except:
            print('Exception occured in divide_const')

    def mult_const(self,spec=False,multiplier=0,outname='MultSpec'):
        """divde a spectrum by a constant"""
        while outname in self.database:
            outname=outname+'1'
        try:
            if spec != False:
                self.fname=spec
            else:
                spec=self.fname
            self.grabSpectra(self.fname)
            if multiplier == 0:
                multiplier, okPressed = QInputDialog.getDouble(self, "Multiplier","Multiply by", 1, -100000, 100000, 3)
            if multiplier != 0:
                newflux=self.flux.value*multiplier

                self.database[outname]={}
                self.database[outname]['wavelength']=self.wavelength
                self.database[outname]['wavelength_orig']=self.wavelength
                self.database[outname]['flux']=newflux*u.flx
                self.database[outname]['flux_orig']=newflux*u.flx
                self.Spectra.updatespectrum()
                self.stackrebuild()
                self.ax.clear()
                self.plotSpectra(spec=outname)
                self.reset()
                self.singleplottoggle()
        except:
            print('Exception occured in mult_const')

    def dividespec(self,stack=False,outname='Divided'):
        '''Divide one spectrum by another. First in stack is numerator.'''
        while outname in self.database:
            outname=outname+'1'
        try:
            if stack == False:
                print('stack error in divide spec')
            wlist=[]
            flist=[]
            for i,f in enumerate(stack):
                self.fname=f
                self.grabSpectra(self.fname)
                wlist.append(self.wavelength.value)
            wflat=[y for x in wlist for y in x]
            winterp=np.arange(min(wflat),max(wflat),abs(wflat[1]-wflat[0]))
            divflux=[]
            for i,f in enumerate(stack):
                self.fname=f
                self.grabSpectra(self.fname)
                newflux=np.interp(winterp,self.wavelength.value,self.flux.value)
                if i == 0:
                    numerator=newflux
                else:
                    divflux=numerator/newflux

            self.database[outname]={}
            self.database[outname]['wavelength']=winterp*self.wavelength.unit
            self.database[outname]['wavelength_orig']=winterp*self.wavelength.unit
            self.database[outname]['flux']=divflux*u.flx
            self.database[outname]['flux_orig']=divflux*u.flx
            self.database[outname]['header']=self.header
            self.Spectra.updatespectrum()
            self.stackrebuild()
            self.ax.clear()
            self.plotSpectra(spec=outname)
            self.reset()
            self.singleplottoggle()
        except:
            print('Exception occured in dividespec')


    def multiplyspec(self,stack=False,outname='Multiplied'):
        '''Multiplied one spectrum by another.'''
        while outname in self.database:
            outname=outname+'1'
        try:
            if stack == False:
                print('stack error in multiply spec')
            wlist=[]
            flist=[]
            for i,f in enumerate(stack):
                self.fname=f
                self.grabSpectra(self.fname)
                wlist.append(self.wavelength.value)
            wflat=[y for x in wlist for y in x]
            winterp=np.arange(min(wflat),max(wflat),abs(wflat[1]-wflat[0]))
            divflux=[]
            for i,f in enumerate(stack):
                self.fname=f
                self.grabSpectra(self.fname)
                newflux=np.interp(winterp,self.wavelength.value,self.flux.value)
                if i == 0:
                    numerator=newflux
                else:
                    divflux=numerator*newflux

            self.database[outname]={}
            self.database[outname]['wavelength']=winterp*self.wavelength.unit
            self.database[outname]['wavelength_orig']=winterp*self.wavelength.unit
            self.database[outname]['flux']=divflux*u.flx
            self.database[outname]['flux_orig']=divflux*u.flx
            self.database[outname]['header']=self.header
            self.Spectra.updatespectrum()
            self.stackrebuild()
            self.ax.clear()
            self.plotSpectra(spec=outname)
            self.reset()
            self.singleplottoggle()
        except:
            print('Exception occured in multiplyspec')


    def linearize(self):
        try:
            wlist=[]
            flist=[]
            # for i,f in enumerate(self.database):

            self.grabSpectra(self.fname)
            wlist.append(self.wavelength.value)
            wflat=[y for x in wlist for y in x]
            winterp=np.arange(min(wflat),max(wflat),abs(wflat[1]-wflat[0]))
            addflux=[]
            # for i,f in enumerate(self.database):
            #     self.fname=f
            # self.grabSpectra(self.fname)
            addflux=np.interp(winterp,self.wavelength.value,self.flux.value)
            # if i == 0:
            #     addflux=newflux
            # else:
            #     addflux=addflux+newflux

            self.database[self.fname]['wavelength']=winterp*self.wavelength.unit
            self.database[self.fname]['flux']=addflux*u.flx
            # try:
            #     self.database[self.fname]['header']=header
            # except:
            #     pass
            self.Spectra.updatespectrum()
            self.stackrebuild()
            self.ax.clear()
            self.plotSpectra(spec=self.fname)
            # self.plotRegions()
        except:
            print('Exception occured in linearize')


    def abort_stack(self):
        #Print a message to user to indicate reason for fault.
        self.message.append("Please select region.")
        self.outputupdate()

    def stackrebuild(self):
        self.stack=[] #rebuild stack everytime so that we have a list of all files open.
        if self.poplast == True:
            self.database.pop(self.fname)
            self.poplast=False
        t=len(self.database)
        for i,row in enumerate(self.database):
            self.database[row]['stacknumber']=i+1
            self.database[row]['plotcolor']=cm.gist_rainbow(i/t)
            self.stack.append(row)

    def stack_startmessage(self):
        self.message.append("Stack Operation Started Please Wait")
        self.outputupdate()

    def stacker(self,func=False):
        """A helper function to automate tasks for many spectra."""
        self.stack_startmessage()
        # self.message="Stack Operation in Progress"
        # self.outputupdate()
        # # self.canvas.draw()
        # self.stackforsaving=[]
        self.stacknum=0
        self.script=True
        self.suffix=False
        if func == False:
            print("No function selected for stacker.")
        else:
            for f in self.database:
                self.fname=f
                # self.measuremode()
                # self.ax.clear()
                # self.selectSpectra(spec=self.fname)
                self.grabSpectra(self.fname)
                if func == "norm":
                    self.normalize()
                elif func == "scopy":
                    self.scopy()
                elif func == "eqw":
                    self.eqw()
                elif func == "gauss":
                    self.fit(func="gauss")
                elif func == "voigt":
                    self.fit(func="voigt")
                elif func == "lorentz":
                    self.fit(func="lorentz")
                elif func == "moffat":
                    self.fit(func="moffat")
                elif func == "bisect":
                    self.BisectLine()
                elif func == "align":
                    self.align()
                    # self.saveFits(extend='-align.fits')
                elif func == "restoreall":
                    self.Spectra.restore()
                elif func == "snr":
                    self.signal2noise()
                elif func == "coord":
                    self.coord()
                elif func == "boxcar":
                    self.smooth(func="boxcar")
                elif func == "gsmooth":
                    self.smooth(func="gaussian")
                elif func == "savepng":
                    self.saveImage()
                elif func == "savefits":
                    self.savefits()
                elif func == "savetext":
                    self.savetext()
                else:
                    print("Stacker cannot handle a function.")
                self.stacknum=self.stacknum+1

        if func == "norm":
            self.reset()
        self.message.append("Stack Operation Complete")
        self.outputupdate()
        self.script=False


    def stackwindowplot(self,spec):
        self.singleplottoggle()
        self.ax.clear()
        # self.measuremode()
        self.plotSpectra(spec=spec)


    def stackdown(self):
        '''Move the display down one spectrum in the stack.'''
        self.singleplottoggle()
        self.ax.clear()
        # self.measuremode()
        spec=self.stack[self.stack.index(self.fname)-1]
        # print('test ',self.database(0))
        self.plotSpectra(spec=spec)
        # self.plotRegions()

    def stackup(self):
        '''Move the display up one spectrum in the stack.'''
        self.singleplottoggle()
        self.ax.clear()
        # self.measuremode()
        try:
            spec=self.stack[self.stack.index(self.fname)+1]
        except:
            spec=self.stack[0]
        self.plotSpectra(spec=spec)
        # self.plotRegions()


    def gotostack(self):
        # self.singleplottoggle()
        # self.measuremode()
        self.ax.clear()
        # self.slist=[]
        selected, okPressed = QInputDialog.getText(self, "Plot Selected Spectra","Spectrum Number(s) (1,4), Ranges (3-4)", QLineEdit.Normal, "")
        # print(list(selected))
        self.selection=list(selected.split(","))
        # sel=list(selected.split(","))
        if okPressed:
            self.overplot=True

            try:
                self.minilist()
                self.plotSpectra(spec=self.slist)
            except:
                pass

            # self.plotRegions()
        else:
            pass


    def minilist(self):
        self.slist=[]
        for s in self.selection:
            if '-' in s:
                x,y=s.split("-")
                srange=np.arange(int(x),int(y)+1)
                for k in srange:
                    # print(k)
                    self.slist.append(self.stack[int(k)-1])
            elif s == '0':
                pass

            else:
                self.slist.append(self.stack[int(s)-1])


    def stackwindowopenspectra(self):
        try:
            self.openSpectra()
        except:
            print('Exception occured in stackwindowopenspectra')

    def delfromstack(self):
        try:
            # self.singleplottoggle()
            # self.measuremode()
            # self.ax.clear()
            # i, okPressed = QInputDialog.getInt(self, "Integer","Plot Spectrum by Stack Number:", 1, 1, 50, 1)
            selected, okPressed = QInputDialog.getText(self, "Delete Selected Spectra","Spectrum Number(s) (1,4), Ranges (3-4)", QLineEdit.Normal, "")
            # print(list(selected))
            self.selection=list(selected.split(","))
            self.minilist()
            if okPressed:
                for s in self.slist:
                    try:
                        self.stack.remove(s)
                        self.database.pop(s)
                    except:
                        pass

                self.stackrebuild()
                self.ax.clear()
                self.singleplottoggle()
                # self.measuremode()
                self.canvas.draw()
                self.stackpane()
        except:
            print('Exception occured in delfromstack')

    def removefromstack(self):
        try:
            self.stack.remove(self.fname)
            self.database.pop(self.fname)
        except:
            self.message.append("Trouble Deleting from Stack.")
            self.outputupdate()
        self.ax.clear()
        self.canvas.draw()
        self.stackpane()
        self.message.append("You may replot the stack (]), overplot ([), or choose a new spectrum from the stack.")
        self.outputupdate()

    def sortby(self,keyword=None):
        pass
        #should have pop up window, choose or enter keyword to sort by
        #sort order, smallest to largest or vice versa
        #
        # if keyword != None:
        #     for i,row in enumerate(self.database):
        #         self.header=self.database[row]['header']
        #         try:
        #             jd.append(self.header['JD'])
        #         except:
        #             jd.append(i)
    def phasefold(self,period=None,T0=None):
        if peirod != None and T0 != None:
            for jd in self.jd:
                phase=(jd-T0)/period
            #add to header, yes
            #again, only a calculation, not a sort.
            #should have a separate sort function that will sort based on header entries
            #this is more general and flexible.
