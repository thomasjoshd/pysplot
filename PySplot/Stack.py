import os
import sys
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMessageBox,QInputDialog, QLineEdit
from PyQt5.QtCore import Qt, QSize


class Stack(QtWidgets.QMainWindow):
    def __init__(self,parent):
        super(Stack, self).__init__(parent)

    # def stackwindowopenspectra(self):
    #     try:
    #         self.openSpectra()
    #     except:
    #         print('Exception occured in stackwindowopenspectra')

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
                        self.parent().database.pop(s)
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
            self.parent().database.pop(self.fname)
        except:
            self.parent().message.append("Trouble Deleting from Stack.")
            self.parent().outputupdate()
        self.ax.clear()
        self.canvas.draw()
        self.stackpane()
        self.parent().message.append("You may replot the stack (]), overplot ([), or choose a new spectrum from the stack.")
        # self.parent().outputupdate()

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
        # print('test ',self.parent().database(0))
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
#
    def abort_stack(self):
        #Print a message to user to indicate reason for fault.
        self.parent().message.append("Please select region.")
        self.parent().outputupdate()

    def stackrebuild(self):
        self.stack=[] #rebuild stack everytime so that we have a list of all files open.
        if self.poplast == True:
            self.parent().database.pop(self.fname)
            self.poplast=False
        t=len(self.parent().database)
        for i,row in enumerate(self.parent().database):
            self.parent().database[row]['stacknumber']=i+1
            self.parent().database[row]['plotcolor']=cm.gist_rainbow(i/t)
            self.stack.append(row)

    def stack_startmessage(self):
        self.parent().message.append("Stack Operation Started Please Wait")
        self.parent().outputupdate()

    def stacker(self,func=False):
        """A helper function to automate tasks for many spectra."""
        self.stack_startmessage()
        # self.message="Stack Operation in Progress"
        # self.parent().outputupdate()
        # # self.canvas.draw()
        # self.stackforsaving=[]
        self.stacknum=0
        self.script=True
        self.suffix=False
        if func == False:
            print("No function selected for stacker.")
        else:
            for f in self.parent().database:
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
        self.parent().message.append("Stack Operation Complete")
        self.parent().outputupdate()
        self.script=False
