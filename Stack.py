import os
import sys
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMessageBox,QInputDialog, QLineEdit
from PyQt5.QtCore import Qt, QSize


class Stack(QtWidgets.QMainWindow):
    def __init__(self,parent):
        super(Stack, self).__init__(parent)

    def stackwindowopenspectra(self):
        try:
            self.openSpectra()
        except:
            pass

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
                self.measuremode()
                self.canvas.draw()
                self.stackpane()
        except:
            pass

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
