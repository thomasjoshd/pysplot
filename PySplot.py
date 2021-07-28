# import os
import sys

from PyQt5 import QtWidgets

from MainWin import MainWin

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    # sys.excepthook = traceback.print_exception()
    # creating main window
    mw = MainWin()
    mw.show()
    sys.exit(app.exec())


#ToDoList
    #add in the telluric-new.py functionality
    #add convoluion by given profile. put this as a function in the sarith window?
    #cross correction
    #open echelle format spectra: hopeless?
    #save images
    #save animated gif
        #https://stackoverflow.com/questions/45852424/how-to-animate-a-function-in-matplotlib
    #set custom titles, labels
    #insert label text annotations
    #sort stack by julian date
    #phase fold data for dynamical?
    #Stack Modify add:
        #linearize
        #linear shift
        #convert wavelength

#Stuff for big rewrite
    #overwrite files doesn't work right, seems to be an astropy issue.
    #split into multiple classes -- in progress
    #finish adding status tips (though these don't seem to show on a mac, or linux.)
    #split fit function into separate functions. (save for big rewrite)
    #better display of line fitting. (longer term)
    #fit double gaussians (longer term, not for this release)
