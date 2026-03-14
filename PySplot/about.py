# from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMessageBox
# from PyQt5.QtCore import Qt, QSize
from version import VERSION,UPDATED

def about():
    """About window"""
    t=QMessageBox()
    t.setWindowTitle("About")
    L1="PySplot Version %s\n"%VERSION
    L2="Last Updated %s\n\n"%UPDATED
    L3="Author: Dr. Joshua Thomas \n Clarkson University \n jthomas@clarkson.edu \n thomas.joshd@gmail.com\n\n"
    L4="This program was designed to emulate some basic IRAF splot functions for 1-D spectra.\n"
    L5="I made use of Astropy, Matplotlib, Scipy, and PyQT5, among others.\n"
    L6="Check for updates at:  https://people.clarkson.edu/~jthomas/pysplot.html\n"
    L7="Known issue: overwriting fits files does not work."
    t.setText(L1+L2+L3+L4+L5+L6+L7)
    # t.setStyleSheet("background-color: white")
    t.setStandardButtons(QMessageBox.Close)
    t.exec_()
