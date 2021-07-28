from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMessageBox,QInputDialog, QLineEdit
from PyQt5.QtCore import Qt, QSize

from functools import partial


class Menu(QtWidgets.QMainWindow):
    def __init__(self,parent):
        super(Menu, self).__init__(parent)

    def menu(self):
        # filling up a menu bar
        bar = self.parent().menuBar()
        filemenu = bar.addMenu('File')

        fileopen = QtWidgets.QAction('Open 1-D Spectra (.fit(s), .txt, .s, .list)', self.parent())
        fileopen.setShortcut('Ctrl+o')
        fileopen.setStatusTip('Open a new Fits, or list of spectra.')
        fileopen.triggered.connect(self.parent().Spectra.open)

        restore=QtWidgets.QAction("Restore Original",self.parent())
        restore.setShortcut('Shift+o')
        restore.setStatusTip('Reloads original opened spectrum')
        restore.triggered.connect(self.parent().Spectra.restore)

        rname=QtWidgets.QAction("Rename Currently Open Spectra",self.parent())
        rname.triggered.connect(self.parent().rname)

        header=QtWidgets.QAction("FITS Header",self.parent())
        header.setShortcut('h')
        header.triggered.connect(self.parent().imhead)

        exit=QtWidgets.QAction("Exit", self.parent())
        exit.setShortcut('Ctrl+q')
        exit.triggered.connect(self.parent().close_out)

        closeall=QtWidgets.QAction("Close All Spectra",self.parent())
        closeall.setShortcut('Ctrl+w')
        closeall.triggered.connect(self.parent().stackreset)

        writelog=QtWidgets.QAction("Write Log/Start New Log ",self.parent())
        # writelog.setShortcut('')
        writelog.triggered.connect(self.parent().log.endoflog)

        savetxt=QtWidgets.QAction("Save Headerless Text Spectrum",self.parent())
        savetxt.setShortcut('Shift+Ctrl+s')
        savetxt.triggered.connect(self.parent().Spectra.save1DText)

        savefits=QtWidgets.QAction("Save New Fits",self.parent())
        savefits.setShortcut('Ctrl+s')
        savefits.triggered.connect(self.parent().Spectra.saveFits)

        #list the actions in the menu
        filemenu.addAction(fileopen)
        filemenu.addAction(restore)
        filemenu.addAction(header)
        filemenu.addAction(rname)
        filemenu.addSeparator()
        filemenu.addAction(savefits)
        filemenu.addAction(savetxt)
        filemenu.addSeparator()
        filemenu.addAction(closeall)
        filemenu.addAction(writelog)
        filemenu.addAction(exit)
        # filemenu.addAction(QtWidgets.QAction("Image Header (h)",self.parent())) #needs popup window

        viewmenu = bar.addMenu('View')
        reset=QtWidgets.QAction("Reset View",self.parent())
        reset.setShortcut('r')
        reset.setStatusTip('Resets the spectrum view.  Removing any fit, or measurement marks.')
        reset.triggered.connect(self.parent().reset)

        gridtoggle=QtWidgets.QAction('Grid Toggle',self.parent())
        gridtoggle.setShortcut('|')
        gridtoggle.setStatusTip('Toggle the gridlines on the plot.')
        gridtoggle.triggered.connect(self.parent().gridtoggle)

        legtog=QtWidgets.QAction('Legend Toggle',self.parent())
        legtog.setShortcut('=')
        legtog.triggered.connect(self.parent().legendtoggle)

        titletog=QtWidgets.QAction('Title Toggle',self.parent())
        titletog.setShortcut('+')
        titletog.triggered.connect(self.parent().titletoggle)

        overplottoggle=QtWidgets.QAction("Over Plot",self.parent())
        overplottoggle.setShortcut('[')
        overplottoggle.setStatusTip('Toggle Over plot mode.')
        overplottoggle.triggered.connect(self.parent().overplottoggle)

        stackplottoggle=QtWidgets.QAction("Stack Plot",self.parent())
        stackplottoggle.setShortcut(']')
        stackplottoggle.setStatusTip('Toggle stack plot mode.')
        stackplottoggle.triggered.connect(self.parent().stackplottoggle)

        singleplottoggle=QtWidgets.QAction("single Plot",self.parent())
        singleplottoggle.setShortcut('\\')
        singleplottoggle.setStatusTip('Toggle single plot mode.')
        singleplottoggle.triggered.connect(self.parent().singleplottoggle)

        style_line=QtWidgets.QAction("Plot as Line",self.parent())
        style_line.setShortcut('i')
        style_line.triggered.connect(self.parent().set_style_line)

        style_point=QtWidgets.QAction("Plot as Points",self.parent())
        style_point.setShortcut('o')
        style_point.triggered.connect(self.parent().set_style_point)

        style_step=QtWidgets.QAction("Plot as Steps 'Histogram'",self.parent())
        style_step.setShortcut('p')
        style_step.triggered.connect(self.parent().set_style_step)


        # viewmenu.addAction(QtWidgets.QAction("Dynamical--WIP (|)",self.parent()))

        viewmenu.addAction(reset)
        viewmenu.addSeparator()
        viewmenu.addAction(gridtoggle)
        viewmenu.addAction(legtog)
        viewmenu.addAction(titletog)
        viewmenu.addSeparator()
        viewmenu.addAction(overplottoggle)
        viewmenu.addAction(stackplottoggle)
        viewmenu.addAction(singleplottoggle)
        viewmenu.addSeparator()
        viewmenu.addAction(style_line)
        viewmenu.addAction(style_point)
        viewmenu.addAction(style_step)

        #viewmenu.addAction(dynamical) #likely to break uses a new window, need to replace tkinter

        modmenu = bar.addMenu('Modify')
        setcontinuum=QtWidgets.QAction("Set Continuum Range", self.parent())
        setcontinuum.setShortcut('c')
        setcontinuum.setStatusTip('Set a range of points to be teated as continuum.')
        setcontinuum.triggered.connect(self.parent().continuum)

        normalize=QtWidgets.QAction("Normalize", self.parent())
        normalize.setShortcut('t')
        normalize.setStatusTip('Fits sets of contiuua with a polynomial of specified order.')
        normalize.triggered.connect(self.parent().normalize)

        normsave=QtWidgets.QAction("Save Norm Parameters",self.parent())
        normsave.setShortcut('Shift+Ctrl+n')
        normsave.triggered.connect(self.parent().SaveNorm)

        normload=QtWidgets.QAction("Load Norm Parameters",self.parent())
        normload.setShortcut('Alt+Ctrl+n')
        normload.triggered.connect(self.parent().LoadNorm)


        smooth=QtWidgets.QAction("Boxcar Smooth", self.parent())
        smooth.setShortcut('Ctrl+b')
        smooth.setStatusTip('Apply Boxcar smoothing with chosen integer.')
        smooth.triggered.connect(partial(self.parent().Spectra.smooth,func="boxcar"))

        gsmooth=QtWidgets.QAction("Gaussian Smooth", self.parent())
        gsmooth.setShortcut('Ctrl+g')
        gsmooth.setStatusTip('Apply Gaussian smoothing with chosen stdev.')
        gsmooth.triggered.connect(partial(self.parent().Spectra.smooth,func="gaussian"))

        align=QtWidgets.QAction("Align Gaussian to a Wavelength", self.parent())
        align.setShortcut('Ctrl+a')
        align.setStatusTip('Linear offset spectra to align features')
        align.triggered.connect(self.parent().align)

        calign=QtWidgets.QAction("Align Click to a Wavelength", self.parent())
        calign.setShortcut('2')
        calign.setStatusTip('Linear offset spectra to align features')
        calign.triggered.connect(self.parent().clickalign)


        crop=QtWidgets.QAction("Crop Spectra", self.parent())
        crop.setShortcut('Ctrl+c')
        crop.setStatusTip('Crop Spectra, for saving to new file.')
        crop.triggered.connect(self.parent().scopy)

        lin=QtWidgets.QAction("Resample to be Linear",self.parent())
        lin.triggered.connect(self.parent().linearize)

        shift=QtWidgets.QAction("Linearly Shift Wavelengths",self.parent())
        shift.triggered.connect(self.parent().linrescale)

        velocity=QtWidgets.QAction("Convert Wavelength <-> Velocity", self.parent())
        velocity.setShortcut('u')
        velocity.setStatusTip('Convert horizontal axis between wavelength and velocity.')
        velocity.triggered.connect(self.parent().velocity)

        arith=QtWidgets.QAction("Spectrum Arithmetic",self.parent())
        arith.setShortcut('a')
        arith.triggered.connect(self.parent().imarith)

        dop=QtWidgets.QAction("Apply Doppler Shift",self.parent())
        # dop.setShortcut('Shift+Ctrl+s')
        dop.triggered.connect(self.parent().Spectra.doppler)


        modmenu.addAction(setcontinuum)
        modmenu.addAction(normalize)
        modmenu.addAction(normsave)
        modmenu.addAction(normload)
        modmenu.addSeparator()
        modmenu.addAction(smooth)
        modmenu.addAction(gsmooth)
        modmenu.addSeparator()
        modmenu.addAction(align)
        modmenu.addAction(calign)
        modmenu.addSeparator()
        modmenu.addAction(crop)
        modmenu.addSeparator()
        modmenu.addAction(lin)
        modmenu.addAction(shift)
        modmenu.addSeparator()
        modmenu.addAction(velocity)
        modmenu.addAction(dop)
        modmenu.addSeparator()
        modmenu.addAction(arith)

        regionmenu = bar.addMenu('Measure')

        defregion=QtWidgets.QAction("Define Region", self.parent())
        defregion.setShortcut('x')
        defregion.setStatusTip('Define a measurement region.')
        defregion.triggered.connect(self.parent().regionload)

        regionclear=QtWidgets.QAction("Clear Region", self.parent())
        regionclear.setShortcut('ctrl+z')
        regionclear.setStatusTip('Clears the region, without resizing the view.')
        regionclear.triggered.connect(self.parent().region_clear)

        equivalent=QtWidgets.QAction("Equivalent Width", self.parent())
        equivalent.setShortcut('e')
        equivalent.setStatusTip('Measure equivalent width, using two points on either side of spectrum.')
        equivalent.triggered.connect(self.parent().eqw)

        snr=QtWidgets.QAction("Signal To Noise", self.parent())
        snr.setShortcut('n')
        snr.setStatusTip('Mean/standard deviation for selected region.')
        snr.triggered.connect(self.parent().signal2noise)

        fitgaus=QtWidgets.QAction("Gaussian Fit", self.parent())
        fitgaus.setShortcut('g')
        fitgaus.setStatusTip('Fit Gaussian between two points.')
        fitgaus.triggered.connect(partial(self.parent().fit,func="gauss"))

        fitvoigt=QtWidgets.QAction("Voigt Fit", self.parent())
        fitvoigt.setShortcut('v')
        fitvoigt.setStatusTip('Fit Voigt between two points.')
        fitvoigt.triggered.connect(partial(self.parent().fit,func="voigt"))

        fitlorentz=QtWidgets.QAction("Lorentzian Fit", self.parent())
        fitlorentz.setShortcut('l')
        fitlorentz.setStatusTip('Fit Lorentzian between two points.')
        fitlorentz.triggered.connect(partial(self.parent().fit,func="lorentz"))

        fitmoffat=QtWidgets.QAction("Moffat Fit", self.parent())
        fitmoffat.setShortcut('m')
        fitmoffat.setStatusTip('Fit Moffat function between two points.')
        fitmoffat.triggered.connect(partial(self.parent().fit,func="moffat"))

        saveregion=QtWidgets.QAction("Save EQW/Fit Region", self.parent())
        saveregion.setShortcut('Shift+Ctrl+x')
        saveregion.triggered.connect(self.parent().SaveRegion)

        loadregion=QtWidgets.QAction("Load EQW/Fit Region", self.parent())
        loadregion.setShortcut('Alt+Ctrl+x')
        loadregion.triggered.connect(self.parent().LoadRegion)

        bisect=QtWidgets.QAction("Bisect Feature",self.parent())
        bisect.setShortcut('w')
        bisect.setStatusTip('Bisect a feature by selecting a range of points on the left and right.')
        bisect.triggered.connect(self.parent().BisectLine)

        bisectsave=QtWidgets.QAction("Save Bisection Regions", self.parent())
        bisectsave.setShortcut('Shift+Ctrl+w')
        bisectsave.triggered.connect(self.parent().SaveBisect)

        bisectload=QtWidgets.QAction("Load Bisection Regions", self.parent())
        bisectload.setShortcut('Alt+Ctrl+w')
        bisectload.triggered.connect(self.parent().LoadBisect)

        coords=QtWidgets.QAction("Print Coordinates of Click", self.parent())
        coords.setShortcut('space')
        coords.triggered.connect(self.parent().coord)

        vert=QtWidgets.QAction("Set Height", self.parent())
        vert.setShortcut('1')
        vert.triggered.connect(self.parent().click_height)

        regionmenu.addAction(defregion)
        regionmenu.addAction(regionclear)
        regionmenu.addSeparator()
        regionmenu.addAction(equivalent)
        regionmenu.addAction(snr)
        regionmenu.addAction(fitgaus)
        regionmenu.addAction(fitvoigt)
        regionmenu.addAction(fitlorentz)
        regionmenu.addAction(fitmoffat)
        regionmenu.addAction(crop)
        regionmenu.addAction(saveregion)
        regionmenu.addAction(loadregion)
        regionmenu.addSeparator()
        regionmenu.addAction(bisect)
        regionmenu.addAction(bisectsave)
        regionmenu.addAction(bisectload)
        regionmenu.addSeparator()
        regionmenu.addAction(coords)
        regionmenu.addAction(vert)


        downstack=QtWidgets.QAction("Display Next in Stack Down",self.parent())
        downstack.setShortcut(',')
        downstack.triggered.connect(self.parent().stackdown)

        upstack=QtWidgets.QAction("Display Next in Stack Up",self.parent())
        upstack.setShortcut('.')
        upstack.triggered.connect(self.parent().stackup)

        stackgoto=QtWidgets.QAction("Plot Spectra by Number or Range of Numbers",self.parent())
        stackgoto.setShortcut('/')
        stackgoto.triggered.connect(self.parent().gotostack)

        savestack=QtWidgets.QAction("Save Stack List to File",self.parent())
        savestack.triggered.connect(self.parent().savestack)

        stacksave=QtWidgets.QAction("Save Stack as New Fits Files",self.parent())
        stacksave.setShortcut('Alt+s')
        stacksave.triggered.connect(self.parent().stacksave)

        stacksave_txt=QtWidgets.QAction("Save Stack as New text Files",self.parent())
        #stacksave_txt.setShortcut('Alt+s')
        stacksave_txt.triggered.connect(self.parent().stacksave_txt)

        stacknorm=QtWidgets.QAction("Normalize",self.parent())
        stacknorm.setShortcut('Alt+t')
        stacknorm.triggered.connect(partial(self.parent().stacker,func="norm"))

        stackeqw=QtWidgets.QAction("Equivalent Width",self.parent())
        stackeqw.setShortcut('Alt+e')
        stackeqw.triggered.connect(partial(self.parent().stacker,func="eqw"))

        stackgaus=QtWidgets.QAction("Gaussian Fit",self.parent())
        stackgaus.setShortcut('Alt+g')
        stackgaus.triggered.connect(partial(self.parent().stacker,func="gauss"))

        stackvoigt=QtWidgets.QAction("Voigt Fit",self.parent())
        stackvoigt.setShortcut('Alt+v')
        stackvoigt.triggered.connect(partial(self.parent().stacker,func="voigt"))

        stacklorentz=QtWidgets.QAction("Lorentzian Fit",self.parent())
        stacklorentz.setShortcut('Alt+l')
        stacklorentz.triggered.connect(partial(self.parent().stacker,func="lorentz"))

        stackmoffat=QtWidgets.QAction("Moffat Fit",self.parent())
        stackmoffat.setShortcut('Alt+m')
        stackmoffat.triggered.connect(partial(self.parent().stacker,func="moffat"))

        stackcrop=QtWidgets.QAction("Crop",self.parent())
        stackcrop.setShortcut('Alt+c')
        stackcrop.triggered.connect(partial(self.parent().stacker,func="scopy"))

        stackbisect=QtWidgets.QAction("Bisect",self.parent())
        stackbisect.setShortcut('Alt+w')
        stackbisect.triggered.connect(partial(self.parent().stacker,func="bisect"))

        stackalign=QtWidgets.QAction("Align",self.parent())
        stackalign.setShortcut('Alt+a')
        stackalign.triggered.connect(partial(self.parent().stacker,func="align"))

        stackrestore=QtWidgets.QAction("Reset to Original Spectra",self.parent())
        stackrestore.setShortcut('Alt+o')
        stackrestore.triggered.connect(partial(self.parent().stacker,func="restoreall"))

        stacksnr=QtWidgets.QAction("Singnal to Noise",self.parent())
        stacksnr.setShortcut('Alt+n')
        stacksnr.triggered.connect(partial(self.parent().stacker,func="snr"))

        dstack=QtWidgets.QAction("Delete From Stack",self.parent())
        dstack.setShortcut('backspace')
        dstack.triggered.connect(self.parent().delfromstack)

        coordstack=QtWidgets.QAction("Print Coordinates of Click",self.parent())
        coordstack.setShortcut('Alt+space')
        coordstack.triggered.connect(partial(self.parent().stacker,func="coord"))

        boxcarstack=QtWidgets.QAction("Boxcar Smooth",self.parent())
        boxcarstack.setShortcut('Alt+Ctrl+b')
        boxcarstack.triggered.connect(partial(self.parent().stacker,func="boxcar"))

        gaussstack=QtWidgets.QAction("Gaussian Smooth",self.parent())
        gaussstack.setShortcut('Alt+Ctrl+g')
        gaussstack.triggered.connect(partial(self.parent().stacker,func="gsmooth"))

        savepng=QtWidgets.QAction("Save Stack to Images",self.parent())
        savepng.triggered.connect(partial(self.parent().stacker,func="savepng"))

        stackview = bar.addMenu('Stack')
        stackview.addAction(stackplottoggle)
        stackview.addAction(downstack)
        stackview.addAction(upstack)
        stackview.addAction(stackgoto)
        stackview.addSeparator()
        stackview.addAction(dstack)
        stackview.addAction(savestack)
        stackview.addSeparator()
        stackview.addAction(stacksave)
        stackview.addAction(stacksave_txt)
        stackview.addAction(savepng)
        stackview.addSeparator()
        stackview.addAction(stackrestore)

        stackmod = bar.addMenu('Stack Modify')
        stackmod.addAction(stacknorm)
        stackmod.addAction(stackalign)
        stackmod.addAction(stackcrop)
        stackmod.addSeparator()
        stackmod.addAction(boxcarstack)
        stackmod.addAction(gaussstack)
        stackmod.addSeparator()
        stackmod.addAction(arith)


        stackmeas = bar.addMenu('Stack Measure')
        stackmeas.addAction(defregion)
        stackmeas.addAction(regionclear)
        stackmeas.addSeparator()
        stackmeas.addAction(stackeqw)
        stackmeas.addAction(stackgaus)
        stackmeas.addAction(stackvoigt)
        stackmeas.addAction(stacklorentz)
        stackmeas.addAction(stackmoffat)
        stackmeas.addAction(stackbisect)
        stackmeas.addAction(stacksnr)
        stackmeas.addSeparator()
        stackmeas.addAction(coordstack)

        aboutmenu = bar.addMenu('About')
        about=QtWidgets.QAction('About', self.parent())
        about.setShortcut('Ctrl+i')
        about.setStatusTip('Click this to find out about PySplot, Nya')
        about.triggered.connect(self.parent().About)
        aboutmenu.addAction(about)
