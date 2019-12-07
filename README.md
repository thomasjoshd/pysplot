# pysplot
An IRAF-like tool for manipulating, and measuring 1-D FITS spectra.

Just download pysplot.py and run from the python3 commandline: pysplot3 pysplot.py

Requires the following libararies:
functools, numpy >= numpy-1.17.4, os, tkinter, matplotlib, astropy >= astropy-3.2.3 , datetime.

shortcut keys in the menus are denoted with ().

Features:

Open fit, fits, and headerless text spectra

Save modifed spectra as fit, fits, headerless txt.

Display image header.

Display 1-D spectra
  single spectrum mode
  overplot mode
  stack plot mode
  (dynamical-time series planned)

Polynomial continuum normalization.
Boxcar smoothing.
Convert display axis from wavelength to velocity and vice versa.

Measure Equivalent Width.
Fit Gaussian, Voigt, and Lorentzian to spectral features. (currently the display doesn't look good.)
Trim or crop spectra, which can be saved to a new fits file.
A spectral feature bisection based on a pair of points on either side of the profile.

Normalization parameters, regions, and bisection regions can all be saved/loaded.

The stack menu allows the application of any of the above features to an entire stack of spectra. (currently only tested on a stack of approximately uniform wavelength coverage.)

(A csv log output is planned.)
