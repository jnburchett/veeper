# -*- coding: utf-8 -*-
"""
Created on Tues Aug 13 14:59:45 2024

@author: burchett
"""
from __future__ import print_function, absolute_import, division, unicode_literals

from PyQt5.uic import loadUiType
from PyQt5 import QtGui
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5 import QtWidgets

from joebvp.atomicdata import atomicdata
import joebvp.joebgoodies as jbg

from joebvp import joebvpfit
from joebvp import utils as jbu
import os
from linetools.spectra.io import readspec
import numpy as np
from astropy.constants import c
try:
    from importlib import reload
except:
    pass
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)

import importlib

modpath=os.path.abspath(os.path.dirname(__file__))
print(os.path.abspath(os.path.dirname(__file__)))
Ui_MainWindow, QMainWindow = loadUiType(modpath+'/mainvpwindow.ui')

#matplotlib.rcParams['font.size']=cfg.general_fontsize

c= c.to('km/s').value

'''try:
    import joebvp_cfg as cfg
except:
    print("joebvp.VPmeasure: No local joebvp_cfg.py found, using default cfg.py file from joebvp.")
    from joebvp import cfg'''

def multispecfit(specfiles,parfile,cfgfiles):
    # Import arbitrary number of spectra and cfg files; append to lists
    cfglist = []
    spectra = []
    for i,cf in enumerate(cfgfiles):
        if cf[-3:]=='.py':
            cf = cf[:-3]
        cfglist.append(importlib.import_module(cf))
        spectra.append(readspec(specfiles[i]))

    # initialize fit parameters
    fitpars, fiterrors, parinfo, linecmts = joebvpfit.readpars(parfile)


    for i,spec in enumerate(spectra):
        thiscfg = cfglist[i]
        thiscfg.lsfs = []
        thiscfg.fgs = []
        thiscfg.wavegroups = []
        thiscfg.wgidxs = []
        thiscfg.uqwgidxs = []

        wave=spec.wavelength.value
        normflux=spec.flux.value/spec.co.value
        normsig=spec.sig.value/spec.co.value
        thiscfg.wave=wave
        thiscfg.spectrum = spec # need this for defining bad pixels later

        #TODO: make joebvpfit_multi that takes lists of wave,flux,sig,cfg
        #TODO: make voigterrfunc_multi (?) that does appropriate evals and convolutions


        #fitpars,fiterrors=joebvpfit.fit_to_convergence(wave,normflux,normsig,fitpars,parinfo, **kwargs)
    import pdb; pdb.set_trace()


