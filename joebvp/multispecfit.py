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
from joebvp import nmpfit
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
        
        thiscfg.wave=spec.wavelength.value
        thiscfg.spectrum = spec # need this for defining bad pixels later
        thiscfg.normflux = spec.flux.value/spec.co.value
        thiscfg.normsig=spec.sig.value/spec.co.value
    
        #fitpars,fiterrors=joebvpfit.fit_to_convergence(wave,normflux,normsig,fitpars,parinfo, **kwargs)
    joebvpfit_multi(cfglist,fitpars,parinfo)
    #import pdb; pdb.set_trace()



def joebvpfit_multi(cfglist,linepars,flags):

    xtol=1e-11
    gtol=1e-11
    # Only feed to the fitter the parameters that go into the model
    partofit=linepars[:5]
    parinfo=joebvpfit.prepparinfo(partofit,flags)
    # Prep parameters for fitter
    partofit=joebvpfit.unfoldpars(partofit)
    # Save the velocity windows to add back to the parameter array
    vlim1=linepars[5] ; vlim2=linepars[6]
    # Set up lists of wavelength, etc., arrays
    xs = []
    ys = []
    errs = []
    for cfg in cfglist:
        # Get atomic data
        lam,fosc,gam=atomicdata.setatomicdata(linepars[0])
        cfg.lams=lam ; cfg.fosc=fosc ; cfg.gam=gam
        # Set fit regions
        cfg.fitidx = joebvpfit.fitpix(cfg.wave, linepars,fitcfg=cfg)
        xs.append(cfg.wave)
        ys.append(cfg.normflux)
        errs.append(cfg.normsig)
    #modelvars={'xs':xs,'ys':ys,'errs':errs,'cfglist':cfglist}
    modelvars = {'cfglist':cfglist}
    
    # Do the fit and translate the parameters back into the received format
    m=nmpfit.mpfit(voigterrfunc_multi,partofit,functkw=modelvars,parinfo=parinfo,nprint=1,quiet=0,fastnorm=1,ftol=1e-10,xtol=xtol,gtol=gtol)
    if m.status <= 0: print('Fitting error:',m.errmsg)
    fitpars=joebvpfit.foldpars(m.params)
    fiterrors = joebvpfit.foldpars(m.perror)
    # Add velocity windows back to parameter array
    fitpars.append(vlim1) ; fitpars.append(vlim2)


    print('\nFit results: \n')
    for i in range(len(fitpars[0])):
        print(jbg.tabdelimrow([round(fitpars[0][i],2),jbg.decimalplaces(fitpars[3][i],5),jbg.roundto(fitpars[1][i],5),jbg.roundto(fitpars[2][i],5),jbg.roundto(fitpars[4][i],5)])[:-2])
        print(jbg.tabdelimrow([' ',' ',' ',round(fiterrors[1][i],3),round(fiterrors[2][i],3),round(fiterrors[4][i],3)]))

    return fitpars,fiterrors

def voigterrfunc_multi(p,cfglist,fjac=None):
    fp=joebvpfit.foldpars(p)
    diffs = np.array([],dtype=np.float32)
    for i,cfg in enumerate(cfglist):
        thismodel = joebvpfit.voigtfunc(cfg.wave,fp,cfg)
        y = cfg.normflux
        sig = cfg.normsig
        thisdiff = ((y[cfg.fitidx] - thismodel[cfg.fitidx]) / sig[cfg.fitidx])
        diffs = np.concatenate([diffs,thisdiff])
    status = 0
    return([status, diffs])

