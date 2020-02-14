#######################################################################
###
###  Main function for generating Voigt profiles.  
###  convenience function: cosvoigt(waves,params)
###  params format (each entry a list): [restwave,coldens,bval,z,vel]
###
#######################################################################
from __future__ import print_function, absolute_import, division, unicode_literals

from PyQt5 import QtGui, QtCore

from joebvp import joebgoodies as jbg
import numpy as np  
from scipy.signal import convolve
from scipy.special import wofz
import scipy.stats
import joebvp.atomicdata as atomicdata
from astropy import constants as const
import sys,os
try:
    import joebvp_cfg as cfg
except:
    print("joebvp.makevoigt: No local joebvp_cfg.py found, using default cfg.py file from joebvp.")
    from joebvp import cfg
from sklearn.cluster import MeanShift, estimate_bandwidth
import astropy.units as u
from astropy.io import ascii
from astropy.table import Table
from linetools.spectra.lsf import LSF

ln2=np.log(2)
c=const.c.value/1e3

modpath=os.path.dirname(__file__)

def Hfunc(x,a):
    z=x+1j*a
    I=wofz(z).real
    return I

def cosvoigt(vwave,vpars):
    from joebvp import joebvpfit
    pars,info=joebvpfit.initlinepars(vpars[3],vpars[0],vpars,initinfo=[0,0,0])
    cfg.fitidx=joebvpfit.fitpix(vwave, pars)
    cfg.wave=vwave
    vflux=np.zeros(len(vwave))+1.
    factor=voigt(vwave,vpars[0],vpars[1],vpars[2],vpars[3],vpars[4])
    convfactor=convolvecos(vwave,factor,vpars[0],vpars[3])
    vflux*=convfactor
    return vflux

def cosvoigt_cont(vwave,cont,vpars):
    from joebvp import joebvpfit
    pars,info=joebvpfit.initlinepars(vpars[3],vpars[0],vpars,initinfo=[0,0,0])
    cfg.fitidx=joebvpfit.fitpix(vwave, pars,find_bad_pixels=False)
    cfg.wave=vwave
    vflux=np.zeros(len(vwave))+1.
    factor=voigt(vwave,vpars[0],vpars[1],vpars[2],vpars[3],vpars[4])
    convfactor=convolvecos(vwave,factor*cont,vpars[0],vpars[3])
    vflux*=convfactor
    return vflux

def voigt(waves,line,coldens,bval,z,vels):
    tautot=np.zeros(len(waves))
    if len(cfg.lams)==0:
        lam,fosc,gam=atomicdata.setatomicdata(line)
        cfg.lams=lam ; cfg.fosc=fosc ; cfg.gam=gam
    for i in range(len(coldens)):
        #thatfactor=(1.+z[i])*(1.+vels[i]/c)
        thatfactor=(1.+z[i])

        lam0=cfg.lams[i]
        gam=cfg.gam[i]
        fosc=cfg.fosc[i]
        lam = waves/thatfactor
        dlam=bval[i]*lam0/c  #Doppler param in wavelength
        x=(lam-lam0-lam0*vels[i]/c)/dlam
        a=gam/(4.*np.pi*(c*1e13/lam0**2*dlam))
        vp=Hfunc(x,a)
        tauval= np.sqrt(np.pi) * cfg.echarge ** 2 / cfg.m_e / cfg.c ** 2 * (lam0 * 1e-8) ** 2 / dlam * 1e8 * (10 ** coldens[i]) * fosc * vp
        tautot+=tauval
    return np.exp(-tautot)

def get_lsfs():

    lsfobjs=[]
    for i,inst in enumerate(cfg.instr):
        if inst in ['COS','STIS']:
            lsfobjs.append(LSF(dict(name=inst, grating=cfg.gratings[i],
                                life_position=cfg.lps[i], cen_wave=cfg.cen_wave[i],
                                slit=cfg.slits[i])))
        elif inst in ['Gaussian','gaussian']:
            lsfobjs.append(LSF(dict(name=inst, pixel_scale=cfg.pixel_scales[i],
                                    FWHM=cfg.fwhms[i])))
    cfg.lsfs=[]
    for fg in cfg.fgs:
        if isinstance(fg,int):
            lamobs=cfg.wave[fg]
            lsfmatch = jbg.wherebetween(lamobs, cfg.lsfranges[:, 0], cfg.lsfranges[:, 1])
            lsf = lsfobjs[lsfmatch[0]].interpolate_to_wv_array(cfg.wave[cfg.fgs] * u.AA, kind='cubic')
            cfg.lsfs.append(lsf['kernel'])
            break
        else:
            lamobs=np.median(cfg.wave[fg])
            lsfmatch = jbg.wherebetween(lamobs, cfg.lsfranges[:, 0], cfg.lsfranges[:, 1])
            lfg = len(fg)
            if len(fg) < 10:
                print("Line at {:.2f} AA is undersampling the LSF. Will increase number of pixels at either side to"
                    " include at least 10.".format(lamobs))
                if np.isnan(lamobs):
                    import pdb;pdb.set_trace()
                    # NT: this could happen when pixels with S/N<0 exist, i.e. pixels where flux is <0 (e.g. black lines).
                    # For this reason, is better to remove the option to eliminate pixels based on S/N.
                n_more_side = int(np.ceil((10 - len(fg))/2 + 1))
                inds_right = [np.max(fg) + ii + 1 for ii in range(n_more_side)]
                inds_left = [np.min(fg) - ii - 1 for ii in range(n_more_side)]
                inds_left.sort()
                fg = inds_left + fg.tolist() + inds_right
                fg = np.array(fg)
                print("New fg is: {}".format(fg))
            try:
                lsf = lsfobjs[lsfmatch[0]].interpolate_to_wv_array(cfg.wave[fg] * u.AA, kind='cubic')
            except:
                import pdb; pdb.set_trace()


            # except:
            # 	QtCore.pyqtRemoveInputHook()
            # 	import pdb; pdb.set_trace()
            # 	QtCore.pyqtRestoreInputHook()

            cfg.lsfs.append(lsf['kernel'])

def convolvecos(wave,profile,lines,zs):
    if len(wave)>len(cfg.fitidx):
        if len(cfg.fitidx) > 0:
            fitwaves=wave[cfg.fitidx]
        else:
            # todo: figure out why are there cfg.fitidx = []?
            #fitwaves = wave
            raise ValueError('No valid pixel ranges for fit, possibly due spectral_gaps setting in cfg.py')
            import pdb; pdb.set_trace()
    else:
        fitwaves=wave
    if cfg.wavegroups==[]:
        X = np.array(list(zip(fitwaves,np.zeros(len(fitwaves)))), dtype=float)
        ms = MeanShift(bandwidth=25.)
        ms.fit(X)
        cfg.wgidxs = ms.labels_
        cfg.wavegroups = ms.cluster_centers_
        cfg.uqwgidxs=np.unique(cfg.wgidxs)
        # Identify groups of fitidxs
        buf=4
        df= cfg.fitidx[1:] - cfg.fitidx[:-1]
        dividers = np.where(df > buf)[0] #These are the last indices of each group
        if len(dividers)==0:
            cfg.fgs=[cfg.fitidx]
        else:
            cfg.fgs=[np.arange(cfg.fitidx[0], cfg.fitidx[dividers[0]])] #1st group
            for i, idx in enumerate(dividers[:-1]):
                newfg = np.arange(cfg.fitidx[idx + 1], cfg.fitidx[dividers[i + 1]]) # 'i+1' b/c 1st group handled separately
                if len(newfg)>1.:
                    cfg.fgs.append(newfg)
            cfg.fgs.append(np.arange(cfg.fitidx[dividers[-1] + 1], cfg.fitidx[-1] + 1))  #last group
            newfgs=[]
            ### Deal with wavegroups that may have huge jumps in wavelength between pixels
            for i,fg in enumerate(cfg.fgs):

                diffs= cfg.wave[fg[1:]] - cfg.wave[fg[:-1]]
                outlier=np.where(diffs>(30.*np.median(diffs)))[0]

                if len(outlier)>0:
                    if len(outlier)==1:
                        cfg.fgs[i]=np.arange(fg[0], fg[outlier + 1])
                        newfgs.append(np.arange(fg[outlier]+1,fg[-1]))
                        if len(newfgs[-1])<15:
                            nummorepix=15-len(newfgs)
                            newfgs[-1]=np.concatenate([newfgs[-1],np.arange(fg[-1],fg[-1]+nummorepix+1)])  #add 1 to 2nd arg to make sure this val is included
                    else:
                        cfg.fgs[i]=np.arange(fg[0], fg[outlier[0]])
                        for j,out in enumerate(outlier):
                            if out==outlier[-1]:
                                newfgs.append(np.arange(fg[out]+1,fg[-1]+1)) #add 1 to 2nd arg to make sure this val is included
                            else:
                                newfgs.append(np.arange(fg[out]+1,fg[outlier[j+1]]+1))   #add 1 to 2nd arg to make sure this val is included
                            if len(newfgs[-1])<15:
                                nummorepix=15-len(newfgs)
                                newfgs[-1]=np.concatenate([newfgs[-1],np.arange(fg[-1],fg[-1]+nummorepix+1)]) #add 1 to 2nd arg to make sure this val is included

            for fg in newfgs:
                cfg.fgs.append(fg)
        get_lsfs()
    convprof=profile
    for i,ll in enumerate(cfg.fgs):
        if isinstance(ll,int):
            lsfwidth= int(np.ceil(len(cfg.fgs) / 2 + 1))
            paddedprof = np.insert(profile[cfg.fgs], 0, [1.] * lsfwidth)
            paddedprof = np.append(paddedprof, [1.] * lsfwidth)
            convprof[cfg.fgs] = convolve(paddedprof, cfg.lsfs[i], mode='same')[lsfwidth:-lsfwidth]
            break
        else:

            lsfwidth=int(np.ceil(len(ll)/2+1))
            paddedprof = np.insert(profile[ll], 0, [1.] * lsfwidth)
            paddedprof=np.append(paddedprof,[1.]*lsfwidth)
            convprof[ll] = convolve(paddedprof, cfg.lsfs[i], mode='same')[lsfwidth:-lsfwidth]


    return convprof