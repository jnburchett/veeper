from __future__ import print_function, absolute_import, division, unicode_literals

import os
import numpy as np
from astropy.constants import c as C

wave=0
flux=0
sigup=0
normwave=0
normflux=0
normsig=0
fitidx=0
fitcoeff=0
fitcovmtx=0
fiterrors=0
origpars=0
origparinfo=0
filename=''
field=''
homedir=os.path.expanduser('~')
spectrum = []  # this will store the spectrum
bad_pixels = []  # this will store bad pixel indices
islocal = False


Todd = False
if Todd:
    instr=['COS','COS','COS','COS','STIS']  # Options in include 'COS', 'STIS', 'Gaussian'
    gratings=['G130M','G160M','G185M','G225M','E230M']
    slits=['NA','NA','NA','NA','0.2x0.2']  # Required when using instr = 'STIS'
    lsfranges=np.array([[1100,1460],[1460,1800],[1800,2100],[2100,2278],[2278,4100]])
    lps=['1','1','1','1','1']  # Lifetime positions; required when using instr = 'COS'
    cen_wave=['1327','1600','1953','2250','2707']
    spectral_gaps = []
    pixel_scales = ['N/A','N/A','N/A','N/A','N/A']  # Required when using instr = 'Gaussian' (in Angstroms)
    fwhms = ['N/A','N/A','N/A','N/A','N/A']  # Required when using instr = 'Gaussian' (in Angstroms)

else:
    # sample parameters for the casual COS user
    instr=['COS','COS']  # Options in include 'COS', 'STIS', 'Gaussian'
    lsfranges=np.array([[1130,1450],[1450,1800]])
    gratings=['G130M','G160M']
    cen_wave=['1291','1611']
    slits=['NA','NA']  # Required when using instr = 'COS'
    lps=['2','2']  # Lifetime positions; required when using instr = 'COS'
    #spectral_gaps = [[0,1162.5], [1198,1201.5], [1213.3, 1217.93], [1299.3,1321.6],[1596,1612.8],[1782,2000]]
    spectral_gaps = []
    pixel_scales = ['N/A','N/A']  # Required when using instr = 'Gaussian' (in Angstroms)
    fwhms = ['N/A','N/A']  # Required when using instr = 'Gaussian' (in Angstroms)
    ### Comment the above and uncomment below for sample Gaussian LSF parameters
    '''instr = ['Gaussian']  # Options in include 'COS', 'STIS', 'Gaussian'
    lsfranges = np.array([[3500, 10000]])
    spectral_gaps = []
    pixel_scales = [0.03]  # Required when using instr = 'Gaussian' (in Angstroms)
    fwhms = [0.1413]  # Required when using instr = 'Gaussian' (in Angstroms)'''

# fundamental constants
echarge = 4.803204505713468e-10
m_e = 9.10938291e-28
c = C.to('cm/s').value
c2 = 29979245800.0

# store LSFs and FGs
lsfs=[]
fgs=[]
wavegroups=[]
wgidxs=[]
uqwgidxs=[]

outputdir = './'

VPparoutfile = outputdir + field + '_VP.dat'
VPmodeloutfile = outputdir + field + 'VPmodel.fits'
contoutfile = outputdir + 'continua.dat'
largeVPparfile = outputdir + '_VP_log.dat'
defaultcol = 13.1
defaultb = 20.0
defaultvlim = 1000.
lowblim = 1.
upperblim = 85.
upperblim_HI = 210.

# Plotting
ylim = (-0.2, 1.4)
general_fontsize = 10.
xtick_fontsize = 'small'
ytick_fontsize = 'small'
xy_fontsize = 'small'
x_labelpad = 0
y_labelpad = 0
label_ypos = 0.2 * (ylim[0]+ylim[1])
label_fontsize = 10.
residual_markersize = 1
spec_linewidth = 1


# Handy definitions
lyseries=np.array([913.8260,914.0390,914.2860,914.5760,914.919,915.329,915.824,916.429,917.1806,918.1294,919.3514,920.9631,923.1504,926.2257,930.7483,937.8035,949.7431,972.5368,1025.7223,1215.6701])
