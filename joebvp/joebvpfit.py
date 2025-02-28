### Fit a voigt profile to a normalized spectrum
### Parameters: [lam_rest,column_density,doppler_param,z,velocity_offset]

import numpy as np
import joebvp.joebgoodies as jbg

# from stsci.tools import nmpfit
from joebvp import nmpfit
from joebvp import utils
from joebvp import makevoigt
try:
	import joebvp_cfg as cfg
except:
	print("joebvp.joebvpfit: No local joebvp_cfg.py found, using default cfg.py file form joebvp.")
	from joebvp import cfg
import joebvp.atomicdata as atomicdata
from astropy.io import ascii
from astropy import units as u
from numpy import random
import warnings
import sys
from linetools import utils as ltu
from linetools.spectra.xspectrum1d import XSpectrum1D

c=299792.458

def foldpars(pars,numpars=5):
	rows=len(pars)/numpars
	fpars=[]
	for i in np.arange(numpars,dtype='int'):
		fparrow=[]
		for j in np.arange(rows,dtype='int'):
			fparrow.append(pars[i+numpars*j])
		fpars.append(fparrow)
	return fpars

def unfoldpars(pars,numpars=5):
	rows=len(pars[0])
	ufpars=[]
	for j in range(rows):
		for i in np.arange(numpars,dtype='int'):
			ufpars.append(pars[i][j])
	return ufpars

def voigtfunc(vwave,vpars,fitcfg=cfg):
	### Check to see if cfg variables are set
	if isinstance(fitcfg.fitidx, int)|isinstance(fitcfg.wave, int):
		fitcfg.fitidx = fitpix(vwave, vpars)
		fitcfg.wave = vwave
	if len(fitcfg.lsfs) == 0:
		makevoigt.get_lsfs(fitcfg)
	vflux=np.zeros(len(vwave))+1.
	factor=makevoigt.voigt(vwave,vpars[0],vpars[1],vpars[2],vpars[3],vpars[4])
	convfactor=makevoigt.convolveprof(vwave,factor,fitcfg)
	vflux*=convfactor
	return vflux

def voigterrfunc(p,x,y,err,fjac=None):
	fp=foldpars(p)
	model=voigtfunc(x,fp)
	status=0
	return([status, (y[cfg.fitidx] - model[cfg.fitidx]) / err[cfg.fitidx]])

def update_bad_pixels(fitcfg=cfg):
	# define bad pixels
	try:
		cond_badpix = (fitcfg.spectrum.wavelength <= fitcfg.spectrum.wvmin) | \
					(fitcfg.spectrum.wavelength >= fitcfg.spectrum.wvmax) | \
					(fitcfg.spectrum.sig <= 0)
					#(cfg.spectrum.flux / cfg.spectrum.sig < cfg.min_sn)  # bad S/N
	except:
		pass
	# spectral gaps
	for gap in fitcfg.spectral_gaps:
		cond_gap = (fitcfg.spectrum.wavelength >= gap[0]*u.AA) & (fitcfg.spectrum.wavelength <= gap[1]*u.AA)
		cond_badpix = cond_badpix | cond_gap
	bad_pixels = np.where(cond_badpix)[0]
	return bad_pixels

def fitpix(wave,pararr,find_bad_pixels=True,fitcfg=cfg):
	if find_bad_pixels:
		# define bad pixels
		fitcfg.bad_pixels = update_bad_pixels(fitcfg) # this variable stores the indices of bad pixels
	else:
		fitcfg.bad_pixels = []
	ll=pararr[0]
	lz=pararr[3]
	lv1=pararr[5]
	lv2=pararr[6]
	relpix=[]
	for i in range(len(ll)):
		vels=jbg.veltrans(lz[i],wave,ll[i])
		p1=jbg.closest(vels,lv1[i])
		p2=jbg.closest(vels,lv2[i])

		if ((p1>=10) & (p2<=(len(wave)-1-10))):
			relpix.extend(range(p1-10,p2+10))
		elif (p1<10):
			relpix.extend(range(0, p2 + 10))
		else:
			relpix.extend(range(p1 - 10, len(wave)-1))
	rp = np.unique(np.array(relpix))
	clean_rp = np.array([i for i in rp if i not in fitcfg.bad_pixels])
	return clean_rp

def prepparinfo(linepars,parflags):
	parinfo=[]
	parflags=np.array(parflags)
	numpars=5
	for i in range(len(parflags[0])):
		parinfo.extend([{'fixed':1},{'fixed':0},{'fixed':0},{'fixed':1},{'fixed':0}])
		for j in range(1,len(parflags)):
			if parflags[j][i]==1: parinfo[i*numpars+j]['fixed']=1
			elif parflags[j][i]<=0: parinfo[i*numpars+j]['fixed']=0
			else:
				matches=np.where(np.array(parflags[j])==parflags[j][i])[0]
				if matches[0]!=i:
					tiedpar=int(matches[0]*numpars+j)
					parinfo[i*numpars+j]['tied']='p['+str(tiedpar)+']'
		col=round(linepars[1][i],2)
		vel=round(linepars[4][i],2)
		bpar=round(linepars[2][i],2)
		parinfo[numpars*i+1]['limited']=[1,1]
		parinfo[numpars*i+1]['limits']=[round(col-5.,2),round(col+5.,2)]
		parinfo[numpars*i+2]['limited']=[1,1]
		### adjust b-value limits to allow for broadened HI features
		lydiff=abs(linepars[0][i] - cfg.lyseries)
		lymatch = np.where(abs(lydiff)<=0.05)[0]
		if lymatch:
			parinfo[numpars*i+2]['limits']=[max([cfg.lowblim, bpar - 10.]), min([bpar + 10, cfg.upperblim_HI])]
		else:
			parinfo[numpars*i+2]['limits']=[max([cfg.lowblim, bpar - 10.]), min([bpar + 10, cfg.upperblim])]
		parinfo[numpars*i+2]['step']=0.5
		parinfo[numpars*i+2]['mpside']=2
		#parinfo[numpars*i+2]['relstep']=0.0001
		parinfo[numpars*i+4]['limited']=[1,1]
		### Allow velocity to flop around
		if parflags[4][i]<0:
			flopamt=abs(parflags[4][i])
			parinfo[numpars*i+4]['limits']=[round(vel-flopamt,2),round(vel+flopamt,2)]
		elif len(linepars)>5:
			v1=round(linepars[5][i],2) ; v2=round(linepars[6][i],2)
			parinfo[numpars*i+4]['limits']=[v1,v2]
		else:
			parinfo[numpars*i+4]['limits']=[round(vel-50.,2),round(vel+50.,2)]
		parinfo[numpars*i+4]['step']=1.
		#parinfo[numpars*i+4]['relstep']=0.01
	return parinfo

def joebvpfit(wave,flux,sig,linepars,flags):

	xtol=1e-11
	gtol=1e-11
	# Only feed to the fitter the parameters that go into the model
	partofit=linepars[:5]
	parinfo=prepparinfo(partofit,flags)
	# Save the velocity windows to add back to the parameter array
	vlim1=linepars[5] ; vlim2=linepars[6]
	# Get atomic data
	lam,fosc,gam=atomicdata.setatomicdata(linepars[0])
	cfg.lams=lam ; cfg.fosc=fosc ; cfg.gam=gam
	# Set fit regions
	cfg.fitidx = fitpix(wave, linepars)
	# Prep parameters for fitter
	partofit=unfoldpars(partofit)
	modelvars={'x':wave,'y':flux,'err':sig}
	# Do the fit and translate the parameters back into the received format
	m=nmpfit.mpfit(voigterrfunc,partofit,functkw=modelvars,parinfo=parinfo,nprint=1,quiet=0,fastnorm=1,ftol=1e-10,xtol=xtol,gtol=gtol)
	if m.status <= 0: print('Fitting error:',m.errmsg)
	fitpars=foldpars(m.params)
	fiterrors = foldpars(m.perror)
	# Add velocity windows back to parameter array
	fitpars.append(vlim1) ; fitpars.append(vlim2)


	print('\nFit results: \n')
	for i in range(len(fitpars[0])):
		print(jbg.tabdelimrow([round(fitpars[0][i],2),jbg.decimalplaces(fitpars[3][i],5),jbg.roundto(fitpars[1][i],5),jbg.roundto(fitpars[2][i],5),jbg.roundto(fitpars[4][i],5)])[:-2])
		print(jbg.tabdelimrow([' ',' ',' ',round(fiterrors[1][i],3),round(fiterrors[2][i],3),round(fiterrors[4][i],3)]))

	return fitpars,fiterrors


def initlinepars(zs,restwaves,initvals=[],initinfo=[]):
	'''

	Parameters
	----------
	zs : numpy vector of floats
		Redshifts of lines (this parameter will be fixed during fitting)
	restwaves : numpy vector of floats
		Rest frame wavelengths of lines to be fitted
	initvals : list of numpy vectors, optional
		Contains the following (in order): [restwaves,linecol,lineb,zs,linevel,linevlim1,linevlim2]
		Will default to values set in cfg.py if not set
	initinfo : list of numpy vectors, optional
		Contains the flags for fitting (in order): [colflag,bflag,velflag]
		Parameters with flags = 0 and 1 will freely value and fixed, respectively
		If 2 or more lines have the same flag value for the same parameter, the parameters will
			be tied to one another.

	Returns
	-------
	initpars : list of lists
		Parameters for fit ready for fitter!
	parinfo : array of arrays
		Flags to be used in fit
	'''

	### Set atomic data for each line
	lam,fosc,gam=atomicdata.setatomicdata(restwaves)
	cfg.lams=lam ; cfg.fosc=fosc ; cfg.gam=gam

	initpars=[[],[],[],[],[],[],[]]
	defaultcol=cfg.defaultcol
	defaultb=cfg.defaultb
	if initvals==[]:
		for i in range(len(restwaves)):
			initpars[0].extend([restwaves[i]])
			initpars[1].extend([defaultcol])
			initpars[2].extend([defaultb])
			initpars[3].extend([zs[i]])
			initpars[4].extend([0.])
			initpars[5].extend([-cfg.defaultvlim])
			initpars[6].extend([cfg.defaultvlim])
	else:
		if len(initvals)==5:
			for i in range(len(restwaves)):
				initpars[0].extend([initvals[0][i]])
				initpars[1].extend([initvals[1][i]])
				initpars[2].extend([initvals[2][i]])
				initpars[3].extend([initvals[3][i]])
				initpars[4].extend([initvals[4][i]])
				initpars[5].extend([-cfg.defaultvlim])
				initpars[6].extend([cfg.defaultvlim])
		else:
			initpars=[[],[],[],[],[],[],[]]
			for i in range(len(restwaves)):
				initpars[0].extend([initvals[0][i]])
				initpars[1].extend([initvals[1][i]])
				initpars[2].extend([initvals[2][i]])
				initpars[3].extend([initvals[3][i]])
				initpars[4].extend([initvals[4][i]])
				initpars[5].extend([initvals[5][i]])
				initpars[6].extend([initvals[6][i]])

	### If hard limits on Doppler b-value are smaller or greater than cfg.lowblim or cfg.upperblim,
	### modify those limits
	maxb=np.max(initpars[2][:])
	minb=np.min(initpars[2][:])
	if maxb>cfg.upperblim:
		cfg.upperblim= maxb + 10.
	if minb<cfg.lowblim: cfg.lowblim= minb - 2.

	parinfo=np.zeros([5,len(restwaves)],dtype=int)
	parinfo[0]=parinfo[0]+1
	parinfo[3]=parinfo[3]+1
	if ((initinfo==[])&(initvals==[])):

		### Look for multiplet membership of each line
		seriesassoc = np.zeros(len(restwaves)) - 99
		for i in range(len(restwaves)):
			for j in range(len(cfg.multiplets)):
				currmult = np.array(cfg.multiplets[j])
				if (abs(restwaves[i] - currmult[jbg.closest(currmult, restwaves[i])]) < 0.01):
					seriesassoc[i] = j

		uqions=np.unique(seriesassoc).tolist()

		if -99 in uqions: uqions.remove(-99)
		flagctr=2
		for uqion in uqions:
			ionmatch=np.where(seriesassoc == uqion)[0]
			uqzs=np.unique(zs[ionmatch])
			for uz in uqzs:
				matchcrit=(zs==uz)&(seriesassoc==uqion)
				rellines=np.where(matchcrit)[0]
				uqlams=np.unique(restwaves[rellines])
				if len(uqlams)>1:
					complist=[]
					numcomps=[]
					for ul in uqlams:
						matchcrit2=matchcrit&(restwaves==ul)
						matches=np.where(matchcrit2)[0]
						complist.append(matches)
						numcomps.append(len(matches))
					numcomps=np.array(numcomps)
					complist=np.array(complist)
					compidxsort=sorted(range(len(numcomps)),key = lambda x: numcomps[x],reverse=True)
					numcompsort=numcomps[compidxsort]
					complistsort=complist[compidxsort]
					maxcomps=numcompsort[0]
					for compidx in range(maxcomps):
						for li in range(len(complistsort)):
							if compidx<numcompsort[li]:
								parinfo[1][complistsort[li][compidx]]=flagctr
								parinfo[2][complistsort[li][compidx]]=flagctr
								parinfo[4][complistsort[li][compidx]]=flagctr
							else: continue
						flagctr+=1

	elif initinfo!=[]:
		parinfo[1]=initinfo[0]
		parinfo[2]=initinfo[1]
		parinfo[4]=initinfo[2]
	elif ((initinfo==[])&(initvals!=[])):

		### Look for multiplet membership of each line
		seriesassoc = np.zeros(len(restwaves)) - 99
		for i in range(len(restwaves)):
			for j in range(len(cfg.multiplets)):
				currmult = np.array(cfg.multiplets[j])
				if (abs(restwaves[i] - currmult[jbg.closest(currmult, restwaves[i])]) < 0.01):
					seriesassoc[i] = j

		### Fix measurements that are imported
		for i in range(len(restwaves)):
			if ((initpars[1][i]!=defaultcol)&(initpars[2][i]!=defaultb)):
				parinfo[1][i]=1 ; parinfo[2][i]=1 ; parinfo[4][i]=1
		uqions=np.unique(seriesassoc).tolist()
		if -99 in uqions: uqions.remove(-99)
		flagctr=2
		for uqion in uqions:
			ionmatch=np.where(seriesassoc == uqion)[0]
			uqzs=np.unique(zs[ionmatch])
			for uz in uqzs:
				matchcrit=(zs==uz)&(seriesassoc==uqion)
				rellines=np.where(matchcrit)[0]
				uqlams=np.unique(restwaves[rellines])
				if len(uqlams)>1:
					complist=[]
					numcomps=[]
					for ul in uqlams:
						matchcrit2=matchcrit&(restwaves==ul)
						matches=np.where(matchcrit2)[0]
						complist.append(matches.tolist())
						numcomps.append(len(matches))
					numcomps=np.array(numcomps)
					complist=np.array(complist)
					compidxsort=sorted(range(len(numcomps)),key = lambda x: numcomps[x])
					complistsort=complist[compidxsort]
					complistsort=complistsort.tolist()
					for i in range(len(complistsort)-1):
						for idx in complistsort[i]:
							parinfo[1][idx]=flagctr
							parinfo[2][idx]=flagctr
							parinfo[4][idx]=flagctr
							for j in range(i+1,len(complistsort)):
								idxidx=jbg.closest(initvals[4][complistsort[j]],initvals[4][idx])
								clidx=complistsort[j][idxidx]
								parinfo[1][clidx]=flagctr
								parinfo[2][clidx]=flagctr
								parinfo[4][clidx]=flagctr
								complistsort[j].remove(clidx)
							flagctr+=1

	return initpars,parinfo

def errors_mc(fitpars,fiterrors,parinfo,wave,flux,err,sig=6,numiter=2000):
	cfg.wave=wave
	cfg.flux=flux
	cfg.sigup=err
	numiter=numiter
	colflags=parinfo[1]
	bflags=parinfo[2]
	velflags=parinfo[4]
	uqcfs=np.unique(colflags)
	uqbfs=np.unique(bflags)
	uqvfs=np.unique(velflags)
	
	colerrs=fiterrors[1]
	berrs = fiterrors[2]
	velerrs = fiterrors[4]
	
	# Set up arrays to hold initial parameters and widths of random deviations
	mcpararr=np.transpose(np.array([fitpars[1],fitpars[2],fitpars[3]]))
	diffarr=np.transpose(np.array([sig*fiterrors[1],sig*fiterrors[2],sig*fiterrors[4]]))
	# Build array of random numbers (final array will be this shape)
	randarr = random.random_sample(size=(numiter,len(fitpars[0]),3))
	# Deviate parameters
	mcpararr=(mcpararr-diffarr)+2.*diffarr*randarr
	
	# Some lines' parameters are tied to other lines'; also fits sometime top or bottom out b-values
	for i,ff in enumerate(uqcfs):
		theselines=np.where(colflags==ff)[0]
		thiswitherr=np.where((colflags==ff)&(colerrs>0))[0]
		if len(thiswitherr)==0:
			diff = 5.
		else:
			diff = sig*colerrs[thiswitherr]
		mcpararr[:,theselines,0]=(mcpararr[:,theselines,0]-diff)+2.*diff*randarr[:,theselines,0]
	for i,ff in enumerate(uqbfs):
		theselines=np.where(bflags==ff)[0]
		thiswitherr=np.where((bflags==ff)&(berrs>0))[0]
		if len(thiswitherr)==0:
			diff = 5.
		else:
			diff = sig*berrs[thiswitherr]
		mcpararr[:,theselines,1]=(mcpararr[:,theselines,1]-diff)+2.*diff*randarr[:,theselines,1]
	for i,ff in enumerate(uqvfs):
		theselines=np.where(velflags==ff)[0]
		thiswitherr=np.where((velflags==ff)&(velerrs>0))[0]
		if len(thiswitherr)==0:
			diff = 5.
		else:
			diff = sig*velerrs[thiswitherr]
		mcpararr[:,theselines,2]=(mcpararr[:,theselines,2]-diff)+2.*diff*randarr[:,theselines,2]


	chisq=np.zeros(numiter)
	for i in range(numiter):
		mcpars=[fitpars[0],mcpararr[i,:,0],mcpararr[i,:,1],fitpars[3],mcpararr[i,:,2],fitpars[5],fitpars[6]]
		cfg.fitidx=fitpix(wave, mcpars)
		vef= voigterrfunc(unfoldpars(mcpars),wave,flux,err,fjac=None)
		chisq[i]=np.sum(vef[1]**2)
	return mcpararr,chisq

def readpars(filename,wave1=None,wave2=None):
	'''

	Parameters
	----------
	filename : str
		Name of parameter input (or joebvp output) file
		File should at least have following columns:
			specfile|restwave|zsys|col|sigcol|bval|sigbval|vel|sigvel|nflag|bflag|vflag|vlim1|vlim2
	wave1 : float, optional
		Beginning of wavelength range over which to load lines (must be set with wave2)
	wave2 : float, optional
		End of wavelength range over which to load lines (must be set with wave1)

	Returns
	-------
	fitpars : list of lists
		Parameters for fit ready for fitter!
	fiterrors : array of numpy vectors
		Error array for the fitting initialized to '0' for each param
	parinfo : array of arrays
		Flags to be used in fit
	linecmts: list of lists
		Reliability and comment flags, e.g., from igmguesses
	'''

	linelist = ascii.read(filename)
	linerestwave = linelist['restwave'].data
	linez = linelist['zsys'].data
	if (wave1 == None)&(wave2 == None):
		lineshere = np.arange(len(linelist))
	elif ((wave1 == None)|(wave2 == None))|(wave1>=wave2):
		lineshere = np.arange(len(linelist))
		warnings.warn('Note that both \'wave1\' and \'wave2\' must be declared or neither must be. \n Loading all lines in list.')
	else:
		lineobswave = linerestwave * (1. + linez)
		lineshere = np.where((lineobswave > wave1) & (lineobswave < wave2))[0]
	linelist=linelist[lineshere]
	linelist['ions']=atomicdata.lam2ion(linelist['restwave'])

	linelist.sort(['ions','zsys','vel','restwave'])

	linerestwave = linelist['restwave']
	zs = linelist['zsys']
	linecol = linelist['col']
	lineb = linelist['bval']
	linevel = linelist['vel']
	linevlim1 = linelist['vlim1']
	linevlim2 = linelist['vlim2']
	colflag = linelist['nflag']
	bflag = linelist['bflag']
	velflag = linelist['vflag']
	restwaves = linerestwave
	if (('rely' in linelist.colnames)&('comment' in linelist.colnames)):
		pass
	elif ('rely' in linelist.colnames):
		linelist['comment']=['none']*len(linelist)
	else:
		linelist['rely'] = ['-'] * len(linelist)
		linelist['comment'] = ['none'] * len(linelist)

	reliability = linelist['rely']
	comment = linelist['comment']
	initinfo = [colflag, bflag, velflag]
	initpars = [restwaves, linecol, lineb, zs, linevel, linevlim1, linevlim2]
	fitpars, parinfo = initlinepars(zs, restwaves, initpars, initinfo=initinfo)
	fiterrors = np.zeros([5, len(fitpars[0])])  # Initialize errors to zero
	linecmts = [reliability,comment]
	#fiterrors[1] = colsig
	#fiterrors[2] = bsig
	#fiterrors[4] = velsig
	return fitpars,fiterrors,parinfo,linecmts


def writelinepars(fitpars,fiterrors,parinfo, specfile, outfilename, linecmts=None):
	'''
	Write fit parameters out to file.

	Parameters
	----------
	fitpars : list of lists
		Parameters for fit ready for fitter!
	fiterrors : array of numpy vectors
		Error array for the fitting initialized to '0' for each param
	parinfo : array of arrays
		Flags to be used in fit
	specfile : str
		Name of the input file containing the spectrum
	outfilename : str
		Parameter output filename
	linecmts : list of lists, optional
		Reliability flags and comments, e.g., from igmguesses

	'''
	import os
	### Set outputs and open files
	bigfiletowrite = cfg.largeVPparfile
	filetowrite = outfilename
	if os.path.isfile(filetowrite):
		VPparfile = open(filetowrite, 'wb')
		bigparfile = open(bigfiletowrite, 'ab') # Append to the running list
	else:
		VPparfile = open(filetowrite, 'wb')
		bigparfile = open(bigfiletowrite, 'wb')

	### Prep header of line parameter file
	if linecmts is not None:
		header = b'specfile|restwave|zsys|col|sigcol|bval|sigbval|vel|sigvel|nflag|bflag|vflag|vlim1|vlim2|wobs1|wobs2|pix1|pix2|z_comp|trans|rely|comment \n'
	else:
		header = b'specfile|restwave|zsys|col|sigcol|bval|sigbval|vel|sigvel|nflag|bflag|vflag|vlim1|vlim2|wobs1|wobs2|pix1|pix2|z_comp|trans \n'
	VPparfile.write(header)
	bigparfile.write(header)

	### Grab parameters/info for each line
	for i in range(len(fitpars[0])):
		zline = fitpars[3][i]
		vlim1 = fitpars[5][i]
		vlim2 = fitpars[6][i]
		restwave = fitpars[0][i]
		wobs1 = restwave * (1 + zline + vlim1 / 299792.458)
		wobs2 = restwave * (1 + zline + vlim2 / 299792.458)
		pix1 = jbg.closest(cfg.wave, wobs1)
		pix2 = jbg.closest(cfg.wave, wobs2)
		trans = atomicdata.lam2ion(fitpars[0][i])
		z_comp = ltu.z_from_dv(fitpars[4][i]*u.km/u.s, zline)
		if linecmts is not None:
			towrite = jbg.pipedelimrow(
				[specfile, restwave, round(zline, 5), round(fitpars[1][i], 3), round(fiterrors[1][i], 3),
				 round(fitpars[2][i], 3), round(fiterrors[2][i], 3), round(fitpars[4][i], 3), round(fiterrors[4][i], 3),
				 parinfo[1][i], parinfo[2][i], parinfo[4][i], vlim1, vlim2, wobs1, wobs2, pix1, pix2,round(z_comp, 5), trans,
				 linecmts[0][i],linecmts[1][i]])
		else:
			towrite = jbg.pipedelimrow(
				[specfile, restwave, round(zline, 5), round(fitpars[1][i], 3), round(fiterrors[1][i], 3),
				 round(fitpars[2][i], 3), round(fiterrors[2][i], 3), round(fitpars[4][i], 3), round(fiterrors[4][i], 3),
				 parinfo[1][i], parinfo[2][i], parinfo[4][i], vlim1, vlim2, wobs1, wobs2, pix1, pix2, round(z_comp, 5),trans])
		VPparfile.write(towrite.encode())
		bigparfile.write(towrite.encode())
	VPparfile.close()
	bigparfile.close()
	print('Line parameters written to:')
	print(filetowrite)


def writeVPmodel(outfile, wave, fitpars, normflux, normsig):
	from astropy.table import Table
	model = voigtfunc(wave, fitpars)
	modeltab = Table([wave, model, normflux, normsig], names=['wavelength', 'model', 'normflux', 'normsig'])
	# modeltab.write(outfile, format='fits', overwrite=True)
	dummycont = np.ones(len(wave))
	spec = XSpectrum1D.from_tuple((modeltab['wavelength'], modeltab['model'], modeltab['normsig'], dummycont))
	spec.write_to_fits(outfile)

	print('Voigt profile model written to:')
	print(outfile)

def writeVPmodelByComp(outdir, spectrum, fitpars):
	import copy
	import os,glob
	if cfg.spectrum == []:
		cfg.spectrum = spectrum
	os.mkdir(outdir)
	linelist = utils.abslines_from_fitpars(fitpars)
	complist = utils.abscomponents_from_abslines(linelist)

	for comp in complist:
		print('\n',comp.name)
		wave,model = modelFromAbsComp(cfg.spectrum, comp)
		spec = copy.deepcopy(cfg.spectrum)
		spec.flux = model
		flist = glob.glob(outdir+'/'+comp.name+'*')
		if len(flist)!=0:
			fname = comp.name+'_'+str(len(flist))+ '_VPmodel.fits'
		else:
			fname = comp.name+'_VPmodel.fits'
		spec.write_to_fits(outdir+'/'+fname)

	print('Voigt profile models written to directory:', outdir)

def modelFromAbsComp(spectrum,abscomp):
	'''
	Generate Voigt profile using
	Parameters
	----------
	spectrum: linetools XSpectrum1D
		Spectrum to provide wavelength input for Voigt profile evaluation
	abscomp: linetools AbsComponent
		Component objects whose AbsLines will be used to generate profile

	Returns
	-------
	spectrum.wavelength: 1D array
		Wavelength vector corresponding to profile
	profile: 1D array
		Model consisting of Voigt profiles across spectrum

	'''
	restwaves = []
	zs = []
	cols = []
	bs = []
	vels = []
	vlim1s = []
	vlim2s = []
	for absline in abscomp._abslines:
		### Grab parameters from AbsLine object
		restwaves.append(absline.wrest.value)
		zs.append(absline.z)
		cols.append(absline.attrib['logN'])
		bs.append(absline.attrib['b'])
		vels.append(np.mean(absline.limits.vlim).value)
		vlim1s.append(absline.limits.vlim[0].value)
		vlim2s.append(absline.limits.vlim[1].value)
	pars = [restwaves,cols,bs,zs,vels,vlim1s,vlim2s]

	### Set atomic data so that the right profiles are evaluated!
	lam,fosc,gam=atomicdata.setatomicdata(restwaves)
	cfg.lams=lam ; cfg.fosc=fosc ; cfg.gam=gam

	profile = voigtfunc(spectrum.wavelength.value, pars)
	### For debugging old feature
	#nonone = np.where(np.abs(profile - 1.) > 0.01)[0]
	#print spectrum.wavelength[nonone]
	return spectrum.wavelength.value,profile


def fit_to_convergence(wave,flux,sig,linepars,parinfo,maxiter=50,itertol=0.0001):
	'''

	Parameters
	----------
	wave
	flux
	sig
	linepars
	parinfo

	maxiter : int
		Maximum number of times to run the fit while striving for convergence

	itertol : float
		Maximum difference in any parameter from one fitting iteration to the next.  Routine will fit again if
		any difference in the measurements exceeds itertol.


	Returns
	-------

	'''
	fitpars = linepars
	oldfitpars = np.zeros([7, len(fitpars[0])]) - 99
	ctr = 0
	okay = 1
	while ((np.max(np.abs(fitpars - oldfitpars)) > itertol) & (ctr < maxiter)):
		ctr += 1

		try:
			oldfitpars = fitpars
			fitpars, fiterrors = joebvpfit(wave, flux, sig,
											 fitpars, parinfo)
			fitpars = np.array(fitpars)
			print('Iteration', ctr, '-')

		except:
			print('Fitting error!')
			print("Unexpected error:", sys.exc_info()[0])
			okay = 0
			raise


			#break

	if okay != 0:
		print('Fit converged after',ctr,'iterations.')
		return fitpars, fiterrors
	else:
		return linepars,fiterrors
