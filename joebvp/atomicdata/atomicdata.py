from __future__ import print_function, absolute_import, division, unicode_literals

from .. import joebgoodies as jbg
import numpy as np
from linetools.spectralline import AbsLine
from linetools.lists import parse as lilp
from linetools.lists.linelist import LineList
import astropy.units as u
import imp

ilist = LineList('ISM')

jbvp_path = imp.find_module('joebvp')[1]

vernerlist=np.genfromtxt(jbvp_path+'/atomicdata/verner6.txt',dtype=None,delimiter=[10,8,3,4,3,2,9,6])
vernlam=jbg.arrfromcol(vernerlist,0)
vernion=jbg.arrfromcol(vernerlist,1)
vernzatom=jbg.arrfromcol(vernerlist,2)
vernnume=jbg.arrfromcol(vernerlist,3)
verngl=jbg.arrfromcol(vernerlist,4)
verngu=jbg.arrfromcol(vernerlist,5)
vernosc=jbg.arrfromcol(vernerlist,6)
vernp=jbg.arrfromcol(vernerlist,7)

# A start to using the linetools atomic data framework
adata=lilp.parse_morton03()
vdata=lilp.parse_verner96()


for i in range(len(vernion)):
    vernion[i]=vernion[i].strip()

def setatomicdata(lines,precise=True,linelist=ilist):
	lam=np.zeros(len(lines)) ; fosc=np.zeros(len(lines)) ; gam=np.zeros(len(lines))
	for i,ll in enumerate(lines):
		try:
			al=AbsLine(ll*u.AA,closest=True,linelist=linelist)
			lam[i]=al.data['wrest'].value ; fosc[i]=al.data['f'] ; gam[i]=al.data['gamma'].value
		except:
			idx=jbg.closest(adata['wrest'],ll)
			lam[i]=adata['wrest'][idx] ; fosc[i]=adata['f'][idx] ; gam[i]=adata['gamma'][idx]
		if ((abs(lam[i]-ll)>0.01)&(precise==True)):
			idx = jbg.closest(vdata['wrest'], ll)
			try:
				lam[i] = vdata['wrest'][idx].value; fosc[i] = vdata['f'][idx]; gam[i] = vdata['gamma'][idx].value
			except:
				lam[i] = vdata['wrest'][idx]; fosc[i] = vdata['f'][idx]; gam[i] = vdata['gamma'][idx]
	return lam,fosc,gam

def closestlam(restwave):
    lam,fosc,gam = setatomicdata(restwave)

    return lam[0]

def lam2ion(restwave):
    if (isinstance(restwave,int))|(isinstance(restwave,float)):
        ionstr = vernion[jbg.closest(vernlam,restwave)].strip()
        import pdb; pdb.set_trace()
        return ionstr
    else: 
        ions=[]
        for rw in restwave: ions.append(vernion[jbg.closest(vernlam,rw)].strip())
        return ions

def lam2osc(restwave):

    return round(vernosc[jbg.closest(vernlam,restwave)],3)

def lam2vernp(restwave):
    if (isinstance(restwave,int)|isinstance(restwave,float)):
        return round(vernp[jbg.closest(vernlam,restwave)],2)
    else:
        return vernp[jbg.closest(vernlam,restwave)].round(2)

def ion2lam(ion):
    return vernlam[np.where(vernion==ion)[0]]

def ion2laminrange(ion,wave1,wave2,z=0.,frame='obs',pthresh=9.5):
    '''
    Usage: ion2laminrange(ion,wave1,wave2,z=0.,frame='obs',pthresh=9.5)
    Given ion (string) and wavelength range observed, return lambdas to search.
    Use frame='rest' to return restframe wavelengths even if supplying a redshift.
    '''
    try: lamidx=np.where(vernion==ion)[0]
    except: 
        splitname=ion[0:2]+' '+ion[2:]
        try: lamidx=np.where(vernion==ion)[0]
        except:
            print('Ion name doesn\'t match database.')
        else:
            restlams=vernlam[lamidx]
            obswaves=restlams*(1.+z)
            linesinrange=lamidx[np.where(jbg.between(obswaves,wave1,wave2)&(lam2vernp(restlams)>pthresh))]
            if frame=='obs':
               return vernlam[linesinrange]*(1.+z)
            if frame=='rest':
               return vernlam[linesinrange]
    else:
        restlams=vernlam[lamidx]
        obswaves=restlams*(1.+z)
        linesinrange=lamidx[np.where(jbg.between(obswaves,wave1,wave2)&(lam2vernp(restlams)>pthresh))]
        if frame=='obs':
           return vernlam[linesinrange]*(1.+z)
        if frame=='rest':
           return vernlam[linesinrange]
