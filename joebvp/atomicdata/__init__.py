from __future__ import print_function, absolute_import, division, unicode_literals

from .. import joebgoodies as jbg
import numpy as np
from linetools.spectralline import AbsLine
from linetools.lists import parse as lilp
import astropy.units as u
from astropy.table import Table
import imp
from linetools.lists.linelist import LineList

llist = LineList('ISM')

jvp = importlib.util.find_spec('joebvp')
jbvp_path = jvp.submodule_search_locations[0]
#jbvp_path = imp.find_module('joebvp')[1]

vernerlist=np.genfromtxt(jbvp_path+'/atomicdata/verner6.txt',dtype=None,delimiter=[10,8,3,4,3,2,9,6])
vernlam=jbg.arrfromcol(vernerlist,0)
vernion=jbg.arrfromcol(vernerlist,1)
vernzatom=jbg.arrfromcol(vernerlist,2)
vernnume=jbg.arrfromcol(vernerlist,3)
verngl=jbg.arrfromcol(vernerlist,4)
verngu=jbg.arrfromcol(vernerlist,5)
vernosc=jbg.arrfromcol(vernerlist,6)
vernp=jbg.arrfromcol(vernerlist,7)

starts = (0,10,18,21,24,27,30,40)
verntab = Table.read(jbvp_path+'/atomicdata/verner6.txt',
                     format='ascii.fixed_width_no_header',
                     col_starts=starts,fill_values=('','0'),
                     names=['lam','ion','Z','nume','gl','gu','fosc','P'])
verntab['fosc'].astype(float)
vernion=verntab['ion']

# A start to using the linetools atomic data framework
adata=lilp.parse_morton03()
vdata=lilp.parse_verner96()

for i in range(len(vernion)):
    vernion[i]=vernion[i].strip()

def setatomicdata(lines,precise=True):

    lam=np.zeros(len(lines)) ; fosc=np.zeros(len(lines)) ; gam=np.zeros(len(lines))
    for i,ll in enumerate(lines):
        try:
            al=AbsLine(ll*u.AA,closest=True,linelist = llist)
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
        ionstr = vernion[jbg.closest(vernlam, restwave)].strip()
        return ionstr
    else:
        ions=[]
        for rw in restwave: ions.append(vernion[jbg.closest(vernlam,rw)].strip())
        return ions

def lam2osc(restwave):
    lam,fosc,gam = setatomicdata(restwave)
    return fosc[0]

def lam2vernp(restwave):
    if (isinstance(restwave,int)|isinstance(restwave,float)):
        return round(vernp[jbg.closest(vernlam,restwave)],2)
    else:
        return vernp[jbg.closest(vernlam,restwave)].round(2)
