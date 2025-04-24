from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
#from angles import *
#from velcorr import *
import sys
import subprocess
import warnings

### Return index of array element closest to value
def closest(arr,value):
    try:
        iter(value)
        idx=[]
        for val in value: idx.append((np.abs(arr-val)).argmin())
    except TypeError:
        idx = (np.abs(arr-value)).argmin()
    return idx
    '''if (isinstance(value,int)|isinstance(value,float)):
        idx = (np.abs(arr-value)).argmin()
    else:
        idx=[]
        for val in value: idx.append((np.abs(arr-val)).argmin())
    return idx'''

### Transform wavelength into velocity space centered on some line
def veltrans(redshift,waves,line):
    c=299792.458
    if (isinstance(line,int)|isinstance(line,float)):
        transline =c*(waves-(1+redshift)*line)/line/(1+redshift)
    else:
        transline=[]
        for ll in line: transline.append(c*(waves-(1+redshift)*ll)/ll/(1+redshift))
    return transline 

### Return row of pipe-delimited data
def pipedelimrow(data):
    numcols=len(data)
    rowstring=str(data[0])
    for i in range(1,numcols):
        rowstring=rowstring+'|'+str(data[i])
    return rowstring +' \n'

### Quickly load pipe-delimited data file
def loadpipe(filename,names=False):
    if names: tablearr=np.genfromtxt(filename,names=True,dtype=None,delimiter='|')
    else: tablearr=np.genfromtxt(filename,dtype=None,delimiter='|')
    return tablearr

### Return row of comma-delimited data
def commadelimrow(data):
    numcols=len(data)
    rowstring=str(data[0])
    for i in range(1,numcols):
        rowstring=rowstring+','+str(data[i])
    return rowstring +' \n'

### Quickly load comma-delimited data file
def loadcomma(filename,names=False):
    if names: tablearr=np.genfromtxt(filename,names=True,dtype=None,delimiter=',')
    else: tablearr=np.genfromtxt(filename,dtype=None,delimiter=',')
    return tablearr

### Return row of tab-delimited data
def tabdelimrow(data):
    numcols=len(data)
    rowstring=str(data[0])
    for i in range(1,numcols):
        rowstring=rowstring+'\t '+str(data[i])
    return rowstring +' \n'

### Make array from table column
def arrfromcol(datatable,col):
    arr=np.array([row[col] for row in datatable])
    return arr

### Return fraction of L_star for abs. mag.
def llstar(mag,lstar=-21.12):
    return 10.**((lstar-mag)/2.5)



### Convert hhmmss to degrees
def hhmmssdeg(coord):
    if coord[2]=='h':
        hours=float(coord.partition('h')[0])
        mins=float(coord.partition('h')[2].partition('m')[0])
        secs=float(coord.partition('h')[2].partition('m')[2].partition('s')[0])
    elif coord[2]==' ':
        parted=coord.partition(' ')
        hours=float(parted[0])
        mins=float(parted[2].partition(' ')[0])
        secs=float(parted[2].partition(' ')[2])
    elif coord[2]==':':
        parted=coord.partition(':')
        hours=float(parted[0])
        mins=float(parted[2].partition(':')[0])
        secs=float(parted[2].partition(':')[2])
    if coord[0]=='-':
        degs=-1.*(hours*15.+15.*mins/60.+15*secs/3600.)

        warnings.warn('RED FLAG! You shouldn\'t have negative hours!')
    else:
        degs=hours*15.+15.*mins/60.+15.*secs/3600.
    return degs

### Convert ddmmss to decimal degrees
def ddmmssdeg(coord):
    if ((coord[2]=='d') | (coord[3]=='d')):
        degrees=np.abs(float(coord.partition('d')[0]))
        mins=float(coord.partition('d')[2].partition('m')[0])
        secs=float(coord.partition('d')[2].partition('m')[2].partition('s')[0])
    elif ((coord[2]==' ') | (coord[3]==' ')):
        parted=coord.partition(' ')
        degrees=abs(float(parted[0]))
        mins=float(parted[2].partition(' ')[0])
        secs=float(parted[2].partition(' ')[2])
    elif ((coord[2]==':') | (coord[3]==':')):
        parted=coord.partition(':')
        degrees=np.abs(float(parted[0]))
        mins=float(parted[2].partition(':')[0])
        secs=float(parted[2].partition(':')[2])
    degs=degrees+mins/60.+secs/3600.
    if (coord[0]=='-'): degs=-1.*degs
    return degs

### Round to n significant digits
def roundto(x,n,retstring=1):
    if np.isnan(x)==True: return '-'
    if (x==0): number=0
    else: number=round(x, -int(np.floor(np.log10(abs(x))))+(n-1))
    numstr=str(number)
    if ((len(numstr)-1)!=n) & (numstr[0]!='0'):
        leftdig=len(numstr.partition('.')[0])
        rightdig=len(numstr.partition('.')[2])
        zerostoadd=n-leftdig-rightdig
        numstr=numstr+'0'*zerostoadd
        return numstr
    else:
        return number

### Round number to n significant figures
def decimalplaces(x,n):
    #numstr=str(round(x,n))
    numstr="{:f}".format(float(round(x,n)))
    lhs=numstr.split('.')[0]
    rhs=numstr.split('.')[1]
    if len(str(rhs))<n: rhs=rhs+'0'*(n-len(str(rhs)))
    return lhs+'.'+rhs
    
### Volume of a sphere
def volsphere(r):
    return 4./3.*np.pi*r**3

### Relativistic velocity separation between z1,z2
def velsep(z1,z2):
    c=2.99792458e5
    #return c*(((1.+z1)**2-1.)/((1.+z1)**2+1.)-(((z2+1)**2-1)/((z2+1)**2+1)))
    return c*((1.+z1)**2-(1.+z2)**2)/((1.+z1)**2+(1.+z2)**2)

### Combined uncertainty of velocity separation between z1,z2
def sigvelsep(z1,z2,sigz1,sigz2):
    c=2.99792458e5
    difz1=2.*c*(1.+z1)*((1.+z1)**2+(1.+z2)**2)**-1 + -2.*c*(1.+z1)*((1.+z1)**2+(1.+z2)**2)**-2
    difz2=-2.*c*(1.+z2)*((1.+z1)**2+(1.+z2)**2)**-1 + -2.*c*(1.+z2)*((1.+z1)**2+(1.+z2)**2)**-2
    return np.sqrt((difz1*sigz1)**2+(difz2*sigz2)**2)

def SDSSnav(ra,dec):
    pre='http://skyserver.sdss3.org/dr9/en/tools/chart/navi.asp?ra='
    rastr=str(ra)
    decstr=str(dec)
    totstr=pre+rastr+'&dec='+decstr
    if 'linux' in sys.platform:  proc=subprocess.Popen(['xdg-open',totstr])
    else: proc=subprocess.Popen(['open',totstr])

def changefocus(window):
    '''
    Change focus between windows.  Accepts: 'terminal' and 'figure'
    '''
    if window=='terminal':
        if 'linux' in sys.platform:
            subprocess.call(['wmctrl','-a','cardinal'])
        else: 
            subprocess.call(['open','-a','terminal'])    
    if window=='figure':
        if 'linux' in sys.platform:
            subprocess.call(['wmctrl','-a','Figure'])
        else: 
            subprocess.call(['open','-a','XQuartz'])

### Test whether a value is two others
def between(testval,lim1,lim2):
    return ((testval>lim1)&(testval<lim2))



def wherebetween(testval,lim1,lim2):
    '''
	Given multiple ranges of values (lim1 and lim2 are same-sized),
	return indices where testval falls between lim1[idx] and lim2[idx].
	'''
    lim1=np.array(lim1);lim2=np.array(lim2)
    result=np.where((testval>lim1)&(testval<lim2))[0]
    return result


### Generate appropriate calls for a grid of subplots
def subplotgrid(num):
  if num==1: calls=[[1,1,1]]
  elif num==2: calls=[[2,1,1],[2,1,2]]
  elif num==3: calls=[[3,1,1],[3,1,2],[3,1,3]]
  elif num==4: calls=[[4,1,1],[4,1,2],[4,1,3],[4,1,4]]
  elif num==5: calls=[[3,2,1],[3,2,3],[3,2,5],[3,2,2],[3,2,4]]
  elif num==6: calls=[[3,2,1],[3,2,3],[3,2,5],[3,2,2],[3,2,4],[3,2,6]]
  elif num==7: calls=[[4,2,1],[4,2,3],[4,2,5],[4,2,2],[4,2,4],[4,2,6],[4,2,7]]
  elif num==8: calls=[[4,2,1],[4,2,3],[4,2,5],[4,2,7],[4,2,2],[4,2,4],[4,2,6],[4,2,8]]
  elif num==9: calls=[[5,2,1],[5,2,3],[5,2,5],[5,2,7],[5,2,2],[5,2,4],[5,2,6],[5,2,8],[5,2,9]]
  elif num==10: calls=[[5,2,1],[5,2,3],[5,2,5],[5,2,7],[5,2,2],[5,2,4],[5,2,6],[5,2,8],[5,2,9],[5,2,10]]
  elif num==11: calls=[[5,3,1],[5,3,4],[5,3,7],[5,3,10],[5,3,13],[5,3,2],[5,3,5],[5,3,8],[5,3,11],[5,3,14],[5,3,3]]
  elif num==12: calls=[[5,3,1],[5,3,4],[5,3,7],[5,3,10],[5,3,13],[5,3,2],[5,3,5],[5,3,8],[5,3,11],[5,3,14],[5,3,3],[5,3,6]]
  elif num==13: calls=[[5,3,1],[5,3,4],[5,3,7],[5,3,10],[5,3,13],[5,3,2],[5,3,5],[5,3,8],[5,3,11],[5,3,14],[5,3,3],[5,3,6],[5,3,9]]
  elif num==14: calls=[[5,3,1],[5,3,4],[5,3,7],[5,3,10],[5,3,13],[5,3,2],[5,3,5],[5,3,8],[5,3,11],[5,3,14],[5,3,3],[5,3,6],[5,3,9],[5,3,12]]
  elif num==15: calls=[[5,3,1],[5,3,4],[5,3,7],[5,3,10],[5,3,13],[5,3,2],[5,3,5],[5,3,8],[5,3,11],[5,3,14],[5,3,3],[5,3,6],[5,3,9],[5,3,12],[5,3,15]]
  return calls

def subplotgridspec(num):
    if num==1: calls=[[0]]
    elif num==2: calls=[[0],[1]]
    elif num==3: calls=[[0],[1],[2]]
    elif num==4: calls=[[0],[1],[2],[3]]
    elif num==5: calls=[[0,0],[0,1],[0,2],[1,0],[1,1]]
    elif num==6: calls=[[0,0],[1,0],[2,0],[0,1],[1,1],[2,1]]
    elif num==7: calls=[[0,0],[1,0],[2,0],[3,0],[0,1],[1,1],[2,1]]
    elif num==8: calls=[[0,0],[1,0],[2,0],[3,0],[0,1],[1,1],[2,1],[3,1]]
    return calls

