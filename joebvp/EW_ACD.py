from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
from joebvp import joebgoodies
from joebvp.atomicdata import atomicdata
import numpy.polynomial.legendre as L
from joebvp import cfg
from scipy import stats

def contFitLegendreAboutLine(wave,flux,err,restlam,z,velfitregions,uniform=True,**kwargs):
    '''
    Fits a continuum over a spectral region using the formalism of Sembach
    & Savage 1992, estimating the uncertainty at each point of the continuum.

    Parameters
    ----------
    wave
    flux
    err
    restlam
    z
    velfitregions
    uniform

    Returns
    -------

    '''
    vels=joebgoodies.veltrans(z,wave,restlam)

    ### Find indices corresponding to velocities defined in velfitregions
    fitidxs=[]
    for reg in velfitregions:
        thesevels=np.where((vels>=reg[0])&(vels<=reg[1]))[0]
        fitidxs.extend(thesevels)

    if uniform != True:
        fitsol,err,parsol,errmtx,scalefac = fitLegendre(vels[fitidxs],flux[fitidxs],
                                           sig=err[fitidxs],**kwargs)
    else:
        fitsol,err,parsol,errmtx,scalefac = fitLegendre(vels[fitidxs],flux[fitidxs],**kwargs)

    ### Normalize x so that domain of fit is -1 to +1
    vscale = vels/scalefac # must scale to fitted domain

    ### Evaluate continuum using all points, not just those in fit regions
    continuum = L.legval(vels/scalefac,parsol)
    err = errorsLegendre(vscale,errmtx)

    return wave,continuum,err

def fitLegendre(x, y, sig=None, minord=1, maxord=8):
    '''
    Fits data with a series of Legendre polynomials and provides errors at each
    point in the fit due to the uncertainty in the fit parameters.  The order
    of the Legendre polynomial series will be determined by an F-test and
    constrained between minord and maxord.

    Parameters
    ----------
    x: 1D array
        Independent variable for fit
    y: 1D array
        Dependent variable for fit
    sig: 1D array, optional
        Measurement errors on y; if not declared, fit will proceed with uniform weighting
    minord: int, optional
         Lowest order of Legendre polynomial series to try for fit
    maxord: int, optional
         Highest order of Legendre polynomial series to try for fit

    Returns
    -------
    fitsol: 1D array
        Fitted curve to data
    err: 1D array
        Uncertainty in fit at each point
    parsol: 1D array
        Coefficients of Legendre polynomial series fit to data
    errmtx: 2D matrix
        Error matrix of fit
    '''

    ### If no data errors are given, use uniform weighting in fit
    if sig==None:
        uniform = True
        sig=np.ones(len(y))
    else:
        uniform = False

    ### Normalize x so that domain of fit is -1 to +1
    scalefac = np.max(np.abs(x))
    xscale = x/scalefac

    ### Construct matrix of Legendre polynomials
    legmtx = basisMatrixLegendre(xscale,maxord)

    ### Prepare lists for comparison in F-test
    sumsqdiffs=[]
    dfs=[]
    errmtxs=[]
    print(minord)
    ### Loop through possible orders between minord and maxord
    for ord in range(minord,maxord+1):

        ### Prepare alpha & beta matrices
        alpha = np.zeros([ord+1,ord+1])
        beta = np.zeros(ord+1)
        for i in range(ord+1):
            beta[i]=np.sum(y*legmtx[i,:]/sig**2)
            for j in range(ord+1):
                alpha[j,i]=np.sum(legmtx[j,:]*legmtx[i,:]/sig**2)
        alphamtx=np.matrix(alpha)
        betamtx=np.matrix(beta)

        ### Evaluate parameters and fit
        errmtx=alphamtx.I # The error matrix is simply the inverse of alpha
        errmtxs.append(errmtx) # Save to possibly adopt later after F-test
        parsol=np.array(betamtx*errmtx)[0] # Parameter solution
        fitsol = L.legval(xscale, parsol) # Evaluate Leg. poly. series for fit solution

        ### Prepare for F-test
        dfs.append(len(xscale) - len(parsol) - 1.) # degress of freedom
        sumsqdiffs.append(np.sum((y-fitsol)**2.))  # numerator of chi-square

        ### Conduct F-test if we've been through loop already
        if ord!=minord:
            fval = sumsqdiffs[-1]/sumsqdiffs[-2] # F = ratio of variances
            ### Get probability that F < fval in random distribution
            fprob = stats.f.cdf(fval,dfs[-1],dfs[-2])
            ### If this variance is not less than previous w/95% confidence, stop
            if ((fprob >= 0.05)|(ord==maxord)):
                ### Use fit of previous order
                ord = ord - 1
                errmtx = errmtxs[-2]
                parsol = parsol[0:ord+1] # keep only parameters up to last fit
                fitsol = L.legval(xscale, parsol) # evaluate one last time
                ssd = sumsqdiffs[-2]
                df = dfs[-2]
                break

    if uniform!=False: # Need to scale reduced error matrix if no errors on data
        errmtx = errmtx * ssd/df

    err = errorsLegendre(xscale,errmtx)

    return fitsol,err,parsol,errmtx,scalefac

def errorsLegendre(x,errmtx):
    '''
    Evaluates errors at each point in fit using the covariance matrix

    Parameters
    ----------
    x: 1D array
        Points at which to evaluate fit errors
    errmtx: numpy 2d array of shape (order+1,order+1)
        Covariance matrix from fit

    Returns
    -------
    err: numpy 1d array
        Error at each point
    '''

    order = len(errmtx)-1
    legmtx = basisMatrixLegendre(x,order)

    numrows = len(legmtx)
    numcols = len(legmtx[0])

    ### Populate error vector
    errsq = np.zeros(numcols)
    for k in range(numcols): # Evaluate including covariances of fit coeffs
        for i in range(numrows):
            for j in range(numrows):
                errsq[k]+=errmtx[i,j]*legmtx[j,k]*legmtx[i,k]
    err = np.sqrt(errsq)

    return err

def basisMatrixLegendre(x,order):
    '''
    Builds (order+1,x) 2-d array of Legendre polynomials evaluated at x

    Parameters
    ----------
    x: numpy 1-d array
        Points at which to evaluate basis functions
    order:
        Maximum order of polynomial set

    Returns
    -------
    legmtx: numpy 2d array of shape (order+1,x)
        Legendre polynomials evaluated at x
    '''

    ### Construct matrix of Legendre polynomials
    legmtx = np.zeros((order+1,len(x)))

    ### Prepare matrix of Leg. poly. basis functions evaluated at x
    for i in range(order+1):
        coeffs=np.zeros(i+1)
        coeffs[i]=1.
        legmtx[i,:]=L.legval(x, coeffs)

    return legmtx

def EW_ACD_array(wave,flux,ferr,cont,conterr,restlam,zabs,vellim=[-50,50],**kwargs):
    '''
    Returns arrays of equivalent width and apparent column density per pixel as
    well as their associated uncertainties due to continuum placement and flux errors.

    Parameters
    ----------
    wave: 1D float array
    flux: 1D float array
    ferr: 1D float array
        Error in flux
    cont: 1D float array
        Continuum fitted to data such as that from contFitLegendreAboutLine()
    conterr: 1D float array
        Uncertainty in continuum placement at each pixel
    restlam: float
        Rest-frame wavelength of transition to measure
    zabs: float
        Redshift of the the absorption system
    vellim: 2-element list
        Lower and upper bounds over which to compute arrays

    Returns
    -------
    EWpix: 1D float array
        Equivalent width in each pixel
    sigEWf: 1D float array
        EW uncertainty in each pixel due to flux errors
    sigEWc: 1D float array
        EW uncertainty in each pixel due to continuum placement errors
    Npix: 1D float array
        Apparent column density * dv (N) in each pixel
    sigNf: 1D float array
        Uncertainty in N due to flux errors
    sigNc: 1D float array
        Uncertainty in N due to continuum placement errors
    '''

    ### Set atomic data and constants
    c = 299792.458
    restlam=atomicdata.closestlam([restlam])
    osc=atomicdata.lam2osc([restlam])

    ### Transform to velocity space
    vel=joebgoodies.veltrans(zabs,wave,restlam)
    velidx=np.where((vel>=vellim[0])&(vel<=vellim[1]))[0]
    velup = vel[velidx+1]
    veldown = vel[velidx]
    dv = np.abs(velup-veldown)

    ### Identify pixels where flux level is below noise level & replace w/noise
    effflux = flux
    belowerr = np.where(flux<ferr)[0]
    effflux[belowerr] = ferr[belowerr]

    ### Calculate EW and errors due to flux and continuum placement
    print(dv,effflux[velidx],cont[velidx],restlam,c)
    EWpix = dv * (1.-effflux[velidx]/cont[velidx])*restlam/c
    sigEWf = dv / cont[velidx] * ferr[velidx] * restlam/c
    sigEWc = dv * restlam/c * effflux[velidx] / cont[velidx]**2 * conterr[velidx]  # Will not be added in quad.

    ### Calculate optical depth and uncertainty due to flux and continuum
    tauv=np.log(cont[velidx]/(effflux[velidx]))
    tauverr_f = ferr[velidx]/effflux[velidx]
    tauverr_c = conterr[velidx]/cont[velidx]
    Npix=1./2.654e-15/restlam/osc*dv*tauv
    sigNf = 1./2.654e-15/restlam/osc*dv*tauverr_f
    sigNc=1./2.654e-15/restlam/osc*dv*tauverr_c

    return EWpix,sigEWf,sigEWc,Npix,sigNf,sigNc

def EW_SS92err(wave,flux,ferr,cont,conterr,restlam,zabs,vellim=[-50,50]):
    '''
    Calculates the equivalent width, apparent column density, and their
    associated errors a la Sembach & Savage 1992.

    '''

    EWpix,sigEWf,sigEWc,Npix,sigNf,sigNc = EW_ACD_array(wave,flux,ferr,cont,conterr,
                                            restlam,zabs,vellim=vellim)

    ### Totals and errors from each contribution
    EW = np.sum(EWpix)
    N = np.sum(Npix)
    ## Sum flux error contributions in quadrature
    sigEWf_tot = np.sqrt(np.sum(sigEWf**2))
    sigNf_tot = np.sqrt(np.sum(sigNf**2))
    ## Sum continuum error contributions not in quadrature
    sigEWc_tot = np.sum(sigEWc)
    sigNc_tot = np.sum(sigNc)

    return EW,sigEWf_tot,sigEWc_tot,N,sigNf_tot,sigNc_tot

def vel_moment_fitcont(wave,flux,sig,restlam,zabs,continuumregions,vellim=[-50,50],**kwargs):
    ### Fit the continuum in the line-free regions
    wave,cont,conterr=contFitLegendreAboutLine(wave,flux,sig,restlam,zabs,continuumregions,**kwargs)
    EWpix,sigEWf,sigEWc,Npix,sigNf,sigNc = EW_ACD_array(wave,flux,sig,cont,conterr,
                                            restlam,zabs,vellim=vellim)
    ### Transform to velocity space
    vel=joebgoodies.veltrans(zabs,wave,restlam)
    velidx=np.where((vel>=vellim[0])&(vel<=vellim[1]))[0]
    velup = vel[velidx+1]
    veldown = vel[velidx]

    ### Calculate moments
    moment0 = np.sum(Npix)
    moment1 = np.sum(vel[velidx]*Npix)

    ### Mean velocity
    meanv = moment1/moment0
    return meanv

def vel_moment(spec,restlam,zabs,vellim=[-50,50]):
    wave=spec.wavelength.value
    flux=spec.flux.value
    sig=spec.sig.value
    cont=spec.co
    conterr=np.zeros_like(sig) # Continuum errors will be garbage
    EWpix,sigEWf,sigEWc,Npix,sigNf,sigNc = EW_ACD_array(wave,flux,sig,cont,conterr,
                                            restlam,zabs,vellim=vellim)
    ### Transform to velocity space
    vel=joebgoodies.veltrans(zabs,wave,restlam)
    velidx=np.where((vel>=vellim[0])&(vel<=vellim[1]))[0]


    ### Calculate moments
    moment0 = np.sum(Npix)
    moment1 = np.sum(vel[velidx]*Npix)


    ### Mean velocity
    meanv = moment1/moment0
    return meanv






