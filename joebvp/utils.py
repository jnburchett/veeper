from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
from matplotlib import pyplot as plt
from joebvp import makevoigt
from joebvp import atomicdata
from linetools.spectra.io import readspec
from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.lists.linelist import LineList
from astropy.table import Table,vstack
from astropy.io import ascii
try:
    import joebvp_cfg
except:
    print("joebvp.utils: No local joebvp_cfg.py found, using default cfg.py file from joebvp.")
    import joebvp.cfg as cfg
import astropy.units as u

def compose_model(spec,filelist,outfile):
    '''
    Generates full-model spectrum (XSpectrum1D) from a set of joebvp output files and writes out a file with the spectrum.

    Parameters
    ----------
    spec : string or XSpectrum1D
        The spectrum to be fitted with the input lines

    filelist : list of strings or str
        This should be a list containing the names of VP input files or a string referring to a file simply
        listing the input files.
        See joebvpfit.readpars for details of file format

    outfile : str
        Name of model spectrum output file

    '''
    from joebvp import joebvpfit

    ### Deal with alternate input types
    if isinstance(spec,str):
        specobj = readspec(spec)
    else:
        specobj = spec

    ### Load essentials from spectrum
    wave = specobj.wavelength.value
    normsig=specobj.sig.value/specobj.co.value
    cfg.wave = wave

    ### Concatenate all the parameter tables
    concatenate_line_tables(filelist,outtablefile='compiledVPinputs.dat')

    ### Make the model!
    fitpars, fiterrors, parinfo, linecmts = joebvpfit.readpars('compiledVPinputs.dat')  # read this back in to load params
    cfg.fitidx = joebvpfit.fitpix(wave, fitpars)  # set the pixels in the line regions
    model = joebvpfit.voigtfunc(wave,fitpars)  # generate the Voigt profiles of all lines

    ### Instantiate XSpectrum1D object and write it out
    outspec=XSpectrum1D.from_tuple((wave,model,normsig))
    outspec.write_to_fits(outfile)

def model_profile(spec,fitpars,instr=None,gratings=None,lsfranges=None,
                  cen_wave=None,lps=None,slits=None):
    '''Produce Voigt profile model from parameters and instrumental setup for
    convolving proper line spread function.  Note

    Parameters
    ----------
    spec : string or XSpectrum1D
        The spectrum to be fitted with the input lines
    fitpars : list of lists
        The joebvp parameter array that includes line measurements
    instr : list of str, optional
        Instruments used to observe the spectral regions defined in lsfranges
    gratings : list of str, optional
        Gratings used for each spectral region defined in lsfranges
    lsfranges : array of size len(instr)x2
        Beginning and ending wavelengths of spectral ranges for each setup
    cen_wave : list of str, optional
        Central wavelength setting for each spectral region defined in lsfranges
    lps : list of str, optional
        Lifetime positions of detector (particularly for COS)
    slits : list of str, optional
        Slit setting for each spectral region defined in lsfranges

    Returns
    -------
    wavelength : Quantity
        Wavelength array from spec
    profile : array
        Voigt profile model

    '''

    ### Setup LSF parameters if provided
    if instr is not None:
        cfg.instr = instr
        cfg.lsfranges = lsfranges
        cfg.gratings = gratings
        cfg.cen_wave = cen_wave
        cfg.lps = lps
        cfg.slits = slits

    ### Deal with alternate input types
    if isinstance(spec,str):
        specobj = readspec(spec)
    else:
        specobj = spec

    makevoigt.get_lsfs()
    #import pdb; pdb.set_trace()
    cfg.wave = specobj.wavelength.value
    cfg.spectrum = specobj
    profile = makevoigt.cosvoigt(cfg.wave, fitpars )
    return specobj.wavelength,profile




def concatenate_line_tables(filelist,outtablefile='compiledVPoutputs.dat'):
    '''
    Compiles the output from several fitting runs into a single table

    Parameters
    ----------
    filelist : list of strings or str
        This should be a list containing the names of VP input files or a string referring to a file simply
        listing the input files.
        See joebvpfit.readpars for details of file format

    outtablefile : str
        Name of compiled model parameter output file

    '''

    if isinstance(filelist, bytes):
        filelist = filelist.decode()
        lstarr = np.genfromtxt(filelist, dtype=None)
        listofiles = lstarr.tolist()
    elif isinstance(filelist, str):
        lstarr=np.genfromtxt(filelist,dtype=None)
        listofiles=lstarr.tolist()
    else:
        listofiles=filelist
    if sum(1 for line in open(filelist)) == 1:
        listofiles = [listofiles]
    tabs = []
    for i, ff in enumerate(listofiles):
        if isinstance(ff,bytes):
            ff = ff.decode()
        #import pdb; pdb.set_trace()
        tabs.append(ascii.read(ff))
    bigpartable = vstack(tabs)
    ascii.write(bigpartable, output=outtablefile, delimiter='|',overwrite=True)  # write out compiled table

def abslines_from_fitpars(fitpars,ra=None,dec=None,linelist=None):
    '''
    Takes a joebvp parameter array and builds a list of linetools AbsLines from the measurements

    Parameters
    ----------
    fitpars : list of lists
        The joebvp parameter array that includes line measurements
    ra : float, optional
        Right Ascension of the QSO in decimal degrees
    dec : float, optional
        Declination of the QSO in decimal degress

    Returns
    -------
    abslinelist: list
        List of AbsLine objects
    '''
    from linetools.spectralline import AbsLine
    import astropy.units as u
    if linelist is None:
        llist = LineList('ISM')
    else:
        llist = linelist

    abslinelist = [] # Initiate list to populate
    for i,rw in enumerate(fitpars[0]):
        vcent = fitpars[4][i]
        v1 = vcent + fitpars[5][i]
        v2 = vcent + fitpars[6][i]

        line=AbsLine(fitpars[0][i]*u.AA, z=fitpars[3][i],linelist=llist)
        vlims=[v1,v2]*u.km/u.s
        line.limits.set(vlims)
        ### Set other parameters
        line.attrib['logN'] = fitpars[1][i]
        #line.attrib['sig_N'] = colerr
        line.attrib['b'] = fitpars[2][i]
        #line.attrib['sig_b'] = berr
        line.analy['spec']=cfg.spectrum
        line.attrib['flag_N'] = 1
        ### Add it to the list and go on
        abslinelist.append(line)
    return abslinelist

def fitpars_from_abslines(abslinelist):
    '''Builds veeper parameter array from list of linetools AbsLines

    Parameters
    ----------
    abslinelist: list
        List of AbsLine objects

    Returns
    -------
    fitpars : list of lists
        The joebvp parameter array that includes line measurements
    '''
    from joebvp import joebvpfit
    restwaves = [al.wrest.value for al in abslinelist]
    linecol = [al.attrib['logN'] for al in abslinelist]
    lineb = [al.attrib['b'].value for al in abslinelist]
    linevel = [al.attrib['vel'].value for al in abslinelist]
    linevlim1 = [al.limits.vlim[0].value for al in abslinelist]
    linevlim2 = [al.limits.vlim[1].value for al in abslinelist]
    zs = [al.z for al in abslinelist]
    initpars = [restwaves, linecol, lineb, zs, linevel, linevlim1, linevlim2]
    initinfo = [[1]*len(abslinelist)]*len(initpars)
    fitpars, parinfo = joebvpfit.initlinepars(zs, restwaves, initpars,
                                              initinfo=initinfo)

    return fitpars

def abslines_from_VPfile(parfile,specfile=None,ra=None,dec=None,linelist=None):
    '''
    Takes a joebvp parameter file and builds a list of linetools AbsLines from the measurements therein.

    Parameters
    ----------
    parfile : str
        Name of the parameter file in the joebvp format
    ra : float, optional
        Right Ascension of the QSO in decimal degrees
    dec : float, optional
        Declination of the QSO in decimal degress

    Returns
    -------
    abslinelist: list
        List of AbsLine objects
    '''
    from linetools.spectralline import AbsLine
    from linetools.analysis import absline as ltaa
    import astropy.units as u
    if linelist is None:
        llist = LineList('ISM')
    else:
        llist = linelist
    if specfile!=None:
        spec=readspec(specfile) # Allow spectrum file to be declared in call
    linetab = ascii.read(parfile) # Read parameters from file
    linetab['restwave']=linetab['restwave']*u.AA
    abslinelist = [] # Initiate list to populate
    for i,row in enumerate(linetab):
        ### Check to see if errors for this line are defined
        colerr,berr,velerr=get_errors(linetab,i)
        ### Adjust velocity limits according to centroid errors and limits from file
        vcentmin = row['vel']-velerr
        vcentmax = row['vel']+velerr
        v1 = vcentmin + row['vlim1']
        v2 = vcentmax + row['vlim2']
        line=AbsLine(row['restwave']*u.AA, z=row['zsys'],closest=True, linelist=llist)
        vlims=[v1,v2]*u.km/u.s
        line.limits.set(vlims)
        ### Set other parameters
        line.attrib['logN'] = row['col']
        line.attrib['sig_logN'] = colerr
        line.attrib['flag_N'] = 1
        lincol = ltaa.linear_clm(line.attrib)
        line.attrib['N']= lincol[0]
        line.attrib['sig_N'] = lincol[1]
        line.attrib['b'] = row['bval'] * u.km/u.s
        line.attrib['sig_b'] = berr * u.km/u.s
        line.attrib['vel'] = row['vel'] * u.km/u.s
        ### Attach the spectrum to this AbsLine but check first to see if this one is same as previous
        if specfile==None:
            if i==0:
                spec=readspec(row['specfile'])
            elif row['specfile']!=linetab['specfile'][i-1]:
                spec=readspec(row['specfile'])
            else:
                pass
        line.analy['spec']=spec
        ### Add it to the list and go on
        abslinelist.append(line)
    return abslinelist

def get_errors(partable,idx2check):
    '''
    Get actual error values from VP fitting.
    When lines are tied in fit, algorithm returns '0' for all lines but 1 in fit.
    This routine will return errors from row in table that actually have errors.

    Parameters
    ----------
    partable : astropy Table
        Imported Table from parameter file in the joebvp format
    idx2check : int
        Index of row to
    Returns
    -------
    colerr : float
        Fitting error of column density
    berr: float
        Fitting error of Doppler b value
    velerr: float
        Fitting error of velocity centroid
    '''

    row = partable[idx2check]

    ### Process column density, Doppler b value, then velocity centroid
    if row['sigcol']==0:  # Check to see if errors for this line are defined
        simulfit=np.where((partable['nflag']==row['nflag'])&(partable['zsys']==row['zsys']))[0] # Find lines tied to this one
        if len(simulfit)>0: # Are there any?
            simulfit_nz=simulfit[np.where(partable['nflag'][simulfit]!=0)[0]]  # Find the one(s) with nonzero errors
            if len(simulfit_nz)>0: # Are there any?
                colerr=partable['sigcol'][simulfit_nz[0]] # Adopt nonzero error of this row
            else:
                colerr=row['sigcol'] # No nonzero errors; accept value from table read in
        else:
            colerr=row['sigcol'] # No rows tied to this one; accept value from table read in
    else:
        colerr=row['sigcol'] # This row has nonzero errors; accept value from table read in

    if row['sigbval']==0:
        simulfit=np.where((partable['bflag']==row['bflag'])&(partable['zsys']!=row['zsys']))[0]
        if len(simulfit)>0:
            simulfit_nz=simulfit[np.where(partable['bflag'][simulfit]!=0)[0]]
            if len(simulfit_nz)>0:
                berr=partable['sigbval'][simulfit_nz[0]]
            else:
                berr=row['sigbval']
        else:
            berr=row['sigbval']
    else:
        berr=row['sigbval']

    if row['sigvel']==0:
        simulfit=np.where((partable['vflag']==row['vflag'])&(partable['zsys']!=row['zsys']))[0]
        if len(simulfit)>0:
            simulfit_nz=simulfit[np.where(partable['vflag'][simulfit]!=0)[0]]
            if len(simulfit_nz)>0:
                velerr=partable['sigvel'][simulfit_nz[0]]
            else:
                velerr=row['sigvel']
        else:
            velerr=row['sigvel']
    else:
        velerr=row['sigvel']

    return colerr, berr, velerr
    
def inspect_fits(parfile,output='FitInspection.pdf',vlim=[-300,300]*u.km/u.s, **kwargs):
    '''
    Produce pdf of fitting results for quality check.

    Parameters
    ----------
    parfile : str
        Name of the parameter file in the joebvp format
    output : str,optional
        Name of file to write stackplots output
    vlim : Quantity array, optional
        Velocity range of stackplots
        e.g.: [-400,

    Returns
    -------

    '''
    from joebvp import joebvpfit
    from matplotlib.backends.backend_pdf import PdfPages
    llist = LineList('ISM')

    pp=PdfPages(output)

    #import pdb;
    #pdb.set_trace()
    all=abslines_from_VPfile(parfile,linelist=llist,**kwargs) # Instantiate AbsLine objects and make list
    acl=abscomponents_from_abslines(all,vtoler=15.)  # Instantiate AbsComponent objects from this list
    fitpars,fiterrors,parinfo,linecmts = joebvpfit.readpars(parfile)


    fullmodel=makevoigt.cosvoigt(all[0].analy['spec'].wavelength.value,fitpars)

    ### Make a stackplot for each component
    for comp in acl:
        fig=comp.stack_plot(ymnx=(-0.1,1.3),show=False,return_fig=True,vlim=vlim, tight_layout=True)
        if (len(comp._abslines)<6):
            numrow = len(comp._abslines)%6
        else:
            numrow = 6
        height = numrow*1.0+0.5
        fig.set_figheight(height)
        if len(comp._abslines)<6:
            fig.set_figwidth(5.)
        else:
            fig.set_figwidth(10.)
        stackaxes = fig.axes

        for i,ax in enumerate(stackaxes):
            line = comp._abslines[i]
            thesepars=[[line.wrest.value],[line.attrib['logN']],[line.attrib['b'].value],
                       [line.z],[line.attrib['vel'].value],[line.limits.vlim[0].value],[line.limits.vlim[1].value]]
            try:
                thismodel=makevoigt.cosvoigt(line.analy['spec'].wavelength.value,thesepars)
            except:
                ax.text(-100,0.05,'These pixels not used to fit this line')
                continue
            axlin=ax.get_lines()
            veldat=axlin[0].get_data()[0]
            ax.plot(veldat,fullmodel,color='red',alpha=0.8)
            ax.plot(veldat,thismodel,color='purple',linestyle='dashed')
        botbuf = 0.5/height
        fig.subplots_adjust(bottom=botbuf,left=0.1,right=0.95,hspace=0.,wspace=0.35)

        title=comp.name
        fig.suptitle(title)
        fig.savefig(pp,format='pdf')
        plt.close(fig)

    ### Write out full model fits file
    normflux = line.analy['spec'].flux/line.analy['spec'].co
    normsig = line.analy['spec'].sig / line.analy['spec'].co
    modeltab = Table([line.analy['spec'].wavelength.value, fullmodel, normflux, normsig], names=['wavelength', 'model', 'normflux', 'normsig'])
    # modeltab.write(outfile, format='fits', overwrite=True)
    dummycont = np.ones(len(normflux))
    spec = XSpectrum1D.from_tuple((modeltab['wavelength'], modeltab['model'], modeltab['normsig'], dummycont))
    spec.write_to_fits(output[:-4]+'.fits')

    pp.close()

def abscomponents_from_abslines(abslinelist, **kwargs):
    '''
    Organizes list of AbsLines into components based on the following order: redshift, velocity, and species

    Parameters
    ----------
    abslinelist : list
        List of AbsLine objects

    Returns
    -------
    complist : list
        List AbsComponent objects
    '''
    ### Import machine learning clustering algorithm for velocity grouping
    from sklearn.cluster import MeanShift, estimate_bandwidth
    from linetools.isgm.abscomponent import AbsComponent

    ### Populate arrays with line redshifts, velocity centroids, and name of species
    zarr=np.zeros(len(abslinelist))
    varr=np.zeros(len(abslinelist))
    sparr=np.chararray(len(abslinelist),itemsize=6)
    for i,absline in enumerate(abslinelist):
        zarr[i]=absline.z
        # import pdb;        pdb.set_trace()
        varr[i]=absline.attrib['vel'].value
        sparr[i]=atomicdata.lam2ion(absline.wrest.value)

    abslinelist=np.array(abslinelist) # Convert to array for the indexing used below

    ### Group lines with the same redshifts and similar velocities
    complines=[]
    uqzs=np.unique(zarr)
    for uqz in uqzs:
        ### Identify velocity groups
        thesez = np.where(zarr==uqz)[0]
        X = np.array(list(zip(varr[thesez],np.zeros(len(varr[thesez])))), dtype=float)
        ms = MeanShift(bandwidth=7.)
        ms.fit(X)
        vidxs = ms.labels_
        vs = ms.cluster_centers_
        uqvidxs=np.unique(vidxs)
        ### Make lists of lines that will belong to each component
        for idx in uqvidxs:
            theseuqvs=thesez[np.where(vidxs==idx)[0]] # isolate velocity-grouped lines with this redshift
            uqsp=np.unique(sparr[theseuqvs]) # identify the unique species with lines in this group
            ### Get lines belonging to each species and add them to become components
            for sp in uqsp:
                spidx=theseuqvs[np.where(sparr[theseuqvs]==sp)]
                complines.append(abslinelist[spidx])
    ### Instantiate the AbsComponents
    comps=[]
    for lst in complines:
        if '*' in lst[0].name:
            linename=lst[0].name
            starct=linename.count('*')
            stars='*'*starct
        else:
            stars = None
        thiscomp=AbsComponent.from_abslines(lst.tolist(), stars=stars,
                                            chk_vel=False, adopt_median=True,
                                            **kwargs)
        comps.append(thiscomp)
    return comps

def add_continuum(specfile,outfilename):
    '''Append dummy continuum of 1s for normalized spectra'''
    from linetools.spectra.io import readspec
    spec = readspec(specfile)
    spec.co = np.ones(len(spec.flux))
    spec.meta['headers'] = [spec.header]
    spec.write_to_fits(outfilename,clobber=True)

