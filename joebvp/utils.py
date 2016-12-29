import numpy as np
import joebvpfit
import makevoigt
import atomicdata
from linetools.spectra.io import readspec
from linetools.spectra.xspectrum1d import XSpectrum1D
from astropy.table import Table,vstack
from astropy.io import ascii
import joebvp.cfg as cfg

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
    fitpars, fiterrors, parinfo = joebvpfit.readpars('compiledVPinputs.dat')  # read this back in to load params
    cfg.fitidx = joebvpfit.fitpix(wave, fitpars)  # set the pixels in the line regions
    model = joebvpfit.voigtfunc(wave,fitpars)  # generate the Voigt profiles of all lines

    ### Instantiate XSpectrum1D object and write it out
    outspec=XSpectrum1D.from_tuple((wave,model,normsig))
    outspec.write_to_fits(outfile)

def concatenate_line_tables(filelist,outtablefile='compiledVPinputs.dat'):
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

    if isinstance(filelist, str):
        lstarr=np.genfromtxt(filelist,dtype=None)
        listofiles=lstarr.tolist()
    else:
        listofiles=filelist

    tabs = []
    for i, ff in enumerate(listofiles):
        tabs.append(ascii.read(ff))
    bigpartable = vstack(tabs)
    ascii.write(bigpartable, output=outtablefile, delimiter='|')  # write out compiled table

def abslines_from_VPfile(parfile,ra=None,dec=None):
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
    import astropy.units as u
    linetab = ascii.read(parfile) # Read parameters from file
    abslinelist = [] # Initiate list to populate
    for i,row in enumerate(linetab):
        ### Check to see if errors for this line are defined
        colerr,berr,velerr=get_errors(linetab,i)
        ### Adjust velocity limits according to centroid errors and limits from file
        vcentmin = row['vel']-velerr
        vcentmax = row['vel']+velerr
        v1 = vcentmin + row['vlim1']
        v2 = vcentmax + row['vlim2']
        line=AbsLine(row['restwave']*u.AA, z=row['zsys'])
        vlims=[v1,v2]*u.km/u.s
        line.limits.set(vlims)
        ### Set other parameters
        line.attrib['logN'] = row['col']
        line.attrib['sig_N'] = colerr
        line.attrib['b'] = row['bval']
        line.attrib['sig_b'] = berr
        ### Attach the spectrum to this AbsLine but check first to see if this one is same as previous
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
    
def inspect_fits(parfile,output='FitInspection.pdf'):
    '''
    Produce pdf of fitting results for quality check.

    Parameters
    ----------
    parfile : str
        Name of the parameter file in the joebvp format

    output : str
        Name of file to write stackplots output

    Returns
    -------

    '''

    from matplotlib.backends.backend_pdf import PdfPages
    pp=PdfPages(output)

    all=abslines_from_VPfile(parfile) # Instantiate AbsLine objects and make list
    acl=abscomponents_from_abslines(all,vtoler=15.)  # Instantiate AbsComponent objects from this list

    ### Make a stackplot for each component
    for comp in acl:
        fig=comp.stack_plot(show=False,return_fig=True, tight_layout=True)
        title=comp.name
        fig.suptitle(title)
        fig.savefig(pp,format='pdf')
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
        varr[i]=np.mean(absline.limits.vlim.value)
        sparr[i]=atomicdata.lam2ion(absline.wrest.value)

    abslinelist=np.array(abslinelist) # Convert to array for the indexing used below

    ### Group lines with the same redshifts and similar velocities
    complines=[]
    uqzs=np.unique(zarr)
    for uqz in uqzs:
        ### Identify velocity groups
        thesez = np.where(zarr==uqz)[0]
        X = np.array(zip(varr[thesez],np.zeros(len(varr[thesez]))), dtype=float)
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
        thiscomp=AbsComponent.from_abslines(lst.tolist(), **kwargs)
        comps.append(thiscomp)
    return comps

