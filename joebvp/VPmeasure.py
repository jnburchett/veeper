# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 22:50:45 2016

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
try:
    import joebvp_cfg as cfg
except:
    print("joebvp.VPmeasure: No local joebvp_cfg.py found, using default cfg.py file from joebvp.")
    from joebvp import cfg
from joebvp import joebvpfit
from joebvp import utils as jbu
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

modpath=os.path.abspath(os.path.dirname(__file__))
print(os.path.abspath(os.path.dirname(__file__)))
Ui_MainWindow, QMainWindow = loadUiType(modpath+'/mainvpwindow.ui')

matplotlib.rcParams['font.size']=cfg.general_fontsize

c= c.to('km/s').value

class LineParTableModel(QAbstractTableModel):
    def __init__(self,fitpars,fiterrors,parinfo,linecmts=None,parent=None):
            QAbstractTableModel.__init__(self,parent)
            self.fitpars=fitpars
            self.fiterrors=fiterrors
            self.parinfo=parinfo
            self.linecmts=linecmts
            self.headers=['Wavelength','Species','z','N','sig(N)','b','sig(b)','v','sig(v)','rely','comment']

    def rowCount(self,parent):
        return len(self.fitpars[0])

    def columnCount(self, parent):
        return 11

    def data(self,index,role):
        if role == Qt.EditRole:
            row = index.row()
            column = index.column()
            if column == 1: toret=atomicdata.lam2ion(self.fitpars[0][row])
            elif column == 2: toret=round(self.fitpars[3][row],5)
            elif column == 3: toret=round(self.fitpars[1][row],3)
            elif column == 4: toret=round(self.fiterrors[1][row],3)
            elif column == 5: toret=round(self.fitpars[2][row],3)
            elif column == 6: toret=round(self.fiterrors[2][row],3)
            elif column == 7: toret=round(self.fitpars[4][row],3)
            elif column == 8: toret=round(self.fiterrors[4][row],3)
            elif column == 9: toret=self.linecmts[0][row]
            elif column == 10: toret=self.linecmts[1][row]
            else: toret=QVariant(round(self.fitpars[column][row],3))
            return toret

        if index.isValid() and role == Qt.DisplayRole:
            row = index.row()
            column = index.column()
            if column == 0: toret=round(self.fitpars[0][row],3)
            elif column == 1:
                toret=atomicdata.lam2ion(self.fitpars[0][row])
            elif column == 2: toret=round(self.fitpars[3][row],5)
            elif column == 3: toret=round(self.fitpars[1][row],3)
            elif column == 4: toret=round(self.fiterrors[1][row],3)
            elif column == 5: toret=round(self.fitpars[2][row],3)
            elif column == 6: toret=round(self.fiterrors[2][row],3)
            elif column == 7: toret=round(self.fitpars[4][row],3)
            elif column == 8: toret=round(self.fiterrors[4][row],3)
            elif column == 9: toret=self.linecmts[0][row]
            elif column == 10: toret=self.linecmts[1][row]
            else: toret=QVariant(round(self.fitpars[column][row],3))
            return str(toret)

        else: return None

    def flags(self,index):
        return Qt.ItemIsEditable | Qt.ItemIsEnabled | Qt.ItemIsSelectable


    def setData(self,index,value,role=Qt.EditRole):
        if role == Qt.EditRole:
            row = index.row()
            column = index.column()
            if value.toFloat()[1]==True:
                val = round(float(value.toFloat()[0]),2)
                if column == 1: pass
                if column == 2: self.fitpars[3][row] = val
                elif column == 3: self.fitpars[1][row] = val
                elif column == 4: self.fiterrors[1][row] = val
                elif column == 5: self.fitpars[2][row] = val
                elif column == 6: self.fiterrors[2][row] = val
                elif column == 7: self.fitpars[4][row] = val
                elif column == 8: self.fiterrors[4][row] = val
                elif column == 9: self.linecmts[0][row] = value
                elif column == 10: self.linecmts[1][row] = value
                else: self.fitpars[column][row] = val
                self.dataChanged.emit(index,index)
                return True

        return False


    def addLine(self,wave,newrestwave,newz,newcol,newb,newvel,newvel1,newvel2,newrely=None,newcmt=None):
        ### Setup new parameters and append them to main arrays of the data model
        newpars=[[newrestwave],[newcol],[newb],[newz],[newvel],[newvel1],[newvel2]]
        fitpars=np.hstack((self.fitpars,newpars))
        newerrors=[[-99],[-99],[-99],[-99],[-99]]
        fiterrors=np.hstack((self.fiterrors,newerrors))
        newindex=np.max(self.parinfo[1])+1
        newinfo=[[1],[newindex],[newindex],[1],[newindex]]
        parinfo=np.hstack((self.parinfo,newinfo))
        newcmts=[newrely,newcmt]
        linecmts=np.hstack((self.linecmts,newcmts))
        ### Call initlinepars to set atomic data in cfg.fosc, etc.
        junk,junk=joebvpfit.initlinepars(fitpars[3],fitpars[0],initvals=fitpars,initinfo=parinfo)
        ### Do the update
        midx=QModelIndex()  # Needed for 'beginInsertRows'
        self.beginInsertRows(midx,len(self.fitpars[0]),len(self.fitpars[0]))
        self.updatedata(fitpars,fiterrors,parinfo,linecmts)
        self.endInsertRows()
        ### Reset pixels for fit and wavegroups for convolution
        cfg.fitidx=joebvpfit.fitpix(wave, fitpars) #Reset pixels for fit
        cfg.wavegroups=[]
    '''
    def insertRows(self,position,rows,parent=QModelIndex()):
        self.beginInsertRows(parent,position,position+rows-1)
        if isinstance(newrestwave,float):
            self.fitpars[0].extend([newrestwave])
            self.fitpars[1].extend([newcol])
            self.fitpars[2].extend([newb])
            self.fitpars[3].extend([newz])
            self.fitpars[4].extend([newvel])
            self.fitpars[5].extend([-cfg.defaultvlim])
            self.fitpars[6].extend([cfg.defaultvlim])
            self.fiterrors[0].extend([0.]) ;  self.fiterrors[1].extend([0.]) ;  self.fiterrors[2].extend([0.])
            self.fiterrors[3].extend([0.]) ;  self.fiterrors[4].extend([0.])
        for i in range(rows):
            defaultFitpars=[[0]*len(self.fitpars)]
            defaultFiterrors = [[0] * len(self.fiterrors)]
            defaultParinfo = [[0] * len(self.parinfo)]
            self.fitpars.insert(position,defaultFitpars)
        self.endInsertRows()
    '''

    def updatedata(self,fitpars,fiterrors,parinfo,linecmts):
        self.fitpars=fitpars
        self.fiterrors=fiterrors
        self.parinfo=parinfo
        self.linecmts=linecmts
        self.dataChanged.emit(self.index(0,0),self.index(self.rowCount(self)-1,self.columnCount(self)-1))


    def headerData(self,section,orientation,role):
        if role == Qt.DisplayRole:
            if orientation==Qt.Horizontal:
                return self.headers[section]

class newLineDialog(QDialog):
    def __init__(self, parent = None):
        super(newLineDialog,self).__init__(parent)

        layout = QGridLayout(self)

        self.lamLabel=QLabel(self)
        self.lamLabel.setText('Rest Wavelength:')
        layout.addWidget(self.lamLabel,0,0)
        self.lamBox = QLineEdit(self)
        #self.lamBox.setMaximumWidth = 30
        self.lamBox.editingFinished.connect(self.validateWavelength)
        layout.addWidget(self.lamBox,1,0)
        self.zLabel=QLabel(self)
        self.zLabel.setText('z:')
        layout.addWidget(self.zLabel,0,1)
        self.zBox = QLineEdit(self)
        #self.lamBox.setMaximumWidth = 30
        layout.addWidget(self.zBox,1,1)
        self.ionLabel = QLabel(self)
        layout.addWidget(self.ionLabel, 1, 2)
        
        self.colLabel=QLabel(self)
        self.colLabel.setText('N:')
        layout.addWidget(self.colLabel,2,0)
        self.colBox = QLineEdit(self)
        self.colBox.setMaximumWidth = 30
        layout.addWidget(self.colBox,3,0)
        self.bLabel=QLabel(self)
        self.bLabel.setText('b:')
        layout.addWidget(self.bLabel,2,1)
        self.bBox = QLineEdit(self)
        layout.addWidget(self.bBox,3,1)
        self.velLabel=QLabel(self)
        self.velLabel.setText('vel:')
        layout.addWidget(self.velLabel,2,2)
        self.velBox = QLineEdit(self)
        layout.addWidget(self.velBox,3,2)
        
        

        buttons = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel,
            Qt.Horizontal, self)
        layout.addWidget(buttons)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)


    def validateWavelength(self):
        ### Handle approximate wavelength entries by looking up precise values on focus exit
        if self.lamBox.text()!='':
            restlam=self.lamBox.text()
            restlam=float(restlam)
            try:
                lam,fosc,gam=atomicdata.setatomicdata([restlam],precise=False)
                self.lamBox.setText(str(lam[0]))
                self.ionLabel.setText(atomicdata.lam2ion(lam[0]))
            except:
                pass

    def lineParams(self):
        vel=float(self.velBox.text())
        vel1= vel - cfg.defaultvlim ; vel2= vel + cfg.defaultvlim
        return self.lamBox.text(),self.zBox.text(),self.colBox.text(),self.bBox.text(),self.velBox.text(),vel1,vel2

    @staticmethod
    def get_newline(parent = None):
        dialog = newLineDialog(parent)
        result = dialog.exec_()
        if result == 1:
            newlam,newz,newcol,newb,newvel,newvel1,newvel2 = dialog.lineParams()
            return newlam,newz,newcol,newb,newvel,newvel1,newvel2
        else:
            return 0

class Main(QMainWindow, Ui_MainWindow):
    def __init__(self,specfilename,parfilename=None,wave1=None,wave2=None,numchunks=8,parent=None):
        QtWidgets.QMainWindow.__init__(self, parent)
        #super(Main,self).__init__()
        self.setupUi(self)

        ### Initialize stuff
        self.line_dict = {}
        self.fitpars = None
        self.parinfo = None
        self.linecmts = None
        self.wave1 = wave1
        self.wave2 = wave2
        self.numchunks = numchunks
        self.spls = []
        self.labeltog = 1
        self.pixtog = 0
        self.restog = 1
        self.fitconvtog = 0
        self.lastclick=1334.

        ### Read in spectrum and list of lines to fit
        self.specfilename=specfilename
        self.spectrum = readspec(specfilename)
        self.wave=self.spectrum.wavelength.value
        self.normflux=self.spectrum.flux/self.spectrum.co
        self.normsig=self.spectrum.sig/self.spectrum.co
        cfg.spectrum = self.spectrum
        cfg.wave=self.wave
        cfg.normflux=self.normflux
        cfg.filename=self.specfilename

        if not parfilename==None:
            self.initialpars(parfilename)


        ### Connect signals to slots
        self.fitButton.clicked.connect(self.fitlines)
        self.fitConvBox.clicked.connect(self.togfitconv)
        self.boxLineLabel.clicked.connect(self.toglabels)
        self.boxFitpix.clicked.connect(self.togfitpix)
        self.boxResiduals.clicked.connect(self.togresiduals)
        self.loadParsButton.clicked.connect(self.openParFileDialog)
        self.addLineButton.clicked.connect(self.addLineDialog)
        self.writeParsButton.clicked.connect(self.writeParFileDialog)
        self.writeModelButton.clicked.connect(self.writeModelFileDialog)
        self.writeModelCompButton.clicked.connect(self.writeModelCompFileDialog)
        self.quitButton.clicked.connect(self.quitGui)

        ### Initialize spectral plots
        fig=Figure(figsize=(5,3))
        self.fig=fig
        self.initplot(fig)

        ### Initialize side plot
        sidefig=Figure(figsize=(5.85,3.75))
        self.sidefig = sidefig
        self.addsidempl(self.sidefig)
        self.sideplot(self.lastclick)  #Dummy initial cenwave setting

    def initplot(self,fig,numchunks=8):
        wlen=len(self.spectrum.wavelength)/numchunks
        self.spls=[]
        if self.wave1==None:  waveidx1=0  # Default to plotting entire spectrum for now
        else: waveidx1=jbg.closest(self.wave,self.wave1)
        if self.fitpars!=None:
                model=joebvpfit.voigtfunc(self.wave,self.datamodel.fitpars)
        sg=jbg.subplotgrid(numchunks)
        for i in range(numchunks):
            self.spls.append(fig.add_subplot(sg[i][0],sg[i][1],sg[i][2]))
            pixs=np.arange(waveidx1+i*wlen,waveidx1+(i+1)*wlen, dtype='int')
            self.spls[i].plot(self.wave[pixs],self.normflux[pixs],
                              linestyle='solid',linewidth=cfg.spec_linewidth)
            if self.fitpars!=None:
                self.spls[i].plot(self.wave,model,'r')
            self.spls[i].set_xlim(self.wave[pixs[0]],self.wave[pixs[-1]])
            self.spls[i].set_ylim(cfg.ylim)
            self.spls[i].set_xlabel('wavelength', fontsize=cfg.xy_fontsize,
                                    labelpad=cfg.x_labelpad)
            self.spls[i].set_ylabel('relative flux', fontsize=cfg.xy_fontsize,
                                    labelpad=cfg.y_labelpad)
            self.spls[i].get_xaxis().get_major_formatter().set_scientific(False)
            self.spls[i].tick_params(axis='both', which='major',direction='in',
                                     pad=2,length=2)
        fig.subplots_adjust(top=0.98,bottom=0.05,left=0.08,right=0.97,
                            wspace=0.15,hspace=0.24)
        self.addmpl(fig)

    def initialpars(self,parfilename):
        ### Deal with initial parameters from line input file
        self.fitpars,self.fiterrors,self.parinfo,self.linecmts = joebvpfit.readpars(parfilename)
        cfg.fitidx=joebvpfit.fitpix(self.wave, self.fitpars) #Set pixels for fit
        cfg.wavegroups=[]
        self.datamodel = LineParTableModel(self.fitpars,self.fiterrors,self.parinfo,linecmts=self.linecmts)
        self.tableView.setModel(self.datamodel)
        self.datamodel.updatedata(self.fitpars,self.fitpars,self.parinfo,self.linecmts)

    def sideplot(self,cenwave,wavebuf=3):
        if len(self.sidefig.axes)==0:
            self.sideax=self.sidefig.add_subplot(111)
        self.sideax.clear()
        self.sideax.plot(self.wave, self.normflux, linestyle='solid',
                         linewidth=cfg.spec_linewidth)
        if self.pixtog == 1:
            self.sideax.plot(self.wave[cfg.fitidx], self.normflux[cfg.fitidx], 'gs', markersize=4, mec='green')
        model = joebvpfit.voigtfunc(self.wave, self.fitpars)
        res = self.normflux - model
        self.sideax.plot(self.wave, model, 'r')
        if self.restog == 1:
            self.sideax.plot(self.wave, -res, '.', color='black',
                             ms=cfg.residual_markersize)
        self.sideax.plot(self.wave, [0] * len(self.wave), color='gray')
        self.sideax.set_xlabel('wavelength', fontsize=cfg.xy_fontsize,
                                labelpad=cfg.x_labelpad)
        self.sideax.set_ylabel('relative flux', fontsize=cfg.xy_fontsize,
                                labelpad=cfg.y_labelpad)
        self.sideax.tick_params(axis='both', which='major', direction='in',
                                pad=1, length=2)
        self.sidefig.subplots_adjust(top=0.98, bottom=0.17, left=0.12,
                                     right=0.97)

        ### label lines we are trying to fit
        if self.labeltog == 1:
            for j in range(len(self.fitpars[0])):
                labelloc = self.fitpars[0][j] * (1. + self.fitpars[3][j]) + self.fitpars[4][j] / c * \
                                                                            self.fitpars[0][j] * (
                                                                            1. + self.fitpars[3][j])
                label = ' {:.1f}_\nz{:.4f}'.format(self.fitpars[0][j], self.fitpars[3][j])
                self.sideax.text(labelloc, cfg.label_ypos, label, rotation=90, ha='center', va='bottom',
                                 clip_on=True, fontsize=cfg.label_fontsize)

        self.sideax.plot(self.wave, self.normsig, linestyle='solid', color='red', lw=0.5)
        self.sideax.plot(self.wave, -self.normsig, linestyle='solid', color='red', lw=0.5)
        self.sideax.get_xaxis().get_major_formatter().set_scientific(False)
        self.sideax.get_xaxis().get_major_formatter().set_useOffset(False)
        try:
            self.sideax.set_xlim(cenwave-wavebuf,cenwave+wavebuf)
            self.sideax.set_ylim(cfg.ylim)
            self.changesidefig(self.sidefig)
        except TypeError:
            pass
    def fitlines(self):
        print('VPmeasure: Fitting line profile(s)...')
        print(len(self.fitpars[0]),'lines loaded for fitting.')
        if self.fitconvtog:
            self.fitpars, self.fiterrors = joebvpfit.fit_to_convergence(self.wave, self.normflux, self.normsig,
                                                               self.datamodel.fitpars, self.datamodel.parinfo)
        else:
            self.fitpars, self.fiterrors = joebvpfit.joebvpfit(self.wave, self.normflux,self.normsig, self.datamodel.fitpars,self.datamodel.parinfo)
        self.datamodel.updatedata(self.fitpars,self.fiterrors,self.parinfo,self.linecmts)
        self.tableView.resizeColumnsToContents()
        self.updateplot()
        self.sideplot(self.lastclick)

    def togfitconv(self):
        if self.fitconvtog==1: self.fitconvtog=0
        else: self.fitconvtog=1

    def quitGui(self):
        self.deleteLater()
        self.close()

    def toglabels(self):
        if self.labeltog==1: self.labeltog=0
        else: self.labeltog=1
        self.updateplot()

    def togfitpix(self):
        if self.pixtog == 1:
            self.pixtog = 0
        else:
            self.pixtog = 1
        self.updateplot()

    def togresiduals(self):
        if self.restog == 1:
            self.restog = 0
        else:
            self.restog = 1
        self.updateplot()
        self.sideplot(self.lastclick)

    def openParFileDialog(self):
        fname = QtWidgets.QFileDialog.getOpenFileName(self, 'Open line parameter file','.')
        fname = str(fname[0])
        if fname != '':
            self.initialpars(fname)

        self.updateplot()

    def writeParFileDialog(self):
        fname = QtWidgets.QFileDialog.getSaveFileName(self, 'Save line parameter file', cfg.VPparoutfile)
        fname = str(fname[0])
        if fname != '':
            joebvpfit.writelinepars(self.datamodel.fitpars, self.datamodel.fiterrors, self.datamodel.parinfo, self.specfilename, fname, self.datamodel.linecmts)

    def writeModelFileDialog(self):
        fname = QtWidgets.QFileDialog.getSaveFileName(self, 'Save model to file', cfg.VPmodeloutfile)
        fname = str(fname[0])
        if fname != '':
            joebvpfit.writeVPmodel(fname, self.wave, self.fitpars, self.normflux, self.normsig)

    def writeModelCompFileDialog(self):
        dirDialog = QtWidgets.QFileDialog(self)
        dirDialog.setFileMode(dirDialog.Directory)
        dirDialog.setOption(dirDialog.ShowDirsOnly, True)
        defDirName = cfg.VPmodeloutfile[:-5]
        dname = dirDialog.getSaveFileName(self, 'Save model to files split by components',defDirName)
        dname = str(dname[0])
        if dname != '':
            joebvpfit.writeVPmodelByComp(dname, self.spectrum,self.fitpars)

    def addLineDialog(self):
        dlgOutput=newLineDialog.get_newline()
        if (dlgOutput != 0):
            if '' not in dlgOutput:
                newlam,newz,newcol,newb,newvel,newvel1,newvel2 = dlgOutput
                self.datamodel.addLine(self.wave,float(newlam), float(newz), float(newcol), float(newb), float(newvel), float(newvel1), float(newvel2))        #dialog=newLineDialog(parent=None)
                self.fitpars = self.datamodel.fitpars
                self.fiterrors = self.datamodel.fiterrors
                self.parinfo = self.datamodel.parinfo
                self.tableView.setModel(self.datamodel)


    def updateplot(self):
        if self.wave1==None:  waveidx1=0  # Default to plotting entire spectrum for now
        else: waveidx1=jbg.closest(self.wave,self.wave1)
        wlen=len(self.spectrum.wavelength)/self.numchunks
        for i,sp in enumerate(self.spls):
                sp.clear()
                prange=np.arange(waveidx1+i*wlen,waveidx1+(i+1)*wlen,dtype='int')
                if ((len(self.fitpars[0])>0)):
    
                    sp.plot(self.wave,self.normflux,linestyle='solid')
                    if self.pixtog==1:
                        sp.plot(self.wave[cfg.fitidx], self.normflux[cfg.fitidx], 'gs', markersize=4, mec='green')
                    model=joebvpfit.voigtfunc(self.wave,self.fitpars)
                    res=self.normflux-model
                    sp.plot(self.wave,model,'r')
                    if self.restog==1:
                        sp.plot(self.wave,-res,'.',color='black', ms=cfg.residual_markersize)
                    sp.plot(self.wave,[0]*len(self.wave),color='gray')
    
                    ### label lines we are trying to fit
                    if self.labeltog==1:
                        for j in range(len(self.fitpars[0])):
                            labelloc=self.fitpars[0][j]*(1.+self.fitpars[3][j])+self.fitpars[4][j]/c*self.fitpars[0][j]*(1.+self.fitpars[3][j])
                            label = ' {:.1f}_\nz{:.4f}'.format(self.fitpars[0][j], self.fitpars[3][j])
                            sp.text(labelloc, cfg.label_ypos, label, rotation=90, withdash=True, ha='center', va='bottom', clip_on=True, fontsize=cfg.label_fontsize)
                
    
                sp.plot(self.wave,self.normsig,linestyle='solid',color='red', lw=0.5)
                sp.plot(self.wave,-self.normsig,linestyle='solid',color='red', lw=0.5)
                sp.set_ylim(cfg.ylim)
                sp.set_xlim(self.wave[prange[0]],self.wave[prange[-1]])
                sp.set_xlabel('wavelength (A)', fontsize=cfg.xy_fontsize, labelpad=cfg.x_labelpad)
                sp.set_ylabel('normalized flux', fontsize=cfg.xy_fontsize, labelpad=cfg.y_labelpad)
                sp.get_xaxis().get_major_formatter().set_scientific(False)
                sp.get_xaxis().get_major_formatter().set_useOffset(False)
        self.changefig(self.fig)

    def changefig(self, item):
        #text = str(item.text())
        self.canvas.draw()
        #self.rmmpl()
        #self.addmpl(self.fig)


    def changesidefig(self, item):
        #text = str(item.text())
        self.sidecanvas.draw()
        #self.rmsidempl()
        #self.addsidempl(self.sidefig)
        
        
    def on_click(self, event):
        self.lastclick=event.xdata
        self.sideplot(self.lastclick)
        
    def addmpl(self, fig):
        self.canvas = FigureCanvas(fig)
        self.canvas.mpl_connect('button_press_event',self.on_click)
        self.mplvl.addWidget(self.canvas)
        self.canvas.draw()
        self.toolbar = NavigationToolbar(self.canvas,self.mplwindow,coordinates=True)
        self.mplvl.addWidget(self.toolbar)

    def addsidempl(self, sidefig):
        self.sidecanvas = FigureCanvas(sidefig)
        self.sidecanvas.setParent(self.sideMplWindow)
        if len(self.sidefig.axes) == 0:
            self.sidemplvl = QVBoxLayout()
        if len(self.sidefig.axes) != 0:
            self.sidemplvl.addWidget(self.sidecanvas)
        if len(self.sidefig.axes) == 0:
            self.sideMplWindow.setLayout(self.sidemplvl)
        self.sidecanvas.draw()

        
    def rmmpl(self,):
        self.canvas.close()
        self.mplvl.removeWidget(self.canvas)
        self.toolbar.close()
        self.mplvl.removeWidget(self.toolbar)

    def rmsidempl(self, ):
        self.sidemplvl.removeWidget(self.sidecanvas)
        self.sidecanvas.close()
        #self.sideMplWindow.setLayout(None)
        #self.mplvl.removeWidget(self.toolbar)
        #self.toolbar.close()

def go(specfilename, parfilename):
    import sys
    import numpy as np
    from astropy.io import fits as pf
    from astropy.io import ascii
    from linetools.spectra.io import readspec

    app = QtWidgets.QApplication.instance()
    if not app:
        app = QtWidgets.QApplication(sys.argv)
        app.aboutToQuit.connect(app.deleteLater)
    main = Main(specfilename, parfilename)
    main.show()
    app.exec_()

def batch_fit(spec, filelist, outparfile='.VP', outmodelfile='_VPmodel.fits', inspect=True, **kwargs):
    """
    Takes a number of input files and fits the lines in them.  The fitting algorithm will
    be run until convergence for each input file. The program will then ask whether to run the failed
    runs in GUI mode.

    Parameters
    ----------
    spec : string or XSpectrum1D
        The spectrum to be fitted with the input lines
    filelist : list of strings or str
        This should be a list containing the names of VP input files or a string referring to a file simply
        listing the input files.
        See joebvpfit.readpars for details of file format
    outparfile : str, optional
        Suffix for output file for fitted line parameters.
    outmodelfile: str, optional
        Suffix for output fits file for fitted model as the 'flux'.
    inspect : bool, optional
        Whether to produce plots and model for inspection of individual parameter input files.
        Default is True.


    **kwargs : maxiter, itertol
        These are fed on to the joebvp.fit_to_convergence() function

    """
    if isinstance(spec,str):
        spectofit = readspec(spec)
        specfile = spectofit.filename
    else:
        spectofit = spec
        specfile = cfg.filename

    if isinstance(filelist, str):
        lstarr=np.genfromtxt(filelist,dtype=None)
        listofiles=lstarr.tolist()
        if sum(1 for line in open(filelist)) == 1:
            listofiles = [listofiles]
    else:
        listofiles=filelist

    wave=spectofit.wavelength.value
    normflux=spectofit.flux.value/spectofit.co.value
    normsig=spectofit.sig.value/spectofit.co.value
    cfg.wave=wave
    cfg.spectrum = spectofit # need this for defining bad pixels later

    q_pass = 0
    q_fail = 0
    fails = [] # store failures filenames
    for i,ff in enumerate(listofiles):
        if isinstance(ff,bytes):
            ff = ff.decode()
        i += 1
        fitpars, fiterrors, parinfo, linecmts = joebvpfit.readpars(ff)

        cfg.lsfs = []
        cfg.fgs = []
        cfg.wavegroups = []
        cfg.wgidxs = []
        cfg.uqwgidxs = []

        try:
            fitpars,fiterrors=joebvpfit.fit_to_convergence(wave,normflux,normsig,fitpars,parinfo, **kwargs)
            print('VPmeasure: Fit converged:', ff)
            paroutfilename = ff.split('.')[0] + outparfile
            modeloutfilename = ff.split('.')[0] + outmodelfile
            joebvpfit.writelinepars(fitpars, fiterrors, parinfo, specfile, paroutfilename, linecmts)
            joebvpfit.writeVPmodel(modeloutfilename, wave, fitpars, normflux, normsig)
            if inspect:
                jbu.inspect_fits(paroutfilename, output=paroutfilename.split('.')[0]+"_inspect.pdf")
            q_pass += 1
        except:
            print('VPmeasure: Fitting failed:', ff)
            #import pdb; pdb.set_trace()
            fails += [ff]
            q_fail += 1
    print("")
    print("VPmeasure: {}/{} fits converged, {}/{} failed (see log for details).\n".format(q_pass, i, q_fail, i))

    # ask whether run GUI mode on failures
    if q_fail > 0:
        print("VPmeasure: Failed runs are: {}\n".format(fails))
        while True:
            answer = input("VPmeasure: Would you like to run the failed fits in GUI mode? (y/n): ")
            if answer in ['y','n', 'yes','no']:
                break
        if answer in ['y', 'yes']:
            for ff in fails:
                try:
                    go(spec, ff)
                except:
                    raise ValueError('Spectrum/VP input failed to load in interactive mode.')

    # Concatenate and create fits inspection files
    concatenate_all(spectofit)
    print("VPmeasure: Done.")


def concatenate_all(spectrum):
    """Takes all *.VP files in the working directory, and concatenates them into a single VP file
    as well as creates a single PDF file for fit inspection.

    spectrum : XSpectrum1D
        Original spectrum
    """
    # concatenate and inspect fit
    print("VPmeasure: concatenating individual outputs and creating figures for inspection.")
    os.system("ls *.VP > all_VP.txt")
    jbu.concatenate_line_tables("all_VP.txt")
    reload(cfg)  # Clear out the LSFs from the last fit
    cfg.spectrum = spectrum
    jbu.inspect_fits("compiledVPoutputs.dat")


if __name__ == '__main__':
        import sys
        specfilename=str(sys.argv[1])
        parfilename=str(sys.argv[5])
        go(specfilename,parfilename)
