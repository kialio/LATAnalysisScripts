import os
import math
import logging
from gt_apps import *

class quickAnalysis:

    def __init__(self,base='MySource',ra=0,dec=0,rad=10,tmin="INDEF",tmax="INDEF",emin=100,emax=300000,zmax=105,binned=False):
        self.base = base
        self.ra = ra
        self.dec = dec
        self.rad = rad
        self.tmin = tmin
        self.tmax = tmax
        self.emin = emin
        self.emax = emax
        self.zmax = zmax
        self.binned = binned

        self.logger = logging.getLogger('quickAnalysis')
        self.logger.setLevel(logging.DEBUG)
        fh = logging.FileHandler(self.base+'_quickAnalysis.log')
        fh.setLevel(logging.DEBUG)
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        self.logger.addHandler(fh)
        self.logger.addHandler(ch)
        self.logger.info("Created quickAnalysis object: base="+self.base+\
                             ",ra="+str(self.ra)+\
                             ",dec="+str(self.dec)+\
                             ",rad="+str(self.rad)+\
                             ",tmin="+str(self.tmin)+\
                             ",tmax="+str(self.tmax)+\
                             ",emin="+str(self.emin)+\
                             ",emax="+str(self.emax)+\
                             ",zmax="+str(self.zmax)+\
                             ",binned="+str(self.binned))
            
    def runCommand(self,AppCommand,run=True):

        if(run):
            AppCommand.run()
            self.logger.info(AppCommand.command())
        else:
            print AppCommand.command()
            

    def runSelect(self,run = True,evclsmin=3,evclsmax=4,convtype=-1):

        if(self.binned):
            filter['rad'] = self.rad * math.sqrt(2)
        else:
            filter['rad'] = self.rad

        filter['evclsmin'] = evclsmin
        filter['evclsmax'] = evclsmax
        filter['infile'] = "@"+self.base+".list"
        filter['outfile'] = self.base+"_filtered.fits"
        filter['ra'] = self.ra
        filter['dec'] = self.dec
        filter['tmin'] = self.tmin
        filter['tmax'] = self.tmax
        filter['emin'] = self.emin
        filter['emax'] = self.emax
        filter['zmax'] = self.zmax
        filter['convtype'] = convtype

        self.runCommand(filter,run)
        
    def runGTI(self, run = True, filterString="DATA_QUAL==1 && LAT_CONFIG==1 && ABS(ROCK_ANGLE)<52",roi = 'yes'):

        maketime['scfile'] = self.base+'_SC.fits'
        maketime['filter'] = filterString
        maketime['roicut'] = roi
        maketime['evfile'] = self.base+'_filtered.fits'
        maketime['outfile'] = self.base+'_filtered_gti.fits'

        self.runCommand(maketime,run)

    def runLTCube(self, run=True, zmax=180):

        expCube['evfile'] = self.base+'_filtered_gti.fits'
        expCube['scfile'] = self.base+'_SC.fits'
        expCube['outfile'] = self.base+'_ltcube.fits'
        expCube['dcostheta'] = 0.025
        expCube['binsz'] = 1
        expCube['zmax'] = zmax

        self.runCommand(expCube,run)

    def runExpMap(self, run=True,irfs="P6_V11_DIFFUSE"):

        expMap['evfile'] = self.base+'_filtered_gti.fits'
        expMap['scfile'] = self.base+'_SC.fits'
        expMap['expcube'] = self.base+'_ltcube.fits'
        expMap['outfile'] = self.base+'_expMap.fits'
        expMap['irfs'] = irfs
        expMap['srcrad'] = self.rad + 10.
        expMap['nlong'] = 120
        expMap['nlat'] = 120
        expMap['nenergies'] = 20

        self.runCommand(expMap,run)

    def runCCUBE(self, run=True,npix=100,nbins=30):

        bin_size = self.rad*2.0 / float(npix)

        evtbin['evfile'] = self.base+'_filtered_gti.fits'
        evtbin['outfile'] = self.base+'_CCUBE.fits'
        evtbin['algorithm'] = 'CCUBE'
        evtbin['nxpix'] = npix
        evtbin['nypix'] = npix
        evtbin['binsz'] = bin_size
        evtbin['coordsys'] = 'CEL'
        evtbin['xref'] = self.ra
        evtbin['yref'] = self.dec
        evtbin['axisrot'] = 0
        evtbin['proj'] = 'AIT'
        evtbin['ebinalg'] = 'LOG'
        evtbin['emin'] = self.emin
        evtbin['emax'] = self.emax
        evtbin['enumbins'] = nbins

        self.runCommand(evtbin,run)

    def runCMAP(self, run=True,npix=100):

        bin_size = self.rad*2.0 / float(npix)

        evtbin['evfile'] = self.base+'_filtered_gti.fits'
        evtbin['outfile'] = self.base+'_CMAP.fits'
        evtbin['algorithm'] = 'CMAP'
        evtbin['nxpix'] = npix
        evtbin['nypix'] = npix
        evtbin['binsz'] = bin_size
        evtbin['coordsys'] = 'CEL'
        evtbin['xref'] = self.ra
        evtbin['yref'] = self.dec
        evtbin['axisrot'] = 0
        evtbin['proj'] = 'AIT'
    
        self.runCommand(evtbin,run)

    def runExpCube(self,run=True,irfs="P6_V11_DIFFUSE"):

        cmd = "gtexpcube2 infile="+self.base+"_ltcube.fits cmap="\
                                  +self.base+"_CCUBE.fits outfile="\
                                  +self.base+"_BinnedExpMap.fits irfs="\
                                  +irfs

        if(run):
            os.system(cmd)
            self.logger.info(cmd)
        else:
            print cmd

    def runSrcMaps(self, run=True, irfs="P6_V11_DIFFUSE"):

        srcMaps['scfile'] = self.base+"_SC.fits"
        srcMaps['expcube'] = self.base+"_ltcube.fits"
        srcMaps['cmap'] = self.base+"_CCUBE.fits"
        srcMaps['srcmdl'] = self.base+"_model.xml"
        srcMaps['bexpmap'] = self.base+"_BinnedExpMap.fits"
        srcMaps['outfile'] = self.base+"_srcMaps.fits"
        srcMaps['irfs'] = irfs
        srcMaps['rfactor'] = 4
        srcMaps['emapbnds'] = "no"

        self.runCommand(srcMaps,run)

    def runModel(self,run=True, irfs="P6_V11_DIFFUSE"):
        
        model_map['srcmaps'] = self.base+"_srcMaps.fits"
        model_map['srcmdl'] = self.base+"_model.xml"
        model_map['outfile'] = self.base+"_modelMap.fits"
        model_map['expcube'] = self.base+"_ltcube.fits"
        model_map['irfs'] = irfs
        model_map['bexpmap'] = self.base+"_BinnedExpMap.fits"

        self.runCommand(model_map,run)

    def runAll(self, run=True, irfs="P6_V11_DIFFUSE"):

        self.logger.info("***Running gtselect***")
        self.runSelect(run)
        self.logger.info("***Running gtmktime***")
        self.runGTI(run)
        self.logger.info("***Running gtltcube***")
        self.runLTCube(run)

        if(self.binned):
            self.logger.info("***Running gtbin***")
            self.runCCUBE(run)
            self.logger.info("***Running gtexpcube2***")
            self.runExpCube(run,irfs)
            self.logger.info("***Running gtsrcMaps***")
            self.runSrcMaps(run,irfs)
        else:
            self.logger.info("***Running gtexpmap***")
            self.runExpMap()
                
