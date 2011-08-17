#!/usr/bin/env python

"""Perform event selections and exposure calculations for Fermi LAT data.

This module prepares Fermi LAT data for a likelihood anlaysis.  The
user supplies a text file with a list of data runs (one per line) and
a spacecraft file.  These should be called <basename>.list and
<basename>_SC.fits respectively where <basename> is a user defined
prefix.  There are two ways to run this module: from within python or
from the command line.  The simplest way to run it is from the command
line.

First, generate a default config file:

> quickAnalysis -c

The, edit the config file to match your specific analysis by filling
out the various options.  Next, run the command again:

> quickAnalysis <basename>

where <basename> is the prefix you've chosen to use; usually the name
of your source of interest but not necissarily so.

If you want to run this from within python, you'll need to first
create a quickAnalysis object and then you can use the various
functions below.  See the documentation for the individual functions
for more details.

This module logs all of the steps to a file called
<basename>_quickAnalysis.log as well as to the screen.

"""

__author__ = 'Jeremy S. Perkins (FSSC)'
__version__ = '0.1'

import sys
import os
import math
import logging
import ConfigParser
from gt_apps import *

class FileNotFound: pass

class quickAnalysis:

    """This is the base class.  If you want to use this, you first
    need to create an object from this method:

    >>> qA = quickAnalysis('example_name', configFile = True)

    will create an object called qA with the <basename> of
    'example_name' and will read in all of the options from the config
    file.  You can create an example config frile via the writeConfig
    function or via the command line with the -c option.  You can also
    pass all of the variables via the intial object initialiation
    function (see below).  Once you've created this object, you can
    just execute the runAll function to execute all of the steps, or
    use the functions individually as needed.
    """

    def __init__(self,
                 base = 'MySource',
                 configFile = False,
                 ra = 0,
                 dec = 0,
                 rad = 10,
                 tmin = "INDEF",
                 tmax = "INDEF",
                 emin = 100,
                 emax = 300000,
                 zmax = 100,
                 irfs = "P7SOURCE_V6",
                 binned = False,
                 verbosity = 0):

        self.base = base

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

        if(configFile):
            self.logger.info('Reading from config file')
            if(not os.path.exists(self.base+".cfg")):
                self.logger.criticial(self.base+".cfg doesn't exist.")
                return
            config = ConfigParser.RawConfigParser()
            config.read(self.base+'.cfg')
            ra = config.getfloat('quickAnalysis', 'ra')
            dec = config.getfloat('quickAnalysis', 'dec')
            rad = config.getfloat('quickAnalysis', 'rad')
            tmin = config.getfloat('quickAnalysis', 'tmin')
            tmax = config.getfloat('quickAnalysis', 'tmax')
            emin = config.getfloat('quickAnalysis', 'emin')
            emax = config.getfloat('quickAnalysis', 'emax')
            zmax = config.getfloat('quickAnalysis', 'zmax')

            irfs = config.get('common','irfs')
            binned = config.getboolean('common', 'binned')
            verbosity = config.getboolean('common', 'verbosity')

        self.ra = ra
        self.dec = dec
        self.rad = rad
        self.tmin = tmin
        self.tmax = tmax
        self.emin = emin
        self.emax = emax
        self.zmax = zmax
        self.binned = binned
        self.verbosity = verbosity
        self.irfs = irfs
        
        self.logger.info("Created quickAnalysis object: base="+self.base+\
                             ",ra="+str(self.ra)+\
                             ",dec="+str(self.dec)+\
                             ",rad="+str(self.rad)+\
                             ",tmin="+str(self.tmin)+\
                             ",tmax="+str(self.tmax)+\
                             ",emin="+str(self.emin)+\
                             ",emax="+str(self.emax)+\
                             ",zmax="+str(self.zmax)+\
                             ",irfs="+self.irfs+\
                             ",binned="+str(self.binned))
            
    def writeConfig(self):

        """Writes all of the initialization variables to the config
        file called <basename>.cfg."""


        config = ConfigParser.RawConfigParser()
        config.read(self.base+'.cfg')
        if(not config.has_section('common')):
            config.add_section('common')
        if(config.has_section('quickAnalysis')):
            self.logger.info("quickAnalysis config exists, overwriting...")
        else:
            config.add_section('quickAnalysis')

        config.set('common', 'base', self.base)
        config.set('common', 'binned', self.binned)
        config.set('common', 'verbosity', self.verbosity)
        config.set('common', 'irfs', self.irfs)
        
        config.set('quickAnalysis', 'ra', self.ra)
        config.set('quickAnalysis', 'dec', self.dec)
        config.set('quickAnalysis', 'rad', self.rad)
        config.set('quickAnalysis', 'tmin', self.tmin)
        config.set('quickAnalysis', 'tmax', self.tmax)
        config.set('quickAnalysis', 'emin', self.emin)
        config.set('quickAnalysis', 'emax', self.emax)
        config.set('quickAnalysis', 'zmax', self.zmax)
 
        with open(self.base+'.cfg', 'wb') as configfile:
            config.write(configfile)

    def runCommand(self,AppCommand,run=True):

        """Runs a giving command if run is True.  If run is False,
        prints out what the function would run."""

        if(run):
            AppCommand.run()
            self.logger.info(AppCommand.command())
        else:
            print AppCommand.command()
            
    def checkForFiles(self):
        """Checks for existence of needed files.  You need to have two
        files to start with: <basename>.list and <basename>_SC.fits.
        The first one is just a text file with a list of the raw data
        files (one per line) and the other is the spacecraft file.
        <basename> is a user defined prefix (usually the source name
        but not necissarily).  Returns an exception if any of the
        files are not found."""

        if(not os.path.exists(self.base+".list")):
            self.logger.critical(self.base+".list doesn't exist.")
            raise FileNotFound
            return
        if(not os.path.exists(self.base+"_SC.fits")):
            self.logger.critical(self.base+"_SC.fits doesn't exist.")
            raise FileNotFound
            return

    def runSelect(self,run = True,evclass=2,convtype=-1):

        """Runs gtselect on the data using the initialization
        parameters. User selected parameters include the conversion
        type and the eventclass."""

        filter['rad'] = self.rad
        filter['evclass'] = evclass
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

        """Executes gtmktime with the given filter"""

        maketime['scfile'] = self.base+'_SC.fits'
        maketime['filter'] = filterString
        maketime['roicut'] = roi
        maketime['evfile'] = self.base+'_filtered.fits'
        maketime['outfile'] = self.base+'_filtered_gti.fits'

        self.runCommand(maketime,run)

    def runLTCube(self, run=True, zmax=180):

        """Generates a livetime cube"""

        expCube['evfile'] = self.base+'_filtered_gti.fits'
        expCube['scfile'] = self.base+'_SC.fits'
        expCube['outfile'] = self.base+'_ltcube.fits'
        expCube['dcostheta'] = 0.025
        expCube['binsz'] = 1
        expCube['zmax'] = zmax

        self.runCommand(expCube,run)

    def runExpMap(self, run=True):

        """Generates an exposure map that is 10 degrees larger than
        the ROI and has 120 pixels in each direction."""

        expMap['evfile'] = self.base+'_filtered_gti.fits'
        expMap['scfile'] = self.base+'_SC.fits'
        expMap['expcube'] = self.base+'_ltcube.fits'
        expMap['outfile'] = self.base+'_expMap.fits'
        expMap['irfs'] = self.irfs
        expMap['srcrad'] = self.rad + 10.
        expMap['nlong'] = 120
        expMap['nlat'] = 120
        expMap['nenergies'] = 20

        self.runCommand(expMap,run)

    def runCCUBE(self, run=True,bin_size=0.1,nbins=30):

        """Generates a counts cube.  The dimensions of which are the
        largest square subtended by the ROI.  Note that if the ROI is
        exceptionally small or the bin size exceptionally large, the
        square might not be the largest posible since the npix
        calculation floors the calculated value.  The number of energy
        bins is logarithmic and is defined by the nbins variable."""

        npix = int((self.rad/math.sqrt(2.0))*2.0 / bin_size)

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

    def runCMAP(self, run=True,bin_size=0.1):
        
        """Generates a counts map.  The dimensions of which are the
        largest square subtended by the ROI.  Note that if the ROI is
        exceptionally small or the bin size exceptionally large, the
        square might not be the largest posible since the npix
        calculation floors the calculated value."""

        npix = int((self.rad/math.sqrt(2.0))*2.0 / bin_size)

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

    def runExpCube(self,run=True,bin_size=0.1,nbins=30):

        """Generates a binned exposure map that is 20 degrees larger
        than the ROI.  The binned exposure map needs to take into
        account the exposure on sources outside of the ROI.  20
        degrees is the size of the PSF at low energies plus an extra
        10 degrees for security.  The energy binning is logarithmic
        and the number of energy bins is defined by the nbins
        variable."""

        npix = int(((self.rad+20.)/math.sqrt(2.0))*2.0 / bin_size)

        cmd = "gtexpcube2 infile="+self.base+"_ltcube.fits"\
            +" cmap=none"\
            +" outfile="+self.base+"_BinnedExpMap.fits"\
            +" irfs="+self.irfs\
            +" xref="+str(self.ra)\
            +" yref="+str(self.dec)\
            +" nxpix="+str(npix)\
            +" nypix="+str(npix)\
            +" binsz="+str(bin_size)\
            +" coordsys=CEL"\
            +" axisrot=0"\
            +" proj=AIT"\
            +" ebinalg=LOG"\
            +" emin="+str(self.emin)\
            +" emax="+str(self.emax)\
            +" enumbins="+str(nbins)
            
        if(run):
            os.system(cmd)
            self.logger.info(cmd)
        else:
            print cmd

    def generateXMLmodel(self, 
                         galactic_file="gal_2yearp7v6_v0.fits", 
                         isotropic_file="iso_p7v6source.txt", 
                         catalog_file="gll_psc_v05.fit"):

        """Checks to see if <basename>_model.xml exists and creates
        one using the make2FGLXML module if it doesn't.
        make2FGLXml.py must be in the working directory or accessable
        to python for this function to work.  The galactic and
        isotropic models plus the Fermi LAT catalog must also be in
        the working directory.  Additionally, if any extended sources
        are in the ROI, the diffuse templates for those sources should
        be in the working directory."""

        if(os.path.exists(self.base+"_model.xml")):
            self.logger.info(self.base+"_model.xml exists, won't create a new one.")
        else:
            self.logger.info(self.base+"_model.xml doesn't exists, will create a new one.")
            import make2FGLxml
        
            mymodel = make2FGLxml.srcList(catalog_file,self.base+"_filtered_gti.fits",self.base+"_model.xml")
            mymodel.makeModel(galactic_file, 'gal_2yearp7v6_v0', isotropic_file, 'iso_p7v6source')

            self.logger.info("NOTE: if there are extended sources in your ROI, make sure the "\
                                 +"correspoinding diffuse template is in the working directory.")
        

    def runSrcMaps(self, run=True):

        """Generates a source map for your region."""

        self.generateXMLmodel()

        srcMaps['scfile'] = self.base+"_SC.fits"
        srcMaps['expcube'] = self.base+"_ltcube.fits"
        srcMaps['cmap'] = self.base+"_CCUBE.fits"
        srcMaps['srcmdl'] = self.base+"_model.xml"
        srcMaps['bexpmap'] = self.base+"_BinnedExpMap.fits"
        srcMaps['outfile'] = self.base+"_srcMaps.fits"
        srcMaps['irfs'] = self.irfs
        srcMaps['rfactor'] = 4
        srcMaps['emapbnds'] = "no"

        self.runCommand(srcMaps,run)

    def runModel(self,run=True):

        """Generates a model map."""
        
        model_map['srcmaps'] = self.base+"_srcMaps.fits"
        model_map['srcmdl'] = self.base+"_model.xml"
        model_map['outfile'] = self.base+"_modelMap.fits"
        model_map['expcube'] = self.base+"_ltcube.fits"
        model_map['irfs'] = self.irfs
        model_map['bexpmap'] = self.base+"_BinnedExpMap.fits"

        self.runCommand(model_map,run)

    def runAll(self, run=True):

        """Does a full event selection and exposure calculation.  This
        is the function called when this module is run from the
        command line."""

        self.logger.info("***Checking for files***")
        try:
            self.checkForFiles()
        except(FileNotFound):
            self.logger.critical("One or more needed files do not exist")
            return
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
            self.runExpCube(run)
            self.logger.info("***Running gtsrcMaps***")
            self.runSrcMaps(run)
        else:
            self.logger.info("***Running gtexpmap***")
            self.runExpMap(run)
   
# Command-line interface             
def cli():
    """Command-line interface.  Call this without any options for usage notes."""
    import getopt
    class BadUsage: pass
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'c')
        
        for opt, val in opts:
            if opt == '-c':
                print "Creating example configuration file called example.cfg"
                qA = quickAnalysis("example")
                qA.writeConfig()
                return

        if not args: raise BadUsage
        for arg in args:
            qA = quickAnalysis(arg, True)
            qA.runAll(False)

    except (getopt.error, BadUsage):
        cmd = os.path.basename(sys.argv[0])
        print """quickAnalysis - Perform event selections and exposure calculations for Fermi LAT data.

%s <basename> ...  Perform an analysis on <basename>.  <basename> is
    the prefix used for this analysis.  You must already have a
    configuration file if using the command line interface.
""" %(cmd)

if __name__ == '__main__': cli()
