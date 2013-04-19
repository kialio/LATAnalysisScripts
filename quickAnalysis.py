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

> quickAnalysis (-i|--initialize)

Then edit the config file to match your specific analysis by filling
out the various options.  Next, run the command again:

> quickAnalysis -a (-n |--basename==)<basename>

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
__version__ = '0.1.11'

import sys
import os
from gt_apps import filter, maketime, expMap, expCube, evtbin, srcMaps
import quickUtils as qU

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
                 base='MySource',
                 configFile = False,
                 analysisConfig = {"ra" : 0,
                                   "dec" : 0,
                                   "rad" : 10,
                                   "tmin" : "INDEF",
                                   "tmax" : "INDEF",
                                   "emin" : 100,
                                   "emax" : 300000,
                                   "zmax" : 100,
                                   "binsize" : 0.1,
                                   "convtype" : -1},
                 commonConfig = {"base" : 'MySource',
                                 "binned" : False,
                                 "eventclass" : 2,
                                 "irfs" : "P7SOURCE_V6",
                                 "verbosity" : 0,
                                 "multicore": 0}):

        commonConfig['base'] = base

        self.logger = qU.initLogger(base, 'quickAnalysis')

        if(configFile):
            try:
                commonConfigRead,analysisConfigRead,likelihoodConfigRead,plotConfigRead,curveConfigRead = qU.readConfig(self.logger,base)
            except(qU.FileNotFound):
                self.logger.critical("One or more needed files do not exist")
                return
            try:
                commonConfig = qU.checkConfig(self.logger,commonConfig,commonConfigRead)
            except(KeyError):
                return
            try:
                analysisConfig = qU.checkConfig(self.logger,analysisConfig,analysisConfigRead)
            except(KeyError):
                return

        self.commonConf = commonConfig
        self.analysisConf = analysisConfig

        logString = "Created quickAnalysis object: "
        for variable, value in commonConfig.iteritems():
            logString += variable+"="+str(value)+","
        for variable, value in analysisConfig.iteritems():
            logString += variable+"="+str(value)+","
        self.logger.info(logString)
            
    def writeConfig(self):

        """Writes all of the initialization variables to the config
        file called <basename>.cfg."""

        qU.writeConfig(quickLogger=self.logger,
                       commonDictionary=self.commonConf,
                       analysisDictionary=self.analysisConf)

    def runSelect(self,run = True):

        """Runs gtselect on the data using the initialization
        parameters. User selected parameters include the conversion
        type and the eventclass."""

        filter['rad'] = self.analysisConf['rad']
        filter['evclass'] = self.commonConf['eventclass']
        filter['evclsmin'] = 0
        filter['evclsmax'] = 10
        filter['infile'] = "@"+self.commonConf['base']+".list"
        filter['outfile'] = self.commonConf['base']+"_filtered.fits"
        filter['ra'] = self.analysisConf['ra']
        filter['dec'] = self.analysisConf['dec']
        filter['tmin'] = self.analysisConf['tmin']
        filter['tmax'] = self.analysisConf['tmax']
        filter['emin'] = self.analysisConf['emin']
        filter['emax'] = self.analysisConf['emax']
        filter['zmax'] = self.analysisConf['zmax']
        filter['convtype'] = self.analysisConf['convtype']

        qU.runCommand(filter,self.logger,run)
        
    def runGTI(self, run = True, filterString="DATA_QUAL==1 && LAT_CONFIG==1 && ABS(ROCK_ANGLE)<52",roi = 'yes'):

        """Executes gtmktime with the given filter"""

        maketime['scfile'] = self.commonConf['base']+'_SC.fits'
        maketime['filter'] = filterString
        maketime['roicut'] = roi
        maketime['evfile'] = self.commonConf['base']+'_filtered.fits'
        maketime['outfile'] = self.commonConf['base']+'_filtered_gti.fits'

        qU.runCommand(maketime,self.logger,run)

    def runLTCube(self, run=True, zmax=180):

        """Generates a livetime cube"""

        expCube['evfile'] = self.commonConf['base']+'_filtered_gti.fits'
        expCube['scfile'] = self.commonConf['base']+'_SC.fits'
        expCube['outfile'] = self.commonConf['base']+'_ltcube.fits'
        expCube['dcostheta'] = 0.025
        expCube['binsz'] = 1
        expCube['zmax'] = zmax

        np = int(self.commonConf['multicore'])
        
        if(np > 0):
            import gtltcube_mp as ltmp
            ltmp.gtltcube_mp(np, expCube['scfile'], expCube['evfile'], expCube['outfile'],False,expCube['zmax'])
        else:
            qU.runCommand(expCube,self.logger,run)

    def runExpMap(self, run=True):

        """Generates an exposure map that is 10 degrees larger than
        the ROI and has 120 pixels in each direction."""

        expMap['evfile'] = self.commonConf['base']+'_filtered_gti.fits'
        expMap['scfile'] = self.commonConf['base']+'_SC.fits'
        expMap['expcube'] = self.commonConf['base']+'_ltcube.fits'
        expMap['outfile'] = self.commonConf['base']+'_expMap.fits'
        expMap['irfs'] = self.commonConf['irfs']
        expMap['srcrad'] = float(self.analysisConf['rad']) + 10.
        expMap['nlong'] = 120
        expMap['nlat'] = 120
        expMap['nenergies'] = 20

        np = int(self.commonConf['multicore'])
        
        if(np > 0):
            import gtexpmap_mp as expmp
            expmp.gtexpmap_mp(expMap['nlong'],expMap['nlat'],np,1,
                              expMap['scfile'],expMap['evfile'],expMap['expcube'], 
                              expMap['irfs'],expMap['srcrad'],expMap['nenergies'],
                              expMap['outfile'],False)
        else:
            qU.runCommand(expMap,self.logger,run)

    def runCCUBE(self, run=True,nbins=30):

        """Generates a counts cube.  The dimensions of which are the
        largest square subtended by the ROI.  Note that if the ROI is
        exceptionally small or the bin size exceptionally large, the
        square might not be the largest posible since the npix
        calculation floors the calculated value.  The number of energy
        bins is logarithmic and is defined by the nbins variable."""

        npix = qU.NumberOfPixels(float(self.analysisConf['rad']), float(self.analysisConf['binsize']))

        evtbin['evfile'] = self.commonConf['base']+'_filtered_gti.fits'
        evtbin['outfile'] = self.commonConf['base']+'_CCUBE.fits'
        evtbin['algorithm'] = 'CCUBE'
        evtbin['nxpix'] = npix
        evtbin['nypix'] = npix
        evtbin['binsz'] = self.analysisConf['binsize']
        evtbin['coordsys'] = 'CEL'
        evtbin['xref'] = self.analysisConf['ra']
        evtbin['yref'] = self.analysisConf['dec']
        evtbin['axisrot'] = 0
        evtbin['proj'] = 'AIT'
        evtbin['ebinalg'] = 'LOG'
        evtbin['emin'] = self.analysisConf['emin']
        evtbin['emax'] = self.analysisConf['emax']
        evtbin['enumbins'] = nbins
        evtbin['scfile'] = self.commonConf['base']+"_SC.fits"

        qU.runCommand(evtbin,self.logger,run)

    def runCMAP(self, run=True):
        
        """Generates a counts map.  The dimensions of which are the
        largest square subtended by the ROI.  Note that if the ROI is
        exceptionally small or the bin size exceptionally large, the
        square might not be the largest posible since the npix
        calculation floors the calculated value."""

        qU.runCMAP(self.logger, 
                   self.commonConf['base'],
                   self.analysisConf['rad'],
                   self.analysisConf['binsize'],
                   self.analysisConf['ra'],
                   self.analysisConf['dec'])

    def runExpCube(self,run=True,nbins=30, ExpBuffer=30):

        """Generates a binned exposure map that is 30 degrees larger
        than the ROI.  The binned exposure map needs to take into
        account the exposure on sources outside of the ROI.  30
        degrees is the size of the PSF at low energies plus some extra
        tens of degrees for security.  You can adjust this buffer by
        setting the ExpBuffer option.  The energy binning is
        logarithmic and the number of energy bins is defined by the
        nbins variable."""

        npix = qU.NumberOfPixels(float(self.analysisConf['rad'])+ExpBuffer, float(self.analysisConf['binsize']))

        cmd = "gtexpcube2 infile="+self.commonConf['base']+"_ltcube.fits"\
            +" cmap=none"\
            +" outfile="+self.commonConf['base']+"_BinnedExpMap.fits"\
            +" irfs="+self.commonConf['irfs']\
            +" xref="+str(self.analysisConf['ra'])\
            +" yref="+str(self.analysisConf['dec'])\
            +" nxpix="+str(npix)\
            +" nypix="+str(npix)\
            +" binsz="+str(self.analysisConf['binsize'])\
            +" coordsys=CEL"\
            +" axisrot=0"\
            +" proj=AIT"\
            +" ebinalg=LOG"\
            +" emin="+str(self.analysisConf['emin'])\
            +" emax="+str(self.analysisConf['emax'])\
            +" enumbins="+str(nbins)
            
        if(run):
            os.system(cmd)
            self.logger.info(cmd)
        else:
            print cmd

    def generateXMLmodel(self):
        
        """Calls the quickUtils function to make an XML model of your
        region based on the 2FGL. make2FGLXml.py needs to be in your
        python path.  This needs to have the galactic and isotropic
        diffuse models in your working directory as well as the 2FGL
        catalog in FITS format.  See the corresponding function in
        quickUtils for more details."""
        
        try:
            qU.generateXMLmodel(self.logger, self.commonConf['base'])
        except(qU.FileNotFound):
            self.logger.critical("One or more needed files do not exist")
            return

    def runSrcMaps(self, run=True):

        """Generates a source map for your region.  Checks to make
        sure that there's an XML model to be had and if not, creates
        one from the 2FGL."""

        self.generateXMLmodel()

        try:
            qU.checkForFiles(self.logger,[self.commonConf['base']+"_CCUBE.fits",
                                       self.commonConf['base']+"_ltcube.fits",
                                       self.commonConf['base']+"_BinnedExpMap.fits",
                                       self.commonConf['base']+"_SC.fits",])
        except(qU.FileNotFound):
            self.logger.critical("One or more needed files do not exist")
            return

        srcMaps['scfile'] = self.commonConf['base']+"_SC.fits"
        srcMaps['expcube'] = self.commonConf['base']+"_ltcube.fits"
        srcMaps['cmap'] = self.commonConf['base']+"_CCUBE.fits"
        srcMaps['srcmdl'] = self.commonConf['base']+"_model.xml"
        srcMaps['bexpmap'] = self.commonConf['base']+"_BinnedExpMap.fits"
        srcMaps['outfile'] = self.commonConf['base']+"_srcMaps.fits"
        srcMaps['irfs'] = self.commonConf['irfs']
        srcMaps['rfactor'] = 4
        srcMaps['emapbnds'] = "yes"

        qU.runCommand(srcMaps,self.logger,run)

    def runModel(self, run=True):

        """Wrapper for the same function in quickUtils"""

        qU.runModel(self.logger,self.commonConf['base'],self.commonConf['irfs'],run)

    def runAll(self, run=True):

        """Does a full event selection and exposure calculation.  This
        is the function called when this module is run from the
        command line.  You need to have two files to start with:
        <basename>.list and <basename>_SC.fits.  The first one is just
        a text file with a list of the raw data files (one per line)
        and the other is the spacecraft file.  <basename> is a user
        defined prefix (usually the source name but not necissarily).
        Returns an exception if any of the files are not found."""

        self.logger.info("***Checking for files***")
        try:
            qU.checkForFiles(self.logger,[self.commonConf['base']+".list",self.commonConf['base']+"_SC.fits"])
        except(qU.FileNotFound):
            self.logger.critical("One or more needed files do not exist")
            return
        self.logger.info("***Running gtselect***")
        self.runSelect(run)
        self.logger.info("***Running gtmktime***")
        self.runGTI(run)
        self.logger.info("***Running gtltcube***")
        self.runLTCube(run)

        if(self.commonConf['binned']):
            self.logger.info("***Running gtbin***")
            self.runCCUBE(run)
            self.logger.info("***Running gtexpcube2***")
            self.runExpCube(run)
            self.logger.info("***Running gtsrcMaps***")
            self.runSrcMaps(run)
        else:
            self.logger.info("***Running gtexpmap***")
            self.runExpMap(run)

def printCLIHelp():
    """This function prints out the help for the CLI."""

    cmd = os.path.basename(sys.argv[0])
    print """
                        - quickAnalysis - 

Perform event selections and exposure calculations for Fermi LAT data.
You can use the command line functions listed below or run this module
from within python. For full documentation on this module execute
'pydoc quickAnalysis'.

%s (-h|--help) ... This help text.

%s (-a|--analyze) (-n |--basename=)<basename> ...  Perform an analysis
    on <basename>.  <basename> is the prefix used for this analysis.
    You must already have a configuration file if using the command
    line interface.

%s (-i|--initialize) ... Generate a default config file called
    example.cfg.  Edit this file and rename it <basename>.cfg for use
    in the quickAnalysis module.

%s (-x|--xml) (-n |--basename=)<basename> ... Generate a model file
    from the 2FGL.  You need to already have
    <basename>_filtered_gti.fits in your working directory.  You can
    get this file by running the functions runSelect and runGTI on
    your data.  You also need to have the galactic and isotropic
    diffuse models in your working directory as well as the 2FGL model
    file.

%s (-m|--model) (-n |--basename=)<basename> ... Generate a model map
    based on the model file in your config file.  You need to have
    several files already computed.  It's best to do the runAll script
    before trying this.

%s (-c|--cmap) (-n |--basename=)<basename> {(-b |--binsize=)<binsz>} ...
    Generate a counts map for the specific analysis defined by
    <basename>.  You must have generated several files already.  It's
    best to do the runAll script (or execute this script without any
    options) first.  You can give it a bin size (in deg/bin) if you
    want.

%s --filter (-n | --basename=)<basename> ... Generated a file that has
    been event and GTI selected.

%s --sourcemap (-n | --basename=)<basename> ... Generate a source map
    based on the model file in your config file.  You need to have
    already produced several other files.

%s --ccube (-n | --basename=)<basename> ... Generate a counts cube for
    your analysis.  You need to have already produced several other
    files.

%s --bexpmap (-n | --basename=)<basename> ... Generate a binned
    exposure map for your analysis.  You need to have already produced
    several other files.

""" %(cmd,cmd,cmd,cmd,cmd,cmd,cmd,cmd,cmd,cmd)

# Command-line interface             
def cli():
    """Command-line interface.  Call this without any options for usage notes."""
    import getopt
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hiamxcb:n:', ['help',
                                                                'analyze',
                                                                'initialize',
                                                                'modelmap',
                                                                'sourcemap',
                                                                'ccube',
                                                                'bexpmap',
                                                                'filter',
                                                                'xml',
                                                                'cmap',
                                                                'basename=',
                                                                'binsize=',
                                                                ])
        #Loop through first and check for the basename
        haveBase = False
        basename = 'example'
        for opt,val in opts:
            if opt in ('-n','--basename'):
                haveBase = True
                basename = val

        for opt, val in opts:
            if opt in ('-h', '--help'):
                printCLIHelp()
                return
            elif opt in ('-a', '--analyze'):
                if not haveBase: raise getopt.GetoptError("Must specify basename, printing help.")
                qA = quickAnalysis(basename, True)
                qA.runAll(True)         
                return
            elif opt in ('-i', '--initialize'):
                print "Creating example configuration file called example.cfg"
                qA = quickAnalysis(basename)
                qA.writeConfig()
                return
            elif opt in ('-m','--modelmap'):
                if not haveBase: raise getopt.GetoptError("Must specify basename, printing help.")
                print "Creating model map"
                qA = quickAnalysis(basename, True)
                qU.runModel(qA.logger,qA.commonConf['base'],qA.commonConf['irfs'])
                return
            elif opt in ('--sourcemap'):
                if not haveBase: raise getopt.GetoptError("Must specify basename, printing help.")
                print "Creating source map"
                qA = quickAnalysis(basename, True)
                qA.runSrcMaps()
                return
            elif opt in ('--ccube'):
                if not haveBase: raise getopt.GetoptError("Must specify basename, printing help.")
                print "Creating counts cube"
                qA = quickAnalysis(basename, True)
                qA.runCCUBE()
                return
            elif opt in ('--bexpmap'):
                if not haveBase: raise getopt.GetoptError("Must specify basename, printing help.")
                print "Creating binned exposure map"
                qA = quickAnalysis(basename, True)
                qA.runExpCube()
                return
            elif opt in ('--filter'):
                if not haveBase: raise getopt.GetoptError("Must specify basename, printing help.")
                print "Creating filtered (selected and GTI'd) file"
                qA = quickAnalysis(basename, True)
                qA.runSelect()
                qA.runGTI()
                return
            elif opt in ('-x', '--xml'):
                if not haveBase: raise getopt.GetoptError("Must specify basename, printing help.")
                print "Creating XML model file from 2FGL"
                qA = quickAnalysis(basename, True)
                qA.generateXMLmodel()
                return
            elif opt in ('-c', '--cmap'):
                if not haveBase: raise getopt.GetoptError("Must specify basename, printing help.")
                print "Creating counts map"
                qA = quickAnalysis(basename, True)
                for option,value in opts:
                    if option in ('-b', '--binsize'): qA.analysisConf['binsize'] = value
                qA.runCMAP()
                return
                
        if not opts: raise getopt.GetoptError("Must specify an option, printing help.")

    except getopt.error as e:
        print "Command Line Error: " + e.msg
        printCLIHelp()


if __name__ == '__main__': cli()
