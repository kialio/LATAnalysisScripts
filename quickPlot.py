#!/usr/bin/env python

"""Plots DS9 field of view of the Fermi LAT data analysis results.

This module generates the residual and significance maps starting
from the count and model maps (created and used in the likelihood 
analysis of the Fermi LAT data). All the maps are then plotted using 
DS9 for an easy comparison.

There are two ways to run this module: from within python or from 
the command line.  The simplest way to run it is from the command line.

First, generate a default config file:

> quickPlot -c

Then edit the config file to match your specific analysis by filling
out the various options.  Next, run the command again:

> quickPlot <basename>

where <basename> is the prefix you've chosen to use; usually the name
of your source of interest but not necessarily so.

If you want to run this from within python, you'll need to first
create a quickPlot object and then you can use the various
functions below.  See the documentation for the individual functions
for more details.

This module logs all of the steps to a file called
<basename>_quickPlot.log as well as to the screen.

"""

__author__ = 'Davide Donato (FSSC)'
__version__ = '0.1.1'

import sys
import os
import pyfits
from math import *
from ds9 import *
from quickUtils import *

class quickPlot:

    """This is the base class.  If you want to use this, you first
    need to create an object from this method:

    >>> qP = quickPlot('example_name', configFile = True)

    will create an object called qP with the <basename> of
    'example_name' and will read in all of the options from the config
    file.  You can create an example config file via the writeConfig
    function or via the command line with the -c option.  You can also
    pass all of the variables via the intial object initialiation
    function (see below).  Once you've created this object, you can
    just execute the runAll function to execute all of the steps, or
    use the functions individually as needed.
    """

    def __init__(self,
                 base='MySource',
                 configFile = False,
                 plotConfig = {"scaletypeframe1"   : "log",
                               "scalemodeframe1"    : "minmax",
			       "colorframe1"       : "aips0",
                               "scaletypeframe2"   : "log",
                               "scalemodeframe2"    : "minmax",
			       "colorframe2"       : "a",
                               "scaletypeframe3"   : "sqrt",
                               "scalemodeframe3"    : "zscale",
			       "colorframe3"       : "aips0",
                               "scaletypeframe4"   : "sqrt",
                               "scalemodeframe4"    : "minmax",
			       "colorframe4"       : "aips0",
                               "sourceregioncolor" : "black",
                               "sourceregionwidth" : 3,
                               "sourceregionfont"  : "helvetica 14 bold",
                               "sourceregiontype"  : "cross 15",
                               "labelcolor"        : "black",
                               "labelwidth"        : 3,
                               "labelfont"         : "helvetica 10 bold",
                               "labelposition"     : "40 15",
                               "grid"              : "yes",
                               "gridcolor"         : "black",
                               "gridfont"          : 14},
                 likelihoodConfig = {"model"      : "MySource_model.xml",
                                     "sourcename" : "Source Name"},
                 commonConfig = {"base"       : 'MySource',
                                 "binned"     : False,
                                 "eventclass" : 2, 
                                 "irfs"       : "P7SOURCE_V6",
                                 "verbosity"  : 0}):
 
        commonConfig['base'] = base

        self.logger = initLogger(base, 'quickPlot')

        if(configFile):
            try:
                commonConfigRead,analysisConfigRead,likelihoodConfigRead,plotConfigRead = readConfig(self.logger,base)
            except(FileNotFound):
                self.logger.critical("One or more needed files do not exist")
                return
            try:
                commonConfig = checkConfig(self.logger,commonConfig,commonConfigRead)
            except(KeyError):
                return
            try:
                plotConfig = checkConfig(self.logger,plotConfig,plotConfigRead)
            except(KeyError):
                return
            try:
                likelihoodConfig = checkConfig(self.logger,likelihoodConfig,likelihoodConfigRead)
            except(KeyError):
                return
                
        self.commonConf = commonConfig
        self.plotConf = plotConfig
        self.likelihoodConf = likelihoodConfig

        logString = "Created quickPlot object: "
        for variable, value in commonConfig.iteritems():
            logString += variable+"="+str(value)+","
        for variable, value in plotConfig.iteritems():
            logString += variable+"="+str(value)+","
        self.logger.info(logString)
            
	    
    def writeConfig(self):

        """Writes all of the initialization variables to the config
        file called <basename>.cfg."""

        writeConfig(quickLogger=self.logger,
                    commonDictionary=self.commonConf,
                    plotDictionary=self.plotConf)


    def createModelMap(self,run=True):

        """Wrapper for the model map routine in quickUtils"""

        runModel(self.logger,self.commonConf['base'],self.commonConf['irfs'],run)


    def createResidMap(self,run=True):

        """Generates a residual map"""

        try:
            checkForFiles(self.logger,
	    		 [self.commonConf['base']+"_CMap.fits",
			  self.commonConf['base']+"_modelMap.fits"])
        except(FileNotFound):
            self.logger.critical("One or more needed files do not exist")
            return

        onImage  = pyfits.open(self.commonConf['base']+"_CMAP.fits")
        onData   = onImage[0].data.copy()
        onHeader = onImage[0].header
        offImage = pyfits.open(self.commonConf['base']+"_modelMap.fits")
        offData  = offImage[0].data.copy()
        resData  = offImage[0].data.copy()

        for x,row in enumerate(resData):
            for y in enumerate(row):
                resData[x,y[0]] = (onData[x,y[0]]-offData[x,y[0]])

        newImage = pyfits.PrimaryHDU(resData)
        newImage.header = onHeader
        newImage.update_header()

        hdulist = pyfits.HDUList([newImage])
        hdulist.writeto(self.commonConf['base']+"_resMap.fits",clobber=True)        
             

    def createSigMap(self,run=True):

        """Generates a significance map."""

        try:
            checkForFiles(self.logger,
	    		 [self.commonConf['base']+"_CMap.fits",
			  self.commonConf['base']+"_modelMap.fits"])
        except(FileNotFound):
            self.logger.critical("One or more needed files do not exist")
            return

        onImage  = pyfits.open(self.commonConf['base']+"_CMAP.fits")
        onData   = onImage[0].data.copy()
        onHeader = onImage[0].header
        offImage = pyfits.open(self.commonConf['base']+"_modelMap.fits")
        offData  = offImage[0].data.copy()
        sigData  = offImage[0].data.copy()

        for x,row in enumerate(sigData):
            for y in enumerate(row):
                sigData[x,y[0]] = (onData[x,y[0]]-offData[x,y[0]])/sqrt(onData[x,y[0]]+offData[x,y[0]])

        newImage = pyfits.PrimaryHDU(sigData)
        newImage.header = onHeader
        newImage.update_header()

        hdulist = pyfits.HDUList([newImage])
        hdulist.writeto(self.commonConf['base']+"_sigMap.fits",clobber=True)

 
    def plotMaps(self,run=True):

	""""Uses ds9 to plot the count, model, residual and significance maps"""

        try:
            checkForFiles(self.logger,
	    		 [self.commonConf['base']+"_CMap.fits",
			  self.commonConf['base']+"_modelMap.fits",
			  self.commonConf['base']+"_residMap.fits",
			  self.commonConf['base']+"_sigMap.fits",
			  self.commonConf['base']+"_model.xml"])
        except(FileNotFound):
            self.logger.critical("One or more needed files do not exist")
            return

        d = ds9()
	
        d.set('file '+self.commonConf['base']+"_CMAP.fits")
        d.set('scale '+self.plotConf['scaletypeframe1'])
        d.set('scale mode '+self.plotConf['scalemodeframe1'])
        d.set('cmap '+self.plotConf['colorframe1'])
        d.set('regions', 'image; text '+self.plotConf['labelposition']+' # color='+self.plotConf['labelcolor']+' width='+self.plotConf['labelwidth']+' font="'+self.plotConf['labelfont']+'" text={Count map}')

        d.set('tile')

        d.set('frame new')
        d.set('file '+self.commonConf['base']+"_modelMap.fits")
        d.set('scale '+self.plotConf['scaletypeframe2'])
        d.set('scale mode '+self.plotConf['scalemodeframe2'])
        d.set('cmap '+self.plotConf['colorframe2'])
        d.set('regions', 'image; text '+self.plotConf['labelposition']+' # color='+self.plotConf['labelcolor']+' width='+self.plotConf['labelwidth']+' font="'+self.plotConf['labelfont']+'" text={Model map}')

        d.set('frame new')
        d.set('file '+self.commonConf['base']+"_residMap.fits")
        d.set('scale '+self.plotConf['scaletypeframe3'])
        d.set('scale mode '+self.plotConf['scalemodeframe3'])
        d.set('cmap '+self.plotConf['colorframe3'])
        d.set('regions', 'image; text '+self.plotConf['labelposition']+' # color='+self.plotConf['labelcolor']+' width='+self.plotConf['labelwidth']+' font="'+self.plotConf['labelfont']+'" text={Residual}')

        d.set('frame new')
        d.set('file '+self.commonConf['base']+"_sigMap.fits")
        d.set('scale '+self.plotConf['scaletypeframe4'])
        d.set('scale mode '+self.plotConf['scalemodeframe4'])
        d.set('cmap '+self.plotConf['colorframe4'])
        d.set('regions', 'image; text '+self.plotConf['labelposition']+' # color='+self.plotConf['labelcolor']+' width='+self.plotConf['labelwidth']+' font="'+self.plotConf['labelfont']+'" text={Significance}')
 

 	""""Searches through the model file for the sources (RA, Dec and name) to plot on the count map"""
 
        xml = open(self.commonConf['base']+"_model.xml")

        keywords  = ['RA', 'DEC', 'PointSource']
        keywords2 = ['value', 'name']
        ra      = []
        dec     = []
        name    = []
        columns = []

        for line in xml:
	
            for word in line.split():
	    
	 	""""Look for the source RA"""

                if keywords[0] in word:
                    for word2 in line.split():
                        if keywords2[0] in word2:
                            columns = word2.split('"')
                            ra.append(columns[1])
                            break
			    
	 	""""Look for the source Dec"""

                if keywords[1] in word:
                    for word2 in line.split():
                        if keywords2[0] in word2:
                            columns = word2.split('"')
                            dec.append(columns[1])
                            break
			    
	 	""""Look for the source name"""

                if keywords[2] in word:
                    for word2 in line.split('"'):
                        columns.append(word2)
                    i = -1
                    for lookname in columns:
                        i += 1
                        if keywords2[1] in lookname:
                            name.append(columns[i+1])
                            break
			    
        d.set('frame 1')    
	
	""""Loads on the count map the regions associated with each source in the xml file """
	        
        i = -1
        for value in ra:
            i += 1
            d.set('regions', 'fk5; point '+ ra[i]+ ' '+ dec[i]+ ' # point='+self.plotConf['sourceregiontype']+' color='+self.plotConf['sourceregioncolor']+' width='+self.plotConf['sourceregionwidth']+' font="'+self.plotConf['sourceregionfont']+'" text={'+ name[i]+ '}')

	""""Makes ds9 look pretty"""       
	     
        d.set('zoom to fit')
        d.set('match frames wcs')

        if(self.plotConf['grid'] == 'yes'):       
	    d.set('grid yes')
            d.set('grid axes type exterior')
            d.set('grid numlab vertical yes')
            d.set('grid skyformat degrees')
            d.set('grid axes color '+self.plotConf['gridcolor'])
            d.set('grid tick color '+self.plotConf['gridcolor'])
            d.set('grid grid color '+self.plotConf['gridcolor'])
            d.set('grid numlab color '+self.plotConf['gridcolor'])
            d.set('grid numlab fontsize '+self.plotConf['gridfont'])



    def runAll(self, run=True):

        """Generates the model, residual and significance maps and
	plot them together with the input count map.  This is the 
	function called when this module is run from the command 
        line.  You need to have two files to start with:
        <basename>_CMAP.fits (count map) and <basename>_model.xml 
	(source model). <basename> is a user defined prefix (usually 
	the source name but not necessarily).
        Returns an exception if any of the files are not found."""

        self.logger.info("***Checking for files***")
	
        try:
            checkForFiles(self.logger,[self.commonConf['base']+"_CMAP.fits",self.commonConf['base']+"_model.xml"])
        except(FileNotFound):
            self.logger.critical("One or more needed files do not exist")
            return

        self.logger.info("***Creating model map***")
        self.createModelMap(run)
        self.logger.info("***Creating residual map***")
        self.createResidMap(run)
        self.logger.info("***Creating significance map***")
        self.createSigMap(run)
        self.logger.info("***Plotting maps with DS9***")
        self.plotMaps(run)



def cli():
    """Command-line interface.  Call this without any options for usage notes."""
    import getopt
    class BadUsage: pass
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'cMRSP')
        
        for opt, val in opts:
            if opt == '-c':
                print "Creating example configuration file called example.cfg"
                qP = quickPlot("example")
                qP.writeConfig()
                return
            elif opt == '-M':
                print "Creating model map"
                if(not args): 
                    raise BadUsage
                else:
                    qP = quickPlot(sys.argv[2],True)
                    runModel(qP.logger,qP.commonConf['base'],qP.commonConf['irfs'])
                    return
            elif opt == '-R':
                print "Creating residual map"
                if(not args):
                    raise BadUsage
                else:
                    qP = quickPlot(sys.argv[2],True)
                    qP.createResMap()
                    return
            elif opt == '-S':
                print "Creating significance map"
                if(not args):
                    raise BadUsage
                else:
                    qP = quickPlot(sys.argv[2],True)
                    qP.createSigMap()
                    return
            elif opt == '-P':
                print "Creating ds9 plots"
                if(not args):
                    raise BadUsage
                else:
                    qP = quickPlot(sys.argv[2], True)
                    qP.plotMaps()
                    return

        if not args: raise BadUsage
        for arg in args:
            qP = quickPlot(arg, True)
            qP.runAll(True)

    except (getopt.error, BadUsage):
        cmd = os.path.basename(sys.argv[0])
        print """
                           - quickPlot - 

Generates the model, residual, and significance maps starting from the
count map and plots them using ds9.  You can use the command line
functions listed below or run this module from within python. For full
documentation on this module execute 'pydoc quickPlot'.

%s <basename> ...  Plots the results of the Fermi LAT analysis for the
    source <basename>. <basename> is the prefix used for this
    analysis.  You must already have a configuration file if using the
    command line interface.

%s -c ... Generate a default config file called example.cfg.  Edit
    this file and rename it <basename>.cfg for use in the quickPlot
    module.

%s -M <basename> ... Generate a model map from the 2FGL.  You need to
    already have <basename>_filtered_gti.fits in your working
    directory.  You can get this file by running the functions
    runSelect and runGTI (within quickAnalysis) on your data.  You
    also need to have the Galactic and isotropic diffuse models in
    your working directory as well as the 2FGL model file.

%s -R <basename> ... Generate a residual map based on the count and
    model map.

%s -S <basename> ... Generate a significance map based on the count
    and model map.

%s -P <basename> ... Plot the count, model, residual, and significance
	maps using ds9.

""" %(cmd,cmd,cmd,cmd,cmd,cmd)

if __name__ == '__main__': cli()
