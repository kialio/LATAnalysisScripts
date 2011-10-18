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
__version__ = '0.1'


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
    #sourcename and model are already in the config file in the likelihood 
    #section.  I suggest using those and not adding two more which will just
    #add to confusion.

    def __init__(self,
                 base='MySource',
                 configFile = False,
                 plotConfig = {"model" : "MySource_model.xml",
                               "sourcename" : "Source Name"},
                 commonConfig = {"base" : 'MySource',
                                 "binned" : False,
                                 "eventclass" : 2, 
                                 "irfs" : "P7SOURCE_V6",
                                 "verbosity" : 0}):

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
                
        self.commonConf = commonConfig
        self.plotConf = plotConfig

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

        """Generates a residual map using the ftool command "farith"."""

        #JSP - Need file check here.
        
        try:
            checkForCommand(self.logger, ["farith"])
        except(CommandNotFound):
            return

        cmd = "farith "+self.commonConf['base']+"_CMAP.fits "+self.commonConf['base']+"_modelMap.fits "+self.commonConf['base']+"_residMap.fits SUB clobber=yes"
            
        if(run):
            os.system(cmd)
            self.logger.info(cmd)
        else:
            print cmd
             

    def createSigMap(self,run=True):

        """Generates a significance map."""

        #JSP - Need file check here.

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

        #Think about what options to put in the config file.  The user should
        #be able to turn labeling on and off and choose which plots they want.

        try:
            checkForFiles(self.logger,
	    		 [self.commonConf['base']+"_CMap.fits",
			  self.commonConf['base']+"_modelMap.fits",
			  self.commonConf['base']+"_residMap.fits",
			  self.commonConf['base']+"_sigMap.fits"])
        except(FileNotFound):
            self.logger.critical("One or more needed files do not exist")
            return

        d = ds9()
        d.set('file '+self.commonConf['base']+"_CMAP.fits")
        d.set('scale log')
        d.set('scale mode minmax')
        d.set('cmap aips0')
        d.set('regions', 'image; text 40 15 # color=black width=3 font="helvetica 10 bold" text={Count map}')
        d.set('tile')
        d.set('frame new')
        d.set('file '+self.commonConf['base']+"_modelMap.fits")
        d.set('cmap a')
        d.set('regions', 'image; text 40 15 # color=black width=3 font="helvetica 10 bold" text={Model map}')
        d.set('frame new')
        d.set('file '+self.commonConf['base']+"_residMap.fits")
        d.set('scale sqrt')
        d.set('scale mode zscale')
        d.set('cmap aips0')
        d.set('regions', 'image; text 40 15 # color=black width=3 font="helvetica 10 bold" text={Residual}')
        d.set('frame new')
        d.set('file '+self.commonConf['base']+"_sigMap.fits")
        d.set('scale mode minmax')
        d.set('regions', 'image; text 40 15 # color=black width=3 font="helvetica 10 bold" text={Significance}')
 

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
            d.set('regions', 'fk5; point '+ ra[i]+ ' '+ dec[i]+ ' # point=cross 15 color=black width=3 font="helvetica 14 bold" text={'+ name[i]+ '}')
	   

	""""Makes ds9 look pretty"""       
	     
        d.set('zoom to fit')
        d.set('match frames wcs')
        d.set('grid yes')
        d.set('grid axes type exterior')
        d.set('grid numlab vertical yes')
        d.set('grid skyformat degrees')
        d.set('grid axes color black')
        d.set('grid tick color black')
        d.set('grid grid color black')
        d.set('grid numlab color black')
        d.set('grid numlab fontsize 14')


# Needs a function that wraps all the needed functions.  Like the 'runAll' function in 
# quickAnalysis.  This would create all the maps and then run the ds9 command.

def cli():
    """Command-line interface.  Call this without any options for usage notes."""
    import getopt
    class BadUsage: pass
    
    #make sure you get the command line interface to work.  You should have access to all of the
    #individual functions from the command line as well as an overall runall command.

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'cm')
        
        for opt, val in opts:
            if opt == '-c':
                print "Creating example configuration file called example.cfg"
                qP = quickPlot("example")
                qP.writeConfig()
                return
            elif opt == '-m':
                print "Creating model file from 2FGL called example_model.xml"
                if(not args): 
                    raise BadUsage
                else:
                    qA = quickAnalysis(sys.argv[2])
                    qA.generateXMLmodel()
                    return
            elif opt == '-C':
                print "Creating count map"
                if(not args):
                    raise BadUsage
                else:
                    qA = quickAnalysis(sys.argv[2])
                    qA.runCMAP()
                    return
            elif opt == '-M':
                print "Creating model map"
                if(not args):
                    raise BadUsage
                else:
                    qA = quickAnalysis(sys.argv[2])
                    qA.runModel()
                    return

        if not args: raise BadUsage
        for arg in args:
            qP = quickPlot(arg, True)
            qP.runAll(True)

    except (getopt.error, BadUsage):
        cmd = os.path.basename(sys.argv[0])
        print """quickAnalysis - Perform event selections and exposure
        calculations for Fermi LAT data.  You can use the command line
        functions listed below or run this module from within
        python. For full documentation on this module execute 'pydoc
        quickAnalysis'.

python %s <basename> ...  Perform an analysis on <basename>.
    <basename> is the prefix used for this analysis.  You must already
    have a configuration file if using the command line interface.

python %s -c ... Generate a default config file called example.cfg.
    Edit this file and rename it <basename>.cfg for use in the
    quickPlot module.

python %s -m <basename> ... Generate a model file from the 2FGL.  You
    need to already have <basename>_filtered_gti.fits in your working
    directory.  You can get this file by running the functions
    runSelect and runGTI on your data.  You also need to have the
    galactic and isotropic diffuse models in your working directory as
    well as the 2FGL model file.

python %s -C <basename> ... Generate a count map based on the model
    file in your config file.  You need to have several files already
    computed.  It's best to do the runAll script before trying this.

python %s -M <basename> ... Generate a model map based on the model
    file in your config file.  You need to have several files already
    computed.  It's best to do the runAll script before trying this.

""" %(cmd,cmd,cmd)

if __name__ == '__main__': cli()
