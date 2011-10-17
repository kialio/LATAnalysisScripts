#!/usr/bin/env python

"""Various common functions and utilities for the quick modules

"""

__author__ = 'Jeremy S. Perkins (FSSC)'
__version__ = '0.1'

import os
import logging
import math
import ConfigParser

class FileNotFound: pass

def checkForFiles(quickLogger, fileList):
    
    """Checks for the existence of needed files in the list."""
    
    for filename in fileList:
        if(not os.path.exists(filename)):
            quickLogger.critical(filename+" doesn't exist.")
            raise FileNotFound
        
def writeConfig(quickLogger, commonDictionary, analysisDictionary = {}, likelihoodDictionary = {}):

    """Writes all of the needed information to the config file called
    <basename>.cfg"""
    
    basename = commonDictionary['base']

    config = ConfigParser.RawConfigParser()
    config.read(basename+'.cfg')
    if(not config.has_section('common')):
        config.add_section('common')

    for variable, value in commonDictionary.iteritems():
        config.set('common', variable, value)

    if(analysisDictionary):
        if(config.has_section('quickAnalysis')):
            print "quickAnalysis config exists, overwriting..."
            quickLogger.info("quickAnalysis config exists, overwriting...")        
        else:
            config.add_section('quickAnalysis')            
        for variable, value in analysisDictionary.iteritems():
            config.set('quickAnalysis', variable, value)

    if(likelihoodDictionary):
        if(config.has_section('quickLike')):
            print "quickLike config exists, overwriting..."
            quickLogger.info("quickLike config exists, overwriting...")        
        else:
            config.add_section('quickLike')            
        for variable, value in likelihoodDictionary.iteritems():
            config.set('quickLike', variable, value)

    with open(basename+'.cfg', 'wb') as configfile:
        config.write(configfile)

def readConfig(quickLogger,basename):

    """Returns all of the needed information from the config file
    called <basename>.cfg.  Also checks to make sure all of the 
    config parameters are set based on the configure dictionaries
    given in the configDictionaryList."""

    commonDictionary = {}
    analysisDictionary = {}
    likelihoodDictionary = {}

    try:
        checkForFiles(quickLogger,[basename+".cfg"])
        quickLogger.info('Reading from config file ('+basename+'.cfg)')            
        config = ConfigParser.RawConfigParser()
        config.read(basename+'.cfg')
        
        if(config.has_section('common')):
            quickLogger.info('Reading common variables...')
            commonDictionary = dict(config.items('common'))
            if( commonDictionary['binned'] in ['True', 'true', '1', 'yes']):
                commonDictionary['binned'] = True
            else:
                commonDictionary['binned'] = False
            
        if(config.has_section('quickAnalysis')):
            quickLogger.info('Reading quickAnalysis variables...')
            analysisDictionary = dict(config.items('quickAnalysis'))

        if(config.has_section('quickLike')):
            quickLogger.info('Reading quickLike variables...')
            likelihoodDictionary = dict(config.items('quickLike'))

        return commonDictionary,analysisDictionary,likelihoodDictionary

    except(FileNotFound):
        raise FileNotFound
        return

def checkConfig(referenceDictionary,testDictionary):

    """Checks a dictionary against a refernece to make sure that all
    of the parameters are there."""

    try:
        for key in referenceDictionary:
            item = testDictionary[key]
        return 0
    except KeyError as inst:
        return inst.args[0]
    
        
def initLogger(base, name):

    """Sets up and returns a properly configured logging object."""

    quickLogger = logging.getLogger(name)
    quickLogger.setLevel(logging.DEBUG)
    #Prevents duuplicate log entries after reinitialization.                                                        
    if(not quickLogger.handlers):
        fh = logging.FileHandler(base+'_'+name+'.log')
        fh.setLevel(logging.DEBUG)
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        quickLogger.addHandler(fh)
        quickLogger.addHandler(ch)

    return quickLogger

def NumberOfPixels(radius, bin_size):

    """Returns the number of pixels needed to fill in the largest
    possible square subtended by the circle with the given radius.
    The size of the bins is also needed."""

    return int((radius/math.sqrt(2.0))*2.0 / bin_size)

def generateXMLmodel(quickLogger,
                     base,
                     galactic_file="gal_2yearp7v6_v0.fits",
                     isotropic_file="iso_p7v6source.txt",
                     catalog_file="gll_psc_v05.fit"):
    
    """Checks to see if <basename>_model.xml exists and creates one
    using the make2FGLXML module if it doesn't.  make2FGLXml.py must
    be in the working directory or accessable to python for this
    function to work.  The galactic and isotropic models plus the
    Fermi LAT catalog must also be in the working directory.
    Additionally, if any extended sources are in the ROI, the diffuse
    templates for those sources should be in the working directory."""


    try:
        checkForFiles(quickLogger,[base+"_model.xml"])
        quickLogger.info(base+"_model.xml exists, won't create a new one.")
    except(FileNotFound):
        quickLogger.info(base+"_model.xml doesn't exist, will create a new one.")        
        try:
            checkForFiles(quickLogger,[base+"_filtered_gti.fits",galactic_file,isotropic_file,catalog_file])
            import make2FGLxml
            mymodel = make2FGLxml.srcList(catalog_file,base+"_filtered_gti.fits",base+"_model.xml")
            mymodel.makeModel(galactic_file, 'gal_2yearp7v6_v0', isotropic_file, 'iso_p7v6source')
            quickLogger.info("NOTE: if there are extended sources in your ROI, make sure the "\
                             +"correspoinding diffuse template is in the working directory.")
        except(FileNotFound):
            raise FileNotFound
        
        
