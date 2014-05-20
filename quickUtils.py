#!/usr/bin/env python

"""Various common functions and utilities for the quick modules

"""

__author__ = 'Jeremy S. Perkins (FSSC)'
__version__ = '0.2.0'

import os
import logging
import math
import ConfigParser
import numpy as np

from gt_apps import *

class FileNotFound: pass
class CommandNotFound: pass

def log_array(npts, xmin, xmax):

    '''This function creates an array with npts-1 logarithmically spaced
    bins borrowed from macro Jim sent to do profile likelihood.'''

    xstep = np.log(xmax/xmin)/(npts - 1)
    return xmin*np.exp(np.arange(npts, dtype=np.float)*xstep)

def checkForFiles(quickLogger, fileList):
    
    """Checks for the existence of needed files in the list."""
    
    for filename in fileList:
        if(not os.path.exists(filename)):
            quickLogger.critical(filename+" doesn't exist.")
            raise FileNotFound

def checkForCommand(quickLogger, commandList):

    """Checks for the existence of a certain command."""

    for command in commandList:

        cmd = "which -s " + command + " > " + os.devnull + " 2>&1"
        retcode = os.system(cmd)
        
        if(retcode):
            quickLogger.critical("unix command "+command+" not found.")
            raise CommandNotFound
        
        
def writeConfig(quickLogger, commonDictionary, analysisDictionary = {}, likelihoodDictionary = {}, plotDictionary = {}, curveDictionary = {}):

    """Writes all of the needed information to the config file called
    <basename>.cfg"""
    
    basename = commonDictionary['base']

    config = ConfigParser.RawConfigParser()
    config.read(basename+'.cfg')
    if(not config.has_section('common')):
        config.add_section('common')

    for variable, value in commonDictionary.iteritems():
        config.set('common', variable, value)
    quickLogger.info("wrote common config to "+basename+".cfg.")

    if(analysisDictionary):
        if(config.has_section('quickAnalysis')):
            quickLogger.info("quickAnalysis config exists, overwriting...")        
        else:
            config.add_section('quickAnalysis')            
        for variable, value in analysisDictionary.iteritems():
            config.set('quickAnalysis', variable, value)
        quickLogger.info("wrote quickAnalysis config to "+basename+".cfg.")

    if(likelihoodDictionary):
        if(config.has_section('quickLike')):
            quickLogger.info("quickLike config exists, overwriting...")        
        else:
            config.add_section('quickLike')            
        for variable, value in likelihoodDictionary.iteritems():
            config.set('quickLike', variable, value)
        quickLogger.info("wrote quickLikeconfig to "+basename+".cfg.")

    if(plotDictionary):
        if(config.has_section('quickPlot')):
            quickLogger.info("quickPlot config exists, overwriting...")        
        else:
            config.add_section('quickPlot')            
        for variable, value in plotDictionary.iteritems():
            config.set('quickPlot', variable, value)
        quickLogger.info("wrote quickPlot config to "+basename+".cfg.")

    if(curveDictionary):
        if(config.has_section('quickCurve')):
            quickLogger.info("quickCurve config exists, overwriting...")        
        else:
            config.add_section('quickCurve')            
        for variable, value in curveDictionary.iteritems():
            config.set('quickCurve', variable, value)
        quickLogger.info("wrote quickCurve config to "+basename+".cfg.")

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
    plotDictionary = {}
    curveDictionary = {}

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

        if(config.has_section('quickPlot')):
            quickLogger.info('Reading quickPlot variables...')
            plotDictionary = dict(config.items('quickPlot'))

        if(config.has_section('quickCurve')):
            quickLogger.info('Reading quickCurve variables...')
            curveDictionary = dict(config.items('quickCurve'))
            if( curveDictionary['sliding'] in ['True', 'true', '1', 'yes']):
                curveDictionary['sliding'] = True
            else:
                curveDictionary['sliding'] = False

        return commonDictionary,analysisDictionary,likelihoodDictionary,plotDictionary,curveDictionary

    except(FileNotFound):
        raise FileNotFound
        return

def checkConfig(quickLogger, referenceDictionary,testDictionary):

    """Checks a dictionary against a refernece to make sure that all
    of the parameters are there.  If all is good, it'll returen the 
    checked dictionary.  If not, it'll return the reference dictionary
    and raise an exception."""

    try:
        for key in referenceDictionary:
            item = testDictionary[key]
        return testDictionary
    except KeyError as inst:
        quickLogger.critical("Cannot find "+inst.args[0]+" in the config file.")
        raise KeyError
        return referenceDictionary

        
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
                     galactic_file="gll_iem_v05.fit",
                     galactic_name="gll_iem_v05",
                     isotropic_file="iso_source_v05.txt",
                     isotropic_name="iso_source_v05",
                     catalog_file="gll_psc_v08.fit"):
    
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
            mymodel.makeModel(galactic_file, galactic_name, isotropic_file, isotropic_name)
            quickLogger.info("NOTE: if there are extended sources in your ROI, make sure the "\
                             +"correspoinding diffuse template is in the working directory.")
        except(FileNotFound):
            raise FileNotFound
      

def runCommand(AppCommand,quickLogger,run=True,printCmd=False):

    """Runs a giving command if run is True.  If run is False,
    prints out what the function would run."""

    if(run):
        AppCommand.run(print_command=printCmd)
        quickLogger.info(AppCommand.command())
    else:
        print AppCommand.command()
            

def runModel(quickLogger,
	     base,
             modelFile="",
	     irfs="P7REP_SOURCE_V15",
             run=True):
	
    """Generates a model map.  You need to have already run
    the general quickAnlysis tool and then fit your model with 
    quickLike so that all of the needed files exsist."""
    
    if(modelFile):
        model = modelFile
    else:
        model = base+"_likeMinuit.xml"


    try:
        checkForFiles(quickLogger,
                      [base+"_srcMaps.fits",
                       model,
                       base+"_ltcube.fits",
                       base+"_BinnedExpMap.fits"])
    except(FileNotFound):
        quickLogger.critical("One or more needed files do not exist.")
        return

    model_map['srcmaps'] = base+"_srcMaps.fits"
    model_map['srcmdl']  = model
    model_map['outfile'] = base+"_modelMap.fits"
    model_map['expcube'] = base+"_ltcube.fits"
    model_map['irfs']    = irfs
    model_map['bexpmap'] = base+"_BinnedExpMap.fits"
  
    runCommand(model_map,quickLogger,run)

def runCMAP(quickLogger,
            base,
            rad,
            binsize,
            ra,
            dec,
            nxpix,
            nypix,
            run=True):
            	
        """Generates a counts map.  The dimensions of which are the
        largest square subtended by the ROI.  Note that if the ROI is
        exceptionally small or the bin size exceptionally large, the
        square might not be the largest posible since the npix
        calculation floors the calculated value."""

        if nxpix < 0 or nypix < 0:
        	nxpix = NumberOfPixels(float(rad),float(binsize))
        	nypix = NumberOfPixels(float(rad),float(binsize))

        evtbin['evfile'] = base+'_filtered_gti.fits'
        evtbin['outfile'] = base+'_CMAP.fits'
        evtbin['scfile'] = base+"_SC.fits"
        evtbin['algorithm'] = 'CMAP'
        evtbin['nxpix'] = nxpix
        evtbin['nypix'] = nypix
        evtbin['binsz'] = binsize
        evtbin['coordsys'] = 'CEL'
        evtbin['xref'] = ra
        evtbin['yref'] = dec
        evtbin['axisrot'] = 0
        evtbin['proj'] = 'AIT'
    
        runCommand(evtbin,quickLogger,run)

class quickMath:

    '''This class is used in the quickCurve script and are various
    statistical functions developed by Stephen Fegan
    <sfegan@llr.in2p3.fr> based on Numerical recipes in C.'''

    _def_itmax  = 100
    _def_reltol = 1e-15

    @staticmethod
    def _gamma_ser(x, a, gammalna = None):

        '''See Numerical recipes in C - eq 6.2.5.'''

        if x<=0.0:
            if x==0.0: return 0.0
            else: raise ValueError("Argument x is negative: x=%f"%x)
        if a<=0.0:
            raise ValueError("Argument a is zero or negative: a=%f"%a)
        if gammalna == None:
            gammalna = math.lgamma(a)
        actr = a
        dsum = 1.0/a
        sum = dsum
        for i in range(1,quickMath._def_itmax):
            actr += 1.0
            dsum *= x/actr
            # print i,dsum,sum,sum+dsum
            sum  += dsum
            if math.fabs(dsum) < quickMath._def_reltol*math.fabs(sum):
                return sum*math.exp(a*math.log(x) - x - math.lgamma(a))
        raise RuntimeError("Maximum number of iterations exceeded")

    @staticmethod
    def _gamma_cfrac(x, a, gammalna = None):
        # See Numerical recipes in C - eq 6.2.7 and section 5.2
        if x<=0.0:
            if x==0.0: return 0.0
            else: raise ValueError("Argument x is negative: x=%f"%x)
        if a<=0.0:
            raise ValueError("Argument a is zero or negative: a=%f"%a)
        if gammalna == None:
            gammalna = math.lgamma(a)
        tiny = 1e-30
        A = 0.0
        B = x+1.0-a
        F = B
        if F == 0.0: F = tiny
        C = F
        D = 0.0
        for i in range(1,quickMath._def_itmax):
            A += a + 1.0 - 2.0*i
            B += 2.0
            D = B + A*D
            if D == 0.0: D = tiny
            C = B + A/C
            if C == 0.0: C = tiny
            D = 1.0/D
            delta = C*D
            F *= delta
            # print i,F
            if math.fabs(delta-1.0) < quickMath._def_reltol:
                return math.exp(a*math.log(x) - x - gammalna)/F
        raise RuntimeError("Maximum number of iterations exceeded")

    @staticmethod
    def gammainc(x, a):
        if x<=0.0:
            if x==0.0: return 0.0
            else: raise ValueError("Argument x is negative: x=%f"%x)
        if a<=0.0:
            raise ValueError("Argument a is zero or negative: a=%f"%a)
        if x<a+1.0:
            return quickMath._gamma_ser(x, a)
        else:
            return 1.0-quickMath._gamma_cfrac(x, a)

    @staticmethod
    def gammaincc(x, a):
        if x<=0.0:
            if x==0.0: return 0.0
            else: raise ValueError("Argument x is negative: x=%f"%x)
        if a<=0.0:
            raise ValueError("Argument a is zero or negative: a=%f"%a)
        if x<a+1.0:
            return 1.0-quickMath._gamma_ser(x, a)
        else:
            return quickMath._gamma_cfrac(x, a)

    @staticmethod
    def gammainv(p, a):
        if p<=0.0:
            if p==0.0: return 0.0
            else: raise ValueError("Argument p is negative: p=%f"%p)
        elif p>=1:
            if p==1: return float('infinity')
            else: raise ValueError("Argument p is greater than unity: p=%f"%p)
        if a<=0:
            raise ValueError("Argument a is zero or negative: a=%f"%a)

        gammalna = math.lgamma(a)
        xtest = a+1.0
        ftest = quickMath._gamma_ser(xtest, a, gammalna)
        ffind = p
        if(ffind<ftest):
            f = lambda x: quickMath._gamma_ser(x, a, gammalna)
            dfdx = lambda x,f: math.exp((a-1.0)*math.log(x) - x - gammalna)
            ffind = p
            if ffind<0.75:
                xtest = math.exp((math.log(ffind) + gammalna + math.log(a))/a)
                if ffind<0.1*quickMath._def_reltol:
                    return xtest
        else:
            f = lambda x: math.log(quickMath._gamma_cfrac(x, a, gammalna))
            dfdx = lambda x,f: -math.exp((a-1.0)*math.log(x) - x - gammalna - f)
            ffind = math.log(1-p)
        ftest = f(xtest)
        dfdxtest = dfdx(xtest, ftest)
        for i in range(1,quickMath._def_itmax):
            xnext = xtest - (ftest-ffind)/dfdxtest
            if xnext <= 0.0:
                xnext = 0.5*(xtest+0.0)
                #print i, xtest, ftest, dfdxtest, xnext, ffind, math.fabs(xnext-xtest), 10*quickMath._def_reltol*(xnext+xtest)
            if math.fabs(xnext-xtest) < 10*quickMath._def_reltol*(xnext+xtest):
                return xnext
            xtest = xnext
            ftest = f(xtest)
            dfdxtest = dfdx(xtest, ftest)
        raise RuntimeError("Maximum number of iterations exceeded")

    @staticmethod
    def gammainvc(p, a):
        if p<=0.0:
            if p==0.0: return 0.0
            else: raise ValueError("Argument p is negative: p=%f"%p)
        elif p>=1:
            if p==1: return float('infinity')
            else: raise ValueError("Argument p is greater than unity: p=%f"%p)
        if a<=0:
            raise ValueError("Argument a is zero or negative: a=%f"%a)

        gammalna = math.lgamma(a)
        xtest = a+1.0
        ftest = quickMath._gamma_cfrac(xtest, a, gammalna)
        ffind = p
        if(ffind<ftest):
            f = lambda x: math.log(quickMath._gamma_cfrac(x, a, gammalna))
            dfdx = lambda x,f: -math.exp((a-1.0)*math.log(x) - x - gammalna - f)
            ffind = math.log(p)
        else:
            f = lambda x: quickMath._gamma_ser(x, a, gammalna)
            dfdx = lambda x,f: math.exp((a-1.0)*math.log(x) - x - gammalna)
            ffind = 1-p
            if ffind<0.75:
                xtest = math.exp((math.log(ffind) + gammalna + math.log(a))/a)
                if ffind<0.1*quickMath._def_reltol:
                    return xtest
        ftest = f(xtest)
        dfdxtest = dfdx(xtest, ftest)
        for i in range(1,quickMath._def_itmax):
            xnext = xtest - (ftest-ffind)/dfdxtest
            if xnext <= 0.0:
                xnext = 0.5*(xtest+0.0)
    #        print i, xtest, ftest, dfdxtest, xnext, ffind, math.fabs(xnext-xtest), 10*quickMath._def_reltol*(xnext+xtest)
            if math.fabs(xnext-xtest) < 0.5*quickMath._def_reltol*(xnext+xtest):
                return xnext
            xtest = xnext
            ftest = f(xtest)
            dfdxtest = dfdx(xtest, ftest)
        raise RuntimeError("Maximum number of iterations exceeded")

    @staticmethod
    def chi2cdf(x, a):
        return quickMath.gammainc(0.5*x, 0.5*a)

    @staticmethod
    def chi2cdfc(x, a):
        return quickMath.gammaincc(0.5*x, 0.5*a)

    @staticmethod
    def chi2inv(p, a):
        return 2.0*quickMath.gammainv(p, 0.5*a)

    @staticmethod
    def chi2invc(p, a):
        return 2.0*quickMath.gammainvc(p, 0.5*a)
