#!/usr/bin/env python

"""Perform a likelihood analysis for Fermi LAT data.  

You should have completed all of the event selection and exposure
calculations in quickAnalysis before using this moduel. 

First, generate a default config file:

> quickLike (-i|--initialize)

Then edit the config file to match your specific analysis by filling
out the various options.  Rename the config file from example.cfg to
<basename>.cfg where <basename> is a user defined prefix; usually the
name of your source of interest but not necisarily so.

If you haven't created a model file (if you have, make sure it's
called <basename>_model.xml before you use this module) you can create
one from the 2FGL by running

> quickAnalysis (-x|--xml=)<basename>

To run this, you need to have all of the relevant diffuse model files
and the 2FGL catalog file in your working directory as well as
<basename>_filtered_gti.fits in your working directory.  See the
documentation for quickAnalysis for more details.

Once you have a config file (<basename>.cfg) and model file
(<basename>_model.xml) you can then run python and execute the other
functions within the python interpreter.

This module logs all of the steps to a file called
<basename>_quickAnalysis.log as well as to the screen.

"""

import LATAnalysisScripts as LAS

__author__ = LAS.__author__
__version__ = LAS.__version__

import re
import os
import sys
import quickUtils as qU
import pyLikelihood as pyLike
import UnbinnedAnalysis as UbAn
import BinnedAnalysis as BAn
import numpy as np

from UpperLimits import UpperLimits
from LikelihoodState import LikelihoodState
from LATAnalysisScripts import defaultConfig as dC
from LATAnalysisScripts.Logger import Logger
from LATAnalysisScripts.Logger import logLevel as ll


class quickLike:

    """ This is the base class.  A usual likelihood analysis will
    consists of running the following functions (assuming you have a
    configuration file):
    
    * qL = quickLike('MySource', True)
    * qL.makeObs()
    * qL.initDRM()
    * qL.fitDRM()
    * qL.initMIN()
    * qL.fitMIN()

    This will set up all of the objects needed for the analysis and do an
    initial fit with one of the DRM optimizers.  It'll save these results
    and use them for the second fit with one of the Minuit optimizers.
    
    If you do not have a configuration file, you'll need to input all
    of the options for this module when you create the quickLike
    object (see the various options below).  You can create a
    configuration file by executing writeConfig().

    * qL.writeConfig()

    This module will catch any failures from the optimizers and will
    report them to the user.  There are a few functions that are useful to
    use in this case:"""

    __version__ = LAS.__version__

    def __init__(self,
                 base = 'MySource',
                 configFile = False,
                 likelihoodConfig = dC.defaultLikelihoodConfig,
                 commonConfig = dC.defaultCommonConfig):
                                
        commonConfig['base'] = base

        self.logger = Logger(base, self.__class__.__name__,ll(commonConfig['verbosity'])).get()
        self.logger.info("This is quickLike version {}.".format(__version__))

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
            #Reset the verboisty from the config file if needed
            if ll(commonConfig['verbosity']) != self.logger.handlers[1].level:
                self.logger.info("Resetting the log level on the console to {}.".format(commonConfig['verbosity']))
                self.logger.handlers[1].setLevel(ll(commonConfig['verbosity']))
            try:
                likelihoodConfig = qU.checkConfig(self.logger,likelihoodConfig,likelihoodConfigRead)
            except(KeyError):
                return

        self.commonConf = commonConfig
        self.likelihoodConf = likelihoodConfig
        
        self.ret = re.compile('\n')
        self.fitbit = False
        self._init_psf()
        self.Print()

    def __saveModel__(self):

        """Private function that stores the current model into a
        dictionary.  Used when restoring a value to its default."""

        self.oldmdl = {}
        for src in self.MIN.sourceNames():
            parNames = pyLike.StringVector()
            self.MIN[src].src.spectrum().getParamNames(parNames)
            spec = {}
            for par in parNames:
                spec[par] = self.MIN.model[src].funcs['Spectrum'].getParam(par).value()
                self.oldmdl[src] = spec

    def writeConfig(self):

        """Writes all of the initialization variables to the config
        file called <basename>.cfg"""

        qU.writeConfig(quickLogger=self.logger,
                       commonDictionary=self.commonConf,
                       likelihoodDictionary=self.likelihoodConf)

    def Print(self):

        """Prints out information about the various objects to the
        terminal and to the log file."""

        logString = "Created quickLike object: "
        for variable, value in self.commonConf.iteritems():
            logString += variable+"="+str(value)+","
        for variable, value in self.likelihoodConf.iteritems():
            logString += variable+"="+str(value)+","
        self.logger.info(logString)

    def makeObs(self):
        
        """Creates either a binned or unbinned observation object for
        use in the likelihood analysis.  This function checks for all
        of the needed files first.  If you do not have a needed file,
        see the quickAnalysis module for creation.  This function
        should be run before any of the init or fit functions."""

        if(self.commonConf['binned']):
            try:
                qU.checkForFiles(self.logger,[self.commonConf['base']+'_srcMaps.fits',
                                          self.commonConf['base']+'_ltcube.fits',
                                          self.commonConf['base']+'_BinnedExpMap.fits'])
                self.obs = BAn.BinnedObs(srcMaps=self.commonConf['base']+'_srcMaps.fits',
                                         expCube=self.commonConf['base']+'_ltcube.fits',
                                         binnedExpMap=self.commonConf['base']+'_BinnedExpMap.fits',
                                         irfs=self.commonConf['irfs'])
            except(qU.FileNotFound):
                self.logger.critical("One or more needed files do not exist")
                sys.exit()
        else:
            try:
                qU.checkForFiles(self.logger,[self.commonConf['base']+'_filtered_gti.fits',
                                              self.commonConf['base']+'_SC.fits',
                                              self.commonConf['base']+'_expMap.fits',
                                              self.commonConf['base']+'_ltcube.fits'])
                self.obs = UbAn.UnbinnedObs(self.commonConf['base']+'_filtered_gti.fits',
                                       self.commonConf['base']+'_SC.fits',
                                       expMap=self.commonConf['base']+'_expMap.fits',
                                       expCube=self.commonConf['base']+'_ltcube.fits',
                                       irfs=self.commonConf['irfs'])
            except(qU.FileNotFound):
                self.logger.critical("One or more needed files do not exist")
                sys.exit()
        self.logger.info(self.ret.subn(', ',str(self.obs))[0])

    def initDRM(self):

        """Initializes the DRM optimizer (either binned or unbinned).
        This is usually the second function that you run when using
        this module.  You need to run makeObs before you run this
        function.  If it hasn't been run, this function will exit."""
        
        try:
            self.obs
        except AttributeError:
            self.logger.critical("Obs object does not exist.  Create it first with the makeObs function")
            return

        try:
            qU.checkForFiles(self.logger,[self.likelihoodConf['model']])
            if(self.commonConf['binned']):
                self.DRM = BAn.BinnedAnalysis(self.obs,self.likelihoodConf['model'],optimizer="DRMNGB")
            else:
                self.DRM = UbAn.UnbinnedAnalysis(self.obs,self.likelihoodConf['model'],optimizer="DRMNGB")
                self.DRM.tol = float(self.likelihoodConf['drmtol'])
                self.logger.info(self.ret.subn(', ',str(self.DRM))[0])
        except(qU.FileNotFound):
            self.logger.critical("One or more needed files do not exist")
            return

    def initAltFit(self,opt="MINUIT"):

        """Initiallizes a minuit optimizer to use as a backup to the
        DRM optimizer.  This function is used internally in the fitDRM
        function so you probably will never use it.  You need to run
        makeObs before you run this function.  If it hasn't been run,
        this function will exit."""

        try:
            self.obs
        except AttributeError:
            self.logger.critical("Obs object does not exist.  Create it first with the makeObs function")
            return

        try:
            qU.checkForFiles(self.logger,[self.likelihoodConf['model']])
            if(self.commonConf['binned']):
                self.ALTFIT = BAn.BinnedAnalysis(self.obs,self.likelihoodConf['model'],optimizer=opt)
            else:
                self.ALTFIT = UbAn.UnbinnedAnalysis(self.obs,self.likelihoodConf['model'],optimizer=opt)
            self.ALTFIT.tol = float(self.likelihoodConf['drmtol'])
            self.ALTFITobj = pyLike.Minuit(self.ALTFIT.logLike)
            self.logger.info(self.ret.subn(', ',str(self.ALTFIT))[0])
        except(qU.FileNotFound):
            self.logger.critical("One or more needed files do not exist")
            return

    def initMIN(self, useBadFit=False, modelFile="", useEdisp=False, new=True):

        """Initiallizes a New Minuit optimizer to use as a backup to
        the DRM optimizer.  This is usually run after you have
        initially run fitDRM and created a <basename>_likeDRM.xml
        model file which is used a seed for the New Minuit optimizer.
        You can skip the DRM process if you like but you need to have
        the proper model file (<basename>_likeDRM.xml) present in the
        working directory. You need to run makeObs before you run this
        function.  If it hasn't been run, this function will exit.  If
        you want to use the non convergant fit from fitDRM, set
        useBadFit to True.  You can also pass a custom model file name
        via the modelFile parameter."""

        try:
            self.obs
        except AttributeError:
            self.logger.critical("Obs object does not exist.  Create it first with the makeObs function.")
            return

        if(useBadFit):
            model = self.commonConf['base']+'_badDRMFit.xml'
        else:
            model = self.commonConf['base']+'_likeDRM.xml'

        if(modelFile):
            model = modelFile

        if(new):
            opt = 'NewMinuit'
        else:
            opt = 'Minuit'

        try:
            qU.checkForFiles(self.logger,[model])
            if(self.commonConf['binned']):
                self.MIN = BAn.BinnedAnalysis(self.obs,model,optimizer=opt)
            else:
                self.MIN = UbAn.UnbinnedAnalysis(self.obs,model,optimizer=opt)
            self.MIN.tol = float(self.likelihoodConf['mintol'])            
            if(new):
                self.MINobj = pyLike.NewMinuit(self.MIN.logLike)
            else:
                self.MINobj = pyLike.Minuit(self.MIN.logLike)    
            self.pristine = LikelihoodState(self.MIN)
            self.__saveModel__()
            self.logger.info(self.ret.subn(', ',str(self.MIN))[0])
            if(useEdisp):
                self.MIN.logLike.set_edisp_flag(useEdisp)
        except(qU.FileNotFound):
            self.logger.critical("One or more needed files do not exist")
            return

    def fitDRM(self):

        """Performs a DRM inital fit on your data using the
        <basename>_model.xml model file.  It tries an intial fit and
        if that fails, tries a tighter tolerance.  If that fails, it
        tries a looser tolerance.  If that fails, it tries to do this
        initial fit with the MINUIT optimizer.  If that fails, this
        function bails.  If the fit converges, it saves the results to
        <basename>_likeDRM.xml which will be used in the NewMinuit
        fit.  If no fit is found, it will save the results to
        <basename>_badDRMFit.xml.  You can use this in the NewMinuit fit
        if you use the useBadFit option in initMIN.  You need to have run
        initDRM before you run this function."""

        try:
            self.DRM
        except AttributeError:
            self.logger.critical("DRM object does not exist.  Create it first with the initDRM function.")
            return

        altfit=False
        try:
            self.DRM.fit(verbosity=int(self.commonConf['verbosity']))
        except:
            self.logger.error("Initial DRM Fit Failed")
            try:
                self.logger.info("Trying tighter tolerance (DRMtol*0.1)")
                self.DRM.tol = float(self.likelihoodConf['drmtol']) * 0.1
                self.DRM.fit(verbosity= int(self.commonConf['verbosity']))
            except:
                self.logger.error("Second DRM Fit Failed")
                try:
                    self.logger.info("Trying looser tolerance (drmtol*10.)")
                    self.DRM.tol = float(self.likelihoodConf['drmtol']) * 10.
                    self.DRM.fit(verbosity= int(self.commonConf['verbosity']))
                except:
                    self.logger.error("Third DRM Fit Failed")
                    self.logger.info("Trying alternate fit algorithm (MINUIT)")
                    self.initAltFit()
                    self.ALTFIT.fit(verbosity=int(self.commonConf['verbosity']),covar=True,optObject=self.ALTFITobj)
                    self.logger.info("Fit quality: {}".format(self.ALTFITobj.getQuality()))
                    if(self.ALTFITobj.getQuality() < 3):
                        self.logger.error("Alternative fit algorithm failed, bailing")
                        self.logger.error(self.decodeRetCode('Minuit',self.ALTFITobj.getRetCode()))
                        self.ALTFIT.logLike.writeXml(self.commonConf['base']+'_badDRMFit.xml')
                        self.logger.info("Saved ALTFIT as "+self.commonConf['base']+"_badDRMFit.xml")
                        return
                    else:
                        self.logger.info(self.decodeRetCode('Minuit',self.ALTFITobj.getRetCode()))
                        altfit = True

        if(altfit):
            self.logger.info("ALTFIT Fit Finished.  -log(likelihood): "+str(self.ALTFIT.logLike.value()))
            self.ALTFIT.logLike.writeXml(self.commonConf['base']+'_likeDRM.xml')
            self.logger.info("Saved ALTFIT as "+self.commonConf['base']+"_likeDRM.xml")
        else:
            self.DRM.logLike.writeXml(self.commonConf['base']+'_likeDRM.xml')
            self.logger.info("DRM Fit Finished.  -log(likelihood): "+str(self.DRM.logLike.value()))
            self.logger.info("Saved DRM as "+self.commonConf['base']+"_likeDRM.xml")

    def fitMIN(self):

        """Does a New Minuit fit on your data based on the model
        output by the fitDRM function.  You need to have run initMIN
        before running this function.  Saves the results to
        <basename>_likeMIN.xml if there is convergence.  If
        convergence is not found, saves the results to
        <basename>_badMINFit.xml."""

        try:
            self.MIN
        except AttributeError:
            self.logger.critical("MIN object does not exist.  Create it first with the initMIN function.")
            return

        if(isinstance(self.MINobj,pyLike.Minuit)):
            opt = 'Minuit'
            self.MINobj = pyLike.Minuit(self.MIN.logLike)
        else:
            self.MINobj = pyLike.NewMinuit(self.MIN.logLike)
            opt = 'NewMinuit'

        self.MIN.fit(covar=True, optObject=self.MINobj,verbosity=int(self.commonConf['verbosity']))
        self.logger.info("{} Fit Finished.  -log(likelihood): {}".format(opt,self.MIN.logLike.value()))
        self.logger.info("{} Fit Status: {}".format(opt,self.MINobj.getRetCode()))
        self.logger.info("{} fit Distance: {}".format(opt,self.MINobj.getDistance()))
        self.fitbit = True
        if(self.MINobj.getRetCode() > 0):
            self.logger.error("{} DID NOT CONVERGE!!!".format(opt))
            self.logger.error("The fit failed the following tests: "+self.decodeRetCode(opt,self.MINobj.getRetCode()))
            self.MIN.logLike.writeXml(self.commonConf['base']+'_badMINFit.xml')
        else:
            self.MIN.logLike.writeXml(self.commonConf['base']+'_likeMinuit.xml')
            
    def pokeSource(self,source,paramName='Prefactor'):

        """This function pokes a paramter of a source to a value that is 10%
        of what it was.  This is a useful function to use when you are trying
        to get convergence on a fit.  Many times, a fit won't converge because
        the initial model is too close to the final answer (ie. the minimizer
        does not have enough flexibility to accurately calculate a correlation
        matrix).  In this case, run this funciton on one of the stronger
        sources in your model and redo the fit.  This function is also useful
        to determine how robust your fit is. Note that it defaults to using
        poking the 'Prefactor' parameter which might not exist for your
        specific source.  In that case, choose a different paramter."""

        previousValue = self.MIN.model[source].funcs['Spectrum'].getParam(paramName).value()
        self.MIN.model[source].funcs['Spectrum'].func.setParam(paramName, 0.1*previousValue)
        self.logger.info("Resetting the {} of {} from {:,.2f} to {:,.2f}".format(paramName,source,previousValue,0.1*previousValue))

    def printSource(self,source,Emin=100,Emax=300000):

        """Prints various details for a source in your model."""
        try:
            self.MIN
        except AttributeError:
            self.logger.critical("MIN object does not exist. "+\
                                     "Create it first with the initMIN function and then fit it with the fitMIN function.")
            return

        if(not self.fitbit):
            self.logger.warn("Fit isn't current, these values might not be correct. Fun fitMIN first.")

            
        logString = source
        TS = self.MIN.Ts(source)
        print "TS: {:,.2f}".format(TS)
        logString += " TS: {:,.2f} ".format(TS)
        NPred = self.MIN.NpredValue(source)
        print "Npred: {:,.2f}".format(NPred)
        logString += " NPred: {:,.2f} ".format(NPred)
        flux = self.MIN.flux(source,emin=Emin,emax=Emax)
        outString = "Flux: {:,.2e}".format(flux)
        logString += " Flux: {:,.2e} ".format(flux)
        if(self.fitbit):
            fluxErr = self.MIN.fluxError(source,emin=Emin,emax=Emax)
            outString += " +- {:,.2e}".format(fluxErr)
            logString += " Flux Error: {:,.2e} ".format(fluxErr)
        print outString
        for paramName in self.MIN.model[source].funcs['Spectrum'].paramNames:
            paramValue = self.MIN.model[source].funcs['Spectrum'].getParam(paramName).value()
            paramError = self.MIN.model[source].funcs['Spectrum'].getParam(paramName).error()
            paramScale = self.MIN.model[source].funcs['Spectrum'].getParam(paramName).parameter.getScale()
            print paramName,": {:,.2f} +- {:,.2f} x {:,.2e}".format(paramValue,paramError,paramScale)
            logString += paramName + ": {:,.2f} +- {:,.2f} x {:,.2e} ".format(paramValue,paramError,paramScale)

        self.logger.info(logString)

    def customERange(self,Emin,Emax):

        """Sets a smaller energy range for the fitting of both the DRM
        and MIN optimization steps."""

        try:
            self.DRM
        except AttributeError:
            self.logger.warn("DRM object doesn't exist.  Energy range not modified.")
        else:
            self.DRM.setEnergyRange(Emin,Emax)
            self.logger.info("Set energy range for DRM to "+str(self.DRM.emin)+","+str(self.DRM.emax))

        try:
            self.MIN
        except AttributeError:
            self.logger.warn("MIN object doesn't exist.  Energy range not modified.")
        else:
            self.MIN.setEnergyRange(Emin,Emax)
            self.logger.info("Set energy range for MIN to "+str(self.MIN.emin)+","+str(self.MIN.emax))
                
    def calcUpper(self,source,Emin=100,Emax=300000):

        """Calculates an upper limit for a source in your model."""

        self.ul = UpperLimits(self.MIN)
        self.ul[source].compute(emin=Emin,emax=Emax)
        print self.ul[source].results
        self.logger.info(source+" UL: "+str(self.ul[source].results[0]))

    def removeWeak(self,mySource = '',tslimit=0,npredlimit=0,distlimit=0,RemoveFree=False,RemoveFixed=False,recalc=True):

        """This function has two main uses: it will print out details
        on all of the sources in your model and it will remove sources
        according to different requirements.  If you just want to
        print out details, execute it this way:

        <obj>.removeWeak(<my_source>)

        Where <obj> is the quickLike object you're using here and
        <my_source> is the name of your source of interest.  You can
        then remove some of these sources from the model if you like.
        For example, if you want to remove all of the fixed sources
        with TS values less than 1, execute it this way:

        <obj>.removeWeak(<my_source>,tslimit=1,RemoveFixed=True)

        You can mix and match any of the options.  You could remove
        all sources (fixed and free) that are below a TS value of 3
        and are 10 degrees from your source of interest by executing:

        <obj>.removeWeak(<my_source>,tslimit=3,distlimit=10,RemoveFree=True,RemoveFixed=True)

        This funcitons precalculates all of the relavant source
        details and saves them to <obj>.sourceDetails.  You can choose
        to not recalculate these values and use the saved ones by
        passing 'recalc=False' to the function."""

        try:
            self.MIN
        except AttributeError:
            self.logger.critical("MIN object does not exist. "+\
                                     "Create it first with the initMIN function and then fit it with the fitMIN function.")
            return

        if(not self.fitbit):
            self.logger.warn("Fit isn't current, these values might not be correct.  Run fitMIN first.")

        if(mySource == ''):
            mySource = self.likelihoodConf['sourcename']

        try:
            self.MIN.model[mySource].src.getName()
        except AttributeError:
            self.logger.critical(mySource+" is not in the model.  Either pass a valid "+
                                 " sourcename to this function or modify your model "+
                                 "and/or config file to indicate which source you are considering.")
            return

        if(recalc):
            self.sourceDetails = {}

            self.logger.info("Calculating TS/NPred values for all sources.  This could take a bit.")
            for name in self.MIN.sourceNames():
                sourceTS = self.MIN.Ts(name)
                sourceNPred = self.MIN.NpredValue(name)
                if(self.MIN.model[name].src.getType() == 'Point'):
                    distance = self.MIN._separation(self.MIN.model[mySource].src,self.MIN.model[name].src)
                if(np.shape(self.MIN.freePars(name))[0] > 0):
                    free = True
                else:
                    free = False
                
                self.sourceDetails[name] = {'TS' : sourceTS,
                                            'NPred' : sourceNPred, 
                                            'dist' : distance,
                                            'free' : free}


        for source,details in self.sourceDetails.iteritems():
            remove = False
            if(details['free']):
                indexFree = "Free"
                if( (details['TS'] < tslimit) and (details['dist'] > distlimit) and RemoveFree ):
                    remove = True
                if( (details['NPred'] < npredlimit) and (details['dist'] > distlimit) and RemoveFree ):
                    remove = True    
            else:
                indexFree = "Fixed"
                if( (details['TS'] < tslimit) and (details['dist'] > distlimit) and RemoveFixed ):
                    remove = True
                if( (details['NPred'] < npredlimit) and (details['dist'] > distlimit) and RemoveFixed ):
                    remove = True   
            logString =  "{}, TS: {:.2f}, NPred: {:.2f}, Frozen?: {}, Distance: {}".format(source,details['TS'],
                                                                                    details['NPred'],
                                                                                    indexFree,
                                                                                    details['dist'])
            if( remove ):
                self.logger.info("Removing " + logString)
                self.MIN.deleteSource(source)
            else:
                self.logger.info("Retaining "+ logString)

    def unLoadSource(self, name):

        """This function removes a source from the model and stores it so that
        you can use it later if you would like.  This is useful if you are
        working on an upper limit and need to get a fit to work before you can
        calculate the upper limit."""
                    
        self.saved_src = self.MIN.deleteSource(name)
        self.logger.info("Removed "+name+" from the model and saved it.")
        
    def reLoadSource(self):
        
        """This function puts the source removed by the unLoadSource function
        back into the model."""

        try:
            self.MIN.addSource(self.saved_src)        
            self.logger.info("Reloaded saved source.")
        except AttributeError:
            self.logger.critical("Saved Source does not exist. "+\
                                     "Make sure that you've run the unLoadSource function.")
            return

    def paramsAtLimit(self, limit = 0.1):

        """This function will print out any sources whoes parameters
        are close to their limits.  You could use this to find sources
        that are having issues being fit.  This function is useful
        when you're having trouble getting convergence from the New
        Minuit fit routine. The limit is in percentage difference of a
        bound.  If one of the bounds is zero it uses the value of the
        parameter to check for closeness (absolute instead of percent
        differenct).  The default is 0.1 (1%) difference for a measure
        of closeness."""

        try:
            self.MIN
        except AttributeError:
            self.logger.critical("MIN object does not exist. "+\
                                     "Create it first with the initMIN function and then fit it with the fitMIN function.")
            return

        if(not self.fitbit):
            self.logger.warn("Fit isn't current, these values might not be correct. Run fitMIN first.")

        for src in self.MIN.sourceNames():
            for name in self.MIN.model[src].funcs['Spectrum'].paramNames:
                bounds = self.MIN.model[src].funcs['Spectrum'].getParam(name).getBounds()
                value  = self.MIN.model[src].funcs['Spectrum'].getParam(name).value()

                try:
                    distToLower = abs((value - bounds[0])/bounds[0])
                except ZeroDivisionError:
                    distToLower = abs(value)

                try:
                    distToUpper = abs((value - bounds[1])/bounds[1])
                except ZeroDivisionError:
                    distToUpper = abs(value)

                if( distToLower < limit ):
                    self.logger.error("The "+name+" ("+str(value)+") of "+src+" is close ("\
                                          +str(distToLower)+") to its lower limit ("+str(bounds[0])+")")
                if( distToUpper < limit):
                    self.logger.error("The "+name+" ("+str(value)+") of "+src+" is close ("\
                                          +str(distToUpper)+") to its upper limit ("+str(bounds[1])+")")

    def calcBowtie(self,srcName,minE,maxE,numBins):
        
        '''This is derived from T. Johnson's likeSED code which was in turn
        derived from D. Sanchez's pyUnfoldPlot code which was probably
        based on some code developed by J. Chiang.  '''

        '''make some energy bounds for the fit, same max and min as for the
        bands before but with more bins.'''

        modEs=qU.log_array(numBins,minE,maxE)
        centEs=[0.5*(e1+e2) for e1,e2 in zip(modEs[0:-1],modEs[1:])]

        '''Get the model.'''
        mysrc=pyLike.PointSource_cast(self.MIN[srcName].src)
        spec=[float(1000.*mysrc.spectrum()(pyLike.dArg(x))) for x in centEs]

        if(self.MIN.covariance is None):
            print "Whoa, you didn't compute the covariance yet..."
            bt=[0]
        else:
            bt=[]
            covArray=np.array(self.MIN.covariance)
            srcCovArray=[]
            par_index_map={}
            indx=0
            for src in self.MIN.sourceNames():
                parNames=pyLike.StringVector()
                self.MIN[src].src.spectrum().getFreeParamNames(parNames)
                for par in parNames:
                    par_index_map['::'.join((src,par))]=indx
                    indx +=1
            srcPars=pyLike.StringVector()
            self.MIN[srcName].src.spectrum().getFreeParamNames(srcPars)
            pars=['::'.join((srcName,x)) for x in srcPars]
            for xpar in pars:
                ix=par_index_map[xpar]
                srcCovArray.append([covArray[ix][par_index_map[ypar]] for ypar in pars])
            cov=np.array(srcCovArray)
            ''' The whole point here is to get the srcCovArray.'''
            for x in centEs:
                arg=pyLike.dArg(x)
                partials=np.array([mysrc.spectrum().derivByParam(arg,y) for y in srcPars])
                val=np.sqrt(np.dot(partials,np.dot(cov,partials)))
                '''These should come out same as the model so convert to ph/cm^2/s/GeV as well.'''
                bt+=[float(1000.*val)]
        return centEs,bt,spec
        
    def plot2D(self):

        """Plots a two dimensinal counts plot with the counts spectrum of each 
            source along with the total counts spectrum and the measured counts.  
            It also plots a residual plot (data-model/model) for the full model.  
            It saves the figure as <basename_2DModel.png>."""

        from LATAnalysisScripts.quickPlot2 import Plot2DModel
        Plot2DModel(self.MIN,filename="{}_2DModel.png".format(self.commonConf['base']))
        self.logger.info('Saved 2D Model as {}_2DModel.png.'.format(self.commonConf['base']))


    def plotSigMap(self):

        if hasattr(self,'sigmaMap'):
            from LATAnalysisScripts.quickPlot2 import PlotSigMap
            self.logger.info('Saved Significance Plots as {}_SigMap.png.'.format(self.commonConf['base']))
            PlotSigMap(self,filename="{}_SigMap.png".format(self.commonConf['base']))
        else:
            self.logger.error('You need to run calcSigMaps first.')


    def plotSigMaps(self):

        if hasattr(self,'sigmaMap'):
            from LATAnalysisScripts.quickPlot2 import PlotSigMaps
            self.logger.info('Saved Significance Plots as {}_SigMaps.png.'.format(self.commonConf['base']))
            PlotSigMaps(self,filename="{}_SigMaps.png".format(self.commonConf['base']))
        else:
            self.logger.error('You need to run calcSigMaps first.')

    def decodeRetCode(self, optimizer, retCode):

        """Decodes the return codes from the Minuit and New Minuit fit
        functions.  Used in the fitting functions in this module.
        You'll probably never use this function."""

        if(optimizer == 'NewMinuit'):
            
            failure = ""
        
            retCode -= 100
            
            if(retCode & 1):
                failure += " IsAboveMaxEdm"
            if(retCode & 2):
                failure += " HasCovariance"
            if(retCode & 4):
                failure += " HesseFailed"
            if(retCode & 8):
                failure += " HasMadePosDefCovar"
            if(retCode & 16):
                failure += " HasPosDefCovar"
            if(retCode & 32):
                failure += " HasAccurateCovar"
            if(retCode & 64):
                failure += " HasValidCovariance"
            if(retCode & 128):
                failure += " HasValidParameters"
            if(retCode & 256):
                failure += " IsValid"
            
            return failure

        if(optimizer == 'Minuit'):

            failure = "Unknown."

            if(retCode == 0):
                failure = "Error matrix not calculated at all"
            if(retCode == 1):
                failure = "Diagonal approximation only, not accurate"
            if(retCode == 2):
                failure = "Full matrix, but forced positive-definite (i.e. not accurate)"
            if(retCode == 3):
                failure = "Full accurate covariance matrix (After MIGRAD, this is the indication of normal convergence.)"

            return failure

    def _init_psf(self):

        '''Eventually, should caulculate the PSF for the IRFs used.'''

        from scipy import interpolate

        psf_data = np.array([[29.378069628779578, 11.974138435550062],
            [53.123991394359386, 7.957256520655659],
            [93.21305548209773, 5.218638956603347],
            [168.5559881858929, 3.333502332614621],
            [295.7537313410467, 2.073929041512339],
            [524.1755314484437, 1.3074112382033989],
            [947.8599776522383, 0.8027497013599825],
            [1679.9281120535993, 0.5127712463425323],
            [2977.4001732388747, 0.3452772300468028],
            [5383.995040496365, 0.2549696862936114],
            [9542.258178354516, 0.19331232401691525],
            [16575.885731063016, 0.15655123158441375],
            [29973.964322973727, 0.13364545640595094],
            [53660.061093243254, 0.1234818289581322],
            [95103.75714811166, 0.12512043461577865],
            [166872.09599510877, 0.13016787553232645],
            [304797.6595808452, 0.13903680232834606]])

        self.psf = interpolate.interp1d(psf_data[:,0],psf_data[:,1])


    def calcSigMaps(self, save=False):

        from quickUtils import poisson_lnl,smooth
        import pyfits

        try:
            qU.checkForFiles(self.logger,[self.commonConf['base']+'_CCUBE.fits',self.commonConf['base']+'_modelCube.fits'])
            ccube_hdu = pyfits.open("{}_CCUBE.fits".format(self.commonConf['base']))
            model_hdu = pyfits.open("{}_modelCube.fits".format(self.commonConf['base']))
            
        except(qU.FileNotFound):
            self.logger.critical("One or more needed files do not exist")
            sys.exit()

        modelSmoothed = np.zeros_like(model_hdu[0].data,dtype=np.float32)
        countSmoothed = np.zeros_like(model_hdu[0].data,dtype=np.float32)
        tsSmoothed = np.zeros_like(model_hdu[0].data,dtype=np.float32)
        residSmoothed = np.zeros_like(model_hdu[0].data,dtype=np.float32)
        sigmaSmoothed = np.zeros_like(model_hdu[0].data,dtype=np.float32)

        for Ebin,Emin in enumerate(ccube_hdu[1].data['E_MIN']):
            Emax = ccube_hdu[1].data['E_MAX'][Ebin]
            psf = self.psf(Emin/1000.)
            if psf > 0.5:
                self.logger.warn("Setting Smoothing to 0.5 Degree.")
                psf = 0.5
            deg_per_pix = 0.1
    
            modelSmoothed[Ebin],var = smooth(model_hdu[0].data[Ebin],psf,deg_per_pix, summed=True)
            countSmoothed[Ebin],var = smooth(ccube_hdu[0].data[Ebin].astype(np.float32),psf, deg_per_pix, summed=True)
            tsSmoothed[Ebin] = 2.0*(poisson_lnl(countSmoothed[Ebin],countSmoothed[Ebin]) - 
                poisson_lnl(countSmoothed[Ebin],modelSmoothed[Ebin]))
            residSmoothed[Ebin] = countSmoothed[Ebin] - modelSmoothed[Ebin]
            sigmaSmoothed[Ebin] = np.sqrt(tsSmoothed[Ebin])
            sigmaSmoothed[Ebin][residSmoothed[Ebin]<0.] *= -1.

        psf = 0.2
        deg_per_pix = 0.1
        modelFullSmoothed,var = smooth(np.sum(model_hdu[0].data,axis=0),psf,deg_per_pix, summed=True)
        countFullSmoothed,var = smooth(np.sum(ccube_hdu[0].data.astype(np.float32),axis=0),psf,deg_per_pix, summed=True)

        tsFullSmoothed = 2.0*(poisson_lnl(countFullSmoothed,countFullSmoothed) - 
            poisson_lnl(countFullSmoothed,modelFullSmoothed))
        residFullSmoothed = countFullSmoothed - modelFullSmoothed
        sigmaFullSmoothed = np.sqrt(tsFullSmoothed)
        sigmaFullSmoothed[residFullSmoothed<0.] *= -1.

        if save:
            ccube_hdu[0].data = sigmaSmoothed
            ccube_hdu.writeto("{}_SigMaps.fits".format(self.commonConf['base']), clobber = True)
            self.logger.info("Wrote {}_SigMaps.fits.".format(self.commonConf['base']))
            ccube_hdu[0].data = sigmaFullSmoothed
            ccube_hdu.writeto("{}_SigMap.fits".format(self.commonConf['base']), clobber = True)
            self.logger.info("Wrote {}_SigMap.fits.".format(self.commonConf['base']))

        self.sigmaMaps = sigmaSmoothed
        self.sigmaMap = sigmaFullSmoothed
        self.energyBins = ccube_hdu[1].data

    def restoreParam(self,source, param, setFree=True):

        """This function restores a paramter to its initial value.
        You can optionally set it to fixed (or free)."""

        oldvalue = self.MIN.model[source].funcs['Spectrum'].getParam(param).value()
        newvalue = self.oldmdl[source][param]

        self.MIN.model[source].funcs['Spectrum'].setParam(param,newvalue)
        self.logger.info('Reset the {} of {} from {} to {}.'.format(param,source,oldvalue,newvalue))
        self.MIN.model[source].funcs['Spectrum'].getParam(param).setFree(setFree)


def printCLIHelp():
    """This function prints out the help for the CLI."""
    
    cmd = os.path.basename(sys.argv[0])
    print """
                        - quickLike - 

Perform a liklihood analysis on Fermi LAT data.  You can use the
command line functins listed below or run this module from withing
python. For full documentation on this module execute 'pydoc
quickLike'.
                        
%s (-h|--help) ... This help text.
                      
%s (-i|--initialize) ... Generate a default config file called
    example.cfg.  Edit this file and rename it <basename>.cfg for use
    in the quickLike module.

""" %(cmd,cmd)

# Command-line interface    
def cli():
    """Command-line interface.  Call this without any options for usage notes."""
    import getopt

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hi', ['help',
                                                        'initialize'])

        for opt, val in opts:
            if opt in ('-h','--help'):
                printCLIHelp()
            elif opt in ('-i','--initialize'):
                print "Creating example configuration file called example.cfg"
                qL = quickLike("example")
                qL.writeConfig()
                return
            
        if not opts: raise getopt.GetoptError("Must specify an option, printing help.")
            
    except getopt.error as e:
        print "Command Line Error: " + e.msg
        printCLIHelp()
                                                                                                                                            
if __name__ == '__main__': cli()
