#!/usr/bin/env python

"""Perform a likelihood analysis for Fermi LAT data.  You should have
completed all of the event selection and exposure calculations in
quickAnalysis before using this moduel.  A usual likelihood analysis
will consists of running the following functions (assuming you have a
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

You can create a configuration file by executing writeConfig().  Or by
running

> python quickLike -c

from the command line.  

This module will catch any failures from the optimizers and will
report them to the user.  There are a few functions that are useful to
use in this case:



"""

__author__ = 'Jeremy S. Perkins (FSSC)'
__version__ = '0.1'

import pyLikelihood
import re
from quickUtils import *
from UnbinnedAnalysis import *
from BinnedAnalysis import *
from UpperLimits import UpperLimits

class quickLike:

    def __init__(self,
                 base = 'MySource',
                 configFile = False,
                 likelihoodConfig = {"model" : "MySource_model.xml",
                                     "sourceName" : "Source Name",
                                     "drmtol" : 0.1,
                                     "mintol" : 1e-4},
                 commonConfig = {"base" : 'MySource',
                                 "eventClass" : 2,
                                 "binned" : False,
                                 "irfs" : "P7SOURCE_V6",
                                 "verbosity" : 0}):
                                  
        commonConfig['base'] = base

        self.logger = initLogger(base, 'quickLike')

        if(configFile):
            try:
                commonConfig,analysisConfig,likelihoodConfig = readConfig(self.logger,base)
            except(FileNotFound):
                self.logger.critical("One or more needed files do not exist")
                return

        self.commonConf = commonConfig
        self.likelihoodConf = likelihoodConfig
        
        self.ret = re.compile('\n')

        self.Print()
        
    def writeConfig(self):

        """Writes all of the initialization variables to the config
        file called <basename>.cfg"""

        writeConfig(quickLogger=self.logger,
                    commonDictionary=self.commonConf,
                    likelihoodDictionary=self.likelihoodConf)

    def Print(self):
        logString = "Created quickLike object: "
        for variable, value in self.commonConf.iteritems():
            logString += variable+"="+str(value)+","
        for variable, value in self.likelihoodConf.iteritems():
            logString += variable+"="+str(value)+","
        self.logger.info(logString)

    def makeObs(self):
        if(self.commonConf['binned']):
            try:
                checkForFiles(self.logger,[self.commonConf['base']+'_srcMaps.fits',
                                          self.commonConf['base']+'_ltcube.fits',
                                          self.commonConf['base']+'_BinnedExpMap.fits'])
                self.obs = BinnedObs(srcMaps=self.commonConf['base']+'_srcMaps.fits',
                                     expCube=self.commonConf['base']+'_ltcube.fits',
                                     binnedExpMap=self.commonConf['base']+'_BinnedExpMap.fits',
                                     irfs=self.commonConf['irfs'])
            except(FileNotFound):
                self.logger.critical("One or more needed files do not exist")
                return
        else:
            try:
                checkForFiles(self.logger,[self.commonConf['base']+'_filtered_gti.fits',
                                           self.commonConf['base']+'_SC.fits',
                                           self.commonConf['base']+'_expMap.fits',
                                           self.commonConf['base']+'_ltCube.fits'])
                self.obs = UnbinnedObs(self.commonConf['base']+'_filtered_gti.fits',
                                       self.commonConf['base']+'_SC.fits',
                                       expMap=self.commonConf['base']+'_expMap.fits',
                                       expCube=self.commonConf['base']+'_ltCube.fits',
                                       irfs=self.commonConf['irfs'])
            except(FileNotFound):
                self.logger.critical("One or more needed files do not exist")
                return
        self.logger.info(self.ret.subn(', ',str(self.obs))[0])

    def initDRM(self):
        
        try:
            checkForFiles(self.logger,[self.likelihoodConf['model']])
            if(self.commonConf['binned']):
                self.DRM = BinnedAnalysis(self.obs,self.likelihoodConf['model'],optimizer="DRMNGB")
            else:
                self.DRM = UnbinnedAnalysis(self.obs,self.likelihoodConf['model'],optimizer="DRMNGB")
                self.DRM.tol = float(self.likelihoodConf['drmtol'])
                self.logger.info(self.ret.subn(', ',str(self.DRM))[0])
        except(FileNotFound):
            self.logger.critical("One or more needed files do not exist")
            return

    def initAltFit(self,opt="MINUIT"):
        try:
            checkForFiles(self.logger,[self.likelihoodConf['model']])
            if(self.commonConf['binned']):
                self.ALTFIT = BinnedAnalysis(self.obs,self.likelihoodConf['model'],optimizer=opt)
            else:
                self.ALTFIT = UnbinnedAnalysis(self.obs,self.likelihoodConf['model'],optimizer=opt)
            self.ALTFIT.tol = float(self.likelihoodConf['drmtol'])
            self.ALTFITobj = pyLike.Minuit(self.ALTFIT.logLike)
            self.logger.info(self.ret.subn(', ',str(self.ALTFIT))[0])
        except(FileNotFound):
            self.logger.critical("One or more needed files do not exist")
            return

    def initMIN(self):
        try:
            checkForFiles(self.logger,[self.commonConf['base']+'_likeDRM.xml'])
            if(self.commonConf['binned']):
                self.MIN = BinnedAnalysis(self.obs,self.commonConf['base']+'_likeDRM.xml',optimizer='NewMinuit')
            else:
                self.MIN = UnbinnedAnalysis(self.obs,self.commonConf['base']+'_likeDRM.xml',optimizer='NewMinuit')
            self.MIN.tol = float(self.likelihoodConf['mintol'])
            self.MINobj = pyLike.NewMinuit(self.MIN.logLike)
            self.logger.info(self.ret.subn(', ',str(self.MIN))[0])
        except(FileNotFound):
            self.logger.critical("One or more needed files do not exist")
            return

    def fitDRM(self):

        altfit=False
        try:
            self.DRM.fit(verbosity=int(self.commonConf['verbosity']))
        except:
            self.logger.error("Initial DRM Fit Failed")
            try:
                self.logger.info("Trying tighter tolerance (DRMtol*0.1)")
                self.DRM.tol = float(self.likelihoodConf['drmtol']) * 0.1
                self.DRM.fit(verbosity= self.commonConf['verbosity'])
            except:
                self.logger.error("Second DRM Fit Failed")
                try:
                    self.logger.info("Trying looser tolerance (drmtol*10.)")
                    self.DRM.tol = float(self.likelihoodConf['drmtol']) * 10.
                    self.DRM.fit(verbosity= self.commonConf['verbosity'])
                except:
                    self.logger.error("Third DRM Fit Failed")
                    try:
                        self.logger.info("Trying alternate fit algorithm (MINUIT)")
                        self.initAltFit()
                        self.ALTFIT.fit(verbosity=self.commonConf['verbosity'],covar=True,optObject=self.ALTFITobj)
                        print self.ALTFITobj.getQuality()
                        altfit = True
                    except:
                        self.logger.error("Alternative fit algorithm failed, bailing")
                        self.logger.error(self.decodeRetCode('Minuit',self.MINobj.getRetCode()))
                        return

        if(altfit):
            self.logger.info("ALTFIT Fit Finished.  Total TS: "+str(self.ALTFIT.logLike.value()))
            self.ALTFIT.logLike.writeXml(self.commonConf['base']+'_likeDRM.xml')
            self.logger.info("Saved ALTFIT as "+self.commonConf['base']+"_likeDRM.xml")
        else:
            self.DRM.logLike.writeXml(self.commonConf['base']+'_likeDRM.xml')
            self.logger.info("DRM Fit Finished.  Total TS: "+str(self.DRM.logLike.value()))
            self.logger.info("Saved DRM as "+self.commonConf['base']+"_likeDRM.xml")

    def fitMIN(self):
        self.MIN.fit(covar=True, optObject=self.MINobj,verbosity=int(self.commonConf['verbosity']))
        self.MIN.logLike.writeXml(self.commonConf['base']+'_likeMinuit.xml')
        self.logger.info("NEWMINUIT Fit Finished.  Total TS: "+str(self.MIN.logLike.value()))
        self.logger.info("NEWMINUIT Fit Status: "+str(self.MINobj.getRetCode()))
        self.logger.info("NEWMINUIT fit Distance: "+str(self.MINobj.getDistance()))
        if(self.MINobj.getRetCode() > 0):
            self.logger.error("NEWMINUIT DID NOT CONVERGE!!!")
            self.logger.error("The fit failed the following tests: "+self.decodeRetCode('NewMinuit',self.MINobj.getRetCode()))

    def printSource(self,source,Emin=100,Emax=300000):
        print "TS: ",self.MIN.Ts(source)
        print "Npred: ",self.MIN.NpredValue(source)
        print "Flux: ",self.MIN.flux(source,emin=Emin,emax=Emax)
        print "Flux Error: ",self.MIN.fluxError(source,emin=Emin,emax=Emax)
        print "Index: ",self.MIN.model[source].funcs['Spectrum'].getParam('Index').value()
        print "Index Error: ",self.MIN.model[source].funcs['Spectrum'].getParam('Index').error()
        self.logger.info(source+" TS: "+str(self.MIN.Ts(source))+\
                             " Npred: "+str(self.MIN.NpredValue(source))+\
                             " Flux: "+str(self.MIN.flux(source,emin=Emin,emax=Emax))+\
                             " Flux Error: "+str(self.MIN.fluxError(source,emin=Emin,emax=Emax))+\
                             " Index: "+str(self.MIN.model[source].funcs['Spectrum'].getParam('Index').value())+\
                             " Index Error: "+str(self.MIN.model[source].funcs['Spectrum'].getParam('Index').error()))
        
    def calcUpper(self,source,Emin=100,Emax=300000):
        self.ul = UpperLimits(self.MIN)
        self.ul[source].compute(emin=Emin,emax=Emax)
        print self.ul[source].results
        self.logger.info(source+" UL: "+str(self.ul[source].results[0]))

    def removeWeak(self,mySource = '',tslimit=0,DistLimit=0,RemoveFree=False,RemoveFixed=False):
        if(mySource == ''):
            mySource = self.likelihoodConf['sourceName']
        for name in self.MIN.sourceNames():
            remove = False
            distance = 0
            sourceTS = self.MIN.Ts(name)
            if(self.MIN.model[name].src.getType() == 'Point'):
                distance = self.MIN._separation(self.MIN.model[mySource].src,self.MIN.model[name].src)
            if(self.MIN.freePars(name).size() > 0):
                indexFree = "Free"
                if( (sourceTS < tslimit) and (distance > DistLimit) and RemoveFree ):
                    remove = True
            else:
                indexFree = "Fixed"
                if( (sourceTS < tslimit) and (distance > DistLimit) and RemoveFixed ):
                    remove = True
            if( remove ):
                self.logger.info("Removing "+name+", TS: "+str(sourceTS)+", Frozen?: "+str(indexFree)+", Distance: "+str(distance))
                self.MIN.deleteSource(name)
            else:
                self.logger.info("Retaining "+name+", TS: "+str(sourceTS)+", Frozen?: "+str(indexFree)+", Distance: "+str(distance))
                    

    def paramsAtLimit(self, limit = 0.01):
        for src in self.MIN.sourceNames():
            for name in self.MIN.model[src].funcs['Spectrum'].paramNames:
                bounds = self.MIN.model[src].funcs['Spectrum'].getParam(name).getBounds()
                value  = self.MIN.model[src].funcs['Spectrum'].getParam(name).value()
                if( abs(bounds[0] - value) < limit ):
                    self.logger.error("The "+name+" ("+str(value)+") of "+src+" is close to its lower limit ("+str(bounds[0])+")")
                if( abs(bounds[1] - value) < limit):
                    self.logger.error("The "+name+" ("+str(value)+") of "+src+" is close to its lower limit ("+str(bounds[1])+")")



    def decodeRetCode(self, optimizer, retCode):

        if(optimizer == 'NewMinuit'):

            retCode -= 100
            
            failure = ""
            
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
            
            if(retCode == 0):
                failure = "Error matrix not calculated at all"
            if(retCode == 1):
                failure = "Diagonal approximation only, not accurate"
            if(retCode == 2):
                failure = "Full matrix, but forced positive-definite (i.e. not accurate)"
            if(retCode == 3):
                failure = "Full accurate covariance matrix (After MIGRAD, this is the indication of normal convergence.)"

            return failure

    def delObs(self):
        del self.obs
