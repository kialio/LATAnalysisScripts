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

__author__ = 'Jeremy S. Perkins (FSSC)'
__version__ = '0.1.2'

import pyLikelihood
import re
from quickUtils import *
from UnbinnedAnalysis import *
from BinnedAnalysis import *
from UpperLimits import UpperLimits

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

    def __init__(self,
                 base = 'MySource',
                 configFile = False,
                 likelihoodConfig = {"model" : "MySource_model.xml",
                                     "sourcename" : "Source Name",
                                     "drmtol" : 0.1,
                                     "mintol" : 1e-4},
                 commonConfig = {"base" : 'MySource',
                                 "eventclass" : 2,
                                 "binned" : False,
                                 "irfs" : "P7SOURCE_V6",
                                 "verbosity" : 0}):
                                  
        commonConfig['base'] = base

        self.logger = initLogger(base, 'quickLike')

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
                likelihoodConfig = checkConfig(self.logger,likelihoodConfig,likelihoodConfigRead)
            except(KeyError):
                return

        self.commonConf = commonConfig
        self.likelihoodConf = likelihoodConfig
        
        self.ret = re.compile('\n')
        self.fitbit = False
        self.Print()
        
    def writeConfig(self):

        """Writes all of the initialization variables to the config
        file called <basename>.cfg"""

        writeConfig(quickLogger=self.logger,
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
                                           self.commonConf['base']+'_ltcube.fits'])
                self.obs = UnbinnedObs(self.commonConf['base']+'_filtered_gti.fits',
                                       self.commonConf['base']+'_SC.fits',
                                       expMap=self.commonConf['base']+'_expMap.fits',
                                       expCube=self.commonConf['base']+'_ltcube.fits',
                                       irfs=self.commonConf['irfs'])
            except(FileNotFound):
                self.logger.critical("One or more needed files do not exist")
                return
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

    def initMIN(self, useBadFit=False):

        """Initiallizes a New Minuit optimizer to use as a backup to
        the DRM optimizer.  This is usually run after you have
        initially run fitDRM and created a <basename>_likeDRM.xml
        model file which is used a seed for the New Minuit optimizer.
        You can skip the DRM process if you like but you need to have
        the proper model file (<basename>_likeDRM.xml) present in the
        working directory. You need to run makeObs before you run this
        function.  If it hasn't been run, this function will exit.  If
        you want to use the non convergant fit from fitDRM, set
        useBadFit to True."""

        try:
            self.obs
        except AttributeError:
            self.logger.critical("Obs object does not exist.  Create it first with the makeObs function.")
            return

        if(useBadFit):
            model = self.commonConf['base']+'_badDRMFit.xml'
        else:
            model = self.commonConf['base']+'_likeDRM.xml'

        try:
            checkForFiles(self.logger,[model])
            if(self.commonConf['binned']):
                self.MIN = BinnedAnalysis(self.obs,model,optimizer='NewMinuit')
            else:
                self.MIN = UnbinnedAnalysis(self.obs,model,optimizer='NewMinuit')
            self.MIN.tol = float(self.likelihoodConf['mintol'])
            self.MINobj = pyLike.NewMinuit(self.MIN.logLike)
            self.logger.info(self.ret.subn(', ',str(self.MIN))[0])
        except(FileNotFound):
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
                    try:
                        self.logger.info("Trying alternate fit algorithm (MINUIT)")
                        self.initAltFit()
                        self.ALTFIT.fit(verbosity=int(self.commonConf['verbosity']),covar=True,optObject=self.ALTFITobj)
                        print self.ALTFITobj.getQuality()
                        altfit = True
                    except:
                        self.logger.error("Alternative fit algorithm failed, bailing")
                        self.logger.error(self.decodeRetCode('Minuit',self.ALTFITobj.getRetCode()))
                        self.ALTFIT.logLike.writeXml(self.commonConf['base']+'_badDRMFit.xml')
                        self.logger.info("Saved ALTFIT as "+self.commonConf['base']+"_badDRMFit.xml")
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

        self.MIN.fit(covar=True, optObject=self.MINobj,verbosity=int(self.commonConf['verbosity']))
        self.MIN.logLike.writeXml(self.commonConf['base']+'_likeMinuit.xml')
        self.logger.info("NEWMINUIT Fit Finished.  Total TS: "+str(self.MIN.logLike.value()))
        self.logger.info("NEWMINUIT Fit Status: "+str(self.MINobj.getRetCode()))
        self.logger.info("NEWMINUIT fit Distance: "+str(self.MINobj.getDistance()))
        self.fitbit = True
        if(self.MINobj.getRetCode() > 0):
            self.logger.error("NEWMINUIT DID NOT CONVERGE!!!")
            self.logger.error("The fit failed the following tests: "+self.decodeRetCode('NewMinuit',self.MINobj.getRetCode()))
            self.MIN.fit(covar=True, optObject=self.MINobj,verbosity=int(self.commonConf['verbosity']))
            self.MIN.logLike.writeXml(self.commonConf['base']+'_badMINFit.xml')
            

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
        print "TS: ",TS
        logString += " TS: " + str(TS)
        NPred = self.MIN.NpredValue(source)
        print "Npred: ",NPred
        logString += " NPred: " + str(NPred)
        flux = self.MIN.flux(source,emin=Emin,emax=Emax)
        print "Flux: ",flux
        logString += "Flux: "+str(flux)
        if(self.fitbit):
            fluxErr = self.MIN.fluxError(source,emin=Emin,emax=Emax)
            print "Flux Error: ",fluxErr
            logString += "Flux Error: "+str(fluxErr)
        for paramName in self.MIN.model[source].funcs['Spectrum'].paramNames:
            paramValue = self.MIN.model[source].funcs['Spectrum'].getParam(paramName).value()
            print paramName,": ",paramValue
            logString += paramName + ": " + str(paramValue) + " "

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

    def removeWeak(self,mySource = '',tslimit=0,distlimit=0,RemoveFree=False,RemoveFixed=False):

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

        <obj>.removeWeak(<my_source>,tslimit=3,distlimit=10,RemoveFree=True,RemoveFixed=True)"""

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
        for name in self.MIN.sourceNames():
            remove = False
            distance = 0
            sourceTS = self.MIN.Ts(name)
            if(self.MIN.model[name].src.getType() == 'Point'):
                distance = self.MIN._separation(self.MIN.model[mySource].src,self.MIN.model[name].src)
            if(self.MIN.freePars(name).size() > 0):
                indexFree = "Free"
                if( (sourceTS < tslimit) and (distance > distlimit) and RemoveFree ):
                    remove = True
            else:
                indexFree = "Fixed"
                if( (sourceTS < tslimit) and (distance > distlimit) and RemoveFixed ):
                    remove = True
            if( remove ):
                self.logger.info("Removing "+name+", TS: "+str(sourceTS)+", Frozen?: "+str(indexFree)+", Distance: "+str(distance))
                self.MIN.deleteSource(name)
            else:
                self.logger.info("Retaining "+name+", TS: "+str(sourceTS)+", Frozen?: "+str(indexFree)+", Distance: "+str(distance))
                    

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

    def decodeRetCode(self, optimizer, retCode):

        """Decodes the return codes from the Minuit and New Minuit fit
        functions.  Used in the fitting functions in this module.
        You'll probably never use this function."""

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

# Command-line interface                                                                                                                      
def cli():
    """Command-line interface.  Call this without any options for usage notes."""
    import getopt

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'i', ['initialize'])

        for opt, val in opts:
            if opt in ('-i','--initialize'):
                print "Creating example configuration file called example.cfg"
                qL = quickLike("example")
                qL.writeConfig()
                return

        if not opts: raise getopt.GetoptError("Must specify an option, printing help.")

    except getopt.error as e:
        print "Command Line Error: " + e.msg
        cmd = os.path.basename(sys.argv[0])
        print """
                        - quickLike - 

Perform a liklihood analysis on Fermi LAT data.  You can use the
command line functins listed below or run this module from withing
python. For full documentation on this module execute 'pydoc
quickLike'.
                                              
%s (-i|--initialize) ... Generate a default config file called
    example.cfg.  Edit this file and rename it <basename>.cfg for use
    in the quickLike module.

""" %(cmd)
                                                                                                                                            
if __name__ == '__main__': cli()
