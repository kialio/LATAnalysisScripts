#!/usr/bin/env python

"""Generate a light curve


"""

__author__ = 'Jeremy S. Perkins (FSSC)'
__version__ = '0.1.12'

import sys
import os
import quickAnalysis as qA
import quickLike as qL
import quickUtils as qU
import numpy as np
import pyLikelihood as pyLike
import UnbinnedAnalysis as UbAn
from UpperLimits import UpperLimits
from multiprocessing import Pool
from gt_apps import *
import glob

def runAnalysisStepMP(bininfo):

    bin = bininfo[0]
    tmin = bininfo[1]
    commonConf = bininfo[2]
    analysisConf = bininfo[3]
    curveConf = bininfo[4]
    tmax = tmin + float(curveConf['tstep'])

    print bin,tmin,tmax
    dir = "quickCurve_bin" + str(bin) 

    if not os.path.isdir(dir):
        os.mkdir(dir)

    analysisConfig = {"ra" : analysisConf["ra"],
                      "dec" : analysisConf["dec"],
                      "rad" : analysisConf["rad"],
                      "tmin" : tmin,
                      "tmax" : tmax,
                      "emin" : analysisConf["emin"],
                      "emax" : analysisConf["emax"],
                      "zmax" : analysisConf["zmax"],
                      "binsize" : analysisConf["binsize"],
                      "convtype" : analysisConf["convtype"]}
    commonConfig = {"base" : commonConf["base"],
                    "eventclass" : commonConf["eventclass"],
                    "binned" : False,
                    "irfs" : commonConf["irfs"],
                    "verbosity" : commonConf["verbosity"],
                    "multicore" : 0}


    qA_bin = qA.quickAnalysis(commonConf['base'],False,analysisConfig,commonConfig)
    qA_bin.analysisConf['tmin'] = tmin
    qA_bin.analysisConf['tmax'] = tmax

    qA_bin.runSelect(True, False,
                     outfile = dir + "/" + commonConf['base'] + "_filtered.fits")
    qA_bin.runGTI(True, 
                  evfile=dir + "/" + commonConf['base']+'_filtered.fits',
                  outfile=dir + "/" + commonConf['base']+'_filtered_gti.fits')
    qA_bin.runLTCube(True,
                     evfile = dir + "/" + commonConf['base']+'_filtered_gti.fits',
                     outfile = dir + "/" + commonConf['base']+'_ltcube.fits')
    qA_bin.runExpMap(True,
                     evfile = dir + "/" + commonConf['base']+'_filtered_gti.fits',
                     expcube = dir + "/" + commonConf['base']+'_ltcube.fits',
                     outfile = dir + "/" + commonConf['base']+'_expMap.fits')

class quickCurve:

    """This is the base class"""
    
    def __init__(self,
                 base = 'MySource',
                 configFile = False,
                 curveConfig = {'tstart' : 0,
                                'tstop' : 0,
                                'tstep' : 86400,
                                'tsmin' : 4,
                                'model' : 'my_model.xml'},
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
        
        self.logger = qU.initLogger(base, 'quickCurve')
        
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
                likelihoodConfig = qU.checkConfig(self.logger,likelihoodConfig,likelihoodConfigRead)
            except(KeyError):
                return
            try:
                analysisConfig = qU.checkConfig(self.logger,analysisConfig,analysisConfigRead)
            except(KeyError):
                return
            try:
                curveConfig = qU.checkConfig(self.logger,curveConfig,curveConfigRead)
            except(KeyError):
                return

        self.commonConf = commonConfig
        self.curveConf = curveConfig
        self.likelihoodConf = likelihoodConfig
        self.analysisConf = analysisConfig

    def writeConfig(self):

        """Writes all of the initialization variables to the config
        file called <basename>.cfg."""

        qU.writeConfig(quickLogger=self.logger,
		       curveDictionary=self.curveConf,
		       likelihoodDictionary=self.likelihoodConf,
                       commonDictionary=self.commonConf,
                       analysisDictionary=self.analysisConf)

    def runLikelihoodStep(self,bininfo = [0,0,0],tslimit=4.0):

        bin = bininfo[0]
        tmin = bininfo[1]
        tmax = bininfo[1] + float(self.curveConf['tstep'])

        dir = "quickCurve_bin" + str(bin) 
        basename = dir + "/" + self.commonConf['base']

        try:
            qU.checkForFiles(self.logger,[basename+'_filtered_gti.fits',
                                          self.commonConf['base']+'_SC.fits',
                                          basename+'_expMap.fits',
                                          basename+'_ltcube.fits'])
            obs = UbAn.UnbinnedObs(basename+'_filtered_gti.fits',
                                   self.commonConf['base']+'_SC.fits',
                                   expMap=basename+'_expMap.fits',
                                   expCube=basename+'_ltcube.fits',
                                   irfs=self.commonConf['irfs'])
        except(qU.FileNotFound):
            self.logger.critical("One or more needed files do not exist")
            sys.exit()
        
        try:
            qU.checkForFiles(self.logger,[self.curveConf["model"]])
            MIN = UbAn.UnbinnedAnalysis(obs, self.curveConf["model"], optimizer='NewMinuit')
            MIN.tol = float(self.likelihoodConf['mintol'])
            MINobj = pyLike.NewMinuit(MIN.logLike)
        except(qU.FileNotFound):
            self.logger.critical("One or more needed files do not exist")
            sys.exit()

        MIN.fit(covar=True, optObject=MINobj, verbosity=int(self.commonConf['verbosity']))
        self.logger.info("NEWMINUIT Fit Finished.  -log(likelihood): "+str(MIN.logLike.value()))
        self.logger.info("NEWMINUIT Fit Status: "+str(MINobj.getRetCode()))
        self.logger.info("NEWMINUIT fit Distance: "+str(MINobj.getDistance()))
        if(MINobj.getRetCode() > 0):
            self.logger.error("NEWMINUIT DID NOT CONVERGE!!!")
            #self.logger.error("The fit failed the following tests: "+self.decodeRetCode('NewMinuit',self.MINobj.getRetCode()))
            MIN.logLike.writeXml(basename+'_badMINFit.xml')
        else:
            MIN.logLike.writeXml(basename+'_likeMinuit.xml')

        sourcename = self.likelihoodConf['sourcename']

        if(MIN.Ts(sourcename) < tslimit):
            print "Will calulate upper limit for bin {}.".format(bin)
            ul = UpperLimits(MIN)
            upper = ul[sourcename].compute(emin=float(self.analysisConf["emin"]),
                                           emax=float(self.analysisConf["emax"]))
        else:
            upper = (0.0,0.0)
            
        return "{} {} {} {:,.2f} {:,.2f} {:,.2e} {:,.2e} {:,.2e} {}".format(bin,tmin,tmax,
                                                                            MIN.Ts(sourcename),
                                                                            MIN.NpredValue(sourcename),
                                                                            MIN.flux(sourcename,
                                                                                            float(self.analysisConf["emin"]),
                                                                                            float(self.analysisConf["emax"])),
                                                                            MIN.fluxError(sourcename,
                                                                                                 float(self.analysisConf["emin"]),
                                                                                                 float(self.analysisConf["emax"])),
                                                                            upper[0],
                                                                            MINobj.getRetCode())


    def runCurve(self,runAnalysis=True,runLike=True, delete = True):

        tbins = np.arange(float(self.curveConf['tstart']),
                          float(self.curveConf['tstop']),
                          float(self.curveConf['tstep']))

        bins = np.arange(0,np.size(tbins))

        binsinfo = zip(np.arange(0,np.size(tbins)),
                       tbins,
                       [self.commonConf for bin in bins],
                       [self.analysisConf for bin in bins],
                       [self.curveConf for bin in bins])

        if(runAnalysis):
            if int(self.commonConf['multicore']) > 1:
                pool = Pool(processes = int(self.commonConf['multicore']))
                pool.map(runAnalysisStepMP,binsinfo)
            else:
                for bininfo in binsinfo:
                    runAnalysisStepMP(bininfo)

        if(runLike):
            filename = self.commonConf['base'] + ".lc"
            f = open(filename, 'w')
            f.write("#bin tmin tmax TS NPred Flux FluxErr Upper FitStatus\n")
            for bininfo in zip(bins,tbins):
                output = self.runLikelihoodStep(bininfo)
                print output
                f.write(output+"\n")

        if(delete):
            templist = glob.glob("*_bin" + str(binnum) + "*")
            for t in templist:
                os.remove(t)
		

     

def printCLIHelp():
    """This function prints out the help for the CLI."""
    
    cmd = os.path.basename(sys.argv[0])
    print """
                        - quickCurve - 

Perform a liklihood analysis on Fermi LAT data.  You can use the
command line functions listed below or run this module from within
python. For full documentation on this module execute 'pydoc
quickCurve'.
                        
%s (-h|--help) ... This help text.
                      
%s (-i|--initialize) ... Generate a default config file called
    example.cfg.  Edit this file and rename it <basename>.cfg for use
    in the quickLike module.

%s (-a|--analyze) (-n |--basename=)<basename> ...  Perform an analysis
    on <basename>.  <basename> is the prefix used for this analysis.
    You must already have a configuration file if using the command
    line interface.

""" %(cmd,cmd,cmd)

# Command-line interface    
def cli():
    """Command-line interface.  Call this without any options for usage notes."""
    import getopt

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hri:n:', ['help',
                                                            'run',
                                                            'initialize'])

        #Loop through first and check for the basename
        haveBase = False
        basename = 'example'
        for opt,val in opts:
            if opt in ('-n','--basename'):
                haveBase = True
                basename = val

        for opt, val in opts:
            if opt in ('-h','--help'):
                printCLIHelp()
                return
            elif opt in ('-r', '--run'):
                if not haveBase: raise getopt.GetoptError("Must specify basename, printing help.")
                qC = quickCurve(basename, True)
                qC.runCurve(True)         
                return
            elif opt in ('-i','--initialize'):
                print "Creating example configuration file called example.cfg"
                qC = quickCurve(basename)
                qC.writeConfig()
                return
            
        if not opts: raise getopt.GetoptError("Must specify an option, printing help.")
            
    except getopt.error as e:
        print "Command Line Error: " + e.msg
        printCLIHelp()

if __name__ == '__main__': cli()

