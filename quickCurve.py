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
from UpperLimits import UpperLimits
from gt_apps import *
import glob

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


        
    def runAnalysisStep(self,bin=0,tmin=0,tmax=0,delete=True):
            
        basename = self.commonConf['base'] + "_bin" + str(bin)

        os.symlink(self.commonConf['base'] + ".list", basename+".list")        
        os.symlink(self.commonConf['base'] + "_SC.fits", basename+"_SC.fits")        

        analysisConfig = {"ra" : self.analysisConf["ra"],
                          "dec" : self.analysisConf["dec"],
                          "rad" : self.analysisConf["rad"],
                          "tmin" : tmin,
                          "tmax" : tmax,
                          "emin" : self.analysisConf["emin"],
                          "emax" : self.analysisConf["emax"],
                          "zmax" : self.analysisConf["zmax"],
                          "binsize" : self.analysisConf["binsize"],
                          "convtype" : self.analysisConf["convtype"]}
        commonConfig = {"base" : basename,
                        "eventclass" : self.commonConf["eventclass"],
                        "binned" : False,
                        "irfs" : self.commonConf["irfs"],
                        "verbosity" : self.commonConf["verbosity"],
                        "multicore" : 0}


        qA_bin = qA.quickAnalysis(basename,False,analysisConfig,commonConfig)
        qA_bin.analysisConf['tmin'] = tmin
        qA_bin.analysisConf['tmax'] = tmax

        qA_bin.runSelect(True)
        qA_bin.runGTI(True)
        qA_bin.runLTCube(True)
        qA_bin.runExpMap(True)

        
    def runLikelihoodStep(self,bin=0,tmin=0,tmax=0,tslimit=4.0):

        basename = self.commonConf['base'] + "_bin" + str(bin)

        likelihoodConfig = {"model" : self.likelihoodConf["model"],
                            "sourcename" : self.likelihoodConf["sourcename"],
                            "drmtol" : self.likelihoodConf["drmtol"],
                            "mintol" : self.likelihoodConf["mintol"]}
        commonConfig = {"base" : basename,
                        "eventclass" : self.commonConf["eventclass"],
                        "binned" : False,
                        "irfs" : self.commonConf["irfs"],
                        "verbosity" : self.commonConf["verbosity"]}

        qL_bin = qL.quickLike(basename,False,likelihoodConfig,commonConfig)

        qL_bin.makeObs()
        qL_bin.initMIN(modelFile=self.curveConf['model'])
        qL_bin.fitMIN()


        sourcename = qL_bin.likelihoodConf['sourcename']

        if(qL_bin.MIN.Ts(sourcename) < tslimit):
            print "Will calulate upper limit for {}".format(sourcename)
            ul = UpperLimits(qL_bin.MIN)
            upper = ul[sourcename].compute(emin=float(self.analysisConf["emin"]),
                                           emax=float(self.analysisConf["emax"]))
        else:
            upper = (0.0,0.0)
            
        return "{} {} {} {:,.2f} {:,.2f} {:,.2e} {:,.2e} {:,.2e} {}".format(bin,tmin,tmax,
                                                                            qL_bin.MIN.Ts(sourcename),
                                                                            qL_bin.MIN.NpredValue(sourcename),
                                                                            qL_bin.MIN.flux(sourcename,
                                                                                            float(self.analysisConf["emin"]),
                                                                                            float(self.analysisConf["emax"])),
                                                                            qL_bin.MIN.fluxError(sourcename,
                                                                                                 float(self.analysisConf["emin"]),
                                                                                                 float(self.analysisConf["emax"])),
                                                                            upper[0],
                                                                            qL_bin.MINobj.getRetCode())


    def runCurve(self,runAnalysis=True,runLike=True, delete = True):

        tbins = np.arange(float(self.curveConf['tstart']),
                          float(self.curveConf['tstop']),
                          float(self.curveConf['tstep']))

        if(runLike):
            filename = self.commonConf['base'] + ".lc"
            f = open(filename, 'w')
            f.write("#bin tmin tmax TS NPred Flux FluxErr Upper FitStatus\n")

        for binnum,t in enumerate(tbins):
            if(runAnalysis):
                self.runAnalysisStep(binnum,t,t+float(self.curveConf['tstep']),delete=True)
            if(runLike):
                output = self.runLikelihoodStep(binnum,t,t+float(self.curveConf['tstep']))
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

