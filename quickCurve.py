#!/usr/bin/env python

"""Generate a light curve from Fermi LAT data.

You should have already generated all of the needed files with
quickAnalysis and quickLike including an XML model which describes
your region.  This XML file must have errors properly calculated at
least for your source of interest.

At any time execute

> quickCurve -h 

for help.

First, generate a default config file

> quickCurve initialize

and then edit the config file to match your specific analysis.  Copy
the quickCurve section from the example.cfg into your <BASENAME>.cfg
file that you used for the quickAnalysis and quickLike steps.  For
more options type

> quickCurve initialize -h

You will then need to perform 3 steps to generate a light curve.  Do
the following

> quickCurve run <BASENAME>
> quickCurve compute <BASENAME>
> quickCurve summary <BASENAME>

The first step will generate all of the needed files, the second does
the likelihood calculation for each bin and the final step merges the
results into a final summary file (usually called lc_summary.dat).
For more options and details on all of them execute

> quickCurve run -h
> quickCurve compute -h
> quickCurve summary -h

This module logs all output to a file called <BASENAME>_quickCurve.log.

This code is based on a script writting by S. Fegan."""

__author__ = 'Jeremy S. Perkins (FSSC)'
__version__ = '0.2.0'

import os
import glob
import pickle

import numpy as np
import scipy as sp
import quickAnalysis as qA
import quickLike as qL
import quickUtils as qU
import SummedLikelihood as SL
import UnbinnedAnalysis as UA
import BinnedAnalysis as BA
import IntegralUpperLimit as IUL

from copy import copy
from math import sqrt
from multiprocessing import Pool
from quickUtils import quickMath as MyMath

def runAnalysisStepMP(bininfo):

    bin = bininfo[0]
    tmin = bininfo[1]
    commonConf = bininfo[2]
    analysisConf = bininfo[3]
    curveConf = bininfo[4]
    tmax = tmin + float(curveConf['tstep'])

    print bin,tmin,tmax
    dir = "quickCurve_"+ str(curveConf['tstep']) + "_bin" + str(bin) 

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
                      "convtype" : analysisConf["convtype"],
                      "filter" : analysisConf["filter"],
                      "roicut" : analysisConf["roicut"],
                      "ltzmax" : analysisConf["ltzmax"]}
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
    def sfinit(self, srcName=None, ft2=None, irfs=None, model=None,
                 optimizer="Minuit"):
        self.lc = []

        if(srcName == None):
            return
        self.srcName = srcName
        if model != None:
            self.model = model
        else:
            self.model = srcName + "_fitmodel.xml"
        self.ft2 = ft2
        self.irfs = irfs
        self.optimizer = optimizer
        self.obsfiles = []
    
    def __init__(self,
                 base = 'MySource',
                 configFile = False,
                 curveConfig = {'tstart'  : 0,
                                'tstop'   : 0,
                                'tstep'   : 86400,
                                'tsmin'   : 1,
                                'model'   : 'my_model.xml',
                                'summary' : 'lc_summary.dat',
                                'output'  : 'lc.pickle',
                                'ulfluxdf' : 2.0,
                                'ulbayes'  : 4.0,
                                'ulchi2'  : 4,
                                'ulcl'    : 0.95,
                                'opt'     : 'MINUIT',
                                'sliding' : False,
                                'rebin'   : 1},
                 analysisConfig = {"ra" : 0,
                                   "dec" : 0,
                                   "rad" : 10,
                                   "tmin" : "INDEF",
                                   "tmax" : "INDEF",
                                   "emin" : 100,
                                   "emax" : 300000,
                                   "zmax" : 100,
                                   "ltzmax" : 180,
                                   "binsize" : 0.1,
                                   "convtype" : -1,
                                   "filter" : "DATA_QUAL==1 && LAT_CONFIG==1",
                                   "roicut" : "yes"},
                 likelihoodConfig = {"model" : "MySource_model.xml",
                                     "sourcename" : "Source Name",
                                     "drmtol" : 0.1,
                                     "mintol" : 1e-4},
                 commonConfig = {"base" : 'MySource',
                                 "eventclass" : 2,
                                 "binned" : False,
                                 "irfs" : "P7REP_SOURCE_V15",
                                 "verbosity" : 0}):
        
        commonConfig['base'] = base
        
        self.logger = qU.initLogger(base, 'quickCurve')
        
        if(configFile):
            try:
                commonConfigRead,analysisConfigRead,\
                    likelihoodConfigRead,plotConfigRead,curveConfigRead = qU.readConfig(self.logger,base)
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

        self.model = self.curveConf['model']
        self.ft2 = self.commonConf['base']+"_SC.fits"
        self.irfs = self.commonConf['irfs']
        self.optimizer = self.curveConf['opt']

        self.lc = []
        self.obsfiles=[]
        

    def writeConfig(self):

        """Writes all of the initialization variables to the config
        file called <basename>.cfg."""

        qU.writeConfig(quickLogger=self.logger,
		       curveDictionary=self.curveConf,
		       likelihoodDictionary=self.likelihoodConf,
                       commonDictionary=self.commonConf,
                       analysisDictionary=self.analysisConf)


    def globStandardObsDir(self, directory_glob, analysis='unbinned',
                           ft2=None, irfs=None, nbin = 1,
                           sliding_window = False):

        directories = glob.glob(directory_glob)
        directories.sort()
        if nbin == 1:
            for d in directories:
                if os.path.isdir(d):
                    self.addStandardObsDir(d, ft2, irfs, self.obsfiles,
                                           analysis=analysis)
        elif sliding_window == False:
            if self.obsfiles:
                obslist = self.obsfiles.pop()
            else:
                obslist = list()
            for d in directories:
                if len(obslist) == nbin:
                    self.obsfiles.append(obslist)
                    obslist = list()
                if os.path.isdir(d):
                    self.addStandardObsDir(d, ft2, irfs, obslist,
                                           analysis=analysis)
            if len(obslist) != 0:
                self.obsfiles.append(obslist)
        else:
            if self.obsfiles:
                obslist = self.obsfiles.pop()
            else:
                obslist = list()
            for d in directories:
                if len(obslist) == nbin:
                    self.obsfiles.append(copy(obslist))
                    obslist.pop(0)
                if os.path.isdir(d):
                    self.addStandardObsDir(d, ft2, irfs, obslist,
                                           analysis=analysis)
            if len(obslist) != 0:
                self.obsfiles.append(obslist)

    def addStandardObsDir(self, directory, ft2=None, irfs=None, obslist=None,
                          analysis='unbinned'):
        prefix = directory+"/"+self.commonConf['base']
        ecube = prefix + "_ltcube.fits"
        if analysis=='unbinned':
            ft1 = prefix + "_filtered_gti.fits"
            emap = prefix + "_expMap.fits"
            self.addUnbinnedObs(ft1, emap, ecube, ft2, irfs, obslist)
        elif analysis=='binned':
            smaps = prefix + "_srcMaps.fits"
            bemap = prefix + "_binExpMap.fits"
            self.addBinnedObs(smaps, bemap, ecube, irfs, obslist)
        else:
            raise NameError("Unknown analysis type: \""+f['analysis']+
                            "\" for directory \""+directory+"\"")

    def addUnbinnedObs(self, ft1, emap, ecube,
                       ft2=None, irfs=None, obslist=None):
        if ft2 != None: _ft2 = ft2
        else: _ft2 = self.ft2
        if irfs != None: _irfs = irfs
        else: _irfs = self.irfs
        if not os.path.isfile(ft1):
            raise IOError('FT1 file not found: '+ft1);
        if not os.path.isfile(_ft2):
            raise IOError('FT2 file not found: '+_ft2);
        if not os.path.isfile(emap):
            raise IOError('ExpMap not found: '+emap);
        if not os.path.isfile(ecube):
            raise IOError('ExpCube not found: '+ecube);
        obsfiles = dict(analysis = 'unbinned',
                        ft1       = ft1,
                        ft2       = _ft2,
                        emap      = emap,
                        ecube     = ecube,
                        irfs      = _irfs)
        if obslist == None:
            self.obsfiles.append(obsfiles)
        else:
            obslist.append(obsfiles)

    def addBinnedObs(self, smaps, bemap, ecube, irfs=None, obslist=None):
        if irfs != None: _irfs = irfs
        else: _irfs = self.irfs
        if not os.path.isfile(smaps):
            raise IOError('SourceMaps file not found: '+smaps);
        if not os.path.isfile(bemap):
            raise IOError('Binned ExpMap not found: '+bemap);
        if not os.path.isfile(ecube):
            raise IOError('ExpCube not found: '+ecube);
        obsfiles = dict(analysis  = 'unbinned',
                        smaps     = smaps,
                        bemap     = bemap,
                        ecube     = ecube,
                        irfs      = _irfs)
        if obslist == None:
            self.obsfiles.append(obsfiles)
        else:
            obslist.append(obsfiles)

    def loadUnbinnedObs(self, f, verbosity=0):
        if verbosity:
            print 'Loading unbinned observation:',f['ft1']
        obs = UA.UnbinnedObs(eventFile=f['ft1'], scFile=f['ft2'],
                                           expMap=f['emap'],expCube=f['ecube'],
                                           irfs=f['irfs'])
        like = UA.UnbinnedAnalysis(obs, srcModel=self.model,
                                                 optimizer=self.optimizer)
        return [ obs, like ]

    def loadBinnedObs(self, f, verbosity=0):
        if verbosity:
            print 'Loading binned observation:',f['smaps']
        obs = BA.BinnedObs(srcMaps=f['smaps'], expCube=f['ecube'],
                                       binnedExpMap=f['bemap'], irfs=f['irfs'])
        like = BA.BinnedAnalysis(obs, srcModel=self.model,
                                             optimizer=self.optimizer)
        return [ obs, like ]

    def loadObs(self, f, verbosity=0):
        if f['analysis'] == 'unbinned':
            return self.loadUnbinnedObs(f, verbosity)
        elif f['analysis'] == 'binned':
            return self.loadBinnedObs(f, verbosity)
        else:
            raise NameError("Unknown analysis type: \""+f['analysis']+"\"")

    def processAllObs(self, fix_shape=True, delete_below_ts=None,
                      ul_flux_dflux=0, ul_chi2_ts=None, ul_bayes_ts=4.0,
                      ul_cl=0.95, verbosity=0, emin=0, emax=0,
                      interim_save_filename=None):

        for f in self.obsfiles:
            lc = dict()
            lc['config'] = dict()
            lc['config']['fix_shape']       = fix_shape
            lc['config']['delete_below_ts'] = delete_below_ts
            lc['config']['ul_flux_dflux']   = ul_flux_dflux
            lc['config']['ul_chi2_ts']      = ul_chi2_ts
            lc['config']['ul_bayes_ts']     = ul_bayes_ts
            lc['config']['ul_cl']           = ul_cl
            lc['config']['emin']            = emin
            lc['config']['emax']            = emax
            lc['config']['files']           = f
            #lc['config']['argv']            = sys.argv

            lc['e_min'] = emin;
            lc['e_max'] = emax;

            if type(f) != list:
                [ obs, like ] = self.loadObs(f,verbosity)
                lc['t_min'] = obs.roiCuts().minTime()
                lc['t_max'] = obs.roiCuts().maxTime()
                if (emin == 0 or emax == 0):
                    lc['e_min'] = obs.roiCuts().getEnergyCuts()[0];
                    lc['e_max'] = obs.roiCuts().getEnergyCuts()[1];

            else:
                lc['t_min'] = None
                lc['t_max'] = None
                like = SL.SummedLikelihood(self.optimizer)
                for ff in f:
                    [ obs, like1 ] = self.loadObs(ff,verbosity)
                    tmin = obs.roiCuts().minTime()
                    tmax = obs.roiCuts().maxTime()
                    if lc['t_min'] == None or tmin<lc['t_min']:
                        lc['t_min'] = tmin
                    if lc['t_max'] == None or tmax>lc['t_max']:
                        lc['t_max'] = tmax
                    if (lc['e_min'] == 0 or lc['e_max'] == 0):
                        ecuts = obs.roiCuts().getEnergyCuts()
                        lc['e_min'] = ecuts[0]
                        lc['e_max'] = ecuts[1]
                    elif (emin == 0 or emax == 0):
                        ecuts = obs.roiCuts().getEnergyCuts()
                        lc['e_min'] = max(lc['e_min'], ecuts[0])
                        lc['e_max'] = min(lc['e_max'], ecuts[1])
                    like.addComponent(like1)

            emin = lc['e_min']
            emax = lc['e_max']

            like.tol = like.tol*0.01;

            if verbosity > 1:
                print '- Time:',lc['t_min'],'to',lc['t_max']

            src = like[self.likelihoodConf['sourcename']]
            if src == None:
                raise NameError("No source \""+self.likelihoodConf['sourcename']+"\" in model "+
                                self.model)
            srcfreepar=like.freePars(self.likelihoodConf['sourcename'])
            srcnormpar=like.normPar(self.likelihoodConf['sourcename'])
            if len(srcfreepar)>0:
                like.setFreeFlag(self.likelihoodConf['sourcename'], srcfreepar, 0)
                like.syncSrcParams(self.likelihoodConf['sourcename'])


            meanvalue = srcnormpar.getValue()
            meanerror = srcnormpar.error()
            if meanerror == 0:
                self.logger.critical("The error on the normalization for your source is 0!  You need to do a global fit first (with quickLike) and provide the final XML file (<basename>_likeMinuit.xml) with errors included before you run compute.")
                return

            lc['original']=dict()
            lc['original']['normpar_init_value'] = meanvalue
            lc['original']['normpar_name'] = srcnormpar.getName()
            lc['original']['nfree'] = len(like.freePars(self.likelihoodConf['sourcename']))
            lc['original']['flux'] = like[self.likelihoodConf['sourcename']].flux(emin, emax)
            lc['original']['logL'] = like.logLike.value()
            if verbosity > 1:
                print '- Original log Like:',lc['original']['logL']

            if fix_shape:
                if verbosity > 1:
                    print '- Fixing spectral shape parameters'
                sync_name = ""
                for p in like.params():
                    if sync_name != "" and sync_name != p.srcName:
                        like.syncSrcParams(sync_name)
                        sync_name = ""
                    if(p.isFree() and p.srcName!=self.likelihoodConf['sourcename'] and
                       p.getName()!=like.normPar(p.srcName).getName()):
                        if verbosity > 2:
                            print '-- '+p.srcName+'.'+p.getName()
                        p.setFree(False)
                        sync_name = p.srcName
                if sync_name != "" and sync_name != p.srcName:
                    like.syncSrcParams(sync_name)
                    sync_name = ""

           # ----------------------------- FIT 1 -----------------------------

            if verbosity > 1:
                print '- Fit 1 - All parameters of',self.likelihoodConf['sourcename'],'fixed'
            like.fit(max(verbosity-3, 0))

            lc['allfixed'] = dict()
            lc['allfixed']['logL'] = like.logLike.value()
            fitstat = like.optObject.getRetCode()
            if verbosity > 1 and fitstat != 0:
                print "- Fit 1 - Minimizer returned with code: ", fitstat
            lc['allfixed']['fitstat'] = fitstat
            if verbosity > 1:
                print '- Fit 1 - log Like:',lc['allfixed']['logL']

            if delete_below_ts:
                frozensrc = []
                if verbosity > 1:
                    print '- Deleting point sources with TS<'+str(delete_below_ts)
                deletesrc = []
                for s in like.sourceNames():
                    freepars = like.freePars(s)
                    if(s!=self.likelihoodConf['sourcename'] and like[s].type == 'PointSource'
                       and len(freepars)>0):
                        ts = like.Ts(s)
                        if ts<delete_below_ts:
                            deletesrc.append(s)
                            if verbosity > 2:
                                print '--',s,'(TS='+str(ts)+')'
                if deletesrc:
                    for s in deletesrc:
                        like.deleteSource(s)
                    if verbosity > 1:
                        print '- Fit 1 - refitting model'
                    like.fit(max(verbosity-3, 0))
                    lc['allfixed']['fitstat_initial'] = \
                        lc['allfixed']['fitstat']
                    fitstat = like.optObject.getRetCode()
                    if verbosity > 1 and fitstat != 0:
                        print "- Fit 1 - Minimizer returned with code: ",\
                            fitstat
                    lc['allfixed']['fitstat'] = fitstat
                    lc['allfixed']['logL'] = like.logLike.value()
                    if verbosity > 1:
                        print '- Fit 1 - log Like:',lc['allfixed']['logL']


            lc['allfixed']['flux']=like[self.likelihoodConf['sourcename']].flux(emin, emax)
            pars = dict()
            for pn in like[self.likelihoodConf['sourcename']].funcs['Spectrum'].paramNames:
                p = like[self.likelihoodConf['sourcename']].funcs['Spectrum'].getParam(pn)
                pars[p.getName()] = dict(name      = p.getName(),
                                         value     = p.getTrueValue(),
                                         error     = p.error()*p.getScale(),
                                         free      = p.isFree())
            lc['allfixed']['pars'] = pars
    

            # ------------------ N SIGMA PROFILE LIKELIHOOD -------------------

            prof_sigma = (-1,-0.5,0,0.5,1.0)
            lc['profile'] = dict();
            lc['profile']['sigma'] = []
            lc['profile']['value'] = []
            lc['profile']['logL'] = []
            lc['profile']['flux'] = []
            lc['profile']['fitstat'] = []

            if verbosity > 1:
                print '- Fit 1 - generating %d point likelihood profile'%\
                      len(prof_sigma)
            for sigma in prof_sigma:
                val = sigma*meanerror+meanvalue
                if val < srcnormpar.getBounds()[0]:
                    val = srcnormpar.getBounds()[0]
                if (lc['profile']['value']
                    and lc['profile']['value'][-1]==val):
                    continue
                lc['profile']['value'].append(val)
                lc['profile']['sigma'].append((val-meanvalue)/meanerror)
                if(val == meanvalue):
                    lc['profile']['logL'].append(lc['allfixed']['logL'])
                    lc['profile']['flux'].append(lc['allfixed']['flux'])
                else:
                    srcnormpar.setValue(val)
                    like.syncSrcParams(self.likelihoodConf['sourcename'])
                    like.fit(max(verbosity-3, 0))
                    fitstat = like.optObject.getRetCode()
                    if verbosity > 2 and fitstat != 0:
                        print "- Fit 1 - profile: Minimizer returned code: ",\
                            fitstat
                    lc['profile']['fitstat'].append(fitstat)
                    lc['profile']['logL'].append(like.logLike.value())
                    lc['profile']['flux'].append(like[self.likelihoodConf['sourcename']].\
                                              flux(emin, emax))
                if verbosity > 2:
                    print '- Fit 1 - profile: %+g, %f -> %f'%\
                          (sigma,lc['profile']['value'][-1],
                           lc['profile']['logL'][-1]-lc['allfixed']['logL'])

            srcnormpar.setValue(meanvalue)
            like.syncSrcParams(self.likelihoodConf['sourcename'])

            # ----------------------------- FIT 2 -----------------------------

            if verbosity > 1:
                print '- Fit 2 - Normalization parameter of',\
                      self.likelihoodConf['sourcename'],'free'
            srcnormpar.setFree(1)
            like.syncSrcParams(self.likelihoodConf['sourcename'])
            like.fit(max(verbosity-3, 0))
            lc['normfree'] = dict()
            fitstat = like.optObject.getRetCode()
            if verbosity > 1 and fitstat != 0:
                print "- Fit 2 - Minimizer returned with code: ", fitstat
            lc['normfree']['fitstat'] = fitstat
            lc['normfree']['logL'] = like.logLike.value()
            lc['normfree']['ts'] = like.Ts(self.likelihoodConf['sourcename'])
            lc['normfree']['flux_dflux'] = \
                srcnormpar.getValue()/srcnormpar.error()
            if verbosity > 1:
                print '- Fit 2 - log Like:',lc['normfree']['logL'],\
                      '(TS='+str(lc['normfree']['ts'])+')'

            lc['normfree']['nfree']=len(like.freePars(self.likelihoodConf['sourcename']))
            lc['normfree']['flux']=like[self.likelihoodConf['sourcename']].flux(emin, emax)
            pars = dict()
            for pn in like[self.likelihoodConf['sourcename']].funcs['Spectrum'].paramNames:
                p = like[self.likelihoodConf['sourcename']].funcs['Spectrum'].getParam(pn)
                pars[p.getName()] = dict(name      = p.getName(),
                                         value     = p.getTrueValue(),
                                         error     = p.error()*p.getScale(),
                                         free      = p.isFree())
            lc['normfree']['pars'] = pars
            ul_type = None
            if ul_bayes_ts != None and lc['normfree']['ts'] < ul_bayes_ts:

                ul_type = 'bayesian'
                [ul_flux, ul_results] = \
                    IUL.calc_int(like,self.likelihoodConf['sourcename'],cl=ul_cl,
                                                skip_global_opt=True,
                                                verbosity = max(verbosity-2,0),
                                                emin=emin, emax=emax,
                                            poi_values = lc['profile']['value'])
            elif ( ul_flux_dflux != None and \
                   lc['normfree']['flux_dflux'] < ul_flux_dflux ) or \
                   ( ul_chi2_ts != None and lc['normfree']['ts'] < ul_chi2_ts):
                ul_type = 'chi2'
                [ul_flux, ul_results] = \
                    IUL.calc_chi2(like,self.likelihoodConf['sourcename'],cl=ul_cl,
                                                 skip_global_opt=True,
                                                 verbosity = max(verbosity-2,0),
                                                 emin=emin, emax=emax)
            if ul_type != None:
                lc['normfree']['ul'] = dict(flux    = ul_flux,
                                            results = ul_results,
                                            type    = ul_type)

            # ----------------------------- FIT 3 -----------------------------

            if verbosity > 1:
                print '- Fit 3 - All parameters of',self.likelihoodConf['sourcename'],'free'
            like.setFreeFlag(self.likelihoodConf['sourcename'], srcfreepar, 1)
            like.syncSrcParams(self.likelihoodConf['sourcename'])
            like.fit(max(verbosity-3, 0))
            lc['allfree'] = dict()
            fitstat = like.optObject.getRetCode()
            if verbosity > 1 and fitstat != 0:
                print "- Fit 3 - Minimizer returned with code: ", fitstat
            lc['allfree']['fitstat'] = fitstat
            lc['allfree']['logL'] = like.logLike.value()
            lc['allfree']['ts'] = like.Ts(self.likelihoodConf['sourcename'])
            if verbosity > 1:
                print '- Fit 3 - log Like:',lc['allfree']['logL'],\
                      '(TS='+str(lc['allfree']['ts'])+')'
            lc['allfree']['nfree']=len(like.freePars(self.likelihoodConf['sourcename']))
            lc['allfree']['flux']=like[self.likelihoodConf['sourcename']].flux(emin, emax)
            pars = dict()
            for pn in like[self.likelihoodConf['sourcename']].funcs['Spectrum'].paramNames:
                p = like[self.likelihoodConf['sourcename']].funcs['Spectrum'].getParam(pn)
                pars[p.getName()] = dict(name      = p.getName(),
                                         value     = p.getTrueValue(),
                                         error     = p.error()*p.getScale(),
                                         free      = p.isFree())
            lc['allfree']['pars'] = pars

            self.lc.append(lc)
            if interim_save_filename != None:
                self.saveProcessedObs(interim_save_filename)

    def saveProcessedObs(self,filename):
        file=open(filename,'w')
        pickle.dump(self.lc,file)

    def loadProcessedObs(self,filename):
        file=open(filename,'r')
        lcs=pickle.load(file)
        for lc in lcs:
            self.lc.append(lc)

    def generateLC(self, verbosity=0):
        # First: calculate logL of fixed flux model at true minimum - hoping
        # it lies somewhere in the profile we computed
        first = True
        profile_x = []
        profile_y = []
        for lc in self.lc:
            if first:
                profile_x = lc['profile']['value']
                profile_y = lc['profile']['logL']
                first = False
            else:
                profile_y = map(lambda x,y:x+y,lc['profile']['logL'],profile_y)

        p = sp.polyfit(profile_x, profile_y, 2);
        prof_max_val = -p[1]/(2*p[0])
        prof_max_logL = p[2]-p[1]*p[1]/(4*p[0])
        if (prof_max_val<min(profile_x)) or (prof_max_val>max(profile_x)):
            print "Warning: corrected minimum %f is outside profile range [%f to %f]" \
                  %(prof_max_val,min(profile_x),max(profile_x))
            print profile_x, profile_y

        profile_fity = sp.polyval(p,profile_x)
        profile_max_diff = max(map(lambda x,y:abs(x-y),profile_y,profile_fity))
        if profile_max_diff>0.5:
            print "Warning: large difference between profile and fit: %f"%profile_max_diff
            print profile_x, profile_y, profile_fity

        if verbosity>0:
            print profile_x, profile_y, profile_fity

        # Second: process data for LC, accumulating required values to
        # allow calculation of variability stats

        vals = []
        pars = []
        dchi2_specfree     = 0
        dchi2_normfree     = 0
        dchi2_normfree_alt = 0
        dchi2_normfree_ul  = 0
        allfixed_logL      = 0
        allfixed_val       = 0
        npar_specfree      = 0
        npar_normfree      = 0
        first = True

        for lc in self.lc:
            np=lc['original']['normpar_name']
            if first:
                pars.append(np)
                allfixed_val = lc['original']['normpar_init_value']
            scale=lc['normfree']['flux']/lc['normfree']['pars'][np]['value']
            val = []
            val.append(lc['t_min']/86400 + 51910)
            val.append(lc['t_max']/86400 + 51910)
            if lc['normfree'].has_key('ul'):
                val.append(lc['normfree']['ul']['flux'])
                val.append(0)
            else:
                val.append(lc['normfree']['flux'])
                val.append(lc['normfree']['pars'][np]['error']*scale)
            val.append(lc['normfree']['ts'])
            val.append(lc['allfree']['flux'])
            val.append(lc['allfree']['pars'][np]['error']*scale)
            for p in lc['allfree']['pars']:
                if p != np and lc['allfree']['pars'][p]['free'] == True:
                    if first: pars.append(p)
                    val.append(lc['allfree']['pars'][p]['value'])
                    val.append(lc['allfree']['pars'][p]['error'])
            val.append(lc['allfree']['ts'])

            allfixed_logL += lc['allfixed']['logL']
            dchi2_specfree += 2*(lc['allfree']['logL']-lc['normfree']['logL'])
            dchi2_normfree_alt += lc['normfree']['logL']

            # Arbitrarily assume a quadratic is an OK fit
            y = lc['profile']['logL']
            p = sp.polyfit(profile_x, y, 2);
            dchi2_normfree += 2*(lc['normfree']['logL']
                                 - sp.polyval(p, prof_max_val))

            if (lc['normfree'].has_key('ul') and
                lc['normfree']['ul']['type'] == 'bayesian'):
                # Arbitrarily assume a quadratic is an OK fit
                y = lc['normfree']['ul']['results']['poi_chi2_equiv']
                p = sp.polyfit(profile_x, y, 2);
                dchi2_normfree_ul += sp.polyval(p, prof_max_val)
            else:
                dchi2_normfree_ul += 2*(lc['normfree']['logL']
                                        - sp.polyval(p, prof_max_val))
            if not first:
                npar_specfree += lc['allfree']['nfree']-lc['normfree']['nfree']
                npar_normfree += lc['normfree']['nfree']
            first = False
            
            vals.append(val)

        dchi2_normfree_alt = 2*(dchi2_normfree_alt-prof_max_logL)
        corr_logL = prof_max_logL - allfixed_logL

        if(abs(dchi2_normfree-dchi2_normfree_alt) > 0.01):
            print "Warning: normfree log likelhood calculations differ by more than 0.01"
            print dchi2_normfree, dchi2_normfree_alt, dchi2_normfree_ul

        stats = dict(dchi2_specfree            = dchi2_specfree,
                     dchi2_normfree            = dchi2_normfree,
                     dchi2_normfree_ul         = dchi2_normfree_ul,
                     npar_specfree             = npar_specfree,
                     npar_normfree             = npar_normfree,
                     pars                      = pars,
                     prof_x                    = profile_x,
                     prof_y                    = profile_y,
                     prof_max_val              = prof_max_val,
                     prof_max_logL             = prof_max_logL,
                     prof_corr_logL            = corr_logL,
                     allfixed_val              = allfixed_val,
                     allfixed_logL             = allfixed_logL)
        return vals, stats

    
    def writeLC(self, filename=None, lc=None, stats=None,
                header=True, headstart='% ', verbosity=0):
        if lc == None or stats == None:
            [lc, stats] = self.generateLC(verbosity=verbosity)
        #file = sys.stdout
        if filename != None:
            file=open(filename,'w')
        if header:
            # print >>file, '%sOptions: %s'%(headstart,' '.join(lc[0]['config']['argv'][1:]))
            chi2 = stats['dchi2_normfree']
            ndof = stats['npar_normfree']
            try:
                prob = MyMath.chi2cdfc(chi2,ndof)
            except ValueError:
                self.logger.critical("Chi^2 Probability not well defined.  Setting to 0.")
                prob = 0.
            sigma = sqrt(MyMath.chi2invc(prob,1))
            print >>file, '%sVariable flux (no UL): chi^2=%.3f (%d DOF) - Pr(>X)=%g (~%g sigma)'%(headstart,chi2,ndof,prob,sigma)
            chi2 = stats['dchi2_normfree_ul']
            ndof = stats['npar_normfree']
            try:
                prob = MyMath.chi2cdfc(chi2,ndof)
            except ValueError:
                self.logger.critical("Chi^2 Probability not well defined.  Setting to 0.")
                prob = 0.
            sigma = sqrt(MyMath.chi2invc(prob,1))
            print >>file, '%sVariable flux (w/UL):  chi^2=%.3f (%d DOF) - Pr(>X)=%g (~%g sigma)'%(headstart,chi2,ndof,prob,sigma)
            chi2 = stats['dchi2_specfree']
            ndof = stats['npar_specfree']
            try:
                prob = MyMath.chi2cdfc(chi2,ndof)
            except ValueError:
                self.logger.critical("Chi^2 Probability not well defined.  Setting to 0.")
                prob = 0.
            sigma = sqrt(MyMath.chi2invc(prob,1))
            print >>file, '%sVariable spectrum:     chi^2=%.3f (%d DOF) - Pr(>X)=%g (~%g sigma)'%(headstart,chi2,ndof,prob,sigma)
            print >>file, '%sProfile minimum: %f (search range: %f to %f)'%(headstart,stats['prof_max_val'],min(stats['prof_x']),max(stats['prof_x']))
            print >>file, '%sLogL correction: %f (WRT logL @ prescribed val of %g)'%(headstart,stats['prof_corr_logL'],stats['allfixed_val'])

            print >>file, '%sColumn 1: Start of time bin [MJD]'%(headstart)
            print >>file, '%sColumn 2: End of time bin [MJD]'%(headstart)
            print >>file, '%sColumn 3: Fixed spectral shape: Flux [ph/cm^2/s]'%(headstart)
            print >>file, '%sColumn 4: Fixed spectral shape: Error on Flux [ph/cm^2/s]'%(headstart)
            print >>file, '%sColumn 5: Fixed spectral shape: TS'%(headstart)
            print >>file, '%sColumn 6: Optimized spectral shape: Flux [ph/cm^2/s]'%(headstart)
            print >>file, '%sColumn 7: Optimized spectral shape: Error on Flux [ph/cm^2/s]'%(headstart)
            nc=8
            for i in range(1,len(stats['pars'])):
                pn = stats['pars'][i]
                print >>file, '%sColumn %d: Optimized spectral shape: %s'%(headstart,nc+i-1,pn)
                print >>file, '%sColumn %d: Optimized spectral shape: Error on %s'%(headstart,nc+i,pn)
                nc+=2
            print >>file, '%sColumn %d: Optimized spectral shape: TS'%(headstart,nc+i-1)
        for p in lc:
            s = '%.3f %.3f %.3e %.3e %7.2f'%(p[0],p[1],p[2],p[3],p[4])
            for i in range(5,len(p)-1):
                s += ' %.3e'%(p[i])
            s += ' %7.2f'%(p[-1])
            print >>file, s

    def runCurve(self, runAnalysis=True, delete = False):

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


        if(delete):
            templist = glob.glob("*_bin" + str(binnum) + "*")
            for t in templist:
                os.remove(t)

def overrideConfig(logger,dictionary,argVars):

    for variable, value in dictionary.iteritems():
        if variable in argVars:
            if argVars[variable] != None:
                dictionary[variable] = argVars[variable]

def cli():

    from argparse import ArgumentParser, RawTextHelpFormatter
    
    parser = ArgumentParser(description = "                 - quickCurve - \n\n"+
                 "Compute lightcurves from Fermi data. The program opeartes in\n"+
                 "three modes:\n\n\tinitialize, run, summary and compute,\n\n"+
                 "specified with the run, summary or compute options.  The run\n"+
                 "mode generates all of the needed files for the next two modes\n"+
                 "and puts them in seperate directories in the working directory\n"+
                 "(named <basename>_binX). In the compute mode one or many Fermi\n"+
                 "observations are analyzed using the pyLikelihood tools to produce\n"+
                 "a summary file. In the summary mode, these summary files are\n"+
                 "read and the lightcurve is produced.  All of the options can\n"+
                 "be stored in a config file which can be read if you use the\n"+
                 "--config option.",
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument("--verbosity", type=int,
                        help="Verbosity (1,2 or 3)")

    subparsers = parser.add_subparsers(dest="mode")

    init_parser = subparsers.add_parser('initialize', 
                                        help= "Generate a default config file called <BASENAME>.cfg.\n"+
                                        "CAREFUL, it will overwrite the current file.")
    init_parser.add_argument("--basename", type=str,
                             help = "Name the example config file <BASENAME>.cfg instead of\n"+
                                    "example.cfg.")

    run_parser = subparsers.add_parser('run', help="Generate all of the needed files for the lightcurve\n"+
                                       "analysis.  You must already have a config file if\n"+
                                       "using the command line interface.")
    run_parser.add_argument("basename", type=str,
                             help = "basename of the analysis.  The config file should be named\n"+
                                    "<BASENAME>.cfg and all output files will have this prefix.\n"+
                                    "All of the parameters will be read from this config file \n"+
                                    "but you can ovveride any of them from the command line.")
    run_parser.add_argument("--tstart", type=float,
                            help = "Lightcurve start time (in MET)")
    run_parser.add_argument("--tstop", type=float,
                            help = "Lightcurve stop time (in MET)")
    run_parser.add_argument("--tstep", type=float,
                            help = "Lightcurve step size (in seconds, default is 86400)")


    compute_parser = subparsers.add_parser('compute', help="The files produced in the run mode reanalyzed using\n"+
                                           "the pyLikelihood tools to produce a summary file.")
    compute_parser.add_argument("basename", type=str,
                             help = "basename of the analysis.  The config file should be named\n"+
                                    "<BASENAME>.cfg and all output files will have this prefix.\n"+
                                    "All of the parameters will be read from this config file \n"+
                                    "but you can ovveride any of them from the command line.")
    compute_parser.add_argument("--tsmin", type=float, 
                                help = "TS value below which background sources \n"+
                                "are deleted from the model (default is 1).")
    compute_parser.add_argument("--ulbayes", type=float,
                                help = "TS value below which the Baysian upper \n"+
                                "limit is computed (default is 4).")
    compute_parser.add_argument("--ulchi2", type=float,
                                help = "TS value below which the Profile Likelihood\n"+
                                "upper limit is computed (default is 4).")
    compute_parser.add_argument("--ulfluxdf", type=float,
                                help = "Set the value of the flux/error below which\n"+
                                "a Profile likelihood upper limit is calculated (unless\n"+
                                "it is preempted by the Bayes method based on TS value)\n"+
                                "(default is 2)")
    compute_parser.add_argument("--ulcl", type=float,
                                help = "Set the confidence limit of upper limits (default\n"+
                                "is 0.95)")
    compute_parser.add_argument("--model", type=str,
                                help = "The filename of the XML model from the full fit.\n"+
                                "Note that this must be the output of a fit.  It needs to \n"+
                                "have error information on at least your source of interest\n"+
                                "(default = <BASENAME>_likeMinuit.xml)")
    compute_parser.add_argument("--rebin", type=int,
                                help = "Combine <REBIN> time bins into one (default is 1)")
    compute_parser.add_argument("--sliding", type=bool,
                                help = "Combine the time bins using a sliding window\n"+
                                "so that they overlap (default is False)")
    compute_parser.add_argument("--output", type=str,
                        help="Output file name (default is 'lc.pickle')")

    summary_parser = subparsers.add_parser('summary', help="Generate a light curve from the likelihood computations\n"+
                                           "performed by the 'compute' method.")
    summary_parser.add_argument("basename", type=str,
                             help = "basename of the analysis.  The config file should be named\n"+
                                    "<BASENAME>.cfg and all output files will have this prefix.\n"+
                                    "All of the parameters will be read from this config file \n"+
                                    "but you can ovveride any of them from the command line.")
    summary_parser.add_argument("--summary", type=str,
                                help="Output file name (default is lc_summary.dat)")


    args = parser.parse_args()

    if args.mode == 'initialize':
        print "Creating example config file named example.cfg..."
        if args.basename:
            qC = quickCurve(args.basename)
        else:
            qC = quickCurve('example')
        qC.writeConfig()
        return
    elif args.mode == 'run':
        qC = quickCurve(args.basename, True)
        qC.logger.info("Generating files...")
        argVars = vars(args)
        overrideConfig(qC.logger,qC.commonConf,argVars)
        overrideConfig(qC.logger,qC.likelihoodConf,argVars)
        overrideConfig(qC.logger,qC.curveConf,argVars)
        qC.runCurve(True,False)
        return
    elif args.mode == 'compute':
        qC = quickCurve(args.basename, True)
        qC.logger.info("Computing likelihoood...")
        argVars = vars(args)
        overrideConfig(qC.logger,qC.commonConf,argVars)
        overrideConfig(qC.logger,qC.likelihoodConf,argVars)
        overrideConfig(qC.logger,qC.curveConf,argVars)
        
        dirs = glob.glob('quickCurve_'+str(qC.curveConf['tstep'])+'_bin*')
        if qC.commonConf['binned']:
            analysis = 'binned'
        else:
            analysis = 'unbinned'
        for d in dirs:
            qC.globStandardObsDir(d, 
                                  nbin=int(qC.curveConf['rebin']), 
                                  analysis=analysis,
                                  sliding_window = qC.curveConf['sliding'])
        if(qC.curveConf['ulchi2']<0): 
            qC.curveConf['ulchi2']=None
            qC.logger.info("Profile likelihood upper limit disabled.")
        if(qC.curveConf['ulbayes']<0): 
            qC.curveConf['ulbayes']=None
            qC.logger.info("Baysian upper limit disabled.")
        qC.processAllObs(verbosity=int(qC.commonConf['verbosity']), 
                         delete_below_ts=float(qC.curveConf['tsmin']),
                         ul_chi2_ts=float(qC.curveConf['ulchi2']), 
                         ul_flux_dflux = float(qC.curveConf['ulfluxdf']),
                         ul_bayes_ts=float(qC.curveConf['ulbayes']), 
                         ul_cl=float(qC.curveConf['ulcl']),
                         interim_save_filename=qC.curveConf['output'])
        return
    elif args.mode == 'summary':
        qC = quickCurve(args.basename, True)
        argVars = vars(args)
        overrideConfig(qC.logger,qC.commonConf,argVars)
        overrideConfig(qC.logger,qC.likelihoodConf,argVars)
        overrideConfig(qC.logger,qC.curveConf,argVars)
        qC.loadProcessedObs(qC.curveConf['output'])
        qC.writeLC(qC.curveConf['summary'],verbosity=qC.commonConf['verbosity'])
                            

if __name__ == '__main__': 
    cli()

