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

        self.lc = []
        self.obsfiles = []

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
                    self.obsfiles.append(copy.copy(obslist))
                    obslist.pop(0)
                if os.path.isdir(d):
                    self.addStandardObsDir(d, ft2, irfs, obslist,
                                           analysis=analysis)
            if len(obslist) != 0:
                self.obsfiles.append(obslist)

    def addStandardObsDir(self, directory, ft2=None, irfs=None, obslist=None,
                          analysis='unbinned'):
        prefix = directory+"/"+self.srcName
        ecube = prefix + "_ltcube.fits"
        if analysis=='unbinned':
            #ft1 = prefix + "_ev_roi.fits"
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
        obs = UnbinnedAnalysis.UnbinnedObs(eventFile=f['ft1'], scFile=f['ft2'],
                                           expMap=f['emap'],expCube=f['ecube'],
                                           irfs=f['irfs'])
        like = UnbinnedAnalysis.UnbinnedAnalysis(obs, srcModel=self.model,
                                                 optimizer=self.optimizer)
        return [ obs, like ]

    def loadBinnedObs(self, f, verbosity=0):
        if verbosity:
            print 'Loading binned observation:',f['smaps']
        obs = BinnedAnalysis.BinnedObs(srcMaps=f['smaps'], expCube=f['ecube'],
                                       binnedExpMap=f['bemap'], irfs=f['irfs'])
        like = BinnedAnalysis.BinnedAnalysis(obs, srcModel=self.model,
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
            lc['version'] = self.ver
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
            lc['config']['argv']            = sys.argv

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
                like = SummedLikelihood.SummedLikelihood(self.optimizer)
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

            src = like[self.srcName]
            if src == None:
                raise NameError("No source \""+self.srcName+"\" in model "+
                                self.model)
            srcfreepar=like.freePars(self.srcName)
            srcnormpar=like.normPar(self.srcName)
            if len(srcfreepar)>0:
                like.setFreeFlag(self.srcName, srcfreepar, 0)
                like.syncSrcParams(self.srcName)


            meanvalue = srcnormpar.getValue()
            meanerror = srcnormpar.error()
            lc['original']=dict()
            lc['original']['normpar_init_value'] = meanvalue
            lc['original']['normpar_name'] = srcnormpar.getName()
            lc['original']['nfree'] = len(like.freePars(self.srcName))
            lc['original']['flux'] = like[self.srcName].flux(emin, emax)
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
                    if(p.isFree() and p.srcName!=self.srcName and
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
                print '- Fit 1 - All parameters of',self.srcName,'fixed'
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
                    if(s!=self.srcName and like[s].type == 'PointSource'
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


            lc['allfixed']['flux']=like[self.srcName].flux(emin, emax)
            pars = dict()
            for pn in like[self.srcName].funcs['Spectrum'].paramNames:
                p = like[self.srcName].funcs['Spectrum'].getParam(pn)
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
                    like.syncSrcParams(self.srcName)
                    like.fit(max(verbosity-3, 0))
                    fitstat = like.optObject.getRetCode()
                    if verbosity > 2 and fitstat != 0:
                        print "- Fit 1 - profile: Minimizer returned code: ",\
                            fitstat
                    lc['profile']['fitstat'].append(fitstat)
                    lc['profile']['logL'].append(like.logLike.value())
                    lc['profile']['flux'].append(like[self.srcName].\
                                              flux(emin, emax))
                if verbosity > 2:
                    print '- Fit 1 - profile: %+g, %f -> %f'%\
                          (sigma,lc['profile']['value'][-1],
                           lc['profile']['logL'][-1]-lc['allfixed']['logL'])

            srcnormpar.setValue(meanvalue)
            like.syncSrcParams(self.srcName)

            # ----------------------------- FIT 2 -----------------------------

            if verbosity > 1:
                print '- Fit 2 - Normalization parameter of',\
                      self.srcName,'free'
            srcnormpar.setFree(1)
            like.syncSrcParams(self.srcName)
            like.fit(max(verbosity-3, 0))
            lc['normfree'] = dict()
            fitstat = like.optObject.getRetCode()
            if verbosity > 1 and fitstat != 0:
                print "- Fit 2 - Minimizer returned with code: ", fitstat
            lc['normfree']['fitstat'] = fitstat
            lc['normfree']['logL'] = like.logLike.value()
            lc['normfree']['ts'] = like.Ts(self.srcName)
            lc['normfree']['flux_dflux'] = \
                srcnormpar.getValue()/srcnormpar.error()
            if verbosity > 1:
                print '- Fit 2 - log Like:',lc['normfree']['logL'],\
                      '(TS='+str(lc['normfree']['ts'])+')'

            lc['normfree']['nfree']=len(like.freePars(self.srcName))
            lc['normfree']['flux']=like[self.srcName].flux(emin, emax)
            pars = dict()
            for pn in like[self.srcName].funcs['Spectrum'].paramNames:
                p = like[self.srcName].funcs['Spectrum'].getParam(pn)
                pars[p.getName()] = dict(name      = p.getName(),
                                         value     = p.getTrueValue(),
                                         error     = p.error()*p.getScale(),
                                         free      = p.isFree())
            lc['normfree']['pars'] = pars
            ul_type = None
            if ul_bayes_ts != None and lc['normfree']['ts'] < ul_bayes_ts:
                ul_type = 'bayesian'
                [ul_flux, ul_results] = \
                    IntegralUpperLimit.calc_int(like,self.srcName,cl=ul_cl,
                                                skip_global_opt=True,
                                                verbosity = max(verbosity-2,0),
                                                emin=emin, emax=emax,
                                            poi_values = lc['profile']['value'])
            elif ( ul_flux_dflux != None and \
                   lc['normfree']['flux_dflux'] < ul_flux_dflux ) or \
                   ( ul_chi2_ts != None and lc['normfree']['ts'] < ul_chi2_ts):
                ul_type = 'chi2'
                [ul_flux, ul_results] = \
                    IntegralUpperLimit.calc_chi2(like,self.srcName,cl=ul_cl,
                                                 skip_global_opt=True,
                                                 verbosity = max(verbosity-2,0),
                                                 emin=emin, emax=emax)
            if ul_type != None:
                lc['normfree']['ul'] = dict(flux    = ul_flux,
                                            results = ul_results,
                                            type    = ul_type)

            # ----------------------------- FIT 3 -----------------------------

            if verbosity > 1:
                print '- Fit 3 - All parameters of',self.srcName,'free'
            like.setFreeFlag(self.srcName, srcfreepar, 1)
            like.syncSrcParams(self.srcName)
            like.fit(max(verbosity-3, 0))
            lc['allfree'] = dict()
            fitstat = like.optObject.getRetCode()
            if verbosity > 1 and fitstat != 0:
                print "- Fit 3 - Minimizer returned with code: ", fitstat
            lc['allfree']['fitstat'] = fitstat
            lc['allfree']['logL'] = like.logLike.value()
            lc['allfree']['ts'] = like.Ts(self.srcName)
            if verbosity > 1:
                print '- Fit 3 - log Like:',lc['allfree']['logL'],\
                      '(TS='+str(lc['allfree']['ts'])+')'
            lc['allfree']['nfree']=len(like.freePars(self.srcName))
            lc['allfree']['flux']=like[self.srcName].flux(emin, emax)
            pars = dict()
            for pn in like[self.srcName].funcs['Spectrum'].paramNames:
                p = like[self.srcName].funcs['Spectrum'].getParam(pn)
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

        p = scipy.polyfit(profile_x, profile_y, 2);
        prof_max_val = -p[1]/(2*p[0])
        prof_max_logL = p[2]-p[1]*p[1]/(4*p[0])
        if (prof_max_val<min(profile_x)) or (prof_max_val>max(profile_x)):
            print "Warning: corrected minimum %f is outside profile range [%f to %f]" \
                  %(prof_max_val,min(profile_x),max(profile_x))
            print profile_x, profile_y

        profile_fity = scipy.polyval(p,profile_x)
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
            p = scipy.polyfit(profile_x, y, 2);
            dchi2_normfree += 2*(lc['normfree']['logL']
                                 - scipy.polyval(p, prof_max_val))

            if (lc['normfree'].has_key('ul') and
                lc['normfree']['ul']['type'] == 'bayesian'):
                # Arbitrarily assume a quadratic is an OK fit
                y = lc['normfree']['ul']['results']['poi_chi2_equiv']
                p = scipy.polyfit(profile_x, y, 2);
                dchi2_normfree_ul += scipy.polyval(p, prof_max_val)
            else:
                dchi2_normfree_ul += 2*(lc['normfree']['logL']
                                        - scipy.polyval(p, prof_max_val))
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
        file = sys.stdout
        if filename != None:
            file=open(filename,'w')
        if header:
            # print >>file, '%sOptions: %s'%(headstart,' '.join(lc[0]['config']['argv'][1:]))
            chi2 = stats['dchi2_normfree']
            ndof = stats['npar_normfree']
            prob = MyMath.chi2cdfc(chi2,ndof)
            sigma = math.sqrt(MyMath.chi2invc(prob,1))
            print >>file, '%sVariable flux (no UL): chi^2=%.3f (%d DOF) - Pr(>X)=%g (~%g sigma)'%(headstart,chi2,ndof,prob,sigma)
            chi2 = stats['dchi2_normfree_ul']
            ndof = stats['npar_normfree']
            prob = MyMath.chi2cdfc(chi2,ndof)
            sigma = math.sqrt(MyMath.chi2invc(prob,1))
            print >>file, '%sVariable flux (w/UL):  chi^2=%.3f (%d DOF) - Pr(>X)=%g (~%g sigma)'%(headstart,chi2,ndof,prob,sigma)
            chi2 = stats['dchi2_specfree']
            ndof = stats['npar_specfree']
            prob = MyMath.chi2cdfc(chi2,ndof)
            sigma = math.sqrt(MyMath.chi2invc(prob,1))
            print >>file, '%sVariable spectrum:     chi^2=%.3f (%d DOF) - Pr(>X)=%g (~%g sigma)'%(headstart,chi2,ndof,prob,sigma)
            print >>file, '%sProfile minimum: %f (search range: %f to %f)'%(headstart,stats['prof_max_val'],min(stats['prof_x']),max(stats['pro\
f_x']))
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

