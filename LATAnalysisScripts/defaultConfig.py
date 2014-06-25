#!/usr/bin/env python

defaultCommonConfig = {	"base" : 'MySource',
                        "eventclass" : 2,
                        "binned" : False,
                        "irfs" : "P7SOURCE_V6",
                        "verbosity" : 4,
                        "multicore" : 0}

defaultLikelihoodConfig = {"model" : "MySource_model.xml",
       	                   "sourcename" : "Source Name",
                           "drmtol" : 0.1,
                           "mintol" : 1e-4}

defaultAnalysisConfig = {"ra" : 0,
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
                         "nxpix" : -1,
                         "nypix" : -1,
                         "filter" : "DATA_QUAL==1 && LAT_CONFIG==1",
                         "roicut" : "yes"}

defaultCurveConfig = {'tstart'  : 0,
                      'tstop'   : 0,
                      'tstep'   : 86400,
                      'tsmin'   : 1,
                      'model'   : 'MySource_model_lc.xml',
                      'summary' : 'lc_summary.dat',
                      'output'  : 'lc.pickle',
                      'ulfluxdf' : 2.0,
                      'ulbayes'  : 4.0,
                      'ulchi2'  : 4,
                      'ulcl'    : 0.95,
                      'opt'     : 'MINUIT',
                      'sliding' : False,
                      'rebin'   : 1}