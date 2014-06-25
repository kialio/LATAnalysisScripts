#!/usr/bin/env python

defaultCommonConfig = {	"base" : 'MySource',
                        "eventclass" : 2,
                        "binned" : False,
                        "irfs" : "P7SOURCE_V6",
                        "verbosity" : 4,
                        "multicore" : 0}

defaultlikelihoodConfig = {"model" : "MySource_model.xml",
       	                   "sourcename" : "Source Name",
                           "drmtol" : 0.1,
                           "mintol" : 1e-4}