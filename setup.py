#!/usr/bin/evn python

import sys
import os
from distutils.core import setup

try:
    import ds9
except ImportError:
    print "pyds9 not found (required for making plots with quickPlot)."
    sys.exit(1)
else:
    print "pyds9 found."


if os.environ.get("FERMI_DIR"):
    print "The Fermi Science tools seem to be set up."
else:
    print "The Fermi Science tools are not set up."
    sys.exit()

setup(name='LATAnalysisScripts',
      version='0.1.0',
      description='Fermi LAT Analysis Scripts (quickScripts)',
      author='Jeremy S. Perkins and Davide Donato (FSSC)',
      author_email='fermihelp@milkyway.gsfc.nasa.gov',
      url='http://fermi.gsfc.nasa.gov/ssc/',
      py_modules=['quickUtils',
                  'quickAnalysis',
                  'quickLike',
                  'quickPlot'],
      )



