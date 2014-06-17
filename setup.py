#!/usr/bin/evn python

import sys
import os
from distutils.core import setup
import LATAnalysisScripts as LAS

#try:
#    import ds9
#except ImportError:
#    print "pyds9 not found (required for making plots with quickPlot)."
#    sys.exit(1)
#else:
#    print "pyds9 found."

#try:
#    import make2FGLxml
#except ImportError:
#    print "make2FGLxml not found (required for making an XML model)."
#    sys.exit(1)
#else:
#    print "make2FGLxml found."

fermi_dir = os.environ.get("FERMI_DIR")

if fermi_dir:
    print "The Fermi Science tools seem to be set up."
else:
    print "The Fermi Science tools are not set up."
    sys.exit()

setup(name='LATAnalysisScripts',
      version=LAS.__version__,
      description='Fermi LAT Analysis Scripts (quickScripts)',
      author=LAS.__author__,
      author_email='fermihelp@milkyway.gsfc.nasa.gov',
      url='http://fermi.gsfc.nasa.gov/ssc/',
      packages=['LATAnalysisScripts'],
      py_modules=['quickUtils',
                  'quickAnalysis',
                  'quickPlot',
                  'quickCurve',
                  'make2FGLxml'],
      data_files=[(fermi_dir+"/bin",['scripts/quickAnalysis',
                                     'scripts/quickLike',
                                     'scripts/quickPlot',
				      'scripts/quickCurve'])],
      )



