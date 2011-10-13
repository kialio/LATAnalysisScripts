import pyfits
from gt_apps import *
from math import *
from ds9 import *

"ciao"

class quickPlot:

    def __init__(self,base='MySource',irfs='P6_V11_DIFFUSE',model='MySource_model.xml'):

        self.base = base
        self.irfs = irfs
        self.model = model

    def createModel(self,run=True):

        model_map['srcmaps'] = self.base+"_srcMaps.fits"
        model_map['srcmdl'] = self.model
        model_map['outfile'] = self.base+"_modelMap.fits"
        model_map['expcube'] = self.base+"_ltcube.fits"
        model_map['irfs'] = self.irfs
        model_map['bexpmap'] = self.base+"_BinnedExpMap.fits"

        model_map.run()

    def createResidMap(self):

        cmd = "farith "+self.base+"_CMAP.fits "+self.base+"_modelMap.fits "+self.base+"_residMap.fits SUB clobber=yes"
        os.system(cmd)
        print cmd

    def createSigMap(self):

        onImage = pyfits.open(self.base+"_CMAP.fits")
        onData=onImage[0].data.copy()
        onHeader=onImage[0].header
        offImage = pyfits.open(self.base+"_modelMap.fits")
        offData=offImage[0].data.copy()
        sigData=offImage[0].data.copy()

        for x,row in enumerate(sigData):
            for y in enumerate(row):
                sigData[x,y[0]] = (onData[x,y[0]]-offData[x,y[0]])/sqrt(onData[x,y[0]]+offData[x,y[0]])

        newImage = pyfits.PrimaryHDU(sigData)
        newImage.header = onHeader
        newImage.update_header()

        hdulist = pyfits.HDUList([newImage])
        hdulist.writeto(self.base+"_sigMap.fits",clobber=True)

    def createMaps(self):
        
        self.createModel()
        self.createResidMap()
        self.createSigMap()

    def plotMaps(self):

        d = ds9()
        d.set('file '+self.base+'_CMAP.fits')
        d.set('scale log')
        d.set('scale mode minmax')
        d.set('cmap aips0')
        d.set('regions', 'image; text 40 15 # color=black width=3 font="helvetica 10 bold" text={Count map}')
        d.set('tile')
        d.set('frame new')
        d.set('file '+self.base+'_modelMap.fits')
        d.set('cmap a')
        d.set('regions', 'image; text 40 15 # color=black width=3 font="helvetica 10 bold" text={Model map}')
        d.set('frame new')
        d.set('file '+self.base+'_residMap.fits')
        d.set('scale sqrt')
        d.set('scale mode zscale')
        d.set('cmap aips0')
        d.set('regions', 'image; text 40 15 # color=black width=3 font="helvetica 10 bold" text={Residual}')
        d.set('frame new')
        d.set('file '+self.base+'_sigMap.fits')
        d.set('scale mode minmax')
        d.set('regions', 'image; text 40 15 # color=black width=3 font="helvetica 10 bold" text={Significance}')
        xml = open(self.model)
        keywords  = ['RA', 'DEC', 'PointSource']
        keywords2 = ['value', 'name']
        ra = []
        dec = []
        name = []
        columns= []
        for line in xml:
            for word in line.split():
                if keywords[0] in word:
                    for word2 in line.split():
                        if keywords2[0] in word2:
                            columns = word2.split('"')
                            ra.append(columns[1])
                            break
                if keywords[1] in word:
                    for word2 in line.split():
                        if keywords2[0] in word2:
                            columns = word2.split('"')
                            dec.append(columns[1])
                            break
                if keywords[2] in word:
                    for word2 in line.split('"'):
                        columns.append(word2)
                    i = -1
                    for lookname in columns:
                        i += 1
                        if keywords2[1] in lookname:
                            name.append(columns[i+1])
                            break
        d.set('frame 1')              
        i = -1
        for value in ra:
            i += 1
            d.set('regions', 'fk5; point '+ ra[i]+ ' '+ dec[i]+ ' # point=cross 15 color=black width=3 font="helvetica 14 bold" text={'+ name[i]+ '}')
            
        d.set('zoom to fit')
        d.set('match frames wcs')
        d.set('grid yes')
        d.set('grid axes type exterior')
        d.set('grid numlab vertical yes')
        d.set('grid skyformat degrees')
        d.set('grid axes color black')
        d.set('grid tick color black')
        d.set('grid grid color black')
        d.set('grid numlab color black')
        d.set('grid numlab fontsize 14')
