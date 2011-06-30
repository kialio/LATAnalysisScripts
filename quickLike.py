import logging
import pyLikelihood
import re
from UnbinnedAnalysis import *
from BinnedAnalysis import *
from UpperLimits import UpperLimits

class quickLike:

    def __init__(self,base='MySource',irfs='P6_V11_DIFFUSE',model='MySource_model.xml',DRMtol = 0.1,binned=False):
        self.base=base
        self.irfs=irfs
        self.model=model
        self.DRMtol=DRMtol
        self.binned=binned
        
        self.ret = re.compile('\n')

        self.logger = logging.getLogger('quickLike')
        self.logger.setLevel(logging.DEBUG)
        fh = logging.FileHandler(self.base+'_quickLike.log')
        fh.setLevel(logging.DEBUG)
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        self.logger.addHandler(fh)
        self.logger.addHandler(ch)
        self.logger.info("Created quickLike object: base="+self.base+\
                             ",irfs="+self.irfs+\
                             ",model="+self.model+\
                             ",DRMtol="+str(self.DRMtol)+\
                             ",binned="+str(self.binned))

    def Print(self):
        self.logger.info("Created quickLike object: base="+self.base+\
                             ",irfs="+self.irfs+\
                             ",model="+self.model+\
                             ",DRMtol="+str(self.DRMtol)+\
                             ",binned="+str(self.binned))
    def makeObs(self):
        if(self.binned):
            self.obs = BinnedObs(srcMaps=self.base+'_srcMaps.fits',
                                 expCube=self.base+'_ltcube.fits',
                                 binnedExpMap=self.base+'_BinnedExpMap.fits',
                                 irfs=self.irfs)
        else:
            self.obs = UnbinnedObs(self.base+'_filtered_gti.fits',
                                   self.base+'_SC.fits',
                                   expMap=self.base+'_expMap.fits',
                                   expCube=self.base+'_ltCube.fits',
                                   irfs=self.irfs)
        self.logger.info(self.ret.subn(', ',str(self.obs))[0])

    def initDRM(self):
        if(self.binned):
            self.DRM = BinnedAnalysis(self.obs,self.model,optimizer="DRMNGB")
        else:
            self.DRM = UnbinnedAnalysis(self.obs,self.model,optimizer="DRMNGB")
        self.DRM.tol = self.DRMtol
        self.logger.info(self.ret.subn(', ',str(self.DRM))[0])

    def initAltFit(self,opt="MINUIT"):
        if(self.binned):
            self.ALTFIT = BinnedAnalysis(self.obs,self.model,optimizer=opt)
        else:
            self.ALTFIT = UnbinnedAnalysis(self.obs,self.model,optimizer=opt)
        self.ALTFIT.tol = self.DRMtol
        self.logger.info(self.ret.subn(', ',str(self.ALTFIT))[0])
                
    def fitDRM(self):

        altfit=False
        try:
            self.DRM.fit(verbosity=3)
        except:
            self.logger.error("Initial DRM Fit Failed")
            try:
                self.logger.info("Trying looser tolerance (DRMtol*10.)")
                self.DRM.tol = self.DRMtol * 10.
                self.DRM.fit(verbosity= 3)
            except:
                self.logger.error("Second DRM Fit Failed")
                try:
                    self.logger.info("Trying tighter tolerance (DRMtol*0.1)")
                    self.DRM.tol = self.DRMtol * 0.1
                    self.DRM.fit(verbosity= 3)
                except:
                    self.logger.error("Third DRM Fit Failed")
                    try:
                        self.logger.info("Trying alternate fit algorithm (MINUIT)")
                        self.initAltFit()
                        self.ALTFIT.fit(verbosity=3)
                        altfit = True
                    except:
                        self.logger.error("Alternative fit algorithm failed, bailing")
                        return

        if(altfit):
            self.logger.info("ALTFIT Fit Finished.  Total TS: "+str(self.ALTFIT.logLike.value()))
            self.ALTFIT.logLike.writeXml(self.base+'_likeDRM.xml')
            self.logger.info("Saved ALTFIT as "+self.base+"_likeDRM.xml")
        else:
            self.DRM.logLike.writeXml(self.base+'_likeDRM.xml')
            self.logger.info("DRM Fit Finished.  Total TS: "+str(self.DRM.logLike.value()))
            self.logger.info("Saved DRM as "+self.base+"_likeDRM.xml")

    def initMIN(self,tol=1e-4):
        if(self.binned):
            self.MIN = BinnedAnalysis(self.obs,self.base+'_likeDRM.xml',optimizer='NewMinuit')
        else:
            self.MIN = UnbinnedAnalysis(self.obs,self.base+'_likeDRM.xml',optimizer='NewMinuit')
        self.MIN.tol = tol
        self.MINobj = pyLike.NewMinuit(self.MIN.logLike)
        self.logger.info(self.ret.subn(', ',str(self.MIN))[0])

    def fitMIN(self):
        self.MIN.fit(covar=True, optObject=self.MINobj,verbosity=4)
        self.MIN.logLike.writeXml(self.base+'_likeMinuit.xml')
        self.logger.info("NEWMINUIT Fit Finished.  Total TS: "+str(self.MIN.logLike.value()))
        self.logger.info("NEWMINUIT Fit Status: "+str(self.MINobj.getRetCode()))
        self.logger.info("NEWMINUIT fit Distance: "+str(self.MINobj.getDistance()))
        if(self.MINobj.getRetCode() > 0):
            self.logger.error("NEWMINUIT DID NOT CONVERGE!!!")
            self.logger.error("The fit failed the following tests: "+self.decodeRetCode(self.MINobj.getRetCode()))

    def printSource(self,source,Emin=100,Emax=300000):
        print "TS: ",self.MIN.Ts(source)
        print "Npred: ",self.MIN.NpredValue(source)
        print "Flux: ",self.MIN.flux(source,emin=Emin,emax=Emax)
        print "Flux Error: ",self.MIN.fluxError(source,emin=Emin,emax=Emax)
        print "Index: ",self.MIN.model[source].funcs['Spectrum'].getParam('Index').value()
        print "Index Error: ",self.MIN.model[source].funcs['Spectrum'].getParam('Index').error()
        self.logger.info(source+" TS: "+str(self.MIN.Ts(source))+\
                             " Npred: "+str(self.MIN.NpredValue(source))+\
                             " Flux: "+str(self.MIN.flux(source,emin=Emin,emax=Emax))+\
                             " Flux Error: "+str(self.MIN.fluxError(source,emin=Emin,emax=Emax))+\
                             " Index: "+str(self.MIN.model[source].funcs['Spectrum'].getParam('Index').value())+\
                             " Index Error: "+str(self.MIN.model[source].funcs['Spectrum'].getParam('Index').error()))
        
    def calcUpper(self,source,Emin=100,Emax=300000):
        self.ul = UpperLimits(self.MIN)
        self.ul[source].compute(emin=Emin,emax=Emax)
        print self.ul[source].results
        self.logger.info(source+" UL: "+str(self.ul[source].results[0]))

    def removeWeak(self,mySource,tslimit=0,DistLimit=0,RemoveFree=False,RemoveFixed=False):
        print "Name, TS Value, Frozen?, Distance"
        for name in self.MIN.sourceNames():
            remove = False
            distance = 0
            sourceTS = self.MIN.Ts(name)
            if(self.MIN.model[name].src.getType() == 'Point'):
                distance = self.MIN._separation(self.MIN.model[mySource].src,self.MIN.model[name].src)
            if(self.MIN.freePars(name).size() > 0):
                indexFree = "Free"
                if( (sourceTS < tslimit) and (distance > DistLimit) and RemoveFree ):
                    remove = True
            else:
                indexFree = "Fixed"
                if( (sourceTS < tslimit) and (distance > DistLimit) and RemoveFixed ):
                    remove = True
            if( remove ):
                self.logger.info("Removing "+name+", TS: "+str(sourceTS)+", Frozen?: "+str(indexFree)+", Distance: "+str(distance))
                self.MIN.deleteSource(name)
            else:
                self.logger.info("Retaining "+name+", TS: "+str(sourceTS)+", Frozen?: "+str(indexFree)+", Distance: "+str(distance))
                    

    def paramsAtLimit(self, limit = 0.01):
        for src in self.MIN.sourceNames():
            for name in self.MIN.model[src].funcs['Spectrum'].paramNames:
                bounds = self.MIN.model[src].funcs['Spectrum'].getParam(name).getBounds()
                value  = self.MIN.model[src].funcs['Spectrum'].getParam(name).value()
                if( abs(bounds[0] - value) < limit ):
                    self.logger.error("The "+name+" ("+str(value)+") of "+src+" is close to its lower limit ("+str(bounds[0])+")")
                if( abs(bounds[1] - value) < limit):
                    self.logger.error("The "+name+" ("+str(value)+") of "+src+" is close to its lower limit ("+str(bounds[1])+")")



    def decodeRetCode(self, retCode):

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

    def delObs(self):
        del self.obs
