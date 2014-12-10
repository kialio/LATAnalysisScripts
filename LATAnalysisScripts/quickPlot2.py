#!/usr/bin/env python

def Plot2DModel(like,filename='2DModel.png'):

    '''This function plots out the model and data counts
    along with the models of each source.  It takes a 
    likelihood object as input.  This is similar to the 
    plot function when running the balistic gtlike tool.'''

    #import matplotlib as mpl
    #mpl.use('Agg')

    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    import numpy as np

    E = (like.energies[:-1] + like.energies[1:])/2.
    plt.figure(figsize=(10,5))
    
    gs = gridspec.GridSpec(3,3)
    gs.update(left=0.1,right=0.9,bottom=0.15, hspace=0.00, wspace=0.05)
    
    ax1 = plt.subplot(gs[0:2,0:2])    
    
    #plt.ylim((0.4,1e4))
    plt.xlim((like.energies[0],like.energies[-1]))
    sum_counts = np.zeros_like(like._srcCnts(like.sourceNames()[0]))
    for sourceName in like.sourceNames():
        sum_counts = sum_counts + like._srcCnts(sourceName)
        if sourceName[0] == '_':
            sN = sourceName.replace('_','',1)
        else:
            sN = sourceName
        ax1.loglog(E,like._srcCnts(sourceName),label='{}'.format(sN))
    plt.ylim(ymax = 2*np.max(sum_counts))
    ax1.loglog(E,sum_counts,label='Total Model')
    #ax1.errorbar(E,like._Nobs(),yerr=np.sqrt(like._Nobs()), fmt='o',label='Counts')
    ax1.errorbar(E,like.nobs,yerr=np.sqrt(like.nobs), fmt='o',label='Counts')
    if len(like.sourceNames()) > 16:
        font_size = 6
    else:
        font_size = 14
    ax1.legend(bbox_to_anchor=(1.05, 1.03), loc=2,prop={'size':font_size})
    plt.ylabel(r'Counts [ph s$^{-1}$ cm$^{-2}$]')
    plt.tick_params(axis='x',labelbottom='off')
    
    
    #resid = (like._Nobs() - sum_counts)/sum_counts
    #resid_err = (np.sqrt(like._Nobs())/sum_counts)
    resid = (like.nobs - sum_counts)/sum_counts
    resid_err = (np.sqrt(like.nobs)/sum_counts)

    ax2 = plt.subplot(gs[2,0:2])
    plt.xscale('log')
    ax2.errorbar(E,resid,yerr=resid_err,fmt='o')
    ax2.axhline(0.0,ls=':')
    plt.ylabel('Residuals')
    plt.xlabel('Energy [MeV]')
    
    plt.savefig(filename)
    plt.show()

def PSF(energy):

    from scipy import interpolate

    psf_data = np.array([[29.378069628779578, 11.974138435550062],
        [53.123991394359386, 7.957256520655659],
        [93.21305548209773, 5.218638956603347],
        [168.5559881858929, 3.333502332614621],
        [295.7537313410467, 2.073929041512339],
        [524.1755314484437, 1.3074112382033989],
        [947.8599776522383, 0.8027497013599825],
        [1679.9281120535993, 0.5127712463425323],
        [2977.4001732388747, 0.3452772300468028],
        [5383.995040496365, 0.2549696862936114],
        [9542.258178354516, 0.19331232401691525],
        [16575.885731063016, 0.15655123158441375],
        [29973.964322973727, 0.13364545640595094],
        [53660.061093243254, 0.1234818289581322],
        [95103.75714811166, 0.12512043461577865],
        [166872.09599510877, 0.13016787553232645],
        [304797.6595808452, 0.13903680232834606]])

    f = interpolate.interp1d(psf_data[:,0],psf_data[:,1])

    return f(energy)

def PlotSigMaps(quickLikeObj,filename="SigMaps.png"):

    import matplotlib.pyplot as plt
    import numpy as np

    fig, axs = plt.subplots(nrows=15, ncols=4,figsize=(20,60))
    for idx,ax in enumerate(axs):
        for offset in [0,1]:    
            Ebin = idx*2 + offset
            if hasattr(quickLikeObj, 'energyBins'):
                Emax = quickLikeObj.energyBins['E_MAX'][Ebin]
                Emin = quickLikeObj.energyBins['E_MIN'][Ebin]
                psf = quickLikeObj.psf(Emin/1000.)    
            else:
                Emax = quickLikeObj.energies[1:][Ebin]
                Emin = quickLikeObj.energies[:-1][Ebin]
                psf = PSF(Emin/1000.)

            if psf > 0.5:
                psf = 0.5
            deg_per_pix = 0.1
    
            column = offset*2
        
            im = ax[column].imshow(quickLikeObj.sigmaMaps[Ebin],vmin=-5,vmax=5,interpolation='none')    
            ax[column].contour(quickLikeObj.sigmaMaps[Ebin],[5],colors=('w'))
            ax[column].set_title('E = {:.2f} - {:.2f} MeV'.format(Emin/1000.,Emax/1000.))
            h = ax[column+1].hist(quickLikeObj.sigmaMaps[Ebin].flatten(),bins=np.linspace(-5,5,30))
            ax[column+1].set_title('Smoothing = {:.2f} Deg.'.format(float(psf)))
            ax[column+1].set_xlim([-6,6])
            x0, x1 = ax[column+1].get_xlim()
            y0, y1 = ax[column+1].get_ylim()
            ax[column+1].set_aspect((x1-x0)/(y1-y0))
    
    cax = fig.add_axes([0.1, 0.05, 0.8, 0.01])
    fig.colorbar(im, cax=cax, orientation='horizontal')
    plt.savefig(filename)
    plt.show()
