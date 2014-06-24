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
    ax1.errorbar(E,like._Nobs(),yerr=np.sqrt(like._Nobs()), fmt='o',label='Counts')
    if len(like.sourceNames()) > 16:
        font_size = 6
    else:
        font_size = 14
    ax1.legend(bbox_to_anchor=(1.05, 1.03), loc=2,prop={'size':font_size})
    plt.ylabel(r'Counts [ph s$^{-1}$ cm$^{-2}$]')
    plt.tick_params(axis='x',labelbottom='off')
    
    
    resid = (like._Nobs() - sum_counts)/sum_counts
    resid_err = (np.sqrt(like._Nobs())/sum_counts)

    ax2 = plt.subplot(gs[2,0:2])
    plt.xscale('log')
    ax2.errorbar(E,resid,yerr=resid_err,fmt='o')
    ax2.axhline(0.0,ls=':')
    plt.ylabel('Residuals')
    plt.xlabel('Energy [MeV]')
    
    plt.savefig(filename)
    plt.show()
