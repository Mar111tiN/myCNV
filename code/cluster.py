import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture as GMM



def get_centers(merge_df, runs=25, comps=3, VAF_limits=(0.05, 0.95), exclude_X=True):
    '''
    use GMM to identify the center cluster and get the means from that
    because GMM occasionally does not identify the center cluster,
    I let the GMM proceed several times and minimize the center cluster
    next, the center cluster can be identified as the maximum center
    '''
    VAFmin, VAFmax = VAF_limits
    # fit the centers to the data 
    if exclude_X:
        merge_df = merge_df.query('Chr != "chrX"')     
    X = merge_df.query('@VAFmin < VAF < @VAFmax and log2ratiomean == log2ratiomean')[['log2ratiomean', 'VAF']]

    gmm = GMM(n_components=comps, covariance_type='diag', n_init=runs).fit(X)
    labels = gmm.predict(X)
    # get the size of the 
    _, counts = np.unique(labels, return_counts=True)
    maxcount = np.max(counts)
    centers = pd.DataFrame(gmm.means_, columns=['log2ratio', 'VAF'])
    # get mean_cov and meanVAF from largest cluster
    meanCov, meanVAF = centers.loc[np.argmax(counts)]
    size = maxcount
            
    print(f'GMM using {runs} inits: center size {size} meanVAF = {round(meanVAF, 2)} meanCov={round(meanCov, 2)}')
    
    return meanCov, meanVAF, centers


def center_data(snp_df, config):
    '''
    retrieve the centers for scaling using GMM
    '''
    
    meanCov, meanVAF, _ = get_centers(snp_df, VAF_limits=config['heteroSNP']['filter']['VAF'])
    # center coverage 
    if config['coverage']['center']:
        print("log2ratio centered around", meanCov)
        snp_df.loc[:, 'log2ratiomean'] = snp_df['log2ratiomean'] - meanCov
    if config['heteroSNP']['center']:
        print("heteroSNP centered around", meanVAF)
        snp_df.loc[:, 'VAF'] = snp_df['VAF'] - meanVAF + 0.5
    return snp_df
