""" This module incorporates several methods to seamlessly implement selection
    in seqbreed software
    Main arguments
        - criterion: Random, Phenotype, BLUP, GBLUP, Sstep, ...
        - generation: discrete(nsize) / continuous
        - intensity: male, female
        - mating: assortative / random
    The main output is an extended pedigree with offspring of selected parents

    ok: gwas option, fdr option
    TODO:
    - optimize reading /writing snp data
    - Add qtn positions in gwas plot
    - add labels pca,
    - add xvalidation = cross_val_predict(lr, boston.data, y, cv=10)
    - ULL: check pedigree order in complex settings (BLUP / GSSTEP)
"""
import numpy as np
import sys
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import fdrcorrection as fdr


# missing code, ULL: it must be float
MISSING = -999999.

#########
class Pca:
#########
    """ plots PCA of selected inds and markers """
    def __init__(self, X):
        self.X = X
        self.p = []

    def fit(self):
        ''' performs PCA on pop and chip snps
            scikit assumes nsamples x nvariables
        '''
        pca = PCA(n_components=2)
        self.p = pca.fit(self.X.T).fit_transform(self.X.T)

    def plot(self, plotFile=None, labels=None):
        plt.title('PCA decomposition')
        if labels is None:
            plt.scatter(self.p[:, 0], self.p[:, 1])
        else:
            for i in range(max(labels)+1):
                plt.scatter(self.p[labels==i, 0], self.p[labels==i, 1], label='Gen '+str(i))
            plt.legend(loc='best')
        if plotFile is None:
            plt.show()
        else:
            plt.savefig(plotFile)
        plt.close('all')

##########
class GWAS:
##########
    """ Implements GWAS
        - X: the genotypes used to perform the gwas
        - b(float array): regression coefficients
        - pvalue(float array)
        - qvalue(float array)
    """
    def __init__(self, X):
        self.X = X
        self.b = []
        self.se = []
        self.pvalue = []
        self.qvalue = []

    def fit(self, y):
        """ carry out GWAS
        """
        nind = len(y)
        # output
        out = np.array([])
        # genotypes provided in X
        for i in range(self.X.shape[0]):
              # trick to sum genotypes over haplotypes
              x = np.split(self.X[i, :], nind)
              x = np.asarray(list(map(sum, x)))
              b, intercept, r_value, p_value, std_err = stats.linregress(x, y)
              out = np.append(out, [b, std_err, p_value])

        out = out.reshape(len(out) // 3, 3)
        self.b, self.se, self.pvalue = out[:,0], out[:,1], out[:,2]
        # FDR obtained from statsmodels package
        self.fdr = fdr(self.pvalue)[1]

    def plot(self, plotFile=None, fdr=None):
        """ GWAS plot """
        # p-values are printed
        if fdr is None:
            plt.scatter(list(range(len(self.pvalue))), -np.log10(self.pvalue), marker='o')
            plt.ylabel('-log10 P-value')
        # FDR values are printed
        else:
            plt.scatter(list(range(len(self.fdr))), -np.log10(self.fdr), marker='o')
            plt.ylabel('-log10 FDR')

        plt.title('GWAS')
        plt.xlabel('SNP')
        if plotFile is not None:
            plt.savefig(plotFile)
        else:
            plt.show()
        plt.close('all')

    def print(self, gwasFile=None):
        """ prints to file or STDOUT
        """
        f = sys.stdout if gwasFile == None else open(qtnFile, 'w')
        f.write('#SNP COEFF SE PVALUE FDR' + '\n')
        for isnp in range(len(self.b)):
            line = str(isnp+1) + ' ' + str(self.b[isnp]) + ' ' + str(self.se[isnp]) \
                   + ' ' + str(self.pvalue[isnp]) + ' ' + str(self.fdr[isnp])
            f.write(line + '\n')
        pass



##########
def doGRM(X, nh=2):
##########
    """ Computes GRM
        - X contains genotypes
        - nh: ploidy
    """
    p = X.mean(axis=1)/nh
    c = sum(nh * p * (1.-p))
    s = StandardScaler(with_std=False)
    X = s.fit_transform(X.T).T
    G = np.matmul(X.T, X) / c
    return G


##############
def doAInverse(ped):
##############
    """
    Returns A-inverse using Henderson's rules, w/o considering inbreeding
    ped is a nind x 3 np.array with id, father and mother ids, 0 for unknown parent
    """
    w = np.array([1., -0.5, -0.5])
    res = np.array([2., 4./3., 1.])
    nind = ped.shape[0]

    # a dummy 0 row and 0 column is created to facilitate addressing posns
    AI = np.zeros(shape=((nind+1),nind+1))

    c = np.zeros(nind, dtype=int)
    c[ped[:,1]==0] += 1
    c[ped[:,2]==0] += 1

    # ULL with id indices
    for i in range(nind):
        for k1 in range(3):
            for k2 in range(3):
                AI[ped[i,k1],ped[i,k2]] += w[k1] * w[k2] * res[c[i]]

    # rm row 0 and col 0
    return AI[1:,1:]

###########
def dogblup(h2, y, grmFile=None, G=None, invert=True):
###########
    """ (G)BLUP evaluation, returns EBVs
        G can be passed as argument or printed in grmFile
        G is assumed to be the inverse if Invert=False
        - h2(float): h2 used [req]
        - y(array float): phenotypes, can contain missing values as coded by MISSING [req]
        - grmFile(str): file containing Cov matrix of breeding values or its inverse [None]
        - G(nind x nind float array): np.array with Cov matrix of breeding values or its inverse [None]
        - invert(bool): True if G should be inverted or False if already inverted [True]
    """
    # G is read from file
    if grmFile!=None:
        G = np.array(pd.read_csv(grmFile, header=None, comment='#', sep='\s+'))

    nind = len(y)
    if nind != G.shape[0]:
        sys.exit('GRM matrix size must be equal to no. inds')

    # builds MME
    if invert: 
        np.fill_diagonal(G, np.diag(G) * 1.05)
        lhs = np.linalg.inv(G)

    else:
        lhs = G

    y[y==MISSING] = 0
    x1 = np.repeat(1.,nind)
    x1[y==0] = 0
    x2 = np.append(x1, len(x1[x1!=0])) # last element is no. of data

    lhs = lhs * ((1.-h2)/h2)
    np.fill_diagonal(lhs, np.diag(lhs)+x1)
    lhs = np.vstack([lhs, x1])  # add row
    lhs = np.hstack((lhs, x2.reshape(nind+1,1))) # add col
    rhs = np.append(y, np.sum(y))

    # solves, all els but last are the ebvs
    ebv = np.linalg.solve(lhs, rhs)[:-1]
    return ebv


###########
def dosstep(h2, y, im, AInv, A, G):
###########
    """ Single step evaluation
        - h2(float): h2 used [req]
        - y(array float): phenotypes, can contain missing values as coded by MISSING [req]
        - im(int array): set of genotyped individuals (indexed starting with 0) [req]
        - ainvFile(str): file name that contains A inverse [req]
        - aFile(str): file name that contains A [req]
        - grmFile(str): file name that contains GRM of genotyped individuals [req]
    """
    nind = len(y)
    if nind != A.shape[0]: sys.exit('NRM matrix size must be equal to # inds')

    ngind = len(im)
    if ngind != G.shape[0]: sys.exit('GRM matrix size must be equal to # genotyped inds')

    # just in case npd
    np.fill_diagonal(G, np.diag(G) * 1.05)

    # builds H-1
    Hinv = AInv
    Ginv = np.linalg.inv(G)
    A22inv = np.linalg.inv(A[im,:][:,im])
    Z = Ginv - A22inv
    j=0
    for i in im:
        Hinv[i,im] += Z[j,:]
        j+=1

    ebv = dogblup(h2, y, G=Hinv, invert=False)
    return ebv

#########
def doEbv0(criterion, X, y, h2, mkrIds=None, yIds=[], nh=2, ped=None):
#########
    """ returns estimated breeding values (EBVs), assume evaluation is on first trait (itrait=0)
    - y
    - criterion(str): evaluation method: 'phenotype', 'blup', 'gblup', 'sstep' ['random']
    - h2(float): used h2
    - X: contains genotypes
    - mkrIds(int array): set of genotyped individuals (indexed starting with 0) [None, req for sstep]
    - yIds(int array): integer array specifying individuals with data (indexed starting with 0) [all]
    """
    criterion = criterion.lower()
    nind = len(y)

    # remove missing phenotypes ( none by default)
    if len(yIds)==0:
        yIds = np.arange(nind)
    y0 = y
    y = np.repeat(MISSING, nind)
    y[yIds] = y0[yIds]

    # computes A inverse & A required for blup or sstep
    if criterion == 'blup' or criterion=='sstep':
        if criterion == 'blup':
            AI = doAInverse(ped) #--> computes A-inv
        elif criterion == 'sstep':
            AI = doAInverse(ped)    # --> computes A-inv
            A =  np.linalg.inv(AI)  #--> computes A

    # computes GRM if undefined, required for glbup or sstep
    if (criterion == 'gblup' or criterion=='sstep'):
        G = doGRM(X, nh)

    # ---> Computes EBVs
    # pure drift
    if criterion == 'blup':
        ebv = dogblup(h2, y, G=AI, invert=False)
    # GBLUP
    elif criterion == 'gblup':
        ebv = dogblup(h2, y, G=G, grmFile=None)
    # single step
    elif criterion == 'sstep':
        ebv = dosstep(h2, y, mkrIds, AInv=AI, A=A, G=G)
    # unknown
    else:
        sys.exit('Unknown selection criterion ' + criterion)

    return ebv
