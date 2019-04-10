"""
    Generic Genomic Selection Module
    Reads plink or vcf file and allows computing EBV with different criteria
    BLUP, GBLUP, SSTEP, GWAS or PCA plots

"""

import numpy as np
import pandas as pd

# specific modules
import gselection as gs

# input file dir
ddir='/home/miguel/PycharmProjects/gs/toy/'

# pedigree and phenotype file
pedfile = ddir + '/toy.pedy'

# genotypes file 
xfile = ddir + 'toy.gen'
Transpose = False
ploidy = 2 # set ploidy

# STEP 1: uploads genotypes
X = np.array(pd.read_csv(xfile, header=None, comment='#', sep='\s+'), dtype=float)
nind = X.shape[1]
print('N markers read: ' + str(X.shape[0]))
print('N inds read: ' + str(nind))
print('If you have a nind x nsnp matrix, set Transpose to True')

if Transpose is True: X = X.T

# STEP 2: uploads pedigree and phenotypes
ped = np.array(pd.read_csv(pedfile, header=None, comment='#', sep='\s+'))
y = ped[:,3:]
ped = np.array(ped[:,:3], int)
ntrait = y.shape[1]
print('N inds read: ' + str(X.shape[0]))
print('N traits read: ' + str(ntrait))

# STEP 3: PCA plot
pca = gs.Pca(X)
pca.fit()
pca.plot()
pca.plot(plotFile='pca.pdf') # PCA in pdf file

# STEP 4: GWAS plot for first trait
itrait = 0
gwas = gs.GWAS(X=X)
gwas.fit(y=y[:,itrait])
gwas.plot() # pvalue
gwas.plot(fdr=True) # FDR
gwas.plot(plotFile='gwas.png',fdr=False) # p-values in png file
gwas.plot(plotFile='gwas.pdf',fdr=True) # FDR in pdf file
gwas.print() # prints gwas results

# STEP 5: Obtain inverse NRM
AI = gs.doAInverse(ped)

# STEP 6: Predicting Breeding Values for first trait
itrait = 0
h2 = 0.3    # assumed h2

# BLUP evaluation when last 10 individuals have no phenotypes
# phenotyped individuals (all but last 10)
yids = np.arange(nind-10, dtype=np.int) # contains ids of phenotyped individuals (0 is first indiv)
ebv_blup10 = gs.doEbv0(criterion='blup', X=X, y=y[:,itrait], yIds=yids, h2=h2, ped=ped)

# GBLUP, all individuals phenotyped
ebv_gblup = gs.doEbv0(criterion='gblup', X=X, y=y[:,itrait], h2=h2, nh=ploidy)

# Single Step evaluation assuming only last half of population is genotyped and all inds phenotyped but last 10 ones
# mkrids contains ids of genotyped individuals (0 is first indiv)
mkrids = np.arange(nind//2,nind, dtype=np.int)
Xss = X[:,mkrids]
ebv_sstep05 = gs.doEbv0(criterion='sstep', X=Xss,y=y[:,itrait],  yIds=yids, mkrIds=mkrids, h2=h2, nh=ploidy, ped=ped)


