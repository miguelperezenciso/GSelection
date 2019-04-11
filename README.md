### Python Genomic Selection functions
#### Miguel Perez-Enciso (miguel.perez@uab.es)

### Installation

Download or clone the repository.

### Quick startup
Check the jupyter notebook (```main.ipynb```) and ```main.py``` script

### Input files
Inputs are:
        
- a nsnp x nind genotype file coded as 0, 1, 2 for each genotype.
- a pedigree and phenotypes file, with inds coded as integers 1,2.. nind, 
0 for unknwon parents, followed by phenotype values:

           id   id_father   id_mother   y1   y2 ...

### Coding for missing values
- No missing values are allowed in genotypes matrix.
- For phenotypes, the missing code is -999999. An alternative is to set the vector 
```yIds```. This vector is an integer vector with indices of individuals containing individuals 
without missing phenotypes. For instance, in calling the main function

            gs.doEbv0(criterion, X, y, yIds, h2, ped)

```criterion``` is the evaluation ciriterion ('blup', 'gblup', 'sstep'), ```X``` is the matrix with 
genotypes, ```y``` is a nind x ntrait matrix with ntrait phenotypes, ```yIds``` contains id indices (eg ```yIds = [0, 3, 9]```)
means that inds 1st, 4th and 10th have phenotypes, and ```ped``` is the pedigree.

### Citation
Please cite this if you use or reuse the code:

M. Perez-Enciso, L.C. Ramirez-Ayala,  L.M. Zingaretti. 
SeqBreed: a python tool to evaluate genomic prediction in complex scenarios.
To be submitted.

### How to contribute
Please send comments, suggestions or report bugs to miguel.perez@uab.es. From 2020 on, I will be working at
INIA in Madrid (Spain), check at www.inia.es or in the internet for my new email.