# Overview #

**SiCloneFit** is a Bayesian method for joint inference of tumor clones, clonal genotypes and clonal phylogeny from noisy mutation profile of single cells. Given an imperfect noisy genotype matrix from single cells, **SiCloneFit** employs a Gibbs sampling algorithm to sample from the posterior distribution of the probabilistic graphical model. Each sample contains the cluster indicator for all cells, genotypes of each clone, a clonal phylogeny that represents the genealogical relationships and error rate values.

# Dependencies #

* PhyloNet ([https://bioinfo.cs.rice.edu/phylonet]())
* The Apache Commons Mathematics Library ([http://commons.apache.org/proper/commons-math]())
* Parallel Colt ([https://mvnrepository.com/artifact/net.sourceforge.parallelcolt/parallelcolt]())
* Habanero-Java Library ([https://wiki.rice.edu/confluence/display/PARPROG/HJ+Library]())

# Installation #

The binary [SiCloneFitComplete.jar](https://bitbucket.org/hamimzafar/siclonefit/src/master/SiCloneFiTComplete.jar) is directly downloadable.

# Singlet Model #
The singlet model of **SiCloneFit** can be run as follows
```
java -jar SiCloneFitComplete.jar [-m][-n][-fp][-fn][-r][-df][-ipMat][-missing]
[-burnin][-iter][-printIter][-treeIter][-cellNames][-trueTree][-outDir]
```
## Input Data Arguments ##
* ```-m <Integer>``` Replace <Integer> with the number of cells in the dataset.

* ```-n <Integer>``` Replace <Integer> with the number of mutations (rows) in the dataset.

* ```-ipMat <filename>``` Replace <filename> with the path to the file containing the genotype matrix.

* ```-fp <Double>``` Set <Double> to the estimated false positive rate of the single-cell sequencing experiment. 

* ```-fn <Double>``` Set <Double> to the estimated allelic dropout rate of the single-cell sequencing experiment. This will be used for setting the prior distribution for estimating false negative rate.

* ```-df <Integer>``` Set <Integer> to 1 if the input matrix is ternary. Set <Integer> to 0 if the input matrix is binary.

* ```-missing <Double>``` Set <Double> to the fraction of missing data in the input genotype matrix.

## Gibbs Sampling Setting Arguments ##
* ```-r <Integer>``` Set <Integer> to the desired number of Markov chains to use.

* ```-burnin <Integer>``` Set <Integer> to the number of iterations for burnin of each Markov chain.

* ```-iter <Integer>``` Set <Integer> to the number of iterations to run for each Markov chain after burnin.

* ``` -printIter <Integer>``` Set <Integer> to the number of iterations after which the likelihood is to be printed to the standard output. This should be smaller than the total number of iterations for each Markov chain. 

* ``` -treeIter <Integer>``` Set <Integer> to the number of iterations to use for the Metropolis-Hastings sampler for sampling the clonal phylogeny and evolution model parameters.

## Other Arguments ##
* ```-cellNames <filename>``` Replace <filename> with the path to the file containing the names of the cells. The cell names should be written in a single row separated by blank spaces. The order of the names of the cells should be same as the order in any row of the file containing the genotype matrix. If no such file is provided, the cells are numbered 1 to m.

* ```-trueTree <filename>``` Replace <filename> with the path to the file containing the ground truth tree. The true phylogenetic tree should be written in Newick format. The names of the leaves in the true tree should match the names of the cells. If this file is provided, then the program outputs tree reconstruction error.

* ```-outDir <dir>``` Replace <dir> with the path to the directory where the posterior samples will be stored.

# Doublet Model #
The doublet model of **SiCloneFit** can be run as follows
```
java -cp SiCloneFitComplete.jar siCloneFiT.main.RunSiCloneFiTDoublet [-m][-n][-fp][-fn]
[-r][-df][-ipMat][-missing][-doublet][-burnin][-iter][-printIter]
[-treeIter][-cellNames][-trueTree][-outDir]
```
## Arguments ##
```RunSiCloneFiTDoublet``` program has all the arguments as used by the default singlet model. In addition, it has the following argument

* ```-doublet <Double>``` Set <Double> to the doublet rate of the single cell sequencing experiment.

