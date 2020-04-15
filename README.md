# Overview #

**siCloneFitIO** is a python wrapper for running an single cell clonal phylogeny inference algorithm siCloneFit written in java, and the visualization of the results.

# Dependencies #
python3
matplotlib

# Installation #

The SiCloneFit is downloadable from https://bitbucket.org/hamimzafar/siclonefit/src/master/
Download and decompress the SiCloneFiTComplete.jar in the top folder (same level as the README.md).  

# Running siCloneFitIO #


**SiCloneFit** is a Bayesian method for joint inference of tumor clones, clonal genotypes and clonal phylogeny from noisy mutation profile of single cells. Given an imperfect noisy genotype matrix from single cells, **SiCloneFit** employs a Gibbs sampling algorithm to sample from the posterior distribution of the probabilistic graphical model. Each sample contains the cluster indicator for all cells, genotypes of each clone, a clonal phylogeny that represents the genealogical relationships and error rate values.

# Dependencies #

* PhyloNet ([https://bioinfo.cs.rice.edu/phylonet]())
* The Apache Commons Mathematics Library ([http://commons.apache.org/proper/commons-math]())
* Parallel Colt ([https://mvnrepository.com/artifact/net.sourceforge.parallelcolt/parallelcolt]())
* Habanero-Java Library ([https://wiki.rice.edu/confluence/display/PARPROG/HJ+Library]())

# Installation #

The binary [SiCloneFitComplete.jar](https://bitbucket.org/hamimzafar/siclonefit/src/master/SiCloneFiTComplete.jar) is directly downloadable.




python wrapper with cmd line tool to covert data matrix for use of siCloneFit and pandas



Features
--------

* TODO

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
