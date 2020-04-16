# Overview #


**SiCloneFit** is a Bayesian method for joint inference of tumor clones, clonal genotypes and clonal phylogeny from 
noisy mutation profile of single cells. Given an imperfect noisy genotype matrix from single cells, **SiCloneFit** 
employs a Gibbs sampling algorithm to sample from the posterior distribution of the probabilistic graphical model. 
Each sample contains the cluster indicator for all cells, genotypes of each clone, a clonal phylogeny that represents 
the genealogical relationships and error rate values.

**SiCloneFitIO** is a python wrapper with cmd line tool to covert data matrix for use of siCloneFit from and to pandas
DataFrame and plot imputed matrix overlay with raw data in clustermap.

#### Dependencies and installation for SiCloneFitIO (python package)
1. (Create an virtual environment of choice)
2. `python setup.py install`
3. Install external package `IWSS`
```
git clone https://github.com/BuysDB/IWSS
pip3 install -e ./IWSS
```
4. The binary [SiCloneFitComplete.jar](https://bitbucket.org/hamimzafar/siclonefit/src/master/SiCloneFiTComplete.jar) is directly downloadable.
Download and decompress the SiCloneFiTComplete.jar in the top folder (same level as the README.md).  

#### Dependencies for SiCloneFit (java script)#

* PhyloNet ([https://bioinfo.cs.rice.edu/phylonet]())
* The Apache Commons Mathematics Library ([http://commons.apache.org/proper/commons-math]())
* Parallel Colt ([https://mvnrepository.com/artifact/net.sourceforge.parallelcolt/parallelcolt]())
* Habanero-Java Library ([https://wiki.rice.edu/confluence/display/PARPROG/HJ+Library]())


#### Running siCloneFitIO 
```
siclonefit -j ../hamimzafar-siclonefit/SiCloneFiTComplete.jar -s ../test_data/test1.pickle \
-cn ../test_data/cnv.pickle.gz -o ../test_out/ -n test1 -mm 1 -mp 1
```


Features
--------

* TODO

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
