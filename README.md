# Overview #

**SiCloneFitIO** is a python wrapper with cmd line tool to covert data matrix for use of siCloneFit from and to pandas
DataFrame and plot imputed matrix overlay with raw data in clustermap.

**SiCloneFit** is a Bayesian method for joint inference of tumor clones, clonal genotypes and clonal phylogeny from
noisy mutation profile of single cells. Given an imperfect noisy genotype matrix from single cells, **SiCloneFit**
employs a Gibbs sampling algorithm to sample from the posterior distribution of the probabilistic graphical model.
Each sample contains the cluster indicator for all cells, genotypes of each clone, a clonal phylogeny that represents
the genealogical relationships and error rate values.


#### Getting started with SiCloneFitIO
1. (Create an virtual environment of choice)
2. Download and install siclonefitio
```
git clone https://github.com/zztin/siCloneFitIO.git
cd ./siclonefit
python setup.py install
```
3. Run tests for siclonefitio  
```
cd tests/

* test007 test the complete functionality, takes around 5 mins
python -m unittest test_siclonefitio.TestSiclonefitio.test_007_siclonefitIO

* If failed, run all tests. This takes around 10 mins.
python -m unittest

```

4. install other dependencies for siclonefit if needed (check https://bitbucket.org/hamimzafar/siclonefit/src/master/)

#### Running siCloneFitIO

```
siclonefit -j ../hamimzafar-siclonefit/SiCloneFiTComplete.jar -s {input_snvmatrix_pd_dataframe_pickle} \
-cn ../test_data/cnv.pickle.gz -o {output_dir} -n test1 -mm {minMeasurement} -mp {minPresence}
```

```
# run test data. the output test_out_local/ should be identical to test_out/ on github repo.

siclonefit -j ../hamimzafar-siclonefit/SiCloneFiTComplete.jar -s ../test_data/test1.pickle \
-cn ../test_data/cnv.pickle.gz -o ../test_out_local/ -n test1 -mm 1 -mp 1
```



#### Credits

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
