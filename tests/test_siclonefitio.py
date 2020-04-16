#!/usr/bin/env python

"""Tests for `siclonefitio` package."""


import unittest
import os
import shutil
import pandas as pd
import shlex, subprocess

from siclonefitio import formatting
import visual.plot_imp_matrix as pm


class TestSiclonefitio(unittest.TestCase):
    """Tests for `siclonefitio` package."""

    @classmethod
    def setUpClass(cls):
        print("setUpClass")

    @classmethod
    def tearDownClass(cls):
        print("tearDownClass")
        shutil.rmtree("../unittest_000_002")
        shutil.rmtree("../unittest_005")
        shutil.rmtree("../unittest_006")

    def setUp(self):
        """Set up test fixtures, if any."""
        print("setUP")


    def tearDown(self):
        """Tear down test fixtures, if any."""
        print("tearDown")


    def test_000_makedir(self):
        """Test something."""
        formatting.make_dir("../unittest_000_002/")
        assert os.path.exists("../unittest_000_002/")

    def test_001_siclonefit_path(self):
        """Test if siclonefit binary is properly installed in designed path"""
        assert os.path.exists("../hamimzafar-siclonefit/SiCloneFiTComplete.jar"), "The siclonefit binary should locate at siclonefitio/hamimzafar-siclonefit/SiCloneFiTComplete.jar"

    def test_002_formatting_to_sifit(self):
        """
        Test formatting.sifit_formatting create correct input for siclonefit & exert correct command line for running siclonefi
        """
        cmd = "java -jar ../hamimzafar-siclonefit/SiCloneFiTComplete.jar -m 27 -n 6 -fp 0.001 -fn 0.0001 " \
              "-r 1 -df 0 -ipMat ../test_out//test1_mp1_mm1_siclonefit_input.txt -missing 0.24691 " \
              "-cellNames ../test_out//test1_mp1_mm1_cellNames.txt " \
              "-geneNames ../test_out//test1_mp1_mm1_geneNames.txt " \
              "-outDir ../test_out/"
        out = formatting.sifit_formatting("../test_data/test1.pickle",
                                          "../hamimzafar-siclonefit/SiCloneFiTComplete.jar",
                                          "../test_out/",
                                          minPresence=1,
                                          minMeasurementsPerCell=1,
                                          path_to_cnv="../test_data/cnv.pickle.gz",
                                          cnv_column_name="state")
        assert os.path.exists("../test_out/test1_mp1_mm1_cellNames.txt")
        assert os.path.exists("../test_out/test1_mp1_mm1_geneNames.txt")
        assert os.path.exists("../test_out/test1_mp1_mm1_siclonefit_input.txt")
        assert os.path.exists("../test_out/test1_mp1_mm1_siclonefit_raw.pickle")
        self.assertEqual(cmd, out)

    def test_003_siclonefit(self):
        """
        Test external downloaded binary siclonefit imputation algorithm. Check output file count.
        """
        cmd = "java -jar ../hamimzafar-siclonefit/SiCloneFiTComplete.jar -m 27 -n 6 -fp 0.001 -fn 0.0001 " \
              "-r 1 -df 0 -ipMat ../test_out//test1_mp1_mm1_siclonefit_input.txt -missing 0.24691 " \
              "-cellNames ../test_out//test1_mp1_mm1_cellNames.txt " \
              "-geneNames ../test_out//test1_mp1_mm1_geneNames.txt " \
              "-outDir ../test_out/"
        print(cmd)
        cmd_list = shlex.split(cmd)
        subprocess.check_call(cmd_list)
        # self.assertEqual(subprocess.check_call(cmd_list), 0)

    def test_004_missing_percentage(self):
        """
        Test missing percentage calculation
        """
        snvmatrix = pd.read_pickle("../test_out/test1_mp1_mm1_siclonefit_raw.pickle")
        missing = formatting.missing_percentage(snvmatrix)
        print(missing)
        self.assertEqual(missing, "24p_missing_samples")

    def test_005_convert_siclonefit_result(self):
        """
        Test missing percentage calculation from snvmatrix.
        Test convert back siclonefit result to pandas dataframe and has same size as input dataframe"""
        formatting.make_dir("../unittest_005/")
        shutil.copyfile("../test_out/test1_mp1_mm1_siclonefit_raw.pickle",
                        "../unittest_005/test1_mp1_mm1_siclonefit_raw.pickle")
        shutil.copytree("../test_out/24p_missing_samples",
                        "../unittest_005/24p_missing_samples")

        snvmatrix, imputed_snvmatrix = formatting.convert_siclonefit_result("../unittest_005/test1_mp1_mm1_siclonefit_raw.pickle",
                                                                            "../unittest_005/",
                                                                            minPresence=1,
                                                                            minMeasurementsPerCell=1)
        self.assertEqual(snvmatrix.shape, imputed_snvmatrix.shape)

    def test_006_plotting(self):
        """
        Test plotting output file count
        """
        prefix = "test1_mp1_mm1"
        snvmatrix = pd.read_pickle("../test_out/test1_mp1_mm1_siclonefit_raw.pickle")
        imputed_snvmatrix = pd.read_pickle("../test_out/test1_mp1_mm1_siclonefit_imputed.pickle")
        formatting.make_dir("../unittest_006/")
        cnv = pd.read_pickle("../test_data/cnv.pickle.gz")
        csplot = pm.CellSnvPlotMatrix(snvmatrix,
                                      imputed_snvmatrix,
                                      1,
                                      1,
                                      "../unittest_006/",
                                      "prefix",
                                      cnv,
                                      None,
                                      0.1,
                                      False)

        csplot.plot_snv_cell(sorted=True)
        csplot.plot_snv_cell(sorted=False)
        csplot.keptsSNVs_for_plotting.to_pickle(f"../unittest_006/{prefix}_keptsSNVs_for_plotting.pickle")
        out = os.listdir("../unittest_006")
        #print(out)
        self.assertEqual(len(out), 8)

    def test_007_siclonefitIO(self):
        """
        Test the whole process. Test cli.py
        """
        cmd = "siclonefit -j ../hamimzafar-siclonefit/SiCloneFiTComplete.jar -s ../test_data/test1.pickle -cn ../test_data/cnv.pickle.gz -o ../unittest_007/ -n test1 -mm 1 -mp 1"
        cmd_list = shlex.split(cmd)
        try:
            subprocess.run(cmd_list, check=True)
        except subprocess.CalledProcessError as cpe:
            print(cpe)
            assert False
        self.assertEqual(os.listdir("../unittest_006"), os.listdir("../test_out"))


