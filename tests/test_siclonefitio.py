#!/usr/bin/env python

"""Tests for `siclonefitio` package."""


import unittest
import os
import shutil

from siclonefitio import formatting
import visual


class TestSiclonefitio(unittest.TestCase):
    """Tests for `siclonefitio` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""
        shutil.rmtree("./test_make_dir/")

    def test_000_installation(self):
        """Test something."""
        formatting.make_dir("./test_make_dir/")
        assert os.path.exists("./test_make_dir/")

