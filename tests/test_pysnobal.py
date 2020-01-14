#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_pysnobal
----------------------------------

Tests for `pysnobal` module.
"""

import unittest
import pandas as pd
import numpy as np
import os

from pysnobal.snobal import PySnobal


class TestPysnobal(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        os.remove('tests/test_data_point/snobal.pysnobal_c')

    def test_pysnobal_run(self):
        """ Test PySnobal and compare with Snobal """

        # run PySnobal
        status = PySnobal().run()
        self.assertTrue(status)

        # load in the outputs
        gold = pd.read_csv('tests/test_data_point/gold_ipw/gold.snobal.out',
                           header=None, index_col=0)
        new = pd.read_csv(
            'tests/test_data_point/snobal.pysnobal_c', header=None, index_col=0)

        self.assertTrue(new.shape[0] == 8758)
        self.assertTrue(new.shape[1] == 25)

        result = np.abs(gold - new)
        self.assertFalse(np.any(result > 0))
