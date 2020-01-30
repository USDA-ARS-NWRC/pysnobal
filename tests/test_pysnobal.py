#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_pysnobal
----------------------------------

Tests for `pysnobal` module.
"""

import os
import unittest

import numpy as np
import pandas as pd

from pysnobal.snobal import PySnobal


class TestPysnobal(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass
        # os.remove('tests/test_data_point/output.csv')

    def test_pysnobal_run(self):
        """ Test PySnobal and compare with Snobal """

        # run PySnobal
        status = PySnobal('tests/test_data_point/gold_csv/config.ini').run()
        self.assertTrue(status)

        # load in the outputs
        gold = pd.read_csv(
            'tests/test_data_point/gold_csv/gold.snobal.out.csv',
            index_col='date_time', parse_dates=True)
        gold.index = gold.index.tz_localize('MST')

        new = pd.read_csv(
            'tests/test_data_point/gold_csv/output.csv',
            index_col='date_time', parse_dates=True)

        self.assertTrue(new.shape == gold.shape)

        result = np.abs(gold - new)
        self.assertFalse(np.any(result > 0))
