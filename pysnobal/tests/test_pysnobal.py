import os
import unittest

import numpy as np
import pandas as pd

from pysnobal.pysnobal import PySnobal


class TestPysnobal(unittest.TestCase):

    def setUp(self):
        os.makedirs('pysnobal/tests/output', exist_ok=True)

    def tearDown(self):
        os.remove('pysnobal/tests/output/pysnobal_output.csv')

    def test_pysnobal_run(self):

        # run PySnobal
        status = PySnobal(
            'pysnobal/tests/pysnobal_config.ini').run()
        self.assertTrue(status)

        # load in the outputs
        gold = pd.read_csv(
            'pysnobal/tests/test_data_point/gold_csv/gold.pysnobal.csv',
            index_col='date_time', parse_dates=True)
        gold.index = gold.index.tz_convert('MST')

        new = pd.read_csv(
            'pysnobal/tests/output/pysnobal_output.csv',
            index_col='date_time', parse_dates=True)
        new.index = new.index.tz_convert('MST')

        pd.testing.assert_frame_equal(gold, new)
