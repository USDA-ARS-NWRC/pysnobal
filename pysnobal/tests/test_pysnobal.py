import os
import unittest
from copy import deepcopy
from pathlib import Path

import pandas as pd
from inicheck.tools import MasterConfig, cast_all_variables, get_user_config

import pysnobal
from pysnobal.pysnobal import PySnobal


class TestPysnobal(unittest.TestCase):

    BASE_INI_FILE_NAME = 'pysnobal_config.ini'

    test_dir = Path(pysnobal.__file__).parent.joinpath('tests')
    config_file = os.path.join(test_dir, BASE_INI_FILE_NAME)

    @property
    def base_config(self):
        return self.base_config_copy()

    @classmethod
    def base_config_copy(cls):
        return deepcopy(cls._base_config)

    @classmethod
    def load_base_config(cls):
        master_config = os.path.join(
            Path(pysnobal.__file__).parent, 'pysnobal_core_config.ini')
        mcfg = MasterConfig(path=master_config)
        cls._base_config = get_user_config(cls.config_file, mcfg=mcfg)

    @classmethod
    def setUpClass(cls):
        cls.load_base_config()

    def setUp(self):
        os.makedirs('pysnobal/tests/output', exist_ok=True)

    def tearDown(self):
        os.remove('pysnobal/tests/output/pysnobal_output.csv')

    def test_pysnobal_output_normal(self):

        # run PySnobal
        status = PySnobal(self.base_config).run()
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

    def test_pysnobal_ouput_all(self):

        config = self.base_config_copy()
        config.raw_cfg['files'].update({'output_mode': 'all'})

        config.apply_recipes()
        config = cast_all_variables(config, config.mcfg)

        # run PySnobal
        status = PySnobal(config).run()
        self.assertTrue(status)

        # load in the outputs
        gold = pd.read_csv(
            'pysnobal/tests/test_data_point/gold_csv/gold.pysnobal.all.csv',
            index_col='date_time', parse_dates=True)
        gold.index = gold.index.tz_convert('MST')

        new = pd.read_csv(
            'pysnobal/tests/output/pysnobal_output.csv',
            index_col='date_time', parse_dates=True)
        new.index = new.index.tz_convert('MST')

        pd.testing.assert_frame_equal(gold, new)
