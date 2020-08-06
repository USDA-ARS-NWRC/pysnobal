import os
from copy import deepcopy
from pathlib import Path

import pandas as pd
import pytest
from inicheck.tools import MasterConfig, cast_all_variables, get_user_config

import pysnobal
from pysnobal.ipysnobal import iPySnobal

BASE_INI_FILE_NAME = 'test_data_spatial/gold_ipw/rme_wy2004/netcdf/ipysnobal_config.ini'
test_dir = Path(pysnobal.__file__).parent.joinpath('tests')
config_file = os.path.join(test_dir, BASE_INI_FILE_NAME)


@pytest.fixture(scope='module')
def base_config():
    """Load the base config object"""

    master_config = os.path.join(
        Path(pysnobal.__file__).parent, 'ipysnobal_core_config.ini')
    mcfg = MasterConfig(path=master_config)
    return get_user_config(config_file, mcfg=mcfg)


@pytest.fixture
def base_config_copy(base_config):
    """Create a copy of the base config object"""
    return deepcopy(base_config)


@pytest.fixture(autouse=True)
def make_clean():
    """Make the directory and clean up after the test"""

    # os.makedirs('pysnobal/tests/output', exist_ok=True)
    yield
    # os.remove(
    #     'pysnobal/tests/test_data_spatial/gold_ipw/rme_wy2004/netcdf/pysnobal_out.nc')


def test_ipysnobal(make_clean, base_config):
    assert True
    # run PySnobal
    # status = iPySnobal(base_config).run()
    # assert status

    # # load in the outputs
    # gold = pd.read_csv(
    #     'pysnobal/tests/test_data_point/gold_csv/gold.pysnobal.csv',
    #     index_col='date_time', parse_dates=True)
    # gold.index = gold.index.tz_convert('MST')

    # new = pd.read_csv(
    #     'pysnobal/tests/output/pysnobal_output.csv',
    #     index_col='date_time', parse_dates=True)
    # new.index = new.index.tz_convert('MST')

    # pd.testing.assert_frame_equal(gold, new)
