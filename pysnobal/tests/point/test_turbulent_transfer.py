import pandas as pd
import pytest
from pytest import approx

from pysnobal.core.constants import (GRAVITY, MOL_AIR, SEA_LEVEL, STD_AIRTMP,
                                     STD_LAPSE)
from pysnobal.point import libsnobal


@pytest.fixture
def gold_data():
    """Load the csv gold file data"""

    return pd.read_csv(
        'pysnobal/tests/test_data_point/libsnobal/gold_turbulent_transfer.csv')


def calculate_hle1(inputs):

    P_a = libsnobal.hysat(
        SEA_LEVEL,
        STD_AIRTMP,
        STD_LAPSE,
        inputs['elevation'] / 1000.0,
        GRAVITY,
        MOL_AIR)

    H, L_v_E, E, status, ustar, factor = libsnobal.hle1(
        P_a,
        inputs['ta'],
        inputs['ts'],
        inputs['za'],
        inputs['ea'],
        inputs['es'],
        inputs['za'],
        inputs['u'],
        inputs['zu'],
        inputs['z0']
    )

    return H, L_v_E, E, status


def test_hle1_gold(gold_data):

    for index, row in gold_data.iterrows():

        H, L_v_E, E, status = calculate_hle1(row.to_dict())

        assert approx(H, row.py_H)
        assert approx(L_v_E, row.py_L_v_E)
        assert approx(E, row.py_E)
        assert approx(status, row.py_status)
