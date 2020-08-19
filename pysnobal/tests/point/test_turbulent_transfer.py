import pandas as pd
import pytest
from pytest import approx

from pysnobal.core.constants import (GRAVITY, MOL_AIR, SEA_LEVEL, STD_AIRTMP,
                                     STD_LAPSE)
from pysnobal.core.functions import hysat
from pysnobal.point import libsnobal


@pytest.fixture
def gold_data():
    """Load the csv gold file data"""

    return pd.read_csv(
        'pysnobal/tests/test_data_point/libsnobal/gold_turbulent_transfer.csv')


def calculate_hle1(inputs):

    P_a = hysat(
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


def test_hle1_stable(gold_data):

    data = {
        'elevation': 1264.0,
        'ta': 256.16,
        'ts': 256.1990314427861,
        'za': 4.0,
        'ea': 20,
        'es': 137,
        'u': 38.0,
        'zu': 4.0,
        'z0': 0.01,
        'py_H': 0.0,
        'py_L_v_E': -510.5932737176731,
        'py_E': -0.00017688648886582492,
        'py_status': 0.0
    }

    H, L_v_E, E, status = calculate_hle1(data)

    assert approx(H, data['py_H'])
    assert approx(L_v_E, data['py_L_v_E'])
    assert approx(E, data['py_E'])
    assert approx(status, data['py_status'])
