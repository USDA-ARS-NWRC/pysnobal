import numpy as np
import pytest

from pysnobal.core.constants import FREEZE
from pysnobal.point import InputData


@pytest.fixture
def data():
    """Initial input data values"""
    return {
        'thermal': np.array([[237.]]),
        'net_solar': np.array([[0.]]),
        'air_temp': np.array([[277.16]]),
        'soil_temp': np.array([[273.16]]),
        'precip_temp': np.array([[270.66]]),
        'vapor_pressure': np.array([[496.16]]),
        'precip_mass': np.array([[0.]]),
        'percent_snow': np.array([[1.]]),
        'rho_snow': np.array([[150.]]),
        'wind_speed': np.array([[2.3]])
    }


@pytest.fixture
def initial_input(data):
    """InputData class for the data"""
    return InputData(data)


def compare_gold(gold, data):

    for key, value in gold.items():
        np.testing.assert_almost_equal(
            getattr(data, key),
            value,
            decimal=7,
            err_msg="{} not equal".format(key)
        )


def test_access(data, initial_input):
    compare_gold(data, initial_input)


def test_precip_cold_snow(data):

    data['precip_mass'] = np.array([[1]])

    gold = {
        'mass_rain': np.array([[0.]]),
        'mass_snow': np.array([[1.]]),
        'z_snow': np.array([[0.00666667]]),
        'temp_rain': np.array([[0]]),
        'temp_snow': np.array([[270.66]]),
        'h2o_sat_snow': np.array([[0.]])
    }

    data = InputData(data)

    compare_gold(gold, data)


def test_precip_mixed(data):

    data['precip_mass'] = np.array([[1]])
    data['precip_temp'] = np.array([[274]])
    data['percent_snow'] = np.array([[0.5]])

    gold = {
        'mass_rain': np.array([[0.5]]),
        'mass_snow': np.array([[0.5]]),
        'z_snow': np.array([[0.00333333]]),
        'temp_rain': np.array([[274]]),
        'temp_snow': np.array([[FREEZE]]),
        'h2o_sat_snow': np.array([[1.]])
    }

    data = InputData(data)

    compare_gold(gold, data)


def test_precip_warmass_snow(data):

    data['precip_mass'] = np.array([[1]])
    data['precip_temp'] = np.array([[274]])
    data['percent_snow'] = np.array([[1.]])

    gold = {
        'mass_rain': np.array([[0.]]),
        'mass_snow': np.array([[1.]]),
        'z_snow': np.array([[0.00666667]]),
        'temp_rain': np.array([[0]]),
        'temp_snow': np.array([[FREEZE]]),
        'h2o_sat_snow': np.array([[1.]])
    }

    data = InputData(data)

    compare_gold(gold, data)


def test_precip_rain(data):

    data['precip_mass'] = np.array([[1]])
    data['precip_temp'] = np.array([[280]])
    data['percent_snow'] = np.array([[0.]])

    gold = {
        'mass_rain': np.array([[1.]]),
        'mass_snow': np.array([[0.]]),
        'z_snow': np.array([[0.]]),
        'temp_rain': np.array([[280]]),
        'temp_snow': np.array([[0]]),
        'h2o_sat_snow': np.array([[0.]])
    }

    data = InputData(data)

    compare_gold(gold, data)
