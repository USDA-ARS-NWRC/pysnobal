import pytest
from pytest import approx

from pysnobal.core.constants import FREEZE
from pysnobal.point import SnowState


@pytest.fixture
def snow_state():
    return SnowState()


def test_t_s_0(snow_state):

    snow_state.t_s_0 = 10
    assert snow_state.t_s_0 == 10


def test_e_s_0(snow_state):

    snow_state.t_s_0 = FREEZE

    # This will call the e_s_0 property method
    assert snow_state.e_s_0 == 610.71
    assert snow_state.t_s_0 == FREEZE

    snow_state.t_s_0 = FREEZE - 10
    assert approx(snow_state.e_s_0, 259.7018533254676)


def test_calc_layers_zero(snow_state):
    snow_state.m_s = 0.1
    snow_state.calc_layers()

    assert snow_state.z_s == 0
    assert snow_state.z_s_0 == 0
    assert snow_state.z_s_l == 0
    assert snow_state.layer_count == 0


def test_calc_layers_one(snow_state):
    snow_state.z_s = 0.15
    snow_state.m_s = 8
    snow_state.calc_layers()

    assert snow_state.z_s == 0.15
    assert snow_state.z_s_0 == 0.15
    assert snow_state.z_s_l == 0
    assert snow_state.layer_count == 1


def test_calc_layers_almost_two(snow_state):
    snow_state.z_s = 0.25001
    snow_state.rho = 300
    snow_state.m_s = snow_state.z_s * snow_state.rho
    snow_state.calc_layers()

    assert snow_state.z_s == 0.25001
    assert snow_state.z_s_0 == 0.25001
    assert snow_state.z_s_l == 0
    assert snow_state.layer_count == 1


def test_calc_layers_two(snow_state):
    snow_state.z_s = 1.5
    snow_state.rho = 300
    snow_state.m_s = snow_state.z_s * snow_state.rho
    snow_state.calc_layers()

    assert snow_state.z_s == 1.5
    assert snow_state.z_s_0 == 0.25
    assert snow_state.z_s_l == 1.25
    assert snow_state.layer_count == 2
