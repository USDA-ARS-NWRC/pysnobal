import pytest
from pytest import approx

from pysnobal.core.constants import FREEZE
from pysnobal.point import SnowState


@pytest.fixture
def snow_state():
    return SnowState()


def test_T_s_0(snow_state):

    snow_state.T_s_0 = 10
    assert snow_state.T_s_0 == 10


def test_e_s_0(snow_state):

    snow_state.T_s_0 = FREEZE

    # This will call the e_s_0 property method
    assert snow_state.e_s_0 == 610.71
    assert snow_state.T_s_0 == FREEZE

    snow_state.T_s_0 = FREEZE - 10
    assert approx(snow_state.e_s_0, 259.7018533254676)
