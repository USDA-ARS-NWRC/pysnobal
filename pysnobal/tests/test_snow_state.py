import unittest

from pysnobal.constants import FREEZE
from pysnobal.snow_state import SnowState


class TestSnowState(unittest.TestCase):

    def test_T_s_0(self):

        ss = SnowState()
        ss.T_s_0 = 10
        self.assertTrue(ss.T_s_0 == 10)

    def test_e_s_0(self):

        ss = SnowState()
        ss.T_s_0 = FREEZE

        # This will call the e_s_0 property method
        self.assertTrue(ss.e_s_0 == 610.71)
        self.assertTrue(ss.T_s_0 == FREEZE)

        ss.T_s_0 = FREEZE - 10
        self.assertAlmostEqual(ss.e_s_0, 259.7018533254676)
