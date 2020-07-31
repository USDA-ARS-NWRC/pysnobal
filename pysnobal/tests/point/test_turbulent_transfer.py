import unittest

import pandas as pd

from pysnobal.core.constants import (GRAVITY, MOL_AIR, SEA_LEVEL, STD_AIRTMP,
                                     STD_LAPSE)
from pysnobal.point import libsnobal


class TestTurbulentTransfer(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.data = pd.read_csv(
            'pysnobal/tests/test_data_point/libsnobal/gold_turbulent_transfer.csv')  # noqa

    def calculate_hle1(self, inputs):

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

    def test_hle1_gold(self):

        for index, row in self.data.iterrows():

            H, L_v_E, E, status = self.calculate_hle1(row.to_dict())

            self.assertAlmostEqual(H, row.py_H)
            self.assertAlmostEqual(L_v_E, row.py_L_v_E)
            self.assertAlmostEqual(E, row.py_E)
            self.assertAlmostEqual(status, row.py_status)
