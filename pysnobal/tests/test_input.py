import unittest

import numpy as np

from pysnobal.input import InputData
from pysnobal.constants import FREEZE


class TestInput(unittest.TestCase):

    def compare_gold(self, gold, data):

        for key, value in gold.items():
            np.testing.assert_almost_equal(
                getattr(data, key),
                value,
                decimal=7,
                err_msg="{} not equal".format(key)
            )

    def setUp(self):
        self.data = {
            'I_lw': np.array([[237.]]),
            'S_n': np.array([[0.]]),
            'T_a': np.array([[277.16]]),
            'T_g': np.array([[273.16]]),
            'T_pp': np.array([[270.66]]),
            'e_a': np.array([[496.16]]),
            'm_pp': np.array([[0.]]),
            'percent_snow': np.array([[1.]]),
            'rho_snow': np.array([[150.]]),
            'u': np.array([[2.3]])
        }

        self.input = InputData(self.data)

    def test_access(self):

        self.compare_gold(self.data, self.input)

    def test_precip_cold_snow(self):

        self.data['m_pp'] = np.array([[1]])

        gold = {
            'm_rain': np.array([[0.]]),
            'm_snow': np.array([[1.]]),
            'z_snow': np.array([[0.00666667]]),
            'T_rain': np.array([[0]]),
            'T_snow': np.array([[270.66]]),
            'h2o_sat_snow': np.array([[0.]])
        }

        data = InputData(self.data)

        self.compare_gold(gold, data)

    def test_precip_mixed(self):

        self.data['m_pp'] = np.array([[1]])
        self.data['T_pp'] = np.array([[274]])
        self.data['percent_snow'] = np.array([[0.5]])

        gold = {
            'm_rain': np.array([[0.5]]),
            'm_snow': np.array([[0.5]]),
            'z_snow': np.array([[0.00333333]]),
            'T_rain': np.array([[274]]),
            'T_snow': np.array([[FREEZE]]),
            'h2o_sat_snow': np.array([[1.]])
        }

        data = InputData(self.data)

        self.compare_gold(gold, data)

    def test_precip_warm_snow(self):

        self.data['m_pp'] = np.array([[1]])
        self.data['T_pp'] = np.array([[274]])
        self.data['percent_snow'] = np.array([[1.]])

        gold = {
            'm_rain': np.array([[0.]]),
            'm_snow': np.array([[1.]]),
            'z_snow': np.array([[0.00666667]]),
            'T_rain': np.array([[0]]),
            'T_snow': np.array([[FREEZE]]),
            'h2o_sat_snow': np.array([[1.]])
        }

        data = InputData(self.data)

        self.compare_gold(gold, data)

    def test_precip_rain(self):

        self.data['m_pp'] = np.array([[1]])
        self.data['T_pp'] = np.array([[280]])
        self.data['percent_snow'] = np.array([[0.]])

        gold = {
            'm_rain': np.array([[1.]]),
            'm_snow': np.array([[0.]]),
            'z_snow': np.array([[0.]]),
            'T_rain': np.array([[280]]),
            'T_snow': np.array([[0]]),
            'h2o_sat_snow': np.array([[0.]])
        }

        data = InputData(self.data)

        self.compare_gold(gold, data)
