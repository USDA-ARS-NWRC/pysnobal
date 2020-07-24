import unittest
import os

import numpy as np
import pandas as pd

from pysnobal import libsnobal
from pysnobal.constants import FREEZE


class TestTurbulentTransfer(unittest.TestCase):    

    def calculate_hle1(self, inputs):

        P_a = libsnobal.hysat(
            libsnobal.SEA_LEVEL,
            libsnobal.STD_AIRTMP,
            libsnobal.STD_LAPSE,
            inputs['elevation'] / 1000.0,
            libsnobal.GRAVITY,
            libsnobal.MOL_AIR)

        H, L_v_E, E, status = libsnobal.hle1(
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

    def test_hle1(self):

        inputs = {
            'elevation': 3000,
            'ta': 10 + FREEZE,
            'ts': -1 + FREEZE,
            'za': 2,
            'ea': 1000,
            'es': 560,
            'u': 4,
            'zu': 5,
            'z0': 0.05
        }

        H, L_v_E, E, status = self.calculate_hle1(inputs)
        H
