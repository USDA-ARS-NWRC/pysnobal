import numpy as np

from pysnobal.constants import FREEZE


class InputData():

    INPUT_VARIABLES = [
        'S_n',
        'I_lw',
        'T_a',
        'e_a',
        'u',
        'T_g',
        'm_pp',
        'percent_snow',
        'rho_snow',
        'T_pp',
        'm_snow',
        'm_rain',
        'z_snow',
        'h2o_sat_snow',
        'T_rain',
        'T_snow'
    ]

    def __init__(self, data):

        self.S_n = data['S_n']
        self.I_lw = data['I_lw']
        self.T_a = data['T_a']
        self.e_a = data['e_a']
        self.u = data['u']
        self.T_g = data['T_g']
        self.m_pp = data['m_pp']
        self.percent_snow = data['percent_snow']
        self.rho_snow = data['rho_snow']
        self.T_pp = data['T_pp']

        # derived precip values
        self.m_snow = self.m_pp * self.percent_snow
        self.m_rain = self.m_pp - self.m_snow

        # initialize the other variables to 0
        init = np.zeros_like(self.m_pp)
        self.z_snow = init
        self.h2o_sat_snow = init
        self.T_rain = init
        self.T_snow = init

        if self.m_pp > 0:
            self.precip_now = True

            # self.m_pp = self.m_pp
            # self.m_snow = self.m_pp * self.percent_snow
            # self.m_rain = self.m_pp - self.m_snow

            if (self.m_snow > 0.0):
                if (self.rho_snow > 0.0):
                    self.z_snow = self.m_snow / self.rho_snow
                else:
                    raise ValueError(
                        'rho_snow is <= 0.0 with percent_snow > 0.0')
            else:
                self.z_snow = np.zeros_like(self.T_pp)

            # check the precip, temp. cannot be below freezing if rain present
            if (self.m_rain > 0) and (self.T_pp < FREEZE):
                self.T_pp = FREEZE

            # Mixed snow and rain
            if (self.m_snow > 0) and (self.m_rain > 0):
                self.T_snow = FREEZE
                self.h2o_sat_snow = np.ones_like(self.T_pp)
                self.T_rain = self.T_pp

            elif (self.m_snow > 0):
                # Snow only
                if (self.T_pp < FREEZE):
                    # cold snow
                    self.T_snow = self.T_pp
                    self.h2o_sat_snow = np.zeros_like(self.T_pp)
                else:
                    # warm snow
                    self.T_snow = FREEZE
                    self.h2o_sat_snow = np.ones_like(self.T_pp)

            elif (self.m_rain > 0):
                # rain only
                self.T_rain = self.T_pp


class InputDeltas():

    def __init__(self, input1, input2, tstep_info):

        self.input1 = input1
        self.input2 = input2
        self.tstep_info = tstep_info

    def calculate(self):

        self.deltas = {}
        for tstep in self.tstep_info:

            tstep_deltas = {}
            for variable in self.input1.INPUT_VARIABLES:
                tstep_deltas[variable] = (
                    getattr(self.input2, variable) -
                    getattr(self.input1, variable)
                ) / tstep['intervals']

            self.deltas[tstep['level']] = InputData(tstep_deltas)

        return self.deltas
