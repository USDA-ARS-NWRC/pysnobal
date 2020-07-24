import numpy as np

from pysnobal.constants import FREEZE


class InputData():

    # These will act like cumulative variables where the
    # input deltas will add to them
    INPUT_VARIABLES = [
        'S_n',
        'I_lw',
        'T_a',
        'e_a',
        'u',
        'T_g',
        'm_pp'
    ]

    # Some of the precipitation variables are handled slightly
    # differently where they are not added to the values before
    # but are just split evenly by the intervals
    PRECIP_VARIABLES = [
        'm_pp',
        'm_snow',
        'm_rain',
        'z_snow'
    ]

    # These precipitation variables are assumed constant over the
    # timestep
    PRECIP_CONSTANT = [
        'percent_snow',
        'rho_snow',
        'T_pp',
        'T_rain',
        'T_snow'
    ]

    # These are derived variables from the precipitation inputs
    # and will be calculated
    PRECIP_DERIVED = [
        'z_snow',
        'h2o_sat_snow',
        'T_rain',
        'T_snow'
    ]

    def __init__(self, data, input_delta=False):

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
        init = 0
        for precip_derived in self.PRECIP_DERIVED:
            setattr(self, precip_derived, init)

        if not input_delta:
            self.precipitation_inputs()

    def precipitation_inputs(self):

        self.precip_now = False

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
                self.z_snow = 0

            # check the precip, temp. cannot be below freezing if rain present
            if (self.m_rain > 0) and (self.T_pp < FREEZE):
                self.T_pp = FREEZE

            # Mixed snow and rain
            if (self.m_snow > 0) and (self.m_rain > 0):
                self.T_snow = FREEZE
                self.h2o_sat_snow = 1
                self.T_rain = self.T_pp

            elif (self.m_snow > 0):
                # Snow only
                if (self.T_pp < FREEZE):
                    # cold snow
                    self.T_snow = self.T_pp
                    self.h2o_sat_snow = 0
                else:
                    # warm snow
                    self.T_snow = FREEZE
                    self.h2o_sat_snow = 1

            elif (self.m_rain > 0):
                # rain only
                self.T_rain = self.T_pp

    def add_deltas(self, input_deltas):

        # Add the input data deltas
        self.S_n = self.S_n + input_deltas.S_n
        self.I_lw = self.I_lw + input_deltas.I_lw
        self.T_a = self.T_a + input_deltas.T_a
        self.e_a = self.e_a + input_deltas.e_a
        self.u = self.u + input_deltas.u
        self.T_g = self.T_g + input_deltas.T_g

        self.update_precip_deltas(input_deltas)

    def update_precip_deltas(self, input_deltas):

        # update the precipitation. Snobal takes the input deltas
        # and divides by the intervals
        for precip_variable in self.PRECIP_VARIABLES:
            setattr(self, precip_variable,
                    getattr(input_deltas, precip_variable))


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
                ) / tstep['intervals_per_timestep']

            for precip_variable in self.input1.PRECIP_VARIABLES:
                tstep_deltas[precip_variable] = \
                    getattr(self.input1, precip_variable) / \
                    tstep['intervals_per_timestep']

            for precip_constant in self.input1.PRECIP_CONSTANT:
                tstep_deltas[precip_constant] = \
                    getattr(self.input1, precip_constant)

            self.deltas[tstep['level']] = InputData(tstep_deltas, True)

        return self.deltas
