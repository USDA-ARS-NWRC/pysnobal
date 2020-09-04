from pysnobal.core.constants import FREEZE
from pysnobal.point.libsnobal import sati
from pysnobal.point.input import InputData, InputDeltas


class InputSpatialData(InputData):

    def __init__(self, data, input_delta=False):

        self.S_n = data['S_n']
        self.I_lw = data['I_lw']
        self.t_a = data['t_a']
        self.e_a = data['e_a']
        self.u = data['u']
        self.t_g = data['t_g']
        self.m_pp = data['m_pp']
        self.percent_snow = data['percent_snow']
        self.rho_snow = data['rho_snow']
        self.t_pp = data['t_pp']

        # derived precip values
        self.m_snow = self.m_pp * self.percent_snow
        self.m_rain = self.m_pp - self.m_snow

        # initialize the other variables to 0
        init = 0
        for precip_derived in self.PRECIP_DERIVED:
            setattr(self, precip_derived, init)

        if not input_delta:
            self.precipitation_inputs()

    @property
    def t_a(self):
        return self._t_a

    @t_a.setter
    def t_a(self, var):
        self._t_a = var
        self.__sat_vp = False

    @property
    def sat_vp(self):
        """Calculate the saturation vapor pressure over ice for
        the air temperature

        Returns:
            float: saturation vapor pressure over ice
        """
        if not self.__sat_vp:
            self._sat_vp = sati(self.t_a)
            self.__sat_vp = True
        return self._sat_vp

    @property
    def t_g(self):
        return self._t_g

    @t_g.setter
    def t_g(self, var):
        self._t_g = var
        self.__e_g = False

    @property
    def e_g(self):
        """Calculate the saturation vapor pressure over ice for
        the soil temperature

        Returns:
            float: saturation vapor pressure over ice
        """
        if not self.__e_g:
            self._e_g = sati(self.t_g)
            self.__e_g = True
        return self._e_g

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
            if (self.m_rain > 0) and (self.t_pp < FREEZE):
                self.t_pp = FREEZE

            # Mixed snow and rain
            if (self.m_snow > 0) and (self.m_rain > 0):
                self.t_snow = FREEZE
                self.h2o_sat_snow = 1
                self.t_rain = self.t_pp

            elif (self.m_snow > 0):
                # Snow only
                if (self.t_pp < FREEZE):
                    # cold snow
                    self.t_snow = self.t_pp
                    self.h2o_sat_snow = 0
                else:
                    # warm snow
                    self.t_snow = FREEZE
                    self.h2o_sat_snow = 1

            elif (self.m_rain > 0):
                # rain only
                self.t_rain = self.t_pp

    def add_deltas(self, input_deltas):

        # Add the input data deltas
        self.S_n = self.S_n + input_deltas.S_n
        self.I_lw = self.I_lw + input_deltas.I_lw
        self.t_a = self.t_a + input_deltas.t_a
        self.e_a = self.e_a + input_deltas.e_a
        self.u = self.u + input_deltas.u
        self.t_g = self.t_g + input_deltas.t_g

        self.update_precip_deltas(input_deltas)

    def update_precip_deltas(self, input_deltas):

        # update the precipitation. Snobal takes the input deltas
        # and divides by the intervals
        for precip_variable in self.PRECIP_VARIABLES:
            setattr(self, precip_variable,
                    getattr(input_deltas, precip_variable))


class InputSpatialDeltas(InputDeltas):

    def __init__(self, input1, input2, tstep_info):

        super(InputSpatialData, self).__init__(
            input1, input2, tstep_info)

    def calculate(self):

        self.deltas = {}
        for tstep in self.tstep_info:

            tstep_deltas = {}
            for variable in INPUT_VARIABLES:
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

            self.deltas[tstep['level']] = InputData(tstep_deltas)

        return self.deltas
