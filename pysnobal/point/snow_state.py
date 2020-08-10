from pysnobal.core.constants import CAL_TO_J, FREEZE, RHO_ICE, RHO_W0, MIN_SNOW_TEMP
from pysnobal.core.functions import cp_ice, time_average, heat_stor
from pysnobal.point.libsnobal import sati


class SnowState():

    # time averaged values
    _energy_state = [
        'R_n',      # net allwave radiation (W/m^2)
        'H',        # sensible heat xfr (W/m^2)
        'L_v_E',    # latent heat xfr (W/m^2)
        'G',        # heat transfer from soil to snowcover (W/m^2)
        'G_0',      # heat transfer soil or lower layer to active layer (W/m^2)
        'M',        # advected heat from precip (W/m^2)
        'delta_Q',  # change in snowcover's energy(W/m ^ 2)
        'delta_Q_0'  # change in active layer's energy(W/m ^ 2)
    ]

    def __init__(self, init=0):

        self.zeros = init

        # snowpack state variables
        self.h2o_sat = init
        self.layer_count = init
        self.max_h2o_vol = init
        self.rho = init
        self.t_s = init
        self.t_s_0 = init
        self.t_s_l = init
        self.z_s = init
        self.z_s_0 = init
        self.z_s_l = init

        # TODO could be moved to a property that is calculated when needed
        # self.cc_s = init
        self.cc_s_0 = init
        self.cc_s_l = init
        self.m_s = init
        self.m_s_0 = init
        self.m_s_l = init
        self.h2o = init
        self.h2o_max = init
        self.h2o_total = init
        self.h2o_vol = init

        # Snow energetics state variables
        for variable in self._energy_state:
            setattr(self, variable, init)

            bar_variable = "{}_bar".format(variable)
            setattr(self, bar_variable, init)

        # mass balance vars for current timestep
        # specific melt (kg/m^2 or mm)
        self.melt = init

        # mass flux by evap into air from active layer (kg/m^2/s)
        self.E = init

        # mass of evap into air & soil from snowcover (kg/m^2)
        self.E_s = init

        # predicted specific runoff (m/sec)
        self.swi = init

        # sums of mass balance vars since last output record
        self.melt_sum = init
        self.E_s_sum = init

        self.z_0 = 0.05
        self._isothermal = False

    def set_zeros(self, fields):

        if isinstance(fields, str):
            fields = fields.split()

        for field in fields:
            setattr(self, field, self.zeros)

    def set_from_dict(self, data_dict):
        """Set snow state variables from a dict object

        Args:
            data_dict (dict): data dict to set in snow state
        """

        for variable, data in data_dict.items():
            setattr(self, variable, data)

    @property
    def t_s_0(self):
        return self._t_s_0

    @t_s_0.setter
    def t_s_0(self, var):
        self._t_s_0 = var
        self.__es_s_0 = False

    @property
    def e_s_0(self):
        """Calculate the saturation vapor pressure over ice for
        the surface layer.

        Returns:
            float: saturation vapor pressure over ice
        """
        if not self.__es_s_0:
            self._es_s_0 = sati(self.t_s_0)
            self.__es_s_0 = True
        return self._es_s_0

    @property
    def t_s_l(self):
        return self._t_s_l

    @t_s_l.setter
    def t_s_l(self, var):
        self._t_s_l = var
        self.__es_s_l = False

    @property
    def e_s_l(self):
        """Calculate the saturation vapor pressure over ice for
        the lower layer.

        Returns:
            float: saturation vapor pressure over ice
        """
        if not self.__es_s_l:
            self._es_s_l = sati(self.t_s_l)
            self.__es_s_l = True
        return self._es_s_l

    @property
    def isothermal(self):

        if (self.layer_count == 2) and (self.cc_s_0 == 0.0) and \
                (self.cc_s_l == 0.0):
            self._isothermal = True
        elif (self.layer_count == 1) and (self.cc_s_0 == 0.0):
            self._isothermal = True
        else:
            self._isothermal = False

        return self._isothermal

    # @property
    # def cc_s_0(self):
    #     if self.t_s_0 == MIN_SNOW_TEMP + FREEZE:
    #         self._cc_s_0 = 0
    #     else:
    #         self._cc_s_0 = self.cold_content(self.t_s_0, self.m_s_0)

    #     return val

    # @cc_s_0.setter
    # def cc_s_0(self, val):
    #     self._cc_s_0 = val

    # @property
    # def cc_s_l(self):
    #     return self.cold_content(self.t_s_l, self.m_s_l)

    @property
    def cc_s(self):
        cc_s = 0
        if self.layer_count == 2:
            cc_s = self.cc_s_0 + self.cc_s_l
        elif self.layer_count == 1:
            cc_s = self.cc_s_0
        return cc_s

    # def cold_content(self, temp, mass):
    #     """
    #     This routine calculates the cold content for a layer (i.e., the
    #     energy required to bring its temperature to freezing) from the
    #     layer's temperature and specific mass.

    #     Args:
    #         temp: temperature of layer
    #         mass: specific mass of layer

    #     Returns:
    #         cc: cold content of layer
    #     """

    #     cc = 0
    #     if temp < FREEZE:
    #         cc = heat_stor(cp_ice(temp), mass, temp-FREEZE)
    #     return cc

    def adjust_layer_temps(self):
        """Adjust the layer temperatures
        """

        if self.layer_count == 1:
            self.t_s_0 = self.new_tsno(
                self.m_s_0,
                self.t_s_0,
                self.cc_s_0)
            self.t_s = self.t_s_0

        elif self.layer_count == 2:
            if self.isothermal:
                self.t_s = FREEZE
                self.t_s_l = FREEZE
                self.t_s_0 = FREEZE
            else:
                self.t_s_0 = self.new_tsno(
                    self.m_s_0,
                    self.t_s_0,
                    self.cc_s_0)
                self.t_s_l = self.new_tsno(
                    self.m_s_l,
                    self.t_s_l,
                    self.cc_s_l)
                self.t_s = self.new_tsno(
                    self.m_s,
                    self.t_s,
                    self.cc_s)

    def new_tsno(self, spm, t0, ccon):
        """
        calculates a new temperature for a snow layer from its
        adjusted cold content, and the layer's last (previous) temperature.

        The layer's specific mass (the argument <I>spm</I>) can be computed by
        multiplying the layer's thickness (m) by its density (kg/m^3).

        Args:
            spm: layer's specific mass (kg/m^2)
            t0: layer's last temperature (K)
            ccon: layer's adjusted cold content (J/m^2)

        Returns:
            tsno: snow layer's new temperature (K)
        """

        cp = cp_ice(t0)
        tdif = ccon / (spm * cp)
        tsno = tdif + FREEZE

        return tsno

    @property
    def dry_snow_density(self):
        """dry snow density (without H2O) at a given total snow density
        (with H2O) and snow saturation

        Args:
            rhos (float or array): total density of snow (kg/m^3)
            sat (float or array): snow saturation (see SNO_SAT)

        Returns:
            float or array: dry snow density
        """
        return (self.rho - self.h2o_vol * RHO_W0) / \
            (1 - self.h2o_vol * RHO_W0 / RHO_ICE)

    @property
    def kts(self):
        """thermal conductivity of snow (J/(m sec K))
        (after Yen, 1965, see Anderson, 1976, pg. 31)

        Args:
            rho (float or array): snow density [kg/m3]

        Returns:
            float or array: thermal conductivity of snow
        """
        return CAL_TO_J * 0.0077 * (self.rho/1000.0) * (self.rho/1000.0)

    def value_to_bar(self):
        """Copy the current snow state variable value into the `_bar` value.
        """

        for variable in self._energy_state:
            setattr(self, "{}_bar".format(variable), getattr(self, variable))

        self.E_s_sum = self.E_s
        self.melt_sum = self.melt
        self.swi_sum = self.swi

    def time_average(self, time_since_out, time_step):
        """Update the time averaged value (`_bar`) for the desired variables.

        Args:
            time_since_out (float): seonds since last output
            time_step (float): model time step
        """

        for variable in self._energy_state:
            bar_variable = "{}_bar".format(variable)
            bar_value = getattr(self, bar_variable)
            ta = time_average(bar_value, time_since_out,
                              getattr(self, variable), time_step)
            setattr(self, bar_variable, ta)

        self.E_s_sum += self.E_s
        self.melt_sum += self.melt
        self.swi_sum += self.swi
