from pysnobal.point.libsnobal import sati
from pysnobal.core.functions import time_average
from pysnobal.core.snow import new_tsno, h2o_left, cold_content
from pysnobal.core.constants import CAL_TO_J, FREEZE, RHO_ICE, RHO_W0


class SnowState():
    """
    Attributes:
        layer_count
        max_h2o_vol: maximum liquid water content as volume ratio, used to
            predict swi V_water/(V_snow-V_ice)
        rho
        t_s
        t_s_0
        t_s_l
        z_s
        z_s_0
        z_s_l
        cc_s_0
        cc_s_l
        m_s
        m_s_0
        m_s_l
        h2o_sat: percent saturation, 0 to 1.0
        h2o: actual liquid water content (as specific mass, h2o_sat*h2o_max)
        h2o_max: maximum liquid water content (as specific mass)
        h2o_total: total amount of liquid water in snowpack, value greater than
            h2o_max will become swi
        h2o_vol: liquid water content as volume ratio (h2o_sat*max_h2o_vol)
    """

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

    def __init__(self, init=0, max_z_s_0=0.25, small_threshold=1):

        self.zeros = init
        self.max_z_s_0 = max_z_s_0
        self.small_threshold = small_threshold

        # snowpack state variables
        self.h2o_sat = init
        # self.layer_count = init
        self.max_h2o_vol = init
        self.rho = init
        self.t_s = init
        self.t_s_0 = init
        self.t_s_l = init
        self.z_s = init
        self.z_s_0 = init
        self.z_s_l = init

        # TODO could be moved to a property that is calculated when needed
        self.cc_s_0 = init
        self.cc_s_l = init
        self.m_s = init
        self.m_s_0 = init
        self.m_s_l = init
        self.h2o = init
        # self.h2o_max = init
        self.h2o_total = init
        # self.h2o_vol = init

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
        self.__layer_count = False

        # While this is a little redundant, Snobal can lose the snowcover during
        # a timestep but the self.snowcover is not updated like the snow_state. So
        # the snowcover will only be updated when it comes through the do_tstep
        # function. The swi needs the tstep_snowcover.
        self.tstep_snowcover = False

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
        return self.cc_s == 0.0

    @property
    def snowcover(self):
        return self.layer_count > 0

    @property
    def cc_s(self):
        return self.cc_s_0 + self.cc_s_l

    @property
    def layer_count(self):
        if not self.__layer_count:
            self._layer_count = sum([self.z_s_0 > 0, self.z_s_l > 0])
            self.__layer_count = True

        return self._layer_count

    def adjust_layer_temps(self):
        """Adjust the layer temperatures
        """

        if self.layer_count == 1:
            self.t_s_0 = new_tsno(
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
                self.t_s_0 = new_tsno(
                    self.m_s_0,
                    self.t_s_0,
                    self.cc_s_0)
                self.t_s_l = new_tsno(
                    self.m_s_l,
                    self.t_s_l,
                    self.cc_s_l)
                self.t_s = new_tsno(
                    self.m_s,
                    self.t_s,
                    self.cc_s)

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

    @property
    def dry_snow_density(self):
        """dry snow density (without H2O) at a given total snow density
        (with H2O) and snow saturation

        Args:
            rhos (float or array): total density of snow (kg/m^3)
            sat (float or array): liquid water content as volume
                ratio (h2o_sat*max_h2o_vol)

        Returns:
            float or array: dry snow density
        """
        return (self.rho - self.h2o_vol * RHO_W0) / \
            (1 - self.h2o_vol * RHO_W0 / RHO_ICE)

    @property
    def h2o_max(self):
        """Maximum water that the snowpack can hold as a specific mass given
        the depth, dry snow density and maximum h2o volume
        """
        return h2o_left(
            self.z_s,
            self.dry_snow_density,
            self.max_h2o_vol)

    @property
    def h2o_vol(self):
        return self.h2o_sat * self.max_h2o_vol

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

    def layer_mass(self):
        """
        This routine computes the specific mass for each snow layer in
        the snowcover.  A layer's mass is based its depth and the
        average snowcover density.
        """

        self.m_s_0 = self.rho * self.z_s_0
        self.m_s_l = self.rho * self.z_s_l

    def calc_layers(self):
        """
        This routine determines the # of layers in the snowcover based its
        depth and mass.  Usually, there are are 2 layers: the surface (active)
        and the lower layer.  The depth of the surface layer is set to the
        maximum depth for the surface layer (variable "max_z_s_0").  The
        remaining depth constitutes the lower layer.  The routine checks
        to see if the mass of this lower layer is above the minimum threshold
        (i.e., the mass threshold for the small run timestep).  If not,
        the surface layer is the whole snowcover, and there's no lower
        layer.

        """

        if self.m_s <= self.small_threshold:
            # less than minimum layer mass, so treat as no snowcover
            # can change to the set zeros
            z_s_0 = z_s_l = 0

        else:
            # Split up the snow depth into 1 or 2 layers
            z_s_l = self.z_s - self.max_z_s_0
            z_s_0 = self.max_z_s_0 if z_s_l > 0 else self.z_s
            z_s_l = z_s_l if z_s_l > 0 else 0

            # Make sure there's enough MASS for the lower
            # layer.  If not, then there's only 1 layer
            if z_s_l * self.rho < self.small_threshold:
                z_s_0 = self.z_s
                z_s_l = 0

        self.z_s = z_s_0 + z_s_l
        self.z_s_0 = z_s_0
        self.z_s_l = z_s_l
        self.__layer_count = False

    def check_no_layer_mass(self):
        """Reset the snowstate to zero's if the layer count is 0
        """

        if self.layer_count == 0:
            # If mass > 0, then it must be below threshold.
            # So turn this little bit of mass into water
            if self.m_s > 0.0:
                self.h2o_total += self.m_s

            self.set_zeros([
                'rho', 'm_s', 'm_s_0', 'cc_s_0', 'm_s_l',
                'cc_s_l', 'h2o', 'h2o_sat'
            ])

            # Note: Snow temperatures are set to MIN_SNOW_TEMP
            # (as degrees K) instead of 0 K to keep quantization
            # range in output image smaller.
            for col in ['t_s', 't_s_0', 't_s_l']:
                setattr(self, col, FREEZE)

    def init_layers(self):
        """If there are layers present, initialize the layers
        """

        if self.layer_count > 0:
            # Compute the specific mass and cold content for each layer
            self.layer_mass()
            self.cc_s_0 = cold_content(self.t_s_0, self.m_s_0)

            if self.layer_count == 2:
                self.cc_s_l = cold_content(self.t_s_l, self.m_s_l)
            else:
                self.t_s_l = FREEZE
                self.cc_s_l = 0

            # Compute liquid water content as volume ratio, and
            # snow density without water
            # self.h2o_vol = self.h2o_sat * self.max_h2o_vol

            # Determine the maximum liquid water content (as specific mass)
            # and the actual liquid water content (as specific mass)
            # self.h2o_max = h2o_left(
            #     self.z_s,
            #     dry_snow_density(self.rho, self.h2o_vol),
            #     self.max_h2o_vol)
            self.h2o = self.h2o_sat * self.h2o_max
