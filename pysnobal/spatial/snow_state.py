import xarray as xr
import numpy as np

from pysnobal.core.constants import CAL_TO_J, FREEZE, RHO_ICE, RHO_W0
from pysnobal.core.functions import cp_ice, time_average
from pysnobal.point.libsnobal import sati
from pysnobal.point import SnowState


class SpatialSnowState(SnowState):

    def __init__(self, init=0, max_z_s_0=0.25, small_threshold=1):

        super(SpatialSnowState, self).__init__(
            init, max_z_s_0, small_threshold)

        self.num_grid = init.shape[0] * init.shape[1]

    def set_zeros(self, fields, idx):

        if isinstance(fields, str):
            fields = fields.split()

        for field in fields:
            setattr(self, field, self.zeros)

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

        # less than minimum layer mass, so treat as no snowcover
        # can change to the set zeros
        idx0 = self.m_s <= self.small_threshold
        self.z_s_0[idx0] = 0
        self.z_s_l[idx0] = 0

        num_vals_changed = np.sum(idx0.values)

        if num_vals_changed != self.num_grid:

            # Split up the snow depth into 1 or 2 layers
            z_s_l = self.z_s - self.max_z_s_0
            z_s_0 = self.max_z_s_0 if z_s_l > 0 else self.z_s
            z_s_l = z_s_l if z_s_l > 0 else 0

            # Make sure there's enough MASS for the lower
            # layer.  If not, then there's only 1 layer
            if z_s_l * self.rho < self.small_threshold:
                z_s_0 = self.z_s
                z_s_l = 0

        self.z_s = self.z_s_0 + self.z_s_l
        # self.z_s_0 = z_s_0
        # self.z_s_l = z_s_l
        self.__layer_count = False

    def check_no_layer_mass(self):
        """Reset the snowstate to zero's if the layer count is 0
        """

        idx = self.layer_count == 0

        # If mass > 0, then it must be below threshold.
        # So turn this little bit of mass into water
        idx_mass = np.where(idx & self.m_s > 0.0)
        if np.any(idx_mass):
            self.h2o_total[idx_mass] += self.m_s[idx_mass]

        self.set_zeros([
            'rho', 'm_s', 'm_s_0', 'cc_s_0', 'm_s_l',
            'cc_s_l', 'h2o', 'h2o_sat'
        ], idx)

        # Note: Snow temperatures are set to MIN_SNOW_TEMP
        # (as degrees K) instead of 0 K to keep quantization
        # range in output image smaller.
        self.t_s[idx] = FREEZE
        self.t_s_0[idx] = FREEZE
        self.t_s_l[idx] = FREEZE
