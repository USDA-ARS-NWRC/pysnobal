import xarray as xr

from pysnobal.core.constants import CAL_TO_J, FREEZE, RHO_ICE, RHO_W0
from pysnobal.core.functions import cp_ice, time_average
from pysnobal.point.libsnobal import sati
from pysnobal.point import SnowState


class SpatialSnowState(SnowState):

    def __init__(self, init):

        super(SpatialSnowState, self).__init__(init)

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
            self.T_s_0 = self.new_tsno(
                self.m_s_0,
                self.T_s_0,
                self.cc_s_0)
            self.T_s = self.T_s_0

        elif self.layer_count == 2:
            if self.isothermal:
                self.T_s = FREEZE
                self.T_s_l = FREEZE
                self.T_s_0 = FREEZE
            else:
                self.T_s_0 = self.new_tsno(
                    self.m_s_0,
                    self.T_s_0,
                    self.cc_s_0)
                self.T_s_l = self.new_tsno(
                    self.m_s_l,
                    self.T_s_l,
                    self.cc_s_l)
                self.T_s = self.new_tsno(
                    self.m_s,
                    self.T_s,
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
        self.ro_pred_sum = self.ro_predict

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
        self.ro_pred_sum += self.ro_predict
