import xarray as xr
import numpy as np

from pysnobal.core.constants import FREEZE
from pysnobal.point.libsnobal import sati
from pysnobal.point import InputData, InputDeltas


class InputSpatialData(InputData):

    def __init__(self, data, input_delta=False, init=0):

        self.init = init * xr.ones_like(data['air_temp'])

        super(InputSpatialData, self).__init__(
            data, input_delta, self.init.copy())

        # The creation of the attributes was a little sloppy with the
        # DataArray naming, fix that here
        da_attr = []
        for attr, value in self.__dict__.items():
            if isinstance(value, xr.DataArray):
                da_attr.append(attr)

        # The init has created a reference for all data arrays, need
        # to change name and copy
        for attr in da_attr:
            self.init.name = attr
            setattr(self, attr, self.init.copy())

    @property
    def sat_vp(self):
        """Calculate the saturation vapor pressure over ice for
        the air temperature

        Returns:
            float: saturation vapor pressure over ice
        """
        if not self.__sat_vp:
            self._sat_vp = sati(self.air_temp)
            self.__sat_vp = True
        return self._sat_vp

    @property
    def soil_vp(self):
        """Calculate the saturation vapor pressure over ice for
        the soil temperature

        Returns:
            float: saturation vapor pressure over ice
        """
        if not self.__soil_vp:
            self._soil_vp = sati(self.soil_temp)
            self.__soil_vp = True
        return self._soil_vp

    def precipitation_inputs(self):

        self.precip_now = False

        if np.any(self.precip_mass.values > 0):
            self.precip_now = True

            # self.precip_mass = self.precip_mass
            # self.mass_snow = self.precip_mass * self.percent_snow
            # self.mass_rain = self.precip_mass - self.mass_snow

            self.z_snow = xr.where(
                (self.mass_snow > 0.0) & (self.rho_snow > 0.0),
                self.mass_snow / self.rho_snow,
                0
            )

            # check the precip, temp. cannot be below freezing if rain present
            self.precip_temp = xr.where(
                (self.mass_rain > 0) & (self.precip_temp < FREEZE),
                FREEZE,
                self.precip_temp
            )

            self.temp_snow = self.init.copy()
            self.temp_rain = self.init.copy()  # values are 0 or precip_temp

            # Mixed snow and rain
            idx_snow = self.mass_snow > 0
            idx_rain = self.mass_rain > 0

            idx_mixed = idx_snow & idx_rain
            self.temp_snow.values[idx_mixed] = FREEZE
            self.h2o_sat_snow.values[idx_mixed] = 1
            self.temp_rain = xr.where(
                idx_rain,
                self.precip_temp,
                self.temp_rain)

            # cold snow
            idx_freezing = self.precip_temp < FREEZE
            idx_cold_snow = idx_snow & idx_freezing
            self.temp_snow = xr.where(
                idx_cold_snow,
                self.precip_temp,
                self.temp_snow)
            self.h2o_sat_snow.values[idx_cold_snow] = 0

            # warm snow
            idx_warm_snow = idx_snow & ~idx_freezing
            self.temp_snow.values[idx_warm_snow] = FREEZE
            self.h2o_sat_snow.values[idx_warm_snow] = 1

            # # if (self.mass_snow > 0) and (self.mass_rain > 0):
            # #     self.temp_snow = FREEZE
            # #     self.h2o_sat_snow = 1
            # #     self.temp_rain = self.precip_temp
            # # elif (self.mass_snow > 0):
            #     # Snow only
            #     # if (self.precip_temp < FREEZE):
            #         # cold snow
            #         # self.temp_snow = self.precip_temp
            #         # self.h2o_sat_snow = 0
            #     # else:
            #         # warm snow
            #         self.temp_snow = FREEZE
            #         self.h2o_sat_snow = 1

            # elif (self.mass_rain > 0):
            #     # rain only
            #     self.temp_rain = self.precip_temp

    # def add_deltas(self, input_deltas):

    #     # Add the input data deltas
    #     self.net_solar = self.net_solar + input_deltas.net_solar
    #     self.thermal = self.thermal + \
    #         input_deltas.thermal
    #     self.air_temp = self.air_temp + input_deltas.air_temp
    #     self.vapor_pressure = self.vapor_pressure + input_deltas.vapor_pressure
    #     self.wind_speed = self.wind_speed + input_deltas.wind_speed
    #     self.soil_temp = self.soil_temp + input_deltas.soil_temp

    #     self.update_precip_deltas(input_deltas)

    # def update_precip_deltas(self, input_deltas):

    #     # update the precipitation. Snobal takes the input deltas
    #     # and divides by the intervals
    #     for precip_variable in self.PRECIP_VARIABLES:
    #         setattr(self, precip_variable,
    #                 getattr(input_deltas, precip_variable))


class InputSpatialDeltas(InputDeltas):

    def __init__(self, input1, input2, tstep_info):

        super(InputSpatialDeltas, self).__init__(
            input1, input2, tstep_info)

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

            self.deltas[tstep['level']] = InputSpatialData(tstep_deltas)

        return self.deltas
