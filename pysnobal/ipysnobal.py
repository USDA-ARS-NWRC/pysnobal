import os
import sys

import numpy as np
import pandas as pd
import xarray as xr
from inicheck.config import MasterConfig, UserConfig
from inicheck.output import print_config_report
from inicheck.tools import check_config, get_user_config

from pysnobal import utils
from pysnobal.core.constants import (FREEZE, MEDIUM_TSTEP, NORMAL_TSTEP,
                                     SMALL_TSTEP)
from pysnobal.spatial import iSnobal, InputSpatialData, InputSpatialDeltas
from pysnobal.pysnobal import PySnobal


class iPySnobal(PySnobal):

    CORE_CONFIG = 'ipysnobal_core_config.ini'

    # initial conditions, {name: unit}
    INITIAL_CONDITIONS = {
        'z_0': 'm',
        'z_s': 'm',
        'rho': 'kg/m^3',
        't_s_0': 'C',
        't_s': 'C',
        'h2o_sat': '%'
    }

    def __init__(self, config_file):

        # read in the config file
        self.read_config_file(config_file)
        self.parse_params()
        self.parse_time_steps()
        self.read_input_data()
        self.read_initial_conditions()

    def parse_elevation(self):
        # using load_dataset will load the topo into memory and close the file
        # this will most likely not work with the current topo.nc files which
        # would be a dataset
        self.topo = xr.load_dataarray(self.config['topo']['filename'])
        return self.topo

    def read_input_data(self):
        """Open the input files for reading"""

        file_names = self.config['files'].copy()
        del file_names['output_file']
        del file_names['initial_conditions']
        file_names = list(file_names.values())

        self.dataset = xr.open_mfdataset(file_names)

        # Assign the time zone, NOTE xarray appears to convert these to
        # UTC when displaying the times
        time = self.dataset['time'].to_index()
        time_tz = time.tz_localize(self.config['time']['time_zone'])
        self.dataset = self.dataset.assign_coords(time=time_tz.to_series())

        # convert to Kelvin
        self.dataset['precip_temp'] += FREEZE
        self.dataset['air_temp'] += FREEZE
        self.dataset['soil_temp'] += FREEZE

        self.dataset = self.dataset.rename_vars({
            'precip': 'precip_mass',
            'snow_density': 'rho_snow'
        })

        # get the time delta of the data
        tdelta = np.unique(
            self.dataset.time.values[1:] - self.dataset.time.values[:-1])[0]

        # clip to the start and end date
        # Add one time step as each Snobal timestep requires two input
        # data classes
        input_data = self.dataset.sel(time=slice(
            self.start_date, self.end_date + tdelta))

    def read_initial_conditions(self):

        if self.config['files']['initial_conditions'] is not None:
            self.initial_conditions = xr.open_dataset(
                self.config['files']['initial_conditions'])

        else:
            # use the initial snowpack conditions section to fill in
            init_array = np.ones_like(self.topo)

            ic = []
            for name, unit in self.INITIAL_CONDITIONS.items():
                constant = self.config['initial_snow_properties'][name.lower()]
                z = self.elevation.copy(deep=True, data=constant * init_array)
                z.attrs['units'] = unit
                z = z.rename(name)
                ic.append(z)

            # ic = xr.merge(ic).to_dict()
            # self.initial_conditions = {key: np.array(data['data'])
            #                            for key, data in ic['data_vars'].items()}

            self.initial_conditions = xr.merge(ic)

    def run(self):

        self.init_output()

        # loop through the input
        # do_data_tstep needs two input records
        input_data1 = InputSpatialData(self.dataset.isel(time=0))

        self.snobal = iSnobal(
            self.params,
            self.tstep_info,
            self.initial_conditions,
            self.measurement_heights
        )

        # index2 is the index for the second input timestep
        for idx, date_time in enumerate(self.dataset.time[1:]):

            # if pd.to_datetime(date_time.values).hour == 0:
            print(date_time.values)

            # call do_data_tstep()
            input_data2 = InputSpatialData(self.dataset.isel(time=idx))

            # skip doing anything if there isn't any snow and no
            # precip in the inputs, could expand this further to lazily load
            # the rest of the data by just loading the precip first
            if not self.snobal.any_snowcover and not input_data1.precip_now:
                self.snobal.proceed_no_snow()

            else:

                input_deltas = InputSpatialDeltas(
                    input_data1,
                    input_data2,
                    self.snobal.tstep_info).calculate()

                self.snobal.do_data_tstep(
                    input_data1,
                    input_data2,
                    input_deltas
                )

            # input2 becomes input1
            input_data1 = input_data2

            # output to file
            self.output_to_file()

        return True

    def init_output(self):
        self.output_file = None

    def output_to_file(self):
        """
        Output the model result to a file
        """
        pass
        # # plumb this later, just keep in a dataset
        # if self.output_file is None:
        #     self.output_file = self.snobal.output_rec.copy()
        # else:
        #     self.output_file = xr.merge(
        #         [self.output_file, self.snobal.output_rec])
