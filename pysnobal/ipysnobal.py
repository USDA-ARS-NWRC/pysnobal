import os
import sys

import numpy as np
import pandas as pd
import netCDF4 as nc
import xarray as xr
from inicheck.config import MasterConfig, UserConfig
from inicheck.output import print_config_report
from inicheck.tools import check_config, get_user_config

from pysnobal import utils
from pysnobal.core.constants import (FREEZE, MEDIUM_TSTEP, NORMAL_TSTEP,
                                     SMALL_TSTEP)
from pysnobal.spatial import iSnobal
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

        # convert to Kelvin
        self.dataset['precip_temp'] += FREEZE
        self.dataset['air_temp'] += FREEZE
        self.dataset['soil_temp'] += FREEZE

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

            self.initial_conditions = xr.merge(ic)

    def run(self):

        # loop through the input
        # do_data_tstep needs two input records
        input_data1 = self.dataset.isel(time=0)

        self.snobal = iSnobal(
            self.params,
            self.tstep_info,
            self.initial_conditions,
            self.measurement_heights
        )

        # index2 is the index for the second input timestep
        for date_time, input2 in input_data.items():

            # if index.hour == 0:
            #     print(index)

            # call do_data_tstep()
            input_data2 = InputData(input2)
            self.snobal.do_data_tstep(
                input_data1,
                input_data2
            )

            # input2 becomes input1
            input_data1 = input_data2

        # output to file
        self.output_to_file()

        return True
