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
from pysnobal.point import InputData, Snobal
from pysnobal.pysnobal import PySnobal


class iPySnobal(PySnobal):

    CORE_CONFIG = 'ipysnobal_core_config.ini'

    def __init__(self, config_file):

        # read in the config file
        self.read_config_file(config_file)

        # parse the input parameters, from PySnobal
        self.parse_params()

        # parse the time step parameters, from PySnobal
        self.parse_time_steps()

        # open the input files
        self.read_input_data()

    def parse_elevation(self):
        # using load_dataset will load the topo into memory and close the file
        self.topo = xr.load_dataarray(self.config['topo']['filename'])
        return self.topo

    def read_input_data(self):
        """Open the input files for reading"""

        file_names = self.config['files'].copy()
        del file_names['output_file']
        file_names = list(file_names.values())

        self.dataset = xr.open_mfdataset(file_names)

        # convert to Kelvin
        self.dataset['precip_temp'] += FREEZE
        input_data.air_temp += FREEZE
        input_data.soil_temp += FREEZE

        # check the precip, temp. cannot be below freezing if rain present
        # This is only present in Snobal and not iSnobal
        mass_rain = input_data.precip_mass * (1 - input_data.percent_snow)
        input_data.loc[(mass_rain > 0.0) & (
            input_data.precip_temp < FREEZE), 'precip_temp'] = FREEZE

    def run(self):

        # loop through the input
        # do_data_tstep needs two input records
        input_data = self.input_data.to_dict('index')
        input_data1 = InputData(input_data[self.input_data.index[0]])
        del input_data[self.input_data.index[0]]

        self.snobal = Snobal(
            self.params,
            self.tstep_info,
            self.config['initial_snow_properties'],
            self.measurement_heights,
            self.output_timesteps
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
