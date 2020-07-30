import os
import sys
import traceback
from copy import deepcopy

import numpy as np
import pandas as pd
from inicheck.config import MasterConfig, UserConfig
from inicheck.output import print_config_report
from inicheck.tools import check_config, get_user_config

from pysnobal import utils
from pysnobal.core.constants import FREEZE
from pysnobal.point import InputData, Snobal

DATA_TSTEP = 0
NORMAL_TSTEP = 1
MEDIUM_TSTEP = 2
SMALL_TSTEP = 3


class PySnobal():

    def __init__(self, config_file):
        """
        PySnobal is a wrapper to the Snobal C code.
        """

        self.output_timesteps = None

        # read in the config file
        self.read_config_file(config_file)

        # parse the input parameters
        self.parse_params()

        # parse the time step parameters
        self.parse_time_steps()

        # open the input files
        self.read_input_data()

    def read_config_file(self, config_file):
        """
        Reads in the user's config file and checks

        Args:
            config_file: either a path or a UserConfig instance
        """

        # read the config file and store
        if isinstance(config_file, str):
            if not os.path.isfile(config_file):
                raise Exception('Configuration file does not exist --> {}'
                                .format(config_file))

            # Get the master config file
            master_config = os.path.abspath(os.path.dirname(
                __file__)) + '/pysnobal_core_config.ini'
            mcfg = MasterConfig(path=master_config)

            # user config file
            ucfg = get_user_config(config_file, mcfg=mcfg)

        elif isinstance(config_file, UserConfig):
            ucfg = config_file
            config_file = config_file.filename

        else:
            raise Exception('Config passed to PySnobal is neither file '
                            'name or UserConfig instance')

        # Check the config file
        warnings, errors = check_config(ucfg)
        print_config_report(warnings, errors)
        self.ucfg = ucfg
        self.config = self.ucfg.cfg

        # Exit Pysnobal if config file has errors
        if len(errors) > 0:
            print("Errors in the config file. See configuration"
                  " status report above.")
            sys.exit()

    def parse_time_steps(self):
        """
        Parse the options dict, set the default values if not specified
        May need to divide tstep_info and params up into different
        functions. Snobal is expecting a dictionary in the following format
        for `tstep_info`:

            {
                0: {
                    'level': 0,
                    'output': {{ config.files.output_mode }},
                    'threshold': None,
                    'time_step': {{ from input data }},
                    'intervals': 1
                    }
                1: {
                    'level': 1,
                    'output': {{ config.files.output_mode }},
                    'threshold': {{ config.model.norm_threshold }},
                    'time_step': {{ config.model.norm_time_step }},
                    'intervals': {{ data_tstep / norm_tstep }}
                    }
                2: {
                    'level': 2,
                    'output': {{ config.files.output_mode }},
                    'threshold': {{ config.model.medium_threshold }},
                    'time_step': {{ config.model.medium_time_step }},
                    'intervals': {{ norm_tstep / med_tstep }}
                    }
                3: {
                    'level': 3,
                    'output': {{ config.files.output_mode }},
                    'threshold': {{ config.model.small_threshold }},
                    'time_step': {{ config.model.small_time_step }},
                    'intervals': {{ med_tstep / small_tstep }}
                    }
            }

        Where
            0 : data timestep
            1 : normal run timestep
            2 : medium run timestep
            3 : small run timestep
        """

        # intialize the time step info

        tstep_info = []
        for i in range(4):
            t = {
                'level': i,
                'output': False,
                'threshold': None,
                'time_step': None,
                'time_step_timedelta': None,
                'intervals': None,
                'intervals_per_timestep': 1
            }
            tstep_info.append(t)

        # The input data's time step must be between 1 minute and 6 hours.
        # If it is greater than 1 hour, it must be a multiple of 1 hour, e.g.
        # 2 hours, 3 hours, etc.

        data_tstep = self.config['time']['time_step']
        utils.check_range(data_tstep, 1.0, utils.hrs2min(60),
                          "input data's timestep")
        if ((data_tstep > 60) and (data_tstep % 60 != 0)):
            raise ValueError(
                "Data timestep > 60 min must be multiple of 60 min")
        tstep_info[DATA_TSTEP]['time_step'] = utils.min2sec(data_tstep)
        tstep_info[DATA_TSTEP]['time_step_timedelta'] = pd.to_timedelta(
            data_tstep, unit='min')

        # time step parameters
        norm_tstep = self.config['model']['norm_time_step']
        med_tstep = self.config['model']['medium_time_step']
        small_tstep = self.config['model']['small_time_step']

        # normal time step
        tstep_info[NORMAL_TSTEP].update({
            'time_step': utils.min2sec(norm_tstep),
            'time_step_timedelta': pd.to_timedelta(norm_tstep, unit='min'),
            'intervals': int(data_tstep / norm_tstep),
            'intervals_per_timestep': int(data_tstep / norm_tstep),
            'threshold': self.config['model']['norm_threshold']
        })

        # medium time step
        tstep_info[MEDIUM_TSTEP].update({
            'time_step': utils.min2sec(med_tstep),
            'time_step_timedelta': pd.to_timedelta(med_tstep, unit='min'),
            'intervals': int(norm_tstep / med_tstep),
            'intervals_per_timestep': int(norm_tstep / med_tstep),
            'threshold': self.config['model']['medium_threshold']
        })

        # small time step
        # Changed the interval meaning to number of intervals over
        # the whole model timestep to allow for up front deltas
        # calculation. The original Snobal would divide the intervals
        # into the medium timestep
        tstep_info[SMALL_TSTEP].update({
            'time_step': utils.min2sec(small_tstep),
            'time_step_timedelta': pd.to_timedelta(small_tstep, unit='min'),
            'intervals': int(med_tstep / small_tstep),
            'intervals_per_timestep': int(norm_tstep / small_tstep),
            'threshold': self.config['model']['small_threshold']
        })

        # output
        self.output_mode = 'normal'
        if self.config['files']['output_mode'] == 'normal':
            tstep_info[NORMAL_TSTEP]['output'] = True
        elif self.config['files']['output_mode'] == 'all':
            self.output_mode = 'all'
            tstep_info[NORMAL_TSTEP]['output'] = True
            tstep_info[MEDIUM_TSTEP]['output'] = True
            tstep_info[SMALL_TSTEP]['output'] = True

        self.tstep_info = tstep_info

    def parse_params(self):
        """
        Create a parameters dictionary to pass to Snobal that
        contain some of the run parameters specified in the
        config file.
        """

        # start and end dates for the model simulation
        self.start_date = pd.to_datetime(self.config['time']['start_date'])
        self.end_date = pd.to_datetime(self.config['time']['end_date'])

        self.start_date = self.start_date.tz_localize(
            self.config['time']['time_zone'])
        self.end_date = self.end_date.tz_localize(
            self.config['time']['time_zone'])

        # get the rest of the parameters
        params = {}
        params['start_date'] = self.start_date
        params['elevation'] = self.config['topo']['elevation']
        params['data_tstep'] = self.config['time']['time_step']
        params['max_h2o_vol'] = self.config['model']['max_h2o']
        params['max_z_s_0'] = self.config['model']['max_active']
        params['relative_heights'] = self.config['measurement_heights']['relative_heights']  # noqa

        self.params = params

        # elevation to a 2D array
        self.elevation = np.atleast_2d(np.array(self.params['elevation']))

        # measurement heights
        self.measurement_heights = self.config['measurement_heights']

    def read_input_data(self):
        """
        Read in the input data from the csv file and prepare
        for passing to Snobal. Maps the more verbose input column names
        to those that are used in Snobal
        """

        # read the input file
        input_data = pd.read_csv(
            self.config['files']['input_csv'],
            index_col='date_time',
            parse_dates=True)

        input_data.index = input_data.index.tz_localize(
            self.config['time']['time_zone'])

        # get the time delta of the data
        tdelta = input_data.index.to_series().diff().dropna().unique()[0]

        # clip to the start and end date
        # Add one time step as each Snobal timestep requires two input
        # data classes
        input_data = input_data.loc[self.start_date:self.end_date + tdelta, :]

        # convert to Kelvin
        input_data.precip_temp += FREEZE
        input_data.air_temp += FREEZE
        input_data.soil_temp += FREEZE

        # check the precip, temp. cannot be below freezing if rain present
        # This is only present in Snobal and not iSnobal
        mass_rain = input_data.precip_mass * (1 - input_data.percent_snow)
        input_data.loc[(mass_rain > 0.0) & (
            input_data.precip_temp < FREEZE), 'precip_temp'] = FREEZE

        self.input_data = input_data

        if self.output_mode == 'normal':
            self.output_timesteps = self.input_data.index

        self.map_input_data()

    def map_input_data(self):
        """
        Maps the input data from a dataframe to that needed by Snobal
        """

        mapper = {
            'net_solar': 'S_n',
            'incoming_thermal': 'I_lw',
            'air_temp': 'T_a',
            'vapor_pressure': 'e_a',
            'wind_speed': 'u',
            'soil_temp': 'T_g',
            'precip_mass': 'm_pp',
            'percent_snow': 'percent_snow',
            'rho_snow': 'rho_snow',
            'precip_temp': 'T_pp'
        }

        self.input_data = self.input_data.rename(columns=mapper)

    def run(self):
        """
        mimic the main.c from the Snobal model
        """

        # loop through the input
        # do_data_tstep needs two input records so only go
        # to the last record-1
        input_data = self.input_data[:-1].iterrows()
        index, input1 = next(input_data)    # this is the first input

        self.snobal = Snobal(
            self.params,
            self.tstep_info,
            self.config['initial_snow_properties'],
            self.measurement_heights,
            self.output_timesteps
        )

        # first_step = 1
        for index, input2 in input_data:

            # if index.hour == 0:
            #     print(index)

            try:
                # call do_data_tstep()
                self.snobal.do_data_tstep(
                    InputData(input1.to_dict()),
                    InputData(input2.to_dict())
                )

            except Exception as e:
                traceback.print_exc()
                print('pysnobal error on time step {}'.format(index))
                print(e)
                return False
    #

            # input2 becomes input1
            input1 = deepcopy(input2)

        # output to file
        self.output_to_file()

        return True

    def output_to_file(self):
        """
        Output the model result to a file
        """

        self.output_df = pd.DataFrame(self.snobal.output_list)
        self.output_df.set_index('date_time', inplace=True)
        self.output_df.sort_index(inplace=True)

        # Kelvin to Celcius
        self.output_df['T_s_0'] -= FREEZE
        self.output_df['T_s_l'] -= FREEZE
        self.output_df['T_s'] -= FREEZE

        # Map to the outputs
        mapper = {
            'R_n_bar': 'R_n',
            'H_bar': 'H',
            'L_v_E_bar': 'L_v_E',
            'G_bar': 'G',
            'G_0_bar': 'G_0',
            'M_bar': 'M',
            'delta_Q_bar': 'delta_Q',
            'delta_Q_0_bar': 'delta_Q_0',
            'E_s_sum': 'E_s',
            'melt_sum': 'melt',
            'ro_pred_sum': 'ro_predict'
        }
        self.output_df = self.output_df.rename(columns=mapper)

        keep_list = ['R_n', 'H', 'L_v_E', 'G', 'M', 'delta_Q',
                     'G_0', 'delta_Q_0', 'cc_s_0',
                     'cc_s_l', 'cc_s', 'E_s', 'melt', 'ro_predict',
                     'z_s_0', 'z_s_l', 'z_s',
                     'rho', 'm_s_0', 'm_s_l', 'm_s', 'h2o',
                     'T_s_0', 'T_s_l', 'T_s']
        self.output_df = self.output_df[keep_list]

        # Change the output format to match the original Snobal
        if self.config['files']['format_output']:

            float_map = {
                'R_n': '%.3f',
                'H': '%.3f',
                'L_v_E': '%.3f',
                'G': '%.3f',
                'M': '%.3f',
                'delta_Q': '%.3f',
                'G_0': '%.3f',
                'delta_Q_0': '%.3f',
                'cc_s_0': '%.9e',
                'cc_s_l': '%.9e',
                'cc_s': '%.9e',
                'E_s': '%.8f',
                'melt': '%.8f',
                'ro_predict': '%.8f',
                'z_s_0': '%.6f',
                'z_s_l': '%.6f',
                'z_s': '%.6f',
                'rho': '%.3f',
                'm_s_0': '%.3f',
                'm_s_l': '%.3f',
                'm_s': '%.3f',
                'h2o': '%.3f',
                'T_s_0': '%.5f',
                'T_s_l': '%.5f',
                'T_s': '%.5f'
            }
            for key, value in float_map.items():
                self.output_df[key] = self.output_df[key].map(
                    lambda x: value % x)

        self.output_df.to_csv(self.config['files']['output_csv'])
