import os
import sys
import traceback

import numpy as np
import pandas as pd

from inicheck.tools import get_user_config, check_config
from inicheck.config import UserConfig, MasterConfig
from inicheck.output import print_config_report, generate_config, print_recipe_summary

from pysnobal.snobal import snobal
from pysnobal import utils


DATA_TSTEP = 0
NORMAL_TSTEP = 1
MEDIUM_TSTEP = 2
SMALL_TSTEP = 3

DEFAULT_NORMAL_THRESHOLD = 60.0
DEFAULT_MEDIUM_THRESHOLD = 10.0
DEFAULT_sMALL_THRESHOLD = 1.0
DEFAULT_MEDIUM_TSTEP = 15.0
DEFAULT_sMALL_TSTEP = 1.0

WHOLE_TSTEP = 0x1  # output when tstep is not divided
DIVIDED_TSTEP = 0x2  # output when timestep is divided


class PySnobal():

    def __init__(self, config_file):
        """
        PySnobal is a wrapper to the Snobal C code.
        """

        # read in the config file
        self.read_config_file(config_file)

        # parse the input parameters
        self.parse_params()

        # parse the time step parameters
        self.parse_time_steps()

        # open the input files
        self.read_input_data()

        # initialize the snowpack
        self.initialize_snowpack()

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
            raise Exception('Config passed to PySnobal is neither file name nor '
                            ' UserConfig instance')

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
                    'intervals': None
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
            t = {'level': i, 'output': False, 'threshold': None,
                 'time_step': None, 'intervals': None}
            tstep_info.append(t)

        # The input data's time step must be between 1 minute and 6 hours.
        # If it is greater than 1 hour, it must be a multiple of 1 hour, e.g.
        # 2 hours, 3 hours, etc.

        data_tstep = self.config['time']['time_step']
        utils.check_range(data_tstep, 1.0, utils.hrs2min(60),
                          "input data's timestep")
        if ((data_tstep > 60) and (data_tstep % 60 != 0)):
            raise ValueError(
                "Data timestep > 60 min must be multiple of 60 min (whole hrs)")
        tstep_info[DATA_TSTEP]['time_step'] = utils.min2sec(data_tstep)

        # time step parameters
        norm_tstep = self.config['model']['norm_time_step']
        med_tstep = self.config['model']['medium_time_step']
        small_tstep = self.config['model']['small_time_step']

        # normal time step
        tstep_info[NORMAL_TSTEP].update({
            'time_step': utils.min2sec(norm_tstep),
            'intervals': int(data_tstep / norm_tstep),
            'threshold': self.config['model']['norm_threshold']
        })

        # medium time step
        tstep_info[MEDIUM_TSTEP].update({
            'time_step': utils.min2sec(med_tstep),
            'intervals': int(norm_tstep / med_tstep),
            'threshold': self.config['model']['medium_threshold']
        })

        # small time step
        tstep_info[SMALL_TSTEP].update({
            'time_step': utils.min2sec(small_tstep),
            'intervals': int(med_tstep / small_tstep),
            'threshold': self.config['model']['small_threshold']
        })

        # output
        if self.config['files']['output_mode'] == 'data':
            tstep_info[DATA_TSTEP]['output'] = DIVIDED_TSTEP
        elif self.config['files']['output_mode'] == 'normal':
            tstep_info[NORMAL_TSTEP]['output'] = WHOLE_TSTEP | DIVIDED_TSTEP
        elif self.config['files']['output_mode'] == 'all':
            tstep_info[NORMAL_TSTEP]['output'] = WHOLE_TSTEP
            tstep_info[MEDIUM_TSTEP]['output'] = WHOLE_TSTEP
            tstep_info[SMALL_TSTEP]['output'] = WHOLE_TSTEP
        else:
            tstep_info[DATA_TSTEP]['output'] = DIVIDED_TSTEP

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

        params['elevation'] = self.config['topo']['elevation']
        params['data_tstep'] = self.config['time']['time_step']
        params['max_h2o_vol'] = self.config['model']['max_h2o']
        params['max_z_s_0'] = self.config['model']['max_active']
        params['relative_heights'] = self.config['measurement_heights']['relative_heights']

        self.params = params

        # elevation to a 2D array
        self.elevation = np.atleast_2d(np.array(self.params['elevation']))

        # measurement heights
        self.measurement_heights = self.config['measurement_heights']

    def initialize_snowpack(self):
        """
        Initialize the snowpack state from the config file for the following
        snowpack state variables:

            z_s: total snowcover depth (m)
            rho: average snowcover density (kg/m^3)
            T_s_0: active snow layer temperature (C)
            T_s: average snowcover temperature (C)
            h2o_sat: % of liquid H2O saturation (relative
                    water content i.e. ratio of water in
                    snowcover to water that snowcover could
                    hold at saturation)
            z_0: roughness length (m) (for snow 0.01 to 0.0001)
        """

        self.snowpack_state = self.config['initial_snow_properties']

        # Celcuis to Kelvin and change to uppercase for snobal
        self.snowpack_state['T_s_0'] = self.snowpack_state['t_s_0'] + utils.C_TO_K
        self.snowpack_state['T_s'] = self.snowpack_state['t_s'] + utils.C_TO_K

        self.snowpack_state = self.dict2np(self.snowpack_state)

    def set_snowpack_state(self):
        """
        Potential future use to implement a setter for the snowpack state
        """
        raise NotImplementedError('set_snowpack_state not implemented yet')

    def read_input_data(self):
        """
        Read in the input data from the csv file and prepare
        for passing to Snobal. Maps the more verbose input column names
        to those that are used in Snobal
        """

        # read the input file
        input_data = pd.read_csv(
            self.config['files']['input_csv'], index_col='date_time', parse_dates=True)

        input_data.index = input_data.index.tz_localize(
            self.config['time']['time_zone'])

        # clip to the start and end date
        input_data = input_data.loc[self.start_date:self.end_date, :]

        # convert to Kelvin
        input_data.precip_temp += utils.C_TO_K
        input_data.air_temp += utils.C_TO_K
        input_data.soil_temp = + utils.C_TO_K

        # check the precip, temp. cannot be below freezing if rain present
        # This is only present in Snobal and not iSnobal
        mass_rain = input_data.precip_mass * (1 - input_data.percent_snow)
        input_data.precip_temp[(mass_rain > 0.0) & (
            input_data.precip_temp < utils.FREEZE)] = utils.FREEZE

        self.input_data = input_data

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

    def dict2np(self, d):
        """
        The at least 2d is to trick snobal into thinking it's an ndarray
        """
        return {k: np.atleast_2d(np.array(v, dtype=float)) for k, v in d.items()}

    def initialize_output(self):
        """
        Initialize the output dictionary for Snobal to write results to, the
        `output_rec` is passed to the C code for modification
        """

        # create the self.output_rec with additional fields and fill
        # There are a lot of additional terms that the original self.output_rec does not
        # have due to the output function being outside the C code which doesn't
        # have access to those variables
        sz = self.elevation.shape
        flds = ['mask', 'elevation', 'z_0', 'rho', 'T_s_0', 'T_s_l', 'T_s',
                'cc_s_0', 'cc_s_l', 'cc_s', 'm_s', 'm_s_0', 'm_s_l', 'z_s', 'z_s_0', 'z_s_l',
                'h2o_sat', 'layer_count', 'h2o', 'h2o_max', 'h2o_vol', 'h2o_total',
                'R_n_bar', 'H_bar', 'L_v_E_bar', 'G_bar', 'G_0_bar',
                'M_bar', 'delta_Q_bar', 'delta_Q_0_bar', 'E_s_sum', 'melt_sum', 'ro_pred_sum',
                'current_time', 'time_since_out']
        s = {key: np.zeros(sz) for key in flds}  # the structure fields

        # update the output rec with the initial snowpack state
        for key, val in self.snowpack_state.items():
            if key in flds:
                s[key] = val

        s['mask'] = np.ones(sz)
        self.output_rec = s

        # Create the output storage
        self.output_list = []

    def run(self):
        """
        mimic the main.c from the Snobal model
        """

        # initialize outputs
        self.initialize_output()

        # loop through the input
        # do_data_tstep needs two input records so only go
        # to the last record-1
        it = self.input_data[:-1].iterrows()
        index, input1 = next(it)    # this is the first input

        data_tstep = self.tstep_info[0]['time_step']
        timeSinceOut = 0.0
        start_step = 0  # if restart then it would be higher if this were iSnobal
        step_time = start_step * data_tstep

        self.output_rec['current_time'] = step_time * \
            np.ones(self.output_rec['elevation'].shape)
        self.output_rec['time_since_out'] = timeSinceOut * \
            np.ones(self.output_rec['elevation'].shape)

        s = snobal(self.params, self.tstep_info,
                   self.snowpack_state, self.measurement_heights)

        first_step = 1
        for index, input2 in it:

            try:
                # call do_data_tstep()
                c_snobal.do_tstep_grid(self.dict2np(input1.to_dict()), self.dict2np(
                    input2.to_dict()), self.output_rec, self.tstep_info, self.measurement_heights, self.params, first_step)

                if first_step == 1:
                    first_step = 0

                # output the results
                self.output_timestep(index)

            except Exception as e:
                traceback.print_exc()
                print('pysnobal error on time step {}'.format(index))
                print(e)
                return False
    #

            # input2 becomes input1
            input1 = input2.copy()

        # output to file
        self.output_to_file()

        return True

    def output_timestep(self, index):
        """
        Add the model outputs to the output dataframe

        **
        This is a departure from Snobal that can print out the
        sub-time steps, this will only print out on the data tstep
        (for now)
        **

        """

        # go from a numpy array to a single value
        c = {key: val[0][0] for key, val in self.output_rec.items()}

        c['date_time'] = index

        self.output_list.append(c)

        # reset the time since out
        self.output_rec['time_since_out'][0] = 0

        # # write out to a file
        # # f = self.params['out_file']
        # n = 0
        # if f is not None:

        #     curr_time_hrs = SEC_TO_HR(self.output_rec['current_time'][n])

        #     # time
        #     f.write('%g,' % curr_time_hrs)

        #     # energy budget terms
        #     f.write("%.3f,%.3f,%.3f,%.3f,%.3f,%.3f," %
        #             (self.output_rec['R_n_bar'][n], self.output_rec['H_bar'][n], self.output_rec['L_v_E_bar'][n],
        #              self.output_rec['G_bar'][n], self.output_rec['M_bar'][n], self.output_rec['delta_Q_bar'][n]))

        #     # layer terms
        #     f.write("%.3f,%.3f," %
        #             (self.output_rec['G_0_bar'][n], self.output_rec['delta_Q_0_bar'][n]))

        #     # heat storage and mass changes
        #     f.write("%.9e,%.9e,%.9e," %
        #             (self.output_rec['cc_s_0'][n], self.output_rec['cc_s_l'][n], self.output_rec['cc_s'][n]))
        #     f.write("%.8f,%.8f,%.8f," %
        #             (self.output_rec['E_s_sum'][n], self.output_rec['melt_sum'][n], self.output_rec['ro_pred_sum'][n]))

        #     #             # runoff error if data included */
        #     #             if (ro_data)
        #     #                 fprintf(out, " %.3f",
        #     #                         (ro_pred_sum - (ro * time_since_out)))

        #     # sno properties */
        #     f.write("%.6f,%.6f,%.6f,%.3f," %
        #             (self.output_rec['z_s_0'][n], self.output_rec['z_s_l'][n], self.output_rec['z_s'][n], self.output_rec['rho'][n]))
        #     f.write("%.3f,%.3f,%.3f,%.3f," %
        #             (self.output_rec['m_s_0'][n], self.output_rec['m_s_l'][n], self.output_rec['m_s'][n], self.output_rec['h2o'][n]))
        #     if self.params['temps_in_C']:
        #         f.write("%.5f,%.5f,%.5f\n" %
        #                 (K_TO_C(self.output_rec['T_s_0'][n]), K_TO_C(self.output_rec['T_s_l'][n]), K_TO_C(self.output_rec['T_s'][n])))
        #     else:
        #         f.write("%.5f,%.5f,%.5f\n" %
        #                 (self.output_rec['T_s_0'][n], self.output_rec['T_s_l'][n], self.output_rec['T_s'][n]))

        #     # reset the time since out
        #     self.output_rec['time_since_out'][n] = 0

    def output_to_file(self):
        """
        Output the model result to a file
        """

        self.output_df = pd.DataFrame(self.output_list)
        self.output_df.set_index('date_time', inplace=True)
        self.output_df.sort_index(inplace=True)

        # Kelvin to Celcius
        self.output_df['T_s_0'] -= utils.C_TO_K
        self.output_df['T_s_l'] -= utils.C_TO_K
        self.output_df['T_s'] -= utils.C_TO_K

        # remove uneeded columns
        # self.output_df.drop(
        #     ['mask', 'elevation', 'current_time',
        #         'time_since_out', 'layer_count', 'z_0',
        #         'h2o_total', 'h2o_vol'],
        #     axis=1,
        #     inplace=True)

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
