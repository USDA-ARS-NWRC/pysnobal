import os
import sys

import numpy as np
import pandas as pd
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

    def read_input_data(self):
        pass
