# isnobal code, similar to snobal.py where it implements the model
# but in a gridded mode. Then the 2 can stay seperated or evolve
# differently as isnobal will be quite different than snobal.
import xarray as xr
import numpy as np

from pysnobal.point import Snobal
from pysnobal.core.functions import hysat
from pysnobal.core.constants import SEA_LEVEL, STD_AIRTMP, STD_LAPSE, \
    GRAVITY, MOL_AIR, NORMAL_TSTEP, MEDIUM_TSTEP, SMALL_TSTEP

from pysnobal.spatial import SpatialSnowState


class iSnobal(Snobal):

    def __init__(self, params, tstep_info, inital_conditions, meas_heights):
        """
        Initialize the iSnobal() class with the parameters,
        time step information,

        This follows the initialize() function in snobal

        Args:
            params: dictionary of parameters to run the model
            tstep_info: list of time step information
            inital_conditions: the initial snow properties
            meas_height: measurement heights
        """

        self.params = params
        self.tstep_info = tstep_info

        self.elevation = params['elevation']
        self.start_date = self.params['start_date']
        self.current_datetime = self.params['start_date']

        self.output_divided = False

        self.P_a = hysat(
            SEA_LEVEL,
            STD_AIRTMP,
            STD_LAPSE,
            self.elevation / 1000.0,
            GRAVITY,
            MOL_AIR)

        # get the intial snowcover properties
        self.snow_records = inital_conditions

        # initialize the snowcover
        init = xr.zeros_like(self.elevation)
        self.create_snow_state(SpatialSnowState, init=init)
        self.init_snow()

        # get measurement-height record
        self.measurement_heights = meas_heights
        self.relative_hts = True

        self.time_since_out = 0

        self.init_output()

        self._hle1_init = {
            'init_ustar': None,
            'init_factor': None
        }

    def divide_tstep(self):
        """
        This routine performs the model's calculations for 1 data timestep
        between 2 input-data records which are in 'input_rec1' and
        'input_rec2'.

        If there's precipitation during the data timestep, the flag
        'precip_now' used be TRUE.  Furthermore, the routine requires
        that the following precipitation variables have been initialized:

            precip_mass
            percent_snow
            rho_snow
            precip_temp

        This routine divides the data timestep into the appropriate number
        of normal run timesteps.  The input values for each normal timestep
        are computed from the two input records by linear interpolation.

        If output is desired for any of the run timesteps (normal, medium,
        or small), the appropriate output flags must be set in the proper
        timestep's record (i.e., the array 'tstep_info').  If any output
        flag is set, the routine requires that the global variable 'out_func'
        point to appropriate output function.

        This routine may return in the middle of a data timestep if:

            a)  the output function pointed to by 'out_func' is called, and
            b)  the flag 'run_no_snow' is FALSE, and
            c)  there is no snow remaining on the ground at the end of
                timestep

        In this happens, the flag 'stop_no_snow' is set to TRUE.

        """

        # Fetch the record for the timestep at the next level.
        self.next_level = self.current_level + 1
        curr_lvl_tstep = self.tstep_info[self.current_level]

        if self.input1.precip_now and curr_lvl_tstep['level'] > 1:
            self.input1.update_precip_deltas(
                self.input_deltas[curr_lvl_tstep['level']])

        # For each the new smaller timestep, either subdivide them if
        # below their mass threshold, or run the model for them.
        interval = curr_lvl_tstep['intervals']
        for i in range(interval):

            if (self.current_level != SMALL_TSTEP) and \
                    (self.below_thold(curr_lvl_tstep['threshold'])):
                # increment the level number
                self.current_level += 1
                self.divided_step = True
                if not self.divide_tstep():
                    return False
            else:
                if not self.do_tstep(curr_lvl_tstep):
                    return False

        if self.current_level == 0 and not self.output_divided:
            self.output()

        self.current_level -= 1
        self.next_level -= 1

        return True

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
        idx = self.snow_state.m_s <= self.tstep_info[SMALL_TSTEP]['threshold']
        self.snow_state.layer_count = xr.where(
            idx, 0, self.snow_state.layer_count)
        self.snow_state.z_s = xr.where(idx, 0, self.snow_state.z_s)
        self.snow_state.z_s_0 = xr.where(idx, 0, self.snow_state.z_s_0)
        self.snow_state.z_s_l = xr.where(idx, 0, self.snow_state.z_s_l)

        # not enough depth for surface layer and the lower layer,
        # so just 1 layer: surface layer
        idx = self.snow_state.z_s < self.params['max_z_s_0']
        self.snow_state.layer_count = xr.where(
            idx, 1, self.snow_state.layer_count)
        self.snow_state.z_s_0 = xr.where(
            idx, self.snow_state.z_s, self.snow_state.z_s_0)
        self.snow_state.z_s_l = xr.where(idx, 0, self.snow_state.z_s_l)
        # self.snow_state.z_s = xr.where(idx, 0, self.snow_state.z_s)

        # enough depth for both layers
        idx = self.snow_state.z_s >= self.params['max_z_s_0']
        self.snow_state.layer_count = xr.where(
            idx, 2, self.snow_state.layer_count)
        self.snow_state.z_s_0 = xr.where(
            idx, self.params['max_z_s_0'], self.snow_state.z_s_0)

        self.snow_state.z_s_l = self.snow_state.z_s - self.snow_state.z_s_0
        # z_s = z_s_0 + z_s_l  # not really needed but needed for below

        # However, make sure there's enough MASS for the lower
        # layer.  If not, then there's only 1 layer
        idx = self.snow_state.z_s_l * \
            self.snow_state.rho < self.tstep_info[SMALL_TSTEP]['threshold']
        self.snow_state.layer_count = xr.where(
            idx, 1, self.snow_state.layer_count)
        self.snow_state.z_s_0 = xr.where(
            idx, self.snow_state.z_s, self.snow_state.z_s_0)
        self.snow_state.z_s_l = xr.where(idx, 0, self.snow_state.z_s_l)

    def init_output(self):
        """Initialize the output, keeping all data in memory will not be efficient.
        Therefore, just have a Dataset that houses the current timestep that will
        be written out to disk at each timestep.

        TODO there may be a way with dask and xarray to link this to a file and
        sync the dataset with the disk
        """

        # sz = self.elevation.shape
        flds = ['rho', 't_s_0', 't_s_l', 't_s',
                'cc_s_0', 'cc_s_l', 'cc_s', 'm_s', 'm_s_0', 'm_s_l', 'z_s',
                'z_s_0', 'z_s_l', 'h2o_sat', 'layer_count', 'h2o', 'h2o_max',
                'h2o_vol', 'h2o_total', 'R_n_bar', 'H_bar', 'L_v_E_bar',
                'G_bar', 'G_0_bar', 'M_bar', 'delta_Q_bar', 'delta_Q_0_bar',
                'E_s_sum', 'melt_sum', 'swi_sum']

        self.output_rec = xr.Dataset({
            key: xr.zeros_like(self.elevation) for key in flds
        })

    def output(self):
        """
        Specify where the model output should go
        """

        c = {key: getattr(self.snow_state, key)
             for key in self.output_rec.keys()}
        c['date_time'] = self.current_datetime

        self.output_list.append(c)
        self.time_since_out = 0.0
