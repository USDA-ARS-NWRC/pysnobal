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
        self.init_snow(True)

        # get measurement-height record
        self.measurement_heights = meas_heights
        self.get_measurement_height_rec(True)
        self.relative_hts = False

        # runoff data
        self.ro_data = False

        self.time_since_out = 0

        self.init_output()

        self._hle1_init = {
            'init_ustar': None,
            'init_factor': None
        }

    def create_snow_state(self):
        """Create the initial snow state
        """

        init_array = np.zeros_like(self.elevation)
        zeros = xr.DataArray(
            init_array,
            coords=self.elevation.coords,
            dims=self.elevation.dims
        )

        self.snow_state = SpatialSnowState(zeros)

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
