"""
Class snobal() that will hold all the modeling components

20160109 Scott Havens
"""

import libsnobal
import numpy as np
import pandas as pd
import warnings

# Some constants and equations
HR_TO_SEC = 60
DATA_TSTEP = 0
NORMAL_TSTEP = 1
MEDIUM_TSTEP = 2
SMALL_TSTEP = 3

MIN_SNOW_TEMP = -75
FREEZE = 273.16

# density of water at 0C (kg/m^3) (from CRC handbook pg F-11)
RHO_W0 = 999.87

# density of ice - no air (kg/m^3) (from CRC handbook pg F-1)
RHO_ICE = 917.0

# specific heat of ice (J/(kg K)) (from CRC table D-159; most accurate from 0 to -10 C)
# t - temperature in K
CP_ICE = lambda t: ( 0.238846 * (0.024928 + (0.00176*(t))) / 0.001 )


# 'dry' snow density (without H2O) at a given total snow density (with H2O) and snow saturation
# rhos = total density of snow (kg/m^3)
# sat  = snow saturation (see SNO_SAT)
DRY_SNO_RHO = lambda (rhos,sat): ( ((rhos) - (sat) * RHO_W0) / (1 - (sat) * RHO_W0 / RHO_ICE) )

# water retained by snow at given saturation (see SNO_SAT)
#    d    = total depth of snow (m)
#    rhos = density of snow (kg/m^3)
#    sat  = snow saturation (see SNO_SAT)
H2O_LEFT = lambda (d,rhos,sat): ( (sat * d * RHO_W0 * (RHO_ICE - rhos)) / RHO_ICE )

class snobal():
    """
    
    To-do
    self.snow is the working data frame only containing one record
    self.output will be a dataframe that self.snow gets added 
    
    """
    
    def __init__(self, params, tstep_info, snow_prop, meas_heights, precip):
        """
        Initialize the snobal() class with the parameters,
        time step information,
        
        This follows the initialize() function in snobal
        
        Args:
            params: dictionary of parameters to run the model
            tstep_info: list of time step information
            snow_prop: the initial snow properties record
            meas_height: measurement heights
            precip: precipitation record
        """
        
        self.params = params
        self.tstep_info = tstep_info
        
        self.elevation = params['elevation']
        self.P_a = libsnobal.hysat(libsnobal.SEA_LEVEL, libsnobal.STD_AIRTMP, 
                                   libsnobal.STD_LAPSE, (self.elevation / 1000.0),
                                   libsnobal.GRAVITY, libsnobal.MOL_AIR)
        
        # get the intial snowcover properties
        self.snow_records = snow_prop
        self.get_sn_rec(True)
        
        # initialize the snowcover
        self.init_snow()
        
        # get measurement-height record
        self.mh_prop = meas_heights
        self.get_mh_rec(True)
        self.relative_hts = False
        
        # runoff data
        self.ro_data = False
        
        # get precipitation record
        self.precip_record = precip
        self.get_pr_rec(True)
        
        self.time_since_out = 0
        
        
    def do_data_tstep(self, input1, input2):
        """
        This routine performs the model's calculations for 1 data timestep
        between 2 input-data records which are in 'input_rec1' and 
        'input_rec2'.
        
        If there's precipitation during the data timestep, the flag
        'precip_now' used be TRUE.  Furthermore, the routine requires
        that the following precipitation variables have been initialized:
        
            m_pp
            percent_snow
            rho_snow
            T_pp
        
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
        
        Args:
            input1: first timestep
            input2: second timestep
            
            inputs contain all forcing data:
                ['S_n', 'I_lw', 'T_a', 'e_a', 'u', 'T_g','m_pp',
                    'percent_snow', 'rho_snow', 'T_pp']
                
            
        """
        
        # Compute deltas for the climate input parameters over the data timestep.
        input_deltas = []
        input_deltas[DATA_TSTEP] = input2 - input1
        
        # If there is precipitation, then compute the amount of rain & snow in it.
        # Look at the first input record
        pp_info = pd.Series(index=['m_pp','m_snow','m_rain','z_snow'])
        if input1.m_pp > 0:
            pp_info
        
    
    def get_sn_rec(self, first_rec=False):
        """
        This routine loads the next snow-properties record into the
        proper snow variables.  Before loading the next record though,
        it computes the difference between the current snow properties
        (predicted) and those in the next record (measured).  It then
        reads next record from either the corresponding input file or
        standard input.  If there are no more records are available, the
        global variable "more_sn_recs" is set to FALSE.
        
        Args:
            first_rec: whether or not it's the first record
        """
        
        if first_rec:
            
            self.snow_prop_index = 0    # keep track of which snow property to read
            
            self.time_s = self.snow_records.iloc[self.snow_prop_index].name * HR_TO_SEC
            self.curr_time_hrs = self.snow_records.iloc[self.snow_prop_index].name * HR_TO_SEC
            self.start_time = self.snow_records.iloc[self.snow_prop_index].name * HR_TO_SEC
            self.more_sn_recs = True
        
        else:
            # haven't ever seen this used so not entirly sure now this works
            # increase the index if there are more than one record
            self.snow_prop_index += 1
            
        self.more_sn_recs = False
        
    
    def get_mh_rec(self, first_rec=False):
        """
        This routine loads the next measurement-heights record into
        the proper mh variables.  It then reads the next record from
        either the corresponding input file or standard input.  If
        there are no more records are available, the global variable
        "more_mh_recs" is set to FALSE.
        
        Args:
            first_rec: whether or not it's the first record
        """
        
        if first_rec:
            self.mh_prop_index = 0    # keep track of which mh property to read
            
            self.time_z = self.mh_prop.iloc[self.mh_prop_index].name * HR_TO_SEC
            self.z_u = self.mh_prop.iloc[self.mh_prop_index].z_u
            self.z_T = self.mh_prop.iloc[self.mh_prop_index].z_T
            self.z_0 = self.mh_prop.iloc[self.mh_prop_index].z_0
            self.z_g = self.mh_prop.iloc[self.mh_prop_index].z_g
            
        else:
            self.mh_prop_index += 1
            
        self.more_mh_recs = False
    
    
    def get_pr_rec(self, first_rec=False):
        """
        This routine reads the next precipitation record into the
        proper precip variables.  This record is read from either
        the corresponding input file or standard input.  If there
        are no more records are available, the global variable
        "more_pr_recs" is set to FALSE.
        
        The precip record file should have one record for each time
        step
        
        Args:
            first_rec: whether or not it's the first record
            
        """
        
        if first_rec:
            self.pr_prop_index = 0    # keep track of which precip property to read
            
        else:
            self.pr_prop_index += 1
            
        self.time_pp = self.precip_record.iloc[self.pr_prop_index].name * HR_TO_SEC
        self.m_pp = self.precip_record.iloc[self.pr_prop_index].m_pp
        self.percent_snow = self.precip_record.iloc[self.pr_prop_index].percent_snow
        self.rho_snow = self.precip_record.iloc[self.pr_prop_index].rho_snow
        self.T_pp = self.precip_record.iloc[self.pr_prop_index].T_pp
            
        # Calculate properties for rain and snow portion of precipitation.
        mass_rain = self.m_pp - self.percent_snow * self.m_pp

        # Precip. temp. cannot be below freezing if rain present.
        if ((mass_rain > 0.0) and (self.T_pp < FREEZE)):
            if (self.T_pp < (FREEZE - 0.5)):
                warnings.warn("T_pp < 0.0 C with rain; setting T_pp to 0.0 C")
            self.T_pp = FREEZE
        
        
        
    
    def init_snow(self):
        """
        This routine initializes the properties for the snowcover.  It
        determines the number of layers, their individual properties,
        the cold content for the snowcover and its layers, and the
        snowcover's water content.
        
        Initialize all the following values
        h2o_sat, layer_count, m_s_0, m_s_l, max_h2o_vol, rho, T_s, T_s_0, T_s_l,
        z_s, cc_s, cc_s_0, cc_s_l, h2o, h2o_max, h2o_total, h2o_vol, m_s, m_s_0, m_s_l
        
        Args:
            z_s: depth of snowcover (m)
            rho: density of snowcover (kg/m^3)
            T_s: average temperature of snowcover (K)
            T_s_0: temperature of surface layer of snowcover (K)
            T_s_l: temperature of lower layer of snowcover (K)
            h2o_sat:  % of liquid h2o saturation (0 to 1.0)
            max_h2o_vol: maximum liquid h2o content as volume ratio:
                        V_water/(V_snow - V_ice) (unitless)
                        
        """
        
        cols = ['h2o_sat', 'layer_count', 'm_s_0', 'm_s_l', 'max_h2o_vol', 'rho', 'T_s', 
                'T_s_0', 'T_s_l', 'z_s', 'cc_s', 'cc_s_0', 'cc_s_l', 'h2o', 'h2o_max', 
                'h2o_total', 'h2o_vol', 'm_s', 'm_s_0', 'm_s_l']
        
        # create an empty dataframe
        self.snow = pd.Series(index=cols)
        
        # set the values from the initial snow properties
        self.snow.z_s = self.snow_records.iloc[self.snow_prop_index].z_s
        self.snow.rho = self.snow_records.iloc[self.snow_prop_index].rho
        self.snow.T_s_0 = self.snow_records.iloc[self.snow_prop_index].T_s_0
        self.snow.T_s = self.snow_records.iloc[self.snow_prop_index].T_s
        self.snow.h2o_sat = self.snow_records.iloc[self.snow_prop_index].h2o_sat
        self.snow.max_h2o_vol = self.params['max_h2o_vol']
        
        # initialize the snowpack
        self.snow.m_s = self.snow.rho * self.snow.z_s
        
        # determine the number of layers
        self.calc_layers()
        
        if self.snow.layer_count == 0:
            # If mass > 0, then it must be below threshold.
            # So turn this little bit of mass into water
            
            if self.snow.m_s > 0.0:
                self.snow.h2o_total += self.snow.m_s
            
            for col in ['rho','m_s','cc_s','m_s_0','cc_s_0','m_s_l','cc_s_l','h2o_vol','h2o','h2o_max','h2o_sat']:
                self.snow.loc[col] = 0        
        
            # Note: Snow temperatures are set to MIN_SNOW_TEMP
            # (as degrees K) instead of 0 K to keep quantization
            # range in output image smaller.
            for col in ['T_s', 'T_s_0', 'T_s_l']:
                self.snow.loc[col] = MIN_SNOW_TEMP + FREEZE;
         
        else:
            # Compute the specific mass and cold content for each layer   
            self.layer_mass()
            self.snow.cc_s_0 = self.cold_content(self.snow.T_s_0, self.snow.m_s_0)
            
            if self.snow.layer_count == 2:
                self.snow.cc_s_l = self.cold_content(self.snow.T_s_l, self.snow.m_s_l)
            else:
                self.snow.T_s_l = MIN_SNOW_TEMP + FREEZE
                self.snow.cc_s_l = 0
            
            # Compute liquid water content as volume ratio, and
            # snow density without water
            self.snow.h2o_vol = self.snow.h2o_sat * self.snow.max_h2o_vol
            rho_dry = DRY_SNO_RHO(self.snow.rho, self.snow.h2o_vol)
            
            # Determine the maximum liquid water content (as specific mass)
            # and the actual liquid water content (as specific mass)
            self.snow.h2o_max = H2O_LEFT(self.snow.z_s, rho_dry, self.snow.max_h2o_vol)
            self.snow.h2o = self.snow.h2o_sat * self.snow.h2o_max


        
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
        
        if self.snow.m_s <= self.tstep_info[SMALL_TSTEP]['threshold']:
            # less than minimum layer mass, so treat as no snowcover
            
            layer_count = 0
            z_s = z_s_0 = z_s_l = 0
            
        elif self.snow.z_s < self.params['max_z_s_0']:
            # not enough depth for surface layer and the lower layer,
            # so just 1 layer: surface layer
            
            layer_count = 1
            z_s_0 = self.snow.z_s
            z_s_l = 0
            
        else:
            # enough depth for both layers
            
            layer_count = 2
            z_s_0 = self.params['max_z_s_0']
            z_s_l = self.snow.z_s - z_s_0
            
            # However, make sure there's enough MASS for the lower
            # layer.  If not, then there's only 1 layer
            if z_s_l * self.snow.rho < self.tstep_info[SMALL_TSTEP]['threshold']:
                layer_count = 1
                z_s_0 = self.snow.z_s
                z_s_l = 0
        
            
        self.snow.layer_count = layer_count
        self.snow.z_s = z_s
        self.snow.z_s_0 = z_s_0
        self.snow.z_s_l = z_s_l
        
        
    def layer_mass(self):
        """
        This routine computes the specific mass for each snow layer in
        the snowcover.  A layer's mass is based its depth and the
        average snowcover density.
        """
        
        if self.snow.layer_count == 0:
            self.snow.m_s_0 = 0
            self.snow.m_s_l = 0
            
        else:
            # layer count is 1 or 2
            self.snow.m_s_0 = self.snow.rho * self.snow.z_s_0
            
            if self.snow.layer_count == 2:
                self.snow.m_s_l = self.snow.rho * self.snow.z_s_l
            else:
                self.snow.m_s_l = 0
                
        
    def cold_content(self, temp, mass):
        """
        This routine calculates the cold content for a layer (i.e., the
        energy required to bring its temperature to freezing) from the
        layer's temperature and specific mass.
        
        Args:
            temp: temperature of layer 
            mass: specific mass of layer
            
        Returns:
            cc: cold content of layer
        """
        
        cc = 0
        if temp < FREEZE:
            self.heat_stor(CP_ICE(temp), mass, temp-FREEZE)
        return cc
    

    def heat_stor(self, cp, spm, tdif):
        """
        Calculate the heat storage
        Args:
            cp: specific heat of layer (J/kg K)
            spm: layer specific mass (kg/m^2)
            tdif: temperature change (K)
        """
        
        return cp * spm * tdif

        