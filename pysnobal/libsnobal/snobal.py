"""
Class snobal() that will hold all the modeling components

20160109 Scott Havens
"""

import libsnobal
import numpy as np

HR_TO_SEC = 60

class snobal():
    """
    """
    
    def __init__(self, params, tstep_info, snow_prop, meas_heights):
        """
        Initialize the snobal() class with the parameters,
        time step information,
        
        Args:
            params: dictionary of parameters to run the model
            tstep_info: list of time step information
            snow_prop: the initial snow properties record
            meas_height: measurement heights
        """
        
        self.params = params
        self.tstep_info = tstep_info
        
        self.elevation = params['elevation']
        self.P_a = libsnobal.hysat(libsnobal.SEA_LEVEL, libsnobal.STD_AIRTMP, 
                                   libsnobal.STD_LAPSE, (self.elevation / 1000.0),
                                   libsnobal.GRAVITY, libsnobal.MOL_AIR)
        
        self.snow_records = snow_prop
        self.get_sn_rec(True)
        
        
    
    def get_sn_rec(self, first_rec=False):
        """
        Get the next snow property record
        
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
        
    
    def init_snow(self):
        """
        This routine initializes the properties for the snowcover.  It
        determines the number of layers, their individual properties,
        the cold content for the snowcover and its layers, and the
        snowcover's water content.
        
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
        
        
        
        
        
        