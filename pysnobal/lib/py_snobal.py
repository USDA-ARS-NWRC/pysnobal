"""
Class snobal() that will hold all the modeling components

20160109 Scott Havens
"""

import py_libsnobal as libsnobal
import core_c

# import libsnobal
import numpy as np
# import numpy.ma as ma
# import pandas as pd
import warnings
from copy import copy, deepcopy
import matplotlib.pyplot as plt


# Some constants and equations
WHOLE_TSTEP = 0x1 # output when tstep is not divided
DIVIDED_TSTEP = 0x2  # output when timestep is divided

HR_TO_SEC = 3600.0
SEC_TO_HR = lambda x: x / 3600.0
DATA_TSTEP = 0
NORMAL_TSTEP = 1
MEDIUM_TSTEP = 2
SMALL_TSTEP = 3

MIN_SNOW_TEMP = -75
FREEZE = 273.16
KT_MOISTSAND = 1.65
MAX_SNOW_DENSITY = 600
MAX_LAYERS = 2

# density of water at 0C (kg/m^3) (from CRC handbook pg F-11)
RHO_W0 = 999.87

# specific heat of water at 0C (J / (kg K))
CP_W0 = 4217.7

# density of ice - no air (kg/m^3) (from CRC handbook pg F-1)
RHO_ICE = 917.0

# ratio vaporization to sublimation
VAP_SUB  = 2.501 / 2.835 

# Kelvin to Celcius
K_TO_C = lambda x: x - FREEZE

# specific heat of ice (J/(kg K)) (from CRC table D-159; most accurate from 0 to -10 C)
# t - temperature in K
CP_ICE = lambda t: CAL_TO_J * (0.024928 + (0.00176 * t)) / 0.001 

# specific heat of water (J/(kg K))
#    (from CRC table D-158; most accurate from 0 to +10 C)
#    (incorrect at temperatures above 25 C)
CP_WATER = lambda t: CP_W0 - 2.55 * (t - FREEZE)


# 'dry' snow density (without H2O) at a given total snow density (with H2O) and snow saturation
# rhos = total density of snow (kg/m^3)
# sat  = snow saturation (see SNO_SAT)
DRY_SNO_RHO = lambda rhos,sat: ( ((rhos) - (sat) * RHO_W0) / (1 - (sat) * RHO_W0 / RHO_ICE) )

# water retained by snow at given saturation (see SNO_SAT)
#    d    = total depth of snow (m)
#    rhos = density of snow (kg/m^3)
#    sat  = snow saturation (see SNO_SAT)
H2O_LEFT = lambda d,rhos,sat: ( (sat * d * RHO_W0 * (RHO_ICE - rhos)) / RHO_ICE )


# Convert calories to Joules
CAL_TO_J = 4.186798188


# melt (kg/m^2), Q = available energy (J/m^2)
MELT = lambda Q: Q / libsnobal.LH_FUS(FREEZE)

SNOW_EMISSIVITY = 0.98
STEF_BOLTZ = 5.67032e-8     # Stefan-Boltzmann constant (W / m^2 / deg^4)


def TIME_AVG(avg,total_time,value,time_incr):
    """
    A macro to update a time-weighted average for a quantity.
    avg        current average
    total_time    the time interval the current average applies to
    value        new value to be averaged in
    time_incr    the time interval the new value applies to
    TIME_AVG = lambda avg,total_time,value,time_incr: (avg * total_time + value * time_incr) / (total_time + time_incr)
    """
    return (avg * total_time + value * time_incr) / (total_time + time_incr)

def KTS(rho):
    """
    thermal conductivity of snow (J/(m sec K))
    (after Yen, 1965, see Anderson, 1976, pg. 31)
    rho = snow density (kg/m^3)
    """
#     return CAL_TO_J * 0.0077 * (rho/1000.0) * (rho/1000.0)
    return CAL_TO_J * 0.0077 * np.power(rho/1000.0, 2)

class snobal(object):
    """
    
    To-do
    self.snow is the working data frame only containing one record
    self.output will be a dataframe that self.snow gets added 
    
    """
    
    # came from self.__dict__.keys()
    # slots will help with memory when multiple instances of snobal are used
#     __slots__ = ['em', 'input2', 'input1', 'z_t', 'computed', 'precip_now', 
#                  'snow_records', 'time_step', 'precip_info', 'tstep_level', 
#                  'current_time', 'z_u', 'time_since_out', 'time_s', 'snow_prop_index', 
#                  'input_deltas', 'ro_data', 'snow', 'time_z', 'P_a', 'params', 'z_g', 
#                  'next_level', 'elevation', 'mh_prop_index', 'start_time', 'tstep_info', 
#                  'curr_time_hrs', 'mh_prop', 'z_0', 'more_sn_recs', 'relative_hts', 
#                  'curr_level', 'precip', 'more_mh_recs', 'snowcover', 'isothermal']
    
    
    def __init__(self, params, tstep_info, snow_prop, meas_heights):
        """
        Initialize the snobal() class with the parameters,
        time step information,
        
        This follows the initialize() function in snobal
        
        Args:
            params: dictionary of parameters to run the model
            tstep_info: list of time step information
            snow_prop: the initial snow properties record
                required field: time_s, z, z_s, rho, T_s, T_s_0, h2o_sat
            meas_height: measurement heights
        """
        
        self.params = params
        self.tstep_info = tstep_info
        
#         self.elevation = ma.masked_array(snow_prop['z'], mask=self.mask)
        self.elevation = snow_prop['z']
        
        self.gridded_P_a = libsnobal.hysat(libsnobal.SEA_LEVEL, libsnobal.STD_AIRTMP, 
                                   libsnobal.STD_LAPSE, self.elevation / 1000.0,
                                   libsnobal.GRAVITY, libsnobal.MOL_AIR)
        
        self.shape = self.elevation.shape
        self.ngrid = np.prod(self.shape)
#         self.zeros = np.zeros(self.shape)

        # check for a mask
        if 'mask' in snow_prop.keys():
            self.mask = snow_prop['mask'].astype(np.bool)
        else:
            self.mask = np.ones(self.shape, dtype=np.bool)
        
        # get the intial snowcover properties
        self.snow_records = snow_prop
        self.get_sn_rec(True)
        
        # intialize em/snow that will be the size of the grid, but will vary based
        # on what time step is being processed
        # initialize the em
        self.init_em()
        
        # initialize the snowcover
        self.init_snow(True)
        
        # now copy these into a new instance that will store the full model states
        # so that self.snow/self.em can be used to work with the data time step
        # subsets
        self.gridded_snow = deepcopy(self.snow)
        self.gridded_em = deepcopy(self.em)
        
        # initialize the 
                
        # initialize the precip for the time step
#         self.init_precip()
        
        # get measurement-height record
        self.mh_prop = meas_heights
        self.get_mh_rec(True)
        self.relative_hts = params['relative_heights']
        
        # runoff data
        self.ro_data = False
        
#         # get precipitation record
#         self.precip_record = precip
#         self.get_pr_rec(True)
        
        # some other variables
        self.input_deltas = dict()
        self.computed = dict()
        self.precip_info = dict()
        for level in range(DATA_TSTEP, SMALL_TSTEP+1):
            self.input_deltas[level] = {}
            self.computed[level] = None
            self.precip_info[level] = {}
        
        
        self.time_since_out = np.zeros(self.shape)
        
        # set some attributes
        self.isothermal = np.zeros(self.shape, dtype=bool)
        
        
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
            input1: first timestep dict
            input2: second timestep sixr
            
            inputs contain all forcing data:
                ['S_n', 'I_lw', 'T_a', 'e_a', 'u', 'T_g','m_pp',
                    'percent_snow', 'rho_snow', 'T_pp']
                
            
        """
        
#         print '%.2f' % (self.current_time/3600.0)

        if np.max(self.current_time)/3600.0 > 953:
            self.curr_level
        
        # store the inputs for later
        self.input1 = {}
        self.input2 = {}
        self.gridded_input1 = input1.copy()
        self.gridded_input2 = input2.copy()
        
        # extract_data.c performs an init_snow() which resets some of the variables
        # this is an ugly copy of that which can be improved on later
#         self.snow = deepcopy(self.gridded_snow)
#         self.em = deepcopy(self.gridded_em)
#         self.init_snow()
#         self.gridded_snow = deepcopy(self.snow)
#         self.gridded_em = deepcopy(self.em)
        self.init_water_content_gridded()
        
        # Compute deltas for the climate input parameters over the data timestep.
#         for i in input1.keys():
#             self.input_deltas[DATA_TSTEP][i] = input2[i] - input1[i] 
#         self.input_deltas[DATA_TSTEP] = input2.subtract(input1)
        
        # If there is precipitation, then compute the amount of rain & snow in it.
        # Look at the first input record
#         self.precip = self.init_precip()
#         self.precip_now = np.zeros_like(self.elevation, dtype=bool)
        
        # this precip will hold the temperatures and such, set to zero initially
#         self.precip = {key: 0.0 for key in ['T_rain','T_snow','h2o_sat_snow']}
#         self.precip = EmptyClass(['T_rain','T_snow','h2o_sat_snow'], self.shape)
        self.precip = self.init_precip()
        
#         # fill precip info with empty classes
#         for level in range(DATA_TSTEP, SMALL_TSTEP+1):
#             self.precip_info[level] = self.init_precip()
        
        # check to see if there is precipitation
        precip = input1['m_pp'] > 0
        if np.any(precip):
            
            self.precip.now = precip
            self.precip.m_pp = input1['m_pp']
            self.precip.m_snow = input1['m_pp'] * input1['percent_snow']
            self.precip.m_rain = input1['m_pp'] - self.precip.m_snow
            
            ind = self.precip.m_snow > 0
            if np.any(ind):
                if np.any(input1['rho_snow'][ind] > 0):
                    self.precip.z_snow[ind] = self.precip.m_snow[ind] / input1['rho_snow'][ind]
                else:
                    raise ValueError('input1["rho_snow"] is <= 0.0 with input1["percent_snow"] > 0.0')
             
            # partion the precip based on the pixels m_snow and m_rain   
            
            rain = self.precip.m_rain > 0 
            snow = self.precip.m_snow > 0
            
            snow_only = snow & ~rain
            rain_only = ~snow & rain
            mixed = snow & rain
            
            # Mixed snow and rain
            if np.any(mixed):
                self.precip.T_snow[mixed] = FREEZE
                self.precip.h2o_sat_snow[mixed] = 1.0
                self.precip.T_rain[mixed] = input1['T_pp'][mixed]
            
            
            # snow only
            if np.any(snow_only):
                self.precip.T_snow[snow_only] = input1['T_pp'][snow_only]
                self.precip.h2o_sat_snow[snow_only] = 0.0
                
                # warm snow
                ind = snow_only & (input1['T_pp'] >= FREEZE)
                self.precip.T_snow[ind] = FREEZE
                self.precip.h2o_sat_snow[ind] = 1.0
                
            # rain only
            if np.any(rain_only):
                self.precip.T_rain[rain_only] = input1['T_pp'][rain_only]
                
            # check the precip, temp. cannot be below freezing if rain present
            # this check was only present in Snobal and not iSnobal
#             self.precip.T_rain[self.precip.T_rain < FREEZE] = FREEZE
                
                    
        self.gridded_precip = deepcopy(self.precip)
#         self.precip_info[DATA_TSTEP] = self.precip
                
        # Clear the 'computed' flag at the other timestep levels.
#         for level in range(DATA_TSTEP, SMALL_TSTEP+1):
#             self.computed[level] = False
        
    
        #Divide the data timestep into normal run timesteps.
        self.curr_level = DATA_TSTEP   # keeps track of what time step level the model is on
        self.divide_tstep() 
     
     
    def below_thold_grid(self, threshold):
        """
        This routine determines if any individual layer's mass is below
        a given threshold for the current timestep.

        Args:
            threshold: current timestep's threshold for a 
                   layer's mass

        Returns:
            True    A layer's mass is less than the threshold.
            False    All layers' masses are greater than the threshold.
        """
        
        # if there is no snow anywhere
        if np.sum(self.snow.layer_count == 0) == self.ngrid:
            return False

        return (self.snow.m_s_0[self.snow.m_s_0 > 0] < threshold) or \
            (self.snow.m_s_l[self.snow.z_s_l > 0] < threshold)
           
    
#     @profile
    def copy_gridded(self, index, steps=None, copy_gridded_flag=True):
        """
        Copy from self.gridded* to self.snow and self.em
        """
        
        if copy_gridded_flag:
            # snow
            for key in self.snow.keys:
                setattr(self.snow, key, getattr(self.gridded_snow, key)[index])
            self.snow.shape = self.snow.m_s.shape
                
            # em
            for key in self.em.keys:
                setattr(self.em, key, getattr(self.gridded_em, key)[index])
            self.em.shape = self.em.cc_s.shape 
                
            # input1
            for key in self.gridded_input1.keys():
                self.input1[key] = self.gridded_input1[key][index]       
            
            # input2
            for key in self.gridded_input2.keys():
                self.input2[key] = self.gridded_input2[key][index]
            
            # precip
            for key in self.gridded_precip.keys:
                setattr(self.precip, key, getattr(self.gridded_precip, key)[index])
                
            # precip is a little funky since part of it needs to be divided up
            if steps is not None:
                for key in ['m_pp','m_snow','m_rain','z_snow']:
                    setattr(self.precip, key, getattr(self.gridded_precip, key)[index]/steps)
            self.precip_now = self.precip.m_pp > 0
            
            # measurement heights
            for key in self.gridded_mh.keys():
                if key != 'time_z':
                    setattr(self, key, self.gridded_mh[key][index])
                    
            # get the P_a
            self.P_a = self.gridded_P_a[index]
        
        
#     @profile
    def copy_working(self, index, copy_gridded_flag=True):
        """
        Copy the working subset back into the gridded
        Have to use some numpy trickery to index back into gridded_*
        """
        
        if copy_gridded_flag:
            # snow
            for key in self.snow.keys:
                v = getattr(self.gridded_snow, key) # this creates a reference to the class attribute
                v[index] = getattr(self.snow, key)  # update the reference
                
            # em
            for key in self.em.keys:
                v = getattr(self.gridded_em, key) # this creates a reference to the class attribute
                v[index] = getattr(self.em, key)  # update the reference
            
            if self.tstep_level > NORMAL_TSTEP:
                # input1
                for key in self.gridded_input1.keys():
                    self.gridded_input1[key][index] = self.input1[key]     
                
                # input2
                for key in self.gridded_input2.keys():
                    self.gridded_input2[key][index] = self.input2[key]
            
            # precip
#             for key in self.gridded_precip.keys:
#                 v = getattr(self.gridded_precip, key) # this creates a reference to the class attribute
#                 v[index] = getattr(self.precip, key)  # update the reference
                            
    
    def divide_inputs(self, index, steps):
        """
        Divide the inputs into the intervals
        """
        
        # Compute deltas for the climate input parameters over the data timestep.
        self.input_deltas = {}
        for i in self.gridded_input1.keys():
            self.input_deltas[i] = (self.gridded_input2[i][index] - self.gridded_input1[i][index]) / steps 
        
#         for key in ['m_pp','m_snow','m_rain','z_snow']:
#             setattr(self.precip, key, getattr(self.gridded_, key)[index]/steps)
        
    
    def calc_levels(self):
        """
        Calculate the levels at which each pixel will be calculated
        """    
        
        bins = [self.tstep_info[1]['threshold'],
                self.tstep_info[2]['threshold']]
        
        level_0 = np.digitize(self.gridded_snow.m_s_0, bins) + 1
        level_l = np.digitize(self.gridded_snow.m_s_l, bins) + 1
        level_l[self.gridded_snow.layer_count < 2] = 1   # if no snow in lower layer, normal timestep
        
        level = np.maximum.reduce([level_0, level_l])
        level[self.gridded_snow.m_s == 0] = 1    # if there is no snow, then normal timestep
        
        return level
            
    @profile
    def divide_tstep(self):
        """
        This routine divides the data timestep into the appropriate number
        of normal run timesteps.  The input values for each normal timestep
        are computed from the two input records by linear interpolation.
        
        self.gridded_snow and self.gridd_em hold the full data that was 
        provided to isnobal and will be used as the location to store the
        model state.  self.snow and self.em are the working datasets that
        will be subsetted from gridded based on the mass thresholds.
        
        Steps:
        1. Determine what level each pixel is at
        2. do_tstep for the normal time step
        3. if non normal time step
        3.1 Each medium time step, calculate the input_deltas for the entire timestep
            i.e. medium input/15 or small input/60
        3.2 update_inputs will add deltas to input to ensure that it's always current
        
        """
        
#         print np.max(self.current_time)/3600.0
<<<<<<< HEAD:pysnobal/libsnobal/py_snobal.py
        if np.max(self.current_time)/3600.0 > 801:
=======
        if np.max(self.current_time)/3600.0 > 711:
>>>>>>> f8eb08efebea59287c24b0243afdc3f5e423abae:pysnobal/lib/py_snobal.py
            self.curr_level
                      
        # determine the levels at which each pixel will be calculated
        # the level number will indicate how many steps to take 
        level = self.calc_levels()
        
        ############################################################
        # no snow and normal tstep
        level1 = (level == NORMAL_TSTEP) & self.mask
#         level1 = level1    # the mask here will ensure that the mask is propigated down
        
        if np.any(level1):
            self.do_tstep(self.tstep_info[NORMAL_TSTEP], level1, None)
           
    
        ############################################################
        # medium and small tstep
        # 
        ind = level >= MEDIUM_TSTEP
            
        if np.any(ind):
#             self.current_time = prev_time # reset the time to the start of the time step
            
            # do the medium time step first
            med_interval = self.tstep_info[MEDIUM_TSTEP]['intervals']
            for i in range(med_interval):
                
                # since gridded_input1 is updated after every time step, have to futher subdivide
                # based on how many steps are left    
                med_steps_left = med_interval - i
                sm_interval = self.tstep_info[SMALL_TSTEP]['intervals']
                sm_steps_left = sm_interval * med_steps_left
                
                # do the medium time step
                med_ind = (level == MEDIUM_TSTEP) & self.mask
                if np.any(med_ind):
                    self.divide_inputs(med_ind, med_steps_left)
                    self.do_tstep(self.tstep_info[MEDIUM_TSTEP], med_ind, med_interval)
                
                # do the small time step
                sm_ind = (level == SMALL_TSTEP) & self.mask
                if np.any(sm_ind):
                    self.divide_inputs(sm_ind, sm_steps_left)                    
                    for j in range(sm_interval):
                        
                        # set the flag for the gridded copying
                        copy_gridded_flag = False
                        copy_working_flag = False
                        if j == (sm_interval -1 ):
                            copy_working_flag = True
                        if j == 0:
                            copy_gridded_flag = True 
                            
                        self.do_tstep(self.tstep_info[SMALL_TSTEP], sm_ind, sm_interval*med_interval, 
                                      [copy_gridded_flag, copy_working_flag])
                 
                # update the levels
                prev_level = level.copy()
                level = self.calc_levels()
                level[level1] = NORMAL_TSTEP    # these have already been calculated at level1
                lvl_ind = (prev_level != 1) & (level == 1)
                if np.any(lvl_ind):
                    # if we are in the middle of the small/med timestep and it goes up
                    # to the normal_tstep, then keep at medium tstep to finsh out the time step
                    level[lvl_ind] = MEDIUM_TSTEP
                            
            # Output if this timestep is divided?
            # does a bitwise AND comparison
#         if self.tstep_info[self.curr_level]['output'] & DIVIDED_TSTEP:
# #             print '%.2f output divided tstep' % (self.current_time/3600.0)
#             self.output()
#         self.time_since_out = np.zeros(self.shape)
        
        
    
        
<<<<<<< HEAD:pysnobal/libsnobal/py_snobal.py
    @profile
    def do_tstep(self, tstep, index, step):
=======
#     @profile
    def do_tstep(self, tstep, index, step, copy_flag=[True,True]):
>>>>>>> f8eb08efebea59287c24b0243afdc3f5e423abae:pysnobal/lib/py_snobal.py
        """
        This routine performs the model's calculations for a single timestep.
        It requires that these climate variables have been initialized:
        
            S_n
            I_lw
            T_a
            e_a
            u
            T_g
            
        The routine also requires the precipitation data have been adjusted
        for the timestep, and have been stored in the array:
        
            precip_info
        
        if the flag 'precip_now' is TRUE.  The routine will set the flag
        'stop_no_snow' to TRUE if
        
            a)  the output function pointed to by 'out_func' is called, and
            b)  the flag 'run_no_snow' is FALSE, and
            c)  there is no snow remaining on the ground at the end of
                timestep
                
        Args:
            tstep: tstep_info for the current time step
            
        """
        
#         print self.current_time/3600.0
#         if self.current_time/3600.0 > 686:
#             self.curr_level
            
        # copy the working subset for    
        self.copy_gridded(index, step, copy_flag[0])
        self.index = index    # for debugging purposes
        
        self.time_step = tstep['time_step']
        self.tstep_level = tstep['level']
                
        # get the current time step precip
#         if self.precip_now:
            
        
        self.snow.reset(['h2o_total'], 0.0)
        
        # is there snowcover?
        self.snowcover = self.snow.layer_count > 0
        self.snowcover_domain = np.any(self.snowcover)
        
        if self.snowcover_domain:
            self.input1['e_g'] = core_c.sati_gridded(self.input1['T_g'])
#             self.em.cc_s_0 = self.cold_content(self.snow.T_s_0, self.snow.m_s_0)
#             self.em.cc_s_l= self.cold_content(self.snow.T_s_l, self.snow.m_s_l)
        
        # Calculate energy transfer terms
        self.e_bal()
        
        # Adjust mass and calculate runoff
        self.mass_bal()
        
        # Update the averages for the energy terms and the totals for mass
        # changes since the last output.
        if np.any(self.time_since_out[index] > 0):
            tindx = self.time_since_out[index]
            self.em.R_n_bar[:] = TIME_AVG(self.em.R_n_bar, tindx, self.em.R_n, self.time_step)
            self.em.H_bar[:] = TIME_AVG(self.em.H_bar, tindx, self.em.H, self.time_step)
            self.em.L_v_E_bar[:] = TIME_AVG(self.em.L_v_E_bar, tindx, self.em.L_v_E, self.time_step)
            self.em.G_bar[:] = TIME_AVG(self.em.G_bar, tindx, self.em.G, self.time_step)
            self.em.M_bar[:] = TIME_AVG(self.em.M_bar, tindx, self.em.M, self.time_step)
            self.em.delta_Q_bar[:] = TIME_AVG(self.em.delta_Q_bar, tindx, self.em.delta_Q, self.time_step)
            self.em.G_0_bar[:] = TIME_AVG(self.em.G_0_bar, tindx, self.em.G_0, self.time_step)
            self.em.delta_Q_0_bar[:] = TIME_AVG(self.em.delta_Q_0_bar, tindx, self.em.delta_Q_0, self.time_step)

            self.em.E_s_sum += self.em.E_s
            self.em.melt_sum += self.em.melt
            self.em.ro_pred_sum += self.em.ro_predict
    
            self.time_since_out[index] += self.time_step
            
        else:
            self.em.R_n_bar[:] = self.em.R_n
            self.em.H_bar[:] = self.em.H
            self.em.L_v_E_bar[:] = self.em.L_v_E
            self.em.G_bar[:] = self.em.G
            self.em.M_bar[:] = self.em.M
            self.em.delta_Q_bar[:] = self.em.delta_Q
            self.em.G_0_bar[:] = self.em.G_0
            self.em.delta_Q_0_bar[:] = self.em.delta_Q_0

            self.em.E_s_sum[:] = self.em.E_s
            self.em.melt_sum[:] = self.em.melt
            self.em.ro_pred_sum[:] = self.em.ro_predict
    
            self.time_since_out[index] = self.time_step
            
#         # ensure that the mask is applied
#         if self.mask is not None:
#             self.em.set_zeros(~self.mask)
#             self.snow.set_zeros(~self.mask)    
        
        # increment time
        self.current_time[index] += self.time_step
        
        # update the inputs
        self.update_inputs()
        
        # put the working results back into gridded
        self.copy_working(index, copy_flag[1])
        
        if tstep['output'] & WHOLE_TSTEP:
#             print '%.2f output stuff here' % (self.current_time/3600.0)
            self.output()
            self.time_since_out[index] = 0.0
        
        return True
        
        
    def update_inputs(self):
        """
        Update the model's input parameters
        """
        
        if self.tstep_level > NORMAL_TSTEP:
            self.input1['S_n'] += self.input_deltas['S_n']
            self.input1['I_lw'] += self.input_deltas['I_lw']
            self.input1['T_a'] += self.input_deltas['T_a']
            self.input1['e_a'] += self.input_deltas['e_a']
            self.input1['u'] += self.input_deltas['u']
            self.input1['T_g'] += self.input_deltas['T_g']
        
        
    @profile    
    def mass_bal(self):
        """
        Calculates the point mass budget for 2-layer energy budget snowmelt
        model.  It then solves for new snow temperatures.
        """
        
        # age snow by compacting snow due to time passing */
        self.time_compact()
        
        # process precipitation event
        self.precip_event()
        
        # calculate melt or freezing and adjust cold content
        self.snowmelt()
        
        # calculate evaporation and adjust snowpack 
        self.evap_cond()
        
        # compact snow due to H2O generated (melt and rain)
        self.h2o_compact()
        
        # calculate runoff, and adjust snowcover
        self.runoff()
        
        # adjust layer temps if there was a snowcover at start of the
        # timestep and there's still snow on the ground
        if self.snowcover_domain:
            
            ind = self.snow.layer_count == 1
            self.snow.T_s_0[ind] = self.new_tsno(self.snow.m_s_0[ind], self.snow.T_s_0[ind], self.em.cc_s_0[ind])
            self.snow.T_s[ind] = self.snow.T_s_0[ind]
            
            ind = self.snow.layer_count == 2
            if np.any(ind):
                self.snow.T_s_0[ind] = self.new_tsno(self.snow.m_s_0[ind], self.snow.T_s_0[ind], self.em.cc_s_0[ind])
                self.snow.T_s_l[ind] = self.new_tsno (self.snow.m_s_l[ind], self.snow.T_s_l[ind], self.em.cc_s_l[ind])
                self.snow.T_s[ind] = self.new_tsno (self.snow.m_s[ind], self.snow.T_s[ind], self.em.cc_s[ind])
             
                if np.any(self.isothermal):       
                    self.snow.T_s[self.isothermal & ind] = FREEZE
                    self.snow.T_s_l[self.isothermal & ind] = FREEZE
                    self.snow.T_s_0[self.isothermal & ind] = FREEZE
            
#             if self.snow.layer_count == 1:
#                 self.snow.T_s_0 = self.new_tsno(self.snow.m_s_0, self.snow.T_s_0, self.em.cc_s_0)
#                 self.snow.T_s = self.snow.T_s_0
#                 
#             elif self.snow.layer_count == 2:
#                 if self.isothermal:
#                     self.snow.T_s = FREEZE
#                     self.snow.T_s_l = FREEZE
#                     self.snow.T_s_0 = FREEZE
#                 else:
#                     self.snow.T_s_0 = self.new_tsno(self.snow.m_s_0, self.snow.T_s_0, self.em.cc_s_0)
#                     self.snow.T_s_l = self.new_tsno (self.snow.m_s_l, self.snow.T_s_l, self.em.cc_s_l)
#                     self.snow.T_s = self.new_tsno (self.snow.m_s, self.snow.T_s, self.em.cc_s)
        
    
    def new_tsno(self, spm, t0, ccon):
        """
        calculates a new temperature for a snow layer from its
        adjusted cold content, and the layer's last (previous) temperature.
        
        The layer's specific mass (the argument <I>spm</I>) can be computed by
        multiplying the layer's thickness (m) by its density (kg/m^3).

        Args:
            spm: layer's specific mass (kg/m^2) 
            t0: layer's last temperature (K) 
            ccon: layer's adjusted cold content (J/m^2)
            
        Returns:
            tsno: snow layer's new temperature (K)
        """
        
        cp = CP_ICE(t0)
        tdif = ccon / (spm * cp)
        tsno = tdif + FREEZE
        
        return tsno
     
    
    def runoff(self):
        """
        Calculates runoff for point energy budget 2-layer snowmelt model
        """
        
        # If no snow on ground at start of timestep or no layers currently,
        # then all water (e.g., rain) is runoff.
        
#         if (not self.snowcover) or (self.snow.layer_count == 0):
        # runoff where there isn't snow
        ind = ~self.snowcover | (self.snow.layer_count == 0)
        self.em.ro_predict[ind] = self.snow.h2o_total[ind]
        
        if np.sum(ind) == self.snow.shape[0]:
            return
        
#         if not self.snowcover_domain:
#             self.em.ro_predict += self.snow.h2o_total
#             return
        
        
        # Determine the snow density without any water, and the maximum
        # liquid water the snow can hold.
        m_s_dry = self.snow.m_s - self.snow.h2o_total
        rho_dry = m_s_dry / self.snow.z_s
        self.snow.h2o_max = H2O_LEFT(self.snow.z_s, rho_dry, self.snow.max_h2o_vol)
        # not necessary since z_s will be zero when no snow
        self.snow.h2o_max[~self.snowcover] = 0  # reset values that don't have snow
        
        # Determine runoff
        runoff = self.snow.h2o_total > self.snow.h2o_max
        ind = runoff & self.snowcover
        if np.any(ind):
            self.em.ro_predict[ind] = self.snow.h2o_total[ind] - self.snow.h2o_max[ind]
            self.snow.h2o[ind] = self.snow.h2o_max[ind]
            self.snow.h2o_sat[ind] = 1.0
            self.snow.h2o_vol[ind] = self.snow.max_h2o_vol[ind]
            
        # water left in the snow
        ind = ~runoff & self.snowcover
        if np.any(ind):
            self.em.ro_predict[ind] = 0.0
            self.snow.h2o[ind] = self.snow.h2o_total[ind]
            self.snow.h2o_sat[ind] = self.snow.h2o[ind] / self.snow.h2o_max[ind]
            self.snow.h2o_vol[ind] = self.snow.h2o_sat[ind] * self.snow.max_h2o_vol[ind]
            
        # Update the snowcover's mass for the loss of runoff.
        self.adj_snow(np.zeros(self.snow.shape), -self.em.ro_predict)
        
        
#         if self.snow.h2o_total > self.snow.h2o_max:
#             self.em.ro_predict = self.snow.h2o_total - self.snow.h2o_max
#             self.snow.h2o = self.snow.h2o_max
#             self.snow.h2o_sat = 1.0
#             self.snow.h2o_vol = self.snow.max_h2o_vol
#     
#             # Update the snowcover's mass for the loss of runoff.
#             self.adj_snow(0.0, -self.em.ro_predict)
# 
#         else:
#             self.em.ro_predict = 0.0
#             self.snow.h2o = self.snow.h2o_total
#             self.snow.h2o_sat = self.snow.h2o / self.snow.h2o_max
#             self.snow.h2o_vol = self.snow.h2o_sat * self.snow.max_h2o_vol
    

    def h2o_compact(self):
        """
        This routine compacts or densifies the snowcover based on the
        amount of liquid H2O that was added to the snowcover from melting
        and rain.  The snowcover's density is increased using the
        following "half-saturation" function:
        
            delta_rho(h2o_added) = A / (1 + B/h2o_added)
        
        A = "saturation-level" or asymtope which is the difference between
            the maximum density due to compaction by liquid H2O
            (approximately 550 kg/m^2) and the current density
        B = the point for half of the saturation level is reached (5 %)
            (h2o_added = ratio of mass of liquid h2o added by melting and
                     rain to the mass of the snowcover)
        
                ^
                |
                  A + = = = = = = = = = = = = = = = = = =
        (550 - current  |            *   *
               density) |           *
                |           *
        delta_rho |        *
        (kg/m^2)  |      *
              A/2 + . . . *
                  |     * .
                  |   *   .
                  |  *     .
                  | *     .
                  |*    .
                    0 +-------+-----------------------------+      h2o_added
                  0    B: 5 %                 1.0
  
        """
        # Maximum density due to compaction by liquid H2O added (kg/m^2)
        MAX_DENSITY = 550.0
        
        # ratio where half the difference between maximum density and
        # current density is reached (ratio from 0.0 to 1.0).
        B = 0.4
        
        if self.snowcover_domain:
        
#             if (not self.snowcover) or (self.snow.rho > MAX_DENSITY):
#                 return
            
            A = MAX_DENSITY - self.snow.rho
#             if self.precip_now:
            # should always have a value for the precip even if it's all zeros
            h2o_added = (self.em.melt + self.precip.m_rain) / self.snow.m_s
                
#             else:
#                 h2o_added = self.em.melt / self.snow.m_s
                
            # avoid dividing by zero and only calculate if below max density
            ind = (h2o_added > 0.000001) & (A > 0)
            if np.any(ind):
                            
                self.snow.rho[ind] += A[ind] / (1 + B/h2o_added[ind])
                
                # adjust the snowcover for this new density
                self.snow.z_s[ind] = self.snow.m_s[ind] / self.snow.rho[ind]
                self.adj_layers()
            
        
        
    def evap_cond(self):
        """
        Calculates mass lost or gained by evaporation/condensation
        at a point for 2-layer energy balance snowmelt model snobal.c;
        Also adjusts the liq h2o, mass and depth of the snow layer;
        Assumes that liq h2o is favored in evap as the ratio of
        vaporization to sublimation (0.882); Half the ice lost as evap
        is assumed to be lost depth; the rest reduces the density;
        
        Variables calculated:
            E_s_0: mass of evaporation to air (kg/m^2) 
            E_s_l: mass of evaporation to soil (kg/m^2) 
            E_l: mass flux by evap/cond to soil (kg/m^2/s) 
            e_g: soil vapor press 
            e_s_l: lower snow layer's vapor press 
            k: soil diffusion coef 
            prev_h2o_tot: previous value of h2o_total variable 
            q_delta: difference between snow & soil spec hum's 
            q_g: soil spec hum 
            q_s_l: lower snow layer's spec hum 
            rho_air: air density 
            T_bar: snow-soil mean temp 
        """
        
        # calculate evaporation or condensation   

        # If no snow on ground at start of timestep, then just exit.
        if not self.snowcover_domain:
            self.em.reset(['E_s'], 0.0)
            return
        
        # Total mass change due to evap/cond at surface during timestep
        E_s_0 = self.em.E * self.time_step

        # Adjust total h2o for evaporative losses
        prev_h2o_tot = self.snow.h2o_total.copy()

        ind = self.snow.h2o_total > 0.0 
        if np.any(ind):
            self.snow.h2o_total[ind] += (E_s_0[ind] * VAP_SUB)
            self.snow.h2o_total[self.snow.h2o_total <= 0] = 0
            

        # Determine total mass change due to evap/cond at soil
#         E_s_l = np.zeros_like(E_s_0)
        
#         if self.snow.layer_count == 0: 
#             E_s_l = 0.0
#         else:

#         T_bar = np.zeros_like(E_s_0)
#         e_s_l = np.zeros_like(E_s_0)
        
        # determine what the snow temperature is
        T_s = self.snow.T_s_0.copy()
        T_bar = (self.input1['T_g'] + self.snow.T_s_0) / 2.0
                
        # layer_count == 2
        ind = self.snow.layer_count == 2
        if np.any(ind):
#             e_s_l[ind] = libsnobal.sati_np(self.snow.T_s_l[ind])
            T_s[ind] = self.snow.T_s_l[ind]
            T_bar[ind] = (self.input1['T_g'][ind] + self.snow.T_s_l[ind]) / 2.0
        
        # now calculate sati
        # sati has already been calculated, the only time the snow temperature is changed
        # is during adj_snow() which will only remove temperatures when it's zero
#         e_s_l = core_c.sati_gridded(T_s)
        e_s_l = self.snow.es_0[:]
        e_s_l[ind] = self.snow.es_l[ind]
        

        q_s_l = libsnobal.spec_hum_np(e_s_l, self.P_a)
#         e_g = libsnobal.sati_np(self.input1['T_g'])
        q_g = libsnobal.spec_hum_np(self.input1['e_g'], self.P_a)
#         q_delta = q_g - q_s_l
        rho_air = libsnobal.GAS_DEN(self.P_a, libsnobal.MOL_AIR, T_bar)
        k = libsnobal.DIFFUS(self.P_a, T_bar)

        E_l = libsnobal.EVAP(rho_air, k, q_g - q_s_l, self.z_g)

        # total mass of evap/cond for time step 
        E_s_l = E_l * self.time_step
        
        # set non snowcovered pixels to zero
        E_s_l[self.snow.z_s == 0] = 0

        # adjust h2o_total for evaporative losses
        ind = self.snow.h2o_total > 0.0 
        if np.any(ind):
            self.snow.h2o_total[ind] += (E_s_l[ind] * VAP_SUB)
            self.snow.h2o_total[self.snow.h2o_total <= 0] = 0
        
    
        self.em.E_s = E_s_0 + E_s_l
    
        # adj mass and depth for evap/cond
        ind = self.snow.layer_count > 0        
        if np.any(ind):
            delta_z = np.zeros(self.em.shape)
            delta_z[ind] = ((self.em.E_s[ind] + (prev_h2o_tot[ind] - self.snow.h2o_total[ind])) / self.snow.rho[ind]) / 2.0
            self.adj_snow( delta_z, self.em.E_s)


    def snowmelt(self):
        """
        Calculates melting or re-freezing for point 2-layer energy balance
        snowmelt model.
        
        Variables calculated:
            Q_0: energy available for surface melt 
            Q_l: energy available for lower layer melt 
            Q_freeze: energy used for re-freezing 
            Q_left: energy left after re_freezing 
            h2o_refrozen: amount of liquid H2O that was refrozen 
        """
        
        if not self.snowcover_domain:
            self.em.reset(['melt'], 0.0)
            return
        
        # calculate melt or freezing, and adjust cold content
        
        # calculate surface melt
        # energy for surface melt
        Q_0 = (self.em.delta_Q_0 * self.time_step) + self.em.cc_s_0
                
        warm = Q_0 > 0
        if np.any(warm):
            self.em.melt[warm] = MELT(Q_0[warm])
            self.em.cc_s_0[warm] = 0
        
#         ind = Q_0 < 0
        if np.any(~warm):
            self.em.melt[~warm] = 0
            self.em.cc_s_0[~warm] = Q_0[~warm]
        
        ind = Q_0 == 0
        if np.any(ind):
            self.em.melt[ind] = 0
            self.em.cc_s_0[ind] = 0
        

        # layer count = 1
        Q_l = np.zeros(self.snow.shape)

        # calculate lower layer melt
        ind = self.snow.layer_count == 2
        if np.any(ind):
            Q_l[ind] = ((self.em.G[ind] - self.em.G_0[ind]) * self.time_step) + self.em.cc_s_l[ind]
            
            warm = Q_l > 0
            if np.any(warm):
                self.em.melt[warm] += MELT(Q_l[warm])
                self.em.cc_s_l[warm] = 0
            
            ind = Q_l == 0
            if np.any(ind):
                self.em.melt[ind] += 0
                self.em.cc_s_l[ind] = 0
            
#             ind = Q_l < 0
            if np.any(~warm):
                self.em.melt[~warm] += 0
                self.em.cc_s_l[~warm] = Q_l[~warm]
                   
        
        self.snow.h2o_total += self.em.melt
     
        # adjust layers for re-freezing
        # adjust surface layer
        h2o_refrozen = np.zeros(self.snow.shape)
        Q_freeze = np.zeros(self.snow.shape)
        Q_left = np.zeros(self.snow.shape) 
        
        cold = self.em.cc_s_0 < 0       # if the snowpack is cold
        water = self.snow.h2o_total > 0 # if there is liquid water
        cold_water = cold & water       # if there is liquied water and it's cold
        
        # adjust surface layer for refreezing
        if np.any(cold_water):
            Q_freeze[cold_water] = self.snow.h2o_total[cold_water] * (self.snow.z_s_0[cold_water]/self.snow.z_s[cold_water]) * libsnobal.LH_FUS(FREEZE)
            Q_left[cold_water] = Q_0[cold_water] + Q_freeze[cold_water]
            
            # if there is CC left where there is snow
            ind = (Q_left < 0) & cold_water & (self.snow.z_s > 0)
            h2o_refrozen[ind] = self.snow.h2o_total[ind] * (self.snow.z_s_0[ind]/self.snow.z_s[ind])
            self.em.cc_s_0[ind] = Q_left[ind]
            
            # else the snowpack is melting
            ind = (Q_left >= 0) & cold_water & (self.snow.z_s > 0)
            h2o_refrozen[ind] = self.snow.h2o_total[ind] * (self.snow.z_s_0[ind]/self.snow.z_s[ind]) - MELT(Q_left[ind])
            self.em.cc_s_0[ind] = 0
        
            
        # adjust lower layer for re-freezing
        cold_lower = (self.snow.layer_count == 2) & (self.em.cc_s_l < 0.0) & water
        if np.any(cold_lower):
            Q_freeze[cold_lower] = self.snow.h2o_total[cold_lower] * (self.snow.z_s_l[cold_lower]/self.snow.z_s[cold_lower]) * libsnobal.LH_FUS(FREEZE)
            Q_left[cold_lower] = Q_l[cold_lower] + Q_freeze[cold_lower]
            
            ind = (Q_left < 0) & cold_lower & (self.snow.z_s > 0)
            h2o_refrozen[ind] += self.snow.h2o_total[ind] * (self.snow.z_s_l[ind]/self.snow.z_s[ind])
            self.em.cc_s_l[ind] = Q_left[ind]
            
            ind = (Q_left >= 0) & cold_lower & (self.snow.z_s > 0)
            h2o_refrozen[ind] += self.snow.h2o_total[ind] * (self.snow.z_s_l[ind]/self.snow.z_s[ind]) - MELT(Q_left[ind])
            self.em.cc_s_l[ind] = 0
            
     
        # Note:  because of rounding errors, h2o_refrozen may not
        # be exactly the same as h2o_total.  Check for this
        # case, and if so, then just zero out h2o_total.
#         ind = np.abs(self.snow.h2o_total - h2o_refrozen) <= 1e-8
        self.snow.h2o_total -= h2o_refrozen
        self.snow.h2o_total[self.snow.h2o_total <= 1e-8] = 0
            
        # determine if snowcover is isothermal
        self.isothermal = np.zeros(self.snow.shape, dtype=bool)
#         iso_1lay = (self.snow.layer_count == 1) & (self.em.cc_s_0 == 0.0)
#         iso_2lay = (self.snow.layer_count == 2) & (self.em.cc_s_0 == 0.0) & (self.em.cc_s_l == 0.0)
#         self.isothermal[iso_1lay | iso_2lay] = True
        
        
        iso = (self.em.cc_s_0 == 0.0) & (self.em.cc_s_l == 0.0)
        self.isothermal[iso] = True
        
            
        # adjust depth and density for melt
        if np.sum(self.em.melt) > 0:
            self.adj_snow((-1)*self.em.melt/self.snow.rho, np.zeros(self.snow.shape))
            
        # set total cold content
#         ind = self.snow.layer_count == 2
        self.em.cc_s = self.em.cc_s_0 + self.em.cc_s_l
        
#         ind = self.snow.layer_count == 1
#         self.em.cc_s[ind] = self.em.cc_s_0[ind]
                
                    
    
    def precip_event(self):
        """
        This routine processes a precipitation event, i.e., the current
        precip record, if there's one for the current timestep.  It
        determines if the precip is rain or snow which increases the
        snowcover.
        """
        
        if np.any(self.precip.now):
            
            # if there is precip, but no snow, then create a snowpack
            if np.any(~self.snowcover):
                
                # set the values from the initial snow properties
                # [:] performs a quick deep copy since snow is already initialized
                no_snow = ~self.snowcover
                tsnow = self.precip.T_snow[no_snow]
                self.snow.z_s[no_snow] = self.precip.z_snow[no_snow]
                self.snow.rho[no_snow] = self.input1['rho_snow'][no_snow]
                self.snow.T_s[no_snow] = tsnow
                self.snow.T_s_0[no_snow] = tsnow
                self.snow.T_s_l[no_snow] = tsnow
                self.snow.h2o_sat[no_snow] = self.precip.h2o_sat_snow[no_snow]
            
                # reset the precip values to zero in case adj_snow needs to be called below
                self.precip.set_value(['z_snow','m_pp','h2o_sat_snow'], no_snow, 0.0)
                
                self.init_snow()
            
            
            
            if self.snowcover_domain:
                # Adjust snowcover's depth and mass by snowfall's
                # depth and the total precipitation mass.
                self.adj_snow(self.precip.z_snow, self.precip.m_pp)
            
                # Determine the additional liquid water that's in
                # the snowfall, and then add its mass to liquid
                # water in the whole snowcover.
                h2o_vol_snow = self.precip.h2o_sat_snow * self.snow.max_h2o_vol
                self.snow.h2o += H2O_LEFT(self.precip.z_snow, 
                                          self.input1['rho_snow'],
                                          h2o_vol_snow)
                
#             else:
#                 
#                 # set the values from the initial snow properties
#                 # [:] performs a quick deep copy since snow is already initialized
#                 self.snow.z_s[:] = self.precip.z_snow
#                 self.snow.rho[:] = self.input1['rho_snow']
#                 self.snow.T_s[:] = self.precip.T_snow
#                 self.snow.T_s_0[:] = self.precip.T_snow
#                 self.snow.T_s_l[:] = self.precip.T_snow
#                 self.snow.h2o_sat[:] = self.precip.h2o_sat_snow
#                 
#                 self.init_snow()
            
            # Add rainfall and water in the snowcover to the total liquid water
            self.snow.h2o_total += self.snow.h2o + self.precip.m_rain
                
        else:
            # Add water in the snowcover to total liquid water
            self.snow.h2o_total += self.snow.h2o
            
#     @profile        
    def adj_snow(self, delta_z_s, delta_m_s):
        """
        This routine adjusts the snowcover for a change in its depth or
        its mass or both.  The snowcover's density is updated.  If there
        is a change in the snowcover's depth, the # of layers is recomputed.
        If there's just a change in the snowcover's mass with no change in
        its depth, then just the specific masses for the layers are updated.
        
        The routine ensures that the snowcover's density does NOT exceed
        a maximum density (currently 750 kg/m^3).  If the adjustments to
        the snowcover, for some reason, lead to an excessive density, the
        density is clipped at the maximum, and the depth re-adjusted
        accordingly.
        
        Args:
            delta_z_s: change in snowcover depth
            delta_m_s: change in snowcover mass
        """
        
        # first check for nan which means rho=0 usually
        # and only adjust the snow where there is snow
        ind = ~np.isnan(delta_z_s) & self.snowcover
        
        # Update depth, mass, and then recompute density
        self.snow.z_s[ind] += delta_z_s[ind]
        self.snow.m_s[ind] += delta_m_s[ind]
        
        ind = self.snow.z_s > 0.0
        self.snow.rho[ind] = self.snow.m_s[ind] / self.snow.z_s[ind]
        
        self.snow.rho[~ind] = 0
            
        # clip density at maximum density if necessary
        mrho = self.snow.rho > MAX_SNOW_DENSITY
        if np.any(mrho):
            self.snow.rho[mrho] = MAX_SNOW_DENSITY
            self.snow.z_s[ind] = self.snow.m_s[ind] / self.snow.rho[ind]
#             
#         if self.snow.rho > MAX_SNOW_DENSITY:
# #             self.snow.rho = MAX_SNOW_DENSITY
#             self.snow.z_s = self.snow.m_s / self.snow.rho
#             self.adj_layers()
            
        
        # If a change in depth, adjust the layers' depths and masses
        # Only need to adjust the layers if there is some change in depth,
        # this will ensure that runoff() will not remove a layer after the
        # runoff has already been calculated
        if np.any(delta_z_s != 0.0):
            self.adj_layers()
#         else:
        # Just change in the snowcover's mass, so update the layer masses
        self.layer_mass()
    

    def time_compact(self):
        """
        This routine "ages" the snowcover by accounting for the compaction
        or densification by gravity as time passes.  The snowcover's
        density is increased using the following "half-saturation" function:
        
            rho(time) = A / (1 + B/time)
        
        A = "saturation-level" or asymtope which is the maximum density
            due to compaction by gravity (approximately 350 kg/m^2)
        B = the point for half of the saturation level is reached (10 days)
        
                ^
                |
             A: 350 + = = = = = = = = = = = = = = = = = =
                |            *   *
                |           *
        rho     |           *
        (kg/m^2)    |        *
                |      *
               A/2: 175 + . . . *
                  |     * .
                  |   *   .
                  |  *    .
                  | *     .
                  |*    .
                    0 +-------+----------------------------->
                  0    B: 10 days        time
        """
        
        # Maximum density due to compaction by gravity (kg/m^2)
        A = 350.0
        
        # Time when half "saturation", i.e., maximum density is reached (seconds)
        # (864000 = 10 days * 24 hours/day * 60 mins/hr * 60 secs/min)
        B = 864000.0
        
        # If the snow is already at or above the maximum density due
        # compaction by gravity, then just leave.
#         if (not self.snowcover) or (self.snow.rho > A):
#             return
        ind = (self.snowcover) & (self.snow.rho < A)
        
        if np.any(ind):
        
            # Given the current density, determine where on the time axis
            # we are (i.e., solve the function above for "time").
            time = B / ((A / self.snow.rho[ind]) - 1)
            
            # Move along the time axis by the time step, and calculate the
            # density at this new time.
            self.snow.rho[ind] = A / (1 + B/(time + self.time_step))
            
            # Adjust the snowcover for this new density
            self.snow.z_s[ind] = self.snow.m_s[ind] / self.snow.rho[ind]
            
            self.adj_layers()
        
        
    def adj_layers(self):
        """
        This routine adjusts the layers of the snowcover because the
        snowcover's depth has changed.  It is assumed that the snowcover's
        density has already been updated.  The # of layers are recomputed
        based on the overall snowcover depth.  Their depths and masses
        are updated as well.  If a layer has been created due to an
        increase in the snowcover's depth, its temperature and cold content
        are initialized. 
        
        Recompute then number of layers and see if there's been
        a change in the # of layers.  Note:  since this routine
        is called to adjust an existing snowcover, the current # of
        layers must be either 1 or 2 while the new # of layers may
        either be 0, 1 or 2.
        
          current #    new #
          of layers    of layers
        
             1       -->       0
             1       -->       1    (no change)
             1       -->       2
             2       -->       0
             2       -->       1
             2       -->       2    (no change)

        """
        
        prev_layer_count = self.snow.layer_count.copy()
        
        self.calc_layers()
        
        # set locations with layer_count back to zeros
        ind = self.snow.layer_count == 0
        if np.any(ind):
            
            # if mass > 0, then it must be below threshold.
            # So turn this little bit of mass into water
            # In Snobal, this would be added to ro_predict, but since this is a little different,
            # this water won't get added correctly, so added it to ro_predict and set h2o_total to zero
#             if self.snow.m_s > 0:
            self.snow.h2o_total[ind] += self.snow.m_s[ind]
#             self.em.ro_predict[ind] = self.snow.h2o_total[ind] + self.snow.m_s[ind]
            
            # reset some values back to zero
            index = ['h2o', 'h2o_sat', 'h2o_max', 'h2o_vol', 'm_s', 'm_s_0', 'rho']
            self.snow.set_value(index, ind, 0)
            
            # Note: Snow temperatures are set to MIN_SNOW_TEMP
            # (as degrees K) instead of 0 K to keep quantization
            # range in output image smaller.
            self.snow.set_value(['T_s','T_s_0','T_s_l'], ind, MIN_SNOW_TEMP + FREEZE)
            
            # reset the cold contents
            self.em.set_value(['cc_s', 'cc_s_0', 'cc_s_l'], ind, 0)
            
            # upadate the snowcover
            self.snowcover[ind] = False
        
        # upadate the other areas with layers
#         if np.any(~ind):
        self.layer_mass()
           
        # 1 layer --> 2 layers, add lower layer
        ind = (prev_layer_count == 1) & (self.snow.layer_count == 2) 
        if np.any(ind):
            self.snow.T_s_l[ind] = self.snow.T_s[ind]
            self.em.cc_s_l[ind] = self.cold_content(self.snow.T_s_l[ind], self.snow.m_s_l[ind])
            
        # 2 layers --> 1 layer, remove lower layer
        ind = (prev_layer_count == 2) & (self.snow.layer_count == 1)
        if np.any(ind):
            self.snow.T_s_l[ind] = MIN_SNOW_TEMP + FREEZE
            self.em.cc_s_l[ind] = 0    
            
            
        
#         if self.snow.layer_count == 0:
#             # 1 or 2 layers ---> 0
#             self.snow.rho = 0
#             
#             # if mass > 0, then it must be velow threshold.
#             # So turn this little bit of mass into water
#             if self.snow.m_s > 0:
#                 self.snow.h2o_total += self.snow.m_s
#             
#             # set a bunch of values to 0
#             index = ['h2o', 'h2o_max', 'h2o_total', 'h2o_vol', 'm_s', 'cc_s', 'm_s_0', 'cc_s_0']
#             for i in index:
#                 setattr(self.snow, i, 0.0)
#                         
#             # Note: Snow temperatures are set to MIN_SNOW_TEMP
#             # (as degrees K) instead of 0 K to keep quantization
#             # range in output image smaller.
#             
#             self.snow.T_s = MIN_SNOW_TEMP + FREEZE
#             self.snow.T_s_0 = MIN_SNOW_TEMP + FREEZE
#             
#             if prev_layer_count == 2:
#                 self.snow.m_s_l = 0
#                 self.em.cc_s_l = 0
#                 self.snow.T_s_l = MIN_SNOW_TEMP + FREEZE
#                 
#             self.snowcover = False
#             
#         else:
#             self.layer_mass()
#             
#             if (prev_layer_count == 1) and (self.snow.layer_count == 2):
#                 # 1 layer --> 2 layers, add lower layer
#                 self.snow.T_s_l = self.snow.T_s
#                 self.em.cc_s_l = self.cold_content(self.snow.T_s_l, self.snow.m_s_l)
#                 
#             elif (prev_layer_count == 2) and (self.snow.layer_count == 1):
#                 # 2 layers --> 1 layer, remove lower layer
#                 self.snow.T_s_l = MIN_SNOW_TEMP + FREEZE
#                 self.em.cc_s_l = 0
            
        
    @profile
    def e_bal(self):
        """
        Calculates point energy budget for 2-layer snowcover.
        """
        
        # if there is a snowcover
#         if self.snow.layer_count > 0:
        if self.snowcover_domain:
             
            # precalculate sati
            self.snow.es_0 = core_c.sati_gridded(self.snow.T_s_0)
            self.snow.es_l = core_c.sati_gridded(self.snow.T_s_l)
             
             
            # Calculates net allwave radiation from the net solar radiation
            #incoming thermal/longwave radiation, and the snow surface
            # temperature.
            # replaces _net_rad()
            self.em.R_n = self.input1['S_n'] + (SNOW_EMISSIVITY * (self.input1['I_lw'] - STEF_BOLTZ * np.power(self.snow.T_s_0, 4)))
             
            # calculate H & L_v_E (and E as well)
            self.h_le()
              
            # calculate G & G_0 (conduction/diffusion heat xfr)
#             if self.snow.layer_count == 1:
            g = self.g_soil('surface')
            self.em.G_0[:] = g
            self.em.G[:] = g
                  
            ind = self.snow.layer_count == 2
            if np.any(ind):
                g = self.g_soil ('lower')
                self.em.G[ind] = g[ind]
                gs = self.g_snow()
                self.em.G_0[ind] = gs[ind]
                  
            # calculate advection
            self.advec()
              
            # sum E.B. terms
            # surface energy budget
            self.em.delta_Q_0 = self.em.R_n + self.em.H + self.em.L_v_E + self.em.G_0 + self.em.M
              
            # total snowpck energy budget
#             if self.snow.layer_count == 1:
            self.em.delta_Q[:] = self.em.delta_Q_0
#             else:
            ind = self.snow.layer_count == 2
            if np.any(ind):
                self.em.delta_Q[ind] += self.em.G[ind] - self.em.G_0[ind]
                 
            # since this is iSNOWbal, remove the energy balance components where
            # there is no snow
            self.em.set_value(['R_n', 'H', 'L_v_E', 'E', 'G', 'G_0', 'delta_Q', 'delta_Q_0'], 
                              ~self.snowcover, 0)
                 
        else:
            self.em.reset(['R_n', 'H', 'L_v_E', 'E', 'G', 'G_0', 'delta_Q', 'delta_Q_0'], 0)
        
               
                       
    def advec(self):
        """
        This routine calculates the advected energy for a 2-layer snowcover
        if there's precipitation for the current timestep.
        """    
        
        M = np.zeros(self.snow.shape)
        
        if np.any(self.precip.now):
                        
            M = self.heat_stor(CP_WATER(self.precip.T_rain), \
                               self.precip.m_rain, \
                               self.precip.T_rain - self.snow.T_s_0) + \
                 self.heat_stor(CP_ICE(self.precip.T_snow), \
                                self.precip.m_snow, \
                                self.precip.T_snow - self.snow.T_s_0)
                                
            M /= self.time_step
            
            # set locations without precip back to zero
            M[~self.precip.now | ~self.snowcover] = 0
            
#         else:
#             M = 0
            
        self.em.M = M
        
        
#     @profile
    def g_soil(self, layer):
        """
        conduction heat flow between snow and soil
        """
        
        if layer == 'surface':
            tsno = self.snow.T_s_0
            ds = self.snow.z_s_0
            es = self.snow.es_0
            
        else:
            tsno = self.snow.T_s_l
            ds = self.snow.z_s_l
            es = self.snow.es_l
            
        if np.any(tsno > FREEZE):
            warnings.warn('g_soil: tsno > %8.2f; set to %8.2f\n' % (FREEZE, FREEZE))
            tsno[tsno > FREEZE] = FREEZE
        
        # set effective soil conductivity
        # /***    changed to KT_MOISTSAND by D. Marks, NWRC, 09/30/2003    ***/
        # /***    based on heat flux data from RMSP            ***/
        # /***    note: Kt should be passed as an argument        ***/
        # /***    k_g = efcon(KT_WETSAND, tg, pa);            ***/
        kts = KT_MOISTSAND * np.ones(self.P_a.shape)
        k_g = core_c.efcon_gridded(kts, self.input1['T_g'], self.P_a, self.input1['e_g'])
        
        # calculate G    
        # set snow conductivity
        kcs = KTS(self.snow.rho)
        k_s = core_c.efcon_gridded(kcs, tsno, self.P_a, es)
        g = libsnobal.ssxfr(k_s, k_g, tsno, self.input1['T_g'], ds, self.z_g)
        
        
        return g
        
    
    def g_snow(self):
        """
        conduction heat flow between snow layers
        """
        
        # calculate g, if the difference in temperature is zero, g will be zero
#         g = np.zeros(self.shape)
#         ind = self.snow.T_s_0 != self.snow.T_s_l
        
#         if self.snow.T_s_0 == self.snow.T_s_l:
#             g = 0
#         else:

        # simplify the equations, since rho is for the whole snowpack, we don't need
        # to calculate KTS for both layers and just use the same for each layer
        kcs = KTS(self.snow.rho)
#         kcs2 = KTS(self.snow.rho)
        k_s1 = core_c.efcon_gridded(kcs, self.snow.T_s_0, self.P_a, self.snow.es_0)
        k_s2 = core_c.efcon_gridded(kcs, self.snow.T_s_l, self.P_a, self.snow.es_l)
        
        g = libsnobal.ssxfr(k_s1, k_s2, self.snow.T_s_0, self.snow.T_s_l, self.snow.z_s_0, self.snow.z_s_l)
        
        return g
            
#     @profile
    def h_le(self):
        """
        Calculates point turbulent transfer (H and L_v_E) for a 2-layer snowcover
        """
        
#         if self.snow.T_s_0 < 0:
#             raise Exception('T_s_0 is below 0 K')
#         if self.input1['T_a'] < 0:
#             raise Exception('T_a is below 0 K')
        
        # calculate saturation vapor pressure
#         e_s = libsnobal.sati_np(self.snow.T_s_0)
        
        # error check for bad vapor pressures
        sat_vp = core_c.sati_gridded(self.input1['T_a'])
        
        ind = self.input1['e_a'] > sat_vp
        if np.any(ind):
            self.input1['e_a'][ind] = sat_vp[ind]
            
        # determine relative measurement heights
        if self.relative_hts:
            rel_z_t = self.z_t
            rel_z_u = self.z_u
        else:
            rel_z_t = self.z_t - self.snow.z_s
            rel_z_u = self.z_u - self.snow.z_s
        
        # calculate H & L_v_E, the three different version all do the exact same thing but at different complexity and speed
        
#         2nd fastest, using numpy
#         H, L_v_E, E, status = libsnobal.hle1_np2(self.P_a, self.input1['T_a'], self.snow.T_s_0, rel_z_t, \
#                                              self.input1['e_a'], self.snow.es_0, rel_z_t, self.input1['u'], rel_z_u, \
#                                              self.z_0)
#         
#         # 3rd, looping in pyx file
#         H2, L_v_E2, E2, status = core_c.hle1_c(self.P_a, self.input1['T_a'], self.snow.T_s_0, rel_z_t, \
#                                              self.input1['e_a'], self.snow.es_0, rel_z_t, self.input1['u'], rel_z_u, \
#                                              self.z_0)
        
        # fastest, passing all variables to C to perform the loop
        H, L_v_E, E, status = core_c.hle1_gridded(self.P_a, self.input1['T_a'], self.snow.T_s_0, rel_z_t, \
                                             self.input1['e_a'], self.snow.es_0, rel_z_t, self.input1['u'], rel_z_u, \
                                             self.z_0, error_check=0)
                
        if status != 0:
            idx = np.where(self.index)
            raise Exception("hle1 did not converge at point y=%i, x=%i) , sorry... :(" % (idx[0][status], idx[1][status]))
            
        self.em.H = H
        self.em.L_v_E = L_v_E
        self.em.E = E
         
        
    def below_thold(self, threshold):
        """
        This routine determines if any individual layer's mass is below
        a given threshold for the current timestep.

        Args:
            threshold: current timestep's threshold for a 
                   layer's mass

        Returns:
            True    A layer's mass is less than the threshold.
            False    All layers' masses are greater than the threshold.
        """
        
        # if there is no snow anywhere
        if np.sum(self.snow.layer_count == 0) == self.ngrid:
            return False
#         if self.snow.layer_count == 1:
#             return self.snow.m_s < threshold
#         else:
        return np.any(self.snow.m_s_0[self.snow.m_s_0 > 0] < threshold) or \
            np.any(self.snow.m_s_l[self.snow.z_s_l > 0] < threshold)
            
           
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
            
            self.time_s = self.snow_records['time_s'] * HR_TO_SEC
            self.curr_time_hrs = self.snow_records['time_s'] * HR_TO_SEC
            self.start_time = self.snow_records['time_s'] * HR_TO_SEC
            self.more_sn_recs = True
            
            self.current_time = self.time_s * np.ones(self.shape)
            
        
        
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
            
            self.gridded_mh = {}
            self.gridded_mh['time_z'] = self.mh_prop['time_z'] * HR_TO_SEC
            self.gridded_mh['z_u'] = self.mh_prop['z_u'] * np.ones_like(self.snow.z_s)
            self.gridded_mh['z_t'] = self.mh_prop['z_t'] * np.ones_like(self.snow.z_s)
            self.gridded_mh['z_0'] = self.snow_records['z_0']
            self.gridded_mh['z_g'] = self.mh_prop['z_g'] * np.ones_like(self.snow.z_s)
            
        else:
            self.mh_prop_index += 1
            
        self.more_mh_recs = False
        
    
    def init_precip(self):
        """
        Returns a Pandas series that will contain all information about
        current time step precip
            
        """
        
        cols = ['m_pp','m_snow','m_rain','z_snow','T_rain','T_snow','h2o_sat_snow','now']
        
        return EmptyClass(cols, self.shape)
        
#         return PrecipClass()  
#         return pd.Series(index=['m_pp','m_snow','m_rain','z_snow'])
#         return {key: 0.0 for key in ['m_pp','m_snow','m_rain','z_snow']}
                
    
    def init_snow(self, from_record=False):
        """
        This routine initializes the properties for the snowcover.  It
        determines the number of layers, their individual properties,
        the cold content for the snowcover and its layers, and the
        snowcover's water content.
        
        Initialize all the following values
        h2o_sat, layer_count, max_h2o_vol, rho, T_s, T_s_0, T_s_l,
        z_s, cc_s, cc_s_0, cc_s_l, h2o, h2o_max, h2o_total, h2o_vol, m_s, m_s_0, m_s_l
        
        Args:
            from_record: whether or not to read from the snow.properties file,
                else it will reinitialize based on 
        
            z_s: depth of snowcover (m)
            rho: density of snowcover (kg/m^3)
            T_s: average temperature of snowcover (K)
            T_s_0: temperature of surface layer of snowcover (K)
            T_s_l: temperature of lower layer of snowcover (K)
            h2o_sat:  % of liquid h2o saturation (0 to 1.0)
            max_h2o_vol: maximum liquid h2o content as volume ratio:
                        V_water/(V_snow - V_ice) (unitless)
                        
        """
        
        cols = ['h2o_sat', 'layer_count', 'max_h2o_vol', 'rho', 
                'T_s', 'T_s_0', 'T_s_l', 
                'z_s', 'z_s_0', 'z_s_l', 
                'm_s', 'm_s_0', 'm_s_l',
                'es_0', 'es_l',
                'h2o', 'h2o_max', 'h2o_total', 'h2o_vol']
        
#                 'cc_s', 'cc_s_0', 'cc_s_l',
               
        # set the values from the initial snow properties
        if from_record:
            
            # create an empty dataframe
#             self.snow = pd.Series(index=cols)
            
            # create a dict instead
#             self.snow = Map({key: 0.0 for key in cols})
            self.snow = EmptyClass(cols, self.shape)
        
            self.snow.z_s[:] = self.snow_records['z_s']
            self.snow.rho[:] = self.snow_records['rho']
            self.snow.T_s_0[:] = self.snow_records['T_s_0']
            self.snow.T_s[:] = self.snow_records['T_s']
            self.snow.h2o_sat[:] = self.snow_records['h2o_sat']
            self.snow.max_h2o_vol[:] = self.params['max_h2o_vol'] * np.ones_like(self.snow.z_s)
        
        else:
            try:
                getattr(self, 'snow')
            except:
                raise Exception('The snow has not been intitialized')
        
        # initialize the snowpack
        self.snow.m_s = self.snow.rho * self.snow.z_s
        
        # determine the number of layers and layer heights
        self.calc_layers()
        
        # If mass > 0, then it must be below threshold.
        # So turn this little bit of mass into water
        # and reset all the point's snow values to zero or minimum value
        ind = self.snow.layer_count == 0
        if np.any(ind):
            self.snow.h2o_total[ind] += self.snow.m_s[ind]
#             self.em.ro_predict[ind] = self.snow.h2o_total[ind] + self.snow.m_s[ind]
            self.snow.set_value(['rho','m_s','m_s_0','m_s_l','h2o_vol','h2o','h2o_max','h2o_sat'], 
                                ind, 0)
            
            # Note: Snow temperatures are set to MIN_SNOW_TEMP
            # (as degrees K) instead of 0 K to keep quantization
            # range in output image smaller.
            self.snow.set_value(['T_s', 'T_s_0', 'T_s_l'], ind, MIN_SNOW_TEMP + FREEZE)
         
        # Compute the specific mass and cold content for each layer
        # this will just reset a layer to zero if layer_count=0
#         ind = ~ind
#         if np.any(ind):
        self.layer_mass()
        self.em.cc_s_0 = self.cold_content(self.snow.T_s_0, self.snow.m_s_0)
        
        lc2 = self.snow.layer_count == 2
#         if np.any(lc2):
        self.em.cc_s_l[lc2] = self.cold_content(self.snow.T_s_l[lc2], self.snow.m_s_l[lc2])
#         else:
        self.snow.set_value(['T_s_l'], ~lc2, MIN_SNOW_TEMP + FREEZE)
        self.em.set_value(['cc_s_l'], ~lc2, 0)
            
        
#         if self.snow.layer_count == 2:
#             self.em.cc_s_l = self.cold_content(self.snow.T_s_l, self.snow.m_s_l)
#         else:
#             self.snow.T_s_l = MIN_SNOW_TEMP + FREEZE
#             self.em.cc_s_l = 0
            
            
        # Compute liquid water content as volume ratio, and
        # snow density without water
        self.snow.h2o_vol = self.snow.h2o_sat * self.snow.max_h2o_vol
        rho_dry = DRY_SNO_RHO(self.snow.rho, self.snow.h2o_vol)
        
        # Determine the maximum liquid water content (as specific mass)
        # and the actual liquid water content (as specific mass)
        self.snow.h2o_max = H2O_LEFT(self.snow.z_s, rho_dry, self.snow.max_h2o_vol)
        self.snow.h2o = self.snow.h2o_sat * self.snow.h2o_max
        
        
#         if self.snow.layer_count == 0:
#             # If mass > 0, then it must be below threshold.
#             # So turn this little bit of mass into water
#             
#             if self.snow.m_s > 0.0:
#                 self.snow.h2o_total += self.snow.m_s
#             
#             for col in ['rho','m_s','cc_s','m_s_0','cc_s_0','m_s_l','cc_s_l','h2o_vol','h2o','h2o_max','h2o_sat']:
#                 setattr(self.snow, col, 0.0)        
#         
#             # Note: Snow temperatures are set to MIN_SNOW_TEMP
#             # (as degrees K) instead of 0 K to keep quantization
#             # range in output image smaller.
#             for col in ['T_s', 'T_s_0', 'T_s_l']:
#                 setattr(self.snow, col, MIN_SNOW_TEMP + FREEZE)
# #                 self.snow[col] = MIN_SNOW_TEMP + FREEZE
#          
#         else:
#             # Compute the specific mass and cold content for each layer   
#             self.layer_mass()
#             self.em.cc_s_0 = self.cold_content(self.snow.T_s_0, self.snow.m_s_0)
#             
#             if self.snow.layer_count == 2:
#                 self.em.cc_s_l = self.cold_content(self.snow.T_s_l, self.snow.m_s_l)
#             else:
#                 self.snow.T_s_l = MIN_SNOW_TEMP + FREEZE
#                 self.em.cc_s_l = 0
#             
#             # Compute liquid water content as volume ratio, and
#             # snow density without water
#             self.snow.h2o_vol = self.snow.h2o_sat * self.snow.max_h2o_vol
#             rho_dry = DRY_SNO_RHO(self.snow.rho, self.snow.h2o_vol)
#             
#             # Determine the maximum liquid water content (as specific mass)
#             # and the actual liquid water content (as specific mass)
#             self.snow.h2o_max = H2O_LEFT(self.snow.z_s, rho_dry, self.snow.max_h2o_vol)
#             self.snow.h2o = self.snow.h2o_sat * self.snow.h2o_max


    def init_water_content_gridded(self):
        """
        iSnobal calls init_snow() at the beginning of every time step.  Since pysnobal
        has a more persistant state variables, we only need to update a couple of the
        snowpack state variables at the start of every time step
        """
        
        # Compute liquid water content as volume ratio, and
        # snow density without water
        self.gridded_snow.h2o_vol = self.gridded_snow.h2o_sat * self.gridded_snow.max_h2o_vol
        rho_dry = DRY_SNO_RHO(self.gridded_snow.rho, self.gridded_snow.h2o_vol)
        
        # Determine the maximum liquid water content (as specific mass)
        # and the actual liquid water content (as specific mass)
        self.gridded_snow.h2o_max = H2O_LEFT(self.gridded_snow.z_s, rho_dry, self.gridded_snow.max_h2o_vol)
        self.gridded_snow.h2o = self.gridded_snow.h2o_sat * self.gridded_snow.h2o_max

    
    def init_em(self):
        """
        Initialize a Pandas Series for the energy and max fluxes
        """
        
        col = [
                # energy balance info for current timestep
                'R_n',            # net allwave radiation (W/m^2) 
                'H',              # sensible heat xfr (W/m^2) 
                'L_v_E',          # latent heat xfr (W/m^2) 
                'G',              # heat xfr by conduction & diffusion from soil to snowcover (W/m^2) 
                'G_0',            # heat xfr by conduction & diffusion from soil or lower layer to active layer (W/m^2) 
                'M',              # advected heat from precip (W/m^2) 
                'delta_Q',        # change in snowcover's energy (W/m^2) 
                'delta_Q_0',      # change in active layer's energy (W/m^2) 
                
                #   averages of energy balance vars since last output record   
                'R_n_bar',
                'H_bar',
                'L_v_E_bar',
                'G_bar',
                'G_0_bar',
                'M_bar',
                'delta_Q_bar',
                'delta_Q_0_bar',
                
                #   mass balance vars for current timestep      
                'melt',           # specific melt (kg/m^2 or m) 
                'E',        # mass flux by evap into air from active layer (kg/m^2/s) 
                'E_s',        # mass of evap into air & soil from snowcover (kg/m^2) 
                'ro_predict',     # predicted specific runoff (m/sec) 
                
                #   sums of mass balance vars since last output record   
                'melt_sum',
                'E_s_sum',
                'ro_pred_sum',
                
                # cold content values
                'cc_s', 
                'cc_s_0', 
                'cc_s_l'
               ]
        
#         self.em = pd.Series(data=np.zeros(len(col)), index=col)
#         self.em = Map({key: 0.0 for key in col})
        self.em = EmptyClass(col, self.shape)
        



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
        

        
#         # assume that there are two layers
#         z_s_0 = self.params['max_z_s_0'] * np.ones(self.shape)
#         z_s_l = self.snow.z_s - z_s_0
#         
#         # correct if there is only one or no layers
#         ind = z_s_l < 0
#         z_s_0[ind] = z_s_l[ind] + self.params['max_z_s_0']
#         
#         # make sure that there is enough mass total
#         ind = self.snow.m_s < self.tstep_info[SMALL_TSTEP]['threshold']
#         z_s_0[ind] = 0
#         z_s_l[ind] = 0

#         layer_count = self.snow.layer_count
#         z_s_0 = self.snow.z_s_0     # these are references to self.snow.z_s_0
#         z_s_l = self.snow.z_s_l
        
        # concatenate into 3D array
#         z = np.dstack((z_s_0, z_s_l))
        
        # check the surface layer
        
        
        
        # preallocate
#         layer_count = np.zeros(self.shape)
#         z_s = np.zeros_like(layer_count)
        layer_count = self.snow.layer_count
        z_s_0 = self.snow.z_s_0     # these are references to self.snow.z_s_0
        z_s_l = self.snow.z_s_l
#         z_s = self.snow.z_s
        
        # less than minimum layer mass, so treat as no snowcover
        # but since everything was preallocated, then it shouldn't matter
        no_mass = self.snow.m_s <= self.tstep_info[SMALL_TSTEP]['threshold']
        layer_count[no_mass] = 0
        z_s_0[no_mass] = 0
        z_s_l[no_mass] = 0

        
#         if np.sum(no_mass) == self.ngrid: 
#          
#             layer_count[:] = 0
#             z_s_0[:] = 0
#             z_s_l[:] = 0
#             self.snow.z_s[:] = 0
#             return
        
        # not enough depth for surface layer and the lower layer,
        # so just 1 layer: surface layer
        surf_only = self.snow.z_s < self.params['max_z_s_0']
        
        ind = surf_only & ~no_mass
        layer_count[ind] = 1
        z_s_0[ind] = self.snow.z_s[ind]
        z_s_l[ind] = 0
        
        # enough depth for both layers
#         ind = self.snow.z_s >= self.params['max_z_s_0']
#         if np.any(ind):
        layer_count[~surf_only] = 2
        z_s_0[~surf_only] = self.params['max_z_s_0']
        z_s_l[~surf_only] = self.snow.z_s[~surf_only] - z_s_0[~surf_only]
#         z_s[ind] = z_s_0[ind] + z_s_l[ind] # not really needed but needed for below
        
        # However, make sure there's enough MASS for the lower
        # layer.  If not, then there's only 1 layer
        ind = (z_s_l * self.snow.rho < self.tstep_info[SMALL_TSTEP]['threshold']) & ~surf_only
#         if np.any(ind):
        layer_count[ind] = 1
        z_s_0[ind] = self.snow.z_s[ind]
        z_s_l[ind] = 0
    
        
#         if self.snow.m_s <= self.tstep_info[SMALL_TSTEP]['threshold']:
#             # less than minimum layer mass, so treat as no snowcover
#             
#             layer_count = 0
#             z_s = z_s_0 = z_s_l = 0
#             
#         elif self.snow.z_s < self.params['max_z_s_0']:
#             # not enough depth for surface layer and the lower layer,
#             # so just 1 layer: surface layer
#             
#             layer_count = 1
#             z_s_0 = self.snow.z_s
#             z_s_l = 0
#             z_s = z_s_0
#             
#         else:
#             # enough depth for both layers
#             
#             layer_count = 2
#             z_s_0 = self.params['max_z_s_0']
#             z_s_l = self.snow.z_s - z_s_0
#             z_s = z_s_0 + z_s_l # not really needed but needed for below
#             
#             # However, make sure there's enough MASS for the lower
#             # layer.  If not, then there's only 1 layer
#             if z_s_l * self.snow.rho < self.tstep_info[SMALL_TSTEP]['threshold']:
#                 layer_count = 1
#                 z_s_0 = self.snow.z_s
#                 z_s_l = 0
        
            
#         self.snow.layer_count = layer_count
        self.snow.z_s = z_s_0 + z_s_l
#         self.snow.z_s_0 = z_s_0
#         self.snow.z_s_l = z_s_l
        
        
#     @profile
    def layer_mass(self):
        """
        This routine computes the specific mass for each snow layer in
        the snowcover.  A layer's mass is based its depth and the
        average snowcover density.
        """
        
        
        # see if I can get away with this as rho and z_s should already
        # be updated
#         ind = self.snow.layer_count == 0
#         self.snow.set_value(['m_s_0','m_s_l'], ind, 0.0)
        
        self.snow.m_s_0 = self.snow.rho * self.snow.z_s_0
        self.snow.m_s_l = self.snow.rho * self.snow.z_s_l
        
        
#         if self.snow.layer_count == 0:
#             self.snow.m_s_0 = 0
#             self.snow.m_s_l = 0
#             
#         else:
#             # layer count is 1 or 2
#             self.snow.m_s_0 = self.snow.rho * self.snow.z_s_0
#             
#             if self.snow.layer_count == 2:
#                 self.snow.m_s_l = self.snow.rho * self.snow.z_s_l
#             else:
#                 self.snow.m_s_l = 0
                
    
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
        
        cc = np.zeros_like(temp)
        ind = temp < FREEZE
        if np.any(ind):
            cc[ind] = self.heat_stor(CP_ICE(temp[ind]), mass[ind], temp[ind]-FREEZE)
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


    def output(self):
        """
        Specify where the model output should go
        """
        
        # write out to a file
        if self.params['out_filename'] is not None:
            
            curr_time_hrs = SEC_TO_HR(self.current_time)
            
#             
#             # time
#             self.params['out_file'].write('%g,' % curr_time_hrs)
#             
#             # energy budget terms
#             self.params['out_file'].write("%.1f,%.1f,%.1f,%.1f,%.1f,%.1f," % \
#                     (self.em.R_n_bar, self.em.H_bar, self.em.L_v_E_bar, \
#                     self.em.G_bar, self.em.M_bar, self.em.delta_Q_bar))
# 
#             # layer terms
#             self.params['out_file'].write("%.1f,%.1f," % \
#                     (self.em.G_0_bar, self.em.delta_Q_0_bar))
#     
#             # heat storage and mass changes
#             self.params['out_file'].write("%.6e,%.6e,%.6e," % \
#                     (self.em.cc_s_0, self.em.cc_s_l, self.em.cc_s))
#             self.params['out_file'].write("%.5f,%.5f,%.5f," % \
#                     (self.em.E_s_sum, self.em.melt_sum, self.em.ro_pred_sum))
#     
# #             # runoff error if data included */
# #             if (ro_data)
# #                 fprintf(out, " %.1f",
# #                         (ro_pred_sum - (ro * time_since_out)))
#     
#             # sno properties */
#             self.params['out_file'].write("%.3f,%.3f,%.3f,%.1f," % \
#                     (self.snow.z_s_0, self.snow.z_s_l, self.snow.z_s, self.snow.rho))
#             self.params['out_file'].write("%.1f,%.1f,%.1f,%.1f," % \
#                     (self.snow.m_s_0, self.snow.m_s_l, self.snow.m_s, self.snow.h2o))
#             if self.params['temps_in_C']:
#                 self.params['out_file'].write("%.2f,%.2f,%.2f\n" % 
#                         (K_TO_C(self.snow.T_s_0), K_TO_C(self.snow.T_s_l), K_TO_C(self.snow.T_s)))
#             else:
#                 self.params['out_file'].write("%.2f,%.2f,%.2f\n" % \
#                         (self.snow.T_s_0, self.snow.T_s_l, self.snow.T_s))
            
            
            # time
            self.params['out_file'].write('%g,' % curr_time_hrs)
             
            # energy budget terms
            self.params['out_file'].write("%.3f,%.3f,%.3f,%.3f,%.3f,%.3f," % \
                    (self.em.R_n_bar, self.em.H_bar, self.em.L_v_E_bar, \
                    self.em.G_bar, self.em.M_bar, self.em.delta_Q_bar))
 
            # layer terms
            self.params['out_file'].write("%.3f,%.3f," % \
                    (self.em.G_0_bar, self.em.delta_Q_0_bar))
     
            # heat storage and mass changes
            self.params['out_file'].write("%.9e,%.9e,%.9e," % \
                    (self.em.cc_s_0, self.em.cc_s_l, self.em.cc_s))
            self.params['out_file'].write("%.8f,%.8f,%.8f," % \
                    (self.em.E_s_sum, self.em.melt_sum, self.em.ro_pred_sum))
     
#             # runoff error if data included */
#             if (ro_data)
#                 fprintf(out, " %.3f",
#                         (ro_pred_sum - (ro * time_since_out)))
     
            # sno properties */
            self.params['out_file'].write("%.6f,%.6f,%.6f,%.3f," % \
                    (self.snow.z_s_0, self.snow.z_s_l, self.snow.z_s, self.snow.rho))
            self.params['out_file'].write("%.3f,%.3f,%.3f,%.3f," % \
                    (self.snow.m_s_0, self.snow.m_s_l, self.snow.m_s, self.snow.h2o))
            if self.params['temps_in_C']:
                self.params['out_file'].write("%.5f,%.5f,%.5f\n" % 
                        (K_TO_C(self.snow.T_s_0), K_TO_C(self.snow.T_s_l), K_TO_C(self.snow.T_s)))
            else:
                self.params['out_file'].write("%.5f,%.5f,%.5f\n" % \
                        (self.snow.T_s_0, self.snow.T_s_l, self.snow.T_s))
    
    
class PrecipClass():
    """
    Simple precip class to hold values that will function as the Map()
    """
    
    def __init__(self):
        self.m_pp = 0.0
        self.m_snow = 0.0
        self.m_rain = 0.0
        self.z_snow = 0.0
        
        
class InputClass():
    """
    Simple class to hold input data
    """            
    
    def _init__(self):
        self.S_n = None
        self.I_lw = None
        self.T_a = None
        self.e_a = None
        self.u = None
        self.T_g = None
       
            
class EmptyClass:
    """
    The class holds the model state variables, either for
    em or snow.
    
    The em_out and snow_out are the 2D arrays
    The em and snow are 1D raveled arrays of the above
    """
    def __init__(self, cols, size, mask=None):
        self.keys = cols
        self.shape = size
        self.mask = mask
        
#         if mask is not None:
#             mask1d = np.where(mask.flatten())
        
#         a = np.zeros(size)
#         a.fill(np.nan)
        for c in cols:
#             if mask is not None:
#                 setattr(self, c, ma.masked_array(np.zeros(size), 
#                                                  mask=mask, hard_mask=True))
#             else:
#             setattr(self, '%s_out' % c, np.zeros(size))
            
#             if mask is not None:
#                 r = np.ravel(getattr(self, '%s_out' % c))
#                 setattr(self, c, r[mask1d])
#             else:

            setattr(self, c, np.zeros(size))
            
            
    def set_zeros(self, index):
        """
        Set the attributes to zeros for the given index
        
        Args:
            index: numpy array
        """

        for c in self.keys:
            self.__dict__[c][index] = 0.0
            
    def set_nan(self, index):
        """
        Set the attributes to nan for the given index
        
        Args:
            index: numpy array
        """

        for c in self.keys:
            self.__dict__[c][index] = np.nan
            
    def reset(self, keys, value):
        """
        Reset specified keys back to a given value
        """
        for c in keys:
#             if self.mask is not None:
#                 arr = ma.masked_array(value * np.ones(self.shape), mask=self.mask)
#                 setattr(self, c, arr)
#             else:
            setattr(self, c, value * np.ones(self.shape))
            
            
    def set_value(self, keys, index, value):
        """
        Set the attributes to zeros for the given index
        
        Args:
            keys: attributes to change
            index: numpy array
            value: value to change to
        """

        if not np.any(index):
            pass

        for c in keys:
            self.__dict__[c][index] = value
            
    def check_nan(self):
        """
        Checks to see if any of the fields have nan values
        """
        
        val = []
        
        for c in self.keys:
            val.append(np.any(np.isnan(self.__dict__[c])))
            
        return val
        
                
# class Map(dict):
#     """
#     Example:
#     m = Map({'first_name': 'Eduardo'}, last_name='Pool', age=24, sports=['Soccer'])
#     """
#     def __init__(self, *args, **kwargs):
#         super(Map, self).__init__(*args, **kwargs)
#         for arg in args:
#             if isinstance(arg, dict):
#                 for k, v in arg.iteritems():
#                     self[k] = v
# 
#         if kwargs:
#             for k, v in kwargs.iteritems():
#                 self[k] = v
# 
#     def __getattr__(self, attr):
#         return self.get(attr)
# 
#     def __setattr__(self, key, value):
#         self.__setitem__(key, value)
# 
#     def __setitem__(self, key, value):
#         super(Map, self).__setitem__(key, value)
#         self.__dict__.update({key: value})
# 
#     def __delattr__(self, item):
#         self.__delitem__(item)
# 
#     def __delitem__(self, key):
#         super(Map, self).__delitem__(key)
#         del self.__dict__[key]            
        