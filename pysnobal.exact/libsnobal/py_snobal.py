"""
Class snobal() that will hold all the modeling components

20160109 Scott Havens
"""

import libsnobal
import numpy as np
# import pandas as pd
import warnings
from copy import copy

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

# thermal conductivity of snow (J/(m sec K))
# (after Yen, 1965, see Anderson, 1976, pg. 31)
#     rho = snow density (kg/m^3)
KTS = lambda rho: CAL_TO_J * 0.0077 * (rho/1000.0) * (rho/1000.0)

# melt (kg/m^2), Q = available energy (J/m^2)
MELT = lambda Q: Q / libsnobal.LH_FUS(FREEZE)

SNOW_EMISSIVITY = 0.98
STEF_BOLTZ = 5.67032e-8     # Stefan-Boltzmann constant (W / m^2 / deg^4)

# A macro to update a time-weighted average for a quantity.
#    avg        current average
#    total_time    the time interval the current average applies to
#    value        new value to be averaged in
#    time_incr    the time interval the new value applies to
TIME_AVG = lambda avg,total_time,value,time_incr: (avg * total_time + value * time_incr) / (total_time + time_incr) 

class snobal(object):
    """
    
    To-do
    self.snow is the working data frame only containing one record
    self.output will be a dataframe that self.snow gets added 
    
    """
    
    # came from self.__dict__.keys()
    # slots will help with memory when multiple instances of snobal are used
    __slots__ = ['em', 'input2', 'input1', 'z_T', 'computed', 'precip_now', 
                 'snow_records', 'time_step', 'precip_info', 'tstep_level', 
                 'current_time', 'z_u', 'time_since_out', 'time_s', 'snow_prop_index', 
                 'input_deltas', 'ro_data', 'snow', 'time_z', 'P_a', 'params', 'z_g', 
                 'next_level', 'elevation', 'mh_prop_index', 'start_time', 'tstep_info', 
                 'curr_time_hrs', 'mh_prop', 'z_0', 'more_sn_recs', 'relative_hts', 
                 'curr_level', 'precip', 'more_mh_recs', 'snowcover', 'isothermal']
    
    
    def __init__(self, params, tstep_info, snow_prop, meas_heights):
        """
        Initialize the snobal() class with the parameters,
        time step information,
        
        This follows the initialize() function in snobal
        
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
        
        # get the intial snowcover properties
        self.snow_records = snow_prop
        self.get_sn_rec(True)
        
        # initialize the snowcover
        self.init_snow(True)
        
        # initialize the em
        self.init_em()
        
        # initialize the precip for the time step
#         self.init_precip()
        
        # get measurement-height record
        self.mh_prop = meas_heights
        self.get_mh_rec(True)
        self.relative_hts = False
        
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
            self.input_deltas[level] = Map({})
            self.computed[level] = None
            self.precip_info[level] = Map({})
        
        
        self.time_since_out = 0
        
        # set some attributes
        self.isothermal = False
        
        
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
        
        input1 = Map(input1)
        input2 = Map(input2)
        
        # store the inputs for later
        self.input1 = input1
        self.input2 = input2
        
        # Compute deltas for the climate input parameters over the data timestep.
        for i in input1.keys():
            self.input_deltas[DATA_TSTEP][i] = input2[i] - input1[i] 
#         self.input_deltas[DATA_TSTEP] = input2.subtract(input1)
        
        # If there is precipitation, then compute the amount of rain & snow in it.
        # Look at the first input record
        pp_info = self.init_precip()
        self.precip_now = False
        
        # this precip will hold the temperatures and such, set to zero initially
#         self.precip = pd.Series(0.0, index=['T_rain','T_snow','h2o_sat_snow'])
        self.precip = Map({key: 0.0 for key in ['T_rain','T_snow','h2o_sat_snow']})
        
        if input1.m_pp > 0:
            self.precip_now = True
            
            pp_info.m_pp = input1.m_pp
            pp_info.m_snow = input1.m_pp * input1.percent_snow
            pp_info.m_rain = input1.m_pp - pp_info.m_snow
            
            if (pp_info.m_snow > 0.0):
                if (input1.rho_snow > 0.0):
                    pp_info.z_snow = pp_info.m_snow / input1.rho_snow
                else:
                    raise ValueError('input1.rho_snow is <= 0.0 with input1.percent_snow > 0.0')
            else:
                pp_info.z_snow = 0
                
            # Mixed snow and rain
            if (pp_info.m_snow > 0) and (pp_info.m_rain > 0):
                self.precip.T_snow = FREEZE
                self.precip.h2o_sat_snow = 1.0
                self.precip.T_rain = input1.T_pp
            
            elif (pp_info.m_snow > 0):
                # Snow only
                if (input1.T_pp < FREEZE):
                    # cold snow
                    self.precip.T_snow = input1.T_pp
                    self.precip.h2o_sat_snow = 0
                else:
                    # warm snow
                    self.precip.T_snow = FREEZE
                    self.precip.h2o_sat_snow = 1
                    
            elif (pp_info.m_rain > 0):
                # rain only
                self.precip.T_rain = input1.T_pp
                    
        self.precip_info[DATA_TSTEP] = pp_info
                
        # Clear the 'computed' flag at the other timestep levels.
        for level in range(DATA_TSTEP, SMALL_TSTEP+1):
            self.computed[level] = False
        
        if self.current_time/3600.0 > 1036.99:
            self.curr_level
    
        #Divide the data timestep into normal run timesteps.
        self.curr_level = DATA_TSTEP   # keeps track of what time step level the model is on
        self.divide_tstep() 
           
    
    def divide_tstep(self):
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
                
        """
           
#         print "Current level --> %i" % self.curr_level   
        
        # Fetch the record for the timestep at the next level.
        self.next_level = self.curr_level + 1
        next_lvl_tstep = self.tstep_info[self.next_level]
        
        # get the input deltas
        curr_lvl_deltas = self.input_deltas[self.curr_level]
        next_lvl_deltas = self.input_deltas[self.next_level]
        
        # get the precip info
        curr_lvl_precip = self.precip_info[self.curr_level]
        next_lvl_precip = self.precip_info[self.next_level]
        
        
        # If this is the first time this new level has been used during
        # the current data timestep, then calculate its input deltas
        # and precipitation values.
        
        # To-do can get rid of computed and use the None
        
        if not self.computed[self.next_level]:
#             next_lvl_deltas = curr_lvl_deltas / next_lvl_tstep['intervals']
#             self.input_deltas[self.next_level] = next_lvl_deltas.copy()
            for k in curr_lvl_deltas.keys():
                self.input_deltas[self.next_level][k] = curr_lvl_deltas[k] / next_lvl_tstep['intervals']
            
            if self.precip_now:
#                 next_lvl_precip = curr_lvl_precip / next_lvl_tstep['intervals']
#                 self.precip_info[self.next_level] = next_lvl_precip.copy()
                for k in curr_lvl_precip.keys():
                    self.precip_info[self.next_level][k] = curr_lvl_precip[k] / next_lvl_tstep['intervals']
                
            
            self.computed[self.next_level] = True
            
        # For each the new smaller timestep, either subdivide them if
        # below their mass threshold, or run the model for them.
        interval = next_lvl_tstep['intervals']
        for i in range(interval):
#             print "Current level --> %i, loop %i" % (next_lvl_tstep['level'], i)
            
            if (self.next_level != SMALL_TSTEP) and (self.below_thold(next_lvl_tstep['threshold'])):
                self.curr_level = copy(self.next_level) # increment the level number
                if not self.divide_tstep():
                    return False
            else:
                if not self.do_tstep(next_lvl_tstep):
                    return False
        
        
        # Output if this timestep is divided?
        # does a bitwise AND comparison
        if self.tstep_info[self.curr_level]['output'] & DIVIDED_TSTEP:
#             print '%.2f output divided tstep' % (self.current_time/3600.0)
            self.output()
            self.time_since_out = 0.0
            
        self.curr_level -= 1
        self.next_level -= 1
        
        return True
        
        
#     @profile
    def do_tstep(self, tstep):
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
        
        if self.current_time/3600.0 > 1068.36:
            self.curr_level
        
        self.time_step = tstep['time_step']
        self.tstep_level = tstep['level']
        
        # get the current time step precip
#         if self.precip_now:
            
        
        self.snow.h2o_total = 0
        
        # is there snowcover?
        self.snowcover = self.snow.layer_count > 0
        
        # Calculate energy transfer terms
        self.e_bal()
        
        # Adjust mass and calculate runoff
        self.mass_bal()
        
        # Update the averages for the energy terms and the totals for mass
        # changes since the last output.
        if self.time_since_out > 0:
            self.em.R_n_bar = TIME_AVG(self.em.R_n_bar, self.time_since_out, self.em.R_n, self.time_step)
            self.em.H_bar = TIME_AVG(self.em.H_bar, self.time_since_out, self.em.H, self.time_step)
            self.em.L_v_E_bar = TIME_AVG(self.em.L_v_E_bar, self.time_since_out, self.em.L_v_E, self.time_step)
            self.em.G_bar = TIME_AVG(self.em.G_bar, self.time_since_out, self.em.G, self.time_step)
            self.em.M_bar = TIME_AVG(self.em.M_bar, self.time_since_out, self.em.M, self.time_step)
            self.em.delta_Q_bar = TIME_AVG(self.em.delta_Q_bar, self.time_since_out, self.em.delta_Q, self.time_step)
            self.em.G_0_bar = TIME_AVG(self.em.G_0_bar, self.time_since_out, self.em.G_0, self.time_step)
            self.em.delta_Q_0_bar = TIME_AVG(self.em.delta_Q_0_bar, self.time_since_out, self.em.delta_Q_0, self.time_step)

            self.em.E_s_sum += self.em.E_s
            self.em.melt_sum += self.em.melt
            self.em.ro_pred_sum += self.em.ro_predict
    
            self.time_since_out += self.time_step
            
        else:
            self.em.R_n_bar = self.em.R_n
            self.em.H_bar = self.em.H
            self.em.L_v_E_bar = self.em.L_v_E
            self.em.G_bar = self.em.G
            self.em.M_bar = self.em.M
            self.em.delta_Q_bar = self.em.delta_Q
            self.em.G_0_bar = self.em.G_0
            self.em.delta_Q_0_bar = self.em.delta_Q_0

            self.em.E_s_sum = self.em.E_s
            self.em.melt_sum = self.em.melt
            self.em.ro_pred_sum = self.em.ro_predict
    
            self.time_since_out = self.time_step
            
        # increment time
        self.current_time += self.time_step
        
        if tstep['output'] & WHOLE_TSTEP:
#             print '%.2f output stuff here' % (self.current_time/3600.0)
            self.output()
            self.time_since_out = 0.0
            
        # Update the model's input parameters
        self.input1.S_n += self.input_deltas[tstep['level']].S_n
        self.input1.I_lw += self.input_deltas[tstep['level']].I_lw
        self.input1.T_a += self.input_deltas[tstep['level']].T_a
        self.input1.e_a += self.input_deltas[tstep['level']].e_a
        self.input1.u += self.input_deltas[tstep['level']].u
        self.input1.T_g += self.input_deltas[tstep['level']].T_g
        
        return True
        
#     @profile    
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
        if self.snowcover:
            if self.snow.layer_count == 1:
                self.snow.T_s_0 = self.new_tsno(self.snow.m_s_0, self.snow.T_s_0, self.snow.cc_s_0)
                self.snow.T_s = self.snow.T_s_0
                
            elif self.snow.layer_count == 2:
                if self.isothermal:
                    self.snow.T_s = FREEZE
                    self.snow.T_s_l = FREEZE
                    self.snow.T_s_0 = FREEZE
                else:
                    self.snow.T_s_0 = self.new_tsno(self.snow.m_s_0, self.snow.T_s_0, self.snow.cc_s_0)
                    self.snow.T_s_l = self.new_tsno (self.snow.m_s_l, self.snow.T_s_l, self.snow.cc_s_l)
                    self.snow.T_s = self.new_tsno (self.snow.m_s, self.snow.T_s, self.snow.cc_s)
        
    
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
        
        if (not self.snowcover) or (self.snow.layer_count == 0):
            self.em.ro_predict = self.snow.h2o_total
            return
        
        # Determine the snow density without any water, and the maximum
        # liquid water the snow can hold.
        m_s_dry = self.snow.m_s - self.snow.h2o_total
        rho_dry = m_s_dry / self.snow.z_s
        self.snow.h2o_max = H2O_LEFT(self.snow.z_s, rho_dry, self.snow.max_h2o_vol)
        
        # Determine runoff, and water left in the snow
        if self.snow.h2o_total > self.snow.h2o_max:
            self.em.ro_predict = self.snow.h2o_total - self.snow.h2o_max
            self.snow.h2o = self.snow.h2o_max
            self.snow.h2o_sat = 1.0
            self.snow.h2o_vol = self.snow.max_h2o_vol
    
            # Update the snowcover's mass for the loss of runoff.
            self.adj_snow(0.0, -self.em.ro_predict)

        else:
            self.em.ro_predict = 0.0
            self.snow.h2o = self.snow.h2o_total
            self.snow.h2o_sat = self.snow.h2o / self.snow.h2o_max
            self.snow.h2o_vol = self.snow.h2o_sat * self.snow.max_h2o_vol
    
     
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
        MAX_DENSITY = 550
        
        # ratio where half the difference between maximum density and
        # current density is reached (ratio from 0.0 to 1.0).
        B = 0.4
        
        if (not self.snowcover) or (self.snow.rho > MAX_DENSITY):
            return
        
        A = MAX_DENSITY - self.snow.rho
        if self.precip_now:
            h2o_added = (self.em.melt + self.precip_info[self.tstep_level].m_rain) / self.snow.m_s
            
        else:
            h2o_added = self.em.melt / self.snow.m_s
            
        if h2o_added > 0.000001:
            self.snow.rho += A / (1 + B/h2o_added)
            
            # adjust the snowcover for this new density
            self.snow.z_s = self.snow.m_s / self.snow.rho
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
        if not self.snowcover:
            self.em.E_s = 0.0
            return
        
        # Total mass change due to evap/cond at surface during timestep
        E_s_0 = self.em.E * self.time_step

        # Adjust total h2o for evaporative losses
        prev_h2o_tot = self.snow.h2o_total

        if self.snow.h2o_total > 0.0:
            self.snow.h2o_total += (E_s_0 * VAP_SUB)
            if self.snow.h2o_total <= 0.0:
                self.snow.h2o_total = 0.0

        # Determine total mass change due to evap/cond at soil
        if self.snow.layer_count == 0: 
            E_s_l = 0.0
        else:
            if self.snow.layer_count == 2:
                e_s_l = libsnobal.sati(self.snow.T_s_l)
                T_bar = (self.input1.T_g + self.snow.T_s_l) / 2.0
            
            else:  # layer_count == 1 
                e_s_l = libsnobal.sati(self.snow.T_s_0)
                T_bar = (self.input1.T_g + self.snow.T_s_0) / 2.0
            
    
            q_s_l = libsnobal.spec_hum(e_s_l, self.P_a)
            e_g = libsnobal.sati(self.input1.T_g)
            q_g = libsnobal.spec_hum(e_g, self.P_a)
            q_delta = q_g - q_s_l
            rho_air = libsnobal.GAS_DEN(self.P_a, libsnobal.MOL_AIR, T_bar)
            k = libsnobal.DIFFUS(self.P_a, T_bar)
    
            E_l = libsnobal.EVAP(rho_air, k, q_delta, self.z_g)
    
            # total mass of evap/cond for time step 
            E_s_l = E_l * self.time_step
    
            # adjust h2o_total for evaporative losses
            if self.snow.h2o_total > 0.0:
                self.snow.h2o_total += (E_s_l * VAP_SUB)
                if self.snow.h2o_total <= 0.0:
                    self.snow.h2o_total = 0.0
            
        
    
        self.em.E_s = E_s_0 + E_s_l
    
        # adj mass and depth for evap/cond        
        if self.snow.layer_count > 0:
            self.adj_snow( ((self.em.E_s + (prev_h2o_tot - self.snow.h2o_total)) / self.snow.rho) / 2.0, self.em.E_s)


    
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
        
        if not self.snowcover:
            self.em.melt = 0
            return
        
        # calculate melt or freezing, and adjust cold content
        
        # calculate surface melt
        # energy for surface melt
        Q_0 = (self.em.delta_Q_0 * self.time_step) + self.snow.cc_s_0
        
        if Q_0 > 0:
            self.em.melt = MELT(Q_0)
            self.snow.cc_s_0 = 0
        elif Q_0 == 0:
            self.em.melt = 0
            self.snow.cc_s_0 = 0
        else:
            self.em.melt = 0
            self.snow.cc_s_0 = Q_0

        # calculate lower layer melt
        if self.snow.layer_count == 2:
            Q_l = ((self.em.G - self.em.G_0) * self.time_step) + self.snow.cc_s_l
            
            if Q_l > 0:
                self.em.melt += MELT(Q_l)
                self.snow.cc_s_l = 0
            elif Q_l == 0:
                self.em.melt += 0
                self.snow.cc_s_l = 0
            else:
                self.em.melt += 0
                self.snow.cc_s_l = Q_l
            
        else:
            # layer count = 1
            Q_l = 0
            
        
        self.snow.h2o_total += self.em.melt
     
        # adjust layers for re-freezing
        # adjust surface layer
        h2o_refrozen = 0
        
        if self.snow.cc_s_0 < 0:
            # if liquid h2o present, calc refreezing and adj cc_s_0
            if self.snow.h2o_total > 0:
                Q_freeze = self.snow.h2o_total * (self.snow.z_s_0/self.snow.z_s) * libsnobal.LH_FUS(FREEZE)
                Q_left = Q_0 + Q_freeze
                
                if Q_left <= 0:
                    h2o_refrozen = self.snow.h2o_total * (self.snow.z_s_0/self.snow.z_s)
                    self.snow.cc_s_0 = Q_left
                else:
                    h2o_refrozen = self.snow.h2o_total * (self.snow.z_s_0/self.snow.z_s) - MELT(Q_left)
                    self.snow.cc_s_0 = 0
                    
        # adjust lower layer for re-freezing
        if (self.snow.layer_count == 2) and (self.snow.cc_s_l < 0.0):
            # if liquid h2o, calc re-freezing and adj cc_s_l
            if self.snow.h2o_total > 0.0:
                Q_freeze = self.snow.h2o_total * (self.snow.z_s_l/self.snow.z_s) * libsnobal.LH_FUS(FREEZE)
                Q_left = Q_l + Q_freeze
            
                if Q_left <= 0.0:
                        h2o_refrozen += self.snow.h2o_total * (self.snow.z_s_l/self.snow.z_s)
                        self.snow.cc_s_l = Q_left
                else:
                    h2o_refrozen += ((self.snow.h2o_total * (self.snow.z_s_l/self.snow.z_s)) - MELT(Q_left))
                    self.snow.cc_s_l = 0.0

     
        # Note:  because of rounding errors, h2o_refrozen may not
        # be exactly the same as h2o_total.  Check for this
        # case, and if so, then just zero out h2o_total.
        if np.abs(self.snow.h2o_total - h2o_refrozen) <= 1e-8:
            self.snow.h2o_total = 0
        else:
            self.snow.h2o_total -= h2o_refrozen
            
        # determine if snowcover is isothermal
        if (self.snow.layer_count == 2) and (self.snow.cc_s_0 == 0.0) and (self.snow.cc_s_l == 0.0):
            self.isothermal = True
        elif (self.snow.layer_count == 1) and (self.snow.cc_s_0 == 0.0):
            self.isothermal = True   
        else:
            self.isothermal = False
            
        # adjust depth and density for melt
        if self.em.melt > 0:
            self.adj_snow((-1)*self.em.melt/self.snow.rho, 0)
            
        # set total cold content
        if self.snow.layer_count == 2:
            self.snow.cc_s = self.snow.cc_s_0 + self.snow.cc_s_l
        elif self.snow.layer_count == 1:
            self.snow.cc_s = self.snow.cc_s_0
         
                    
    
    def precip_event(self):
        """
        This routine processes a precipitation event, i.e., the current
        precip record, if there's one for the current timestep.  It
        determines if the precip is rain or snow which increases the
        snowcover.
        """
        
        if self.precip_now:
            if self.snowcover:
                # Adjust snowcover's depth and mass by snowfall's
                # depth and the total precipitation mass.
                self.adj_snow(self.precip_info[self.tstep_level].z_snow, self.precip_info[self.tstep_level].m_pp)
            
                # Determine the additional liquid water that's in
                # the snowfall, and then add its mass to liquid
                # water in the whole snowcover.
                h2o_vol_snow = self.precip.h2o_sat_snow * self.snow.max_h2o_vol
                self.snow.h2o += H2O_LEFT(self.precip_info[self.tstep_level].z_snow, 
                                          self.input1.rho_snow,
                                          h2o_vol_snow)
                
            else:
                
                # set the values from the initial snow properties
                self.snow.z_s = self.precip_info[self.tstep_level].z_snow
                self.snow.rho = self.input1.rho_snow
                self.snow.T_s = self.precip.T_snow
                self.snow.T_s_0 = self.precip.T_snow
                self.snow.T_s_l = self.precip.T_snow
                self.snow.h2o_sat = self.precip.h2o_sat_snow
                
                self.init_snow()
            
            # Add rainfall and water in the snowcover to the total liquid water
            self.snow.h2o_total += self.snow.h2o + self.precip_info[self.tstep_level].m_rain
                
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
        
        # Update depth, mass, and then recompute density
        self.snow.z_s += delta_z_s
        self.snow.m_s += delta_m_s
        
        if self.snow.z_s != 0.0:
            self.snow.rho = self.snow.m_s / self.snow.z_s
        else:
            self.snow.rho = 0
            
        # clip ensity at maximum density if necessary
        if self.snow.rho > MAX_SNOW_DENSITY:
            self.snow.rho = MAX_SNOW_DENSITY
            self.snow.z_s = self.snow.m_s / self.snow.rho
            self.adj_layers()
            
        else:
            # If a change in depth, adjust the layers' depths and masses
            if delta_z_s != 0.0:
                self.adj_layers()
            else:
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
        A = 350
        
        # Time when half "saturation", i.e., maximum density is reached (seconds)
        # (864000 = 10 days * 24 hours/day * 60 mins/hr * 60 secs/min)
        B = 864000
        
        # If the snow is already at or above the maximum density due
        # compaction by gravity, then just leave.
        if (not self.snowcover) or (self.snow.rho > A):
            return
        
        # Given the current density, determine where on the time axis
        # we are (i.e., solve the function above for "time").
        time = B / ((A / self.snow.rho) - 1)
        
        # Move along the time axis by the time step, and calculate the
        # density at this new time.
        self.snow.rho = A / (1 + B/(time + self.time_step))
        
        # Adjust the snowcover for this new density
        self.snow.z_s = self.snow.m_s / self.snow.rho
        
        self.adj_layers()
        
        
#     @profile    
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
        
        prev_layer_count = self.snow.layer_count
        
        self.calc_layers()
        
        if self.snow.layer_count == 0:
            # 1 or 2 layers ---> 0
            self.snow.rho = 0
            
            # if mass > 0, then it must be velow threshold.
            # So turn this little bit of mass into water
            if self.snow.m_s > 0:
                self.snow.h2o_total += self.snow.m_s
            
            # set a bunch of values to 0
            index = ['h2o', 'h2o_max', 'h2o_total', 'h2o_vol', 'm_s', 'cc_s', 'm_s_0', 'cc_s_0']
            for i in index:
                self.snow[i] = 0
                        
            # Note: Snow temperatures are set to MIN_SNOW_TEMP
            # (as degrees K) instead of 0 K to keep quantization
            # range in output image smaller.
            
            self.snow.T_s = MIN_SNOW_TEMP + FREEZE
            self.snow.T_s_0 = MIN_SNOW_TEMP + FREEZE
            
            if prev_layer_count == 2:
                self.snow.m_s_l = 0
                self.snow.cc_s_l = 0
                self.snow.T_s_l = MIN_SNOW_TEMP + FREEZE
                
            self.snowcover = False
            
        else:
            self.layer_mass()
            
            if (prev_layer_count == 1) and (self.snow.layer_count == 2):
                # 1 layer --> 2 layers, add lower layer
                self.snow.T_s_l = self.snow.T_s
                self.snow.cc_s_l = self.cold_content(self.snow.T_s_l, self.snow.m_s_l)
                
            elif (prev_layer_count == 2) and (self.snow.layer_count == 1):
                # 2 layers --> 1 layer, remove lower layer
                self.snow.T_s_l = MIN_SNOW_TEMP + FREEZE
                self.snow.cc_s_l = 0
            
        
#     @profile
    def e_bal(self):
        """
        Calculates point energy budget for 2-layer snowcover.
        """
        
        # if there is a snowcover
        if self.snow.layer_count > 0:
            
            # Calculates net allwave radiation from the net solar radiation
            #incoming thermal/longwave radiation, and the snow surface
            # temperature.
            # replaces _net_rad()
            self.em.R_n = self.input1.S_n + (SNOW_EMISSIVITY * (self.input1.I_lw - STEF_BOLTZ * np.power(self.snow.T_s_0, 4)))
            
            # calculate H & L_v_E (and E as well)
            self.h_le()
            
            # calculate G & G_0 (conduction/diffusion heat xfr)
            if self.snow.layer_count == 1:
                self.g_soil('surface')
                self.em.G_0 = self.em.G
            else:  # layer_count == 2
                self.g_soil ('lower')
                self.g_snow()
                
            # calculate advection
            self.advec()
            
            # sum E.B. terms
            # surface energy budget
            self.em.delta_Q_0 = self.em.R_n + self.em.H + self.em.L_v_E + self.em.G_0 + self.em.M
            
            # total snowpck energy budget
            if self.snow.layer_count == 1:
                self.em.delta_Q = self.em.delta_Q_0
            else:
                self.em.delta_Q = self.em.delta_Q_0 + self.em.G - self.em.G_0
                
        else:
            self.em.R_n = 0
            self.em.H = 0
            self.em.L_v_E = 0
            self.em.E = 0
            self.em.G = 0
            self.em.G_0 = 0
            self.em.delta_Q = 0
            self.em.delta_Q_0 = 0
        
               
                       
    def advec(self):
        """
        This routine calculates the advected energy for a 2-layer snowcover
        if there's precipitation for the current timestep.
        """    
        
        if self.precip_now:
            
            M = self.heat_stor(CP_WATER(self.precip.T_rain), \
                               self.precip_info[self.tstep_level].m_rain, \
                               self.precip.T_rain - self.snow.T_s_0) + \
                 self.heat_stor(CP_ICE(self.precip.T_snow), \
                                self.precip_info[self.tstep_level].m_snow, \
                                self.precip.T_snow - self.snow.T_s_0)
                                
            M /= self.time_step
            
        else:
            M = 0
            
        self.em.M = M
        
        

    def g_soil(self, layer):
        """
        conduction heat flow between snow and soil
        """
        
        if layer == 'surface':
            tsno = self.snow.T_s_0
            ds = self.snow.z_s_0
            
        else:
            tsno = self.snow.T_s_l
            ds = self.snow.z_s_l
            
        if tsno > FREEZE:
            warnings.warn('g_soil: tsno = %8.2f; set to %8.2f\n' % (tsno, FREEZE))
            tsno = FREEZE
        
        # set effective soil conductivity
        # /***    changed to KT_MOISTSAND by D. Marks, NWRC, 09/30/2003    ***/
        # /***    based on heat flux data from RMSP            ***/
        # /***    note: Kt should be passed as an argument        ***/
        # /***    k_g = efcon(KT_WETSAND, tg, pa);            ***/
        k_g = libsnobal.efcon(KT_MOISTSAND, self.input1.T_g, self.P_a)
        
        # calculate G    
        # set snow conductivity
        kcs = KTS(self.snow.rho)
        k_s = libsnobal.efcon(kcs, tsno, self.P_a)
        g = libsnobal.ssxfr(k_s, k_g, tsno, self.input1.T_g, ds, self.z_g)
        
        self.em.G = g
        
    def g_snow(self):
        """
        conduction heat flow between snow layers
        """
        
        # calculate g
        if self.snow.T_s_0 == self.snow.T_s_l:
            g = 0
        else:
            kcs1 = KTS(self.snow.rho)
            kcs2 = KTS(self.snow.rho)
            k_s1 = libsnobal.efcon(kcs1, self.snow.T_s_0, self.P_a)
            k_s2 = libsnobal.efcon(kcs2, self.snow.T_s_l, self.P_a)
            
            g = libsnobal.ssxfr(k_s1, k_s2, self.snow.T_s_0, self.snow.T_s_l, self.snow.z_s_0, self.snow.z_s_l)
        
        self.em.G_0 = g
            
        
    def h_le(self):
        """
        Calculates point turbulent transfer (H and L_v_E) for a 2-layer snowcover
        """
        
#         if self.snow.T_s_0 < 0:
#             raise Exception('T_s_0 is below 0 K')
#         if self.input1.T_a < 0:
#             raise Exception('T_a is below 0 K')
        
        # calculate saturation vapor pressure
        e_s = libsnobal.sati(self.snow.T_s_0)
        
        # error check for bad vapor pressures
        sat_vp = libsnobal.sati(self.input1.T_a)
        
        if self.input1.e_a > sat_vp:
            self.input1.e_a = sat_vp
            
        # determine relative measurement heights
        if self.relative_hts:
            rel_z_T = self.z_T
            rel_z_u = self.z_u
        else:
            rel_z_T = self.z_T - self.snow.z_s
            rel_z_u = self.z_u - self.snow.z_s
        
        # calculate H & L_v_E
        H, L_v_E, E, status = libsnobal.hle1(self.P_a, self.input1.T_a, self.snow.T_s_0, rel_z_T, \
                                             self.input1.e_a, e_s, rel_z_T, self.input1.u, rel_z_u, self.z_0)
        if status != 0:
            raise Exception("hle1 did not converge\nP_a %f, T_a %f, T_s_0 %f\nrelative z_T %f, e_a %f, e_s %f\nu %f, relative z_u %f, z_0 %f\n" % \
                            (self.P_a, self.input1.T_a, self.snow.T_s_0, rel_z_T, \
                            self.input1.e_a, e_s, self.input1.u, rel_z_u, self.z_0))
            
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
        
        if self.snow.layer_count == 0:
            return False
        if self.snow.layer_count == 1:
            return self.snow.m_s < threshold
        else:
            return (self.snow.m_s_0 < threshold) or (self.snow.m_s_l < threshold)
            
           
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
            
            self.current_time = self.time_s
        
        
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
    
    
    def init_precip(self):
        """
        Returns a Pandas series that will contain all information about
        current time step precip
            
        """
                
#         return pd.Series(index=['m_pp','m_snow','m_rain','z_snow'])
        return Map({key: 0.0 for key in ['m_pp','m_snow','m_rain','z_snow']})
                
    
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
                'cc_s', 'cc_s_0', 'cc_s_l',
                'm_s', 'm_s_0', 'm_s_l',
                'h2o', 'h2o_max', 'h2o_total', 'h2o_vol']
               
        # set the values from the initial snow properties
        if from_record:
            
            # create an empty dataframe
#             self.snow = pd.Series(index=cols)
            
            # create a dict instead
            self.snow = Map({key: 0.0 for key in cols})
        
            self.snow.z_s = self.snow_records.iloc[self.snow_prop_index].z_s
            self.snow.rho = self.snow_records.iloc[self.snow_prop_index].rho
            self.snow.T_s_0 = self.snow_records.iloc[self.snow_prop_index].T_s_0
            self.snow.T_s = self.snow_records.iloc[self.snow_prop_index].T_s
            self.snow.h2o_sat = self.snow_records.iloc[self.snow_prop_index].h2o_sat
            self.snow.max_h2o_vol = self.params['max_h2o_vol']
        
        else:
            try:
                getattr(self, 'snow')
            except:
                raise Exception('The snow has not been intitialized')
        
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
                self.snow[col] = 0        
        
            # Note: Snow temperatures are set to MIN_SNOW_TEMP
            # (as degrees K) instead of 0 K to keep quantization
            # range in output image smaller.
            for col in ['T_s', 'T_s_0', 'T_s_l']:
                self.snow[col] = MIN_SNOW_TEMP + FREEZE
         
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
                'E_s_sum'
               ]
        
#         self.em = pd.Series(data=np.zeros(len(col)), index=col)
        self.em = Map({key: 0.0 for key in col})
        



#     @profile
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
            z_s = z_s_0
            
        else:
            # enough depth for both layers
            
            layer_count = 2
            z_s_0 = self.params['max_z_s_0']
            z_s_l = self.snow.z_s - z_s_0
            z_s = z_s_0 + z_s_l # not really needed but needed for below
            
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
        
        
#     @profile
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
            cc = self.heat_stor(CP_ICE(temp), mass, temp-FREEZE)
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
            
            # time
            self.params['out_file'].write('%g' % curr_time_hrs)
            
            # energy budget terms
            self.params['out_file'].write(" %.1f %.1f %.1f %.1f %.1f %.1f" % \
                    (self.em.R_n_bar, self.em.H_bar, self.em.L_v_E_bar, \
                    self.em.G_bar, self.em.M_bar, self.em.delta_Q_bar))

            # layer terms
            self.params['out_file'].write(" %.1f %.1f" % \
                    (self.em.G_0_bar, self.em.delta_Q_0_bar))
    
            # heat storage and mass changes
            self.params['out_file'].write(" %.6e %.6e %.6e" % \
                    (self.snow.cc_s_0, self.snow.cc_s_l, self.snow.cc_s))
            self.params['out_file'].write(" %.5f %.5f %.5f" % \
                    (self.em.E_s_sum, self.em.melt_sum, self.em.ro_pred_sum))
    
#             # runoff error if data included */
#             if (ro_data)
#                 fprintf(out, " %.1f",
#                         (ro_pred_sum - (ro * time_since_out)))
    
            # sno properties */
            self.params['out_file'].write(" %.3f %.3f %.3f %.1f" % \
                    (self.snow.z_s_0, self.snow.z_s_l, self.snow.z_s, self.snow.rho))
            self.params['out_file'].write(" %.1f %.1f %.1f %.1f" % \
                    (self.snow.m_s_0, self.snow.m_s_l, self.snow.m_s, self.snow.h2o))
            if self.params['temps_in_C']:
                self.params['out_file'].write(" %.2f %.2f %.2f\n" % 
                        (K_TO_C(self.snow.T_s_0), K_TO_C(self.snow.T_s_l), K_TO_C(self.snow.T_s)))
            else:
                self.params['out_file'].write(" %.2f %.2f %.2f\n" % \
                        (self.snow.T_s_0, self.snow.T_s_l, self.snow.T_s))
                
                
class Map(dict):
    """
    Example:
    m = Map({'first_name': 'Eduardo'}, last_name='Pool', age=24, sports=['Soccer'])
    """
    def __init__(self, *args, **kwargs):
        super(Map, self).__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for k, v in arg.iteritems():
                    self[k] = v

        if kwargs:
            for k, v in kwargs.iteritems():
                self[k] = v

    def __getattr__(self, attr):
        return self.get(attr)

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super(Map, self).__setitem__(key, value)
        self.__dict__.update({key: value})

    def __delattr__(self, item):
        self.__delitem__(item)

    def __delitem__(self, key):
        super(Map, self).__delitem__(key)
        del self.__dict__[key]            
        