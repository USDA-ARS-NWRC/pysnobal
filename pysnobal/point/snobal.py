"""
Class snobal() that will hold all the modeling components

20160109 Scott Havens
"""

import warnings
from copy import copy

import numpy as np

from pysnobal.core.constants import (
    DATA_TSTEP, FREEZE, GRAVITY, KT_MOISTSAND, MAX_SNOW_DENSITY, MIN_SNOW_TEMP,
    MOL_AIR, RHO_W0, SEA_LEVEL, SMALL_TSTEP, SNOW_EMISSIVITY, STD_AIRTMP,
    STD_LAPSE, STEF_BOLTZ, VAP_SUB)
from pysnobal.core.functions import (
    cp_ice, cp_water, diffusion_coef, gas_density, h2o_left, melt,
    time_average, vapor_flux)
from pysnobal.point import InputDeltas, SnowState, libsnobal


class Snobal(object):
    """

    To-do
    self.snow is the working data frame only containing one record
    self.output will be a dataframe that self.snow gets added

    """

    def __init__(self, params, tstep_info, snow_prop, meas_heights,
                 output_timesteps=None):
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
        self.start_date = self.params['start_date']
        self.current_datetime = self.params['start_date']

        self.output_timesteps = output_timesteps
        self.output_divided = False
        if self.output_timesteps is None:
            self.output_divided = True

        self.P_a = libsnobal.hysat(
            SEA_LEVEL,
            STD_AIRTMP,
            STD_LAPSE,
            self.elevation / 1000.0,
            GRAVITY,
            MOL_AIR)

        # get the intial snowcover properties
        self.snow_records = snow_prop
        self.get_sn_rec(True)

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
            input2: second timestep dict

            inputs contain all forcing data:
                ['S_n', 'I_lw', 'T_a', 'e_a', 'u', 'T_g','m_pp',
                    'percent_snow', 'rho_snow', 'T_pp']


        """

        # store the inputs for later
        self.input1 = input1
        self.input2 = input2

        # Compute deltas for the climate input parameters over the
        # data timestep.
        self.input_deltas = InputDeltas(
            self.input1, self.input2, self.tstep_info).calculate()

        # Divide the data timestep into normal run timesteps.
        # keeps track of what time step level the model is on
        self.current_level = DATA_TSTEP
        self.divided_step = False
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

        # Fetch the record for the timestep at the next level.
        self.next_level = self.current_level + 1
        next_lvl_tstep = self.tstep_info[self.next_level]

        if self.input1.precip_now and next_lvl_tstep['level'] > 1:
            self.input1.update_precip_deltas(
                self.input_deltas[next_lvl_tstep['level']])

        # For each the new smaller timestep, either subdivide them if
        # below their mass threshold, or run the model for them.
        interval = next_lvl_tstep['intervals']
        for i in range(interval):

            if (self.next_level != SMALL_TSTEP) and \
                    (self.below_thold(next_lvl_tstep['threshold'])):
                # increment the level number
                self.current_level = copy(self.next_level)
                self.divided_step = True
                if not self.divide_tstep():
                    return False
            else:
                if not self.do_tstep(next_lvl_tstep):
                    return False

        self.current_level -= 1
        self.next_level -= 1

        return True

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

        self.time_step = tstep['time_step']
        self.tstep_level = tstep['level']

        # get the current time step precip
        if self.input1.precip_now:
            self.input1.update_precip_deltas(
                self.input_deltas[self.tstep_level])

        self.snow_state.set_zeros('h2o_total')

        # is there snowcover?
        self.snowcover = self.snow_state.layer_count > 0

        # Calculate energy transfer terms
        self.energy_balance()

        # Adjust mass and calculate runoff
        self.mass_bal()

        # Update the averages for the energy terms and the totals for mass
        # changes since the last output.
        # TODO move this to SnowState
        if self.time_since_out > 0:
            self.snow_state.R_n_bar = time_average(
                self.snow_state.R_n_bar,
                self.time_since_out,
                self.snow_state.R_n,
                self.time_step)
            self.snow_state.H_bar = time_average(
                self.snow_state.H_bar,
                self.time_since_out,
                self.snow_state.H,
                self.time_step)
            self.snow_state.L_v_E_bar = time_average(
                self.snow_state.L_v_E_bar,
                self.time_since_out,
                self.snow_state.L_v_E,
                self.time_step)
            self.snow_state.G_bar = time_average(
                self.snow_state.G_bar,
                self.time_since_out,
                self.snow_state.G,
                self.time_step)
            self.snow_state.M_bar = time_average(
                self.snow_state.M_bar,
                self.time_since_out,
                self.snow_state.M,
                self.time_step)
            self.snow_state.delta_Q_bar = time_average(
                self.snow_state.delta_Q_bar,
                self.time_since_out,
                self.snow_state.delta_Q,
                self.time_step)
            self.snow_state.G_0_bar = time_average(
                self.snow_state.G_0_bar,
                self.time_since_out,
                self.snow_state.G_0,
                self.time_step)
            self.snow_state.delta_Q_0_bar = time_average(
                self.snow_state.delta_Q_0_bar,
                self.time_since_out,
                self.snow_state.delta_Q_0,
                self.time_step)

            self.snow_state.E_s_sum += self.snow_state.E_s
            self.snow_state.melt_sum += self.snow_state.melt
            self.snow_state.ro_pred_sum += self.snow_state.ro_predict

            self.time_since_out += self.time_step

        else:
            self.snow_state.R_n_bar = self.snow_state.R_n
            self.snow_state.H_bar = self.snow_state.H
            self.snow_state.L_v_E_bar = self.snow_state.L_v_E
            self.snow_state.G_bar = self.snow_state.G
            self.snow_state.M_bar = self.snow_state.M
            self.snow_state.delta_Q_bar = self.snow_state.delta_Q
            self.snow_state.G_0_bar = self.snow_state.G_0
            self.snow_state.delta_Q_0_bar = self.snow_state.delta_Q_0

            self.snow_state.E_s_sum = self.snow_state.E_s
            self.snow_state.melt_sum = self.snow_state.melt
            self.snow_state.ro_pred_sum = self.snow_state.ro_predict

            self.time_since_out = self.time_step

        # increment time
        self.current_time = self.current_time + self.time_step
        self.current_datetime = self.current_datetime + \
            tstep['time_step_timedelta']

        # # output if on the whole timestep
        if self.output_divided or \
                self.current_datetime in self.output_timesteps:
            # print('output whole timestep {}'.format(self.current_datetime))
            self.output()

        # Update the model's input parameters
        # TODO move this to the input delta where it's just accessing that
        # interval value, then we're not adding on the fly which is known to
        # have complex issues in python
        # Also this doesn't need to happen for the normal timestep
        self.input1.add_deltas(self.input_deltas[tstep['level']])

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
            self.snow_state.adjust_layer_temps()

    def runoff(self):
        """
        Calculates runoff for point energy budget 2-layer snowmelt model
        """

        # If no snow on ground at start of timestep or no layers currently,
        # then all water (e.g., rain) is runoff.

        if (not self.snowcover) or (self.snow_state.layer_count == 0):
            self.snow_state.ro_predict = self.snow_state.h2o_total
            return

        # Determine the snow density without any water, and the maximum
        # liquid water the snow can hold.
        m_s_dry = self.snow_state.m_s - self.snow_state.h2o_total
        rho_dry = m_s_dry / self.snow_state.z_s
        self.snow_state.h2o_max = h2o_left(
            self.snow_state.z_s, rho_dry, self.snow_state.max_h2o_vol)

        # Determine runoff, and water left in the snow
        if self.snow_state.h2o_total > self.snow_state.h2o_max:
            self.snow_state.ro_predict = self.snow_state.h2o_total - \
                self.snow_state.h2o_max
            self.snow_state.h2o = self.snow_state.h2o_max
            self.snow_state.h2o_sat = 1.0
            self.snow_state.h2o_vol = self.snow_state.max_h2o_vol

            # Update the snowcover's mass for the loss of runoff.
            self.adj_snow(0.0, -self.snow_state.ro_predict)

        else:
            self.snow_state.ro_predict = 0.0
            self.snow_state.h2o = self.snow_state.h2o_total
            self.snow_state.h2o_sat = self.snow_state.h2o / \
                self.snow_state.h2o_max
            self.snow_state.h2o_vol = self.snow_state.h2o_sat * \
                self.snow_state.max_h2o_vol

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
        (max - current|                     *   *
                      |                *
                      |             *
            delta_rho |          *
            (kg/m^2)  |        *
                  A/2 + . . . *
                      |     * .
                      |   *   .
                      |  *    .
                      | *     .
                      |*      .
                    0 +-------+-----------------------------+      h2o_added
                      0    B: 5 %                 1.0

        """
        # Maximum density due to compaction by liquid H2O added (kg/m^2)
        MAX_DENSITY = 550

        # ratio where half the difference between maximum density and
        # current density is reached (ratio from 0.0 to 1.0).
        B = 0.4

        if (not self.snowcover) or (self.snow_state.rho > MAX_DENSITY):
            return

        A = MAX_DENSITY - self.snow_state.rho
        if self.input1.precip_now:
            h2o_added = (self.snow_state.melt + self.input1.m_rain) / \
                self.snow_state.m_s

        else:
            h2o_added = self.snow_state.melt / self.snow_state.m_s

        if h2o_added > 0.000001:
            self.snow_state.rho += A / (1 + B/h2o_added)

            # adjust the snowcover for this new density (_new_density function)
            self.snow_state.z_s = self.snow_state.m_s / self.snow_state.rho
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
            self.snow_state.E_s = 0.0
            return

        # Total mass change due to evap/cond at surface during timestep
        E_s_0 = self.snow_state.E * self.time_step

        # Adjust total h2o for evaporative losses
        prev_h2o_tot = self.snow_state.h2o_total

        if self.snow_state.h2o_total > 0.0:
            self.snow_state.h2o_total += (E_s_0 * VAP_SUB)
            if self.snow_state.h2o_total <= 0.0:
                self.snow_state.h2o_total = 0.0

        # Determine total mass change due to evap/cond at soil
        if self.snow_state.layer_count == 0:
            E_s_l = 0.0
        else:
            if self.snow_state.layer_count == 2:
                e_s_l = self.snow_state.e_s_l
                T_bar = (self.input1.T_g + self.snow_state.T_s_l) / 2.0

            else:  # layer_count == 1
                e_s_l = self.snow_state.e_s_0
                T_bar = (self.input1.T_g + self.snow_state.T_s_0) / 2.0

            q_s_l = libsnobal.spec_hum(e_s_l, self.P_a)
            q_g = libsnobal.spec_hum(self.input1.e_g, self.P_a)
            q_delta = q_g - q_s_l
            rho_air = gas_density(self.P_a, libsnobal.MOL_AIR, T_bar)
            k = diffusion_coef(self.P_a, T_bar)

            E_l = vapor_flux(rho_air, k, q_delta, self.z_g)

            # total mass of evap/cond for time step
            E_s_l = E_l * self.time_step

            # adjust h2o_total for evaporative losses
            if self.snow_state.h2o_total > 0.0:
                self.snow_state.h2o_total += (E_s_l * VAP_SUB)
                if self.snow_state.h2o_total <= 0.0:
                    self.snow_state.h2o_total = 0.0

        self.snow_state.E_s = E_s_0 + E_s_l

        # adj mass and depth for evap/cond
        if self.snow_state.layer_count > 0:
            delta_z = ((self.snow_state.E_s +
                        (prev_h2o_tot - self.snow_state.h2o_total)) /
                       self.snow_state.rho) / 2.0
            self.adj_snow(delta_z, self.snow_state.E_s)

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
            self.snow_state.melt = 0
            return

        # calculate melt or freezing, and adjust cold content

        # calculate surface melt
        # energy for surface melt
        Q_0 = (self.snow_state.delta_Q_0 * self.time_step) + \
            self.snow_state.cc_s_0

        if Q_0 > 0:
            self.snow_state.melt = melt(Q_0)
            self.snow_state.cc_s_0 = 0
        elif Q_0 == 0:
            self.snow_state.melt = 0
            self.snow_state.cc_s_0 = 0
        else:
            self.snow_state.melt = 0
            self.snow_state.cc_s_0 = Q_0

        # calculate lower layer melt
        if self.snow_state.layer_count == 2:
            Q_l = ((self.snow_state.G - self.snow_state.G_0) *
                   self.time_step) + self.snow_state.cc_s_l

            if Q_l > 0:
                self.snow_state.melt += melt(Q_l)
                self.snow_state.cc_s_l = 0
            elif Q_l == 0:
                self.snow_state.cc_s_l = 0
            else:
                self.snow_state.cc_s_l = Q_l

        else:
            # layer count = 1
            Q_l = 0

        self.snow_state.h2o_total += self.snow_state.melt

        # adjust layers for re-freezing
        # adjust surface layer
        h2o_refrozen = 0

        if self.snow_state.cc_s_0 < 0:
            # if liquid h2o present, calc refreezing and adj cc_s_0
            if self.snow_state.h2o_total > 0:
                Q_freeze = self.snow_state.h2o_total * \
                    (self.snow_state.z_s_0/self.snow_state.z_s) * \
                    libsnobal.lh_fus(FREEZE)
                Q_left = Q_0 + Q_freeze

                if Q_left <= 0:
                    h2o_refrozen = self.snow_state.h2o_total * \
                        (self.snow_state.z_s_0/self.snow_state.z_s)
                    self.snow_state.cc_s_0 = Q_left
                else:
                    h2o_refrozen = self.snow_state.h2o_total * \
                        (self.snow_state.z_s_0/self.snow_state.z_s) - \
                        melt(Q_left)
                    self.snow_state.cc_s_0 = 0

        # adjust lower layer for re-freezing
        if (self.snow_state.layer_count == 2) and \
                (self.snow_state.cc_s_l < 0.0):
            # if liquid h2o, calc re-freezing and adj cc_s_l
            if self.snow_state.h2o_total > 0.0:
                Q_freeze = self.snow_state.h2o_total * \
                    (self.snow_state.z_s_l/self.snow_state.z_s) * \
                    libsnobal.lh_fus(FREEZE)
                Q_left = Q_l + Q_freeze

                if Q_left <= 0.0:
                    h2o_refrozen += self.snow_state.h2o_total * \
                        (self.snow_state.z_s_l/self.snow_state.z_s)
                    self.snow_state.cc_s_l = Q_left
                else:
                    h2o_refrozen += (
                        (self.snow_state.h2o_total *
                         (self.snow_state.z_s_l/self.snow_state.z_s)) -
                        melt(Q_left))
                    self.snow_state.cc_s_l = 0.0

        # Note:  because of rounding errors, h2o_refrozen may not
        # be exactly the same as h2o_total.  Check for this
        # case, and if so, then just zero out h2o_total.
        if np.abs(self.snow_state.h2o_total - h2o_refrozen) <= 1e-8:
            self.snow_state.h2o_total = 0
        else:
            self.snow_state.h2o_total -= h2o_refrozen

        # adjust depth and density for melt
        if self.snow_state.melt > 0:
            self.adj_snow((-1)*self.snow_state.melt/self.snow_state.rho, 0)

        # set total cold content
        # TODO move this to snow state
        if self.snow_state.layer_count == 2:
            self.snow_state.cc_s = self.snow_state.cc_s_0 + \
                self.snow_state.cc_s_l
        elif self.snow_state.layer_count == 1:
            self.snow_state.cc_s = self.snow_state.cc_s_0

    def precip_event(self):
        """
        This routine processes a precipitation event, i.e., the current
        precip record, if there's one for the current timestep.  It
        determines if the precip is rain or snow which increases the
        snowcover.
        """

        if self.input1.precip_now:
            if self.snowcover:
                # Adjust snowcover's depth and mass by snowfall's
                # depth and the total precipitation mass.
                self.adj_snow(
                    self.input1.z_snow, self.input1.m_pp)

                # Determine the additional liquid water that's in
                # the snowfall, and then add its mass to liquid
                # water in the whole snowcover.
                h2o_vol_snow = self.input1.h2o_sat_snow * \
                    self.snow_state.max_h2o_vol
                self.snow_state.h2o += h2o_left(self.input1.z_snow,
                                                self.input1.rho_snow,
                                                h2o_vol_snow)

            elif self.input1.m_snow > 0:

                # set the values from the initial snow properties
                self.snow_state.z_s = self.input1.z_snow
                self.snow_state.rho = self.input1.rho_snow
                self.snow_state.T_s = self.input1.T_snow
                self.snow_state.T_s_0 = self.input1.T_snow
                self.snow_state.T_s_l = self.input1.T_snow
                self.snow_state.h2o_sat = self.input1.h2o_sat_snow

                self.init_snow()

            # Add rainfall and water in the snowcover to the total liquid water
            self.snow_state.h2o_total += self.snow_state.h2o + \
                self.input1.m_rain

        else:
            # Add water in the snowcover to total liquid water
            self.snow_state.h2o_total += self.snow_state.h2o

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
        self.snow_state.z_s += delta_z_s
        self.snow_state.m_s += delta_m_s

        if self.snow_state.z_s != 0.0:
            self.snow_state.rho = self.snow_state.m_s / self.snow_state.z_s
        else:
            self.snow_state.rho = 0

        # clip density at maximum density if necessary
        if self.snow_state.rho > MAX_SNOW_DENSITY:
            self.snow_state.rho = MAX_SNOW_DENSITY
            self.snow_state.z_s = self.snow_state.m_s / self.snow_state.rho
            self.adj_layers()

        else:
            # If a change in depth, adjust the layers' depths and masses
            if delta_z_s != 0.0:
                self.adj_layers()
            else:
                # Just change in the snowcover's mass, update the layer masses
                self.layer_mass()

    def time_compact(self):
        """
        This routine replaces the original simple gravity compaction routine
        which "aged" the snowcover by accounting for the compaction
        or densification by gravity as time passes.  The original routine
        relied on an approximation based on equations in Anderson (1976), which
        increased the snowcover's density using the following "half-saturation"
        function:

            rho(time) = A / (1 + B/time)

        Historically, precise snow density was not a concern, as long as mass
        and SWE were correct.  With the development of the ASO program
        providing time-series lidar snow depth images, and the 2017 snow season
        in the southern Sierra Nevada, with individual storms approaching 500mm
        of deposition in a few days, and upper elevation snow depths of greater
        than 10m, it became clear that a more robust density model was
        required.

        Snow Density:  Snow density is initially a function of the temperature
        of the ice particles (snow flakes) as they fall to the ground during a
        storm.  Under very cold conditions (-10 to -15 C), we can get snow as
        light as 50 kg m-3, which is really light powder and great skiing.
        As the ice particle temperature approaches 0C, new snow densities can
        be as high as 200 kg m-3, which is not light powder, but still pretty
        good skiing, unless there is a lot of it. There are several (4)
        processes that can increase snow density during and after the storm.
        Note that the largest and most rapid changes occur during or just
        following the storm.  Compaction (an increase in density without a
        change in mass) can be caused by:

        1) Destructive mechanical metamorphism (compaction due to wind -
        affects mainly low density near-surface or new snow)
        2) Destructive temperature metamorphism
        3) Pressure - compaction due to snow load or overburden (affects
        both new snow deposition, snow on the ground, including drifting)
        4) Addition of liquid water due to melt or
        Of these compaction processes, iSnobal accounts three - 2,3 & 4.
        This routine addresses 2, temperature metamorphism, and 3., overburden
        metamorphism. The 4th, addition of liquid water is accounted for in the
        routine _h2o_compact.c. We are using equations found in Anderson (1976)
        and Oleson, et al. (2013), who got his equations from Anderson in the
        first place.  (Oleson made it easier to figure out the units...)

        Though many dedicated individuals have worked on the density issue over
        over the last 60+ years, we have, in general, a group working in
        Saporro Japan to thank for most of what we know about snow density.
        Anderson got the data and basic fitting equations from careful field
        and cold room measurements made by Yosida (1963), Mellor (1964) and
        Kojima (1967). The equations we are using are based on those that they
        derived from data that they carefully collected over a number of years.
        It is noteworthy that while rapidly changing computers have made the
        kind of spatial modeling we are attempting possible, snow physics
        remains unchanged – and the field and labortory efforts and data
        fitting equations from more than half a century ago represent our
        best understanding of those physics.

        Tz = 0.0 (freezing temperature, C or K)
        Ts = snow or precipitation temperature (C or K)
        rho = intital snow density (kg/(m^3))
        SWE = snow mass (mm/(m^2))
        K = temperature metamorphism coef.
        rho_n = new snow density (kg/(m^3))
        zs = new snow depth (mm)

        Proportional Destructive Temperature Metamorphism (PTM):

        if (rho < 100)
            K = 1.0
        else
            K = exp(-0.046 * (rho - 100))

        PTM = 0.01 * K * exp(-0.04 * (Tz - Ts))

        Proportional Overburden Compaction (POC):

        POC = (0.026 * exp(-0.08 * (Tz - Ts)) * SWE * exp(-21.0 * rho))

        New snow density and depth

        rho_n = rho + ((PTM + POC) * rho)
        zs_n = SWE / rho_n
        """

        # Maximum density due to compaction by gravity (kg/m^2)
        RHO_MAX = 550

        # R = 48 # in the original but not used?
        R1 = 23.5
        R2 = 24.5

        SWE_MAX = 2000.0

        # If the snow is already at or above the maximum density due
        # compaction by gravity, then just leave.
        if (not self.snowcover) or (self.snow_state.rho >= RHO_MAX):
            return

        # Calculate rate which compaction will be applied per time step.
        # Rate will be adjusted as time step varies.
        if self.snow_state.m_s >= SWE_MAX:
            rate = 1.0
        else:
            rate = R1 * np.cos(np.pi * self.snow_state.m_s / SWE_MAX) + R2
            rate = rate / (self.time_step / 3600.0)

        # Proportional Destructive Temperature Metamorphism (d_rho_m)
        if self.snow_state.rho < 100:
            K = 1.0
        else:
            K = np.exp(-0.046 * (self.snow_state.rho - 100))

        d_rho_m = 0.01 * K * np.exp(-0.04 * (FREEZE - self.snow_state.T_s))
        d_rho_m = d_rho_m / rate

        # Proportional Overburden Compaction (d_rho_c)
        d_rho_c = 0.026 * np.exp(-0.08 * (FREEZE - self.snow_state.T_s)) * \
            self.snow_state.m_s * \
            np.exp(-21.0 * (self.snow_state.rho / RHO_W0))
        d_rho_c = d_rho_c / rate

        self.snow_state.rho = self.snow_state.rho + \
            ((d_rho_m + d_rho_c) * self.snow_state.rho)

        # Adjust the snowcover for this new density
        self.snow_state.z_s = self.snow_state.m_s / self.snow_state.rho
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

        prev_layer_count = self.snow_state.layer_count

        self.calc_layers()

        if self.snow_state.layer_count == 0:
            # 1 or 2 layers ---> 0
            self.snow_state.rho = 0

            # if mass > 0, then it must be velow threshold.
            # So turn this little bit of mass into water
            if self.snow_state.m_s > 0:
                self.snow_state.h2o_total += self.snow_state.m_s

            # set a bunch of values to 0
            index = ['h2o_vol', 'h2o', 'h2o_max', 'h2o_sat',
                     'm_s', 'cc_s', 'm_s_0', 'cc_s_0']
            self.snow_state.set_zeros(index)

            # Note: Snow temperatures are set to MIN_SNOW_TEMP
            # (as degrees K) instead of 0 K to keep quantization
            # range in output image smaller.

            self.snow_state.T_s = MIN_SNOW_TEMP + FREEZE
            self.snow_state.T_s_0 = MIN_SNOW_TEMP + FREEZE

            if prev_layer_count == 2:
                self.snow_state.m_s_l = 0
                self.snow_state.cc_s_l = 0
                self.snow_state.T_s_l = MIN_SNOW_TEMP + FREEZE

            self.snowcover = False

        else:
            self.layer_mass()

            if (prev_layer_count == 1) and (self.snow_state.layer_count == 2):
                # 1 layer --> 2 layers, add lower layer
                self.snow_state.T_s_l = self.snow_state.T_s
                self.snow_state.cc_s_l = self.cold_content(
                    self.snow_state.T_s_l, self.snow_state.m_s_l)

            elif (prev_layer_count == 2) and \
                    (self.snow_state.layer_count == 1):
                # 2 layers --> 1 layer, remove lower layer
                self.snow_state.T_s_l = MIN_SNOW_TEMP + FREEZE
                self.snow_state.cc_s_l = 0

    # @profile
    def energy_balance(self):
        """
        Calculates point energy budget for 2-layer snowcover.
        """

        # if there is a snowcover
        if self.snow_state.layer_count > 0:

            # Calculates net allwave radiation from the net solar radiation
            # incoming thermal/longwave radiation, and the snow surface
            # temperature.
            # replaces _net_rad()
            self.snow_state.R_n = self.input1.S_n + \
                (SNOW_EMISSIVITY * (
                    self.input1.I_lw -
                    STEF_BOLTZ * np.power(self.snow_state.T_s_0, 4)))

            # calculate H & L_v_E (and E as well)
            self.h_le()

            # calculate G & G_0 (conduction/diffusion heat xfr)
            if self.snow_state.layer_count == 1:
                self.g_soil('surface')
                self.snow_state.G_0 = self.snow_state.G
            else:  # layer_count == 2
                self.g_soil('lower')
                self.g_snow()

            # calculate advection
            self.advec()

            # sum E.B. terms
            # surface energy budget
            self.snow_state.delta_Q_0 = self.snow_state.R_n + \
                self.snow_state.H + \
                self.snow_state.L_v_E + \
                self.snow_state.G_0 + \
                self.snow_state.M

            # total snowpack energy budget
            if self.snow_state.layer_count == 1:
                self.snow_state.delta_Q = self.snow_state.delta_Q_0
            else:
                self.snow_state.delta_Q = self.snow_state.delta_Q_0 + \
                    self.snow_state.G - self.snow_state.G_0

        else:
            self.snow_state.R_n = 0
            self.snow_state.H = 0
            self.snow_state.L_v_E = 0
            self.snow_state.E = 0
            self.snow_state.G = 0
            self.snow_state.G_0 = 0
            self.snow_state.delta_Q = 0
            self.snow_state.delta_Q_0 = 0

    def advec(self):
        """
        This routine calculates the advected energy for a 2-layer snowcover
        if there's precipitation for the current timestep.
        """

        if self.input1.precip_now:

            M = self.heat_stor(cp_water(self.input1.T_rain),
                               self.input1.m_rain,
                               self.input1.T_rain - self.snow_state.T_s_0) + \
                self.heat_stor(cp_ice(self.input1.T_snow),
                               self.input1.m_snow,
                               self.input1.T_snow - self.snow_state.T_s_0)

            M /= self.time_step

        else:
            M = 0

        self.snow_state.M = M

    def g_soil(self, layer):
        """
        conduction heat flow between snow and soil
        """

        if layer == 'surface':
            tsno = self.snow_state.T_s_0
            ds = self.snow_state.z_s_0
            es_layer = self.snow_state.e_s_0

        else:
            tsno = self.snow_state.T_s_l
            ds = self.snow_state.z_s_l
            es_layer = self.snow_state.e_s_l

        if tsno > FREEZE:
            warnings.warn('g_soil: tsno = %8.2f; set to %8.2f\n' %
                          (tsno, FREEZE))
            tsno = FREEZE

        # set effective soil conductivity
        # /***    changed to KT_MOISTSAND by D. Marks, NWRC, 09/30/2003    ***/
        # /***    based on heat flux data from RMSP            ***/
        # /***    note: Kt should be passed as an argument        ***/
        # /***    k_g = efcon(KT_WETSAND, tg, pa);            ***/
        k_g = libsnobal.efcon(
            KT_MOISTSAND,
            self.input1.T_g,
            self.P_a,
            es_layer=self.input1.e_g)

        # calculate G
        # set snow conductivity
        kcs = self.snow_state.kts
        k_s = libsnobal.efcon(
            kcs,
            tsno,
            self.P_a,
            es_layer=es_layer)
        g = libsnobal.ssxfr(k_s, k_g, tsno, self.input1.T_g, ds, self.z_g)

        self.snow_state.G = g

    def g_snow(self):
        """
        conduction heat flow between snow layers
        """

        # calculate g
        if self.snow_state.T_s_0 == self.snow_state.T_s_l:
            g = 0
        else:
            kcs1 = self.snow_state.kts
            kcs2 = self.snow_state.kts
            k_s1 = libsnobal.efcon(
                kcs1,
                self.snow_state.T_s_0,
                self.P_a,
                es_layer=self.snow_state.e_s_0)

            k_s2 = libsnobal.efcon(
                kcs2,
                self.snow_state.T_s_l,
                self.P_a,
                es_layer=self.snow_state.e_s_l)

            g = libsnobal.ssxfr(
                k_s1,
                k_s2,
                self.snow_state.T_s_0,
                self.snow_state.T_s_l,
                self.snow_state.z_s_0,
                self.snow_state.z_s_l)

        self.snow_state.G_0 = g

    # @profile
    def h_le(self):
        """
        Calculates point turbulent transfer (H and L_v_E) for a
        2-layer snowcover
        """

        # error check for bad vapor pressures
        sat_vp = self.input1.sat_vp

        if self.input1.e_a > sat_vp:
            warnings.warn(
                """{} h_le: input vapor pressure is greater than """
                """saturation vapor pressure""".format(self.current_datetime))
            self.input1.e_a = sat_vp

        # determine relative measurement heights
        if self.relative_hts:
            rel_z_T = self.z_T
            rel_z_u = self.z_u
        else:
            rel_z_T = self.z_T - self.snow_state.z_s
            rel_z_u = self.z_u - self.snow_state.z_s

        # calculate H & L_v_E
        H, L_v_E, E, status = libsnobal.hle1(
            self.P_a,
            self.input1.T_a,
            self.snow_state.T_s_0,
            rel_z_T,
            self.input1.e_a,
            self.snow_state.e_s_0,
            rel_z_T,
            self.input1.u,
            rel_z_u,
            self.z_0
        )

        if status != 0:
            raise Exception(
                """hle1 did not converge \n"""
                """P_a={}, T_a={}, T_s_0={}\n"""
                """relative z_T={}, e_a={}, e_s={}\n"""
                """u={}, relative z_u={}, z_0={}\n""".format(
                    self.P_a, self.input1.T_a, self.snow_state.T_s_0, rel_z_T,
                    self.input1.e_a, self.snow_state.e_s_0, self.input1.u,
                    rel_z_u, self.z_0))

        self.snow_state.H = H
        self.snow_state.L_v_E = L_v_E
        self.snow_state.E = E

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

        if self.snow_state.layer_count == 0:
            return False
        if self.snow_state.layer_count == 1:
            return self.snow_state.m_s < threshold
        else:
            return (self.snow_state.m_s_0 < threshold) or \
                (self.snow_state.m_s_l < threshold)

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

            # keep track of which snow property to read
            self.snow_prop_index = 0

            self.time_s = 0
            self.curr_time_hrs = 0
            self.start_time = 0
            self.more_sn_recs = True
            self.z_0 = self.snow_records['z_0']

            self.current_time = self.time_s

        else:
            # haven't ever seen this used so not entirly sure now this works
            # increase the index if there are more than one record
            self.snow_prop_index += 1

        self.more_sn_recs = False

    def get_measurement_height_rec(self, first_rec=False):
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
            self.measurement_heights_index = 0

            # self.time_z = self.measurement_heights['name * HR_TO_SEC
            self.z_u = self.measurement_heights['z_u']
            self.z_T = self.measurement_heights['z_t']
            self.z_g = self.measurement_heights['z_g']

        else:
            self.measurement_heights_index += 1

        self.more_mh_recs = False

    def init_snow(self, from_record=False):
        """init_sno
        This routine initializes the properties for the snowcover.  It
        determines the number of layers, their individual properties,
        the cold content for the snowcover and its layers, and the
        snowcover's water content.

        Initialize all the following values
        h2o_sat, layer_count, max_h2o_vol, rho, T_s, T_s_0, T_s_l,
        z_s, cc_s, cc_s_0, cc_s_l, h2o, h2o_max, h2o_total, h2o_vol,
        m_s, m_s_0, m_s_l

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

        # set the values from the initial snow properties
        if from_record:

            self.snow_state = SnowState()

            self.snow_state.z_s = self.snow_records['z_s']
            self.snow_state.rho = self.snow_records['rho']
            self.snow_state.T_s_0 = self.snow_records['t_s_0']
            self.snow_state.T_s = self.snow_records['t_s']
            self.snow_state.h2o_sat = self.snow_records['h2o_sat']
            self.snow_state.max_h2o_vol = self.params['max_h2o_vol']

        elif not hasattr(self, 'snow_state'):
            raise Exception('The snow_state has not been intitialized')

        # initialize the snowpack
        self.snow_state.m_s = self.snow_state.rho * self.snow_state.z_s

        # determine the number of layers
        self.calc_layers()

        if self.snow_state.layer_count == 0:
            # If mass > 0, then it must be below threshold.
            # So turn this little bit of mass into water

            if self.snow_state.m_s > 0.0:
                self.snow_state.h2o_total += self.snow_state.m_s

            self.snow_state.set_zeros([
                'rho', 'm_s', 'cc_s', 'm_s_0', 'cc_s_0', 'm_s_l',
                'cc_s_l', 'h2o_vol', 'h2o', 'h2o_max', 'h2o_sat'
            ])

            # Note: Snow temperatures are set to MIN_SNOW_TEMP
            # (as degrees K) instead of 0 K to keep quantization
            # range in output image smaller.
            for col in ['T_s', 'T_s_0', 'T_s_l']:
                setattr(self.snow_state, col, MIN_SNOW_TEMP + FREEZE)

        else:
            # Compute the specific mass and cold content for each layer
            self.layer_mass()
            self.snow_state.cc_s_0 = self.cold_content(
                self.snow_state.T_s_0,
                self.snow_state.m_s_0)

            if self.snow_state.layer_count == 2:
                self.snow_state.cc_s_l = self.cold_content(
                    self.snow_state.T_s_l,
                    self.snow_state.m_s_l)
            else:
                self.snow_state.T_s_l = MIN_SNOW_TEMP + FREEZE
                self.snow_state.cc_s_l = 0

            # Compute liquid water content as volume ratio, and
            # snow density without water
            self.snow_state.h2o_vol = self.snow_state.h2o_sat * \
                self.snow_state.max_h2o_vol

            # Determine the maximum liquid water content (as specific mass)
            # and the actual liquid water content (as specific mass)
            self.snow_state.h2o_max = h2o_left(
                self.snow_state.z_s,
                self.snow_state.dry_snow_density,
                self.snow_state.max_h2o_vol)
            self.snow_state.h2o = self.snow_state.h2o_sat * \
                self.snow_state.h2o_max

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

        if self.snow_state.m_s <= self.tstep_info[SMALL_TSTEP]['threshold']:
            # less than minimum layer mass, so treat as no snowcover

            layer_count = 0
            z_s = z_s_0 = z_s_l = 0

        elif self.snow_state.z_s < self.params['max_z_s_0']:
            # not enough depth for surface layer and the lower layer,
            # so just 1 layer: surface layer

            layer_count = 1
            z_s_0 = self.snow_state.z_s
            z_s_l = 0
            z_s = z_s_0

        else:
            # enough depth for both layers

            layer_count = 2
            z_s_0 = self.params['max_z_s_0']
            z_s_l = self.snow_state.z_s - z_s_0
            z_s = z_s_0 + z_s_l  # not really needed but needed for below

            # However, make sure there's enough MASS for the lower
            # layer.  If not, then there's only 1 layer
            if z_s_l * self.snow_state.rho < self.tstep_info[SMALL_TSTEP]['threshold']:  # noqa
                layer_count = 1
                z_s_0 = self.snow_state.z_s
                z_s_l = 0

        self.snow_state.layer_count = layer_count
        self.snow_state.z_s = z_s
        self.snow_state.z_s_0 = z_s_0
        self.snow_state.z_s_l = z_s_l

    def layer_mass(self):
        """
        This routine computes the specific mass for each snow layer in
        the snowcover.  A layer's mass is based its depth and the
        average snowcover density.
        """

        if self.snow_state.layer_count == 0:
            self.snow_state.m_s_0 = 0
            self.snow_state.m_s_l = 0

        else:
            # layer count is 1 or 2
            self.snow_state.m_s_0 = self.snow_state.rho * self.snow_state.z_s_0

            if self.snow_state.layer_count == 2:
                self.snow_state.m_s_l = self.snow_state.rho * \
                    self.snow_state.z_s_l
            else:
                self.snow_state.m_s_l = 0

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
            cc = self.heat_stor(cp_ice(temp), mass, temp-FREEZE)
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

    def init_output(self):
        """Initialize the output, it will be a list of dictionaries that
        will be converted to a dataframe at the end
        """

        # sz = self.elevation.shape
        flds = ['rho', 'T_s_0', 'T_s_l', 'T_s',
                'cc_s_0', 'cc_s_l', 'cc_s', 'm_s', 'm_s_0', 'm_s_l', 'z_s',
                'z_s_0', 'z_s_l', 'h2o_sat', 'layer_count', 'h2o', 'h2o_max',
                'h2o_vol', 'h2o_total', 'R_n_bar', 'H_bar', 'L_v_E_bar',
                'G_bar', 'G_0_bar', 'M_bar', 'delta_Q_bar', 'delta_Q_0_bar',
                'E_s_sum', 'melt_sum', 'ro_pred_sum']
        s = {key: 0.0 for key in flds}  # the structure fields

        # # update the output rec with the initial snowpack state
        # for key, val in self.snowpack_state.items():
        #     if key in flds:
        #         s[key] = val

        # s['mask'] = np.ones(sz)
        self.output_rec = s

        # Create the output storage
        self.output_list = []

    def output(self):
        """
        Specify where the model output should go
        """

        c = {key: getattr(self.snow_state, key)
             for key in self.output_rec.keys()}
        c['date_time'] = self.current_datetime

        self.output_list.append(c)
        self.time_since_out = 0.0
