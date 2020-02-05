def SEC_TO_HR(x): return x / 3600.0


def someFunction(text):
    print('You passed this Python program {} from C! Congratulations!'.format(text))
    return 12345


def output_timestep_to_file(output_dict):
    """
    Output the model results to a file
    ** 
    This is a departure from Snobal that can print out the
    sub-time steps, this will only print out on the data tstep
    (for now) 
    **

    """

    # write out to a file
    with open('tests/test_data_point/snobal.pysnobal_c', 'a') as f:

        curr_time_hrs = SEC_TO_HR(output_dict['current_time'])

        # time
        f.write('%g,'.format(curr_time_hrs))

        # # energy budget terms
        # f.write("%.3f,%.3f,%.3f,%.3f,%.3f,%.3f," %
        #         (self.output_rec['R_n_bar'], self.output_rec['H_bar'], self.output_rec['L_v_E_bar'],
        #             self.output_rec['G_bar'], self.output_rec['M_bar'], self.output_rec['delta_Q_bar']))

        # # layer terms
        # f.write("%.3f,%.3f," %
        #         (self.output_rec['G_0_bar'], self.output_rec['delta_Q_0_bar']))

        # # heat storage and mass changes
        # f.write("%.9e,%.9e,%.9e," %
        #         (self.output_rec['cc_s_0'], self.output_rec['cc_s_l'], self.output_rec['cc_s']))
        # f.write("%.8f,%.8f,%.8f," %
        #         (self.output_rec['E_s_sum'], self.output_rec['melt_sum'], self.output_rec['ro_pred_sum']))

        # #             # runoff error if data included */
        # #             if (ro_data)
        # #                 fprintf(out, " %.3f",
        # #                         (ro_pred_sum - (ro * time_since_out)))

        # # sno properties */
        # f.write("%.6f,%.6f,%.6f,%.3f," %
        #         (self.output_rec['z_s_0'], self.output_rec['z_s_l'], self.output_rec['z_s'], self.output_rec['rho']))
        # f.write("%.3f,%.3f,%.3f,%.3f," %
        #         (self.output_rec['m_s_0'], self.output_rec['m_s_l'], self.output_rec['m_s'], self.output_rec['h2o']))
        # if self.params['temps_in_C']:
        #     f.write("%.5f,%.5f,%.5f\n" %
        #             (K_TO_C(self.output_rec['T_s_0']), K_TO_C(self.output_rec['T_s_l']), K_TO_C(self.output_rec['T_s'])))
        # else:
        #     f.write("%.5f,%.5f,%.5f\n" %
        #             (self.output_rec['T_s_0'], self.output_rec['T_s_l'], self.output_rec['T_s']))

        f.write("\n")
        # # reset the time since out
        # self.output_rec['time_since_out'] = 0
