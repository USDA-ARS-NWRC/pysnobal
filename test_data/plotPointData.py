'''
20160112 Scott Havens

Plot the original Snobal vs pySnobal
'''


import numpy as np
# from mpl_toolkits.axes_grid1 import host_subplot
import matplotlib.pyplot as plt
import os


#------------------------------------------------------------------------------ 
# read the input and output files

org = np.loadtxt('snobal.output.all')
new = np.loadtxt('snobal.out')

new_time = new[:,0]
org_time = org[:,0]

ind = 

#------------------------------------------------------------------------------ 
# labels and such

data_label = np.array(['S_n', 'I_lw', 'T_a', 'e_a', 'u', 'T_g'])
ppt_label = np.array(['m_ppt','%_snow','rho_snow','T_pp'])
output_label = np.array(['time_s','R_n','H','L_v_E','G','M','delta_Q','G_0','delta_Q_0',
            'cc_s_0','cc_s_l','cc_s','E_s','melt','ro_predict','z_s_0','z_s_l',
            'z_s','rho','m_s_0','m_s_l','m_s','h2o','T_s_0','T_s_l','T_s'])

output_em = [1,2,3,4,5,6,7,8,11,12,13]
output_snow = [17,18,21,23,24,25,16,22]

#------------------------------------------------------------------------------ 
# plot the original outputs
f, axo = plt.subplots(2, 2, sharex=True)

axo[0][0].plot(org_time, org[:,output_em])
# axo[0][0].legend(output_label[output_em], loc=2)
axo[0][0].set_ylim([-1500,1500])
axo[0][0].set_title('EM original')

axo[0][1].plot(org_time, org[:,output_snow])
# axo[0][1].legend(output_label[output_snow], loc=2)
axo[0][1].set_title('SNOW original')

#------------------------------------------------------------------------------ 
# plot the new outputs

axo[1][0].plot(new_time, new[:,output_em])
# axo[1][0].legend(output_label[output_em], loc=2)
axo[1][0].set_ylim([-1500,1500])
axo[1][0].set_title('EM new')

axo[1][1].plot(new_time, new[:,output_snow])
# axo[1][1].legend(output_label[output_snow], loc=2)
axo[1][1].set_title('SNOW new')


#------------------------------------------------------------------------------ 
# plot the difference



# axo[2][0].plot(org[:,output_em]-new[:,output_em])
# # axo[2][0].legend(output_label[output_em], loc=2)
# axo[2][0].set_ylim([-1500,1500])
# axo[2][0].set_title('EM difference')
# 
# axo[2][1].plot(org[:,output_snow]-new[:,output_snow])
# # axo[2][1].legend(output_label[output_snow], loc=2)
# axo[2][1].set_title('SNOW difference')


plt.show()


