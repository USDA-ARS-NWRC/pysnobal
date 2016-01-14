'''
20160112 Scott Havens

Plot the original Snobal vs pySnobal
'''


import numpy as np
import pandas as pd
# from mpl_toolkits.axes_grid1 import host_subplot
import matplotlib.pyplot as plt
import os


#------------------------------------------------------------------------------ 
# read the input and output files

output_label = np.array(['time_s','R_n','H','L_v_E','G','M','delta_Q','G_0','delta_Q_0',
            'cc_s_0','cc_s_l','cc_s','E_s','melt','ro_predict','z_s_0','z_s_l',
            'z_s','rho','m_s_0','m_s_l','m_s','h2o','T_s_0','T_s_l','T_s'])

# org = np.loadtxt('snobal.output.all')
# new = np.loadtxt('snobal.out')
org = pd.read_csv('snobal.v1', sep=' ', index_col=[0], names=output_label)
new = pd.read_csv('snobal.out', sep=' ', index_col=[0], names=output_label)

d = org - new

#------------------------------------------------------------------------------ 
# labels and such

data_label = np.array(['S_n', 'I_lw', 'T_a', 'e_a', 'u', 'T_g'])
ppt_label = np.array(['m_ppt','%_snow','rho_snow','T_pp'])


output_em = np.array([1,2,3,4,5,6,7,8,11,12,13])-1
output_snow = np.array([17,18,21,23,24,25,16,22])-1

#------------------------------------------------------------------------------ 
# plot the original outputs
f, axo = plt.subplots(3, 2, sharex=True)

org.plot(y=output_em, ax=axo[0][0], ylim=(-1500,1000))
axo[0][0].set_title('EM original')

org.plot(y=output_snow, ax=axo[0][1])
axo[0][1].set_title('SNOW original')

new.plot(y=output_em, ax=axo[1][0], ylim=(-1500,1000))
axo[1][0].set_title('EM new')

new.plot(y=output_snow, ax=axo[1][1])
axo[1][1].set_title('SNOW new')

d.plot(y=output_em, ax=axo[2][0], ylim=(-1500,1000))
axo[2][0].set_title('EM diff')

d.plot(y=output_snow, ax=axo[2][1])
axo[2][1].set_title('SNOW diff')

plt.show()


# # axo[0][0].plot(org_time, org[:,output_em])
# # axo[0][0].legend(output_label[output_em], loc=2)
# # axo[0][0].set_ylim([-1500,1500])
# # axo[0][0].set_title('EM original')
# 
# axo[0][1].plot(org_time, org[:,output_snow])
# # axo[0][1].legend(output_label[output_snow], loc=2)
# axo[0][1].set_title('SNOW original')
# 
# #------------------------------------------------------------------------------ 
# # plot the new outputs
# 
# axo[1][0].plot(new_time, new[:,output_em])
# # axo[1][0].legend(output_label[output_em], loc=2)
# axo[1][0].set_ylim([-1500,1500])
# axo[1][0].set_title('EM new')
# 
# axo[1][1].plot(new_time, new[:,output_snow])
# # axo[1][1].legend(output_label[output_snow], loc=2)
# axo[1][1].set_title('SNOW new')


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


