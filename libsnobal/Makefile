CATEGORY := snobal
DESCRIPTION := energy-balance snowmelt model

FUNCS := init_snow  do_data_tstep

PRIVATE_SRCS := _adj_layers.c _adj_snow.c _advec.c _below_thold.c \
		_calc_layers.c _cold_content.c _divide_tstep.c _do_tstep.c \
		_e_bal.c _evap_cond.c _h2o_compact.c _h_le.c _layer_mass.c \
		_mass_bal.c _net_rad.c _new_density.c _precip.c _runoff.c \
		_snowmelt.c _time_compact.c _vars.c \
		vars.c

include $(IPW)/make/funcCategory
