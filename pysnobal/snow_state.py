import numpy as np


class SnowState():

    def __init__(self):

        init = 0
        self.zeros = init

        # snowpack state variables
        self.h2o_sat = init
        self.layer_count = init
        self.max_h2o_vol = init
        self.rho = init
        self.T_s = init
        self.T_s_0 = init
        self.T_s_l = init
        self.z_s = init
        self.z_s_0 = init
        self.z_s_l = init
        self.cc_s = init
        self.cc_s_0 = init
        self.cc_s_l = init
        self.m_s = init
        self.m_s_0 = init
        self.m_s_l = init
        self.h2o = init
        self.h2o_max = init
        self.h2o_total = init
        self.h2o_vol = init

        # Snow energetics state variables
        # net allwave radiation (W/m^2)
        self.R_n = init

        # sensible heat xfr (W/m^2)
        self.H = init

        # latent heat xfr (W/m^2)
        self.L_v_E = init

        # heat xfr by conduction & diffusion from soil to snowcover (W/m^2)
        self.G = init

        # heat xfr by conduction & diffusion from soil or lower layer
        # to active layer (W/m^2)
        self.G_0 = init

        # advected heat from precip (W/m^2)
        self.M = init

        # change in snowcover's energy(W/m ^ 2)
        self.delta_Q = init

        # change in active layer's energy(W/m ^ 2)
        self.delta_Q_0 = init

        #   averages of energy balance vars since last output record
        self.R_n_bar = init
        self.H_bar = init
        self.L_v_E_bar = init
        self.G_bar = init
        self.G_0_bar = init
        self.M_bar = init
        self.delta_Q_bar = init
        self.delta_Q_0_bar = init

        # mass balance vars for current timestep
        # specific melt (kg/m^2 or m)
        self.melt = init

        # mass flux by evap into air from active layer (kg/m^2/s)
        self.E = init

        # mass of evap into air & soil from snowcover (kg/m^2)
        self.E_s = init

        # predicted specific runoff (m/sec)
        self.ro_predict = init

        # sums of mass balance vars since last output record
        self.melt_sum = init
        self.E_s_sum = init

    def set_zeros(self, fields):

        if isinstance(fields, str):
            fields = fields.split()

        for field in fields:
            setattr(self, field, self.zeros)

    # def calc_layers(self):
    #     """
    #     This routine determines the # of layers in the snowcover based its
    #     depth and mass.  Usually, there are are 2 layers: the surface (active)
    #     and the lower layer.  The depth of the surface layer is set to the
    #     maximum depth for the surface layer (variable "max_z_s_0").  The
    #     remaining depth constitutes the lower layer.  The routine checks
    #     to see if the mass of this lower layer is above the minimum threshold
    #     (i.e., the mass threshold for the small run timestep).  If not,
    #     the surface layer is the whole snowcover, and there's no lower
    #     layer.

    #     """

    #     if self.m_s <= self.tstep_info[SMALL_TSTEP]['threshold']:
    #         # less than minimum layer mass, so treat as no snowcover

    #         layer_count = 0
    #         z_s = z_s_0 = z_s_l = 0

    #     elif self.z_s < self.params['max_z_s_0']:
    #         # not enough depth for surface layer and the lower layer,
    #         # so just 1 layer: surface layer

    #         layer_count = 1
    #         z_s_0 = self.z_s
    #         z_s_l = 0
    #         z_s = z_s_0

    #     else:
    #         # enough depth for both layers

    #         layer_count = 2
    #         z_s_0 = self.params['max_z_s_0']
    #         z_s_l = self.z_s - z_s_0
    #         z_s = z_s_0 + z_s_l  # not really needed but needed for below

    #         # However, make sure there's enough MASS for the lower
    #         # layer.  If not, then there's only 1 layer
    #         if z_s_l * self.rho < self.tstep_info[SMALL_TSTEP]['threshold']:
    #             layer_count = 1
    #             z_s_0 = self.z_s
    #             z_s_l = 0

    #     self.layer_count = layer_count
    #     self.z_s = z_s
    #     self.z_s_0 = z_s_0
    #     self.z_s_l = z_s_l
