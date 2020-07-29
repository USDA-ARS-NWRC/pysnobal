from pysnobal.point.libsnobal import sati


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
        self.cc_s = init  # could be moved to a property that is calculated when needed
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

    @property
    def T_s_0(self):
        return self._T_s_0

    @T_s_0.setter
    def T_s_0(self, var):
        self._T_s_0 = var
        self.__es_s_0 = False

    @property
    def e_s_0(self):
        """Calculate the saturation vapor pressure over ice for
        the surface layer.

        Returns:
            float: saturation vapor pressure over ice
        """
        if not self.__es_s_0:
            self._es_s_0 = sati(self.T_s_0)
            self.__es_s_0 = True
        return self._es_s_0

    @property
    def T_s_l(self):
        return self._T_s_l

    @T_s_l.setter
    def T_s_l(self, var):
        self._T_s_l = var
        self.__es_s_l = False

    @property
    def e_s_l(self):
        """Calculate the saturation vapor pressure over ice for
        the lower layer.

        Returns:
            float: saturation vapor pressure over ice
        """
        if not self.__es_s_l:
            self._es_s_l = sati(self.T_s_l)
            self.__es_s_l = True
        return self._es_s_l
