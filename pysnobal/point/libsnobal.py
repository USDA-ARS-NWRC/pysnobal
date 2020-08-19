import math

import numpy as np

from pysnobal.core.constants import (BOIL, FREEZE, GRAVITY, LOG_SEA_LEVEL,
                                     MOL_AIR, VON_KARMAN)
from pysnobal.core.functions import (diffusion_coef, gas_density, lh_fus,
                                     lh_sub, lh_vap, mix_ratio, spec_hum,
                                     virtual_temperature)

# specific heat of air at constant pressure (J / kg / deg)
CP_AIR = 1.005e3
DALR = GRAVITY / CP_AIR  # dry adiabatic lapse rate (deg / m)

AH = 1.0                # ratio sensible/momentum phi func
AV = 1.0                # ratio latent/momentum phi func
ITMAX = 50              # max # iterations allowed
PAESCHKE = 7.35         # Paeschke's const (eq. 5.3)
THRESH = 1.e-5          # convergence threshold

SM = 0
SH = 1
SV = 2
BETA_S = 5.2
BETA_U = 16.0

LOG_10 = math.log(10.0)  # log(10)
PI_2 = np.pi / 2.


def satw_np(tk):
    """
    Saturation vapor pressure of water. from IPW satw
    Args:
        tk: temperature in Kelvin
    Returns:
        saturated vapor pressure over water
    """

    # remove bad values
    tk[tk < 0] = np.nan

    btk = BOIL/tk
    x = -7.90298*(btk - 1.0) + 5.02808*np.log(btk)/LOG_10 - \
        1.3816e-7*(np.power(10.0, 1.1344e1*(1.0 - tk/BOIL))-1.) + \
        8.1328e-3*(np.power(10.0, -3.49149*(btk - 1.0)) - 1.0) + \
        LOG_SEA_LEVEL/LOG_10

    x = np.power(10.0, x)

    return x


def sati_np(tk):
    """
    saturation vapor pressure over ice. From IPW sati
    Args:
        tk: temperature in Kelvin
    Returns:
        saturated vapor pressure over ice
    20151027 Scott Havens
    """

    # remove bad values
    tk[tk < 0] = np.nan

    # preallocate
    x = np.empty(tk.shape)

    # vapor above freezing
    ind = tk > FREEZE
    x[ind] = satw_np(tk[ind])

    # vapor below freezing
    x[~ind] = 100.0 * np.power(10.0, -9.09718*((FREEZE/tk[~ind]) - 1.0) -
                               3.56654*np.log(FREEZE/tk[~ind])/LOG_10 +
                               8.76793e-1*(1.0 - (tk[~ind]/FREEZE)) +
                               np.log(6.1071)/LOG_10)

    return x


# @profile
def satw(tk):
    '''
    Saturation vapor pressure of water. from IPW satw but for a single value
    20160112 Scott Havens
    '''

    if tk < 0:
        raise ValueError('tk < 0')

    btk = BOIL/tk
    x = -7.90298 * (btk - 1.0) + \
        5.02808 * math.log(btk) / LOG_10 - \
        1.3816e-7 * (math.pow(10.0, 11.344 * (1.0 - tk/BOIL)) - 1.) + \
        8.1328e-3 * (math.pow(10.0, -3.49149 * (btk - 1.0)) - 1.0) + \
        LOG_SEA_LEVEL / LOG_10

    x = math.pow(10.0, x)

    return x


# @profile
def sati(tk):
    '''
    saturation vapor pressure over ice. From IPW sati but for a single value
    20151027 Scott Havens
    '''

    if tk < 0:
        raise ValueError('tk < 0')

    # vapor above freezing
    if tk > FREEZE:
        x = satw(tk)

    else:
        # vapor below freezing
        ftk = FREEZE/tk
        x = 100.0 * math.pow(
            10.0,
            -9.09718 * (ftk - 1.0) -
            3.56654 * math.log(ftk) / LOG_10 +
            8.76793e-1 * (1.0 - tk/FREEZE) +
            math.log(6.1071) / LOG_10)

    return x


# @profile
def psi(zeta, code):
    """
    psi-functions

    The original clone from the IPW psi function

    code =   SM    momentum
             SH    sensible heat flux
             SV    latent heat flux
    """

    if code not in ['SM', 'SH', 'SV']:
        raise ValueError("psi-function code not of these: SM, SH, SV")

    if zeta > 0.0:        # stable
        result = -BETA_S * min(zeta, 1.0)

    elif zeta < 0.0:  # unstable

        x = math.sqrt(math.sqrt(1.0 - BETA_U * zeta))

        if code == 'SM':
            result = 2. * math.log((1. + x)/2.) + \
                math.log((1. + x * x)/2.) - \
                2. * math.atan(x) + \
                PI_2

        else:
            result = 2. * math.log((1. + x * x)/2.)

    else:  # neutral
        result = 0

    return result


def psi_momentum(x):
    """psi function for momentum (SM)

    Args:
        x (float): x value

    Returns:
        float: momentum value
    """
    return 2. * math.log((1. + x)/2.) + \
        math.log((1. + x * x)/2.) - \
        2. * math.atan(x) + \
        PI_2


def psi_heat_flux(x):
    """psi function for sensible heat flux (SH) and
    latent heat flux (SV)

    Args:
        x (float): x value

    Returns:
        float: heat flux value
    """
    return 2. * math.log((1. + x * x)/2.)


# @profile
def psi_func(zeta, func):
    """Calculate the psi value

    Args:
        zeta (float): zeta value
        func (instance): function name to call

    Returns:
        float: psi value
    """

    if zeta > 0.0:        # stable
        result = -BETA_S * min(zeta, 1.0)

    elif zeta < 0.0:  # unstable
        x = math.sqrt(math.sqrt(1.0 - BETA_U * zeta))
        result = func(x)

    else:  # neutral
        result = 0

    return result


# @profile
def hle1(press, air_temp, surface_temp, za, ea, es, zq, wind_speed, zu, z0,
         init_ustar=None, init_factor=None):
    """
    computes sensible and latent heat flux and mass flux given
    measurements of temperature and specific humidity at surface
    and one height, wind speed at one height, and roughness
    length.  The temperature, humidity, and wind speed measurements
    need not all be at the same height.

    Args
        press air pressure (Pa)
        air_temp air temperature (K) at height za
        surface_temp surface temperature (K)
        za height of air temp measurement (m)
        ea vapor pressure (Pa) at height zq
        es vapor pressure (Pa) at surface
        zq height of spec hum measurement (m)
        wind_speed wind speed (m/s) at height zu
        zu height of wind speed measurement (m)
        z0 roughness length (m)

    Outputs
        h sens heat flux (+ to surf) (W/m^2)
        le latent heat flux (+ to surf) (W/m^2)
        e mass flux (+ to surf) (kg/m^2/s)
        status status of convergence
            0      successful calculation
            -1      no convergence

    """

    # Check for bad inputs
    # heights must be positive
    if (z0 <= 0) or (zq <= z0) or (zu <= z0) or (za <= z0):
        raise ValueError(
            "height not positive z0={}, zq={}, zu={}, za={}".format(
                z0, zq, zu, za))

    # temperatures are Kelvin
    if (air_temp <= 0) or (surface_temp <= 0):
        raise ValueError("temps not K air_temp=%f\tsurface_temp=%f" %
                         (air_temp, surface_temp))

    # pressures must be positive
    if (ea <= 0) or (es <= 0) or (press <= 0) \
            or (ea >= press) or (es >= press):
        raise ValueError("press < 0 ea=%f\tes=%f\tpress=%f" % (ea, es, press))

    # vapor pressures can't exceed saturation
    # this should never error if vapor pressure was derived from air temp
    ea_w = satw(air_temp)
    if ea - 25.0 > ea_w:
        raise ValueError(
            """vapor pressure exceeds air saturation: """
            """air_temp={}, ea={}, ea_sat={}""".format(
                air_temp, ea, sati(air_temp)
            ))

    if ea > ea_w:
        ea = ea_w

    # displacement plane height, eq. 5.3 & 5.4
    d0 = 2 * PAESCHKE * z0 / 3

    # constant log expressions
    ltsm = math.log((zu - d0) / z0)
    ltsh = math.log((za - d0) / z0)
    ltsv = math.log((zq - d0) / z0)

    # convert vapor pressures to specific humidities
    q_diff = spec_hum(ea, press) - spec_hum(es, press)

    # convert temperature to potential temperature
    air_temp += DALR * za
    t_diff = air_temp - surface_temp

    # air density at press, virtual temp of geometric mean
    # of air and surface
    dens = gas_density(
        press,
        MOL_AIR,
        virtual_temperature(
            math.sqrt(air_temp * surface_temp),
            math.sqrt(ea * es),
            press)
    )

    # if neutral stability ignore starting values and calculate
    # when would this happen? With floating point precision, this will almost
    # never happen. Even if it did, the solution would be found quickly?
    if air_temp == surface_temp:

        # starting value, assume neutral stability, so psi-functions
        # are all zero
        ustar = VON_KARMAN * wind_speed / ltsm
        factor = VON_KARMAN * ustar * dens
        e = q_diff * factor * AV / ltsv
        h = t_diff * factor * CP_AIR * AH / ltsh
        it = 0

    # if not neutral stability, iterate on Obukhov stability
    # length to find solution
    else:
        # it = 0
        lo = 1e500

        # initialize ustar and factor
        if init_ustar is None:
            ustar = VON_KARMAN * wind_speed / ltsm
        else:
            ustar = init_ustar

        if init_factor is None:
            factor = VON_KARMAN * ustar * dens
        else:
            factor = init_factor

        e = q_diff * factor * AV / ltsv
        h = t_diff * factor * CP_AIR * AH / ltsh

        for it in range(ITMAX):
            last = lo

            # Eq 4.25, but no minus sign as we define
            # positive H as toward surface
            lo = ustar * ustar * ustar * dens / \
                (VON_KARMAN * GRAVITY * (h/(air_temp * CP_AIR) + 0.61 * e))

            # friction velocity, eq. 4.34'
            ustar = VON_KARMAN * wind_speed / \
                (ltsm - psi_func(zu/lo, psi_momentum))

            # evaporative flux, eq. 4.33'
            factor = VON_KARMAN * ustar * dens
            e = q_diff * factor * AV / (ltsv - psi_func(zq/lo, psi_heat_flux))

            # sensible heat flux, eq. 4.35'
            # with sign reversed
            h = t_diff * factor * AH * CP_AIR / \
                (ltsh - psi_func(za/lo, psi_heat_flux))

            diff = last - lo

            if (abs(diff) < THRESH) and (abs(diff/lo) < THRESH):
                break

    ier = -1 if (it >= ITMAX) else 0

    xlh = lh_vap(surface_temp)
    if surface_temp <= FREEZE:
        xlh += lh_fus(surface_temp)

    # latent heat flux (- away from surf)
    le = xlh * e

    return h, le, e, ier, ustar, factor


def efcon(k, layer_temp, p_a, es_layer=None):
    """
    calculates the effective thermal conductivity for a layer
    accounting for both conduction and vapor diffusion.
    Saturation within the layer is assumed.

    Args:
        k: layer thermal conductivity (J/(m K sec))
        layer_temp: layer temperature (K)
        p_a: air pressure (Pa)
        es_layer: saturation vapor pressure over ice of layer, Optional
            will be calculated if not provided

    Returns:
        etc: effective thermal conductivity (J/(m K sec))
    """

    # calculate effective layer diffusion (see Anderson, 1976, pg. 32)
    de = diffusion_coef(p_a, layer_temp)

    # set latent heat from layer temp.
    if layer_temp > FREEZE:
        lh = lh_vap(layer_temp)
    elif layer_temp == FREEZE:
        lh = (lh_vap(layer_temp) + lh_sub(layer_temp)) / 2.0
    else:
        lh = lh_sub(layer_temp)

    # set mixing ratio from layer temp.
    if es_layer is None:
        es_layer = sati(layer_temp)
    q = mix_ratio(es_layer, p_a)

    # calculate effective layer conductivity
    return k + (lh * de * q)


def ssxfr(k1, k2, t1, t2, d1, d2):
    """
    calculates the steady state heat transfer between two layers.

    Args:
        k1: layer 1's thermal conductivity (J / (m K sec))
        k2: layer 2's    "         "
        t1: layer 1's average layer temperature (K)
        t2: layer 2's    "      "        "
        d1: layer 1's thickness (m)
        d2: layer 2's    "       "

    Returns:
        g - heat transfer between layers (W/m^2)

    """

    return 2.0 * (k1 * k2 * (t2 - t1)) / ((k2 * d1) + (k1 * d2))
