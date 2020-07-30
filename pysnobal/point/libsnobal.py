import math

import numpy as np

from pysnobal.core.constants import (BOIL, FREEZE, GRAVITY,
                                     LOG_SEA_LEVEL, MOL_AIR, MOL_H2O,
                                     RGAS, SEA_LEVEL, VON_KARMAN)
from pysnobal.core.functions import (gas_density, lh_fus, lh_sub, lh_vap,
                                     mix_ratio, virtual_temperature)

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
BETA_U = 16

ITMAX = 50
LOG_10 = math.log(10.0)  # log(10)


def hysat(pb, tb, L, h, g, m):
    '''
    integral of hydrostatic equation over layer with linear
    temperature variation

        pb = base level pressure
        tb = base level temp (K)
        L  = lapse rate (deg/km)
        h  = layer thickness (km)
        g  = grav accel (m/s^2)
        m  = molec wt (kg/kmole)

     (the factors 1.e-3 and 1.e3 are for units conversion)
     20151027 Scott Havens
     '''

    if L == 0:
        return pb * np.exp(-g * m * h * 1.e3/(RGAS * tb))
    else:
        return pb * np.power(tb/(tb + L * h), g * m/(RGAS * L * 1.e-3))


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


def spec_hum(e, P):
    """
    specific humidity from vapor pressure

    e = vapor pressure
    P = pressure (same units as e)
    """

    return e * MOL_H2O / (MOL_AIR * P + e * (MOL_H2O - MOL_AIR))


# @profile
def psi(zeta, code):
    """
    psi-functions
    code =   SM    momentum
             SH    sensible heat flux
             SV    latent heat flux
    """

    if zeta > 0:        # stable
        if zeta > 1:
            zeta = 1
        result = -BETA_S * zeta

    elif zeta < 0:  # unstable

        x = math.sqrt(math.sqrt(1 - BETA_U * zeta))

        if code == 'SM':
            result = 2 * math.log((1 + x)/2) + math.log((1 + x * x)/2) - \
                2 * math.atan(x) + np.pi/2

        elif (code == 'SH') or (code == 'SV'):
            result = 2 * math.log((1 + x * x)/2)

        else:  # shouldn't reach
            raise ValueError("psi-function code not of these: SM, SH, SV")

    else:  # neutral
        result = 0

    return result


# @profile
def hle1(press, ta, ts, za, ea, es, zq, u, zu, z0):
    """
    computes sensible and latent heat flux and mass flux given
    measurements of temperature and specific humidity at surface
    and one height, wind speed at one height, and roughness
    length.  The temperature, humidity, and wind speed measurements
    need not all be at the same height.

    Args
        press air pressure (Pa)
        ta air temperature (K) at height za
        ts surface temperature (K)
        za height of air temp measurement (m)
        ea vapor pressure (Pa) at height zq
        es vapor pressure (Pa) at surface
        zq height of spec hum measurement (m)
        u wind speed (m/s) at height zu
        zu height of wind speed measurement (m)
        z0 roughness length (m)

    Outputs
        h sens heat flux (+ to surf) (W/m^2)
        le latent heat flux (+ to surf) (W/m^2)
        e mass flux (+ to surf) (kg/m^2/s)
        status status of convergence
            0      successful calculation
            -1      no convergence

    20160111 Scott Havens
    """

    # Check for bad inputs
    # heights must be positive
    if (z0 <= 0) or (zq <= z0) or (zu <= z0) or (za <= z0):
        raise ValueError(
            "height not positive z0={}, zq={}, zu={}, za={}".format(
                z0, zq, zu, za))

    # temperatures are Kelvin
    if (ta <= 0) or (ts <= 0):
        raise ValueError("temps not K ta=%f\tts=%f" % (ta, ts))

    # pressures must be positive
    if (ea <= 0) or (es <= 0) or (press <= 0) \
            or (ea >= press) or (es >= press):
        raise ValueError("press < 0 ea=%f\tes=%f\tpress=%f" % (ea, es, press))

    # vapor pressures can't exceed saturation
    # if way off stop
    es_sat = sati(ts)
    ea_w = satw(ta)
    if ((es - 25.0) > es_sat) or ((ea - 25.0) > ea_w):
        raise ValueError(
            "vp > sat es={}, ea_sat={}, ea={}, ea_sat={}".format(
                es, es_sat, ea, sati(ta)))

    # else fix them up
    if es > es_sat:
        es = es_sat

    if ea > ea_w:
        ea = ea_w

    # displacement plane height, eq. 5.3 & 5.4
    d0 = 2 * PAESCHKE * z0 / 3

    # constant log expressions
    ltsm = math.log((zu - d0) / z0)
    ltsh = math.log((za - d0) / z0)
    ltsv = math.log((zq - d0) / z0)

    # convert vapor pressures to specific humidities
    qa = spec_hum(ea, press)
    qs = spec_hum(es, press)
    q_diff = qa - qs

    # convert temperature to potential temperature
    ta += DALR * za
    t_diff = ta - ts

    # air density at press, virtual temp of geometric mean
    # of air and surface
    dens = gas_density(press, MOL_AIR, virtual_temperature(
        math.sqrt(ta * ts), math.sqrt(ea * es), press))

    # starting value, assume neutral stability, so psi-functions
    # are all zero
    ustar = VON_KARMAN * u / ltsm
    factor = VON_KARMAN * ustar * dens
    e = q_diff * factor * AV / ltsv
    h = t_diff * factor * CP_AIR * AH / ltsh

    # if not neutral stability, iterate on Obukhov stability
    # length to find solution
    it = 0
    if ta != ts:

        lo = 1e500

        while True:
            last = lo

            # Eq 4.25, but no minus sign as we define
            # positive H as toward surface
            lo = ustar * ustar * ustar * dens / \
                (VON_KARMAN * GRAVITY * (h/(ta * CP_AIR) + 0.61 * e))

            # friction velocity, eq. 4.34'
            ustar = VON_KARMAN * u / (ltsm - psi(zu/lo, 'SM'))

            # evaporative flux, eq. 4.33'
            factor = VON_KARMAN * ustar * dens
            e = q_diff * factor * AV / (ltsv - psi(zq/lo, 'SV'))

            # sensible heat flus, eq. 4.35'
            # with sign reversed
            h = t_diff * factor * AH * CP_AIR / (ltsh - psi(za/lo, 'SH'))

            diff = last - lo

            it += 1
            if (np.abs(diff) < THRESH) and (np.abs(diff/lo) < THRESH):
                break
            if it > ITMAX:
                break

    ier = -1 if (it >= ITMAX) else 0
#     print 'iterations: %i' % it
    xlh = lh_vap(ts)
    if ts <= FREEZE:
        xlh += lh_fus(ts)

    # latent heat flux (- away from surf)
    le = xlh * e

    return h, le, e, ier


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
#     de = diffusion_coef(p, t)
    de = 0.65 * (SEA_LEVEL / p_a) * \
        math.pow(layer_temp/FREEZE, 14.0) * (0.01 * 0.01)

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
