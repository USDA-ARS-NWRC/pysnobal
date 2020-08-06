import numpy as np

from pysnobal.core.constants import (CAL_TO_J, CP_W0, FREEZE, MOL_AIR, MOL_H2O,
                                     RGAS, RHO_ICE, RHO_W0, SEA_LEVEL)


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


def k_to_c(x):
    """Kelvin to Celcius

    Args:
        x (float or array): temperature in Kelvin

    Returns:
        float or array: temperature in Celcius
    """
    return x - FREEZE


def cp_ice(t):
    """specific heat of ice (J/(kg K)) (from CRC table D-159;
    most accurate from 0 to -10 C) t - temperature in K

    Args:
        t (float or array): temperature in Kelvin

    Returns:
        [float or array]: specific heat of ice
    """
    return CAL_TO_J * (0.024928 + (0.00176 * t)) / 0.001


def cp_water(t):
    """specific heat of water (J/(kg K))
    (from CRC table D-158; most accurate from 0 to +10 C)
    (incorrect at temperatures above 25 C)

    Args:
        t (float or array): temperature in Kelvin

    Returns:
        float or array: specific heat of water
    """
    return CP_W0 - 2.55 * (t - FREEZE)


def h2o_left(d, rhos, sat):
    """water retained by snow at given saturation (see SNO_SAT)

    Args:
        d (float or array): total depth of snow [m]
        rhos (float or array): density of snow [kg/m3]
        sat (float or array): snow saturation

    Returns:
        float or array: water held within the snowpack
    """
    return (sat * d * RHO_W0 * (RHO_ICE - rhos)) / RHO_ICE


def melt(Q):
    """snow melt

    Args:
        Q (float or array): available energy [J/m2]

    Returns:
        float or array: snow melt [kg/m2]
    """
    return Q / lh_fus(FREEZE)


def time_average(avg, total_time, value, time_incr):
    """A macro to update a time-weighted average for a quantity.

    Args:
        avg (float or array): current average
        total_time (float): the time interval the current average applies to
        value (float or array): new value to be averaged in
        time_incr (float): the time interval the new value applies to

    Returns:
        float or array: time averaged value
    """
    return (avg * total_time + value * time_incr) / \
        (total_time + time_incr)


def gas_density(pressure, mole_air, temp):
    """equation of state, to give density of a gas (kg/m^3)

    Args:
        pressure (float or array): pressure
        mole_air (float or array): molecular weight of air
        temp (float or array): temperature

    Returns:
        float or array: density of gas [kg/m3]
    """

    return pressure * mole_air/(RGAS * temp)


def virtual_temperature(temp, vapor_pressure, pressure):
    """virtual temperature, i.e. the fictitious temperature that air must
    have at the given pressure to have the same density as a water vapor
    and air mixture at the same pressure with the given temperature and
    vapor pressure.

    Args:
        temp (float or array): temperature [K]
        vapor_pressure (float or array): vapor pressure [Pa]
        pressure (float or array): pressure [Pa]

    Returns:
        float or array: virtual temperature [K]
    """
    return temp/(1. - (1. - MOL_H2O/MOL_AIR) * (vapor_pressure/pressure))


def lh_vap(temp):
    """latent heat of vaporization

    Args:
        temp (float or array): temperature [K]

    Returns:
        float or array: latent heat of vaporization
    """
    return 2.5e6 - 2.95573e3 * (temp - FREEZE)


def lh_fus(temp):
    """latent heat of fusion

    Args:
        temp (float or array): temperature [K]

    Returns:
        float or array: latent heat of fusion
    """
    return 3.336e5 + 1.6667e2 * (FREEZE - temp)


def lh_sub(temp):
    """latent heat of sublimination (J/kg)

    Args:
        temp (float or array): temperature [K]

    Returns:
        float or array: latent heat of sublimination (J/kg)
    """
    return lh_vap(temp) + lh_fus(temp)


def mix_ratio(vapor_pressure, pressure):
    """mixing ratio

    Args:
        vapor_pressure (float or array): vapor pressure [Pa]
        pressure (float or array): pressure [Pa]

    Returns:
        float or array: mixing ratio
    """
    return (MOL_H2O/MOL_AIR) * vapor_pressure/(pressure - vapor_pressure)


def diffusion_coef(pressure, temp):
    """effective diffusion coefficient(m ^ 2/sec) for saturated porous layer
    (like snow...).  See Anderson, 1976, pg. 32, eq. 3.13.

    Args:
        pressure (float or array): pressure [Pa]
        temp (float or array): temperature [K]

    Returns:
        float or array: effective diffusion coefficient
    """

    return 0.65 * (SEA_LEVEL / pressure) * \
        np.power(temp/FREEZE, 14.0) * (0.01 * 0.01)


def vapor_flux(air_density, k, q_dif, z_dif):
    """water vapor flux (kg/(m^2 sec)) between two layers

    Args:
        air_density (float or array): air density (kg/m^3)
        k (float or array): diffusion coef. (m^2/sec)
        q_dif (float or array): specific hum. diff between layers (kg/kg)
        z_dif (float or array): absolute distance between layers (m)

    Returns:
        float or array: water vapor flux
    """
    return air_density * k * (q_dif/z_dif)
