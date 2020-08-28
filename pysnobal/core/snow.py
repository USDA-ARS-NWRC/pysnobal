from pysnobal.core.constants import RHO_ICE, RHO_W0, FREEZE
from pysnobal.core.functions import lh_fus, cp_ice


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


def new_tsno(spm, t0, ccon):
    """
    calculates a new temperature for a snow layer from its
    adjusted cold content, and the layer's last (previous) temperature.

    The layer's specific mass (the argument <I>spm</I>) can be computed by
    multiplying the layer's thickness (m) by its density (kg/m^3).

    Args:
        spm: layer's specific mass (kg/m^2)
        t0: layer's last temperature (K)
        ccon: layer's adjusted cold content (J/m^2)

    Returns:
        tsno: snow layer's new temperature (K)
    """

    cp = cp_ice(t0)
    tdif = ccon / (spm * cp)
    tsno = tdif + FREEZE

    return tsno


def heat_stor(cp, spm, tdif):
    """
    Calculate the heat storage
    Args:
        cp: specific heat of layer (J/kg K)
        spm: layer specific mass (kg/m^2)
        tdif: temperature change (K)
    """

    return cp * spm * tdif


def cold_content(temp, mass):
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

    return heat_stor(cp_ice(temp), mass, temp - FREEZE)
