"""
The Snobal 'library' which is a collection of functions to run the model

20160108 Scott Havens
"""

import numpy as np

FREEZE = 273.16         # freezing temp K
BOIL = 373.15           # boiling temperature K
STD_LAPSE_M = -0.0065   # lapse rate (K/m)
STD_LAPSE = -6.5        # lapse rate (K/km)
STD_AIRTMP = 2.88e2     # standard sea level air temp (K)
SEA_LEVEL = 1.013246e5  # sea level pressure
RGAS = 8.31432e3        # gas constant (J / kmole / deg)
GRAVITY = 9.80665       # gravity (m/s^2)
MOL_AIR = 28.9644       # molecular weight of air (kg / kmole)


def hysat(pb, tb, L, h, g, m):        
    '''
    integral of hydrostatic equation over layer with linear temperature variation
    
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
       

def satw(tk):
    '''
    Saturation vapor pressure of water. from IPW satw
    20151027 Scott Havens
    '''
    
    # remove bad values
    tk[tk < 0] = np.nan

    l10 = np.log(10.0)

    btk = BOIL/tk
    x = -7.90298*(btk- 1.0) + 5.02808*np.log(btk)/l10 - \
            1.3816e-7*(np.power(10.0,1.1344e1*(1.0 - tk/BOIL))-1.) + \
            8.1328e-3*(np.power(10.0,-3.49149*(btk - 1.0)) - 1.0) + \
            np.log(SEA_LEVEL)/l10;

    x = np.power(10.0,x);

    return x


def sati(tk):
    '''
    saturation vapor pressure over ice. From IPW sati
    20151027 Scott Havens
    '''
    
    # remove bad values
    tk[tk < 0] = np.nan
    
    # preallocate
    x = np.empty(tk.shape)
    
    # vapor above freezing
    ind = tk > FREEZE
    x[ind] = satw(tk[ind])
    
    # vapor below freezing
    l10 = np.log(10.0)
    x[~ind] = 100.0 * np.power(10.0, -9.09718*((FREEZE/tk[~ind]) - 1.0) - 3.56654*np.log(FREEZE/tk[~ind])/l10 + \
            8.76793e-1*(1.0 - (tk[~ind]/FREEZE)) + np.log(6.1071)/l10)


    return x


def brutsaert(ta, l, ea, z, pa):
    '''
    Calculate atmosphere emissivity
    
    ta - air temp (K)
    l - temperature lapse rate (deg/m)
    ea - vapor pressure (Pa)
    z - elevation (z)
    pa - air pressure (Pa)
    
    20151027 Scott Havens
    '''
    
    t_prime = ta - (l * z)
    rh = ea / sati(ta)
    rh[rh > 1] = 1
    
    e_prime = (rh * sati(t_prime))/100.0

    air_emiss = (1.24*np.power(e_prime/t_prime, 1./7.0))*pa/SEA_LEVEL

    air_emiss[air_emiss > 1.0] = 1.0

    return air_emiss