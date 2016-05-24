"""
The Snobal 'library' which is a collection of functions to run the model

20160108 Scott Havens
"""

import numpy as np
cimport numpy as np
# from math import log, pow
from libc.math cimport *
import cython

# We now need to fix a datatype for our arrays. I've used the variable
# DTYPE for this, which is assigned to the usual NumPy runtime
# type info object.
DTYPE = np.float64

# "ctypedef" assigns a corresponding compile-time type to DTYPE_t. For
# every type in the numpy module there's a corresponding compile-time
# type with a _t-suffix.
ctypedef np.float64_t DTYPE_t

cdef double FREEZE = 273.16         # freezing temp K
cdef double BOIL = 373.15           # boiling temperature K
STD_LAPSE_M = -0.0065   # lapse rate (K/m)
STD_LAPSE = -6.5        # lapse rate (K/km)
STD_AIRTMP = 2.88e2     # standard sea level air temp (K)
SEA_LEVEL = 1.013246e5  # sea level pressure
cdef double RGAS = 8.31432e3        # gas constant (J / kmole / deg)
GRAVITY = 9.80665       # gravity (m/s^2)
MOL_AIR = 28.9644       # molecular weight of air (kg / kmole)
MOL_H2O = 18.0153       # molecular weight of water vapor (kg / kmole)
cdef double c_MOL_AIR = 28.9644       # molecular weight of air (kg / kmole)
cdef double c_MOL_H2O = 18.0153       # molecular weight of water vapor (kg / kmole)
cdef double VON_KARMAN = 0.41       # Von Karman constant

cdef double CP_AIR = 1.005e3        # specific heat of air at constant pressure (J / kg / deg)
cdef double DALR = GRAVITY / CP_AIR # dry adiabatic lapse rate (deg / m)

cdef double AH = 1.0                # ratio sensible/momentum phi func    
cdef double AV = 1.0                # ratio latent/momentum phi func    
cdef int ITMAX = 50              # max # iterations allowed        
cdef double PAESCHKE = 7.35         # Paeschke's const (eq. 5.3)        
cdef double THRESH = 1.e-5          # convergence threshold        

cdef int SM = 0
cdef int SH = 1
cdef int SV = 2
cdef double BETA_S = 5.2
cdef double BETA_U = 16


# equation of state, to give density of a gas (kg/m^3)
 
cpdef GAS_DEN (p, m, t): 
    return p * m/(RGAS * t)

cdef double GAS_DEN_c (double p, double m, double t): 
    return p * m/(RGAS * t)

# virtual temperature, i.e. the fictitious temperature that air must
# have at the given pressure to have the same density as a water vapor
# and air mixture at the same pressure with the given temperature and
# vapor pressure.
@cython.cdivision(True)
cdef double VIR_TEMP (double t, double e, double P): 
    return t/(1. - (1. - MOL_H2O/MOL_AIR) * (e/P))

# latent heat of vaporization, t = temperature (K)
cdef double LH_VAP (double t): 
    return 2.5e6 - 2.95573e3 * (t - FREEZE)

# latent heat of fusion, t = temperature (K)
cpdef double LH_FUS (double t): 
    return 3.336e5 + 1.6667e2 * (FREEZE - t)

# latent heat of sublimination (J/kg), t = temperature (K)
cdef double LH_SUB (double t):
    return LH_VAP(t) + LH_FUS(t)

# mixing ratio
cdef MIX_RATIO (np.ndarray[DTYPE_t, ndim=2] e, np.ndarray[DTYPE_t, ndim=2] P): 
    return (c_MOL_H2O/c_MOL_AIR) * e/(P - e)

# effectuve diffusion coefficient (m^2/sec) for saturated porous layer
# (like snow...).  See Anderson, 1976, pg. 32, eq. 3.13.
#    pa = air pressure (Pa)
#    ts = layer temperature (K)
cpdef DIFFUS(pa, ts):
    return 0.65 * (SEA_LEVEL / pa) * pow(ts/FREEZE, 14.0) * (0.01*0.01)

# DIFFUS = lambda pa, ts:  0.65 * (SEA_LEVEL / pa) * np.power(ts/FREEZE,14.0) * (0.01*0.01)

# water vapor flux (kg/(m^2 sec)) between two layers
#   air_d = air density (kg/m^3)
#   k     = diffusion coef. (m^2/sec)
#   q_dif = specific hum. diff between layers (kg/kg)
#   z_dif = absolute distance between layers (m)
#
#   note:   q_dif controls the sign of the computed flux
# EVAP = lambda double air_d, double k, double q_dif, double z_dif: air_d * k * (q_dif/z_dif)
EVAP = lambda air_d,k,q_dif,z_dif: air_d * k * (q_dif/z_dif)


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
       

def satw_np(np.ndarray tk):
    '''
    Saturation vapor pressure of water. from IPW satw but for a numpy array
    20151027 Scott Havens
    '''
    
    # remove bad values
    tk[tk < 0] = np.nan

    l10 = np.log(10.0)

    btk = BOIL/tk
    x = -7.90298*(btk- 1.0) + 5.02808*np.log(btk)/l10 - \
            1.3816e-7*(np.power(10.0,1.1344e1*(1.0 - tk/BOIL))-1.) + \
            8.1328e-3*(np.power(10.0,-3.49149*(btk - 1.0)) - 1.0) + \
            np.log(SEA_LEVEL)/l10

    x = np.power(10.0,x)

    return x

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def sati_np(np.ndarray tk):
    '''
    saturation vapor pressure over ice. From IPW sati but for a numpy array
    20151027 Scott Havens
    '''
    
    cdef int tk_size = tk.shape[0] * tk.shape[1]
    cdef double l10
    cdef np.ndarray ind = np.zeros([tk.shape[0], tk.shape[1]], dtype=np.uint8)
    
    # remove bad values
    tk[tk < 0] = np.nan
    
    if np.sum(np.isnan(tk)) == tk_size:
        raise ValueError('sati_np: All values of tk < 0')
    
    # preallocate
    cdef np.ndarray x = np.zeros([tk.shape[0], tk.shape[1]], dtype=DTYPE)
#     x = np.empty_like(tk.shape)
    
    # vapor below freezing
    l10 = np.log(10.0)
#     x[ind] = 100.0 * np.power(10.0,-9.09718*((FREEZE/tk[ind]) - 1.0) - 3.56654*np.log(FREEZE/tk[ind])/l10 + \
#             8.76793e-1*(1.0 - (tk[ind]/FREEZE)) + np.log(6.1071)/l10)
    x = 100.0 * np.power(10.0,-9.09718*((FREEZE/tk) - 1.0) - 3.56654*np.log(FREEZE/tk)/l10 + \
            8.76793e-1*(1.0 - (tk/FREEZE)) + np.log(6.1071)/l10)

    # vapor above freezing
    ind = tk > FREEZE
    if np.any(ind):
        x[ind] = satw_np(tk[ind])

    return x


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def sati_2d(np.ndarray[DTYPE_t, ndim=2] tk):
    '''
    saturation vapor pressure over ice. From IPW sati but for a 2D numpy array
    20160523 Scott Havens
    '''
    
    cdef int tk_size = tk.shape[0] * tk.shape[1]
    cdef double l10
    cdef np.ndarray[np.uint8_t, ndim=2, cast=True] ind = np.zeros([tk.shape[0], tk.shape[1]], dtype=np.uint8)
    cdef np.ndarray[DTYPE_t, ndim=2] x = np.zeros([tk.shape[0], tk.shape[1]], dtype=DTYPE)
    
    # remove bad values
    tk[tk < 0] = np.nan
    
    if np.sum(np.isnan(tk)) == tk_size:
        raise ValueError('sati_np: All values of tk < 0')
    
    # vapor below freezing
    l10 = np.log(10.0)
    x = 100.0 * np.power(10.0,-9.09718*((FREEZE/tk) - 1.0) - 3.56654*np.log(FREEZE/tk)/l10 + \
            8.76793e-1*(1.0 - (tk/FREEZE)) + np.log(6.1071)/l10)

    # vapor above freezing
    ind = tk > FREEZE
    if np.any(ind):
        x[ind] = satw_np(tk[ind])

    return x




@cython.cdivision(True) 
cdef double satw(double tk):
    '''
    Saturation vapor pressure of water. from IPW satw but for a single value
    20160112 Scott Havens
    '''
    
    cdef double btk, x
    
    if tk < 0:
        raise ValueError('tk < 0')
    
#     l10 = log(10.0)

    btk = BOIL/tk
    x = -7.90298*(btk- 1.0) + 5.02808*log(btk)/log(10.0) - \
            1.3816e-7*(pow(10.0,1.1344e1*(1.0 - tk/BOIL))-1.) + \
            8.1328e-3*(pow(10.0,-3.49149*(btk - 1.0)) - 1.0) + \
            np.log(SEA_LEVEL)/log(10.0)

    x = pow(10.0,x)

    return x


# cpdef double sati(double tk):

@cython.cdivision(True)
cdef double sati(double tk):
    '''
    saturation vapor pressure over ice. From IPW sati but for a single value
    20151027 Scott Havens
    '''
    
    cdef double x
#     cdef double l10
#     cdef double f = FREEZE
    
    if tk < 0:
        raise ValueError('tk < 0')
        
    # vapor above freezing
    if tk > FREEZE:
        x = satw(tk)
        
    else:
        # vapor below freezing
#         l10 = log(10.0)
        x = 100.0 * pow(10.0, -9.09718*((FREEZE/tk) - 1.0) - 3.56654*log(FREEZE/tk)/log(10.0) + \
                8.76793e-1*(1.0 - (tk/FREEZE)) + log(6.1071)/log(10.0))
        

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
    rh = ea / sati_np(ta)
    rh[rh > 1] = 1
    
    e_prime = (rh * sati_np(t_prime))/100.0

    air_emiss = (1.24*np.power(e_prime/t_prime, 1./7.0))*pa/SEA_LEVEL

    air_emiss[air_emiss > 1.0] = 1.0

    return air_emiss


cpdef double spec_hum(double e, double P):
    """
    specific humidity from vapor pressure
 
    e = vapor pressure
    P = pressure (same units as e)
    """
    
    return e * c_MOL_H2O / (c_MOL_AIR * P + e * (c_MOL_H2O - c_MOL_AIR))


def spec_hum_np(np.ndarray[DTYPE_t, ndim=2] e, np.ndarray[DTYPE_t, ndim=2] P):
    """
    specific humidity from vapor pressure
 
    e = vapor pressure
    P = pressure (same units as e)
    """
    
    return e * MOL_H2O / (MOL_AIR * P + e * (MOL_H2O - MOL_AIR))


cdef double psi(double zeta, int code):
    """
    psi-functions
    code =   SM    momentum
             SH    sensible heat flux
             SV    latent heat flux
    """
    
    cdef double result, x
    
    if zeta > 0:        # stable
        if zeta > 1:
            zeta = 1
        result = -BETA_S * zeta
    

    elif zeta < 0:    #unstable

        x = sqrt(sqrt(1 - BETA_U * zeta));

        if code == SM:
            result = 2 * log((1 + x)/2) + log((1 + x * x)/2) - \
                2 * atan(x) + M_PI_2
        
        elif (code == SH) or (code == SV):
            result = 2 * log((1 + x * x)/2)
        
        else: # shouldn't reach 
            raise ValueError("psi-function code not of these: SM, SH, SV")
    
    else:   #neutral
        result = 0
    
    return result


# @profile
def hle1_grid(np.ndarray[DTYPE_t, ndim=2] press, np.ndarray[DTYPE_t, ndim=2] ta, \
              np.ndarray[DTYPE_t, ndim=2] ts, np.ndarray[DTYPE_t, ndim=2] za, \
              np.ndarray[DTYPE_t, ndim=2] ea, np.ndarray[DTYPE_t, ndim=2] es, \
              np.ndarray[DTYPE_t, ndim=2] zq, np.ndarray[DTYPE_t, ndim=2] u, \
              np.ndarray[DTYPE_t, ndim=2] zu, np.ndarray[DTYPE_t, ndim=2] z0, mask=None):
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
        mask only calculate where this is true        

    Outputs
        h sens heat flux (+ to surf) (W/m^2)    
        le latent heat flux (+ to surf) (W/m^2)    
        e mass flux (+ to surf) (kg/m^2/s)
        status status of convergence
            0      successful calculation
            -1      no convergence
        
    20160111 Scott Havens
    """
    
    cdef int ny = press.shape[0]
    cdef int nx = press.shape[1]
    cdef int ier
    
    cdef np.ndarray H = np.zeros([ny, nx], dtype=DTYPE)
    cdef np.ndarray LE = np.zeros([ny, nx], dtype=DTYPE)
    cdef np.ndarray E = np.zeros([ny, nx], dtype=DTYPE)
    cdef np.ndarray M = np.zeros([ny, nx], dtype=bool)
    
#     H = np.zeros_like(press)
#     LE = np.zeros_like(press)
#     E = np.zeros_like(press)
     
    # apply the mask
    if mask is None:
        M = np.zeros([ny, nx], dtype=bool)
    else:
        M = mask
                 
    # iterate over the array and apply hle1()
    for index,val in np.ndenumerate(ta):
        if M[index]:
            h, le, e, ier = hle1 (press[index], ta[index], ts[index], 
                                  za[index], ea[index], es[index], 
                                  zq[index], u[index], zu[index], 
                                  z0[index])
             
            if ier == -1:
                break
            H[index] = h
            LE[index] = le
            E[index] = e
             
    return H, LE, E, ier


# @profile
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def hle1 (double press, double ta, double ts, double za, double ea, \
          double es, double zq, double u, double zu, double z0):
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
    
    # define some constants to keep constant with hle1.c
    k = VON_KARMAN
    av = AV
    ah = AH
    cp = CP_AIR
    g = GRAVITY
    
    cdef double d0, ltsm, ltsh, ltsv, qa, qs, dens
    cdef double ustar, factor
    cdef diff, last, lo
#     cdef double d0, ltsm, ltsh, ltsv, qa, qs, dens
#     cdef double ustar, factor, lo, last, diff
    cdef double e, h
    cdef int it, ier
    
    # Check for bad inputs
    
    # heights must be positive 
    if (z0 <= 0) or (zq <= z0) or (zu <= z0) or (za <= z0):
        raise ValueError("height not positive z0=%f\tzq=%f\tzu=%\tza=%f" % \
                (z0, zq, zu, za))

    # temperatures are Kelvin 
    if (ta <= 0) or (ts <= 0):
        raise ValueError("temps not K ta=%f\tts=%f" % (ta, ts))

    # pressures must be positive 
    if (ea <= 0) or (es <= 0) or (press <= 0) or (ea >= press) or (es >= press):
        raise ValueError("press < 0 ea=%f\tes=%f\tpress=%f" % (ea, es, press))

    # vapor pressures can't exceed saturation 
    # if way off stop 
    if ((es - 25.0) > sati(ts)) or ((ea - 25.0) > satw(ta)):
        raise ValueError("vp > sat es=%f\tessat=%f\tea=%f\teasat=%f" % \
                (es, sati(ts), ea, sati(ta)))
    
    # else fix them up 
    if es > sati(ts):
        es = sati(ts)
        
    if ea > satw(ta):
        ea = satw(ta)
    
    #displacement plane height, eq. 5.3 & 5.4
    d0 = 2 * PAESCHKE * z0 / 3

    # constant log expressions
    ltsm = log((zu - d0) / z0)
    ltsh = log((za - d0) / z0)
    ltsv = log((zq - d0) / z0)
    
    # convert vapor pressures to specific humidities
    qa = spec_hum(ea, press)
    qs = spec_hum(es, press)
    
    # convert temperature to potential temperature
    ta += DALR * za
    
    # air density at press, virtual temp of geometric mean
    # of air and surface
    dens = GAS_DEN_c(press, MOL_AIR, VIR_TEMP(sqrt(ta * ts), sqrt(ea * es), press))
    
    # starting value, assume neutral stability, so psi-functions
    # are all zero
    ustar = k * u / ltsm
    factor = k * ustar * dens
    e = (qa - qs) * factor * av / ltsv
    h = (ta - ts) * factor * cp * ah / ltsh
    
    
    # if not neutral stability, iterate on Obukhov stability
    # length to find solution
    it = 0
    if ta != ts:

        lo = 1e500

        while True:
            last = lo

            # Eq 4.25, but no minus sign as we define
            # positive H as toward surface
            lo = ustar * ustar * ustar * dens / (k * g * (h/(ta * cp) + 0.61 * e))
                                
            # friction velocity, eq. 4.34'
            ustar = k * u / (ltsm - psi(zu/lo, SM))

            # evaporative flux, eq. 4.33'
            factor = k * ustar * dens
            e = (qa - qs) * factor * av / (ltsv - psi(zq/lo, SV))

            
            # sensible heat flus, eq. 4.35'
            # with sign reversed
            h = (ta - ts) * factor * ah * cp / (ltsh - psi(za/lo, SH))

            diff = last - lo

            it+=1
            if (fabs(diff) < THRESH) and (fabs(diff/lo) < THRESH):
                break
            if it > ITMAX:
                break

    ier = -1 if (it >= ITMAX) else 0
#     print 'iterations: %i' % it
    xlh = LH_VAP(ts)
    if ts <= FREEZE:
        xlh += LH_FUS(ts)
    
    # latent heat flux (- away from surf)
    le = xlh * e

    
    return h, le, e, ier
 
    
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef np.ndarray efcon(float k, np.ndarray[DTYPE_t, ndim=2] t, \
                       np.ndarray[DTYPE_t, ndim=2] p, \
                       np.ndarray[DTYPE_t, ndim=2] ea):
    """
    calculates the effective thermal conductivity for a layer
    accounting for both conduction and vapor diffusion.
    Saturation within the layer is assumed.
    
    Args:
        k: layer thermal conductivity (J/(m K sec)) 
        t: layer temperature (K)                    
        p: air pressure (Pa)    
        
    Returns:
        etc: effective thermal conductivity (J/(m K sec))
    """
    
    # calculate effective layer diffusion (see Anderson, 1976, pg. 32)
    #     de = DIFFUS(p, t)
#     cdef double sl = SEA_LEVEL
    cdef int ny = t.shape[0]
    cdef int nx = t.shape[1]
    cdef int i, j
    cdef double val
    
    cdef np.ndarray[DTYPE_t, ndim=2] de = np.zeros([ny, nx], dtype=DTYPE)
#     cdef np.ndarray[np.uint8_t, ndim=2, cast=True] ind = np.zeros([ny, nx], dtype=np.uint8)
    cdef np.ndarray[DTYPE_t, ndim=2] lh = np.zeros([ny, nx], dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] e = np.zeros([ny, nx], dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] q = np.zeros([ny, nx], dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] K = np.ones([ny, nx], dtype=DTYPE)
#     cdef np.ndarray SL = np.ones([ny, nx], dtype=DTYPE)
    
    K = k * K
#     SL = SEA_LEVEL * SL
    
    de = DIFFUS(p, t)
#     de = 0.65 * (SEA_LEVEL / p) * np.power(t/FREEZE, 14.0) * (0.01 * 0.01)

    # set latent heat from layer temp.
#     lh = np.zeros_like(t)
#     if np.any(t > FREEZE):
#         ind = t > FREEZE
#         lh[ind] = LH_VAP(t[ind])
#     if np.any(t == FREEZE):
#         ind = t == FREEZE
#         lh[ind] = (LH_VAP(t[ind]) + LH_SUB(t[ind])) / 2.0
#     if np.any(t < FREEZE):
#         ind = t < FREEZE
#         lh[ind] = LH_SUB(t[ind])
        
#     for index,val in np.ndenumerate(t):
#         i = index[0]
#         j = index[1]
    for i in range(ny):
        for j in range(nx):
            val = t[i,j]
            if val > FREEZE:
                lh[i,j] = LH_VAP(val)
            elif val < FREEZE:
                lh[i,j] = LH_SUB(val)
            else:
                lh[i,j] = (LH_VAP(val) + LH_SUB(val)) / 2.0

    # set mixing ratio from layer temp.
#     e = sati_2d(t)
    q = MIX_RATIO(ea, p)

    # calculate effective layer conductivity
    return K + (lh * de * q)

    
@cython.cdivision(True)
cpdef ssxfr(k1, k2, t1, t2, d1, d2):
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
