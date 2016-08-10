"""
Collection of functions that are better off in C for speed purposes

20160810
"""

import cython
import numpy as np
cimport numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

cdef extern from "envphys.h":
    int hle1(double press, double ta, double ts, double za, double ea, double es, 
             double zq, double u, double zu, double z0, int error_check, 
             double *h, double *le, double *e);
    int hle1_grid(int ngrid, double *press, double *ta, double *ts, double *za, double *ea, 
                   double *es, double *zq, double *u, double *zu, double *z0, int error_check, 
                   double *h, double *le, double *e);
    double sati_grid(int ngrid, double *tk, double *es);
    void efcon_grid(int ngrid, double *k, double *t, double *p, double *e, double *etc);



@cython.boundscheck(False)
@cython.wraparound(False)
def hle1_c(press, ta, ts,  za, ea, es, zq, u, zu, z0):
    """
    computes sensible and latent heat flux and mass flux given
    measurements of temperature and specific humidity at surface
    and one height, wind speed at one height, and roughness
    length.  The temperature, humidity, and wind speed measurements
    need not all be at the same height.
    
    Args:
        press: air pressure (Pa)            
        ta: air temperature (K) at height za    
        ts: surface temperature (K)        
        za: height of air temp measurement (m)    
        ea: vapor pressure (Pa) at height zq    
        es: vapor pressure (Pa) at surface    
        zq: height of spec hum measurement (m)    
        u: wind speed (m/s) at height zu    
        zu: height of wind speed measurement (m)    
        z0: roughness length (m)      

    Outputs
        h: sens heat flux (+ to surf) (W/m^2)    
        le: latent heat flux (+ to surf) (W/m^2)    
        e: mass flux (+ to surf) (kg/m^2/s)
        status: status of convergence
            0      successful calculation
            -1      no convergence
        
    20160111 Scott Havens
    """
    
#     cdef int ny = press.shape[0]
#     cdef int nx = press.shape[1]
    cdef int ier
    cdef int i
    
#     cdef np.ndarray H = np.zeros_like(press, dtype=DTYPE)
#     cdef np.ndarray LE = np.zeros_like(press, dtype=DTYPE)
#     cdef np.ndarray E = np.zeros_like(press, dtype=DTYPE)
#     cdef np.ndarray M = np.zeros_like(press, dtype=bool)
    
    H = np.zeros_like(press)
    LE = np.zeros_like(press)
    E = np.zeros_like(press)
            
        
    cdef double h, le, e
    
                 
    # iterate over the array and apply hle1()
    i = 0
    for index,val in np.ndenumerate(ta):
        ier = hle1 (press[index], ta[index], ts[index], 
                              za[index], ea[index], es[index], 
                              zq[index], u[index], zu[index], 
                              z0[index], 0, &h, &le, &e)
          
        if ier == -1:
            ier = i
            break
        H[index] = h
        LE[index] = le
        E[index] = e
             
        i += 1
              
    return H, LE, E, ier


@cython.boundscheck(False)
@cython.wraparound(False)
# https://github.com/cython/cython/wiki/tutorials-NumpyPointerToC
def hle1_gridded(press, ta, ts,  za, ea, es, zq, u, zu, z0, int error_check=1):
    """
    computes sensible and latent heat flux and mass flux given
    measurements of temperature and specific humidity at surface
    and one height, wind speed at one height, and roughness
    length.  The temperature, humidity, and wind speed measurements
    need not all be at the same height.
    
    Args:
        press: air pressure (Pa)            
        ta: air temperature (K) at height za    
        ts: surface temperature (K)        
        za: height of air temp measurement (m)    
        ea: vapor pressure (Pa) at height zq    
        es: vapor pressure (Pa) at surface    
        zq: height of spec hum measurement (m)    
        u: wind speed (m/s) at height zu    
        zu: height of wind speed measurement (m)    
        z0: roughness length (m)
        error_check: 1 (True) or 0 (False) if the inputs should be error checked (slower)      

    Returns
        h: sens heat flux (+ to surf) (W/m^2)    
        le: latent heat flux (+ to surf) (W/m^2)    
        e: mass flux (+ to surf) (kg/m^2/s)
        status: status of convergence
            0      successful calculation
            -1      no convergence
        
    20160810 Scott Havens
    """
    
    cdef int ngrid, ier
    
    ngrid = press.shape[0]
    
    # convert the press array to C
    cdef np.ndarray[double, mode="c"] press_arr
    press_arr = np.ascontiguousarray(press, dtype=np.float64)
    
    # convert the ta to C
    cdef np.ndarray[double, mode="c"] ta_arr
    ta_arr = np.ascontiguousarray(ta, dtype=np.float64)
    
    # convert the ts to C
    cdef np.ndarray[double, mode="c"] ts_arr
    ts_arr = np.ascontiguousarray(ts, dtype=np.float64)
    
    # convert the za to C
    cdef np.ndarray[double, mode="c"] za_arr
    za_arr = np.ascontiguousarray(za, dtype=np.float64)        
    
    # convert the ea to C
    cdef np.ndarray[double, mode="c"] ea_arr
    ea_arr = np.ascontiguousarray(ea, dtype=np.float64)
    
    # convert the es to C
    cdef np.ndarray[double, mode="c"]es_arr
    es_arr = np.ascontiguousarray(es, dtype=np.float64)

    # convert the zq to C
    cdef np.ndarray[double, mode="c"] zq_arr
    zq_arr = np.ascontiguousarray(zq, dtype=np.float64)
            
    # convert the u to C
    cdef np.ndarray[double, mode="c"] u_arr
    u_arr = np.ascontiguousarray(u, dtype=np.float64)
    
    # convert the zu to C
    cdef np.ndarray[double, mode="c"] zu_arr
    zu_arr = np.ascontiguousarray(zu, dtype=np.float64)
    
    # convert the z0 to C
    cdef np.ndarray[double, mode="c"] z0_arr
    z0_arr = np.ascontiguousarray(z0, dtype=np.float64)
            
    # output variables
    h = np.zeros_like(press)
    le = np.zeros_like(press)
    e = np.zeros_like(press)
    cdef np.ndarray[double, mode="c"] H 
    H = np.ascontiguousarray(h, dtype=np.float64)
    cdef np.ndarray[double, mode="c"] LE 
    LE = np.ascontiguousarray(le, dtype=np.float64)
    cdef np.ndarray[double, mode="c"] E
    E = np.ascontiguousarray(e, dtype=np.float64)
      
    # call the C function
    ier = hle1_grid(ngrid, &press_arr[0], &ta_arr[0], &ts_arr[0], &za_arr[0], 
              &ea_arr[0], &es_arr[0], &zq_arr[0], &u_arr[0], &zu_arr[0], 
              &z0_arr[0], error_check, &H[0], &LE[0], &E[0])

    return H, LE, E, ier


@cython.boundscheck(False)
@cython.wraparound(False)
# https://github.com/cython/cython/wiki/tutorials-NumpyPointerToC
def sati_gridded(t):
    """
    Calculates the saturation vapor pressure over ice
    
    Args:
        t: temperature [K]
        
    Returns
        es: saturation vapor pressure [Pa]
    """
    
    ngrid = t.shape[0]
    
    # convert the t to C
    cdef np.ndarray[double, mode="c"] t_arr
    t_arr = np.ascontiguousarray(t, dtype=np.float64)

    e = np.zeros_like(t)
    cdef np.ndarray[double, mode="c"] es 
    es = np.ascontiguousarray(e, dtype=np.float64)
    
    # calculate the saturation vapor pressure
    sati_grid(ngrid, &t_arr[0], &es[0])
    
    return es


@cython.boundscheck(False)
@cython.wraparound(False)
# https://github.com/cython/cython/wiki/tutorials-NumpyPointerToC
def efcon_gridded(k, t, p, e):
    """
    calculates the effective thermal conductivity for a layer
    accounting for both conduction and vapor diffusion.
    Saturation within the layer is assumed.
    
    Args:
        k: layer thermal conductivity (J/(m K sec)) 
        t: layer temperature (K)                    
        p: air pressure (Pa)
        e: vapor pressure (Pa)
        
    Returns:
        etc: effective thermal conductivity (J/(m K sec))
    """
    
    ngrid = k.shape[0]
    
    # convert the k to C
    cdef np.ndarray[double, mode="c"] k_arr
    k_arr = np.ascontiguousarray(k, dtype=np.float64)

    # convert the t to C
    cdef np.ndarray[double, mode="c"] t_arr
    t_arr = np.ascontiguousarray(t, dtype=np.float64)
    
    # convert the t to C
    cdef np.ndarray[double, mode="c"] p_arr
    p_arr = np.ascontiguousarray(p, dtype=np.float64)
    
    # convert the t to C
    cdef np.ndarray[double, mode="c"] e_arr
    e_arr = np.ascontiguousarray(e, dtype=np.float64)


    # output variable
    tmp = np.zeros_like(t)
    cdef np.ndarray[double, mode="c"] etc 
    etc = np.ascontiguousarray(tmp, dtype=np.float64)
    
    # calculate the saturation vapor pressure
    efcon_grid(ngrid, &k_arr[0], &t_arr[0], &p_arr[0], &e_arr[0], &etc[0])
    
    return etc
