# -*- coding: utf-8 -*-
"""
ipysnobal: the Python implementation of iSnobal

This is not a replica of iSnobal but my interpretation and
porting to Python.  See pysnobal.exact for more direct
interpretation

20160118 Scott Havens
"""

# from libsnobal import libsnobal
from libsnobal.py_snobal import snobal

import ConfigParser
import sys, os
import numpy as np
import pandas as pd
from datetime import timedelta
import netCDF4 as nc
import matplotlib.pyplot as plt
import progressbar
from multiprocessing import Pool
# from functools import partial
import itertools

# os.system("taskset -p 0xff %d" % os.getpid())

DEFAULT_MAX_Z_S_0 = 0.25
DEFAULT_MAX_H2O_VOL = 0.01

DATA_TSTEP = 0
NORMAL_TSTEP = 1
MEDIUM_TSTEP = 2
SMALL_TSTEP = 3

DEFAULT_NORMAL_THRESHOLD = 60.0
DEFAULT_MEDIUM_TSTEP = 15.0
DEFAULT_SMALL_TSTEP = 1.0

WHOLE_TSTEP = 0x1 # output when tstep is not divided
DIVIDED_TSTEP = 0x2  # output when timestep is divided

hrs2min = lambda x: x * 60
min2sec = lambda x: x * 60
C_TO_K = 273.16
FREEZE = C_TO_K


# def _pickle_method(method):
#     func_name = method.im_func.__name__
#     obj = method.im_self
#     cls = method.im_class
#     if func_name.startswith('__') and not func_name.endswith('__'): #deal with mangled names
#         cls_name = cls.__name__.lstrip('_')
#         func_name = '_' + cls_name + func_name
#     return _unpickle_method, (func_name, obj, cls)
# 
# def _unpickle_method(func_name, obj, cls):
#     for cls in cls.__mro__:
#         try:
#             func = cls.__dict__[func_name]
#         except KeyError:
#             pass
#         else:
#             break
#     return func.__get__(obj, cls)
# 
# import copy_reg
# import types
# copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)

# parse configuration file
class MyParser(ConfigParser.ConfigParser):
    def as_dict(self):
        d = dict(self._sections)
        for k in d:
            d[k] = dict(self._defaults, **d[k])
            d[k].pop('__name__', None)
        d = self._make_lowercase(d)
        return d
    
    def _make_lowercase(self, obj):
        if hasattr(obj,'iteritems'):
            # dictionary
            ret = {}
            for k,v in obj.iteritems():
                ret[self._make_lowercase(k)] = v
            return ret
        elif isinstance(obj,basestring):
            # string
            return obj.lower()
        elif hasattr(obj,'__iter__'):
            # list (or the like)
            ret = []
            for item in obj:
                ret.append(self._make_lowercase(item))
            return ret
        else:
            # anything else
            return obj


def check_range(value, min_val, max_val, descrip):
    """
    Check the range of the value
    Args:
        value: value to check
        min_val: minimum value
        max_val: maximum value
        descrip: short description of input
    Returns:
        True if within range
    """
    if (value < min_val) or (value > max_val):
        raise ValueError("%s (%f) out of range: %f to %f", descrip, value, min_val, max_val);
    pass


def date_range(start_date, end_date, increment):
    '''
    Calculate a list between start and end date with
    an increment
    '''
    result = []
    nxt = start_date
    while nxt <= end_date:
        result.append(nxt)
        nxt += increment
    return np.array(result)


def get_args(configFile):
    """
    Parse the configuration file
    
    Args:
        configFile: configuration file for ipysnobal
        
    Returns:
        options: options structure with defaults if not set
        
        options = {
            z: site elevation (m),
            t: time steps: data [normal, [,medium [,small]]] (minutes),
            m: snowcover's maximum h2o content as volume ratio,
            d: maximum depth for active layer (m),
            
            s: snow properties input data file,
            h: measurement heights input data file,
            p: precipitation input data file,
            i: input data file,
            I: initial conditions 
            o: optional output data file,
            O: how often output records written (data, normal, all),
            c: continue run even when no snowcover,
            K: accept temperatures in degrees K,
            T: run timesteps' thresholds for a layer's mass (kg/m^2),
        }
        
    To-do: take all the rest of the defualt and check ranges for the
    input arguements, i.e. rewrite the rest of getargs.c
    """

    # read the config file and store
    if not os.path.isfile(configFile):
        raise Exception('Configuration file does not exist --> %s' % configFile)
        

    f = MyParser()
    f.read(configFile)
    config = f.as_dict()
    
    #------------------------------------------------------------------------------
    # read in the constants
    c = {}
    for v in config['constants']:
        c[v] = float(config['constants'][v])

    # these are the default options
    options = {
        'time_step': 60,
        'max-h2o': 0.01,
#         'max_z0': DEFAULT_MAX_Z_S_0,
        'c': True,
        'K': True,
        'mass_threshold': DEFAULT_NORMAL_THRESHOLD,
        'time_z': 0,
        'max_z_s_0': DEFAULT_MAX_Z_S_0,
        'z_u': 5.0,
        'z_t': 5.0,
        'z_g': 0.5,
        'relative_heights': True
    }
    options.update(c) # update the defult with any user values
    
    config['constants'] = options
        
    #------------------------------------------------------------------------------ 
    # read in the time and ensure a few things
    # nsteps will only be used if end_date is not specified
    data_tstep_min = int(config['time']['time_step'])
    check_range (data_tstep_min, 1.0, hrs2min(60),"input data's timestep")
    if ((data_tstep_min > 60) and (data_tstep_min % 60 != 0)):
        raise ValueError("Data timestep > 60 min must be multiple of 60 min (whole hrs)")
    config['time']['time_step'] = data_tstep_min
    
    # add to constant sections for tstep_info calculation
    config['constants']['time_step'] = config['time']['time_step']  
    
    
    # read in the start date and end date
    start_date = pd.to_datetime(config['time']['start_date'])
    
    if 'end_date' in config['time']:
        end_date = pd.to_datetime(config['time']['end_date'])
        if end_date < start_date:
            raise ValueError('end_date is before start_date')
        nsteps = (end_date-start_date).total_seconds()/60  # elapsed time in minutes
        nsteps = int(nsteps / config['time']['time_step'])
        
    elif 'nsteps' in config['time']:
        nsteps = int(config['time']['nsteps'])
        
        end_date = start_date + timedelta(minutes=nsteps*config['time']['time_step'])
        
    else:
        raise Exception('end_date or nsteps must be specified')
    
    # create a date time vector
    dv = date_range(start_date, end_date, 
                    timedelta(minutes=config['constants']['time_step']))
    
    if len(dv) != nsteps + 1:
        raise Exception('nsteps does not work with selected start and end dates')
    
    config['time']['start_date'] = start_date
    config['time']['end_date'] = end_date
    config['time']['nsteps'] = nsteps
    config['time']['date_time'] = dv
    
    # check the output section
    config['output']['frequency'] = int(config['output']['frequency'])
    
    try:
        config['output']['nthreads'] = int(config['output']['nthreads'])
    except:
        config['output']['nthreads'] = None
        
    return config


def get_tstep_info(options):
    """
    Parse the options dict, set the default values if not specified
    May need to divide tstep_info and params up into different
    functions
    """

    # intialize the time step info
    # 0 : data timestep
    # 1 : normal run timestep
    # 2 : medium  "     "
    # 3 : small   "     "
    
    tstep_info = []
    for i in range(4):
        t = {}
        t['level'] = i;
        t['output'] = False;
        tstep_info.append(t)
    

    # The input data's time step must be between 1 minute and 6 hours.
    # If it is greater than 1 hour, it must be a multiple of 1 hour, e.g.
    # 2 hours, 3 hours, etc.
     
    data_tstep_min = float(options['time_step'])
    tstep_info[DATA_TSTEP]['time_step'] = min2sec(data_tstep_min)
    
    norm_tstep_min = 60.0
    tstep_info[NORMAL_TSTEP]['time_step'] = min2sec(norm_tstep_min)
    tstep_info[NORMAL_TSTEP]['intervals'] = int(data_tstep_min / norm_tstep_min)
    
    med_tstep_min = DEFAULT_MEDIUM_TSTEP
    tstep_info[MEDIUM_TSTEP]['time_step'] = min2sec(med_tstep_min)
    tstep_info[MEDIUM_TSTEP]['intervals'] = int(norm_tstep_min / med_tstep_min)
    
    small_tstep_min = DEFAULT_SMALL_TSTEP
    tstep_info[SMALL_TSTEP]['time_step'] = min2sec(small_tstep_min)
    tstep_info[SMALL_TSTEP]['intervals'] = int(med_tstep_min / small_tstep_min)
    
    # output
#     if options['O'] == 'data':
#         tstep_info[DATA_TSTEP]['output'] = DIVIDED_TSTEP
#     elif options['O'] == 'normal':
#         tstep_info[NORMAL_TSTEP]['output'] = WHOLE_TSTEP | DIVIDED_TSTEP
#     elif options['O'] == 'all':
#         tstep_info[NORMAL_TSTEP]['output'] = WHOLE_TSTEP
#         tstep_info[MEDIUM_TSTEP]['output'] = WHOLE_TSTEP
#         tstep_info[SMALL_TSTEP]['output'] = WHOLE_TSTEP
#     else:
#         tstep_info[DATA_TSTEP]['output'] = DIVIDED_TSTEP
    tstep_info[DATA_TSTEP]['output'] = DIVIDED_TSTEP
    
    # mas thresholds for run timesteps
    threshold = DEFAULT_NORMAL_THRESHOLD
    tstep_info[NORMAL_TSTEP]['threshold'] = threshold
    
    threshold = DEFAULT_MEDIUM_TSTEP
    tstep_info[MEDIUM_TSTEP]['threshold'] = threshold
    
    threshold = DEFAULT_SMALL_TSTEP
    tstep_info[SMALL_TSTEP]['threshold'] = threshold
    
    
    # get the rest of the parameters
    params = {}
    
#     params['elevation'] = options['z']
    params['data_tstep'] = data_tstep_min
    params['max_h2o_vol'] = options['max-h2o']
    params['max_z_s_0'] = options['max_z_s_0']
#     params['sn_filename'] = options['s']
#     params['mh_filename'] = options['h']
#     params['in_filename'] = options['i']
#     params['pr_filename'] = options['p']
    params['out_filename'] = None
#     params['out_file'] = open(params['out_filename'], 'w')
    params['stop_no_snow'] = options['c']
    params['temps_in_C'] = options['K']
    params['relative_heights'] = options['relative_heights']

    return params, tstep_info


def open_files(options):
    """
    Open the netCDF files for initial conditions and inputs
    - Reads in the initial_conditions file
        Required variables are x,y,z,z_0
        The others z_s, rho, T_s_0, T_s, h2o_sat, mask can be specified
        but will be set to default of 0's or 1's for mask
    
    - Open the files for the inputs and store the file identifier
        
    """
    
    #------------------------------------------------------------------------------ 
    # get the initial conditions
    i = nc.Dataset(options['initial_conditions']['file'])
    
    # read the required variables in
    init = {}
    init['x'] = i.variables['x'][:]         # get the x coordinates
    init['y'] = i.variables['y'][:]         # get the y coordinates
    init['z'] = i.variables['z'][:]         # get the elevation
    init['z_0'] = i.variables['z_0'][:]     # get the roughness length
    
    # All other variables will be assumed zero if not present
    all_zeros = np.zeros_like(init['z'])
    flds = ['z_s', 'rho', 'T_s_0', 'T_s', 'h2o_sat', 'mask']
    
    for f in flds:  
        if i.variables.has_key(f):
            init[f] = i.variables[f][:]         # read in the variables
        elif f == 'mask':
            init[f] = np.ones_like(init['z'])   # if no mask set all to ones so all will be ran
        else:            
            init[f] = all_zeros                 # default is set to zeros
    
    i.close()
    
    for key in init.keys():
        init[key] = init[key].astype(np.float64)
    
    # convert temperatures to K
    init['T_s'] += FREEZE
    init['T_s_0'] += FREEZE
    
    
    #------------------------------------------------------------------------------ 
    # get the forcing data and open the file
    force = {}
    force['thermal'] = nc.Dataset(options['inputs']['thermal'], 'r')
    force['air_temp'] = nc.Dataset(options['inputs']['air_temp'], 'r')
    force['vapor_pressure'] = nc.Dataset(options['inputs']['vapor_pressure'], 'r')
    force['wind_speed'] = nc.Dataset(options['inputs']['wind_speed'], 'r')
    force['net_solar'] = nc.Dataset(options['inputs']['net_solar'], 'r')
    
    # soil temp can either be distributed for set to a constant
    try:
        force['soil_temp'] = nc.Dataset(options['inputs']['soil_temp'], 'r')
    except:
        force['soil_temp'] = float(options['inputs']['soil_temp']) * np.ones_like(init['z'])
        
    force['precip_mass'] = nc.Dataset(options['inputs']['precip_mass'], 'r')
    force['percent_snow'] = nc.Dataset(options['inputs']['percent_snow'], 'r')
    force['snow_density'] = nc.Dataset(options['inputs']['snow_density'], 'r')
    force['precip_temp'] = nc.Dataset(options['inputs']['precip_temp'], 'r')
    
    return init, force
    

def close_files(force):
    
    for f in force.keys():
        force[f].close()    
        

def output_files(options, init):
    """
    Create the snow and em output netCDF file
    """
    
    #------------------------------------------------------------------------------
    # EM netCDF 
    m = {}
    m['name'] = ['net_rad','sensible_heat','latent_heat','snow_soil','precip_advected','sum_EB','evaporation','snowmelt','SWI','cold_content']
    m['units'] = ['W m-2','W m-2','W m-2','W m-2','W m-2','W m-2','kg m-2','kg m-2','kg or mm m-2','J m-2']
    m['description'] =['Average net all-wave radiation','Average sensible heat transfer','Average latent heat exchange','Average snow/soil heat exchange',
                     'Average advected heat from precipitation','Average sum of EB terms for snowcover','Total evaporation',
                     'Total snowmelt','Total runoff','Snowcover cold content']
    
    netcdfFile = os.path.join(options['output']['location'], 'em.nc')
    dimensions = ('time','y','x')
    
    em = nc.Dataset(netcdfFile, 'w')
    
    # create the dimensions
    em.createDimension('time',None)
    em.createDimension('y',len(init['y']))
    em.createDimension('x',len(init['x']))
    
    # create some variables
    em.createVariable('time', 'f', dimensions[0])
    em.createVariable('y', 'f', dimensions[1])
    em.createVariable('x', 'f', dimensions[2])
    
    setattr(em.variables['time'], 'units', 'hours since %s' % options['time']['start_date'])
    setattr(em.variables['time'], 'calendar', 'standard')
#     setattr(em.variables['time'], 'time_zone', time_zone)
    em.variables['x'][:] = init['x']
    em.variables['y'][:] = init['y']
    
    # em image
    for i,v in enumerate(m['name']):
        
#         em.createVariable(v, 'f', dimensions[:3], chunksizes=(6,10,10))
        em.createVariable(v, 'f', dimensions[:3])
        setattr(em.variables[v], 'units', m['units'][i])
        setattr(em.variables[v], 'description', m['description'][i])
        
    options['output']['em'] = em
    
    #------------------------------------------------------------------------------
    # SNOW netCDF 
    
    s = {}
    s['name'] = ['thickness','snow_density','specific_mass','liquid_water','temp_surf','temp_lower','temp_snowcover','thickness_lower','water_saturation']
    s['units'] = ['m','kg m-3','kg m-2','kg m-2','C','C','C','m','percent']
    s['description'] =['Predicted thickness of the snowcover','Predicted average snow density','Predicted specific mass of the snowcover',
                       'Predicted mass of liquid water in the snowcover','Predicted temperature of the surface layer',
                       'Predicted temperature of the lower layer','Predicted temperature of the snowcover',
                       'Predicted thickness of the lower layer', 'Predicted percentage of liquid water saturation of the snowcover']
    
    netcdfFile = os.path.join(options['output']['location'], 'snow.nc')
    dimensions = ('time','y','x')
    
    snow = nc.Dataset(netcdfFile, 'w')
    
    # create the dimensions
    snow.createDimension('time',None)
    snow.createDimension('y',len(init['y']))
    snow.createDimension('x',len(init['x']))
    
    # create some variables
    snow.createVariable('time', 'f', dimensions[0])
    snow.createVariable('y', 'f', dimensions[1])
    snow.createVariable('x', 'f', dimensions[2])
    
    setattr(snow.variables['time'], 'units', 'hours since %s' % options['time']['start_date'])
    setattr(snow.variables['time'], 'calendar', 'standard')
#     setattr(snow.variables['time'], 'time_zone', time_zone)
    snow.variables['x'][:] = init['x']
    snow.variables['y'][:] = init['y']
    
    # snow image
    for i,v in enumerate(s['name']):
        
#         snow.createVariable(v, 'f', dimensions[:3], chunksizes=(6,10,10))
        snow.createVariable(v, 'f', dimensions[:3])
        setattr(snow.variables[v], 'units', s['units'][i])
        setattr(snow.variables[v], 'description', s['description'][i])
    
    
    options['output']['snow'] = snow
    
    
def output_timestep(s, tstep, options):
    """
    Output the model results for the current time step
    """
    
    em_out = {'net_rad':'R_n_bar', 'sensible_heat':'H_bar', 'latent_heat': 'L_v_E_bar', 
              'snow_soil': 'G_bar', 'precip_advected': 'M_bar', 
              'sum_EB': 'delta_Q_bar', 'evaporation': 'E_s_sum', 
              'snowmelt': 'melt_sum', 'SWI': 'ro_pred_sum', 
              'cold_content': 'cc_s'}
    snow_out = {'thickness': 'z_s', 'snow_density': 'rho', 
                'specific_mass': 'm_s', 'liquid_water': 'h2o', 
                'temp_surf': 'T_s_0', 'temp_lower': 'T_s_l', 
                'temp_snowcover': 'T_s', 'thickness_lower': 'z_s_l', 
                'water_saturation': 'h2o_sat'}
    
#             sbuf[sbuf_start++] = T_s_0 - FREEZE;
#             sbuf[sbuf_start++] = T_s_l - FREEZE;
#             sbuf[sbuf_start++] = T_s   - FREEZE;
    
    # preallocate
    em = {key: np.zeros(s.shape) for key in em_out.keys()}
    snow = {key: np.zeros(s.shape) for key in snow_out.keys()}
    
    # gather all the data together
    for index, si in np.ndenumerate(s):
        
        if si is not None:
            for key,value in em_out.iteritems():
                em[key][index] = getattr(si.em, value)
                
            for key,value in snow_out.iteritems():
                snow[key][index] = getattr(si.snow, value)
    
    # convert from K to C
    snow['temp_snowcover'] -= FREEZE
    snow['temp_surf'] -= FREEZE
    snow['temp_lower'] -= FREEZE
        
            
    # now find the correct index
    # the current time integer
    times = options['output']['snow'].variables['time']
    t = nc.date2num(tstep.replace(tzinfo=None), times.units, times.calendar)
            
    if len(times) != 0:
        index = np.where(times[:] == t)[0]
        if index.size == 0:
            index = len(times)
        else:
            index = index[0]
    else:
        index = len(times)
        
        
    # insert the time
    options['output']['snow'].variables['time'][index] = t
    options['output']['em'].variables['time'][index] = t
        
    # insert the data
    for key in em_out:
        options['output']['em'].variables[key][index,:] = em[key]
    for key in snow_out:
        options['output']['snow'].variables[key][index,:] = snow[key]
        
    # sync to disk
    options['output']['snow'].sync()
    options['output']['em'].sync()
    
    
def get_timestep(force, tstep):
    """
    Pull out a time step from the forcing files and 
    place that time step into a dict
    """
    
    inpt = {}
    
    # map function from these values to the ones requried by snobal
    map_val = {'air_temp': 'T_a', 'net_solar': 'S_n', 'thermal': 'I_lw',
               'vapor_pressure': 'e_a', 'wind_speed': 'u',
               'soil_temp': 'T_g', 'precip_mass': 'm_pp',
               'percent_snow': 'percent_snow', 'snow_density': 'rho_snow',
               'precip_temp': 'T_pp'}
   
    
    for f in force.keys():
        
        # determine the index
        t = nc.date2index(tstep, force[f].variables['time'], 
                          calendar=force[f].variables['time'].calendar,
                          select='exact')
        
        # pull out the value        
        inpt[map_val[f]] = force[f].variables[f][t,:].astype(np.float64)
    
    # convert from C to K
    inpt['T_a'] += FREEZE
    inpt['T_pp'] += FREEZE
    inpt['T_g'] += FREEZE
    
    return inpt
    
        

def initialize(params, tstep_info, mh, init):
    """
    Initialize pysnobal over the grid
    
    Args:
        params: parameters from get_tstep_info
        tstep_info: time step information
        mh: measurement height dict
        init: initial conditions dictionary
        
    Outputs:
        s: array of pysnobal classes
    """
    
    # variables needed for the snow properties
    # time_s will always be zero since the indicies will always
    # start at zero for ipysnobal
    v = ['time_s', 'z_s', 'rho', 'T_s', 'T_s_0', 'h2o_sat']
    sn = {key: 0.0 for key in v}
    
    # allocate an empty numpy array to hold all the snobal objects
#     s = np.empty_like(init['z'], dtype=object)
    
    
    # add to the parameters
#     params['elevation'] = init['z']
    
    # add to the measurement heights
#     mh['z_0'] = init['z_0']
    
    init['time_s'] = 0.0
    
    s = snobal(params, tstep_info, init, mh)
    
#     for index, x in np.ndenumerate(init['z']):
#         
#         if init['mask'][index]:
#         
#             # fill the initial snow record properties
#             for vi in v[1:]:
#                 sn[vi] = init[vi][index]
#             
#             # add to the parameters
#             params['elevation'] = init['z'][index]
#             
#             # add to the measurement heights
#             mh['z_0'] = init['z_0'][index]
#             
#             # initialize snobal
#             s[index] = snobal(params, tstep_info, sn, mh)
            
    return s
    
    
# class parallel_helper(object):
#     """
#     Simple class to aid in parallelizing the loop over the
#     snobal instances
#     """
#     def __init__(self, input1, input2, s):
#         self.input1 = input1
#         self.input2 = input2
#         self.s = s
#         
#     def run(self, index=None):
#         """
#         Args: index - tuple for the index to be ran, from iterator object
#         """
#         
#         if self.s[index] is not None:
#         
#             in1 = {key: self.input1[key][index] for key in self.input1.keys()}
#             in2 = {key: self.input2[key][index] for key in self.input2.keys()}
#             
#             self.s[index].do_data_tstep(in1, in2)
#             
#     def go(self):
#         
#         pool = Pool(8)
# #         it = np.nditer(self.input1, flags=['multi_index','refs_ok'])
#         it = np.ndenumerate(self.s)
#         pool.map(self.run, it)

    

    
# @profile
def run(s, input1, input2):
    """
    Acutally run the model for a single processor
    """
    
    for index, si in np.ndenumerate(s):
        
        if si is not None:
        
            in1 = {key: input1[key][index] for key in input1.keys()}
            in2 = {key: input2[key][index] for key in input2.keys()}
            
            s[index].do_data_tstep(in1, in2)


def run_map(inpt):
    """
    Run a single point, input is a list
    [snobal instance, input1, input2]
    """    
    if inpt[1] is not None:
        inpt[1].do_data_tstep(inpt[2], inpt[3])
        
    return inpt[1]
        
        

class MyIterator:
    """
    Simple iterator class that will iterate through
    the snobal classes and inputs
    """
    def __init__(self, obj, max_value):
        self.obj = obj
        self.cnt = 0
        self.max_value = max_value
       
    def __iter__(self):
        return self
    
    def next(self):
        """
        Return the next value from the object
        """
        while self.cnt < self.max_value:
            result = self.obj.get(self.cnt)
            self.cnt += 1
            return result
        
        raise StopIteration
    
    
      
class SnobalIterator:
    def __init__(self, s, input1, input2):
        self.s = s
        self.input1 = input1
        self.input2 = input2
        
        self.max_value = s.size
        self.cnt = 0

    def __iter__(self):
        return self  #MyIterator(self, self.s.size)
    
    def next(self):
        """
        Return the next value from the object
        """
        while self.cnt < self.max_value:
            result = self.get(self.cnt)
            self.cnt += 1
            return result
        
        raise StopIteration
    

    def get(self, index):
        # most likely will have to return the values for s,input1/2
        # then will have to pass to another function
       
        # get the index to subindex
        i = np.unravel_index(index, self.s.shape)

        if self.s[i] is not None:
       
            # get the input values
            in1 = {key: self.input1[key][i] for key in self.input1.keys()}
            in2 = {key: self.input2[key][i] for key in self.input2.keys()}
            
            # run the model     
#             self.s[i].do_data_tstep(in1, in2)
            
            return [index, self.s[i], in1, in2]
        
        else:
            return [index, None]
        


#@profile
def main(configFile):
    """
    mimic the main.c from the Snobal model
    
    Args:
        configFile: path to configuration file
    """
    
    # parse the input arguments
    options = get_args(configFile)
    
    # get the timestep info
    params, tstep_info = get_tstep_info(options['constants'])

    # open the files and read in data
    init, force = open_files(options)
    
    # initialize
    s = initialize(params, tstep_info, options['constants'], init)
    
    # create the output files
    output_files(options, init)
    
    # create a pool if needed
#     pool = None
#     if options['output']['nthreads'] is not None:
#         pool = Pool(processes=options['output']['nthreads'])
    
    
    # loop through the input
    # do_data_tstep needs two input records so only go 
    # to the last record-1
    
    input1 = get_timestep(force, options['time']['date_time'][0])
    
    pbar = progressbar.ProgressBar(max_value=len(options['time']['date_time'])-1)
    j = 0
    for tstep in options['time']['date_time'][1:-1]:
        
        input2 = get_timestep(force, tstep)
    
        s.do_data_tstep(input1, input2)
    
        
#         isnobal(s, input1, input2)

    
#         if pool is not None:
#             m = list(pool.imap(run_map, SnobalIterator(s, input1, input2), chunksize=100))
#             s = np.array(m)
#             s = s.reshape(init['z'].shape)
#             
#         else:
        
#         m = list(itertools.imap(run_map, SnobalIterator(s, input1, input2)))
# #             run(s, input1, input2)
        
        input1 = input2.copy()
        
        # output at the frequency and the last time step
        if (j % options['output']['frequency'] == 0) or (j == len(options['time']['date_time'])):
            output_timestep(s, tstep, options)
        
        j += 1
        pbar.update(j)
        
    pbar.finish()
          
    
    # output
#     params['out_file'].close()
    close_files(force)
#     app = MyApplication()
#     app.run()




if __name__ == "__main__":
        
    if len(sys.argv) == 1:
        raise Exception('Configuration file must be specified')
    else:
        configFile = sys.argv[1]
    
    main(configFile)

    
    

