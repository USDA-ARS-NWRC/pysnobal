import math

DATA_TSTEP = 0
NORMAL_TSTEP = 1
MEDIUM_TSTEP = 2
SMALL_TSTEP = 3

FREEZE = 273.16         # freezing temp K
BOIL = 373.15           # boiling temperature K
STD_LAPSE_M = -0.0065   # lapse rate (K/m)
STD_LAPSE = -6.5        # lapse rate (K/km)
STD_AIRTMP = 2.88e2     # standard sea level air temp (K)
SEA_LEVEL = 1.013246e5  # sea level pressure
LOG_SEA_LEVEL = math.log(SEA_LEVEL)  # log of sea level pressure
RGAS = 8.31432e3        # gas constant (J / kmole / deg)
GRAVITY = 9.80665       # gravity (m/s^2)
MOL_AIR = 28.9644       # molecular weight of air (kg / kmole)
MOL_H2O = 18.0153       # molecular weight of water vapor (kg / kmole)
VON_KARMAN = 0.41       # Von Karman constant

HR_TO_SEC = 3600.0

MIN_SNOW_TEMP = -75
KT_MOISTSAND = 1.65
MAX_SNOW_DENSITY = 600

# density of water at 0C (kg/m^3) (from CRC handbook pg F-11)
RHO_W0 = 999.87

# specific heat of water at 0C (J / (kg K))
CP_W0 = 4217.7

# density of ice - no air (kg/m^3) (from CRC handbook pg F-1)
RHO_ICE = 917.0

# ratio vaporization to sublimation
VAP_SUB = 2.501 / 2.835

# Convert calories to Joules
CAL_TO_J = 4.186798188

SNOW_EMISSIVITY = 0.98
STEF_BOLTZ = 5.67032e-8     # Stefan-Boltzmann constant (W / m^2 / deg^4)


# Maximum density due to compaction by gravity (kg/m^2)
RHO_MAX = 550

# R = 48 # in the original but not used?
R1 = 23.5
R2 = 24.5

SWE_MAX = 2000.0
