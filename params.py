'''
A parameter file for the Transient Sky Simulator
'''

# Parameters for the observation

# size of the window (in JD2000)
RA_lo  = "14:10:00" 
RA_hi  = "14:30:00"
DEC_lo = "20:00:00"
DEC_hi = "25:00:00"
"""

RA_lo  = "14:10:00" 
RA_hi  = "14:30:00"
DEC_lo = "20:00:00"
DEC_hi = "25:00:00"

RA_lo  = "17:43:40" 
RA_hi  = "17:47:40"
DEC_lo = "-27:00:28"
DEC_hi = "-31:00:28"

RA_lo  = "17:43:40" 
RA_hi  = "17:47:40"
DEC_lo = "-37:00:28"
DEC_hi = "-41:00:28"
"""


nCells_D = 1000          # number of cells in the distance direction (galactic)
nCells_xgal = int(1e4)   # number of cells in the distance directions (extragalactic)

# time of observation

# eventually: change this to a date so we can track the moon
start_obs = 0.0
dur_obs = 80             # duration of observation in days 
nObs_per_day = 1         # number of observations per day 

# telescope capabilities

mag_limit = {'U':23.0, 'B':23.0, 'V':23.0, 'R':23.0, 'I':23.0, 'J':23.0, 'H':23.0, 'K':23.0, 'u':23.0, 'g':23.0, 'r':23.0, 'i':23.0, 'z':23.0}           # limiting magnitude in each band
mag_resolution = 0.1	   #Resolution of the telescope in magnitudes (in most sensitive band)
color_system = 'UBVRI'#'ugriz'


# Transients to use. Set to zero to leave out of calculation.

use_nova     = 0
use_CVpolar  = 0
use_CVIP     = 0
use_CVdwarf  = 0
use_AMCVn    = 0
use_SNIa     = 0
use_SNIb     = 1
use_SNIc     = 1
use_SNIIL    = 0
use_SNIIP    = 0
use_M3       = 0
use_M3_5     = 0
use_M4       = 0
use_kilonova = 0	#Should have no entry in dataTransients.dat or peakmagsTransients.dat

#Data for the kilonova injection
maxLIGODist = 5.25e26    #cm = 170 Mpc
k_tMerge    = 0	         #s	Time of merger relative to observation start
k_dist	    = 5.25e24    #cm	Can be set to None to set the distance randomly within LIGO limits
k_mvRed	    = [0.04,0.15]#[Msun,c] Ejecta mass and velocity of red component
k_mvBlue    = [0.025,0.3]#[Msun,c] Ejecta mass and velocity of red component
			 #To only have a red(blue) component, set k_mvRed(k_mvBlue)=None

# File holding data for the transients:
transientFile = 'dataTransients.dat'
#transientFile = 'test.dat'

#File holding the peak magnitudes for the transients:
PeakMagFile = 'peakmagsTransients.dat'

# File holding data for M-dwarfs
MDwarfFile = 'data_MDwarfs.dat'


#Place to save the animation
Animation_loc = 'Leuven.mp4'

# parameters for the data display and output
outfile = 'test_SNIa.dat'
# band to show in animation. Must be consistent with choice of color system above.
showbands = 'UBVRI'  
# should we use dust extinction in the calculations?
use_dust = True
