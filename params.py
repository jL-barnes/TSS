'''
A parameter file for the Transient Sky Simulator
'''

# Parameters for the observation

# size of the window

RA_lo  = "14:10:00" 
RA_hi  = "14:30:00"
DEC_lo = "20:00:00"
DEC_hi = "25:00:00"
nCells_D = 1000          # number of cells in the distance direction (galactic)
nCells_xgal = int(1e4)   # number of cells in the distance directions (extragalactic)

# time of observation

# eventually: change this to a date so we can track the moon
start_obs = 0.0
dur_obs = 20             # duration of observation in days 
nObs_per_day = 6         # number of observations per day 

# telescope capabilities

mag_limit = 23.0           # limiting magnitude in most sensitive band
color_system = 'UBVRI'

# factor (Figure out what this is...)
factor = 1.0

# Transients to use. Set to zero to leave out of calculation.

use_nova     = 1
use_CVpolar  = 0
use_CVIP     = 0
use_CVdwarf  = 0
use_AMCVn    = 0
use_SNIa     = 1
use_M3       = 0
use_M3_5     = 0
use_M4       = 0
use_kilonova = 0

# File holding data for the transients:
transientFile = 'dataTransients.dat'





# parameters for the data display and output
outfile = 'test_SNIa.dat'
# band to show in animation. Must be consistent with choice of color system above.
showbands = 'BRI'  
