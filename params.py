'''
A parameter file for the Transient Sky Simulator
'''

# Parameters for the observation

nCells_D = 1000          # number of cells in the distance direction (galactic)
nCells_xgal = int(1e4)   # number of cells in the distance directions (extragalactic)


"""
Telescope capabilities

mag_limit:      limiting magnitude in each filter band
mag_resolution: Resolution of the telescope in magnitudes (in most sensitive band)
aperture_DEC:   Telescope aperture in DEC direction (in degrees)
aperture_RA:    Telescope aperture in RA direction (in degrees)
color_system:   The color system corresponding to ugriz(y)(q)
                Choose from 'sdss', 'lsst' or 'blackgem'
                The UBVRI color scheme is always loaded
"""

mag_limit      = {'U':21.0, 'B':21.0, 'V':21.0, 'R':21.0, 'I':21.0, 'u':20.9, 'g':22.9, 'r':22.3, 'i':21.7, 'z':21.4, 'q': 23.2, 'y':21.0} 
#mag_limit      = {'U':21.0, 'B':21.0, 'V':21.0, 'R':21.0, 'I':21.0, 'u':18.9, 'g':19.9, 'r':19.3, 'i':19.7, 'z':19.4} 
mag_resolution = 0.1	
aperture_DEC   = 1.64317
aperture_RA    = 1.64317
color_system   = 'blackgem'



# Transients to use. Set to zero to leave out of calculation.

use_nova     = 0
use_CVdwarf  = 0
use_AMCVn    = 0
use_SNIa     = 0
use_SNIb     = 0
use_SNIc     = 1
use_SNIIP    = 0
use_SNIIL    = 0
use_SNIInL   = 0
use_SNIInP   = 0
use_M3       = 0
use_M3_5     = 0
use_M4       = 0
use_kilonova = 1	 #Should have no entry in dataTransients.dat or peakmagsTransients.dat

"""
Other transients that you have added yourself can be put in this array.
They will ALWAYS be used
Please use the same name (with the same capitalization) in:
 - dataTransients.dat
 - peakmagsTransients.dat
"""
Extra_transients = []

#Data for the kilonova injection
maxLIGODist = 5.25e26    #cm = 170 Mpc
k_tMerge    = 0	         #s	Time of merger relative to observation start
k_dist	    = 5.25e24    #cm	Can be set to None to set the distance randomly within LIGO limits
k_mvRed	    = [0.04,0.15]#[Msun,c] Ejecta mass and velocity of red component
k_mvBlue    = [0.025,0.3]#[Msun,c] Ejecta mass and velocity of blue component
			 #To only have a red(blue) component, set k_mvRed(k_mvBlue)=None

# File holding general data for the transients:
transientFile = 'dataTransients.dat'
#transientFile = 'test.dat'

#File holding the peak magnitudes for the transients:
PeakMagFile = 'peakmagsTransients.dat'

# File holding data for M-dwarfs
MDwarfFile = 'data_MDwarfs.dat'


#Place to save the animation
Animation_loc = 'Animation.mp4'

#Output file
outfile = 'long.dat'
