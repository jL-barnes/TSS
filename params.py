'''
A parameter file for the Transient Sky Simulator
'''

# Parameters for the observation

n_RA = 6                 # number of cells in the RA direction
n_DEC = 6                # number of cells in the DEC direction
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

mag_limit      = {'U':23.0, 'B':23.0, 'V':21.1, 'R':20.9, 'I':20.2, 'u':20.9, 'g':22.9, 'r':22.3, 'i':21.7, 'z':21.4, 'q': 23.2, 'y':21.0} 
mag_resolution = 0.1	
aperture_DEC   = 1.64317
aperture_RA    = 1.64317
color_system   = 'sdss'



# Transients to use. Set to zero to leave out of calculation.

use_nova     = 1
use_UGem     = 1
use_SUUMa    = 1
use_ZCam     = 1
use_SNIa     = 1
use_SNIb     = 1
use_SNIc     = 1
use_SNIIP    = 1
use_SNIIL    = 1
use_SNIInL   = 1
use_SNIInP   = 1
use_M3       = 1
use_M3_5     = 1
use_M4       = 1
use_kilonova = 1	 #Should have no entry in dataTransients.dat or peakmagsTransients.dat

"""
Other transients that you have added yourself can be put in this array.
They will ALWAYS be used
Please use the same name (with the same capitalization) in dataTransients.dat
Example: Extra_transients = ['AM CVn']
"""
Extra_transients = []

#Data for the kilonova injection
maxLIGODist = 5.25e26    #cm = 170 Mpc
k_tMerge    = -3600      #s	Time of merger relative to observation start
k_dist	    = None    #cm	Can be set to None to set the distance randomly within LIGO limits
k_mvRed	    = [0.04,0.15]#[Msun,c] Ejecta mass and velocity of red component
k_mvBlue    = [0.025,0.3]#[Msun,c] Ejecta mass and velocity of blue component
			 #To only have a red(blue) component, set k_mvRed(k_mvBlue)=None
k_HGE       = 'no'       #The host galaxy extinction probability distribution scale
                         #no = no host galaxy extinction
                         #G%f (with %f a float) is a Gaussian with sigma = %f
                         #E%f (with %f a float) is an exponential with sigma = %f

# File holding general data for the transients:
transientFile = 'dataTransients.dat'
#transientFile = 'test.dat'

# File holding data for M-dwarfs
MDwarfFile = 'data_MDwarfs.dat'


#The following arguments can be overridden by parsing arguments to the python 
# programs on the command line. These are taken as defaults.

#The file with the dates and colorfilters of the observations.
Obstimes = 'Obstimes.dat'

#Colors used by Color-color diagram.
#If you fill in ['A', 'B', 'C'], this will create a plot for A-B against B-C
#If you enter ['A', 'B', 'C', 'D'], this will create a plot for A-B against C-D
#One should always enter 3 or 4 colors. 
#Please enter them in the order of observation.
CC_bands = ['g', 'r', 'i']

#Colors used for the animation
#An animation will be created with a separate subfigure for each color
Ani_bands = ['g', 'r']

#Output file for TSS data
outfile = 'output.dat'

#Output file for Color-Color diagram
CC_outfile = 'CC.png'

#Output file for the photometry animation
Ani_outfile = 'Animation.mp4'

Ani_samefr = False #If True the framerate between two epochs that are far apart in the animation does not change. Should usually be set to False.
