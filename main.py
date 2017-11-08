import numpy as np
#import params as pm
import testparams as pm
import observ as ob
import transient as trns
from Animation import AnimateSky
#import MilkyWay as MW



###################################
# Set up the observation:
###################################

#trns.SetUpTransientData()

###################################
# Set up the observation:
###################################
#%%
observingRun = ob.observation( pm.RA_lo, pm.RA_hi, pm.DEC_lo, pm.DEC_hi,\
                                pm.nCells_D, pm.nCells_xgal, pm.mag_limit,\
                                pm.start_obs, pm.dur_obs, pm.nObs_per_day )
#%%
observingRun.setUpColors( pm.color_system, pm.showbands )                               


###################################
# Set up transients
###################################

transientZoo = trns.TransientSet( observingRun )
for templ in transientZoo.transient_templates:
    print templ.tag

transientZoo.populate()

print templ.tag, templ.N_trans

observingRun.take_data( transientZoo, pm.outfile )


for temp in transientZoo.transient_templates:
    print temp.tag, temp.N_trans, temp.N_transnonsee
print "To do:\n\tPhillips relation\n\tDust Extinction\n\tGive each Mdwarf a good quiescent luminosity in ergs"

#AnimateSky()


#Construct_Animation()
#Construct_ColorColor()

#Show_Animation()
#Show_ColorColor()

####################################
# Get the light curves
####################################

#takeData( transientZoo )

####################################
# Plot the light curves
####################################

