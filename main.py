import time
import sys
import numpy as np
import params as pm
import observ as ob
import transient as trns
from colorcolor import ColorColorSky
from Animation import AnimateSky



def main():
    tic = time.time()
    ###################################
    # Set up the observation:
    ###################################
    observingRun = ob.observation( pm.RA_lo, pm.RA_hi, pm.DEC_lo, pm.DEC_hi,\
                                    pm.nCells_D, pm.nCells_xgal, pm.mag_limit,\
                                    pm.start_obs, pm.dur_obs, pm.nObs_per_day )

    ###################################
    # Set up transient templates
    ###################################

    transientZoo = trns.TransientSet( observingRun )
    for templ in transientZoo.transient_templates:
        print templ.tag

    ###################################
    # Generate transients
    ###################################
    transientZoo.populate()

    if pm.use_kilonova:
        transientZoo.inject_kilonova( pm.k_tMerge, pm.k_dist, pm.k_mvRed,\
                                       pm.k_mvBlue )
        print "Injecting a kilonova"
    tuc = time.time()

    ###################################
    # Make observations
    ###################################
    observingRun.take_data( transientZoo, pm.outfile )
    tec = time.time()
    print "Time for observations", tec - tuc

    TotN_trans = 0
    for temp in transientZoo.transient_templates:
        print temp.tag, temp.N_trans, temp.N_transnonsee
        TotN_trans += temp.N_trans

    ###################################
    # Do the animation
    ###################################
    toc = time.time()
    print "Total time elapsed", toc - tic

    Opts = getopts(sys.argv[1:])
    if TotN_trans > 0 and Opts != {}:
        if 'Animate' in Opts['-o']:
            AnimateSky()
        if 'ColorColor' in Opts['-o']:
            ColorColorSky()
        else:
            print "Neither 'Animate' or 'ColorColor' was given as an option. Saving the results in ", pm.outfile
    else:
        print "No transients were found, exiting simulation..."

def getopts(argv):
    opts = {}
    while argv:
        if argv[0][0] == '-':
            if argv[0] in opts:
                opts[argv[0]].append(argv[1])
            else:
                opts[argv[0]] = [argv[1]]
        argv = argv[1:] 
    return opts

if __name__ == "__main__": 
    """
    Run as: python main.py [-o Animate] [-o ColorColor]
    """
    main()

