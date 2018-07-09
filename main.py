import time
import sys
import os
import shutil
import numpy as np
import params as pm	
import observ as ob
import argparse
import transient as trns
from colorcolor import ColorColorSky
from Animation import AnimateSky

def main():
    tic = time.time()
    Opts = getOpts()
    if Opts.params:
        CreateNewParams(Opts.params)
    import params as pm      #Import the new paramaters as pm

    ###################################
    # Set up the observation:
    ###################################
    observingRun = ob.observation( pm.RA_lo, pm.RA_hi, pm.DEC_lo, pm.DEC_hi,\
                                    pm.nCells_D, pm.nCells_xgal, pm.mag_limit,\
                                    pm.ObsFile, Opts )

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
                                       pm.k_mvBlue, observingRun )
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

    if TotN_trans > 0:
        if Opts.option:
            if 'Animate' in Opts.option:
                AnimateSky()
            if 'ColorColor' in Opts.option:
                ColorColorSky()
        else:
            print "Neither 'Animate' or 'ColorColor' was given as an option. Saving the results to ", pm.outfile
    else:
        print "No transients were found, exiting simulation..."



def getOpts():
    parser = argparse.ArgumentParser(description='Simulate the transient sky')
    parser.add_argument("-o", "--option", nargs = '*', help="Execute an extra function after the generation of transients. Choose either 'Animate' or 'ColorColor'")
    parser.add_argument("-p", "--params", help="Define a file with observation parameters. If a previous params.py file exists, it'll be copied to params_old.py and the new file will overwrite params.py. If params_old.py exists, it will be overwritten.")
    parser.add_argument("-f", "--file", help="File with observation times. This can also be entered in params.py")
    args = parser.parse_args()
    return args

def CreateNewParams(Paramfile):
    if Paramfile != 'params_old.py':
        if os.path.exists('params_old.py'):
            print "Overwriting params_old.py with params.py"
        shutil.copy('params.py', 'params_old.py')
        print "Copying " + Paramfile + " into params.py"
        shutil.copy(Paramfile, 'params.py')
    else: 
        Message = "You're trying to use the params_old.py file for your observation parameters. Please rename this file and try again."
        raise Exception(Message)
        #sys.exit(0)


if __name__ == "__main__": 
    """
    Run as: python main.py [Arguments]
    Optional arguments:
    [-o] [--option]  'Animate' and/or 'ColorColor' 
    [-f] [--file]    The Obstimes file e.g. 'Obstimes.txt' This can also be entered in params.py
    [-p] [--params]  Use a different params file.
    [-h] [--help]    Print help function
    """
    print getOpts()
    main()

