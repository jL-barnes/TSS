import time
import observ as ob
import numpy as np
import argparse
import transient as trns
from colorcolor import ColorColorSky
from Animation import AnimateSky

def main(Opts):
    tic = time.time()

    #Import the paramaters into a pm dict
    pm = {}
    execfile(Opts.params, pm)

    ###################################
    # Set up the observation:
    ###################################
    observingRun = ob.observation( pm, Opts )

    #pm['RA_lo'], pm['RA_hi'], pm['DEC_lo'], 
                                    #pm['DEC_hi'], pm.nCells_D, pm.nCells_xgal, pm.mag_limit,
    ###################################
    # Set up transient templates
    ###################################
    Set_of_transients = []
    for framenr in range(observingRun.Nframes):
        transientZoo = trns.TransientSet( observingRun, framenr )
        for templ in transientZoo.transient_templates:
            print templ.tag

        ###################################
        # Generate transients
        ###################################
        transientZoo.populate()
        Set_of_transients.append(transientZoo)

    if pm['use_kilonova']:
        SKN = np.random.randint(0,len(Set_of_transients))
        Set_of_transients[SKN].inject_kilonova( pm['k_tMerge'], pm['k_dist'],\
                                                pm['k_mvRed'], pm['k_mvBlue'],\
                                                observingRun, SKN )
        print "Injecting a kilonova"
    tuc = time.time()

    ###################################
    # Make observations
    ###################################
    TotalData = []
    TrTypes   = []
    TrTypenrs = []
    No_of_trs = 0
    for T_zoo in Set_of_transients:
        Data, TrTs, TrTnrs, NrTrs = observingRun.take_data( T_zoo, pm['outfile'], TrTypes, TrTypenrs, No_of_trs )
        #print Data
        TotalData.append(Data)
        TrTypes = TrTs
        TrTypenrs = TrTnrs
        No_of_trs = NrTrs
    #print TotalData
    print TrTypes, "trtypes"
    tec = time.time()
    print "Time for observations", tec - tuc

    TotN_trans = 0
    writefile = observingRun.OpenFile(pm['outfile'], TrTypes, observingRun.Trlegend)
    for i,T_zoo in enumerate(Set_of_transients):
        observingRun.WriteToFile(writefile, TotalData[i])
        #for temp in transientZoo.transient_templates:
        #    print temp.tag, temp.N_trans, temp.N_transnonsee
        #    TotN_trans += temp.N_trans
    for tag in TrTypes:
        N_tr_seen = 0
        N_tr_not_seen = 0
        for T_zoo in Set_of_transients:
            for templ in T_zoo.transient_templates:
                if templ.tag == tag:
                    N_tr_seen += templ.N_trans    
                    N_tr_not_seen += templ.N_transnonsee -templ.N_trans
        print tag, "Seen: ", N_tr_seen, " Not seen: ", N_tr_not_seen
        TotN_trans += N_tr_seen

    ###################################
    # Do the animation
    ###################################
    toc = time.time()
    print "Total time elapsed", toc - tic

    if TotN_trans > 0:
        if Opts.option:
            if 'Animate' in Opts.option:
                AnimateSky(observingRun.bands ,pm['outfile'])
            if 'ColorColor' in Opts.option:
                ColorColorSky(observingRun.bands, pm['outfile'])
        else:
            print "Neither 'Animate' or 'ColorColor' was given as an option. Saving the results to ", pm['outfile']
    else:
        print "No transients were found, exiting simulation..."



def getOpts():
    parser = argparse.ArgumentParser(description='Simulate the transient sky')
    parser.add_argument("-o", "--option", nargs = '*', help="Execute an extra function after the generation of transients. Choose either 'Animate' or 'ColorColor'")
    parser.add_argument("-p", "--params", default = 'params.py', help="Define a file with observation parameters. default:params.py")
    parser.add_argument("-f", "--file", default='Obstimes.dat', help="File with observation times. default:Obstimes.dat")
    parser.add_argument("-l", "--offline", action='store_true', help ="Execute the calculations fully offline. This usually takes a lot longer!")
    parser.add_argument("-d", "--nodust", action='store_true', help ="Exclude dust in the calculations")
    parser.add_argument("-c", "--colorsys", choices = ["AB", "Vega", "ABVega"], default = 'AB', help="Color system to use. ABVega will return UBVRI measurements in Vega and ugriz measurements in AB. Default: AB")
    parser.set_defaults(feature=True)
    args = parser.parse_args()
    return args

def printOpts(Opts):
    print "Running TSS with options:"
    if Opts.option:
        print "[-o] [--option]  ", Opts.option 
    if Opts.file:
        print "[-f] [--file]    ", Opts.file
    if Opts.params:
        print "[-p] [--params]  ", Opts.params
    if Opts.offline:
        print "[-l] [--offline] ", Opts.offline
    if Opts.nodust:
        print "[-d] [--nodust]  ", Opts.nodust
    if Opts.colorsys:
        print "[-c] [--colorsys]", Opts.colorsys

if __name__ == "__main__": 
    """
    Run as: python main.py [Arguments]
    Optional arguments:
    [-o] [--option]   'Animate' and/or 'ColorColor' 
    [-f] [--file]     The Obstimes file default: 'Obstimes.dat' 
    [-p] [--params]   Params file to use. default:params.py
    [-l] [--offline]  Execute offline
    [-d] [--nodust]   Exlude dust
    [-c] [--colorsys] Color system to use. ABVega will return UBVRI measurements in Vega and ugriz measurements in AB. Default='AB'
    [-h] [--help]     Print help function
    """
    Opts = getOpts()
    printOpts(Opts)
    
    main(Opts)

