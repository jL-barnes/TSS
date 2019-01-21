import time
import observ as ob
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
        Set_of_transients[0].inject_kilonova( pm['k_tMerge'], pm['k_dist'],\
                                              pm['k_mvRed'], pm['k_mvBlue'],\
                                              observingRun )
        print "Injecting a kilonova"
    tuc = time.time()

    ###################################
    # Make observations
    ###################################
    TotalData = []
    for T_zoo in Set_of_transients:
        Data, TrTypes = observingRun.take_data( T_zoo, pm['outfile'] )
        TotalData.append(Data)
    print TrTypes, "trtypes"
    tec = time.time()
    print "Time for observations", tec - tuc

    TotN_trans = 0
    writefile = observingRun.OpenFile(pm['outfile'], TrTypes)
    for i,T_zoo in enumerate(Set_of_transients):
        observingRun.WriteToFile(writefile, TotalData[i])
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


"""
def CreateNewParams(Paramfile):
    if Paramfile != 'params.py':
        i = 0
        while os.path.exists('params_old%d.py' % (i)):
            i+=1
        print "Copying params.py into params_old%d.py" % (i)
        shutil.copy('params.py', 'params_old%d.py' % (i))
        print "Copying " + Paramfile + " into params.py"
        shutil.copy(Paramfile, 'params.py')
"""

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

