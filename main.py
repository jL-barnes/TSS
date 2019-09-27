import time
import observ as ob
import numpy as np
import argparse
import transient as trns
from colorcolor import ColorColorSky
from Animation import AnimateSky

def main(Opts, pm):
    tic = time.time()

    ###################################
    # Set up the observation:
    ###################################
    observingRun = ob.observation(pm, Opts)

    ###################################
    # Set up transient templates
    ###################################

    Set_of_transients = []
    for framenr in range(observingRun.Nframes):
        transientZoo = trns.TransientSet(observingRun, framenr)
        for templ in transientZoo.transient_templates:
            print(templ.tag)

        ###################################
        # Generate transients
        ###################################
        transientZoo.populate()
        Set_of_transients.append(transientZoo)
    
    #Generate a kilonova
    if pm['use_kilonova']:
        SKN = np.random.randint(0,len(Set_of_transients))
        Set_of_transients[SKN].inject_kilonova(pm['k_tMerge'], pm['k_dist'],
                                               pm['k_mvRed'], pm['k_mvBlue'],
                                               observingRun, SKN )
        print("Injecting a kilonova")
    tuc = time.time()
    print("Time for generation", tuc - tic)

    ###################################
    # Make observations
    ###################################
    TotalData = []
    TrTypes   = []
    TrTypenrs = []
    No_of_trs = 0
    for T_zoo in Set_of_transients:
        Data, TrTs, TrTnrs, NrTrs = observingRun.take_data(T_zoo, 
            pm['outfile'], TrTypes, TrTypenrs, No_of_trs)
        TotalData.append(Data)
        TrTypes = TrTs
        TrTypenrs = TrTnrs
        No_of_trs = NrTrs

    tec = time.time()
    print("Time for observations", tec - tuc)

    TotN_trans = 0
    writefile = observingRun.OpenFile(Opts.output, TrTypes, 
                                      observingRun.Trlegend)
    for i,T_zoo in enumerate(Set_of_transients):
        observingRun.WriteToFile(writefile, TotalData[i])
    observingRun.CloseFile()

    for tag in TrTypes:
        N_tr_seen = 0
        N_tr_not_seen = 0
        for T_zoo in Set_of_transients:
            for templ in T_zoo.transient_templates:
                if templ.tag == tag:
                    N_tr_seen += templ.N_trans    
                    N_tr_not_seen += templ.N_transnonsee -templ.N_trans
        print(tag, "Seen: ", N_tr_seen, " Not seen: ", N_tr_not_seen)
        TotN_trans += N_tr_seen

    ###################################
    # Do the animation
    ###################################
    toc = time.time()
    print("Total time elapsed", toc - tic)

    if TotN_trans > 0:
        if Opts.option:
            if 'Animate' in Opts.option:
                AnimateSky(observingRun.bands, Opts, pm)
            if 'ColorColor' in Opts.option:
                ColorColorSky(observingRun.bands, Opts, pm)
            else:
                print("Warning: A wrong option was given with the argument",\
                      " --option.")
                print("It should be either 'Animate' or 'ColorColor' ")
        else:
            print("Neither 'Animate' or 'ColorColor' was given as an option.",\
                  " Saving the results to ", pm['outfile'])
    else:
        print("No transients were found, exiting simulation...")

def getOpts():
    """ Get the running arguments
    
    Returns
    -------
    A parser object with the arguments with which to run Animation.py
    """
    paramshelp = "Define a file with observation parameters. default:params.py"
    filehelp = ("File with observation times. defaults to the obstimes in"
                + " params.py.")
    offlinehelp = ("Execute the calculations fully offline. This usually takes"
                   + " a lot longer!"
    nodusthelp = "Exclude dust in the calculations"
    colorsyshelp = ("Color system to use. ABVega will return UBVRI"
                    + " measurements in Vega and ugriz measurements in AB."
                    + " The other options return all measurements in either AB" 
                    + "or Vega Default: ABVega")
    outputhelp = "Output file. Defaults to the outputfile in params.py"
    optionhelp = ("Execute an extra function after the generation of"
                  + " transients. Choose either 'Animate' or 'ColorColor'")

    parser = argparse.ArgumentParser(description='Simulate the transient sky')
    parser.add_argument("-p", "--params", default='params.py', 
                        help=paramshelp)
    parser.add_argument("-f", "--file", help=filehelp)
    parser.add_argument("-l", "--offline", action='store_true',
                        help=offlinehelp)
    parser.add_argument("-d", "--nodust", action='store_true', help=nodusthelp)
    parser.add_argument("-c", "--colorsys", choices = ["AB", "Vega", "ABVega"], 
                        default = 'ABVega', help=colorsyshelp)
    parser.add_argument("-o", "--output", help=outputhelp)
    parser.add_argument("-O", "--option", nargs = '*', help=optionhelp)
    parser.set_defaults(feature=True)
    args = parser.parse_args()
    return args

def printOpts(Opts):
    """ Prints the parser arguments with which to run Animation.py

    Parameters
    ----------
    Opts : argparse object
        The arguments parsed to the Animation.py program
    """
    print("Running TSS with options:")
    if Opts.file:
        print("[-f] [--file]    ", Opts.file)
    if Opts.params:
        print("[-p] [--params]  ", Opts.params)
    if Opts.offline:
        print("[-l] [--offline] ", Opts.offline)
    if Opts.nodust:
        print("[-d] [--nodust]  ", Opts.nodust)
    if Opts.colorsys:
        print("[-c] [--colorsys]", Opts.colorsys)
    if Opts.output:
        print("[-o] [--output]  ", Opts.output)
    if Opts.option:
        print("[-O] [--option]  ", Opts.option)

if __name__ == "__main__": 
    """
    Run as: python3 main.py [Arguments]

    Optional arguments:
    [-f] [--file]     The Obstimes file default: 'Obstimes.dat' 
    [-p] [--params]   Params file to use. default:params.py
    [-l] [--offline]  Execute offline
    [-d] [--nodust]   Exlude dust
    [-c] [--colorsys] Color system to use. ABVega will return UBVRI 
                       measurements in Vega and ugriz measurements in 
                       AB. Default='ABVega'
    [-h] [--help]     Print help function
    [-O] [--output]   Output file. Defaults to the outputfile in params.py
    [-o] [--option]   'Animate' and/or 'ColorColor' 
    These arguments will override any arguments in params.py
    """
    Opts = getOpts()
    
    #Import the paramaters into a pm dict
    pm = {}
    with open(Opts.params) as f:
       code = compile(f.read(), Opts.params, 'exec')
       exec(code, pm)
    pm['filename'] = Opts.params
    if not Opts.output:
        Opts.output = pm['outfile']
    if not Opts.file:
        Opts.file = pm['Obstimes']
    printOpts(Opts)
    
    main(Opts, pm)

