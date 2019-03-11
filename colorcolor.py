import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from astLib import astCoords
import argparse

class CCData:
    def __init__(self, datafile, bandList):
        #global markerList
        f = open(datafile, 'r')
        typesList = f.readline().split()
        hdr = f.readline().split()
        f.close()
        #data = np.genfromtxt( datafile, skip_header =2 )
        Data = np.genfromtxt( datafile, skip_header =2 , dtype = [int, int, float, float, float, float, ('string','S5')])
        self.bandObs = np.array([row[6] for row in Data])
        data = np.array( [[row[0], row[1], row[2], row[3], row[4], row[5]] 
                          for row in Data] )
        if np.array(data).shape == (6,):
            raise Exception("Too little data to make an animation!")
        else:
            self.mag_norm = np.min( data[:,5:])
        # set the window for the observation
        tm, self.IDs, ras, decs = np.unique(data[:,2]), np.unique(data[:,0]), data[:,3], data[:,4]
        self.ra_min = np.min(ras)
        self.ra_max = np.max(ras)
        self.dec_min = np.min(decs)
        self.dec_max = np.max(decs)          
        self.nBands = 3
        #col_idstypes = np.array([1,2])
        #col_idsbands = np.linspace( 5 , 5 + self.nBands - 1 , self.nBands ).astype(int) 
        #col_ids = np.append(col_idstypes, col_idsbands)
        #col_ids = np.array([0,1,3,4,5])
        #self.TRdata = [ data[data[:,0] == ID][:,col_ids] for ID in IDs ]
        self.time_data  = []
        self.mag_data   = []
        self.band_data  = []
        self.trtype_data= []
        for ID in self.IDs:
            Bool = data[:,0] == ID
            self.time_data.append( data[Bool, 2] )
            self.mag_data.append( data[Bool, 5] )
            self.band_data.append( self.bandObs[Bool] )
            self.trtype_data.append( data[Bool][0, 1])
        for C in Opts.colors:
            uniq_bands = np.unique(self.bandObs)
            if C not in uniq_bands:
                raise ValueError("The wrong colors were requested. %s was not observed in %s which only contains observations in %s, %s and %s" % (C, Opts.file, uniq_bands[0], uniq_bands[1], uniq_bands[2]) )
        self.times = tm
        #dt = (tm[1] - tm[0])/8.64e4
        self.n_frames = len(self.times)
        #nBands = len(bandList)
        # figure out the marker stuff
        self.trtypes = [int(ty) for ty  in typesList[::2]]
        trNames = typesList[1::2]
        self.nTrTypes = len(self.trtypes)
        self.trNames = trNames
        self.trtypenr = {self.trtypes[i]:self.trNames[i] for i in range(self.nTrTypes)}

class ColorColor:
    def __init__(self, Opts, pm):
        self.data = CCData(Opts.file, Opts.colors)
        self.mag_lim = pm['mag_limit']
        self.bands = Opts.colors
        self.Obstimes = self.data.times
        self.Alphas = np.linspace(0.85,0.1,len(self.Obstimes))
    def plotCC(self):
        cmap = plt.get_cmap('Accent')
        DotColor = {i:cmap(i//1./self.data.nTrTypes) for i in self.data.trtypes}
        Size = 15
        """
        SNIa = []
        SNIb = []
        SNIc = []
        SNIIL = []
        SNIIP = []
        SNIInP = []
        SNIInL = []
        CVdwarf = []
        """
        #print np.array(self.data.TRdata[1]).shape
        #for i,ID in enumerate(self.data.TRdata):
        #print np.unique(self.data.trtype_data)
        #print self.data.trtypes
        #print self.data.IDs
        self.Nr_seen = {i : 0 for i in self.data.trtypes}
        for i,ID in enumerate(self.data.IDs):
            #print i
            #tt = int(ID[0,0])
            tt_index = int(self.data.trtype_data[i])
            #for l,k in enumerate(self.data.trtypes):
            #    if k == tt: tt_index = l
            B1_B0 = []
            B2_B1 = []
            times = []
            j0 = 0
            found0 = False
            obs_bands = self.data.band_data[i]
            obs_mags = self.data.mag_data[i]
            while not found0:
                if obs_bands[j0] == self.bands[0]:
                    B0 = obs_mags[j0]
                    entry = j0
                    if B0 < self.mag_lim[ self.bands[0] ]:
                        found1 = False
                        while not found1:
                            j0+=1
                            if j0 == len(obs_bands):
                                found1 = True
                                break
                            if obs_bands[j0] == self.bands[1]:
                                B1 = obs_mags[j0]
                                if B1 > self.mag_lim[ self.bands[1] ]:
                                    found1 = True
                                else: 
                                    found2 = False
                                    while not found2:
                                        j0+=1
                                        if j0 == len(obs_bands):
                                            found2 = True
                                            break
                                        if obs_bands[j0] == self.bands[2]:
                                            B2 = obs_mags[j0]
                                            if B2 > self.mag_lim[ self.bands[2] ]:
                                                found1 = True
                                            else:
                                                B2_B1.append(B2 - B1)
                                                B1_B0.append(B1 - B0)
                                                times.append(self.data.time_data[i][entry])
                                                found2 = True
                                                found1 = True
                                if j0 == len(obs_bands):
                                    found2 = True
                                    break
                    j0 = entry
                j0 += 1
                if j0 == len(obs_bands):
                    found0 = True
                    break
                
            rgbacolors = np.zeros((len(B1_B0),4))
            rgbacolors[:,:3] = DotColor[tt_index][:3]

            for m in range(len(B1_B0)):
                for a in range(len(self.Alphas)):
                    if self.Obstimes[a] == times[m]: 
                        rgbacolors[m,3] = self.Alphas[a]
            if self.data.trtypenr[tt_index] == 'kilonova': Size = 125
            else: Size = 10
            if self.data.trtypenr[tt_index] == 'SNIa': 
                Zor = 0      #Put SNIa in the background: they are often very numerous
            else: Zor = tt_index + 1
            if len(B2_B1) > 0:
                self.Nr_seen[tt_index] += 1   #One transient of this type was seen.
            #plt.scatter(BandM1[See] - BandM2[See], BandM0[See] - BandM1[See], color = rgbacolors, edgecolors = (1,1,1,0), s = Size, zorder = Zor)  
            plt.scatter(B2_B1, B1_B0, color = rgbacolors, edgecolors = (1,1,1,0), s = Size, zorder = Zor)  

        legendhandles = []
        
        print "Nr of transients seen:"
        for tt_index in self.data.trtypes:
            print self.data.trtypenr[tt_index], ":", self.Nr_seen[tt_index]

        for i, Transient in enumerate(self.data.trNames):
            tt_index = int(self.data.trtypes[i])
            C = DotColor[tt_index]
            legendhandles.append( mpatches.Patch(color=C, label=Transient))
        plt.xlabel('%s - %s' %(self.bands[1], self.bands[2]))
        plt.ylabel('%s - %s' %(self.bands[0], self.bands[1]))
        plt.legend(handles=legendhandles)
        plt.savefig(Opts.output)



def ColorColorSky(Opts, pm):
    """
    Make a color-color plot using the data that is generated with the main.py script 
    """
    print "Drawing Color-color plot"
    cc = ColorColor( Opts, pm )
    cc.plotCC()

def getOpts():
    parser = argparse.ArgumentParser(description='Animate an output of TSS')
    parser.add_argument("-c", "--colors", nargs = '*', default = ['g', 'r', 'i'], help="The colors/passbands to create a color-color plot for. If you fill in A B C, this will create a plot for A-B against B-C. One should always enter 3 colors. Please enter them in the order of observation.")
    parser.add_argument("-p", "--params", default = 'params.py', help="Define a file with observation parameters. default:params.py")
    parser.add_argument("-f", "--file", help="The file (that was output by TSS) to make a color-color diagram of. default: outfile in the params file")
    parser.add_argument("-o", "--output", default = 'CC.png', help ="The output file. default: CC.png")
    parser.set_defaults(feature=True)
    args = parser.parse_args()
    return args

def checkOpts(Opts):
    if len(Opts.colors) != 3:
        raise ValueError('The colors argument did not have 3 passbands entered, but ', len(Opts.colors))

def printOpts(Opts):
    print "Running TSS with options:"
    if Opts.colors:
        print "[-c] [--colors]  ", Opts.colors
    if Opts.params:
        print "[-p] [--params]  ", Opts.params
    if Opts.file:
        print "[-f] [--file]    ", Opts.file
    if Opts.output:
        print "[-o] [--output] ", Opts.output
    
if __name__ == "__main__": 
    """
    Run as: python main.py [Arguments]
    Optional arguments:
    [-c] [--colors]    The colors/passbands to create a color-color plot for. If you fill in A B C, this will create a plot for A-B against B-C. Please enter them in the order of observation.
    [-p] [--params]    Params file to use. default:params.py
    [-f] [--file]      The file (that was output by TSS) to make a color-color plot from. default: outfile in the params file
    [-o] [--output]    The output file. default: CC.png
    """
    Opts = getOpts()
    printOpts(Opts)
    
    pm = {}
    execfile(Opts.params, pm)
    if not Opts.file:
        Opts.file = pm['outfile']
    checkOpts(Opts)
        
    ColorColorSky(Opts, pm)


    #ColorColorSky(['g','r', 'i'], 'long2.dat')

