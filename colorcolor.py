import params as pm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from astLib import astCoords

class CCData:
    def __init__(self, datafile, bandList, fromfile):
        #global markerList
        f = open(datafile, 'r')
        typesList = f.readline().split()
        hdr = f.readline().split()
        f.close()
        data = np.genfromtxt( datafile, skip_header =2 )
        if np.array(data).shape == (6,):
            raise Exception("Too little data to make an animation!")
        else:
            self.mag_norm = np.min( data[:,5:])
        # set the window for the observation
        tm, self.IDs, ras, decs, self.bandObs = np.unique(data[:,2]), np.unique(data[:,0]), data[:,3], data[:,4], data[:,6]
        if not fromfile:
            self.ra_min   = astCoords.hms2decimal(pm.RA_lo,":")
            self.ra_max   = astCoords.hms2decimal(pm.RA_hi,":")
            self.dec_min  = astCoords.dms2decimal(pm.DEC_lo,":")
            self.dec_max     = astCoords.dms2decimal(pm.DEC_hi,":")
            if self.ra_max < self.ra_min: 
                self.ra_min, self.ra_max = self.ra_max, self.ra_min
            if self.dec_max < self.dec_min: 
                self.dec_min, self.dec_max = self.dec_max, self.dec_min
        else:	#pm.RA_lo may not refer to the correct coordinates if read from file
            self.ra_min = np.min(ras)
            self.ra_max = np.max(ras)
            self.dec_min = np.min(decs)
            self.dec_max = np.max(decs)          
        self.nBands = 3
        col_idstypes = np.array([1,2])
        col_idsbands = np.linspace( 5 , 5 + self.nBands - 1 , self.nBands ).astype(int) 
        #col_ids = np.append(col_idstypes, col_idsbands)
        col_ids = np.array([0,1,3,4,5])
        #self.TRdata = [ data[data[:,0] == ID][:,col_ids] for ID in IDs ]
        self.time_data  = []
        self.mag_data   = []
        self.band_data  = []
        self.trtype_data= []
        for ID in IDs:
            Bool = data[:,0] == ID
            self.time_data.append( data[Bool, 2] )
            self.mag_data.append( data[Bool, 5] )
            self.band_data.append( data[Bool, 6] )
            self.trtype_data.append( data[0, 1])
        self.times = tm
        dt = (tm[1] - tm[0])/8.64e4
        self.n_frames = len(self.times)
        nBands = len(bandList)
        # figure out the marker stuff
        self.trtypes = [int(ty) for ty  in typesList[::2]]
        trNames = typesList[1::2]
        #if len(trtypes) > len(markerList):
        #    markerList += markerList[0: len(trtypes)-len(markerList)]
        self.nTrTypes = len(self.trtypes)
        self.trNames = trNames

class ColorColor:
    def __init__(self, dataFile, bandlist, fromfile):
        self.data = CCData(dataFile, bandlist, fromfile)
        self.mag_lim = pm.mag_limit
        self.bands = bandlist
        self.Obstimes = self.data.times
        self.Alphas = np.linspace(0.85,0.2,len(self.Obstimes))
    def plotCC(self):
        cmap = plt.get_cmap('Accent')
        DotColor = [cmap(i//1./self.data.nTrTypes) for i in self.data.trtypes]
        Size = 15
        IDs_seen = []
        SNIa = []
        SNIb = []
        SNIc = []
        SNIIL = []
        SNIIP = []
        SNIInP = []
        SNIInL = []
        CVdwarf = []
        #print np.array(self.data.TRdata[1]).shape
        #for i,ID in enumerate(self.data.TRdata):
        for i,ID in enumerate(self.data.IDs):
            #print i
            #print ID[:,2]	#[i,j]   i=Obs j=0,1,2,3,4: type, time,U,B,V 
            tt = int(ID[0,0])
            tt_index = 0
            #for l,k in enumerate(self.data.trtypes):
            #    if k == tt: tt_index = l
            times  = np.array(ID[:,1])
            BandM0 = np.array(ID[:,2])
            BandM1 = np.array(ID[:,3])
            BandM2 = np.array(ID[:,4])
            See0 = BandM0 < self.mag_lim[self.bands[0]]	#What can be seen in color 0?
            See1 = BandM1 < self.mag_lim[self.bands[1]]
            See2 = BandM2 < self.mag_lim[self.bands[2]]
            See = See0 * See1 * See2	#What can be seen in all colors?
            if [0.] in times: IDs_seen.append(i)
            elif [259200.0] in times:
                if self.data.trNames[tt_index] == 'SNIa': SNIa.append(i)
                if self.data.trNames[tt_index] == 'SNIb': SNIb.append(i)
                if self.data.trNames[tt_index] == 'SNIc': SNIc.append(i)
                if self.data.trNames[tt_index] == 'SNIIL': SNIIL.append(i)
                if self.data.trNames[tt_index] == 'SNIIP': SNIIP.append(i)
                if self.data.trNames[tt_index] == 'SNIInL': SNIInL.append(i)
                if self.data.trNames[tt_index] == 'SNIInP': SNIInP.append(i)
                if self.data.trNames[tt_index] == 'CVdwarf': CVdwarf.append(i)
            if np.any(See):		#The transient can be seen in all colors
                rgbacolors = np.zeros((len(times[See]),4))
                rgbacolors[:,:3] = DotColor[tt_index][:3]
                alphas = np.zeros(len(self.Obstimes), dtype = bool)
                for i in range(len(alphas)):
                    if self.Obstimes[i] in times[See]: alphas[i] = True
                rgbacolors[:,3] = self.Alphas[alphas]
                if self.data.trNames[tt_index] == 'kilonova': Size = 125
                else: Size = 10
                if self.data.trNames[tt_index] == 'SNIa': 
                    Zor = 0
                else: Zor = tt_index + 1
                plt.scatter(BandM1[See] - BandM2[See], BandM0[See] - BandM1[See], color = rgbacolors, edgecolors = (1,1,1,0), s = Size, zorder = Zor)  
        legendhandles = []
        print "Nr of SNIa discovered: ", len(SNIa)
        print "Nr of SNIb discovered: ", len(SNIb)
        print "Nr of SNIc discovered: ", len(SNIc)
        print "Nr of SNIIL discovered: ", len(SNIIL)
        print "Nr of SNIIP discovered: ", len(SNIIP)
        print "Nr of SNIInL discovered: ", len(SNIInL)
        print "Nr of SNIInP discovered: ", len(SNIInP)
        print "Nr of CVdwarf discovered: ", len(CVdwarf)
        print "Nr of transients that were already visible: ", len(IDs_seen)
        for i, Transient in enumerate(self.data.trNames):
            legendhandles.append( mpatches.Patch(color=DotColor[i], label=Transient))
        plt.xlabel('%s - %s' %(self.bands[1], self.bands[2]))
        plt.ylabel('%s - %s' %(self.bands[0], self.bands[1]))
        plt.legend(handles=legendhandles)
        plt.show()



def ColorColorSky(bList, skyfile):
    """
    Make a color-color plot using the data that is generated with the main.py script 
    We take the the first three colors in showbands in params.py
    """
    print "Drawing Color-color plot"
    cc = ColorColor( skyfile, bList, fromfile=False )
    cc.plotCC()

def ColorColorFile( bList, skyfile ):
    """
    Make a color-color plot using the data from the 'outfile' in params.py
    We take the the first three colors in showbands in params.py. These have to correspond to
     the colors used in the data file.
    """
    cc = ColorColor( skyfile, bList, fromfile = True )
    cc.plotCC()

