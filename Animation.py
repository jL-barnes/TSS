#matplotlib.use( 'macosx' )
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import ColorConverter
from matplotlib import rcParams
import argparse

rcParams.update({'figure.autolayout':True})
#rc('font',**{'family':'serif', 'serif':['Palatino']})
#rc('text', usetex=True)


#mag_lim = pm.mag_limit
band_colors = {'U':'Violet', 'B':'Blue', 'V':'Green', 'R':'Red', 'I':'Orange', 'J':'Pink',\
               'H':'Brown', 'K':'Gray', 'u': 'Violet', 'g':'Green', 'r': 'Orange', 'i':'Red', 'z': 'Brown'}

markerList = ['o','s','p','*','^', 'v', 'D','x', '<', '>', 8, '+','|']

nights_per_minute = 20   #nr. nights animated in one minute
ms_per_night = 60*1000 / nights_per_minute  #nr. of miliseconds that each night lasts in the animations
one_night = 24*60*60. #nr. of seconds in a night
Bigjumpfactor = 3.    #The factor by which large jumps in frame rate should be divided to keep the movie length reasonable


myCC = ColorConverter()

class AnData:
    def __init__(self, datafile, bandList):
        global markerList
        f = open(datafile, 'r')
        typesList = f.readline().split()
        hdr = f.readline().split()
        f.close()
        Data = np.genfromtxt( datafile, skip_header =2 , dtype = [int, int, float, float, float, float, float, ('string','S5')])
        self.obs_bands = np.array([row[7] for row in Data])
        data = np.array( [[row[0], row[1], row[2], row[3], row[4], row[5]] 
                          for row in Data] )
        if np.array(data).shape == (6,):
            raise Exception("Too little data to make an animation!")
        else:
            self.mag_norm = np.min( data[:][5:])
        for c in bandList:
            if c not in self.obs_bands:
                raise Exception( "The data does not contain any observations in the %s filter band. Please enter the correct colors in the color option: --colors u g r i z" % (c) )
        # set the window for the observation
        tm, ras, decs = np.unique(data[:,2]), data[:,3], data[:,4]
        self.ra_min = np.min(ras)
        self.ra_max = np.max(ras)
        self.dec_min = np.min(decs)
        self.dec_max = np.max(decs)          
        self.nBands = len(bandList)
        #col_idstypes = np.array([0,1,3,4])
        #col_idsbands = np.linspace( 5 , 5 + self.nBands - 1 , self.nBands ).astype(int) 
        col_ids = np.array([0,1,3,4,5])#np.append(col_idstypes, col_idsbands)
        self.timeFrames = []
        self.tF_bands   = []
        for t in tm:
            self.timeFrames.append(data[data[:,2] == t][:,col_ids])
            obsbands = self.obs_bands[data[:,2] == t]
            self.tF_bands.append(obsbands)
        #self.timeFrames = [ data[data[:,2] == t][:,col_ids] for t in tm ]
        self.times = tm
        dt = (tm[1] - tm[0])/8.64e4
        self.n_frames = len(self.times)
        nBands = len(bandList)
        # figure out the marker stuff
        trtypes = [int(ty) for ty  in typesList[::2]]
        trNames = typesList[1::2]
        if len(trtypes) > len(markerList):
            markerList += markerList[0: len(trtypes)-len(markerList)]
        self.nTrTypes = len(trtypes)
        self.trNames = trNames
        # add information about rates of change
        for i in range( self.n_frames ):
            tF_now = self.timeFrames[i]
            if i < self.n_frames-1:
                tF_nxt = self.timeFrames[i+1]
            else:
                tF_nxt = self.timeFrames[i-1]#[:,:-2]
                #print "hi"
                #tF_nxt = self.timeFrames[i-1][:,:-nBands]
            dmdtDat = np.zeros( np.shape(tF_now[:,4]) )
            for j, iD in enumerate(tF_now[:,0]):
                #print i,tF_nxt[tF_nxt[:,0]==iD][0]#.shape
                mags_now = tF_now[tF_now[:,0]==iD][0,4]
                if iD in tF_nxt[:,0]:
                    mags_nxt = tF_nxt[tF_nxt[:,0]==iD][0,4]
                else: mags_nxt = mags_now
                dMdt = np.abs( mags_now - mags_nxt )
                dMdt /= dt
                dmdtDat[j] = dMdt
            self.timeFrames[i] = np.column_stack( (self.timeFrames[i], dmdtDat ) )
        self.timeFrames = np.array(self.timeFrames)
        np.save('timeframes.npy', np.array(self.timeFrames))
        self.dmdt_max = max( np.max( tF[:,-nBands:] ) for tF in self.timeFrames )

        
        
class Animated:
    def __init__(self, Opts, pm):
        self.Opts = Opts
        bandlist = self.Opts.colors
        self.mag_lim = pm['mag_limit']
        if len(bandlist) <= 3:
            n_row, n_col = 1, len(bandlist)
        else:
            n_col = 2 
            n_row = (len(bandlist)+1)/2
        fig_x, fig_y = n_col*4.5 , n_row*4.5 + 1 

        self.fig, self.axs = plt.subplots( n_row, n_col, figsize = (fig_x, fig_y), sharey=True )
        self.bot = 1.0/fig_y
        self.fig.subplots_adjust( bottom = self.bot )
        if len(bandlist) > 1:
            self.axs = self.axs.flatten()
        else:
            self.axs = [self.axs]
        self.data = AnData(Opts.file, bandlist)
        for i, bnd in enumerate(bandlist):
            self.axs[i].set_xlim(self.data.ra_min, self.data.ra_max)
            self.axs[i].set_ylim(self.data.dec_min, self.data.dec_max)
            self.axs[i].set_title( bnd )
            self.axs[i].set_xlabel( 'RA' )
            self.axs[i].set_ylabel( 'DEC' )
        for j in range( self.data.nTrTypes):
            self.axs[0].scatter([],[],color='k',alpha = 0.5, s = 60,  marker = markerList[j], label=self.data.trNames[j] )
        self.axs[0].legend( frameon=True, loc=3, bbox_to_anchor = [0.1,self.bot],\
                                handletextpad=0.5, scatterpoints=1, bbox_transform=self.fig.transFigure )
        self.fig.set_tight_layout(True)
        self.bands = bandlist
        self.nBands = len(bandlist)
        self.size_min = 5.0
        size_max = 200.0
        self.scatters = []
        self.timetext = ''
        self.size_grad = {}
        for band in self.bands:
            self.size_grad[band] = (size_max - self.size_min)/(self.data.mag_norm -self.mag_lim[band])
        self.alpha_grad = 0.85/self.data.dmdt_max
        #self.plot_init()
        #Frames = np.arange(self.data.n_frames)
        self.init_framelength()
        Frames = [0]
        for i in range(self.data.n_frames - 1):
            for j in range(self.framelength[i]):
                Frames.append(i+1)
        self.anim = animation.FuncAnimation(  self.fig, self.update, Frames,\
                                              init_func=self.plot_init, repeat=False)
        
    def init_framelength( self ):
        """
        This function calculates the different time intervals between each
         obsrvation. This time interval is used to calculate the number of 
         frames each time stamp should encompass. This is done to make the 
         animation more correct with time. For example: The time between the 
         first two observations may take 5 minutes, while the time between the 
         next two is 23 hours. It would look strange if the time between each 
         frame were the same.
        
        self.framelength gives the number of frames each observation should 
         last
        self.FPN is the number of frames a night lasts: frames per night.
        """
        self.timeintervals = np.zeros(self.data.n_frames - 1)
        for i in range(len(self.data.times) - 1):
            self.timeintervals[i] = self.data.times[i+1] - self.data.times[i]
        #self.framelength = np.ones(len(self.data.timeintervals))
        minframelength = min(self.timeintervals)
        self.framelength = (self.timeintervals / minframelength).astype(int)
        if self.Opts.samefr:
            if len(self.framelength) > 4:
                uniq_jumps =  np.sort(np.unique(self.framelength))
                if len(uniq_jumps) > 2:
                    largestjumps = uniq_jumps[:-3]
                    if largestjumps[1] > 2 * Bigjumpfactor * largestjumps[0]:
                        Bool_largest = np.array(self.framelength in largestjumps[:2])
                        self.framelength[Bool_largest] = self.framelength[Bool_largest] /Bigjumpfactor
                    if largestjumps[2] > 2 * Bigjumpfactor * Bigjumpfactor * largestjumps[1]:
                        Bool_largest = np.array(self.framelength == largestjumps[2])
                        self.framelength[Bool_largest] = self.framelength[Bool_largest] / Bigjumpfactor
                else:
                    largestjumps = uniq_jumps
                    if largestjumps[1] > 2 * Bigjumpfactor * largestjumps[0]:
                        Bool_largest = np.array(self.framelength == largestjumps[1])
                        self.framelength[Bool_largest] = self.framelength[Bool_largest] / Bigjumpfactor
        if sum(self.timeintervals) > one_night:   #The animation comprises more than one night:
            #fullframelength = self.framelength * minframelength
            timecumsum      = np.cumsum(self.timeintervals)
            end_night_one   = np.argwhere(timecumsum)[0]
            self.FPN        = np.cumsum(self.framelength)[end_night_one]   #The approx nr. of frames for one night.
        else:
            fractionofnight = np.sum(self.timeintervals) / one_night
            self.FPN        = np.sum(self.framelength) / fractionofnight

    def plot_init( self ):
        self.init_bands = []
        all_bands_initiated = False
        i = 0
        not_initiated_bands = self.bands
        colornr = 0
        #for i,band in enumerate( self.bands ):
        #print "hello", not_initiated_bands
        while not all_bands_initiated:
            #print i, len(self.data.tF_bands)
            #print self.data.tF_bands[i], not_initiated_bands
            uniqbands = np.unique(self.data.tF_bands[i] )
            innotinitb = np.array([b in not_initiated_bands for b in uniqbands] )
            if any(innotinitb):
                for band in uniqbands[innotinitb]:
                    #print np.shape(self.data.timeFrames[i,:]), np.shape(self.data.tF_bands), np.shape(self.data.times[i])
                    init_dat, init_t = self.data.timeFrames[i], self.data.times[i]
                    Boolean1 = np.array(init_dat[:,4] <= self.mag_lim[band])
                    Boolean2 = np.array(self.data.tF_bands[i] == band)
                    Boolean = Boolean1 * Boolean2
                    this_type, this_ra, this_dec, this_mag, this_dmdt = init_dat[Boolean].T[[1,2,3,4,5]]
                    this_type = np.array( this_type, dtype = int )
                    # alpha => rate of change
                    # size => brightness
                    sizes = self.size_min + (this_mag - self.mag_lim[band])*self.size_grad[band]
                    rgba_cols = np.zeros( [len(this_ra), 4] )
                    rgba_cols[:,0:3] = myCC.to_rgb(band_colors[band])
                    rgba_cols[:,-1] = 0.15 + self.alpha_grad * this_dmdt
                    for k in range( self.data.nTrTypes ):
                        ids = np.where( this_type == k )[0]
                        if k == 1:
                            self.idsblablab = ids
                        sc = self.axs[colornr].scatter( this_ra[ids], this_dec[ids], s = sizes[ids], c = rgba_cols[ids],\
                                                       marker=markerList[k], figure=self.fig)  
                        #sc = ax1[0].scatter( this_ra[ids], this_dec[ids])
                        #print (this_ra[ids], this_dec[ids])
                        sc.set_edgecolor('none')
                        self.scatters.append(sc)
                    #plt.show()

                    self.init_bands.append(band)
                    not_initiated_bands = [x for x in not_initiated_bands if x != band]
                    colornr += 1
                i += 1
            else:
                i += 1
            if not_initiated_bands == [] or i > len(self.data.tF_bands):
                all_bands_initiated = True
        self.init_bands = np.array(self.init_bands)
        self.time_text = self.axs[0].text( 0.0, 1.03, "t = %.2f days" % (init_t/8.64e4), transform=self.axs[0].transAxes)
        #self.anim.event_source.interval =  self.data.timeintervals[0] / one_night * ms_per_night   #Sets the interval between two frames. 
        return self.scatters, self.timetext   
    def update( self, frame_no ):
        #print "hi", self.anim.event_source.interval
        dat, t = self.data.timeFrames[frame_no], self.data.times[frame_no]
        uniqbands = np.unique(self.data.tF_bands[frame_no])
        for band in uniqbands:
            if band in self.bands:
                i = np.where(self.init_bands == band)[0][0] 
                for blabla in range(1): #i,band in enumerate( self.bands ):
                    this_type, this_ra, this_dec, this_mag, this_dmdt = dat[dat[:,4] <= self.mag_lim[band]].T[[1,2,3,4,5]]
                    this_type = np.array( this_type, dtype=int )
                    sizes = self.size_min + (this_mag - self.mag_lim[band])*self.size_grad[band]
                    rgba_cols = np.zeros( [len(this_ra), 4] )
                    rgba_cols[:,0:3] = myCC.to_rgb(band_colors[self.bands[i]])
                    rgba_cols[:,-1] = 0.15 + self.alpha_grad * this_dmdt
                    krng = self.data.nTrTypes
                    for k in range( krng ):
                        ids = np.where( this_type == k )[0]
                        self.scatters[i*krng+k].set_offsets( np.column_stack( (this_ra[ids], this_dec[ids]) ))
                        self.scatters[i*krng+k].set_facecolor( rgba_cols[ids])
                        self.scatters[i*krng+k].set_edgecolor('none')
                        self.scatters[i*krng+k].set_sizes( sizes[ids] )
        self.time_text.set_text( 't = %.2f days' % (t/8.64e4))
        #self.anim.event_source.interval =  (self.data.timeintervals[frame_no-1] 
        #                                    / one_night * ms_per_night )  #Sets the interval between two frames. 
        return self.scatters,self.timetext    
    
class pmOpts:
    """
    Objects of this class should only be generated if colorcolor.py is not run
     directly.
    Objects of this class contain the options that would otherwise be generated 
     by getOpts().
    """
    def __init__(self, bands, Opts_TSS, pm):
        self.params = pm['filename']
        self.file   = Opts_TSS.output
        self.output = pm['Ani_outfile']
        self.samefr = pm['Ani_samefr']
        self.colors = self.selectcolors( bands, Opts_TSS, pm )
        
    def selectcolors(self, bands, Opts_TSS, pm ):
        correct_colors = np.all([C in bands for C in pm['Ani_bands']])
        if correct_colors:
            return pm['Ani_bands']
        else:
            print "Warning: the wrong colors were requested. Please check that Ani_bands in params.py contains only colors that were observed."
            print "Therefore, an animation is being created based on the colors that were observed."
            return bands


def _AnimateSky(Opts, pm):
    print "beginning animation of colors ", Opts.colors
    an = Animated( Opts, pm )
    Writer = animation.writers['ffmpeg']
    FPS = int(an.FPN * nights_per_minute)
    writer = Writer(fps = FPS, metadata=dict(artist='me'), bitrate=1800)
    an.anim.save(Opts.output, writer = writer) 
    print "Animation saved to ", Opts.output
    
def AnimateSky(bands, Opts_TSS, pm):
    """
    A function that calls _AnimateSky with parameters that would usually be
     in Opts taken from pm.
    This function should be called by external python files.
    """
    Opts = pmOpts(bands, Opts_TSS, pm)
    _AnimateSky(Opts, pm)

def getOpts():
    parser = argparse.ArgumentParser(description='Animate an output of TSS')
    parser.add_argument("-c", "--colors", nargs = '*', default = ['g'], help="The colors/passbands to animate. Can be multiple. Parse them as 'u g r'")
    parser.add_argument("-p", "--params", default = 'params.py', help="Define a file with observation parameters. default:params.py")
    parser.add_argument("-f", "--file", help="The file (that was output by TSS) to animate. default: outfile in the params file")
    parser.add_argument("-r", "--samefr", action='store_true', help ="Do not change the framerate for two epochs that are far apart")
    parser.add_argument("-o", "--output", default = 'Animation.mp4', help ="The output file. default: Animation.mp4")
    parser.set_defaults(feature=True)
    args = parser.parse_args()
    return args

def printOpts(Opts):
    print "Running Animation.py with options:"
    if Opts.colors:
        print "[-c] [--colors]  ", Opts.colors
    if Opts.params:
        print "[-p] [--params]  ", Opts.params
    if Opts.file:
        print "[-f] [--file]    ", Opts.file
    if Opts.samefr:
        print "[-r] [--samefr] ", Opts.samefr
    if Opts.output:
        print "[-o] [--output] ", Opts.output


if __name__ == "__main__": 
    """
    Run as: python Animation.py [Arguments]
    Optional arguments:
    [-c] [--colors]    The colors/passbands to animate. Can be multiple. 
    [-p] [--params]    Params file to use. default:params.py
    [-f] [--file]      The file (that was output by TSS) to animate. default: outfile in the params file
    [-r] [--samefr]    Do not change the framerate for two epochs that are far apart
    [-o] [--output]    The output file. default: Animation.mp4
    These arguments will override any arguments in params.py
    """
    Opts = getOpts()
    printOpts(Opts)
    
    pm = {}
    execfile(Opts.params, pm)
    pm['filename'] = Opts.params
    if not Opts.file:
        Opts.file = pm['outfile']
        
    #AnimateSky(['u', 'g'], 'long.dat')
    _AnimateSky(Opts, pm)

