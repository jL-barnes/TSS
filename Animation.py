import matplotlib
#matplotlib.use( 'macosx' )
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import ColorConverter
from matplotlib import rcParams
from matplotlib import rc
from astLib import astCoords
import params as pm

rcParams.update({'figure.autolayout':True})
#rc('font',**{'family':'serif', 'serif':['Palatino']})
#rc('text', usetex=True)


mag_lim = pm.mag_limit
band_colors = {'U':'Violet', 'B':'Blue', 'V':'Green', 'R':'Red', 'I':'Orange', 'J':'Pink',\
               'H':'Brown', 'K':'Gray', 'u': 'Violet', 'g':'Green', 'r': 'Orange', 'i':'Red', 'z': 'Brown'}

markerList = ['o','s','p','*','^', 'v', 'D','x', '<', '>', 8, '+','|']


myCC = ColorConverter()

class AnData:
    def __init__(self, datafile, bandList):
        global markerList
        f = open(datafile, 'r')
        typesList = f.readline().split()
        hdr = f.readline().split()
        f.close()
        Data = np.genfromtxt( datafile, skip_header =2 , dtype = [int, int, float, float, float, float, ('string','S5')])
        self.obs_bands = np.array([row[6] for row in Data])
        data = np.array( [[row[0], row[1], row[2], row[3], row[4], row[5]] 
                          for row in Data] )
        if np.array(data).shape == (6,):
            raise Exception("Too little data to make an animation!")
        else:
            self.mag_norm = np.min( data[:][5:])
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
            self.tF_bands.append(self.obs_bands[data[:,2] == t][0])
        
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
        np.save('timeframes.npy', np.array(self.timeFrames))
        self.dmdt_max = max( np.max( tF[:,-nBands:] ) for tF in self.timeFrames )
        
        
fig1, ax1 = plt.subplots(4,1)
class Animated:
    def __init__(self, dataFile, bandlist):
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
        self.data = AnData(dataFile, bandlist)
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
            self.size_grad[band] = (size_max - self.size_min)/(self.data.mag_norm -mag_lim[band])
        self.alpha_grad = 0.85/self.data.dmdt_max
        #self.plot_init()
        Frames = np.arange(20)#np.arange(self.data.n_frames)
        self.anim = animation.FuncAnimation(  self.fig, self.update, Frames,\
                                              init_func=self.plot_init, repeat=False)
        #self.anim = animation.FuncAnimation( fig1, self.update, Frames, \
        #                                        init_func= self.plot_init, repeat=False)
    def plot_init( self ):
        self.init_bands = []
        for i,band in enumerate( self.bands ):
            init_dat, init_t = self.data.timeFrames[i], self.data.times[i]
            this_type, this_ra, this_dec, this_mag, this_dmdt = init_dat[init_dat[:,4] <= mag_lim].T[[1,2,3,4,5]]
            this_pos = np.column_stack( (this_ra, this_dec) )
            this_type = np.array( this_type, dtype = int )
            # alpha => rate of change
            # size => brightness
            sizes = self.size_min + (this_mag-mag_lim[band])*self.size_grad[band]
            rgba_cols = np.zeros( [len(this_ra), 4] )
            rgba_cols[:,0:3] = myCC.to_rgb(band_colors[band])
            rgba_cols[:,-1] = 0.15 + self.alpha_grad * this_dmdt
            for k in range( self.data.nTrTypes ):
                ids = np.where( this_type == k )[0]
                sc = self.axs[i].scatter( this_ra[ids], this_dec[ids], s = sizes[ids], c = rgba_cols[ids],\
                                              marker=markerList[k], figure=self.fig)  
                #sc = ax1[0].scatter( this_ra[ids], this_dec[ids])
                #print (this_ra[ids], this_dec[ids])
                #sc = ax1[i].scatter(this_ra[ids], this_dec[ids], s=sizes[ids], c = rgba_cols[ids], marker=markerList[k], figure=fig1)
                plt.show()

                sc.set_edgecolor('none')
                self.scatters.append(sc)
            self.init_bands.append(band)
        self.init_bands = np.array(self.init_bands)
        self.time_text = self.axs[0].text( 0.0, 1.03, "t = %.2f days" % (init_t/8.64e4), transform=self.axs[0].transAxes)
        return self.scatters, self.timetext   
    def update( self, frame_no ):
        dat, t = self.data.timeFrames[frame_no], self.data.times[frame_no]
        band = self.data.tF_bands[frame_no]
        i = np.where(self.init_bands == band)[0][0] 
        for blabla in range(1): #i,band in enumerate( self.bands ):
            this_type, this_ra, this_dec, this_mag, this_dmdt = dat[dat[:,4] <= mag_lim[band]].T[[1,2,3,4,5]]
            this_pos = np.column_stack( (this_ra, this_dec) )
            this_type = np.array( this_type, dtype=int )
            sizes = self.size_min + (this_mag-mag_lim[band])*self.size_grad[band]
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
        return self.scatters,self.timetext     

def AnimateSky(bList, skyfile):
    #bList = [b for b in pm.showbands]
    print "beginning animation"
    an = Animated( skyfile, bList)
    Writer = animation.writers['ffmpeg']
    #FFwriter = animation.ImageMagickWriter(fps=5)
    writer = Writer(fps=5, metadata=dict(artist='me'), bitrate=1800)
    an.anim.save('Animation.mp4', writer = writer) 
    #an.anim.save('Animation.mp4', writer=FFwriter)
    #plt.show()

def AnimateFile( bList, skyfile ):
    an = Animated( skyfile, bList, fromfile = True )
    an.show()

