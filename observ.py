import numpy as np
import csv
import MilkyWay as MW
import cosmology as cm
import dust
import h5py
from astLib import astCoords
from utils import getFileLine

class cell:
    """ A cell object is a cell in a grid.

        idx: the id/number of the cell
        raMid: the center RA coordinate of the cell in degrees
        decMid: the center DEC coordinate of the cell in degrees
        DMid: the center luminosity distance coordinate of the cell in kpc
        hD: size of the cell in the distance direction (differs between 
         galactic and non-galactic) in kpc
        hDEC: size of the cell in the DEC direction in degrees
        vol: The volume of this cell (t.b. changed later)
        z: The redshift of the center of this cell
        Grid: The parent grid object.
        cellnr_from_left: There are self.Grid.N_RA cells in a row, this 
         is the cfl'th
    """
    def __init__( self, idx, cfl, raMid, decMid, distMid, hdist, hDEC, z, Grid):
        self.idx    = idx
        self.raMid  = raMid
        self.decMid = decMid
        self.DMid   = distMid
        self.hD     = hdist
        self.hDEC   = hDEC
        self.vol    = 0.0
        self.z      = z
        self.Grid   = Grid
        self.cellnr_from_left = cfl # 

    def setDensity( self, Hthin = 0.3):
        """ Set the Galactic transient densities of this cell.

            This sets rho, the transient density of the cell.
        """
        self.rho = MW.get_MW_dens([self.raMid, self.decMid, self.DMid], Hthin)
    def setVolume( self ):
        """ Calculate the cell volume.

            This is done by subtracting the volume of a sphere with 
             r=Dmax from a sphere with r=Dmin. Then we multiply this by 
             the fraction that the FoV is of the full sky and the 
             fraction that the cell is of the FoV.
            cosm is the cosmological factor that converts a volume to a comoving volume
        """
        if self.DMid > self.Grid.Dmax_MK:
            z1 = self.Grid.Cosmo.get_redshift(  (self.DMid - 0.5*self.hD) / 1.e3 )	#convert to Mpc
            z2 = self.Grid.Cosmo.get_redshift(  (self.DMid + 0.5*self.hD) / 1.e3 )	#convert to Mpc
        else:
            z1, z2 = 0,0
        cosm1 = np.power( 1 + z1, 3.)
        cosm2 = np.power( 1 + z2, 3.)
        Sphere1 = 4./3. * np.pi * np.power( self.DMid - 0.5*self.hD, 3.0 ) / cosm1
        Sphere2 = 4./3. * np.pi * np.power( self.DMid + 0.5*self.hD, 3.0 ) / cosm2
        Total_vol = Sphere2 - Sphere1
        Fullsky = 129600 / np.pi
        cell_fraction = 1. / (self.Grid.N_RA * self.Grid.N_DEC)
        self.vol = (self.Grid.FoV / Fullsky) * cell_fraction * Total_vol
    def sampleCoordsold( self ):
        """ Sample a random RA, DEC and Distance coordinate within the 
             cell volume.

            ---------
            returns: a list with the sampled RA, DEC and Distance 
                     coordinates. The distance coordinate is in kpc
        """
        z1 = np.random.random()
        z2 = np.random.random()
        z3 = np.random.random()
        randRA   = self.raMid   + (z1 - 0.5)* self.hRA
        randDEC  = self.decMid  + (z2 - 0.5)* self.hDEC
        randDist = self.DMid + (z3 - 0.5)* self.hD
        return [randRA, randDEC, randDist]   
    def sampleCoords( self ):
        """ Sample a random RA, DEC and Distance coordinate within the 
             cell volume.

            The cell is NOT rectangular. 
            Therefore we first sample a declination.
            At that declination we calculate the width of the cell in 
             RA, and sample the RA coordinate.
            We also sample a random distance coordinate within the cell
            Ultimately we need to know whether this particular 
             coordinate overlaps with other cells. 
            The transients in the overlapping part of the cell need to 
             be distributed accordingly over all overlapping grids.
            We therefore inject this transient in THIS cell/grid with a 
             probability of 1/(nr. of overlapping cells + 1)

            ---------
            returns: a list with the sampled RA, DEC and Distance 
                      coordinates. The distance coordinate is in kpc.
                     If the transient is rejected, an empty list is
                      returned.
        """

        z1 = np.random.random()
        z2 = np.random.random()
        z3 = np.random.random()
        randDEC  = self.decMid  + (z1 - 0.5)* self.hDEC
        width_at_DEC = self.Grid.obRun.aperture_RA / np.cos(randDEC * (np.pi/180.)) 
        self.hRA = width_at_DEC / self.Grid.N_RA
        startRA = ( self.Grid.RA_c - 0.5 * width_at_DEC + 
                    (self.cellnr_from_left + 0.5) * self.hRA )
        randRA   = startRA + (z2 - 0.5)* self.hRA
        randDist = self.DMid + (z3 - 0.5)* self.hD
        
        N_overlapping_grids = 0
        for i,grid in enumerate(self.Grid.obRun.SkyGrids):
            if i != self.Grid.gridnr:
                if grid.WithinFrame(randRA, randDEC):
                    N_overlapping_grids +=1
        rand_nr = np.random.random_sample()
        if rand_nr < 1. / (N_overlapping_grids + 1):
            return [randRA, randDEC, randDist]
        else:
            return []


class grid:
    """ A grid instance is the grid-ified version of a frame. 

        A field-of-view (or frame) is split into a grid of N_RA cells
         in the RA direction, N_DEC cells into the DEC direction and 
         N_DMW + N_xgal cells in the distance direction.
        A grid instance has information on the coordinates of the grid
         and the accompanying grid of dust extinction.
    """
    def __init__(self, pm, Dmax_xgal, gridnr, obRun):
        """
            pm: the parameter file
            Dmax_xgal: the maximum distance up to which transients are 
                       visible in kpc.
            gridnr: the frame number/identifying number of the grid
            obRun: the parent Observation instance
        """
        self.gridnr      = gridnr
        self.obRun       = obRun
        self.SetUpFramecoords()
        self.Dmax_MK     = MW.DMax_MW_kpc
        self.Dmax_xgal   = Dmax_xgal
        self.N_RA        = pm['n_RA']
        self.N_DEC       = pm['n_DEC']
        self.N_cellsRADEC= self.N_RA * self.N_DEC
        self.FoV         = obRun.aperture_DEC * obRun.aperture_RA
        self.N_DMW       = pm['nCells_D']
        self.N_xgal      = pm['nCells_xgal']
        self.h_DMW       = self.Dmax_MK / float(self.N_DMW)
        self.h_xgal      = (Dmax_xgal - self.Dmax_MK) / float(self.N_xgal)
        self.Gal_trans	 = obRun.Gal_trans    #Whether there are Gal transients
        self.Xgal_trans  = obRun.Xgal_trans	#Whether there are Xgal transients
        self.bands       = self.obRun.bands
        if self.Xgal_trans:
            self.Cosmo = cm.Cosmology(Dmax_xgal)
        self.cellGrid    = self.makeCellGrid( self.h_DMW, self.h_xgal )
        self.choose_dust_grid()

    def SetUpFramecoords(self):
        """ Set up the boundary and center coordinates of the grid.

            Keep in mind that in RA/DEC space the frame is rectangular.
            As such the "left" RA_low coordinate is not equal to the 
             "left" RA_high coordinate.
        """
        centercoords     = self.obRun.obFrameCoords[self.gridnr]
        RA_c             = centercoords[0]
        DEC_c            = centercoords[1]
        self.RA_c        = astCoords.hms2decimal(RA_c,":")
        self.DEC_c       = astCoords.dms2decimal(DEC_c,":")
        DEC_lo           = self.DEC_c - 0.5 * self.obRun.aperture_DEC
        DEC_hi           = self.DEC_c + 0.5 * self.obRun.aperture_DEC
        if DEC_hi < DEC_lo: 
            DEC_lo, DEC_hi = DEC_hi, DEC_lo        
        RA_lo_l          = self.RA_c - (0.5 * self.obRun.aperture_RA /
                            np.cos(DEC_lo * (np.pi/180.)) )
        RA_lo_r          = self.RA_c + (0.5 * self.obRun.aperture_RA /
                            np.cos(DEC_lo * (np.pi/180.)) )
        RA_hi_l          = self.RA_c - (0.5 * self.obRun.aperture_RA /
                            np.cos(DEC_hi * (np.pi/180.)) )
        RA_hi_r          = self.RA_c + (0.5 * self.obRun.aperture_RA /
                            np.cos(DEC_hi * (np.pi/180.)) )
        if RA_lo_r < RA_lo_l: 
            RA_lo_r, RA_lo_l = RA_lo_l, RA_lo_r  
        if RA_hi_r < RA_hi_l: 
            RA_hi_r, RA_hi_l = RA_hi_l, RA_hi_r  
        if RA_hi_l < RA_lo_l: 
            RA_lo_r, RA_hi_r, RA_lo_l, RA_hi_l = \
                RA_hi_r, RA_lo_r, RA_hi_l, RA_lo_l
        self.DEC_lo  = DEC_lo
        self.DEC_hi  = DEC_hi
        self.RA_lo_l = RA_lo_l
        self.RA_lo_r = RA_lo_r
        self.RA_hi_l = RA_hi_l
        self.RA_hi_r = RA_hi_r
        
    def makeCellGrid( self, h_DMW, h_xgal ):
        """ Creates the cells inside a grid

            In this function we first loop over all cells in the
             distance direction.
            For each of thos distances we determine the redshift.
            And for each distance we generate the cells in the RA/DEC
             direction and calculate the center RA and DEC coordinates
             for each cell
            hDec is the width of a cell in the DEC-direction. This is 
             constant throughout the frame as the frame is aligned with
             this direction.

            -------
            returns: a list of the created cells
        """
        hDec             = (self.DEC_hi-self.DEC_lo) / float(self.N_DEC)
        cells = []
        for k in range( self.N_DMW + self.N_xgal ):
            if k < self.N_DMW: 
                mid_D = ( float(k) + 0.5) * h_DMW
                dh = self.h_DMW
                z = 0		#Redshift
            else:          
                if self.Xgal_trans:
                    mid_D = self.Dmax_MK + ( float(k) - self.N_DMW + 0.5 ) * h_xgal
                    dh    = h_xgal
                    z     = self.Cosmo.get_redshift( mid_D / 1.e3 )	#convert to Mpc
            for j in range( self.N_DEC ):
                mid_dec      = self.DEC_lo + (float(j) + 0.5) * hDec
                width_at_DEC = ( self.obRun.aperture_RA / 
                                np.cos(mid_dec * (np.pi/180.)) )
                Cwidth       = width_at_DEC / self.N_RA
                for i in range( self.N_RA ):
                    mid_ra  = ( self.RA_c - 0.5 * width_at_DEC + 
                               (float(i) + 0.5) * Cwidth )
                    idx     = i + self.N_RA*j + self.N_RA*self.N_DEC*k
                    newcell = cell( idx, i, mid_ra, mid_dec, mid_D, dh, hDec, z , self )
                    cells.append( newcell )
        return cells   
    def setCellVolumes( self ):
        """ Sets the volume and stellar density of every cell in this Grid
        """
        for cell in self.cellGrid:
            cell.setVolume()
            cell.setDensity()
    def choose_dust_grid( self ):
        """ Choose the right dust grid.

            If the --nodust option was chosen when running main.py, 
             we exclude dust extinction from this calculation.

            We have three different dust maps:
             -The Schlegel, 2D dust map that covers the entire sky
             -The Green map, a 3D dust map of the Milky Way that only
               covers the Northern hemisphere
             -The Schultheis 3D dust map that covers the field towards
               the Galactic bulge.

            For extragalactic transients (and galactic transients with 
             D>60kpc) we choose the Schlegel map.
            If the FOV is included in the Schultheis map, we choose
             that map, if it is included in the Green map, we pick the
             Green 3D dust map.
            For galactic transients (if they are used) it is possible 
             that the FOV is included in the Schultheis as well as the 
             Green map. Then Schultheis gets priority
            If the FOV is neither included in the Schultheis or the
             Green map, we use the 2D dust map from Schlegel. This does
             affect the accuracy, but as of V1.0 there were no 3D dust 
             maps available for the Southern Hemisphere.
        """
        RA_lo = min(self.RA_lo_l, self.RA_lo_r, self.RA_hi_l, self.RA_hi_r)
        RA_hi = max(self.RA_lo_l, self.RA_lo_r, self.RA_hi_l, self.RA_hi_r)
        coordboundaries = (RA_lo, RA_hi, self.DEC_lo, self.DEC_hi)
        offline = self.obRun.Opts.offline
        
        if not self.obRun.nodust:
            self.Xtr_dust = dust.Schlegel_Extinction(self.bands, self.obRun.colorScheme, offline, *coordboundaries)
            if self.Gal_trans:  #There are galactic transients
                if dust.InSchultheisBoundary(*coordboundaries):
                    self.Gal_dust = dust.Schultheis_Extinction(self.Xtr_dust, self.bands, self.obRun.colorScheme, offline, *coordboundaries)
                elif dust.InGreenBoundary(*coordboundaries):
                    self.Gal_dust = dust.Green_Extinction(self.Xtr_dust, self.bands, self.obRun.colorScheme, offline, *coordboundaries)
                else:
                    self.Gal_dust = self.Xtr_dust
        else:
            print("Excluding dust extinction in the simulation")
            self.Xtr_dust = dust.No_dust(self.bands)
            self.Gal_dust = dust.No_dust(self.bands)
    def WithinFrame( self, RA, DEC ):
        """ Check whether the object at coordinates (RA, DEC) is within 
             the frame of this grid instance.

            --------
            returns: a boolean whether the object is in the frame
        """
        if DEC < self.DEC_lo or DEC > self.DEC_hi:
            return False
        width_at_DEC = self.obRun.aperture_RA / np.cos(DEC * (np.pi/180.)) 
        RA_l = self.RA_c - (0.5 * width_at_DEC )
        RA_r = self.RA_c + (0.5 * width_at_DEC )
        if RA < RA_l or RA > RA_r:
            return False
        return True
            

class timeFrame:
    """ A timeFrame instance contains all the information on the
         observation schedule.

        It basically reads out the Obstimes file.
    """
    def __init__( self, File ):
        self.JD        = []
        self.Band      = []
        self.RAc       = []
        self.DECc      = []
        self.read_Obsfile(File)
        self.t_start = self.JD[0]
        self.t_end = self.JD[-1]
        self.Frames    = self.SetUpFrames()
    def read_Obsfile(self, File):
        """ Read the Obstimes file

            This file contains a list of observations. Each observation
             may have its own coordinates and color filter.

            The data is then stored into four lists:
            JD: the date/moment of observation in Julian Date
            Band: the color filter of the observation
            RAc: the center RA coordinate of the observation
            DECc: the center DEC coordinate of the observation
        """
        with open(File) as Obsfile:
            data = csv.reader(Obsfile, delimiter='\t')
            for i,row in enumerate(data):
                if i not in [0,1]:
                    self.JD.append(float(row[0]) * 8.64e4)
                    self.Band.append(row[1])
                    self.RAc.append(row[2])
                    self.DECc.append(row[3])
    def SetUpFrames( self ):
        """ Sets up a few frame parameters
  
            FrameCoords: a list of unique RA and DEC coordinates of this
             observation schedule. These are called frames.
            Nr_Frames: The number of unique frames
            Framenrs: A list of frame numbers
            Frames: A list where each observation has received a number
             that corresponds to the observation's framenr

            -------
            returns: Frames
            
        """
        listcoords       = list(zip(self.RAc, self.DECc))
        self.FrameCoords = list(set(listcoords))
        self.Nr_Frames   = len(self.FrameCoords)
        self.Framenrs    = np.arange(self.Nr_Frames)
        Fcoords          = np.array([''.join(x)  for x in self.FrameCoords])
        Frames           = np.zeros(len(self.RAc))
        for i in range(len(self.RAc)):
            Coord = ''.join(listcoords[i])
            Frames[i] = np.argwhere(Fcoords == Coord).ravel()[0]
        return Frames
    def get_Dates( self ):
        """ Creates a list of the observation dates for each frame

            -------
            returns: this list
        """
        obTimes = np.array(self.JD) - self.t_start
        Dates = []
        for frame in self.Framenrs:
            Dates.append(obTimes[self.Frames == frame])
        return Dates
    def get_ObBands( self ):
        """ Creates a list of the color filters for each frame

            -------
            returns: this list
        """
        Bands = []
        _Band = np.array(self.Band)
        for frame in self.Framenrs:
            Bands.append(_Band[self.Frames == frame])
        return Bands
    def get_Bands( self ):
        """ Return an array of all the unique color filters used
        """
        return np.unique(self.Band)
    def get_Frames( self ):
        """ Return the Frames list
        """
        return self.Frames

class observation:
    """ An observation instance contains all the information on this
         particular observation sequence.
    """
    def __init__(self, pm, Opts): 
        """
            pm: The parameters file
            Opts: The options whith which the program as run

            At the end for each frame a corresponding grid (instance) 
             is generated.
        """
        self.Opts           = Opts
        self.pm             = pm
        self.transientFile  = pm['transientFile']
        self.MdwarfFile     = pm['MDwarfFile']
        self.mag_resolution = pm['mag_resolution']
        self.aperture_DEC   = pm['aperture_DEC']
        self.aperture_RA    = pm['aperture_RA']
        self.maxLIGODist    = pm['maxLIGODist']
        self.Gal_trans      = False
        self.Xgal_trans     = False
        self.set_TransientsList()
        self.TimeParams = timeFrame( self.Opts.file )
        self.nodust         = Opts.nodust
        self.set_Colors( pm['color_system'] )
        self.set_Mag_lim(pm['mag_limit'])
        # determine max distance (extra-galactic)
        Dmax_xgal = self.Get_Max_observing_D()
        # set the properties of the sky window
        self.obTimes       = self.TimeParams.get_Dates()
        self.obBands       = self.TimeParams.get_ObBands()
        self.obFrames      = self.TimeParams.get_Frames()
        self.Nframes       = self.TimeParams.Nr_Frames
        self.obFrameCoords = self.TimeParams.FrameCoords   #Center coordinates of the frames
        self.SkyGrids      = []
        for i in range(self.Nframes):
            self.SkyGrids.append( grid( self.pm, Dmax_xgal, i, self ))
            self.SkyGrids[i].setCellVolumes()
    def set_Colors( self, col ):
        """ Set the color system of this observation

            The Bessell UBVRI system is always used as a base color
             system. On top of this another system can be chosen. 
            This can still be the UBVRI system. Usually the secondary
             system should be a ugriz system of which we need to
             transform the name to lowercase.
        """
        if col.upper() == 'UBVRI': colr = col.upper()
        elif col.lower() == 'sdss': colr = col.lower()
        elif col.lower() == 'ztf': colr = col.lower()
        elif col.lower() == 'blackgem': colr = col.lower()
        elif col.lower() == 'lsst': colr = col.lower()
        else: raise Exception(col, "is an invalid color system.")
        self.colorScheme = colr
        self.bands = self.TimeParams.get_Bands()
        self.nBands = len(self.bands)
    def set_Mag_lim(self, mag_limit):
        """ Sets the magnitude limit in each color band
        """
        self.threshold = mag_limit
        for c in self.bands:
            if c not in self.threshold:
                raise Exception("There is no magnitude limit defined for ", c, " in params.py")
        self.mag_lim        = [ self.threshold[ band ] for band in self.bands]
    def set_TransientsList( self ):
        """ Sets up a list of transients. These are obtained from the 
             parameters file.
            
            It also creates a legend where a transient is coupled to an 
             integer.
            It also decides whether there are any Galactic and any 
             extragalactic transients.
        """
        self.transientsList = []
        if self.pm['use_nova']:     self.transientsList.append( 'nova' )
        if self.pm['use_UGem']:     self.transientsList.append( 'UGem' )
        if self.pm['use_SUUMa']:    self.transientsList.append( 'SUUMa' )
        if self.pm['use_ZCam']:     self.transientsList.append( 'ZCam' )
        if self.pm['use_SNIa']:     self.transientsList.append( 'SNIa' )
        if self.pm['use_SNIb']:     self.transientsList.append( 'SNIb' )
        if self.pm['use_SNIc']:     self.transientsList.append( 'SNIc' )
        if self.pm['use_SNIIL']:    self.transientsList.append( 'SNIIL' )
        if self.pm['use_SNIIP']:    self.transientsList.append( 'SNIIP' )
        if self.pm['use_SNIInL']:   self.transientsList.append( 'SNIInL' )
        if self.pm['use_SNIInP']:   self.transientsList.append( 'SNIInP' )
        if self.pm['use_M3']:       self.transientsList.append( 'M3' )
        if self.pm['use_M3_5']:     self.transientsList.append( 'M3_5' )
        if self.pm['use_M4']:       self.transientsList.append( 'M4' )
        for transient in self.pm['Extra_transients']:
            self.transientsList.append( transient )
        self.Trlegend = {}
        for i,transient in enumerate(self.transientsList):
            self.Trlegend[transient] = i
            trData = getFileLine( self.transientFile, transient )
            Galtype = int( trData[0] )
            if Galtype in [0,1,2]:
                self.Gal_trans = True	#There are galactic transients
            if Galtype in [3,4]:
                self.Xgal_trans = True	#There are extragalactic transients
        if self.pm['use_kilonova']: self.Xgal_trans = True
    def Get_Max_observing_D( self ):
        """ Finds the maximum observing distance to which we need to 
             generate extragalactic transients.

            For each transient we first open the light curve files
             and determine their peak magnitude. If possible we also
             use Kcorrections to determine this.
            The brightest transient is updated with each transient we
             go through.

            -------
            returns: The maximum distance in kpc
        """
        brightest = 100.0 * np.ones( len( self.bands ) )
        for transient in self.transientsList:
            trData = getFileLine( self.transientFile, transient )
            if not int(trData[0]) in [3,4]:
                continue   #We're only calculating for extragalactic transients
            std_mag_R = float(trData[12])
            LCFile    = 'LightCurveFiles/' + trData[5]
            lcdata_up = h5py.File( LCFile + '_UBVRI.hdf5', 'r' )
            lcdata_lo = h5py.File( LCFile + '_%s.hdf5' % self.colorScheme, 'r')
            if self.colorScheme == 'UBVRI':
                for band in self.bands:
                    if band not in ['U','B','V','R','I','J','H','K']:
                        raise TypeError("You're trying to observe in bands that are not in the color_system in %s. Please check your %s and %s" % (self.Opts.params, self.Opts.params, self.Opts.file))
            if int( trData[0] ) == 3: #Has K-corrections
                Kcorfile_up = h5py.File('LightCurveFiles/Kcorrections/%s_%s.hdf5' % (trData[5], self.colorScheme),'r')
                Kcorfile_lo = h5py.File('LightCurveFiles/Kcorrections/%s_%s.hdf5' % (trData[5], self.colorScheme),'r')
                for i,band in enumerate(self.bands):
                    if band.islower():
                        Kcorfile = Kcorfile_lo
                        lcdata = lcdata_lo
                    elif band.isupper():
                        Kcorfile = Kcorfile_up
                        lcdata = lcdata_up
                    if not band in lcdata.keys():
                        print("no data for %s band for %s" % (band, transient))
                    minimum_band = brightest[i]
                    for z in range(len(Kcorfile[band][0,:])):
                        new_minimum = min(lcdata[band][:] + Kcorfile[band][:,z])
                        minimum_band = min(minimum_band, new_minimum)
                    brightest[i] = minimum_band - 3 * std_mag_R
            elif int( trData[0] ) == 4: #Doesn't have K-corrections
                for i,band in enumerate(self.bands):
                    if band.islower():
                        lcdata = lcdata_lo
                    elif band.isupper():
                        lcdata = lcdata_up
                    if not band in lcdata.keys():
                        print("no data for %s band for %s" % (band, transient))
                    minimum_band = min(lcdata[band][:])
                    brightest[i] = min(minimum_band, brightest[i]) -3*std_mag_R
        exp_term = 0.2*( max(self.mag_lim - brightest) ) + 1.0
        Dmax_xgal = np.power( 10.0, exp_term ) # in parsecs
        Dmax_xgal /= 1.0e3 # in kiloparsecs
        Dmax_xgal = max(Dmax_xgal, MW.DMax_MW_kpc)
        print("Maximum distance = %.2e Mpc\n" % (Dmax_xgal/1e3))
        return Dmax_xgal
            
    def set_Framenrlist( self, trtemp ):
        """ Finds in what frames each transient is visible

            Then it stores these framenrs in visibleframes
        """
        if trtemp.transients != []:
            for lc in trtemp.transients:
                for i,grid in enumerate(self.SkyGrids):
                    if grid.WithinFrame(lc.LC_RA, lc.LC_DEC):
                        lc.visibleframes.append(i)
            
    def take_data( self, transientDataSet, outFile, TrTypes, TrTypenrs, nrTrs ):
        """ Takes all the transient observation data and stores it into
             a single array.


            transientDataSet: Set of transient templates with the
                               transients
            outFile: The name of the output file
            TrTypes: A list with transient types that have been 
                      observed in previous observation frames
            TrTypenrs: The corresponding identification numbers of the
                        transient types that have already been observed
            nrTrs: Total number of transients already observed

            In this function a table called allData is created. This 
             table contains a row for every observation of a transient.
            So, if there are three Supernovae Iae observed, each with 
             five observations, 3*5=15 rows are added to allData.
            For each row the following data is noted:
            -trIds: The ID (number) of this transient
            -types: The transient type identification number
            -TimeDat: The time of this observation in JD
            -RADat: The RA coordinate of the transient
            -DECDat: The DEC coordinate of the transient
            -bandDat: The magnitude of the transient in this observation 
            -QmagDat: The quiescence magntiude of this transient in the
                       same color band as in which the observation was 
                       made
            -Obs_Bands: The color band of the observation

            -------
            returns: allData, an updated list of TrTypes, TrTypenrs and 
                      nrTrs
        """
        nTimeSteps = len( np.concatenate(self.obTimes))
        below_threshold = np.max(self.mag_lim) + 0.5 * self.mag_resolution

            
        # how many transients do we have total? This assumes all are visible
        nTr_tot = 0
        for trTemp in transientDataSet.transient_templates:
            self.set_Framenrlist(trTemp)
            if trTemp.N_trans > 0:
                trTemp.sample_all_LCs( np.array(self.obTimes),
                                       np.array(self.obBands), self.threshold )
            nTr_tot += trTemp.N_trans

        # an array for all the magnitude data
        mags = 1000. * np.ones( [nTr_tot, nTimeSteps]  )
        Qmags = 1000. * np.ones( [nTr_tot, nTimeSteps] )
        radec_coords = np.zeros( [nTr_tot,2] )
        Obs_Times = 1.e19 * np.ones( [nTr_tot, nTimeSteps]  )
        Obs_Bands = np.ones((nTr_tot, nTimeSteps), dtype=np.unicode_ )
        # loop over all transients
        i = 0
        it = 0
        imax = 0
        trIds = []
        types = []
        for trTemp in transientDataSet.transient_templates:
            if trTemp.N_trans > 0:
                
                if trTemp.tag not in TrTypes:
                    TrTypes.append(trTemp.tag)
                
                it = self.Trlegend[trTemp.tag]
                radec, bandSets, Quiesmags, Passbands, obTimes = trTemp.take_data(
                                                      self.colorScheme,
                                                      self.Opts.colorsys)
                imax = i + len( radec )
                trIds += range(nrTrs + i, nrTrs + imax)
                radec_coords[i:imax] = radec
                for x,b in enumerate(bandSets):
                    #Check whether the observations make the magnitude limit
                    maglimarray = np.array([self.threshold[band] for band 
                                            in Passbands[x]])
                    b[ b > maglimarray] = below_threshold 
                k = 0
                for j in range(i,imax):
                    len_obs               = len(bandSets[k])
                    mags[j,:len_obs]      = bandSets[k]
                    Qmags[j,:len_obs]     = Quiesmags[k]
                    Obs_Times[j,:len_obs] = obTimes[k]
                    Obs_Bands[j,:len_obs] = Passbands[k]
                    k+=1
                types.extend([it for z in range(i, imax)])
                i = imax
                it += 1
        if i == 0: 	#No transients at all found
            return [], TrTypes, TrTypenrs, nrTrs 
        nTr_tot = nrTrs + imax
        # write out data
        allData = []
        for j in range( nTr_tot - nrTrs ):
            nr_obs = len(Obs_Times[j][Obs_Times[j] < 1.e19])
            for i in range( nr_obs ):
                bandDat = '{0:0.3f}'.format(mags[j,i])
                QmagDat = '{0:0.3f}'.format(Qmags[j,i])
                TimeDat = '{0:0.2f}'.format(Obs_Times[j][i])
                RADEC = radec_coords[j].tolist()
                RADat = '{0:0.4f}'.format(RADEC[0])
                DECDat = '{0:0.4f}'.format(RADEC[1])
                Onerow = ([trIds[j], types[j]] +  [TimeDat] + [RADat] + 
                          [DECDat] + [bandDat] + [QmagDat] + [Obs_Bands[j][i]] )
                allData.append(Onerow)

        return allData,  TrTypes, TrTypenrs, nTr_tot

    def OpenFile(self, outFile, trTypes, trLegend):
        """ Set up the output file.

            First a file with the name outFile is created.
            Then we create a line in that file with a list of the 
             transient types and corresponding identification numbers
             of the types that have been observed.
            Lastly we write the header for the table that will be 
             filled in later.

            -------
            returns: a csv file object in which to write the data
        """
        self.f = open( outFile, 'w' )
        lineTR = []
        for i,Tr in enumerate(trTypes):
            lineTR.append(trLegend[Tr])
            lineTR.append(Tr)
        cwr = csv.writer( self.f, delimiter = '\t' )
        cwr.writerow( lineTR)
        hdrs = ['iD','type','time', 'RA', 'DEC', 'mag', 'Qmag', 'band']
        cwr.writerow(hdrs)
        return cwr

    def WriteToFile(self, cwr, Data):
        """ Write the table Data to the cwr file object
        """
        for row in Data:
            cwr.writerow(row)

    def CloseFile(self):
        """ Close the output file
        """
        self.f.close()

