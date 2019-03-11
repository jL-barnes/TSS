import numpy as np
import csv
from astLib import astCoords
import MilkyWay as MW
import cosmology as cm
import dust


class cell:
    def __init__( self, idx, cfl, raMid, decMid, distMid, hdist, hDEC, z, Grid):
        self.idx    = idx
        self.raMid  = raMid
        self.decMid = decMid
        self.DMid   = distMid      # Luminosity distance of middle of cell
        self.hD     = hdist        # size of the cell in the distance direction
                                   # (differs between galactic and non-galactic)
        #self.hRA    = hRA
        self.hDEC   = hDEC
        self.rho    = []
        self.vol    = 0.0
        self.z      = z            # Median redshift of the cell
        self.Grid   = Grid	   # The grid to which this cell belongs to
        self.cellnr_from_left = cfl # There are self.Grid.N_RA cells in a row, this is the cfl'th
    # a function to get the maximum variability in density over a cell
    """
    def get_rho_var( self,raDec_key, midEdge_key ):
        # define coordinates: ra/dec
        ralo, ramid, rahi = self.raMid-0.5*self.hRA, self.raMid, self.raMid+0.5*self.hRA
        declo, decmid, dechi = self.decMid-0.5*self.hDEC, self.decMid, self.decMid+0.5*self.hDEC
        lolo =   [ralo,  declo,  self.DMid]
        lomid =  [ralo,  decmid, self.DMid]
        lohi =   [ralo,  decmid, self.DMid]
        midlo =  [ramid, declo,  self.DMid]
        midmid = [ramid, decmid, self.DMid]
        midhi =  [ramid, dechi,  self.DMid]
        hilo =   [rahi,  declo,  self.DMid]
        himid =  [rahi,  decmid, self.DMid]
        hihi =   [rahi,  dechi,  self.DMid]
        ratio     = 0.0
        ratio_inv = 0.0
        coordPair = []
        if raDec_key == 'RA':
            if midEdge_key == 'mid':
                coordPair = [ [ lolo, midlo ], [ midlo, hilo ], [ lomid, midmid ], \
                              [ midmid, himid ],[ lohi, midhi ], [ midhi, hihi ] ]
            elif midEdge_key == 'edge':
                coordPair = [ [ lolo, hilo ], [ midlo, midhi ], [ hilo, hihi ] ]
        elif raDec_key == 'DEC':
            if midEdge_key == 'mid':
                coordPair = [ [ lolo, lomid ], [ lomid, lohi ], [ midlo, midmid ],\
                              [ midmid, midhi ], [ hilo, himid ], [ himid, hihi ] ]
            elif midEdge_key == 'edge':
                coordPair = [ [ lolo, lohi ], [ midlo, midhi ], [ hilo, hihi ] ]   
        if len(coordPair) < 1: print 'key error in call to to cell::get_rho_var()'
        for pair in coordPair:
            rhoRatio = MW.get_MW_dens( pair[0], piecewise=False )/MW.get_MW_dens( pair[1], piecewise=False )
            ratio = max( rhoRatio, ratio )
            ratio_inv = max( 1.0/rhoRatio, ratio_inv )
        maxRatio = max(ratio, ratio_inv)
        return maxRatio
    """
    def setDensity( self ):
        self.rho = MW.get_MW_dens( [self.raMid, self.decMid, self.DMid], piecewise = True )
    def setVolume2( self ):
        """
        Calculate the (comoving) volume of the cell
        V_C = (1+z)^3 \int D^2 sin(delta) dD dalpha ddelta
        t1: integral over alpha / 3.0
        t2: integral over D * 3.0
        t3: integral over delta
        t4: cosmological factor (=1 if z=0)
        """
        
        t1 = np.deg2rad( self.hRA )/3.0
        t2 = np.power( self.DMid + 0.5*self.hD, 3.0 ) - np.power( self.DMid - 0.5*self.hD, 3.0 );
        cos_DEChi = np.cos( np.deg2rad( self.decMid + 0.5*self.hDEC ) )
        cos_DEClo = np.cos( np.deg2rad( self.decMid - 0.5*self.hDEC ) )
        t3 = np.abs( cos_DEChi - cos_DEClo )
        t4 = np.power( 1 + self.z, 3.)
        self.vol = t1*t2*t3*t4
    def setVolume( self ):
        """
        Calculate the cell volume.
        This is done by subtracting the volume of a sphere with r=Dmax from a
         sphere with r=Dmin. Then we multiply this by the fraction that the FoV
         is of the full sky and the fraction that the cell is of the FoV
        """
        Sphere1 = 4./3. * np.pi * np.power( self.DMid - 0.5*self.hD, 3.0 )
        Sphere2 = 4./3. * np.pi * np.power( self.DMid + 0.5*self.hD, 3.0 )
        Total_vol = Sphere2 - Sphere1
        Fullsky = 129600 / np.pi
        cell_fraction = 1. / (self.Grid.N_RA * self.Grid.N_DEC)
        self.vol = (self.Grid.FoV / Fullsky) * cell_fraction * Total_vol
    def sampleCoordsold( self ):
        z1 = np.random.random()
        z2 = np.random.random()
        z3 = np.random.random()
        randRA   = self.raMid   + (z1 - 0.5)* self.hRA
        randDEC  = self.decMid  + (z2 - 0.5)* self.hDEC
        randDist = self.DMid + (z3 - 0.5)* self.hD
        return [randRA, randDEC, randDist]   
    def sampleCoords( self ):
        z1 = np.random.random()
        z2 = np.random.random()
        z3 = np.random.random()
        randDEC  = self.decMid  + (z1 - 0.5)* self.hDEC
        width_at_DEC = self.Grid.obRun.aperture_RA / np.cos(randDEC * (np.pi/180.)) 
        self.hRA = width_at_DEC / self.Grid.N_RA
        startRA = ( self.Grid.RA_c - 0.5 * width_at_DEC + 
                    self.cellnr_from_left * self.hRA )
        randRA   = startRA + (z2 - 0.5)* self.hRA
        randDist = self.DMid + (z3 - 0.5)* self.hD
        return [randRA, randDEC, randDist]


class grid:
    def __init__(self, pm, Dmax_xgal, Gal_trans, Xgal_trans, gridnr, obRun):
        """
        rA_lo      = astCoords.hms2decimal(pm['RA_lo'],":")
        rA_hi      = astCoords.hms2decimal(pm['RA_hi'],":")
        dEC_lo     = astCoords.dms2decimal(pm['DEC_lo'],":")
        dEC_hi     = astCoords.dms2decimal(pm['DEC_hi'],":")
        if rA_hi < rA_lo: 
            rA_lo, rA_hi = rA_hi, rA_lo
        if dEC_hi < dEC_lo: 
            dEC_lo, dEC_hi = dEC_hi, dEC_lo
        """
        self.gridnr      = gridnr
        self.obRun       = obRun
        self.SetUpFramecoords()
        """
        self.RA_lo       = rA_lo
        self.RA_hi       = rA_hi
        self.DEC_lo      = dEC_lo
        self.DEC_hi      = dEC_hi
        """
        self.Dmax_MK     = MW.DMax_MW_kpc
        self.Dmax_xgal   = Dmax_xgal
        self.N_RA        = 6
        self.N_DEC       = 6
        self.N_cellsRADEC= self.N_RA * self.N_DEC
        self.FoV         = obRun.aperture_DEC * obRun.aperture_RA
        self.N_DMW       = pm['nCells_D']
        self.N_xgal      = pm['nCells_xgal']
        self.h_DMW       = self.Dmax_MK / float(self.N_DMW)
        self.h_xgal      = (Dmax_xgal - self.Dmax_MK) / float(self.N_xgal)
        self.Gal_trans	 = Gal_trans    #Whether there are Gal transients
        self.Xgal_trans  = Xgal_trans	#Whether there are Xgal transients
        self.bands       = self.obRun.bands
        if self.Xgal_trans:
            self.Cosmo = cm.Cosmology(Dmax_xgal)
        #hRA              = (rA_hi-rA_lo) / float(self.N_RA)
        #hDec             = (dEC_hi-dEC_lo) / float(self.N_DEC)
        self.cellGrid    = self.makeCellGrid( self.h_DMW, self.h_xgal )
        self.choose_dust_grid()

    def SetUpFramecoords(self):
        centercoords     = self.obRun.obFrameCoords[self.gridnr]
        RA_c             = centercoords[0]
        DEC_c            = centercoords[1]
        self.RA_c        = astCoords.hms2decimal(RA_c,":")
        self.DEC_c       = astCoords.dms2decimal(DEC_c,":")
        print self.RA_c, self.DEC_c, "RA/DEC"
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
 
        ###Just for now:
        #self.RA_lo = self.RA_lo_l
        #self.RA_up = self.RA_up_l
        
    def makeCellGrid( self, h_DMW, h_xgal ):
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
                    #mid_ra  = self.RA_lo + ( float(i)+0.5 ) * hRA   
                    mid_ra  = ( self.RA_c - 0.5 * width_at_DEC + 
                               (float(i) + 0.5) * Cwidth )
                    idx     = i + self.N_RA*j + self.N_RA*self.N_DEC*k
                    newcell = cell( idx, i, mid_ra, mid_dec, mid_D, dh, hDec, z , self )
                    cells.append( newcell )
        return cells   
    def resize( self, deltaRhoMax, maxIter ):
        """
        This function is deprecated since version 0.9 when cells became non-rectangular
        """
        ri = 0
        while ri < maxIter:
            # isolate the part of the grid that corresponds to the Milky Way
            idxMax = self.N_RA * self.N_DEC * self.N_DMW
            cells_galactic = self.cellGrid[ 0:idxMax ]
            # get the variation in density
            rho_RA_mid   = [c.get_rho_var('RA', 'mid')   for c in cells_galactic]
            rho_RA_edge  = [c.get_rho_var('RA', 'edge')  for c in cells_galactic]
            rho_DEC_mid  = [c.get_rho_var('DEC', 'mid')  for c in cells_galactic]
            rho_DEC_edge = [c.get_rho_var('DEC', 'edge') for c in cells_galactic]
            rho_RA_mid   = max( rho_RA_mid )
            rho_RA_edge  = max( rho_RA_edge )
            rho_DEC_mid  = max( rho_DEC_mid )
            rho_DEC_edge = max( rho_DEC_edge )
            if rho_RA_mid >= deltaRhoMax:
                self.N_RA *= 2
            elif rho_RA_edge >= deltaRhoMax:
                self.N_RA += 1
            if rho_DEC_mid >= deltaRhoMax:
                self.N_DEC *= 2
            elif rho_DEC_edge >= deltaRhoMax:
                self.N_DEC += 1
            if max( rho_RA_mid, rho_RA_edge, rho_DEC_mid, rho_DEC_edge ) < deltaRhoMax:
                print "Grid successfully resized"
                break
            else:
                hRA  = (self.RA_hi - self.RA_lo)/float(self.N_RA)
                hDEC = (self.DEC_hi - self.DEC_lo)/float(self.N_DEC)
                self.cellGrid = self.makeCellGrid( hRA, hDEC, self.h_DMW, self.h_xgal )
                ri += 1
            print "resizing, N_RA = %d, N_DEC = %d" % (self.N_RA, self.N_DEC )
    def setCellValues( self ):
        galIdx = self.N_RA * self.N_DEC * self.N_DMW
        for cell in self.cellGrid:
            cell.setVolume()
            if cell.idx < galIdx:
                cell.setDensity()
    def choose_dust_grid( self ):
        """
        Choose the right dust grid
        For extragalactic transients (and galactic transients with D>60kpc)
         we choose the Schlegel map
        For galactic transients it is possible that the FOV is included in
         the Schultheis as well as the Green map. Then Schultheis gets
         priority
        """
        RA_lo = min(self.RA_lo_l, self.RA_lo_r, self.RA_hi_l, self.RA_hi_r)
        RA_hi = max(self.RA_lo_l, self.RA_lo_r, self.RA_hi_l, self.RA_hi_r)
        coordboundaries = (RA_lo, RA_hi, self.DEC_lo, self.DEC_hi)
        offline = self.obRun.Opts.offline
        
        if not self.obRun.nodust:
            self.Xtr_dust = dust.Schlegel_Extinction(self.bands, self.obRun.colorScheme, offline, *coordboundaries)
            if self.Gal_trans:  #There are galactic transients
                if dust.InSchultheisBoundary(*coordboundaries):
                    self.Gal_dust = dust.Schultheis_Extinction(self.Xtr_dust, self.bands, self.obRun.colorScheme, *coordboundaries)
                elif dust.InGreenBoundary(*coordboundaries):
                    self.Gal_dust = dust.Green_Extinction(self.Xtr_dust, self.bands, self.obRun.colorScheme, offline, *coordboundaries)
                else:
                    self.Gal_dust = self.Xtr_dust
        else:
            print "Excluding dust extinction in the simulation"
            self.Xtr_dust = dust.No_dust(self.bands)
            self.Gal_dust = dust.No_dust(self.bands)
            

class timeFrame:
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
        with open(File) as Obsfile:
            data = csv.reader(Obsfile, delimiter='\t')
            for i,row in enumerate(data):
                if i not in [0,1]:
                    self.JD.append(float(row[0]) * 8.64e4)
                    self.Band.append(row[1])
                    self.RAc.append(row[2])
                    self.DECc.append(row[3])
    def SetUpFrames( self ):
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
        """
        Return an array of the observation dates for each frame
        """
        #return np.linspace( self.t_start, self.t_end, self.n_steps )
        obTimes = np.array(self.JD) - self.t_start
        Dates = []
        for frame in self.Framenrs:
            Dates.append(obTimes[self.Frames == frame])
        return Dates
        #return np.array(self.JD) - self.t_start
    def get_ObBands( self ):
        """
        Return an array of the observation bands for each frame
        """
        Bands = []
        _Band = np.array(self.Band)
        for frame in self.Framenrs:
            Bands.append(_Band[self.Frames == frame])
        return Bands
        #return self.Band
    def get_Bands( self ):
        return np.unique(self.Band)
    def get_relative_times( self ):
        dates = self.get_Dates()
        return dates[-1] - dates
    def get_Frames( self ):
        return self.Frames

class observation:
    def __init__(self, pm, Opts):
        self.Opts           = Opts
        self.pm             = pm
        self.setTransientsList()
        self.Peak_mag_f     = pm['PeakMagFile']
        self.transientFile  = pm['transientFile']
        self.MdwarfFile     = pm['MDwarfFile']
        self.mag_resolution = pm['mag_resolution']
        self.aperture_DEC   = pm['aperture_DEC']
        self.aperture_RA    = pm['aperture_RA']
        self.maxLIGODist    = pm['maxLIGODist']
        self.TimeParams = timeFrame( self.Opts.file )
        #self.NrFrames       = len(np.unique(self.TimeParams.RAc))
        self.nodust         = Opts.nodust
        self.setUpColors( pm['color_system'] )
        self.setUpMag_lim(pm['mag_limit'])
        # determine max distance (extra-galactic)
        brightest = 100.0 * np.ones( len( self.bands ) )
        Gal_trans, Xgal_trans = False, False 
        f = open( self.Peak_mag_f, 'r' )
        while f.readline() != '\n':
            pass
        if not self.transientsList == []:
            for line in f:
                fields = line.split()
                if fields[0] in self.transientsList:
                    pkmag = {'U':fields[1], 'B':fields[2], 'V':fields[3], 'R':fields[4], 'I':fields[5], 'J':fields[6], 'H':fields[7], 'K':fields[8], 'u':fields[9], 'g':fields[10], 'r':fields[11], 'i':fields[12], 'z':fields[13], 'y':fields[14], 'q':fields[15]}
                    pkmag = np.array([ pkmag[ band ] for band in self.bands]).astype(float)
                    devmag = float(fields[16])
                    brightest = np.minimum(brightest, pkmag - 3.0*devmag)
                    trData = getFileLine( self.transientFile, fields[0] )
                    Galtype = int( trData[0] )
                    if Galtype in [0,1,2]:
                        Gal_trans = True	#There are galactic transients
                    if Galtype == 3:
                        Xgal_trans = True	#There are extragalactic transients
        f.close()
        exp_term = 0.2*( max(self.mag_lim - brightest) ) + 1.0
        Dmax_xgal = np.power( 10.0, exp_term ) # parsecs
        Dmax_xgal /= 1.0e3 # kiloparsecs
        print "Dmax_xgal = %.2e Mpc\n" % (Dmax_xgal/1e3)
        # set the properties of the sky window
        self.obTimes       = self.TimeParams.get_Dates()
        self.obBands       = self.TimeParams.get_ObBands()
        self.obFrames      = self.TimeParams.get_Frames()
        self.Nframes       = self.TimeParams.Nr_Frames
        self.obFrameCoords = self.TimeParams.FrameCoords   #Center coordinates of the frames
        self.SkyGrids      = []
        for i in range(self.Nframes):
            self.SkyGrids.append( grid( self.pm, Dmax_xgal, Gal_trans, 
                                       Xgal_trans, i, self ))
            #self.SkyGrids[i].resize( 2.0, 1 )
            self.SkyGrids[i].setCellValues()
    def setUpColors( self, col ):
        if col.upper() == 'UBVRI': colr = col.upper()
        elif col.lower() == 'sdss': colr = col.lower()
        elif col.lower() == 'blackgem': colr = col.lower()
        elif col.lower() == 'lsst': colr = col.lower()
        else: raise Exception(col, "is an invalid color system.")
        self.colorScheme = colr
        self.bands = self.TimeParams.get_Bands()
        self.nBands = len(self.bands)
    def setUpMag_lim(self, mag_limit):
        self.threshold = mag_limit
        for c in self.bands:
            if c not in self.threshold:
                raise Exception("There is no magnitude limit defined for ", c, " in params.py")
        self.mag_lim        = [ self.threshold[ band ] for band in self.bands]
    def setTransientsList( self):
        self.transientsList = []
        if self.pm['use_nova']:     self.transientsList.append( 'nova' )
        if self.pm['use_CVdwarf']:  self.transientsList.append( 'CVdwarf' )
        if self.pm['use_AMCVn']:    self.transientsList.append( 'AMCVn' )
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
    def take_data( self, transientDataSet, outFile, TrTypes, TrTypenrs, nrTrs ):
        """
        nrTrs:  The number of transients that have already been saved
        """
        #tEnd = self.obTimes[-1]/8.64e4
        framenr = transientDataSet.Framenr
        #nTimeSteps = len( self.obTimes )
        nTimeSteps = len( self.obTimes[framenr])
        below_threshold = np.max(self.mag_lim) + 0.5 * self.mag_resolution
        self.convert_toAB   = np.zeros(nTimeSteps, dtype=bool)
        self.convert_toVega = np.zeros(nTimeSteps, dtype=bool)
        #self.convert_toAB   = np.zeros(len(self.obBands), dtype=bool)
        #self.convert_toVega = np.zeros(len(self.obBands), dtype=bool)
        #for i,band in enumerate(self.obBands):
        for i,band in enumerate(self.obBands[framenr]):
            if self.Opts.colorsys == 'AB':
                if band.isupper():
                    self.convert_toAB[i] = True
            if self.Opts.colorsys == 'Vega':
                if band.islower():
                    self.convert_toVega[i] = True
            
        # how many transients do we have total? This assumes all are visible
        """
        nTr_tot = 0
        for trTemp in transientDataSet.transient_templates:
            print "Nrtrans", trTemp.N_trans
            nTr_tot += trTemp.N_trans
        # an array for all the magnitude data
        mags = np.zeros( [nTr_tot, nTimeSteps]  )
        radec_coords = np.zeros( [nTr_tot,2] )
        # loop over all transients
        i = 0
        it = 0
        trIds = []
        trTypes = []
        types = []
        for trTemp in transientDataSet.transient_templates:
            if trTemp.N_trans > 0:
                trTypes +=  [it, trTemp.tag]
                radec, bandSets, bandbands = trTemp.sample_all_LCs( self.obTimes, self.obBands, self.threshold )
                if len(radec) == 0: 	#No transients for this type
                    it +=1
                    continue
                imax = i + len( radec )
                trIds += range(i, imax)
                radec_coords[i:imax] = radec
                #Nrbands = len(bandSets[0])
                for b in bandSets:
                    print b
                mags[ i:imax,:] = np.array( [np.array(bS).T for bS in bandSets] )
                print mags
                types.extend([it for i in range(i, imax)])
                i = imax
                it += 1
        """
        nTr_tot = 0
        for trTemp in transientDataSet.transient_templates:
            if trTemp.N_trans > 0:
                trTemp.sample_all_LCs( self.obTimes[framenr],
                                      self.obBands[framenr], self.threshold )
            nTr_tot += trTemp.N_trans
            print "Nrtrans", trTemp.tag, trTemp.N_trans

        # an array for all the magnitude data
        #mags = np.zeros( [nrTrs + nTr_tot, nTimeSteps]  )
        #radec_coords = np.zeros( [nrTrs + nTr_tot,2] )
        mags = np.zeros( [nTr_tot, nTimeSteps]  )
        radec_coords = np.zeros( [nTr_tot,2] )
        # loop over all transients
        #i = nrTrs
        i = 0
        it = 0
        imax = 0
        trIds = []
        types = []
        for trTemp in transientDataSet.transient_templates:
            if trTemp.N_trans > 0:
                
                if trTemp.tag not in TrTypes:
                    TrTypes.append(trTemp.tag)
                    """
                    if TrTypenrs == []:
                        it = 0
                    else:
                        it = max(TrTypenrs) + 1
                    TrTypenrs.append(it)
                    """
                
                it = self.Trlegend[trTemp.tag]
                radec, bandSets, Passbands = trTemp.take_data(self.colorScheme,
                                        self.convert_toAB, self.convert_toVega)
                imax = i + len( radec )
                #trIds += range(i, imax)
                trIds += range(nrTrs + i, nrTrs + imax)
                radec_coords[i:imax] = radec
                #Nrbands = len(bandSets[0])
                for x,b in enumerate(bandSets):
                    maglimarray = np.array([self.threshold[band] for band in Passbands[x]])
                    b[ b > maglimarray] = below_threshold 
                mags[ i:imax,:] = np.array( [np.array(bS).T for bS in bandSets] )
                types.extend([it for z in range(i, imax)])
                i = imax
                it += 1
        if i == 0: 	#No transients at all found
            return [], TrTypes, TrTypenrs, nrTrs
        """
        j=0
        for i,trTemp in enumerate(transientDataSet.transient_templates):
            if trTemp.N_trans > 0:
                Bandsets = mags[ j : j + trTemp.N_trans]
                j = trTemp.N_trans
        """
        #nTr_tot = imax  
        nTr_tot = nrTrs + imax
        """
        # set up output file
        f = open( outFile, 'wb' )
        # set up csv
        cwr = csv.writer( f, delimiter = '\t' )
        cwr.writerow( trTypes)# + [1, 'thing_1', 2, 'alien_flashlights'] )
        hdrs = ['iD','type','time', 'RA', 'DEC', 'band']
        cwr.writerow(hdrs)
        """
        # write out data
        allData = []
        for i in range( nTimeSteps ):
            #print "hello", nTimeSteps
            for j in range( nTr_tot - nrTrs ):
                #radec = radec_coords[j]
                #bandDat = mags[j,i]
                #print mags[j,i], j, i
                bandDat = '{0:0.3f}'.format(mags[j,i])
                #TimeDat = '{0:0.2f}'.format(self.obTimes[i])
                TimeDat = '{0:0.2f}'.format(self.obTimes[framenr][i])
                RADEC = radec_coords[j].tolist()
                RADat = '{0:0.4f}'.format(RADEC[0])
                DECDat = '{0:0.4f}'.format(RADEC[1])
                #passband = Passbands[j][i]
                #Onerow = [trIds[j], types[j]] +  [self.obTimes[i]] + radec.tolist() + [bandDat] + [passband]
                #Onerow = [trIds[j], types[j]] +  [self.obTimes[i]] + radec.tolist() + [bandDat] + [self.obBands[i]]
                Onerow = [trIds[j], types[j]] +  [TimeDat] + [RADat] + [DECDat] + [bandDat] + [self.obBands[framenr][i]]
                #print Onerow
                allData.append(Onerow)
                #cwr.writerow( allDat )
        #f.close()
        return allData,  TrTypes, TrTypenrs, nTr_tot
    """
    def writetofile(self, trTemp, bandSets):
        f = open( 'lc_data_test%s.dat' % trTemp.tag, 'wb' )
        clc = csv.writer( f, delimiter = '\t' )
        for k, bS in enumerate(bandSets):
            for j in range( len(bS[0]) ):
                print len(self.obTimes), len(bS), k, j
                data = [k] + [self.obTimes[j]] + [band[j] for band in bS]
                clc.writerow(data)
        f.close()
    """
    def OpenFile(self, outFile, trTypes, trLegend):
        #set up output file
        f = open( outFile, 'wb' )
        lineTR = []
        for i,Tr in enumerate(trTypes):
            #lineTR.append(i)
            lineTR.append(trLegend[Tr])
            lineTR.append(Tr)
        # set up csv
        cwr = csv.writer( f, delimiter = '\t' )
        cwr.writerow( lineTR)# + [1, 'thing_1', 2, 'alien_flashlights'] )
        hdrs = ['iD','type','time', 'RA', 'DEC', 'band']
        cwr.writerow(hdrs)
        return cwr

    def WriteToFile(self, cwr, Data):
        for row in Data:
            cwr.writerow(row)


def getFileLine( file, tag ):
    f = open( file, 'r' )
    while f.readline() != '\n':
        pass
    data = f.readline()
    while data.split()[0] != tag:
        data = f.readline()
    restOfData = data.split()[1:]
    f.close()
    return restOfData
