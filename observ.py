import numpy as np
import csv
from astLib import astCoords
import MilkyWay as MW
import params as pm
import cosmology as cm
import dust


class cell:
    def __init__( self, idx, raMid, decMid, distMid, hdist, hRA, hDEC, z, Grid):
        self.idx    = idx
        self.raMid  = raMid
        self.decMid = decMid
        self.DMid   = distMid      # Luminosity distance of middle of cell
        self.hD     = hdist        # size of the cell in the distance direction
                                   # (differs between galactic and non-galactic)
        self.hRA    = hRA
        self.hDEC   = hDEC
        self.rho    = []
        self.vol    = 0.0
        self.z      = z            # Median redshift of the cell
        self.Grid   = Grid	   # The grid to which this cell belongs to
    # a function to get the maximum variability in density over a cell
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
    def setDensity( self ):
        self.rho = MW.get_MW_dens( [self.raMid, self.decMid, self.DMid], piecewise = True )
    def setVolume( self ):
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
    def sampleCoords( self ):
        z1 = np.random.random()
        z2 = np.random.random()
        z3 = np.random.random()
        randRA   = self.raMid   + (z1 - 0.5)* self.hRA
        randDEC  = self.decMid  + (z2 - 0.5)* self.hDEC
        randDist = self.DMid + (z3 - 0.5)* self.hD
        return [randRA, randDEC, randDist]

class grid:
    def __init__(self, RA_lo, RA_hi, DEC_lo, DEC_hi, Dmax_xgal, N_dist, 
                 N_dist_xgal, Gal_trans, Xgal_trans, obRun):
        rA_lo      = astCoords.hms2decimal(RA_lo,":")
        rA_hi      = astCoords.hms2decimal(RA_hi,":")
        dEC_lo     = astCoords.dms2decimal(DEC_lo,":")
        dEC_hi     = astCoords.dms2decimal(DEC_hi,":")
        if rA_hi < rA_lo: 
            rA_lo, rA_hi = rA_hi, rA_lo
        if dEC_hi < dEC_lo: 
            dEC_lo, dEC_hi = dEC_hi, dEC_lo
        self.RA_lo       = rA_lo
        self.RA_hi       = rA_hi
        self.DEC_lo      = dEC_lo
        self.DEC_hi      = dEC_hi
        self.Dmax_MK     = MW.DMax_MW_kpc
        self.Dmax_xgal   = Dmax_xgal
        self.N_RA        = 6
        self.N_DEC       = 6
        self.N_DMW       = N_dist
        self.N_xgal      = N_dist_xgal
        hRA, hDec        = (rA_hi-rA_lo)/float(self.N_RA), (dEC_hi-dEC_lo)/float(self.N_DEC)
        self.h_DMW       = self.Dmax_MK/float(N_dist)
        self.h_xgal      = (Dmax_xgal - self.Dmax_MK)/float(N_dist_xgal)
        self.Gal_trans	 = Gal_trans    #Whether there are Gal transients
        self.Xgal_trans  = Xgal_trans	#Whether there are Xgal transients
        self.obRun       = obRun
        self.bands       = self.obRun.bands
        if self.Xgal_trans:
            self.Cosmo = cm.Cosmology(Dmax_xgal)
        self.cellGrid    = self.makeCellGrid( hRA, hDec, self.h_DMW, self.h_xgal )
        self.choose_dust_grid()

    def makeCellGrid( self, hRA, hDec, h_DMW, h_xgal ):
        cells = []
        for k in range( self.N_DMW + self.N_xgal ):
            if k < self.N_DMW: 
                mid_D = ( float(k) + 0.5) * h_DMW
                dh = self.h_DMW
                z = 0		#Redshift
            else:          
                if self.Xgal_trans:
                    mid_D = self.Dmax_MK + ( float(k) - self.N_DMW + 0.5 ) * h_xgal
                    dh = h_xgal
                    z = self.Cosmo.get_redshift( mid_D / 1.e3 )	#convert to Mpc
            for j in range( self.N_DEC ):
                mid_dec = self.DEC_lo + (float(j) + 0.5) * hDec
                for i in range( self.N_RA ):
                    mid_ra = self.RA_lo + ( float(i)+0.5 ) * hRA                    
                    idx =  i + self.N_RA*j + self.N_RA*self.N_DEC*k
                    newcell = cell( idx, mid_ra, mid_dec, mid_D, dh, hRA, hDec, z , self )
                    cells.append( newcell )
        return cells   
    def resize( self, deltaRhoMax, maxIter ):
        ri = 0
        while ri < maxIter:
            # isolate the part of the grid that corresponds to the Milky Way
            idxMax = self.N_RA * self.N_DEC * self.N_DMW
            cells_galactic = self.cellGrid[ 0:idxMax ]
            # get the variation in density
            rho_RA_mid  = [c.get_rho_var('RA', 'mid')   for c in cells_galactic]
            rho_RA_edge = [c.get_rho_var('RA', 'edge')  for c in cells_galactic]
            rho_DEC_mid = [c.get_rho_var('DEC', 'mid')  for c in cells_galactic]
            rho_DEC_edge = [c.get_rho_var('DEC', 'edge') for c in cells_galactic]
            rho_RA_mid  = max( rho_RA_mid )
            rho_RA_edge = max( rho_RA_edge )
            rho_DEC_mid = max( rho_DEC_mid )
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
        coordboundaries = (self.RA_lo, self.RA_hi, self.DEC_lo, self.DEC_hi)
        offline = self.obRun.Opts.offline
        
        if pm.use_dust:
            self.Xtr_dust = dust.Schlegel_Extinction(self.bands, self.obRun.colorScheme, offline, *coordboundaries)
            if self.Gal_trans:  #There are galactic transients
                if dust.InSchultheisBoundary(*coordboundaries):
                    self.Gal_dust = dust.Schultheis_Extinction(self.Xtr_dust, self.bands, self.obRun.colorScheme, *coordboundaries)
                elif dust.InGreenBoundary(*coordboundaries):
                    self.Gal_dust = dust.Green_Extinction(self.Xtr_dust, self.bands, offline, self.obRun.colorScheme, *coordboundaries)
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
        self.read_Obsfile(File)
        self.t_start = self.JD[0]
        self.t_end = self.JD[-1]
    def read_Obsfile(self, File):
        with open(File) as Obsfile:
            data = csv.reader(Obsfile, delimiter='\t')
            for i,row in enumerate(data):
                if i not in [0,1]:
                    self.JD.append(float(row[0]) * 8.64e4)
                    self.Band.append(row[1])
    def get_Dates( self ):
        #return np.linspace( self.t_start, self.t_end, self.n_steps )
        return np.array(self.JD) - self.t_start
    def get_ObBands( self ):
        return self.Band
    def get_Bands( self ):
        return np.unique(self.Band)
    def get_relative_times( self ):
        dates = self.get_Dates()
        return dates[-1] - dates

class observation:
    def __init__( self, RA_lo, RA_hi, DEC_lo, DEC_hi, N_dist_MW, N_dist_xgal,\
                     mag_lim, Opts ):
        # set up the transients list
        self.setTransientsList()
        self.Opts = Opts
        # determine max distance (extra-galactic)
        if self.Opts.file:
            self.TimeParams = timeFrame( self.Opts.file )
        else:
            self.TimeParams = timeFrame( pm.ObsFile )
        self.setUpColors( pm.color_system )
        self.threshold = mag_lim
        self.mag_lim = [ mag_lim[ band ] for band in self.bands]
        brightest = 100.0 * np.ones( len( self.bands ) )
        Gal_trans, Xgal_trans = False, False 
        f = open( pm.PeakMagFile, 'r' )
        while f.readline() != '\n':
            pass
        if not self.transientsList == []:
            for line in f:
                fields = line.split()
                if fields[0] in self.transientsList:
                    pkmag = {'U':fields[1], 'B':fields[2], 'V':fields[3], 'R':fields[4], 'I':fields[5], 'J':fields[6], 'H':fields[7], 'K':fields[8], 'u':fields[9], 'g':fields[10], 'r':fields[11], 'i':fields[12], 'z':fields[13]}
                    pkmag = np.array([ pkmag[ band ] for band in self.bands]).astype(float)
                    devmag = float(fields[14])
                    brightest = np.minimum(brightest, pkmag - 3.0*devmag)
                    trData = getFileLine( pm.transientFile, fields[0] )
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
        self.SkyGrid = grid( RA_lo, RA_hi, DEC_lo, DEC_hi, Dmax_xgal, N_dist_MW,\
                                 N_dist_xgal, Gal_trans, Xgal_trans, self )
        #self.SkyGrid.resize( 2.0, 1 )
        self.SkyGrid.setCellValues()
        self.obTimes = self.TimeParams.get_Dates()
        self.obBands = self.TimeParams.get_ObBands()
    def setUpColors( self, col ):
        if col.upper() == 'UBVRI': colr = col.upper()
        elif col.lower() == 'sdss': colr = col.lower()
        elif col.lower() == 'blackgem': colr = col.lower()
        elif col.lower() == 'lsst': colr = col.lower()
        else: print "Invalid color system selected."
        self.colorScheme = colr
        self.bands = self.TimeParams.get_Bands()
        self.nBands = len(self.bands)
    def setTransientsList( self ):
        self.transientsList = []
        if pm.use_nova:     self.transientsList.append( 'nova' )
        if pm.use_CVpolar:  self.transientsList.append( 'CVpolar' )
        if pm.use_CVIP:     self.transientsList.append( 'CVIP' )
        if pm.use_CVdwarf:  self.transientsList.append( 'CVdwarf' )
        if pm.use_AMCVn:    self.transientsList.append( 'AMCVn' )
        if pm.use_SNIa:     self.transientsList.append( 'SNIa' )
        if pm.use_SNIb:     self.transientsList.append( 'SNIb' )
        if pm.use_SNIc:     self.transientsList.append( 'SNIc' )
        if pm.use_SNIIL:    self.transientsList.append( 'SNIIL' )
        if pm.use_SNIIP:    self.transientsList.append( 'SNIIP' )
        if pm.use_SNIInL:   self.transientsList.append( 'SNIInL' )
        if pm.use_SNIInP:   self.transientsList.append( 'SNIInP' )
        if pm.use_M3:       self.transientsList.append( 'M3' )
        if pm.use_M3_5:     self.transientsList.append( 'M3_5' )
        if pm.use_M4:       self.transientsList.append( 'M4' )
    def take_data( self, transientDataSet, outFile ):
        #tEnd = self.obTimes[-1]/8.64e4
        nTimeSteps = len( self.obTimes )
        below_threshold = np.max(self.mag_lim) + 0.5 * pm.mag_resolution
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
                trTemp.sample_all_LCs( self.obTimes, self.obBands, self.threshold )
            nTr_tot += trTemp.N_trans
            print "Nrtrans", trTemp.N_trans

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
                radec, bandSets, Passbands = trTemp.take_data()
                imax = i + len( radec )
                trIds += range(i, imax)
                radec_coords[i:imax] = radec
                #Nrbands = len(bandSets[0])
                for x,b in enumerate(bandSets):
                    maglimarray = np.array([self.threshold[band] for band in Passbands[x]])
                    b[ b > maglimarray] = below_threshold 
                mags[ i:imax,:] = np.array( [np.array(bS).T for bS in bandSets] )
                types.extend([it for i in range(i, imax)])
                i = imax
                it += 1

        if i == 0: 	#No transients at all found
            return
        """
        j=0
        for i,trTemp in enumerate(transientDataSet.transient_templates):
            if trTemp.N_trans > 0:
                Bandsets = mags[ j : j + trTemp.N_trans]
                j = trTemp.N_trans
        """
        nTr_tot = imax
        # set up output file
        f = open( outFile, 'wb' )
        # set up csv
        cwr = csv.writer( f, delimiter = '\t' )
        cwr.writerow( trTypes)# + [1, 'thing_1', 2, 'alien_flashlights'] )
        hdrs = ['iD','type','time', 'RA', 'DEC', 'band']
        cwr.writerow(hdrs)
        # write out data
        for i in range( nTimeSteps ):
            for j in range( nTr_tot ):
                radec = radec_coords[j]
                bandDat = mags[j,i]
                passband = Passbands[j][i]
                allDat = [trIds[j], types[j]] +  [self.obTimes[i]] + radec.tolist() + [bandDat] + [passband]
                cwr.writerow( allDat )
        f.close()

    def writetofile(self, trTemp, bandSets):
        f = open( 'lc_data_test%s.dat' % trTemp.tag, 'wb' )
        clc = csv.writer( f, delimiter = '\t' )
        for k, bS in enumerate(bandSets):
            print bS
            for j in range( len(bS[0]) ):
                print len(self.obTimes), len(bS), k, j
                data = [k] + [self.obTimes[j]] + [band[j] for band in bS]
                clc.writerow(data)
        f.close()

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
