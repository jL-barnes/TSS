import numpy as np
from scipy.interpolate import interp1d
import h5py
import params as pm
import MilkyWay as MW

year_to_seconds = 8.64e4 * 365.25
emptyval = pm.mag_limit + 10.0

f0_UBVRI = {'U':417.5, 'B':632.0, 'V':363.1, 'R':217.7, 'I':112.6, 'J':31.47, 'H':11.38, 'K':3.961}
f0_ugriz = {'u':859.5, 'g':466.9, 'r':278.0, 'i':185.2, 'z':131.5}

class TransientSet:
    def __init__ ( self, obRun ):
        self.grid = obRun.SkyGrid
        self.transient_templates = []
        for tr in obRun.transientsList:
            new_template = makeTemplate( tr )
            self.transient_templates.append( new_template )
        self.mag_lim = obRun.mag_lim
        self.deltaT_obs  = (obRun.TimeParams.t_end - obRun.TimeParams.t_start)/8.64e4
        print 'setting up functional broadband forms...'
        self.setUpLC_funs( obRun.colorScheme, obRun.bands)
        print 'done!'
    def setUpLC_funs(self, colorchoice, bands):
        for trtemp in self.transient_templates:
            trtemp.setUpLC( colorchoice, bands )
    def populate(self):
        # builds lists of light curve blue prints that are used to
        # generate observations when take_data() is called
        galIdx = self.grid.N_RA * self.grid.N_DEC * self.grid.N_DMW
        self.populate_galactic( galIdx )
        self.populate_xgal( galIdx )
    def populate_galactic( self, galIdx ):
        galTrans = [ tr for tr in self.transient_templates if tr.galType != 3 ]
        cells_galactic = self.grid.cellGrid[0:galIdx]
        for c in cells_galactic:
            for gtr in galTrans: 
                gtr.get_blueprints( c, self.deltaT_obs, self.mag_lim )
    def populate_xgal( self, galIdx ):
        xgalTrans = [ tr for tr in self.transient_templates if tr.galType == 3 ]
        cells_xgal = self.grid.cellGrid[galIdx:]
        for c in cells_xgal:
            for xtr in xgalTrans:
                xtr.get_blueprints( c, self.deltaT_obs, self.mag_lim )
            
def makeTemplate( tag ):
    transientsData = getFileLine( pm.transientFile, tag )
    type = int( transientsData[0] )
    if type == 0 :
        temp = galactic_nRec_template( tag, transientsData )
    elif type == 1:
        temp = galactic_recur_template( tag, transientsData )
    elif type == 2:
        temp = Mdwarf_template( tag, transientsData )
    else:
        temp = xgal_template( tag, transientsData )
    return temp

class transientTemplate:
    def __init__( self, tag, transientsData ):
        self.tag         = tag
        self.galType     = int( transientsData[0] )
        self.galkey      = [int(tD) for tD in transientsData[1:5]]
        self.peak_mag    = float( transientsData[5] )
        self.std_mag     = float( transientsData[6] )
        self.LCFile      = 'LightCurveFiles/' + transientsData[7]
        self.broadbands  = {} 
        self.colors      = ''
        self.stellarNrm  = float( transientsData[8 ] )/MW.rho_stellar_sun
        self.scaleH      = float( transientsData[9 ] )
        self.deltaT_LC    = float( transientsData[10] )
        self.N_trans     = 0    # how many transients each template has
        self.transients  = []   # list of light curves blueprints here
        self.N_transnonsee = 0
    def setUpLC( self, colorchoice, bands ):
        print 'for %s' % self.tag
#        if colorchoice == 'UBVRI': 
#            UBRVI, ugriz = True, False
#            self.colors = 'UBVRI'
#        elif colorchoice == 'ugriz': 
#            UBRVI, ugriz = False, True
#            self.colors = 'ugriz'
#        self.allBands = list(bands)
#        else:
        if colorchoice not in ['ugriz', 'UBVRI']:
            print 'Invalid color system.'
            return 0
        self.allBands = list(bands)
        lcdata = h5py.File( self.LCFile + '_%s.hdf5' % colorchoice, 'r' )
        tms = lcdata['times'][:]
        self.broadbands['t'] = tms;
        for band in self.allBands:
            print "%s band" % band
            if not band in lcdata.keys():
                print "no data for %s band" % band
                bb = emptyval * np.ones( len(tms) )
            else:
                bb = lcdata[band][:]
            ff = interp1d( tms, bb, bounds_error = False, fill_value = emptyval )
            def fun( times, fcn = ff ):
                return fcn(times)
            self.broadbands[band] = fun
        lcdata.close()
    def get_all_bbs( self, tsample ):
        bnds = []
        for bandLabel in self.allBands:
            this_bb = self.broadbands[bandLabel]( tsample )
            bnds.append(this_bb)
        return bnds
#        if self.colors == 'UBVRI':
#            uband = self.broadbands['U']( tsample )
#            bband = self.broadbands['B']( tsample )
#            vband = self.broadbands['V']( tsample )
#            rband = self.broadbands['R']( tsample )
#            iband = self.broadbands['I']( tsample )
#            jband = self.broadbands['J']( tsample )
#            hband = self.broadbands['H']( tsample )
#            kband = self.broadbands['K']( tsample )
#            bnds = [uband, bband, vband, rband, iband, jband, hband, kband]
#        elif self.colors == 'ugriz':
#            uband = self.broadbands['u']( tsample )
#            gband = self.broadbands['g']( tsample )
#            rband = self.broadbands['r']( tsample )
#            iband = self.broadbands['i']( tsample )
#            zband = self.broadbands['z']( tsample )
#            bnds  = [uband, gband, rband, iband, zband]
#        return bnds
    def get_Ntransients( self, cell ):
        if abs(self.scaleH - MW.MW_Hthin) < 0.001:
            N_per_y = self.stellarNrm * cell.vol * np.sum( np.array(cell.rho)[ np.where( self.galkey )] )
        else:
            this_rho = MW.get_MW_dens( [cell.raMid, cell.decMid, cell.DMid],\
                                           piecewise=True, Hthin = self.scaleH )
            N_per_y = self.stellarNrm * cell.vol * np.sum( np.array(this_rho)[ np.where( self.galkey )] )
        if N_per_y > 2.0:
            NTr = int( N_per_y )
        else: NTr = geometricP( N_per_y )
        return NTr   
    def get_blueprints( self, c, dtObs, mag_lim ):
        N_transients = self.get_Ntransients( c )
        self.N_transnonsee += N_transients
        for i in range( N_transients ):
            tOutburst = -365.25 * np.random.random()
            if tOutburst + self.deltaT_LC + dtObs > 0.0:
                magDelta = np.random.normal( scale = self.std_mag )
                peak_M = self.peak_mag + magDelta
                newLC = LC_blueprint( c.sampleCoords(), peak_M, magDelta, tOutburst )
                if canSee_LC( newLC, mag_lim ):
                    self.transients.append( newLC )
                    self.N_trans += 1
    def sample_all_LCs( self, obTimes, threshold ):
        radec_list, bandmags_list = [], []
        obT_f = obTimes[-1]
        dts = obT_f - obTimes
        for lc in self.transients: 
            tSamp = -lc.tExplode * 8.64e4 - dts
            dval = lc.LC_Dist
            dist_mod = 5.0 * np.log10(dval*1.0e3) - 5.0
            magDev = lc.magDelta
            bandMags = self.get_all_bbs( tSamp )
            for band in bandMags:
                band[band < emptyval] += dist_mod + magDev
            bbExt = min( [np.min( band ) for band in bandMags] )
            if bbExt > threshold:
                continue
            else:
                radec_list.append( [ lc.LC_RA, lc.LC_DEC ] )
                bandmags_list.append( bandMags )
        return radec_list, bandmags_list

class galactic_nRec_template( transientTemplate ):
    def __init__( self, tag, trData ):
        transientTemplate.__init__( self, tag, trData  )

class galactic_recur_template( transientTemplate ):
    def __init__( self, tag, trData ):
        transientTemplate.__init__( self, tag, trData )
        self.freqDay = float( trData[11] )
        self.freqDev = float( trData[12] )
    def get_blueprints( self, c, dtObs, mag_lim ):
        N_transients = self.get_Ntransients( c )
        for i in range( N_transients ):
            tOutburst = -365.25 * np.random.random()
            magDelta = self.std_mag * 3.0
            peak_M = self.peak_mag + magDelta
            newLC = LC_blueprint( c.sampleCoords(), peak_M, magDelta, tOutburst )
            if canSee_LC( newLC, mag_lim ):
                self.transients.append( newLC )
                self.N_trans += 1
    def sample_all_LCs( self, obTimes, threshold ):
        radec_list, bandmags_list = [], []
        obT_f = obTimes[-1]
        dts = obT_f - obTimes
        for lc in self.transients:
            radec, mags = self.sample_single_LC( lc, dts )
            if len(radec) > 0:
                radec_list.append( radec )
                bandmags_list.append( mags )
        return radec_list, bandmags_list
    def sample_single_LC( lc, tvals ):
        # figure out the times for each outburst
        tOuts = []
        t_OutBurst = lc.tExplode
        # step forward in time
        while tOutBurst < 0.0:
            tOuts.append( t_OutBurst )
            t_OutBurst += np.random.normal( loc = self.freqDay, scale = self.freqDev)
        # step backward in time
        t_OutBurst = lc.tExplode - np.random.normal( loc = self.freqDay, scale = self.freqDev)
        while t_OutBurst + self.deltaT_LC > -dts[0]: 
            tOuts.insert(0, t_OutBurst)
            t_OutBurst -= np.random.normal( loc = self.freqDay, scale = self.freqDev)
        # construct the aggregate light curve
        if len( tOuts > 0 ):
            tO = tOuts[0]
            tSamp = -tvals - tO*8.64e4
            # deviation for this occurrence of the transient
            magDev = np.random.normal( scale = self.std_mag )
            bandMags = self.get_all_bbs( tSamp )
            for band in bandMags:
                band[ band < emptyval ] +=  magDev
            for tO in tOuts[1:]:
                tSamp = -tvals - t0*8.64e4
                magDev = np.random.normal( scale = self.std_mag )
                these_bandMags = self.get_all_bbs( tSamp )
                for bi, band in enumerate( these_bandMags ):
                    band[ band < emptyval ] += magDev
                    oband = bandMags[bi]
                    i3 = ( band != emptyVal ) & ( oband == emptyVal )
                    i4 = ( band != emptyVal ) & ( oband != emptyVal )
                    bandMags[bi][i3] = band[i3]
                    bandMags[bi][i4] = -2.5 * np.log10( np.power(10, -.4*band[i4]) +\
                                                        np.power(10, -.4*oband[i4])  )
        # apply distance modulus
        dval = lc.LC_Dist
        dist_mod = 5.0 * np.log10( dval*1.0e3 ) - 5.0
        for band in bandMags:
            band += dist_mod
        bbExt = min([np.min(band) for band in bandMags])
        if bbExt > threshold:
            return [],[]
        else:
            return [lc.LC_RA, lc.LC_DEC], bandMags

class Mdwarf_template( transientTemplate ):
    def __init__( self, tag, trData ):
        transientTemplate.__init__( self, tag, trData )
        vals = np.genfromtxt( pm.MDwarfFile, names=True ).T[tag]
        self.aval, self.bval, self.alf_ac, self.alf_in,\
            self.bet_ac, self.bet_in, self.da_dz = vals[0:7]
        logE_min, logE_max = vals[7:9]
        self.dlogE_act = np.arange( logE_min, logE_max, .1 )
        logE_min, logE_max = vals[9:11]
        self.dlogE_ina = np.arange( logE_min, logE_max, .1 )
        qs = vals[11:]
        self.Lq = {'U':qs[0],'B':qs[1],'V':qs[2],'R':qs[3],'I':qs[4],'J':qs[5],'H':qs[6],'K':qs[7]}
        if self.colors == 'ugriz': self.Lq = get_ugriz_quiescents( self.Lq )
        self.Epk = {}
    def get_blueprints( self, c, dtObs, mag_lim ):
        N_transients = self.get_Ntransients( c )
        for i in range( N_transients ):
            t_dummy = -365.25 * np.random.random()
            newLC = LC_blueprint( c.sampleCoords(), self.peak_mag, 0.0, t_dummy )
            if canSee_LC( newLC, mag_lim ):
                self.transients.append( newLC )
                self.N_trans += 1
    def setUpLC( self, colorchoice, bands ):
        print 'for %s' % self.tag
        if colorchoice == 'UBVRI': 
            self.flux_0 = f0_UBVRI
        elif colorchoice == 'ugriz': 
            self.flux_0 = f0_ugriz
        else:
            print 'Invalid color system.'
            return 0
        lcdata = h5py.File( self.LCFile + '_%s.hdf5' % colorchoice, 'r' )
        tms = lcdata['times'][:]
        self.broadbands['t'] = tms;
        for band in allBands:
            print "%s band" % band
            if not band in lcdata.keys():
                print "no data for %s band" % band
                bb = np.zeros( len(tms) )
            else:
                bb = lcdata[band][:]
            ff = interp1d( tms, bb, bounds_error = False, fill_value = 0.0 )
            def fun( times, fcn = ff ):
                return fcn(times)
            self.broadbands[band] = fun
            self.Epk[band] = np.max(bb)
        lcdata.close()
    def sample_all_LCs( self, obTimes, threshold ):
        # total time for observation in hours
        tWindow = ( obTimes[-1] - obtimes[0] )/3600.0 
        # lists to hold the output
        radecs, bandmags_list = [], []
        for mdwarf in self.transients:
            # determine activity status and set constants
            active = False
            # galactic height in parsecs
            hZ = MW.get_galactice_height( mdwarf.LC_DEC, mdwarf.LC_RA, mdwarf.LC_Dist ) * 1.0e3 # parsecs
            P_A = self.aval * np.exp( self.bval * hZ )  
            ztest = np.random.random()
            if ztest <= P_A: active = True
            if active: 
                alf = self.alf_ac + self.da_dz * hZ
                bet = self.bet_ac
                eBins = self.dlogE_act
            else: 
                alf = self.alf_in + self.da_dz * hZ
                bet = self.bet_in
                eBins = self.dlogE_ina
            #Determine Flare Frequency
            ne = tWindow * np.power(10.0, alf) * np.power(10.0, bet*eBins)
            z = np.random.random(len(eBins))  #Round probabilities off to integers
            cond = z < ne - ne.astype(int)
            ne[cond] +=1
            nFlares = ne.astype(int)

            if max( nFlares ) < 1: continue
            else:
                lumlist = self.get_all_lums( nFlares, eBins, obTimes )
                for j, lum in enumerate(lumlist):
                    area = 4.0 * np.pi * np.power( mdwarf.LC_dist*1e3, 2.0 ) # pc^2
                    lum = -2.5 * np.log10( lum/area/self.flux_0[ self.allBands[j] ] )
                    lumlist[j] = lum
                if min( np.min(lumin) for lumin in lumlist ) < mag_lim:
                    radec.append([ mdwarf.LC_RA, mdwarf.LC_DEC ])
                    bandmags_list.append( lumlist )
        return radec, bandmags_list
    def get_all_lums( nFlares, eBins, obT ):
        lumlist = []
        window = obT[-1] - obT[0]
        nBands = len(self.allBands)
        lum_matrix = np.zeros( [ nBands, len(obT) ] )
        eF, nF = eBins[nFlares>0], nFlares[nFlares>0]
        for i, n in enumerate(nF):
            thisE = np.power(10.0, eF[i])
            dts = np.random.random( n ) * window
            for dt in dts:
                tSamp = obT-dt
                for j, bnd in enumerate(nBands):
                    new_flux = self.broadbands[bnd]( tSamp )
                    new_flux[ new_flux > 0.0 ] += (thisE - self.Epk[bnd])
                    lum_matrix[j] += new_flux
        for i,row in enumerate(lum_matrix):
            lum_matrix[i] += self.Lq[ self.allBands[i] ]
        lumlist = [row for row in lum_matrix]
        return lumlist
                
    #def setUpLCs():
        # get unscaled light curves 
    #   return #-

class xgal_template( transientTemplate ):
    def __init__( self, tag, trData ):
        transientTemplate.__init__( self, tag, trData )
        self.deltaT_LC = 400.0
        self.stellarNrm = float(trData[8])
    def get_Ntransients( self, cell ):
        kpc3 = cell.vol
        N_per_y = kpc3 * self.stellarNrm
        if N_per_y > 2.0:
            NTr = int( N_per_y )
        else: NTr = geometricP( N_per_y )
        return NTr
#    def get_blueprints( self, c, dtObs, mag_lim ):
#        N_transients = self.get_Ntransients( c )
#        for i in range( N_transients ):
#            tOutBurst = -365.25 * np.random.random()
#            if tOutBurst + self.deltaT_LC + dtObs > 0.0:
#                magDelta = np.random.normal( scale = self.std_mag )
#                peak_M = self.peak_mag + magDelta
#                newLC = LC_blueprint( c.sampleCoords(), peak_M, tOutBurst )
#                if canSee_LC( newLC, mag_lim ):
#                    self.transients.append( newLC )


class LC_blueprint:
    def __init__ (self, coords, magPk, magDelta, tBurst ):
        # explosion time: randomly within the year preceding 
        # the end of the observation
        self.tExplode   = tBurst
        self.LC_RA      = coords[0]
        self.LC_DEC     = coords[1]
        self.LC_Dist    = coords[2]
        self.magDelta   = magDelta
        self.PeakMag    = magPk



def canSee_LC( lightcurve, magLimit ):
    # a crude way to estimate whether the light curve should still be visible
    seen = True
    d_10pc = 0.01
    dist_ratio = lightcurve.LC_Dist/d_10pc
    apparent_mag = lightcurve.PeakMag + 5.0 * np.log10( dist_ratio )
    if apparent_mag > magLimit:
        seen = False
    return seen

def geometricP( nt ):
    Nt = 0
    if nt > 1.0:
        p = 0.5*nt
        z1 = np.random.geometric( p, size=2 )
        if   (z1 == 1).sum() == 2: Nt = 2
        elif (z1 == 1).sum() == 1: Nt = 1
    else:
        z1 = np.random.geometric( nt, size = 1 )
        if z1 == 1: Nt = 1
    return Nt 


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
