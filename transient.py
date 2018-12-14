import numpy as np
import h5py
import coordTransforms as cT
import params as pm
import MilkyWay as MW
import cosmology as cm
import glob
from scipy.interpolate import interp1d


year_to_seconds = 8.64e4 * 365.25
day_to_seconds = 8.64e4
pc_to_cm = 3.08567758e18
emptyval = 1000.0

#Flux zero-points in erg cm^-2 s^-1 A^-1
flux_0 = {'UBVRI': {'U':417.5e-11, 'B':632.0e-11, 'V':363.1e-11, 'R':217.7e-11, 'I':112.6e-11},
          'sdss':  {'u':847.e-11, 'g':490.e-11, 'r':287.e-11, 'i':195.e-11, 'z':137.e-11},
          'blackgem': {'u':754.e-11, 'g':472.e-11, 'r':279.e-11, 'i':186.e-11, 'z':130.e-11, 'q':324.e-11},
          'lsst': {'u':799.e-11, 'g':473.e-11, 'r':284.e-11, 'i':193.e-11, 'z':145.e-11, 'y':115.e-11}
         }
#http://www.astronomy.ohio-state.edu/~martini/usefuldata.html


class TransientSet:
    def __init__ (self, obRun):
        self.grid = obRun.SkyGrid
        self.transient_templates = []
        for tr in obRun.transientsList:
            new_template = makeTemplate( tr, obRun )
            self.transient_templates.append( new_template )
        self.deltaT_obs  = (obRun.TimeParams.t_end - obRun.TimeParams.t_start)/8.64e4
        self.setUpLC_funs( obRun.colorScheme )
    def setUpLC_funs(self, colorchoice):
        print 'setting up functional broadband forms...'
        for trtemp in self.transient_templates:
            trtemp.setUpLC( colorchoice )
        print 'done!'
    def populate(self):
        # builds lists of light curve blue prints that are used to
        # generate observations when take_data() is called
        galIdx = self.grid.N_RA * self.grid.N_DEC * self.grid.N_DMW
        self.populate_galactic(galIdx)
        self.populate_xgal(galIdx)
    def populate_galactic(self, galIdx):
        galTrans = [ tr for tr in self.transient_templates if tr.galType != 3 ]
        cells_galactic = self.grid.cellGrid[0:galIdx]
        for gtr in galTrans:
            D_gal_max = gtr.Max_obs_D * 1.e-3 #kpc
            for c in cells_galactic:
                if D_gal_max < c.DMid + c.hD: break	#A transient past this distance won't be visible anyway
                gtr.get_blueprints( c, self.deltaT_obs)
    def populate_xgal( self, galIdx ):
        xgalTrans = [ tr for tr in self.transient_templates if tr.galType == 3 ]
        cells_xgal = self.grid.cellGrid[galIdx:]
        for xtr in xgalTrans:
            self.grid.Cosmo.create_Nr_density_funct(xtr.NM, xtr.alpha) 
            for c in cells_xgal:
                xtr.get_blueprints( c, self.deltaT_obs )
    def inject_kilonova( self, tMerge, dist, mvRed, mvBlue, obRun ):
        self.transient_templates.append( kilonovaTemplate( tMerge, dist, mvRed, mvBlue, obRun ))
            
def makeTemplate( tag, obRun ):
    transientsData = getFileLine( pm.transientFile, tag )
    PD = getFileLine( pm.PeakMagFile, tag )
    PeakMagData = {'U':PD[0], 'B':PD[1], 'V':PD[2], 'R':PD[3], 'I':PD[4], 
                   'J':PD[5], 'H':PD[6], 'K':PD[7], 'u':PD[8], 'g':PD[9], 
                   'r':PD[10], 'i':PD[11], 'z':PD[12], 'std':PD[13]}
    Type = int( transientsData[0] )
    if Type == 0 :
        temp = galactic_nRec_template( tag, transientsData, PeakMagData, obRun )
    elif Type == 1:
        temp = galactic_recur_template( tag, transientsData, PeakMagData, obRun )
    elif Type == 2:
        temp = Mdwarf_template( tag, transientsData, PeakMagData, obRun)
    else:
        temp = xgal_template( tag, transientsData, PeakMagData, obRun)
    return temp

class transientTemplate:
    def __init__( self, tag, transientsData, PeakMagData, obRun ):
        self.tag         = tag
        self.galType     = int( transientsData[0] )
        self.galkey      = [int(tD) for tD in transientsData[1:5]]
        self.bands       = obRun.bands
        self.mag_lim     = obRun.mag_lim
        self.peak_mag    = np.array([ PeakMagData[ band ] for band in self.bands]).astype(float)
        self.std_mag     = float( PeakMagData['std'] )
        self.LCFile      = 'LightCurveFiles/' + transientsData[5]
        self.broadbands  = {} 
        self.stellarNrm  = float( transientsData[6] )/MW.rho_stellar_sun
        self.scaleH      = float( transientsData[7] )
        self.deltaT_LC    = float( transientsData[8] )
        self.N_trans     = 0    # how many transients each template has
        self.transients  = []   # list of light curves blueprints here
        self.N_transnonsee = 0  #Nr of generated transients, regardless of visibility
        self.radec_list  = []
        self.bandmags_list= []
        self.bandbands_list=[]
    def setUpLC( self, colorchoice ):
        print 'for %s' % self.tag
        if colorchoice not in ['sdss', 'UBVRI', 'blackgem', 'lsst']:
            print 'Invalid color system.'
            return 0
        lcdata = h5py.File( self.LCFile + '_%s.hdf5' % colorchoice, 'r' )
        tms = lcdata['times'][:]
        self.broadbands['t'] = tms;
        for band in self.bands:
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
    def get_all_bbs( self, tsample, obbands ):
        bnds = []
        for i,bandLabel in enumerate(obbands):
            this_bb = self.broadbands[bandLabel]( tsample[i] )
            bnds.append(this_bb)
        return bnds
    def get_Ntransients( self, cell, Nr_of_years ):
        if abs(self.scaleH - MW.MW_Hthin) < 0.001:
            N_per_y = self.stellarNrm * cell.vol * np.sum( np.array(cell.rho)[ np.where( self.galkey )] )
        else:
            this_rho = MW.get_MW_dens( [cell.raMid, cell.decMid, cell.DMid],\
                                           piecewise=True, Hthin = self.scaleH )
            N_per_y = self.stellarNrm * cell.vol * np.sum( np.array(this_rho)[ np.where( self.galkey )] )
        N_per_y = N_per_y * Nr_of_years
        if N_per_y > 2.0:
            NTr = int( N_per_y )
        else: NTr = geometricP( N_per_y )
        return NTr   
    def get_blueprints( self, c, dtObs ):
        Nr_of_years = (self.deltaT_LC + dtObs) / 365.25   #At first observation the first LC
                                                          #should have just ended
        N_transients = self.get_Ntransients( c, Nr_of_years )
        self.N_transnonsee += N_transients
        for i in range( N_transients ):
            tOutburst = -365.25 * Nr_of_years * np.random.random()
            if tOutburst + self.deltaT_LC + dtObs > 0.0:
                magDelta = np.random.normal( scale = self.std_mag )
                peak_M = self.peak_mag + magDelta
                newLC = LC_blueprint( c.sampleCoords(), peak_M, magDelta, tOutburst )
                if self.canSee_LC( newLC ):
                    coords = (newLC.LC_RA, newLC.LC_DEC, newLC.LC_Dist)
                    newLC.Extinction = c.Grid.Gal_dust.Sample_extinction(*coords)
                    if self.canSee_LC( newLC ):
                        self.transients.append(newLC)
                        self.N_trans += 1
    def sample_all_LCs( self, obTimes, obBands, threshold ):
        obT_f = obTimes[-1]
        dts = obT_f - obTimes
        for lc in self.transients: 
            visible = False
            tSamp = -lc.tExplode * 8.64e4 - dts
            dval = lc.LC_Dist	#is in kpc
            dist_mod = 5.0 * np.log10(dval*1.0e3) - 5.0
            magDev = lc.magDelta
            bandMags = self.get_all_bbs( tSamp, obBands )
            for i,label, in enumerate(obBands):
                bandMags[i] += dist_mod + lc.Extinction[label] + magDev
                if bandMags[i] < threshold[label]: visible = True
            if visible:
                self.radec_list.append( [ lc.LC_RA, lc.LC_DEC ] )
                self.bandmags_list.append( bandMags )
                self.bandbands_list.append( obBands )
            else:
                self.N_trans -=1
        self.bandmags_list = np.array(self.bandmags_list)
    def take_data(self):
        """
        Return the observational data to observ.py
        """
        return self.radec_list, self.bandmags_list, self.bandbands_list
    def canSee_LC(self, lightcurve):
        """
        A crude way to estimate whether the light curve should still be visible
        If the Extinction is not zero, it has been adjusted to a real value.
        That lets us calculate the visibility of the transient WITH extinction.
        """
        seen = False
        d_10pc = 0.01
        dist_ratio = lightcurve.LC_Dist/d_10pc
        apparent_mag = lightcurve.PeakMag + 5.0 * np.log10( dist_ratio )
        if lightcurve.Extinction == 0:
            for band, x in enumerate(self.mag_lim):
                if apparent_mag[ band ] < x: 
                    seen = True
                    return seen
        else:
            for band, x in enumerate(self.mag_lim):
                color = self.bands[band]
                if apparent_mag[ band ] + lightcurve.Extinction[color] < x: 
                    seen = True
                    return seen
        return seen

class galactic_nRec_template( transientTemplate ):
    def __init__( self, tag, trData, PkMData, obRun ):
        transientTemplate.__init__( self, tag, trData, PkMData, obRun  )
        self.Max_obs_D     =  Maximum_observing_distance(self.mag_lim, self.peak_mag, self.std_mag)

class galactic_recur_template( transientTemplate ):
    def __init__( self, tag, trData, PkMData, obRun ):
        transientTemplate.__init__( self, tag, trData, PkMData, obRun)
        self.freqDay = float( trData[9] )
        self.freqDev = float( trData[10] )
        self.Max_obs_D     =  Maximum_observing_distance(self.mag_lim, self.peak_mag, self.std_mag)
    def get_blueprints( self, c, dtObs ):
        Nr_of_years = 1.	#Recurrent transient => Nr. of transients does not depend on 
                                #observation time
        N_transients = self.get_Ntransients( c, Nr_of_years )
        self.N_transnonsee += N_transients
        for i in range( N_transients ):
            tOutburst = -365.25 * Nr_of_years * np.random.random()
            magDelta = self.std_mag * 3.0
            peak_M = self.peak_mag + magDelta
            newLC = LC_blueprint( c.sampleCoords(), peak_M, magDelta, tOutburst )
            if self.canSee_LC( newLC):
                coords = (newLC.LC_RA, newLC.LC_DEC, newLC.LC_Dist)
                newLC.Extinction = c.Grid.Gal_dust.Sample_extinction(*coords)
                if self.canSee_LC( newLC ):
                    self.transients.append( newLC )
                    self.N_trans += 1
    def sample_all_LCs( self, obTimes, obBands, threshold ):
        obT_f = obTimes[-1]
        dts = obT_f - obTimes
        for lc in self.transients:
            radec, mags = self.sample_single_LC( lc, dts, obBands, threshold )
            if len(radec) > 0:
                self.radec_list.append( radec )
                self.bandmags_list.append( mags )
                self.bandbands_list.append( obBands )
        self.bandmags_list = np.array(self.bandmags_list)
    def sample_single_LC(self, lc, tvals, obBands, threshold ):
        # figure out the times for each outburst
        tOuts = []
        t_OutBurst = lc.tExplode
        # step forward in time
        while t_OutBurst < self.deltaT_LC + (tvals[-1] / 8.64e4):
            tOuts.append( t_OutBurst )
            t_OutBurst += np.random.normal( loc = self.freqDay, scale = self.freqDev)
        # step backward in time
        t_OutBurst = lc.tExplode - np.random.normal( loc = self.freqDay, scale = self.freqDev)
        while t_OutBurst + self.deltaT_LC > (-tvals[0] / 8.64e4): 
            tOuts.insert(0, t_OutBurst)
            t_OutBurst -= np.random.normal( loc = self.freqDay, scale = self.freqDev)
        # construct the aggregate light curve
        if len( tOuts ) > 0:
            tOuts = tOuts[::-1]	#flip the array to start with the closest outburst
            tO = tOuts[0]
            tSamp = -tvals - tO*8.64e4
            # deviation for this occurrence of the transient
            magDev = np.random.normal( scale = self.std_mag )
            bandMags = self.get_all_bbs( tSamp, obBands )
            for band in bandMags:
                band[ band < emptyval ] +=  magDev
            for tO in tOuts[1:]:
                if tO > (-tvals[0] / 8.64e4) - self.deltaT_LC:
                    tSamp = -tvals - tO*8.64e4
                    magDev = np.random.normal( scale = self.std_mag )
                    these_bandMags = self.get_all_bbs( tSamp, obBands )
                    for bi, band in enumerate( these_bandMags ):
                        band[ band < emptyval ] += magDev
                        oband = bandMags[bi]
                        i3 = ( band != emptyval ) & ( oband == emptyval )
                        i4 = ( band != emptyval ) & ( oband != emptyval )
                        bandMags[bi][i3] = band[i3]
                        bandMags[bi][i4] = -2.5 * np.log10( np.power(10, -.4*band[i4]) +\
                                                        np.power(10, -.4*oband[i4])  )
        # apply distance modulus
        dval = lc.LC_Dist
        dist_mod = 5.0 * np.log10( dval*1.0e3 ) - 5.0
        visible = False
        for i,label, in enumerate(obBands):
            bandMags[i] += dist_mod + lc.Extinction[label]
            if bandMags[i] < threshold[label]: visible = True
        if visible :
            return [lc.LC_RA, lc.LC_DEC], bandMags
        else:
            self.N_trans -= 1
            return [],[]

class Mdwarf_template( transientTemplate ):
    def __init__( self, tag, trData, PkMData, obRun ):
        transientTemplate.__init__( self, tag, trData, PkMData, obRun )
        vals               = np.genfromtxt( pm.MDwarfFile, names=True ).T[tag]
        self.aval, self.bval, self.alf_ac, self.alf_in,\
            self.bet_ac, self.bet_in, self.da_dz = vals[0:7]
        self.Ebins         = 100		#results in 100 bins
        
        logE_min_ac, logE_max_ac = vals[7:9]
        logE_min_in, logE_max_in = vals[9:11]
        qs                 = vals[11:]
        last_logE_ac       = logE_max_ac - (logE_max_ac - logE_min_ac)/self.Ebins
                             #last_logE is the the start of the last logE-bin
        last_logE_in       = logE_max_in - (logE_max_in - logE_min_in)/self.Ebins
        
        self.dlogE_act     = np.linspace(logE_min_ac, last_logE_ac, self.Ebins +2)
        self.dlogE_ina     = np.linspace(logE_min_in, last_logE_in, self.Ebins +2)
        self.Lq            = {'U':qs[0],'B':qs[1],'V':qs[2],'R':qs[3],'I':qs[4],'J':qs[5],'H':qs[6],'K':qs[7],'u':qs[8],'g':qs[9],'r':qs[10],'i':qs[11],'z':qs[12]}
        self.Flare_energy  = float( trData[11] )
        self.Epk           = {}
        self.Max_obs_D     =  Maximum_observing_distance(self.mag_lim, self.peak_mag, self.std_mag)
        self.broadbands_rise  = {}
        self.broadbands_decay = {}
    def get_blueprints( self, c, dtObs ):
        Nr_of_years = 1.	#Recurrent transient => Nr. of transients does not depend on 
                                #observation time
        N_transients = self.get_Ntransients( c, Nr_of_years )
        self.N_transnonsee += N_transients
        for i in range( N_transients ):
            t_dummy = -365.25 * Nr_of_years * np.random.random()
            newLC = LC_blueprint( c.sampleCoords(), self.peak_mag, 0.0, t_dummy )
            if self.canSee_LC( newLC ):
                coords = (newLC.LC_RA, newLC.LC_DEC, newLC.LC_Dist)
                newLC.Extinction = c.Grid.Gal_dust.Sample_extinction(*coords)
                if self.canSee_LC( newLC):	#Test WITH dust
                    self.transients.append( newLC )
                    self.N_trans += 1
    def setUpLC( self, colorchoice ):
        print 'for %s' % self.tag
        self.flux_0 = flux_0[colorchoice]
        Peak1, Peak2 = (0,0)
        self.LC_time = {}
        for TYPE in ['rise', 'decay']:
            lcdata = h5py.File( self.LCFile + '_%s_%s.hdf5' % (colorchoice, TYPE), 'r' )
            tms = lcdata['times'][:]
            if TYPE == 'rise':self.broadbands_rise['t'] = tms;
            else: self.broadbands_decay['t'] = tms
            for band in self.bands:
                print "%s band" % band
                if not band in lcdata.keys():
                    print "no data for %s band" % band
                    bb = np.zeros( len(tms) )
                else:
                    bb = lcdata[band][:]
                ff = interp1d( tms, bb, bounds_error = False, fill_value = 0.0 )
                def fun( times, fcn = ff ):
                    return fcn(times)
                if TYPE == 'rise':
                    self.broadbands_rise[band] = fun
                    Peak1 = np.max(bb)
                    self.LC_time['rise'] = tms[-1]
                elif TYPE == 'decay':
                    self.broadbands_decay[band] = fun
                    Peak2 = np.max(bb)
                    self.LC_time['decay'] = tms[-1]
                self.Epk[band] = max(Peak1, Peak2)
        lcdata.close()
    def sample_all_LCs( self, obTimes, obBands, threshold ):
        # total time for observation in hours
        tWindow = ( obTimes[-1] - obTimes[0] )/3600.0 
        for mdwarf in self.transients:
            visible = False
            # determine activity status and set constants
            active = False
            # galactic height in parsecs
            hZ = cT.get_galactic_height( mdwarf.LC_RA, mdwarf.LC_DEC, mdwarf.LC_Dist ) * 1.0e3 # parsecs
            P_A = self.aval * np.exp( -self.bval * abs(hZ) )  
            ztest = np.random.random()
            if ztest <= P_A: active = True
            if active: 
                alf = self.alf_ac + self.da_dz * abs(hZ) / 1.e3	#da_dz is in kpc^-1
                bet = self.bet_ac
                eBins_l = self.dlogE_act	#lowest logE value of the bins
                eBins_m = self.dlogE_act + 0.5 * (eBins_l[1] - eBins_l[0])
                                                #Middle value of the logE bins
            else: 
                alf = self.alf_in + self.da_dz * abs(hZ) / 1.e3
                bet = self.bet_in
                eBins_l = self.dlogE_ina
                eBins_m = self.dlogE_act + 0.5 * (eBins_l[1] - eBins_l[0])
            #Determine Flare Frequency
            cum_ne = tWindow * np.power(10.0, alf) * np.power(10.0, bet*eBins_l)
            #print alf, bet, "cum_ne", cum_ne, "eBins_l", eBins_l
            #cum_ne is cumulative. We need to decumulate it:
            ne = [cum_ne[i] - cum_ne[i+1] for i in range(self.Ebins +1)]
            ne = np.array(ne)
            z = np.random.random(self.Ebins + 1)  #Round probabilities off to integers
            cond = z < ne - ne.astype(int)
            ne[cond] +=1
            nFlares = ne.astype(int)
            #print nFlares
            
            if max( nFlares ) < 1: 
                self.N_trans -=1
            else:
                lumlist = self.get_all_lums( nFlares, eBins_m[:-1], obTimes, obBands )
                maglist = np.zeros( len(obBands) )
                for j, lum in enumerate(lumlist):
                    bnd = obBands[j]
                    area = 4.0 * np.pi * np.power( mdwarf.LC_Dist * 1e3 * pc_to_cm, 2.0 ) # cm^2
                    mag = -2.5 * np.log10( lum/area/self.flux_0[ bnd ] )
                    maglist[j] = mag
                    #the luminosity (in mags) of the quiescent star
                    mag_non_outb = -2.5 * np.log10( self.Lq[bnd]/area/self.flux_0[ bnd ] )
                    maglist[j] += mdwarf.Extinction[ bnd ]
                    if mag < threshold[bnd] and mag < mag_non_outb - pm.mag_resolution:
                        visible = True
                if visible:
                    """
                    Boolarray = np.array(obBands) == 'u'
                    ubands = maglist[Boolarray]
                    import matplotlib.pyplot as plt	#plot a single Mdwarf lightcurve
                    plt.plot(obTimes[Boolarray], ubands)
                    plt.gca().invert_yaxis()
                    plt.show()
                    """
                    self.radec_list.append([ mdwarf.LC_RA, mdwarf.LC_DEC ])
                    self.bandmags_list.append( maglist )   
                    self.bandbands_list.append( obBands )                 
                else:
                    self.N_trans -=1
        self.bandmags_list = np.array(self.bandmags_list)
    def get_all_lums(self, nFlares, eBins, obT, obbands ):
        window = obT[-1] - obT[0]
        lumlist = np.zeros( len(obT) )
        eF, nF = eBins[nFlares>0], nFlares[nFlares>0]
        for i, n in enumerate(nF):
            logEU = eF[i]
            trise  = np.power(10.0, np.random.normal(0.31*logEU - 7.6, 0.32))
            tdecay = np.power(10.0, np.random.normal(0.51*logEU - 13., 0.22))
            thisE = np.power(10.0, logEU)
            dts = np.random.random( n ) * window
            #generate outbursts that started after t0
            for j,dt in enumerate(dts):
                tSamp = obT-dt
                bnd = obbands[j]
                if self.Epk[bnd] != 0:
                    tSamp_rise  = tSamp * self.LC_time['rise'] / trise
                    tSamp_decay = (tSamp - trise) *\
                                   self.LC_time['decay'] / tdecay
                    new_flux_rise  = ( thisE / self.Flare_energy ) *\
                                      self.broadbands_rise[bnd]( tSamp_rise )
                    new_flux_decay = ( thisE / self.Flare_energy ) *\
                                      self.broadbands_decay[bnd]( tSamp_decay )
                    lumlist += new_flux_rise + new_flux_decay
            #generate outbursts that started before t0
            a = 1
            for j,x in enumerate(obT):	#How many observations can a burst comprise
                if x > self.deltaT_LC:
                    a = j
                    break
            dts = np.random.random( n ) * window
            dts = dts[ dts < self.deltaT_LC * 24.*60.*60.]    #dts in seconds, deltaT_LC in days
            for j,dt in enumerate(dts):
                tSamp = obT[:a]+dt
                bnd = obbands[j]
                if self.Epk[bnd] != 0:
                    tSamp_rise  = tSamp * self.LC_time['rise'] / trise
                    tSamp_decay = (tSamp - trise) *\
                                   self.LC_time['decay'] / tdecay
                    new_flux_rise  = ( thisE / self.Flare_energy ) *\
                                      self.broadbands_rise[bnd]( tSamp_rise )
                    new_flux_decay = ( thisE / self.Flare_energy ) *\
                                      self.broadbands_decay[bnd]( tSamp_decay )
                    lumlist += new_flux_rise + new_flux_decay
        for i,bnd in enumerate(obbands):
            lumlist[i] += self.Lq[ bnd ]
        return lumlist
                

class xgal_template( transientTemplate  ):
    def __init__( self, tag, trData, PkMData, obRun ):
        transientTemplate.__init__( self, tag, trData, PkMData, obRun )
        #self.stellarNrm = 2.47e-14	#for snIa
        self.NM = float(trData[6])
        self.alpha = float(trData[12])
        self.Kcormodel = trData[5]
        #The following are only applicable to SNIa
        self.scale = 1		
        self.scale_0 = 1.09
        self.M_B_max = -18.89
    def Phillips_scale( self, magDev ):
        Mag_zero = self.M_B_max + magDev - 5 * np.log10(cm.H0) + 25
        scale = 1./(0.59 * Mag_zero + 3.0)	#The real scale factor
        self.scale = self.scale_0 / scale	#Adjusted for scale factor of template LC
    def get_all_bbs( self, tsample , obbands, redshift):
        bnds = []
        for i, bandLabel in enumerate(obbands):
            Kcor = self.Kcorfuns.Sample_Kcor(bandLabel, tsample[i] * self.scale,  
                                             redshift)
            this_bb = self.broadbands[bandLabel]( tsample[i] * self.scale ) + float(Kcor)
            bnds.append(float(this_bb))
        return bnds
    def get_Ntransients( self, cell, Nr_of_years ):
        kpc3 = cell.vol
        N_per_y = kpc3 * cell.Grid.Cosmo.Nr_density_xtr( cell.z ) 
        #N_per_y = kpc3 * self.stellarNrm	#Non-cosmological version
        N_per_y = N_per_y * Nr_of_years
        if N_per_y > 2.0:
            NTr = int( N_per_y )
        else: NTr = geometricP( N_per_y )
        return NTr
    def get_blueprints( self, c, dtObs):
        Nr_of_years = (self.deltaT_LC + dtObs) / 365.25   #At first observation the first LC
                                                          #should have just ended
        N_transients = self.get_Ntransients( c, Nr_of_years )
        self.N_transnonsee += N_transients
        for i in range( N_transients ):
            tOutburst = -365.25 * Nr_of_years * np.random.random() 
            if tOutburst + self.deltaT_LC + dtObs > 0.0:
                magDelta = np.random.normal( scale = self.std_mag )
                peak_M = self.peak_mag + magDelta
                newLC = LC_blueprint( c.sampleCoords(), peak_M, magDelta, tOutburst )
                if self.canSee_LC( newLC ):
                    coords = (newLC.LC_RA, newLC.LC_DEC)
                    newLC.Extinction = c.Grid.Xtr_dust.Sample_extinction(*coords)
                    if self.canSee_LC( newLC ):
                        newLC.Redshift = c.Grid.Cosmo.get_redshift(newLC.LC_Dist / 1.e3)	#convert to Mpc
                        self.transients.append(newLC)
                        self.N_trans += 1
    def sample_all_LCs( self, obTimes, obBands, threshold ):
        obT_f = obTimes[-1]
        dts = obT_f - obTimes
        self.Kcorfuns = cm.Kcor(self.Kcormodel, obBands)
        for lc in self.transients: 
            visible = False
            tSamp = -lc.tExplode * 8.64e4 - dts
            dval = lc.LC_Dist	#is in kpc
            dist_mod = 5.0 * np.log10(dval*1.0e3) - 5.0
            magDev = lc.magDelta
            if self.tag == 'SNIa':
                self.Phillips_scale(magDev)
            bandMags = self.get_all_bbs( tSamp, obBands, lc.Redshift )
            for i,label in enumerate(obBands):
                bandMags[i] += dist_mod + lc.Extinction[label] + magDev
                if bandMags[i] < threshold[label]: visible = True
            #for i,label, in enumerate(bandlabels):
            #    bandMags[i][bandMags[i] < emptyval] += dist_mod + lc.Extinction[label] + magDev
            #bbExt = np.array([np.min( band ) for band in bandMags])
            #if np.any( bbExt < threshold ) :

            if visible: 
                self.radec_list.append( [ lc.LC_RA, lc.LC_DEC ] )
                self.bandmags_list.append( bandMags )
                self.bandbands_list.append( obBands )
            else:
                self.N_trans -=1
        self.bandmags_list = np.array(self.bandmags_list)


class kilonovaTemplate( transientTemplate ):
    def __init__(self, tMerge, dist, mvRed, mvBlue, obRun):
        self.tag        = 'kilonova'
        self.bands      = obRun.bands
        self.t_merge    = tMerge
        if dist != None:
            self.dist = dist
        else: self.dist = sampleLIGODistance( pm.maxLIGODist)
        self.RA         = sampleRA( pm.RA_lo, pm.RA_hi )
        self.DEC        = sampleDEC( pm.DEC_lo, pm.DEC_hi )
        self.N_trans    = 1
        self.N_transnonsee = 1
        self.mvRed      = mvRed
        self.mvBlue     = mvBlue
        self.f_band     = {}
        self.netLC      = {}
        self.bandmags_list = np.ones(len(obRun.obBands)) * emptyval

    def sample_all_LCs( self, obTimes, obBands, threshold ):
        if self.mvRed != None:
            mred, vred = get_best_match( 'red', *self.mvRed)
            rlcf = 'knRed_m%.4f_vk%.2f_Xlan1e-1.0.h5' % (mred, vred)
            rlcf = 'kilonovaLib/' + rlcf
            redLC = h5py.File(rlcf, 'r' )
            if self.mvBlue == None: 
                self.netLC['times'] = redLC['days'][:]
                for band in self.bands:
                    self.netLC[band] =  redLC[band][:]
                redLC.close()
        if self.mvBlue != None:
            mblu, vblu = get_best_match( 'blue', *self.mvBlue)
            blcf = 'knBlu_m%.4f_vk%.2f_Xlan1e-5.0.h5' % (mblu, vblu)
            blcf = 'kilonovaLib/' + blcf
            bluLC = h5py.File(blcf, 'r')
            if self.mvRed == None:
                self.netLC['times'] = bluLC['days'][:] 
                for band in self.bands:
                    self.netLC[band] =  bluLC[band][:]
                bluLC.close()
        if not self.netLC: # combine kilonova light curves
            self.netLC['times'] = redLC['days'][:]
            for band in self.bands:
                m1, m2 = redLC[band][:], bluLC[band][:]
                m12 = -2.5*np.log10( np.power(10.0, -.4*m1) + np.power(10.0, -.4*m2) )
                self.netLC[band] = m12
            redLC.close()
            bluLC.close()
        # now that we have the light curve, sample it at the appropriate times
        tSample = (obTimes-self.t_merge)/8.64e4
        dist_mod = 5.0 * np.log10( self.dist/pc_to_cm  ) - 5.0
        for j, band in enumerate(self.bands):
            bandLC = self.netLC[band]
            self.f_band[band] = interp1d( self.netLC['times'][ np.isinf(bandLC) == False],
                                   bandLC[np.isinf(bandLC) == False], 
                                   bounds_error=False, fill_value=np.inf )
        visible = False
        for j,band in enumerate(obBands):
            self.bandmags_list[j] = self.f_band[band]( tSample[j] ) + dist_mod #+ Extinction
            if self.bandmags_list[j] < threshold[band]:visible = True
        #The following three items should be put into an array like with all other transients
        visible = False
        if visible:
            self.radec_list = np.array([[self.RA, self.DEC]])
            self.bandmags_list = np.array([self.bandmags_list])
            self.bandbands_list =  [obBands]
        else:	
            self.N_trans = 0


class LC_blueprint:
    def __init__ (self, coords, magPk, magDelta, tBurst):
        # explosion time: randomly within the year preceding 
        # the end of the observation
        self.tExplode   = tBurst
        self.LC_RA      = coords[0]
        self.LC_DEC     = coords[1]
        self.LC_Dist    = coords[2]
        self.magDelta   = magDelta
        self.PeakMag    = magPk
        self.Extinction = 0
        self.Redshift   = 0	#preliminarily






def Maximum_observing_distance(mag_lim, pkmag, devmag):
    exp_term = 0.2*( max( mag_lim - ( pkmag - 3.0*devmag ) ) ) + 1.0
    Max_D = np.power( 10.0, exp_term ) # parsecs
    return Max_D

def Minimum_Flare_E(D, Lq,  obRun):
    """
    Calculates the minimum flare energy that one could observe of an M dwarf
    This is an absolute lower limit. It is very likely that some brighter flares
     can't be seen either. However, better be safe than sorry.
    The minimum flare energy depends on the quiescence energy, the magnitude 
     limit in the U-band, the magnitude resolution and the distance of the M dwarf.
    This can only be calculated for the U-band. We apply this to all other
     bands. That's okay, because M-dwarfs flare the brightest in the U-band

    D: Distance to transient in kpc
    Lq: quiescent luminosities per colorband in erg/s/Ang
    """
    colorbands = obRun.bands
    #########################################################
    #This is not finished. First we need to model the rise/decay time
    #####################################################
    D = D * 1.e3
    Area = 4.0 * np.pi * np.power( D * pc_to_cm, 2.0 )
    #Quiescent fluxes
    Qflux = { band : Lq[band] / Area for band in colorbands}
    #Quiescent magnitudes
    f_0 = flux_0[obRun.colorScheme]
    Qmag = { band : flux2mag(Qflux[band], f_0[band]) for band in colorbands}

    #The difference between the magnitude limit and quiescent magnitude of the M dwarf
    difmag = {band : Qmag[band] - obRun.mag_limit[band] for band in colorbands}
    mindifmag, mindifband = min( zip(difmag.values(), difmag.keys()) )
    if mindifmag < 0:	#The M dwarf is visible in quiescence
        logminE = 28
        #####################################
        #Here a function should come to investigate the minimum flare required for an M dwarf that is visible. I.e. work with the sensitivity of the telescope: pm.mag_resolution
        ######################33
    else:		#The M dwarf is not visible in quiescence
        #Now calculate how large a flare would be in the U-band to increase 
        # the M dwarf's brightness by mindifmag magnitudes
        Qflux_U = Lq['U'] / Area
        Qmag_U = flux2mag( Qflux_U, f_0['U'] )
        Outburstflux_U = mag2flux( Qmag_U - mindifmag, f_0['U'])
        L_outburst_U = Outburstflux_U * Area
        print Qmag_U, mindifmag, mag2flux(Outburstflux_U, f_0['U']), L_outburst_U
        logminE = np.log10(L_outburst_U - Lq['U']) #####################################
        ############The above function now gives an erg/s: it must be changed to erg with the rise/decay time of a flare
        ###################
    return logminE

def mag2flux(mag,f_0):
    """
    converts magnitude to flux in Jy
    """
    flux = f_0 * 10**(-0.4*mag)
    return flux

def flux2mag(flux,f_0):
    """
    converts flux (Jy) to magnitudes
    """
    mag = -2.5*np.log10(flux/f_0)
    return mag


def geometricP( nt ):
    """
    Use Poisson statistics to obtain a number of transients inside a cell
     if the average is 2. Because we're using Poisson statistics, there is
     a non-zero probability that the number of transients is higher than 2.
    nt: the average number of transients in this cell (float)
    returns: the actual number of transients in this cell (int)
    """
    return np.random.poisson(nt)
    """
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
    """



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



def order_mv():
    ms, vs = [],[]
    for f in glob.glob( 'kilonovaLib/knBlu*.h5' ) + glob.glob( 'kilonovaLib/knRed*.h5'):
        this_m = f[f.find('_m')+2:f.find('_v')]
        this_v = f[f.find('_v')+3:f.find('_X')]
        if this_m not in ms: ms.append(this_m)
        if this_v not in vs: vs.append(this_v)
    ms = np.sort(np.array( [float(m) for m in ms] ))
    vs = np.sort(np.array( [float(v) for v in vs] ))
    return ms, vs

def get_best_match( br, m_msun, v_c ):
    ms, vs = order_mv()
    m_match = ms[np.argmin( np.abs(ms-m_msun) )]
    v_match = vs[np.argmin( np.abs(vs-v_c) )]
    print( 'Best match for %s kilonova part: M = %.2e Msun, v = %.2f c.' % (br, m_match, v_match) )
    return m_match, v_match
    
def sampleRA(ra_lo, ra_hi):
    ra_conv = np.array([15.0, 0.25, 1.0/240.])
    hmslow = np.array([float(val) for val in ra_lo.split(':')])
    hmshi= np.array([float(val) for val in ra_hi.split(':')])
    ralo = np.sum( hmslow*ra_conv )
    rahi = np.sum( hmshi*ra_conv )
    if ralo > rahi: ralo, rahi = rahi, ralo
    RA_sampled = np.random.random() * (rahi - ralo ) + ralo
    return RA_sampled

def sampleDEC( dec_lo, dec_hi ):
    dec_conv = np.array([1.0, 1.0/60.0, 1.0/3600.0])
    dmslow = np.array([float(val) for val in dec_lo.split(':')])
    dmshi= np.array([float(val) for val in dec_hi.split(':')])
    declo = np.sum( dmslow*dec_conv )
    dechi = np.sum( dmshi*dec_conv )
    if declo > dechi: declo, dechi = dechi, declo
    DEC_sampled = np.random.random() * (dechi - declo) + declo
    return DEC_sampled

def sampleLIGODistance(  dmax ):
    while True:
        x = np.random.random()*dmax
        y = np.random.random()*dmax
        z = np.random.random()*dmax
        r = np.sqrt( x*x + y*y + z*z )
        if r <= dmax: break
    return r

