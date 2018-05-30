import numpy as np
import h5py
import coordTransforms as cT
import params as pm
import MilkyWay as MW
import cosmology as cm
import time
import glob
from scipy.interpolate import interp1d


year_to_seconds = 8.64e4 * 365.25
day_to_seconds = 8.64e4
pc_to_cm = 3.08567758e18
emptyval = 1000.0

flux_0 = {'U':417.5e-11, 'B':632.0e-11, 'V':363.1e-11, 'R':217.7e-11, 'I':112.6e-11, 'J':31.47e-11, 'H':11.38e-11, 'K':3.961e-11, 'u':859.5e-11, 'g':466.9e-11, 'r':278.0e-11, 'i':185.2e-11, 'z':131.5e-11}
#http://www.astronomy.ohio-state.edu/~martini/usefuldata.html


class TransientSet:
    def __init__ (self, obRun):
        self.grid = obRun.SkyGrid
        self.transient_templates = []
        self.mag_lim = obRun.mag_lim
        for tr in obRun.transientsList:
            new_template = makeTemplate( tr, self.mag_lim )
            self.transient_templates.append( new_template )
        self.deltaT_obs  = (obRun.TimeParams.t_end - obRun.TimeParams.t_start)/8.64e4
        print 'setting up functional broadband forms...'
        self.setUpLC_funs( obRun.colorScheme )
        print 'done!'
    def setUpLC_funs(self, colorchoice):
        for trtemp in self.transient_templates:
            trtemp.setUpLC( colorchoice )
    def populate(self):
        # builds lists of light curve blue prints that are used to
        # generate observations when take_data() is called
        galIdx = self.grid.N_RA * self.grid.N_DEC * self.grid.N_DMW
        self.populate_galactic(galIdx)
        self.populate_xgal(galIdx)
    def populate_galactic(self, galIdx):
        galTrans = [ tr for tr in self.transient_templates if tr.galType != 3 ]
        cells_galactic = self.grid.cellGrid[0:galIdx]
        i = 0
        tic = time.time()
        for gtr in galTrans:
            D_gal_max = gtr.Max_obs_D * 1.e-3 #kpc
            for c in cells_galactic:
                if D_gal_max < c.DMid + c.hD: break	#A transient past this distance won't be visible anyway
                gtr.get_blueprints( c, self.deltaT_obs, self.mag_lim )
        print "time", time.time() - tic
    def populate_xgal( self, galIdx ):
        xgalTrans = [ tr for tr in self.transient_templates if tr.galType == 3 ]
        cells_xgal = self.grid.cellGrid[galIdx:]
        for xtr in xgalTrans:
            self.grid.Cosmo.create_Nr_density_funct(xtr.NM, xtr.alpha) 
            for c in cells_xgal:
                xtr.get_blueprints( c, self.deltaT_obs, self.mag_lim )
    def inject_kilonova( self, tMerge, dist, mvRed, mvBlue ):
        self.transient_templates.append( kilonovaTemplate( tMerge, dist, mvRed, mvBlue ))
            
def makeTemplate( tag, mag_lim ):
    transientsData = getFileLine( pm.transientFile, tag )
    PD = getFileLine( pm.PeakMagFile, tag )
    PeakMagData = {'U':PD[0], 'B':PD[1], 'V':PD[2], 'R':PD[3], 'I':PD[4], 'J':PD[5], 'H':PD[6], 'K':PD[7], 'u':PD[8], 'g':PD[9], 'r':PD[10], 'i':PD[11], 'z':PD[12], 'std':PD[13]}
    Type = int( transientsData[0] )
    if Type == 0 :
        temp = galactic_nRec_template( tag, transientsData, PeakMagData, mag_lim )
    elif Type == 1:
        temp = galactic_recur_template( tag, transientsData, PeakMagData, mag_lim )
    elif Type == 2:
        temp = Mdwarf_template( tag, transientsData, PeakMagData, mag_lim )
    else:
        temp = xgal_template( tag, transientsData, PeakMagData )
    return temp

class transientTemplate:
    def __init__( self, tag, transientsData, PeakMagData ):
        self.tag         = tag
        self.galType     = int( transientsData[0] )
        self.galkey      = [int(tD) for tD in transientsData[1:5]]
        self.bands       = list( pm.showbands )
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
    def setUpLC( self, colorchoice ):
        print 'for %s' % self.tag
        if colorchoice not in ['ugriz', 'UBVRI']:
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
    def get_all_bbs( self, tsample ):
        bnds = []
        bandlabels = []
        for bandLabel in self.bands:
            this_bb = self.broadbands[bandLabel]( tsample )
            bnds.append(this_bb)
            bandlabels.append(bandLabel)
        return bnds, bandlabels
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
    def get_blueprints( self, c, dtObs, mag_lim ):
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
                if canSee_LC( newLC, mag_lim ):
                    coords = (newLC.LC_RA, newLC.LC_DEC, newLC.LC_Dist)
                    newLC.Extinction = c.Grid.Gal_dust.Sample_extinction(*coords)
                    if canSee_LC( newLC, mag_lim ):
                        self.transients.append(newLC)
                        self.N_trans += 1
    def sample_all_LCs( self, obTimes, threshold ):
        radec_list, bandmags_list = [], []
        obT_f = obTimes[-1]
        dts = obT_f - obTimes
        for lc in self.transients: 
            tSamp = -lc.tExplode * 8.64e4 - dts
            dval = lc.LC_Dist	#is in kpc
            dist_mod = 5.0 * np.log10(dval*1.0e3) - 5.0
            magDev = lc.magDelta
            bandMags, bandlabels = self.get_all_bbs( tSamp )
            for i,label, in enumerate(bandlabels):
                bandMags[i][bandMags[i] < emptyval] += dist_mod + lc.Extinction[label] + magDev
            bbExt = np.array([np.min( band ) for band in bandMags])
            if np.any( bbExt < threshold ) :
                radec_list.append( [ lc.LC_RA, lc.LC_DEC ] )
                bandmags_list.append( bandMags )
            else:
                self.N_trans -=1
        return radec_list, bandmags_list

class galactic_nRec_template( transientTemplate ):
    def __init__( self, tag, trData, PkMData, mag_lim ):
        transientTemplate.__init__( self, tag, trData, PkMData  )
        self.Max_obs_D     =  Maximum_observing_distance(mag_lim, self.peak_mag, self.std_mag)

class galactic_recur_template( transientTemplate ):
    def __init__( self, tag, trData, PkMData, mag_lim ):
        transientTemplate.__init__( self, tag, trData, PkMData)
        self.freqDay = float( trData[9] )
        self.freqDev = float( trData[10] )
        self.Max_obs_D     =  Maximum_observing_distance(mag_lim, self.peak_mag, self.std_mag)
    def get_blueprints( self, c, dtObs, mag_lim ):
        Nr_of_years = 1.	#Recurrent transient => Nr. of transients does not depend on 
                                #observation time
        N_transients = self.get_Ntransients( c, Nr_of_years )
        self.N_transnonsee += N_transients
        for i in range( N_transients ):
            tOutburst = -365.25 * Nr_of_years * np.random.random() 
            magDelta = self.std_mag * 3.0
            peak_M = self.peak_mag + magDelta
            newLC = LC_blueprint( c.sampleCoords(), peak_M, magDelta, tOutburst )
            if canSee_LC( newLC, mag_lim ):
                coords = (newLC.LC_RA, newLC.LC_DEC, newLC.LC_Dist)
                newLC.Extinction = c.Grid.Gal_dust.Sample_extinction(*coords)
                if canSee_LC( newLC, mag_lim ):
                    self.transients.append( newLC )
                    self.N_trans += 1
    def sample_all_LCs( self, obTimes, threshold ):
        radec_list, bandmags_list = [], []
        obT_f = obTimes[-1]
        dts = obT_f - obTimes
        for lc in self.transients:
            radec, mags = self.sample_single_LC( lc, dts, threshold )
            if len(radec) > 0:
                radec_list.append( radec )
                bandmags_list.append( mags )
        return radec_list, bandmags_list
    def sample_single_LC(self, lc, tvals, threshold ):
        # figure out the times for each outburst
        tOuts = []
        t_OutBurst = lc.tExplode
        # step forward in time
        while t_OutBurst < 0.0:
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
            bandMags, bandlabels = self.get_all_bbs( tSamp )
            for band in bandMags:
                band[ band < emptyval ] +=  magDev
            for tO in tOuts[1:]:
                if tO > (-tvals[0] / 8.64e4) - self.deltaT_LC:
                    tSamp = -tvals - tO*8.64e4
                    magDev = np.random.normal( scale = self.std_mag )
                    these_bandMags, bandlabels = self.get_all_bbs( tSamp )
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
        for i,label, in enumerate(bandlabels):
            bandMags[i] += dist_mod - lc.Extinction[label]
        bbExt = np.array([np.min( band ) for band in bandMags])
        if np.any( bbExt < threshold ) :
            return [lc.LC_RA, lc.LC_DEC], bandMags
        else:
            self.N_trans -= 1
            return [],[]

class Mdwarf_template( transientTemplate ):
    def __init__( self, tag, trData, PkMData, mag_lim ):
        transientTemplate.__init__( self, tag, trData, PkMData )
        vals               = np.genfromtxt( pm.MDwarfFile, names=True ).T[tag]
        self.aval, self.bval, self.alf_ac, self.alf_in,\
            self.bet_ac, self.bet_in, self.da_dz = vals[0:7]
        self.Ebins         = 100		#results in 100 bins
        logE_min, logE_max = vals[7:9]
        last_logE          = logE_max - (logE_max - logE_min)/self.Ebins
                             #last_logE is the the start of the last logE-bin
        self.dlogE_act     = np.linspace(logE_min, last_logE, self.Ebins +2 )
        logE_min, logE_max = vals[9:11]
        last_logE          = logE_max - (logE_max - logE_min)/self.Ebins
        self.dlogE_ina     = np.linspace( logE_min, last_logE, self.Ebins +2 )
        qs                 = vals[11:]
        self.Lq            = {'U':qs[0],'B':qs[1],'V':qs[2],'R':qs[3],'I':qs[4],'J':qs[5],'H':qs[6],'K':qs[7],'u':qs[8],'g':qs[9],'r':qs[10],'i':qs[11],'z':qs[12]}
        self.Flare_energy  = float( trData[11] )
        self.Epk           = {}
        self.Max_obs_D     =  Maximum_observing_distance(mag_lim, self.peak_mag, self.std_mag)
    def get_blueprints( self, c, dtObs, mag_lim ):
        Nr_of_years = 1.	#Recurrent transient => Nr. of transients does not depend on 
                                #observation time
        N_transients = self.get_Ntransients( c, Nr_of_years )
        self.N_transnonsee += N_transients
        for i in range( N_transients ):
            t_dummy = -365.25 * Nr_of_years * np.random.random()
            newLC = LC_blueprint( c.sampleCoords(), self.peak_mag, 0.0, t_dummy )
            if canSee_LC( newLC, mag_lim ):
                coords = (newLC.LC_RA, newLC.LC_DEC, newLC.LC_Dist)
                newLC.Extinction = c.Grid.Gal_dust.Sample_extinction(*coords)
                if canSee_LC( newLC, mag_lim ):	#Test WITH dust
                    self.transients.append( newLC )
                    self.N_trans += 1
                """
                else:
                    d_10pc = 0.01
                    dist_ratio = newLC.LC_Dist/d_10pc
                    apparent_mag = newLC.PeakMag + 5.0 * np.log10( dist_ratio ) + newLC.Extinction['B']
                    print newLC.LC_Dist, apparent_mag, newLC.Extinction
                """
    def setUpLC( self, colorchoice ):
        print 'for %s' % self.tag
        self.flux_0 = flux_0
        print self.LCFile + '_%s.hdf5' % colorchoice
        lcdata = h5py.File( self.LCFile + '_%s_rise.hdf5' % colorchoice, 'r' )
        tms = lcdata['times'][:]
        self.broadbands['t'] = tms;
        print tms, "tms"
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
            self.broadbands[band] = fun
            self.Epk[band] = np.max(bb)
        lcdata.close()
    def sample_all_LCs( self, obTimes, threshold ):
        print "nr. of Ms", self.N_trans
        # total time for observation in hours
        tWindow = ( obTimes[-1] - obTimes[0] )/3600.0 
        # lists to hold the output
        radec, bandmags_list = [], []
        for mdwarf in self.transients:
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
                lumlist = self.get_all_lums( nFlares, eBins_m[:-1], obTimes )
                lum_non_outb = np.zeros( len(self.bands) )
                for j, lum in enumerate(lumlist):
                    area = 4.0 * np.pi * np.power( mdwarf.LC_Dist * 1e3 * pc_to_cm, 2.0 ) # cm^2
                    lum = -2.5 * np.log10( lum/area/self.flux_0[ self.bands[j] ] )
                    lumlist[j] = lum
                    #the luminosity (in mags) of the quiescent star
                    lum_non_outb[j] = -2.5 * np.log10( self.Lq[self.bands[j]]/area/self.flux_0[ self.bands[j] ] )
                    lumlist[j] += mdwarf.Extinction[ self.bands[j] ]
                LCminima = np.array( [ np.min(lumin) for lumin in lumlist ])
                if np.any( LCminima < threshold ) and np.any( LCminima < lum_non_outb - pm.mag_resolution ):
                    #import matplotlib.pyplot as plt	#plot a single Mdwarf lightcurve
                    #plt.plot(np.linspace(0,10, len(lumlist[0])), -1.*lumlist[0])
                    #plt.show()
                    radec.append([ mdwarf.LC_RA, mdwarf.LC_DEC ])
                    bandmags_list.append( lumlist )
                else:
                    self.N_trans -=1
        return radec, bandmags_list
    def get_all_lums(self, nFlares, eBins, obT ):
        lumlist = []
        window = obT[-1] - obT[0]
        nBands = len(self.bands)
        lum_matrix = np.zeros( [ nBands, len(obT) ] )
        eF, nF = eBins[nFlares>0], nFlares[nFlares>0]
        for i, n in enumerate(nF):
            logEU = eF[i]
            trise  = np.power(10.0, np.random.normal(0.31*logEU - 7.6, 0.32))
            tdecay = np.power(10.0, np.random.normal(0.51*logEU - 13., 0.22))
            thisE = np.power(10.0, logEU)
            dts = np.random.random( n ) * window
            #generate outbursts that started after t0
            for dt in dts:
                tSamp = obT-dt
                for j, bnd in enumerate(self.bands):
                    if self.Epk[bnd] != 0:
                        new_flux = ( thisE / self.Flare_energy ) * self.broadbands[bnd]( tSamp )
                        lum_matrix[j] += new_flux
            #generate outbursts that started before t0
            a = 1
            for i,x in enumerate(obT):	#How many observations can a burst comprise
                if x > self.deltaT_LC:
                    a = i
                    break
            dts = np.random.random( n ) * window
            dts = dts[ dts < self.deltaT_LC * 24.*60.*60.]    #dts in seconds, deltaT_LC in days
            for dt in dts:
                tSamp = obT[:a]+dt
                for j, bnd in enumerate(self.bands):
                    if self.Epk[bnd] != 0:
                        new_flux = ( thisE / self.Flare_energy ) * self.broadbands[bnd]( tSamp )
                        lum_matrix[j] += new_flux
        for i,row in enumerate(lum_matrix):
            lum_matrix[i] += self.Lq[ self.bands[i] ]
        lumlist = [row for row in lum_matrix]
        return lumlist
                

class xgal_template( transientTemplate  ):
    def __init__( self, tag, trData, PkMData ):
        transientTemplate.__init__( self, tag, trData, PkMData )
        #self.stellarNrm = 4.*4.0e-14
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
    def get_all_bbs( self, tsample , redshift):
        bnds = []
        bandlabels = []
        for bandLabel in self.bands:
            Kcor = np.array(cm.Sample_Kcor(self.Kcorfuns, bandLabel, 
                                           tsample, redshift))
            #Kcor now needs to be "ravelled" because it's a 2D array [[x]] -> [x]
            this_bb = self.broadbands[bandLabel]( tsample * self.scale ) + Kcor.ravel()
            bnds.append(this_bb)
            bandlabels.append(bandLabel)
        return bnds, bandlabels
    def get_Ntransients( self, cell, Nr_of_years ):
        kpc3 = cell.vol
        N_per_y = kpc3 * cell.Grid.Cosmo.Nr_density_xtr( cell.z ) 
        #N_per_y = kpc3 * self.stellarNrm	#Non-cosmological
        N_per_y = N_per_y * Nr_of_years
        if N_per_y > 2.0:
            NTr = int( N_per_y )
        else: NTr = geometricP( N_per_y )
        return NTr
    def get_blueprints( self, c, dtObs, mag_lim ):
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
                if canSee_LC( newLC, mag_lim ):
                    coords = (newLC.LC_RA, newLC.LC_DEC)
                    newLC.Extinction = c.Grid.Xtr_dust.Sample_extinction(*coords)
                    if canSee_LC( newLC, mag_lim ):
                        newLC.Redshift = c.Grid.Cosmo.get_redshift(newLC.LC_Dist / 1.e3)	#convert to Mpc
                        self.transients.append(newLC)
                        self.N_trans += 1
    def sample_all_LCs( self, obTimes, threshold ):
        radec_list, bandmags_list = [], []
        obT_f = obTimes[-1]
        dts = obT_f - obTimes
        self.Kcorfuns = cm.Setup_Kcor(self.Kcormodel, pm.color_system)
        import matplotlib.pyplot as plt
        CM = cm.Cosmology(8.67e6)
        for lc in self.transients: 
            tSamp = -lc.tExplode * 8.64e4 - dts
            dval = lc.LC_Dist	#is in kpc
            dist_mod = 5.0 * np.log10(dval*1.0e3) - 5.0
            magDev = lc.magDelta
            if self.tag == 'SNIa':
                self.Phillips_scale(magDev)
            bandMags, bandlabels = self.get_all_bbs( tSamp, lc.Redshift )
            for i,label, in enumerate(bandlabels):
                bandMags[i][bandMags[i] < emptyval] += dist_mod + lc.Extinction[label] + magDev
            bbExt = np.array([np.min( band ) for band in bandMags])
            if np.any( bbExt < threshold ) :
                radec_list.append( [ lc.LC_RA, lc.LC_DEC ] )
                bandmags_list.append( bandMags )
            else:
                self.N_trans -=1
        """
            if np.min(bandMags) < 21.:
                print CM.get_redshift(lc.LC_Dist/1.e3), lc.tExplode, dist_mod, magDev
                plt.plot(np.linspace(0,80,len(bandMags[0])),bandMags[0])
                plt.plot(np.linspace(0,80,len(bandMags[0])),bandMags[1])
                plt.plot(np.linspace(0,80,len(bandMags[0])),bandMags[2])
                plt.plot(np.linspace(0,80,len(bandMags[0])),bandMags[3])
                plt.plot(np.linspace(0,80,len(bandMags[0])),bandMags[4])
        plt.show()
        """
        return radec_list, bandmags_list

class kilonovaTemplate( transientTemplate ):
    def __init__(self, tMerge, dist, mvRed, mvBlue):
        self.tag = 'kilonova'
        self.bands = list(pm.showbands)
        self.t_merge = tMerge
        if dist != None:
            self.dist = dist
        else: self.dist = sampleLIGODistance( pm.maxLIGODist)
        self.RA = sampleRA( pm.RA_lo, pm.RA_hi )
        self.DEC = sampleDEC( pm.DEC_lo, pm.DEC_hi )
        self.N_trans = 1
        self.N_transnonsee = 1
        self.mvRed = mvRed
        self.mvBlue = mvBlue

    def sample_all_LCs( self, obTimes, threshold ):
        netLC = []
        if self.mvRed != None:
            mred, vred = get_best_match( 'red', *self.mvRed)
            rlcf = 'knRed_m%.4f_vk%.2f_Xlan1e-1.0.h5' % (mred, vred)
            rlcf = 'kilonovaLib/' + rlcf
            redLC = h5py.File(rlcf, 'r' )
            if self.mvBlue == None: 
                netLC.append( redLC['days'][:] )
                for band in self.bands:
                    netLC.append( redLC[band][:] )
                redLC.close()
        if self.mvBlue != None:
            mblu, vblu = get_best_match( 'blue', *self.mvBlue)
            blcf = 'knBlu_m%.4f_vk%.2f_Xlan1e-5.0.h5' % (mblu, vblu)
            blcf = 'kilonovaLib/' + blcf
            bluLC = h5py.File(blcf, 'r')
            if self.mvRed == None:
                netLC.append( bluLC['days'][:] )
                for band in self.bands:
                    netLC.append( bluLC[band][:] )
                bluLC.close()
        if not netLC: # combine kilonova light curves
            times = redLC['days'][:]
            netLC.append(times)
            for band in self.bands:
                m1, m2 = redLC[band][:], bluLC[band][:]
                m12 = -2.5*np.log10( np.power(10.0, -.4*m1) + np.power(10.0, -.4*m2) )
                netLC.append( m12 )
            redLC.close()
            bluLC.close()
        # now that we have the light curve, sample it at the appropriate times
        tSample = (obTimes-self.t_merge)/8.64e4
        dist_mod = 5.0 * np.log10( self.dist/pc_to_cm  ) - 5.0
        for j in range(len(netLC)-1):
            bandLC = netLC[j+1]
            f_band = interp1d( netLC[0][ np.isinf(bandLC) == False],
                               bandLC[np.isinf(bandLC) == False], 
                               bounds_error=False, fill_value=np.inf )
            netLC[j+1] = f_band( tSample ) + dist_mod
        #The following two items should be put into an array like with all other transients
        radec_list = [[self.RA, self.DEC]]
        bandmags_list = [netLC[1:]]	
        return radec_list, bandmags_list



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



def canSee_LC( lightcurve, magLimit ):
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
        for band, x in enumerate(magLimit):
            if apparent_mag[ band ] < x: 
                seen = True
                return seen
    else:
        Bands = list(pm.showbands)
        for band, x in enumerate(magLimit):
            color = Bands[band]
            if apparent_mag[ band ] + lightcurve.Extinction[color] < x: 
                seen = True
                return seen
    return seen


def Maximum_observing_distance(mag_lim, pkmag, devmag):
    exp_term = 0.2*( max( mag_lim - ( pkmag - 3.0*devmag ) ) ) + 1.0
    Max_D = np.power( 10.0, exp_term ) # parsecs
    return Max_D

def Minimum_Flare_E(D, Lq, colorbands):
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
    #########################################################
    #This is not finished. First we need to model the rise/decay time
    #####################################################
    D = D * 1.e3
    Area = 4.0 * np.pi * np.power( D * pc_to_cm, 2.0 )
    #Quiescent fluxes
    Qflux = { band : Lq[band] / Area for band in colorbands}
    #Quiescent magnitudes
    Qmag = { band : flux2mag(Qflux[band], flux_0[band]) for band in colorbands}

    #The difference between the magnitude limit and quiescent magnitude of the M dwarf
    difmag = {band : Qmag[band] - pm.mag_limit[band] for band in colorbands}
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
        Qmag_U = flux2mag( Qflux_U, flux_0['U'] )
        Outburstflux_U = mag2flux( Qmag_U - mindifmag, flux_0['U'])
        L_outburst_U = Outburstflux_U * Area
        print Qmag_U, mindifmag, mag2flux(Outburstflux_U, flux_0['U']), L_outburst_U
        logminE = np.log10(L_outburst_U - Lq['U']) #####################################
        ############The above function now gives an erg/s: it must be changed to erg with the rise/decay time of a flare
        ###################
    return logminE

def mag2flux(mag,flux_0):
    """
    converts magnitude to flux in Jy
    """
    flux = flux_0 * 10**(-0.4*mag)
    return flux

def flux2mag(flux,flux_0):
    """
    converts flux (Jy) to magnitudes
    """
    mag = -2.5*np.log10(flux/flux_0)
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

