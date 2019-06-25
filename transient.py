import numpy as np
import h5py
import coordTransforms as cT
import MilkyWay as MW
import cosmology as cm
import Filterchar as fc
import dust
import glob
from scipy.stats import lognorm
from scipy.interpolate import interp1d


year_to_seconds = 8.64e4 * 365.25
day_to_seconds = 8.64e4
pc_to_cm = 3.08567758e18
emptyval = 1000.0



class TransientSet:
    def __init__ (self, obRun, Framenr):
        self.Framenr = Framenr
        self.grid = obRun.SkyGrids[self.Framenr]
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
        if self.grid.Gal_trans:
            self.populate_galactic(galIdx)
        if self.grid.Xgal_trans:
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
        xgalTrans3 = [ tr for tr in self.transient_templates if tr.galType == 3 ]
        xgalTrans4 = [ tr for tr in self.transient_templates if tr.galType == 4 ]
        cells_xgal = self.grid.cellGrid[galIdx:]
        for xtr in xgalTrans3:
            self.grid.Cosmo.create_Nr_density_funct(xtr.NM, xtr.alpha) 
            for c in cells_xgal:
                xtr.get_blueprints( c, self.deltaT_obs )
        for xtr in xgalTrans4:
            for c in cells_xgal:
                xtr.get_blueprints( c, self.deltaT_obs )
    def inject_kilonova( self, tMerge, dist, mvRed, mvBlue, obRun, Framenr ):
        self.transient_templates.append( kilonovaTemplate( tMerge, dist, mvRed, mvBlue, obRun, Framenr ))
        #print "translegend", obRun.Trlegend
        if obRun.Trlegend == {}:
            obRun.Trlegend['kilonova'] = 0
        else:
            max_val = max(obRun.Trlegend, key=obRun.Trlegend.get)
            #print max_val
            obRun.Trlegend['kilonova'] = obRun.Trlegend[max_val] + 1
            
def makeTemplate( tag, obRun ):
    transientsData = getFileLine( obRun.transientFile, tag )
    Type = int( transientsData[0] )
    if Type == 0 :
        temp = galactic_nRec_template( tag, transientsData, obRun )
    elif Type == 1:
        if tag == 'UGem':
            temp = UGem_template( tag, transientsData, obRun )
        elif tag == 'SUUMa':
            temp = SUUMa_template( tag, transientsData, obRun )
        elif tag == 'ZCam':
            temp = ZCam_template( tag, transientsData, obRun )
        else:
            temp = galactic_recur_template( tag, transientsData, obRun )
    elif Type == 2:
        temp = Mdwarf_template( tag, transientsData, obRun)
    elif Type == 3:
        temp = xgal_template( tag, transientsData, obRun)
    elif Type == 4:
        temp = xgal_no_Kcor_template( tag, transientsData, obRun)
    return temp

class transientTemplate:
    def __init__( self, tag, transientsData, obRun ):
        self.tag         = tag
        self.galType     = int( transientsData[0] )
        self.galkey      = [int(tD) for tD in transientsData[1:5]]
        self.bands       = obRun.bands
        self.mag_lim     = obRun.threshold
        self.broadbands  = {} 
        self.LCFile      = 'LightCurveFiles/' + transientsData[5]
        self.stellarNrm  = float( transientsData[6] )/MW.rho_stellar_sun
        self.scaleH      = float( transientsData[7] )
        self.N_trans     = 0    # how many transients each template has
        self.transients  = []   # list of light curves blueprints here
        self.N_transnonsee = 0  #Nr of generated transients, regardless of visibility
        #self.peak_mag, self.std_mag = PeakMagData( obRun.Peak_mag_f, tag, self.bands )
                                      #Mean and standard deviation of peak magnitude
        self.peak_mag_R    = float( transientsData[10] )
        self.std_mag_R     = float( transientsData[11] )
        self.radec_list    = []
        self.bandmags_list = []
        self.bandbands_list= []
        self.Qmag_list     = []
        self.obT_list      = []
        if self.tag in ['UGem', 'SUUMa', 'ZCam']:
            self.deltaT_LC = transientsData[8].split(',')
        else:
            self.deltaT_LC = float( transientsData[8] )
    def setUpLC( self, colorchoice ):
        print 'for %s' % self.tag
        lcdata_up = h5py.File( self.LCFile + '_UBVRI.hdf5', 'r' )
        lcdata_lo = h5py.File( self.LCFile + '_%s.hdf5' % colorchoice, 'r' )
        tms = lcdata_up['times'][:]    #times(UBVRI) = times(ugriz)
        self.broadbands['t'] = tms
        self.peak_mags = {}   #Empty dict to be filled in later
        for band in self.bands:
            print "%s band" % band
            if band.islower():
                lcdata = lcdata_lo
            else:
                lcdata = lcdata_up
            if not band in lcdata.keys():
                print "no data for %s band" % band
                bb = emptyval * np.ones( len(tms) )
            else:
                scale = self.peak_mag_R - min(lcdata_up['R'][:])
                bb = lcdata[band][:] + scale
            ff = interp1d( tms, bb, bounds_error = False, fill_value = emptyval )
            def fun( times, fcn = ff ):
                return fcn(times)
            self.broadbands[band] = fun
            #Now we need to determine the brightness maximum in each color band
            #This is used in get_blueprints to determine whether the transient is
            # visible or not
            self.peak_mags[band] = min(lcdata[band][:])
        lcdata_lo.close()
        lcdata_up.close()
        self.Max_obs_D     =  Maximum_observing_distance(self.mag_lim, self.peak_mags, self.std_mag_R)
    def get_all_bbs( self, tsample, obbands ):
        bnds = []
        for i,bandLabel in enumerate(obbands):
            this_bb = self.broadbands[bandLabel]( tsample[i] )
            bnds.append(this_bb)
        return np.array(bnds)
    def get_Ntransients( self, cell, Nr_of_years ):
        if abs(self.scaleH - MW.MW_Hthin) < 0.001:
            N_per_y = self.stellarNrm * cell.vol * np.sum( np.array(cell.rho)[ np.where( self.galkey )] )
        else:
            this_rho = MW.get_MW_dens( [cell.raMid, cell.decMid, cell.DMid],\
                                           piecewise=True, Hthin = self.scaleH )
            N_per_y = self.stellarNrm * cell.vol * np.sum( np.array(this_rho)[ np.where( self.galkey )] )
        N_per_y = N_per_y * Nr_of_years
        #if N_per_y > 2.0:
        #    NTr = int( N_per_y )
        #else: NTr = geometricP( N_per_y )
        return geometricP( N_per_y )   
    def get_blueprints( self, c, dtObs ):
        Nr_of_years = (self.deltaT_LC + dtObs) / 365.25   #At first observation the first LC
                                                          #should have just ended
        N_transients = self.get_Ntransients( c, Nr_of_years )
        self.N_transnonsee += N_transients
        for i in range( N_transients ):
            coords = c.sampleCoords()
            if coords != []:
                tOutburst = -365.25 * Nr_of_years * np.random.random()
                if tOutburst + self.deltaT_LC + dtObs > 0.0:
                    magDelta = np.random.normal( scale = self.std_mag_R )
                    peak_M = {band: self.peak_mags[band] + magDelta for band in self.bands}
                    newLC = LC_blueprint( coords, peak_M, magDelta, tOutburst )
                    if self.canSee_LC( newLC ):
                        _Coords = (newLC.LC_RA, newLC.LC_DEC, newLC.LC_Dist)
                        newLC.Extinction = c.Grid.Gal_dust.Sample_extinction(*_Coords)
                        if self.canSee_LC( newLC ):
                            self.transients.append(newLC)
                            self.N_trans += 1
    def sample_all_LCs( self, obTimes, obBands, threshold ):
        for lc in self.transients: 
            vis_frames = np.array(lc.visibleframes)
            obT = np.concatenate(obTimes[vis_frames])
            obB = np.concatenate(obBands[vis_frames])
            obT_f = max(obT)
            dts = obT_f - obT
            visible = False
            tSamp = -lc.tExplode * 8.64e4 - dts
            dval = lc.LC_Dist	#is in kpc
            dist_mod = 5.0 * np.log10(dval*1.0e3) - 5.0
            magDev = lc.magDelta
            bandMags = self.get_all_bbs( tSamp, obB )
            for i,label, in enumerate(obB):
                bandMags[i] += dist_mod + lc.Extinction[label] + magDev
                if bandMags[i] < threshold[label]: visible = True
            if visible:
                self.radec_list.append( [ lc.LC_RA, lc.LC_DEC ] )
                self.bandmags_list.append( bandMags )
                self.bandbands_list.append( obB )
                self.Qmag_list.append(emptyval * np.ones(len(bandMags)))
                self.obT_list.append( obT )
            else:
                self.N_trans -=1
        self.bandmags_list = np.array(self.bandmags_list)
    def take_data(self, colorchoice, colorsys):
        """
        Converts Vega magnitudes to AB magnitudes or vice versa if needed
        Returns the observational data to observ.py
        """
        if colorchoice != 'UBVRI':
            AB_to_Vega = fc.AB_to_Vega[colorchoice]
        Vega_to_AB = fc.Vega_to_AB['UBVRI']
        for i,mags in enumerate(self.bandmags_list):
            for j,band in enumerate(self.bandbands_list[i]):
                if colorsys == 'AB':
                    if band.isupper():
                        self.bandmags_list[i][j] = mags[j] * Vega_to_AB[band]
                if colorsys == 'Vega':
                    if band.islower():
                        self.bandmags_list[i][j] = mags[j] * AB_to_Vega[band]
        return self.radec_list, self.bandmags_list, self.Qmag_list, self.bandbands_list, self.obT_list
        """
        bands = self.bandbands_list[0]
        if np.any(convert_toAB):
            C = [convert_toAB[i] * Vega_to_AB[b] for i,b in enumerate(bands)]
        elif np.any(convert_toVega):
            C = [convert_toVega[i] * AB_to_Vega[b] for i,b in enumerate(bands)]
        else:
            C = np.zeros(len(bands))
        Conv = np.tile(C, (len(self.bandmags_list),1))
        self.bandmags_list += Conv
        return self.radec_list, self.bandmags_list, self.bandbands_list
        """
    def canSee_LC(self, lightcurve):
        """
        A crude way to estimate whether the light curve should still be visible
        If the Extinction is not zero, it has been adjusted to a real value.
        That lets us calculate the visibility of the transient WITH extinction.
        """
        seen = False
        d_10pc = 0.01
        dist_ratio = lightcurve.LC_Dist/d_10pc
        apparent_mag = {band: lightcurve.PeakMag[band]
                        + 5.0 * np.log10( dist_ratio ) for band in self.bands}
        if lightcurve.Extinction == 0:
            for band in self.bands:
                if apparent_mag[band] < self.mag_lim[band]: 
                    seen = True
                    return seen
        else:
            for band in self.bands:
                if apparent_mag[band] + lightcurve.Extinction[band] < self.mag_lim[band]: 
                    seen = True
                    return seen
        return seen

class galactic_nRec_template( transientTemplate ):
    def __init__( self, tag, trData, obRun ):
        transientTemplate.__init__( self, tag, trData, obRun  )

class galactic_recur_template( transientTemplate ):
    def __init__( self, tag, trData, obRun ):
        transientTemplate.__init__( self, tag, trData, obRun)
        freq_normal = trData[9].split(',')
        self.P_c_s    = float( freq_normal[0].replace("(", "") )
        self.P_c_loc  = float( freq_normal[1] )
        self.P_c_scale= float( freq_normal[2].replace(")", "") )
        #self.P_c_mu    = float( trData[9] )#float(freq_normal[0].replace("(", ""))
        #self.P_c_sigma = float( trData[10] )#float(freq_normal[1].replace(")", ""))
    def get_blueprints( self, c, dtObs ):
        Nr_of_years = 1.	#Recurrent transient => Nr. of transients does not depend on 
                                #observation time
        N_transients = self.get_Ntransients( c, Nr_of_years )
        self.N_transnonsee += N_transients
        for i in range( N_transients ):
            coords = c.sampleCoords()
            if coords != []:
                tOutburst = -365.25 * Nr_of_years * np.random.random()
                magDelta = self.std_mag_R * 3.0
                peak_M = {band: self.peak_mags[band] + magDelta for band in self.bands}
                newLC = LC_blueprint( coords, peak_M, magDelta, tOutburst )
                if self.canSee_LC( newLC):
                    _Coords = (newLC.LC_RA, newLC.LC_DEC, newLC.LC_Dist)
                    newLC.Extinction = c.Grid.Gal_dust.Sample_extinction(*_Coords)
                    if self.canSee_LC( newLC ):
                        self.transients.append( newLC )
                        self.N_trans += 1
    def sample_all_LCs( self, obTimes, obBands, threshold ):
        for lc in self.transients:
            vis_frames = np.array(lc.visibleframes)
            obT = np.concatenate(obTimes[vis_frames])
            #print "obT", type(obT), obT
            obB = np.concatenate(obBands[vis_frames])
            obT_f = max(obT)
            dts = obT_f - obT
            radec, mags, Qmags = self.sample_single_LC( lc, dts, obB, threshold )
            if len(radec) > 0:
                self.radec_list.append( radec )
                self.bandmags_list.append( mags )
                self.Qmag_list.append(Qmags)
                self.bandbands_list.append( obB )
                self.obT_list.append(obT)
        self.bandmags_list = np.array(self.bandmags_list)
    def sample_single_LC(self, lc, tvals, obBands, threshold ):
        #Sample outburst frequency
        #self.P_c       = np.random.lognormal(self.P_c_mu, self.P_c_sigma, 1)
        self.P_c = lognorm.rvs(self.P_c_s, self.P_c_loc, self.P_c_scale)
        # figure out the times for each outburst
        timeOuts = []
        t_OutBurst = lc.tExplode
        # step forward in time
        while t_OutBurst < self.deltaT_LC + (max(tvals) / 8.64e4):
            timeOuts.append( t_OutBurst )
            t_OutBurst += self.P_c#np.random.normal( loc = self.freq_normalDay, scale = self.freq_normalDev)
        # step backward in time
        t_OutBurst = lc.tExplode -self.P_c# np.random.normal( loc = self.freq_normalDay, scale = self.freq_normalDev)
        while t_OutBurst + self.deltaT_LC > (-min(tvals) / 8.64e4): 
            timeOuts.insert(0, t_OutBurst)
            t_OutBurst -= self.P_c#np.random.normal( loc = self.freq_normalDay, scale = self.freq_normalDev)
        # construct the aggregate light curve
        if len( timeOuts ) > 0:
            timeOuts = timeOuts[::-1]	#flip the array to start with the closest outburst
            tO = timeOuts[0]
            tSamp = -tvals - tO*8.64e4
            # deviation for this occurrence of the transient
            magDev = np.random.normal( scale = self.std_mag_R )
            bandMags = self.get_all_bbs( tSamp, obBands )
            for band in bandMags:
                band[ band < emptyval ] +=  magDev
            for tO in timeOuts[1:]:
                if tO > (-tvals[0] / 8.64e4) - self.deltaT_LC:
                    tSamp = -tvals - tO*8.64e4
                    magDev = np.random.normal( scale = self.std_mag_R )
                    these_bandMags = self.get_all_bbs( tSamp, obBands )
                    for bi, band in enumerate( these_bandMags ):    #Check for overlap with previous outbursts and add if necessary
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
        Qmags = emptyval * np.ones(len(bandMags))
        for i,label, in enumerate(obBands):
            bandMags[i] += dist_mod + lc.Extinction[label]
            if bandMags[i] < threshold[label]: visible = True
        if visible :
            return [lc.LC_RA, lc.LC_DEC], bandMags, Qmags
        else:
            self.N_trans -= 1
            return [],[],[]

class SUUMa_template( galactic_recur_template ):
    def __init__( self, tag, trData, obRun ):
        galactic_recur_template.__init__( self, tag, trData, obRun  )
        self.type              = 'SUUMa'
        self.LCFile            = 'LightCurveFiles/Dwarfnovae/' + trData[5]
        self.broadbands_normal = {}
        self.broadbands_super  = {}
        self.quiescentmag      = {}
        self.BP_Q_mag          = {}
    def setUpLC( self, colorchoice ):
        print 'for %s' % self.tag
        self.BP_outb_amp = {}
        self.BP_D_no = float( self.deltaT_LC[0].replace("(", "") )
        self.BP_D_so = float( self.deltaT_LC[1].replace(")", "") )
        self.BP_t_s  = {}
        self.BP_t_e  = {}
        for TYPE in ['normal', 'super']:
            lcdata_up = h5py.File( self.LCFile + '_UBVRI_%s.hdf5' % (TYPE), 'r' )
            lcdata_lo = h5py.File( self.LCFile + '_%s_%s.hdf5' % (colorchoice, TYPE), 'r' )
            self.BP_outb_amp[TYPE], self.BP_Q_mag[TYPE] = self.Get_outb_amp(lcdata_up, lcdata_lo, np.append('V',self.bands))
            tms = lcdata_up['times'][:]    #times(UBVRI) = times(ugriz)
            self.BP_t_s[TYPE] = tms[0]
            self.BP_t_e[TYPE] = tms[-1]
            if TYPE == 'normal':
                self.broadbands_normal['t'] = tms
            if TYPE == 'super':
                self.broadbands_super['t'] = tms
            self.peak_mags = {}   #Empty dict to be filled in later
            for band in np.append('V',self.bands):
                print "%s band" % band
                if band.islower():
                    lcdata = lcdata_lo
                else:
                    lcdata = lcdata_up
                if not band in lcdata.keys():
                    print "no data for %s band" % band
                    bb = emptyval * np.ones( len(tms) )
                else:
                    #scale = self.peak_mag_R - min(lcdata_up['R'][:])
                    bb = lcdata[band][:] - max(lcdata[band][:]) #Only record relative change
                ff = interp1d( tms, bb, bounds_error = False, fill_value = emptyval )
               
                this_Q_mag = max(lcdata[band][:]) - self.BP_Q_mag[TYPE]
                if band not in self.quiescentmag.keys():
                    self.quiescentmag[band] = this_Q_mag
                else:
                    self.quiescentmag[band] = max(self.quiescentmag[band], this_Q_mag)
                def fun( times, fcn = ff ):
                    return fcn(times)
                if TYPE == 'normal':
                    self.broadbands_normal[band] = fun
                elif TYPE == 'super':
                    self.broadbands_super[band] = fun
                #Now we need to determine the brightness maximum in each color band
                #This is used in get_blueprints to determine whether the transient is
                # visible or not
                self.peak_mags[band] = min(lcdata[band][:])
                
            lcdata_lo.close()
            lcdata_up.close()
        self.Max_obs_D     =  Maximum_observing_distance(self.mag_lim, self.peak_mags, self.std_mag_R)
    def get_all_normal_bbs( self, tsample, obbands, DN ):
        bnds = []
        for i,bandLabel in enumerate(obbands):
            outb_amp_old = self.BP_outb_amp['normal']
            scale = ( (DN.A_no + outb_amp_old[bandLabel] - outb_amp_old['V']) /
                      (outb_amp_old[bandLabel]) )
            this_bb = ( self.broadbands_normal[bandLabel]( tsample[i] ) * scale
                       + DN.Q_mag[bandLabel] )
            if this_bb > DN.Q_mag[bandLabel]:
                this_bb = DN.Q_mag[bandLabel]
            bnds.append(this_bb)
        return np.array(bnds)
    def get_all_super_bbs( self, tsample, obbands, DN ):
        bnds = []
        for i,bandLabel in enumerate(obbands):
            outb_amp_old = self.BP_outb_amp['super']
            scale = ( (DN.A_no + outb_amp_old[bandLabel] - outb_amp_old['V']) /
                      (outb_amp_old[bandLabel]) )
            this_bb = ( self.broadbands_super[bandLabel]( tsample[i] ) * scale
                       + DN.Q_mag[bandLabel] )
            if this_bb > DN.Q_mag[bandLabel]:
                this_bb = DN.Q_mag[bandLabel]
            bnds.append(this_bb)
        return np.array(bnds)
    """temp:::::   (entire function is temporary)
    def sample_all_LCs( self, obTimes, obBands, threshold ):
        #Temppppp
        #print obTimes, "obtimes"
        #obBands = np.array(['V' for i in range(3650)])
        #obBands = [['u', 'g', 'r', 'i'] for i in range(3650/4)]
        #obBands = np.array(obBands).ravel()
        #obBands = np.append(obBands, ['g', 'r'])
        #obTimes = np.linspace(0,365. * 24. * 60. * 60., 3650)

        for lc in self.transients:
            vis_frames = np.array(lc.visibleframes)
            obT = np.concatenate(obTimes[vis_frames])
            #print "obT", type(obT), obT
            obB = np.concatenate(obBands[vis_frames])
            obT_f = max(obT)
            dts = obT_f - obT
            radec, mags, Qmags = self.sample_single_LC( lc, dts, obB, threshold )
            if len(radec) > 0:
                self.radec_list.append( radec )
                self.bandmags_list.append( mags )
                self.Qmag_list.append(Qmags)
                self.bandbands_list.append( obB )
                self.obT_list.append(obT)
        self.bandmags_list = np.array(self.bandmags_list)
    """
    def sample_single_LC(self, lc, tvals, obBands, threshold ):
        DN = DwarfNova(self.type, self.P_c_s, self.P_c_loc, self.P_c_scale,
                       self.quiescentmag, self.bands)
        max2 = DN.Get_max_lookbacktime()#DN.maxoutblen + DN.maxrecurt
        t_OutBurst = -max2 * np.random.random_sample() - max2   #the time of the first outburst
        if DN.nr_no_per_so > 1.:
            init_outbtype = np.random.choice([0,1], p=[1 - 1./DN.nr_no_per_so, 
                                                       1./DN.nr_no_per_so])
        else:
            init_outbtype = np.random.choice([0,1], p=[1./DN.nr_so_per_no, 
                                                       1 - 1./DN.nr_so_per_no])

        outb_interval = {0 : self.Get_outbinterval( DN, DN.minrecurt, 0), 
                         1 : self.Get_outbinterval( DN, DN.minrecurt, 1)}
        
        
        t_OutBurst += outb_interval[ init_outbtype ]
        timeOuts = [t_OutBurst]
        typeOuts = [init_outbtype]
        

                    
        if DN.P_sc > 2 * DN.P_c:
            nr_no_outb = 0
            if init_outbtype == 0:
                nr_no_outb = 1
             # step forward in time
            while t_OutBurst < max(tvals) / 8.64e4:
                outbtype = self.Get_outbtype(nr_no_outb, DN.nr_no_per_so)
                if outbtype == 0:   nr_no_outb += 1
                elif outbtype == 1: nr_no_outb = 0
                t_OutBurst += outb_interval[ outbtype ]
                typeOuts.append(outbtype)
                timeOuts.append( t_OutBurst )
        elif DN.P_c < 2 * DN.P_sc:
            nr_so_outb = 0
            if init_outbtype == 1:
                nr_so_outb = 1
             # step forward in time
            while t_OutBurst < 1.1 * max(tvals) / 8.64e4:
                outbtype_so = self.Get_outbtype(nr_so_outb, DN.nr_so_per_no)
                if outbtype_so == 1: outbtype = 0
                elif outbtype_so == 0: 
                    outbtype = 1
                    nr_so_outb = 0
                t_OutBurst += outb_interval[ outbtype ]
                typeOuts.append(outbtype)
                timeOuts.append( t_OutBurst )
        else:   #every time there's either a normal or super outburst
            while t_OutBurst < 1.1 * max(tvals) / 8.64e4:
                if typeOuts[-1] == 0:
                    outbtype = 1
                else:
                    outbtype = 0
                if outbtype == 0:
                    if DN.D_no > DN.P_c:
                        t_OutBurst += DN.D_no
                    else:
                        t_OutBurst += DN.P_c
                elif outbtype == 1:
                    if DN.D_so > DN.P_sc:
                        t_OutBurst += DN.D_so
                    else:
                        t_OutBurst += DN.P_sc     
                typeOuts.append(outbtype)
                timeOuts.append( t_OutBurst )
                        
                        
        # construct the aggregate light curve
        timeOuts = -np.array(timeOuts)

        bandMags = []
        Q_bandMags = []
        for i,band in enumerate(obBands):
            this_bb = DN.Q_mag[band]
            bandMags.append(this_bb)
            Q_bandMags.append(this_bb)
        bandMags = np.array(bandMags)
        Q_bandMags = np.array(Q_bandMags)
    
        for i, tO in enumerate(timeOuts):
            if typeOuts[i] == 0:
                tscale = DN.D_no / self.BP_D_no
                tSamp = np.array(-tvals - tO*8.64e4) / tscale #Scale duration
                Bool1 = tSamp < self.BP_t_e['normal'] #self.BP_D_no *8.64e4
                Bool2 = tSamp > self.BP_t_s['normal']
            if typeOuts[i] == 1:
                tscale =  DN.D_so / self.BP_D_so
                tSamp = np.array(-tvals - tO*8.64e4) / tscale #Scale duration
                Bool1 = tSamp < self.BP_t_e['super'] # self.BP_D_so *8.64e4
                Bool2 = tSamp > self.BP_t_s['super']
            if any(Bool1 & Bool2):   #i.e. the outburst is *in* the observations.
                if typeOuts[i] == 0:
                    these_bandMags = np.array( self.get_all_normal_bbs( tSamp, 
                                                             obBands, DN ) )
                if typeOuts[i] == 1:
                    these_bandMags = np.array( self.get_all_super_bbs( tSamp, 
                                                             obBands, DN ) )
                    
                i5 = (these_bandMags != Q_bandMags)
                bandMags[i5] = np.min([these_bandMags[i5],bandMags[i5]], axis=0)
        """temp:::::
        print "P_c, P_sc", DN.P_c, DN.P_sc
        print "A_no, A_so", DN.A_no, DN.A_so
        print "D_no, D_so", DN.D_no, DN.D_so
        print "peakmags", self.peak_mags['V']
        print self.BP_outb_amp
        print timeOuts
        print typeOuts
        import matplotlib.pyplot as plt
        with open('temp', 'w') as file:
            for A in bandMags:
                file.write(str(A))
                file.write('\n')
        file.close()
        BoolB = obBands == 'u'
        BoolV = obBands == 'g'
        BoolR = obBands == 'r'
        BoolI = obBands == 'i'
        tvals = np.array(tvals)
        plt.plot(365 - tvals[BoolB] / 8.64e4, bandMags[BoolB], color = '#0028ff')
        plt.plot(365 - tvals[BoolV] / 8.64e4, bandMags[BoolV], color = '#a6ff00')
        plt.plot(365 - tvals[BoolR] / 8.64e4, bandMags[BoolR], color = '#ff0000')
        plt.plot(365 - tvals[BoolI] / 8.64e4, bandMags[BoolI], color = '#c20000')
        plt.gca().invert_yaxis()
        plt.show()
        """
        
        # apply distance modulus
        dval = lc.LC_Dist
        dist_mod = 5.0 * np.log10( dval*1.0e3 ) - 5.0
        visible = False
        for i,label, in enumerate(obBands):
            bandMags[i] += dist_mod + lc.Extinction[label]
            Q_bandMags[i] += dist_mod + lc.Extinction[label]
            if bandMags[i] < threshold[label]: visible = True
        
        #print "bandmags", self.tag, type(bandMags)
        if visible :
            return [lc.LC_RA, lc.LC_DEC], bandMags, Q_bandMags
        else:
            self.N_trans -= 1
            return [],[], []
                        
    def Get_outbtype( self, nr_no_outb, nr_no_per_so):
        """
        What type of outburst do we expect? if e.g. nr_no_per_so=3.2
         then we expect to have 3 normal outbursts and then a super
         outburst with a probability of 0.2. If there have already 
         been 4 normal outbursts in a row, we get a super outburst.
        """
        if nr_no_outb > nr_no_per_so:
            return 1
        else:
            if nr_no_per_so - nr_no_outb > 1:
                return 0
            else:
                p = np.random.random_sample()
                if p > nr_no_per_so - nr_no_outb:
                    return 0
                else:
                    return 1
    def Get_outbinterval( self, DN, outb_interval, outbtype, D_ss = 0):
        """
        Calculates the interval in between two outbursts.
        This is generally longer for super outbursts.
        returns: the interval between two outbursts.
        """
        if DN.P_c < DN.P_sc:
            if outbtype == 1:
                new_outb_interval = outb_interval + (DN.D_so - DN.D_no)  #Add duration of superoutburst
                if DN.D_so > new_outb_interval:
                    return DN.D_so
                else:
                    return new_outb_interval
            elif outbtype == 2:
                return outb_interval + (D_ss - DN.D_no)
            elif outbtype == 0 and DN.D_no > outb_interval:
                return DN.D_no
            else:
                return outb_interval
        else:
            if outbtype == 0:
                new_outb_interval = outb_interval + (DN.D_no - DN.D_so)  #Add duration of superoutburst
                if DN.D_no > new_outb_interval:
                    return DN.D_no
                else:
                    return new_outb_interval
            elif outbtype == 1 and DN.D_so > outb_interval:
                return DN.D_so
            elif outbtype == 2:
                return outb_interval + (D_ss - DN.D_so)
            else:
                return outb_interval

    def Get_outb_amp( self, lcdata_up, lcdata_lo, bands):
        """
        A little function that finds the outburst amplitude for every band
        It also calculates the quiescence magnitude in the V-band of the light
         curve data.
        returns: the outburst amplitude in every band
                  and the quiescence mag in V
        """
        if 'V' not in lcdata_up.keys(): 
            raise TypeError("There is no V-band data in your UBVRI light curve file for %s" % (self.tag))
        Q_V_mag = max(lcdata_up['V'][:])
        amp = {}
        for band in bands:
            if band.islower():
                lcdata = lcdata_lo
            else:
                lcdata = lcdata_up
            Q_mag = max(lcdata[band][:])
            amp[band] = -min(lcdata[band][:] - Q_mag)
        return amp, Q_V_mag
        
        
class UGem_template( SUUMa_template ):
    """
    The U Gem template is the same as for SU UMa.
    Just like with SU UMa there are two different types of outbursts:
      - The normal one
      - A longer one, but that isn't denoted as "super" in the outburst 
         because it is not brighter than its 'normal' counterpart.
    """
    def __init__( self, tag, trData, obRun ):
        SUUMa_template.__init__( self, tag, trData, obRun  )
        self.type = 'UGem'
        
class ZCam_template( SUUMa_template ):
    def __init__( self, tag, trData, obRun ):
        SUUMa_template.__init__( self, tag, trData, obRun  )
        self.type = 'ZCam'
        self.broadbands_rise   = {}   #for stillstand
        self.broadbands_decay  = {}   #for stillstand
    def setUpLC( self, colorchoice ):
        print 'for %s' % self.tag
        self.BP_outb_amp = {}
        self.BP_D_no = float( self.deltaT_LC[0].replace("(", "") )
        self.BP_D_so = float( self.deltaT_LC[1].replace(")", "") )
        self.BP_t_s  = {}
        self.BP_t_e  = {}
        for TYPE in ['normal', 'super', 'rise', 'decay']:
            print "file", self.LCFile + '_UBVRI_%s.hdf5' % (TYPE), 'r'
            lcdata_up = h5py.File( self.LCFile + '_UBVRI_%s.hdf5' % (TYPE), 'r' )
            lcdata_lo = h5py.File( self.LCFile + '_%s_%s.hdf5' % (colorchoice, TYPE), 'r' )
            self.BP_outb_amp[TYPE], self.BP_Q_mag[TYPE] = self.Get_outb_amp(lcdata_up, lcdata_lo, np.append('V',self.bands))
            tms = lcdata_up['times'][:]    #times(UBVRI) = times(ugriz)
            self.BP_t_s[TYPE] = tms[0]
            self.BP_t_e[TYPE] = tms[-1]
            self.BP_A_ss = {band: 0 for band in np.append('V',self.bands)}
            if TYPE == 'normal':
                self.broadbands_normal['t'] = tms
            elif TYPE == 'super':
                self.broadbands_super['t'] = tms
            elif TYPE == 'rise':
                self.broadbands_rise['t'] = tms
            elif TYPE == 'decay':
                self.broadbands_decay['t'] = tms
            self.peak_mags = {}   #Empty dict to be filled in later
            for band in np.append('V',self.bands):
                print "%s band" % band
                if band.islower():
                    lcdata = lcdata_lo
                else:
                    lcdata = lcdata_up
                if not band in lcdata.keys():
                    print "no data for %s band" % band
                    bb = emptyval * np.ones( len(tms) )
                else:
                    #scale = self.peak_mag_R - min(lcdata_up['R'][:])
                    bb = lcdata[band][:] - max(lcdata[band][:]) #Only record relative change
                ff = interp1d( tms, bb, bounds_error = False, fill_value = emptyval )
                
                this_Q_mag = max(lcdata[band][:]) - self.BP_Q_mag[TYPE]
                if band not in self.quiescentmag.keys():
                    self.quiescentmag[band] = this_Q_mag
                else:
                    self.quiescentmag[band] = max(self.quiescentmag[band], this_Q_mag)
                def fun( times, fcn = ff ):
                    return fcn(times)
                if TYPE == 'normal':
                    self.broadbands_normal[band] = fun
                elif TYPE == 'super':
                    self.broadbands_super[band] = fun
                elif TYPE == 'rise':
                    self.broadbands_rise[band] = fun
                    self.BP_A_ss[band] = min(self.BP_A_ss[band], lcdata[band][-1] - max(lcdata[band][:]))
                    #print lcdata[band][-1], "lcdata"
                elif TYPE == 'decay':
                    self.broadbands_decay[band] = fun
                    self.BP_A_ss[band] = min(self.BP_A_ss[band], lcdata[band][0] - max(lcdata[band][:]))
                #Now we need to determine the brightness maximum in each color band
                #This is used in get_blueprints to determine whether the transient is
                # visible or not
                self.peak_mags[band] = min(lcdata[band][:])
                
            lcdata_lo.close()
            lcdata_up.close()
        self.BP_d_rise = self.BP_t_e['rise'] - self.BP_t_s['rise'] #in s
        self.BP_d_decay = self.BP_t_e['decay'] - self.BP_t_s['decay'] #in s
        self.Max_obs_D     =  Maximum_observing_distance(self.mag_lim, self.peak_mags, self.std_mag_R)
    def get_all_standstill_bbs( self, tsample, obbands, DN, D_ss, tO):
        #print "hallo", D_ss, DN.D_ss
        bnds = []
        rise =  ( self.broadbands_rise['V']( tsample ) * 
                       DN.A_no / self.BP_outb_amp['normal']['V']  
                       + DN.Q_mag['V'] )
        lowvalue = ( emptyval * DN.A_no / self.BP_outb_amp['normal']['V']  
                       + DN.Q_mag['V'] )
        if np.any(rise < lowvalue):
            endrise = tsample[rise < lowvalue][-1]
        else:
            endrise = tO + self.BP_t_e['rise']
            d_till_decay_i = self.BP_D_no / DN.D_no * D_ss
            bool1 = np.any(tsample > endrise)
            bool2 = np.any(tsample < d_till_decay_i * 8.64e4)
            bool3 = np.any(tsample > d_till_decay_i * 8.64e4)
            if not bool1 and bool2 or not bool3:
                return emptyval * np.ones(len(tsample))
        #d_till_decay is the duration of the standstill up till the decay
        # phase in the blueprint ``time-frame''
        d_till_decay = self.BP_D_no / DN.D_no * D_ss + (endrise / 8.64e4)
        for i,bandLabel in enumerate(obbands):
            outb_amp_old = self.BP_outb_amp['normal']
            
            scale = ( (DN.A_no + outb_amp_old[bandLabel] - outb_amp_old['V']) /
                      (outb_amp_old[bandLabel]) )
            this_bb = ( self.broadbands_rise[bandLabel]( tsample[i] ) * scale
                       + DN.Q_mag[bandLabel] )
            if tsample[i] > endrise and tsample[i] < d_till_decay * 8.64e4:
                this_bb = self.BP_A_ss[bandLabel] * scale + DN.Q_mag[bandLabel]
            if tsample[i] > d_till_decay * 8.64e4 :
                time = tsample[i] - d_till_decay * 8.64e4
                this_bb = ( self.broadbands_decay[bandLabel]( time ) * scale
                           + DN.Q_mag[bandLabel] )
            if this_bb > DN.Q_mag[bandLabel]:
                this_bb = DN.Q_mag[bandLabel]
            bnds.append(this_bb)
        #print "SSmag", self.BP_A_ss['I'] * scale + DN.Q_mag['I']
        #print obbands[np.array(bnds) < -20.]
        return np.array(bnds)
    def sample_single_LC(self, lc, tvals, obBands, threshold ):
        """
        #obBands = np.array(['V' for i in range(3650)])
        obBands = [['U','B', 'V', 'R', 'I'] for i in range(3650/5)]
        obBands = np.array(obBands).ravel()
        #obBands = np.append(obBands, ['R', 'I'])
        obTimes = np.linspace(0,365. * 24. * 60. * 60., 3650)
        tvals = max(obTimes) - obTimes
        """
        """
        """
        
        
        DN = DwarfNova(self.type, self.P_c_s, self.P_c_loc, self.P_c_scale,
                       self.quiescentmag, self.bands)
        max2 = DN.Get_max_lookbacktime(D_rise = self.BP_d_rise / 8.64e4, 
                                       D_decay = self.BP_d_decay / 8.64e4)
        t_OutBurst = -max2 * np.random.random_sample() - max2   #the time of the first outburst
        init_outbtype = 2  #a standstill
        DN.D_rise  = (self.BP_d_rise / 8.64e4) * DN.D_no / self.BP_D_no    #scale rise in length
        DN.D_decay = (self.BP_d_decay / 8.64e4) * DN.D_no / self.BP_D_no    #scale decay in length
        outblen_ss = DN.D_rise + DN.D_decay + DN.D_ss 
        #print "outblen_ss", outblen_ss

        outb_interval = {0 : self.Get_outbinterval( DN, DN.minrecurt, 0), 
                         1 : self.Get_outbinterval( DN, DN.minrecurt, 1),
                         2 : self.Get_outbinterval( DN, DN.minrecurt, 2, 
                                                    D_ss = outblen_ss  )}
        
        D_ss     = [DN.D_ss]
        DN.Change_standstill()
        next_ss_outb = t_OutBurst + DN.P_ss
        
        t_OutBurst += outb_interval[ init_outbtype ]
        timeOuts = [t_OutBurst]
        typeOuts = [init_outbtype]
        t_since_so = 0
        
        if DN.P_sc > 2 * DN.P_c:
            nr_no_outb = 0
            if init_outbtype == 0:
                nr_no_outb = 1
            # step forward in time
            while t_OutBurst < 1.1 * max(tvals) / 8.64e4:
                outbtype = self.Get_outbtype(nr_no_outb, DN.nr_no_per_so)
                if t_since_so > DN.P_sc: outbtype = 1
                if t_OutBurst + outb_interval[ outbtype ] > next_ss_outb:
                    outbtype = 2
                if outbtype in [1,2]:
                    t_since_so = 0
                if outbtype == 0:  
                    nr_no_outb += 1
                    t_since_so += outb_interval[ outbtype ]
                elif outbtype == 1: nr_no_outb = 0
                t_OutBurst += outb_interval[ outbtype ]
                D_ss.append(DN.D_ss)
                if outbtype == 2:
                    DN.Change_standstill()
                    outblen_ss = DN.D_rise + DN.D_decay + DN.D_ss
                    outb_interval[2] = self.Get_outbinterval( DN, DN.minrecurt, 
                                                              2, D_ss = outblen_ss)
                    next_ss_outb = t_OutBurst + DN.P_ss
                #print outbtype, outb_interval[ outbtype ], "outbinterval"
                typeOuts.append(outbtype)
                timeOuts.append( t_OutBurst )
        elif DN.P_c < 2 * DN.P_sc:
            nr_so_outb = 0
            if init_outbtype == 1:
                nr_so_outb = 1
             # step forward in time
            while t_OutBurst < 1.1 * max(tvals) / 8.64e4:
                outbtype_so = self.Get_outbtype(nr_so_outb, DN.nr_so_per_no)
                if t_since_so > DN.P_sc: outbtype_so = 0
                if t_OutBurst + outb_interval[ outbtype_so ] > next_ss_outb:
                    outbtype_so = 2
                if outbtype in [0,2]:
                    t_since_so = 0
                if outbtype_so == 1: 
                    outbtype = 0
                elif outbtype_so == 0: 
                    outbtype = 1
                    nr_so_outb = 0
                    t_since_so += outb_interval[ outbtype ]
                elif outbtype_so == 2:
                    outbtype = 2
                    nr_so_outb = 0
                t_OutBurst += outb_interval[ outbtype ]
                D_ss.append(DN.D_ss)
                if outbtype == 2:
                    DN.Change_standstill()
                    outblen_ss = DN.D_rise + DN.D_decay + DN.D_ss
                    outb_interval[2] = self.Get_outbinterval( DN, DN.minrecurt, 
                                                              2, D_ss = outblen_ss)
                    next_ss_outb = t_OutBurst + DN.P_ss
                #print outbtype, outb_interval[ outbtype ],"outbinterval"
                typeOuts.append(outbtype)
                timeOuts.append( t_OutBurst )
        else:   #every time there's either a normal or super outburst
            while t_OutBurst < 1.1 * max(tvals) / 8.64e4:
                if typeOuts[-1] == 0:
                    outbtype = 1
                elif typeOuts[-1] == 1:
                    outbtype = 0
                else:
                    p = np.random.random_sample()
                    if p < 0.5: outbtype = 0
                    else: outbtype = 1
                if t_OutBurst + outb_interval[ outbtype ] > next_ss_outb:
                    outbtype = 2
                D_ss.append(DN.D_ss)
                if outbtype == 0:
                    if DN.D_no > DN.P_c:
                        t_OutBurst += DN.D_no
                    else:
                        t_OutBurst += DN.P_c
                elif outbtype == 1:
                    if DN.D_so > DN.P_sc:
                        t_OutBurst += DN.D_so
                    else:
                        t_OutBurst += DN.P_sc   
                elif outbtype == 2:
                    t_OutBurst += outb_interval[2]
                    DN.Change_standstill()
                    outblen_ss = DN.D_rise + DN.D_decay + DN.D_ss
                    outb_interval[2] = self.Get_outbinterval( DN, DN.minrecurt, 
                                                              2, D_ss = outblen_ss)
                    next_ss_outb = t_OutBurst + DN.P_ss
                typeOuts.append(outbtype)
                timeOuts.append( t_OutBurst )
        # construct the aggregate light curve
        
                      
        # construct the aggregate light curve
        timeOuts = -np.array(timeOuts)

        bandMags = []
        Q_bandMags = []
        for i,band in enumerate(obBands):
            this_bb = DN.Q_mag[band]
            bandMags.append(this_bb)
            Q_bandMags.append(this_bb)
        bandMags = np.array(bandMags)
        Q_bandMags = np.array(Q_bandMags)
    
        for i, tO in enumerate(timeOuts):
            if typeOuts[i] == 0:
                tscale = DN.D_no / self.BP_D_no
                tSamp = np.array(-tvals - tO*8.64e4) / tscale #Scale duration
                Bool1 = tSamp < self.BP_t_e['normal'] #self.BP_D_no *8.64e4
                Bool2 = tSamp > self.BP_t_s['normal']
            if typeOuts[i] == 1:
                tscale =  DN.D_so / self.BP_D_so
                tSamp = np.array(-tvals - tO*8.64e4) / tscale #Scale duration
                Bool1 = tSamp < self.BP_t_e['super'] # self.BP_D_so *8.64e4
                Bool2 = tSamp > self.BP_t_s['super']
            if typeOuts[i] == 2:
                tscale =  DN.D_no / self.BP_D_no
                tSamp = np.array(-tvals - tO*8.64e4) / tscale #Scale duration
                Bool1 = np.ones(len(tSamp), dtype=bool)   #always sample standstills
                Bool2 = np.ones(len(tSamp), dtype=bool)   #always sample standstills
            if any(Bool1 & Bool2):   #i.e. the outburst is *in* the observations.
                if typeOuts[i] == 0:
                    these_bandMags = np.array( self.get_all_normal_bbs( tSamp, 
                                                             obBands, DN ) )
                if typeOuts[i] == 1:
                    these_bandMags = np.array( self.get_all_super_bbs( tSamp, 
                                                             obBands, DN ) )
                if typeOuts[i] == 2:
                    these_bandMags = np.array( self.get_all_standstill_bbs( tSamp, 
                                                             obBands, DN, D_ss[i], tO ) )
                    
                    
                i5 = (these_bandMags != Q_bandMags)
                bandMags[i5] = np.min([these_bandMags[i5],bandMags[i5]], axis=0)
                #bandMags[i5] = these_bandMags[i5]
        """temp:::::
        print "A_no, A_so", DN.A_no, DN.A_so
        print "peakmags", self.peak_mags['g']
        print self.BP_outb_amp
        with open('temp', 'w') as file:
            for A in bandMags:
                file.write(str(A))
                file.write('\n')
        file.close()
        #print 365 - timeOuts
        #print typeOuts
        import matplotlib.pyplot as plt
        print "Distance, Qmag", lc.LC_Dist, DN.quiescent_mag
        print "P_c, P_sc, P_ss", DN.P_c, DN.P_sc, DN.P_ss
        print "D_no, D_so, D_ss", DN.D_no, DN.D_so, DN.D_ss
        print DN.Q_mag
        BoolU = obBands == 'U'
        BoolB = obBands == 'B'
        BoolV = obBands == 'V'
        BoolR = obBands == 'R'
        BoolI = obBands == 'I'
        tvals = np.array(tvals)
        dist_mod = 5.0 * np.log10( lc.LC_Dist*1.0e3 ) - 5.0
        if np.any(bandMags[BoolI] + dist_mod < -20.):
            #print bandMags[BoolI] + dist_mod
            #print 365 - tvals[BoolB] / 8.64e4, bandMags[BoolB]
            plt.plot(365 - tvals[BoolU] / 8.64e4, bandMags[BoolU] + dist_mod, color = '#610061', label='U')
            plt.plot(365 - tvals[BoolB] / 8.64e4, bandMags[BoolB] + dist_mod, color = '#0028ff', label='B')
            plt.plot(365 - tvals[BoolV] / 8.64e4, bandMags[BoolV] + dist_mod, color = '#a6ff00', label='V')
            plt.plot(365 - tvals[BoolR] / 8.64e4, bandMags[BoolR] + dist_mod, color = '#ff0000', label='R')
            plt.plot(365 - tvals[BoolI] / 8.64e4, bandMags[BoolI] + dist_mod, color = '#c20000', label='I')
            plt.gca().invert_yaxis()
            plt.show()
                    
        """

        
        # apply distance modulus
        dval = lc.LC_Dist
        dist_mod = 5.0 * np.log10( dval*1.0e3 ) - 5.0
        visible = False
        for i,label, in enumerate(obBands):
            bandMags[i] += dist_mod + lc.Extinction[label]
            if bandMags[i] < threshold[label]: visible = True
        
        if visible :
            return [lc.LC_RA, lc.LC_DEC], bandMags, Q_bandMags
        else:
            self.N_trans -= 1
            return [],[], []




class Mdwarf_template( transientTemplate ):
    def __init__( self, tag, trData, obRun ):
        transientTemplate.__init__( self, tag, trData, obRun )
        vals               = np.genfromtxt( obRun.MdwarfFile, names=True ).T[tag]
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
        self.Lq            = {'U':qs[0],'B':qs[1],'V':qs[2],'R':qs[3],'I':qs[4],'J':qs[5],'H':qs[6],'K':qs[7],'u':qs[8],'g':qs[9],'r':qs[10],'i':qs[11],'z':qs[12],'q':qs[13],'y':qs[14]}
        self.Epk           = {}
        self.mag_resolution= obRun.mag_resolution
        self.broadbands_rise  = {}
        self.broadbands_decay = {}
        self.LC_time = {'decay': self.deltaT_LC * 24.*3600.}
    def get_blueprints( self, c, dtObs ):
        Nr_of_years = 1.	#Recurrent transient => Nr. of transients does not depend on 
                                #observation time
        N_transients = self.get_Ntransients( c, Nr_of_years )
        self.N_transnonsee += N_transients
        for i in range( N_transients ):
            coords = c.sampleCoords()
            if coords != []:
                t_dummy = -365.25 * Nr_of_years * np.random.random()
                newLC = LC_blueprint( coords, self.peak_mags, 0.0, t_dummy )
                if self.canSee_LC( newLC ):
                    _Coords = (newLC.LC_RA, newLC.LC_DEC, newLC.LC_Dist)
                    newLC.Extinction = c.Grid.Gal_dust.Sample_extinction(*_Coords)
                    if self.canSee_LC( newLC):
                        self.transients.append( newLC )
                        self.N_trans += 1
    def setUpLC( self, colorchoice ):
        print 'for %s' % self.tag
        self.flux_0 = fc.flux_0[colorchoice]
        self.flux_0.update(fc.flux_0['UBVRI'])
        self.Pk_rise = {band:0 for band in np.append('U',self.bands)}
        self.Pk_decay = {band:0 for band in np.append('U',self.bands)}
        for TYPE in ['rise', 'decay']:
            lcdata_lo = h5py.File( self.LCFile + '_%s_%s.hdf5' % (colorchoice, TYPE), 'r' )
            lcdata_up= h5py.File( self.LCFile + '_UBVRI_%s.hdf5' % (TYPE), 'r' )
            tms = lcdata_up['times'][:]
            if TYPE == 'rise':self.broadbands_rise['t'] = tms;
            else: self.broadbands_decay['t'] = tms
            for band in np.append('U',self.bands):
                print "%s band" % band
                if band.islower():
                    lcdata = lcdata_lo
                else:
                    lcdata = lcdata_up                
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
                    self.Pk_rise[band] = np.max(bb)
                    self.LC_time['rise'] = tms[-1]
                elif TYPE == 'decay':
                    self.broadbands_decay[band] = fun
                    self.Pk_decay[band] = np.max(bb)
                    #self.LC_time['decay'] = tms[-1]
                self.Epk[band] = max(self.Pk_rise[band], self.Pk_decay[band])
            if TYPE == 'rise':
                E_rise = np.trapz(lcdata_up['U'][:], lcdata_up['times'][:])
                self.U_rise = lcdata_up['U'][:]
            if TYPE == 'decay':
                E_decay = np.trapz(lcdata_up['U'][:], lcdata_up['times'][:])
                self.U_decay = lcdata_up['U'][:]
        self.E_LC = (E_rise + E_decay)  * 640.4
        lcdata.close()
        
        #self.Max_obs_D     =  Maximum_observing_distance(self.mag_lim, self.peak_mags, self.std_mag_R)
        self.Max_obs_D      = self.Max_Mdwarf_distance()
    def Max_Mdwarf_distance( self ):
        """
        Calculates the maximum distance up to which an Mdwarf can be seen.
        This is done by looking at the maximum energy an Mdwarf can have in
         each filter band and taking the maximum over them.
        """
        max_energy = max(self.dlogE_act[-1], self.dlogE_ina[-1])
        scale      = 10**(max_energy) / self.E_LC
        max_lum = {band:(self.Epk[band] * scale + self.Lq[band]) for band in self.bands}
        maxdist = 0
        self.peak_mags = {}    #Empty dict to be filled in later
        #area in cm^2 of sphere of 10 pc (to calculate absolute magnitudes)
        area_abs_pc = 4 * np.pi * (10 * pc_to_cm) ** 2.    
        for band in self.bands:
            #Calculate the area of the sphere with the distance up to which we
            # can see the M dwarf in case of the largest outburst.            
            area = ( max_lum[band] / self.flux_0[band] * 
                     10 ** (0.4 * self.mag_lim[band]) )    #in cm^2
            dist = np.sqrt( area / (4 * np.pi) ) / ( pc_to_cm )
            if dist > maxdist: maxdist = dist     #in pc
            self.peak_mags[band] = -2.5 * np.log10( max_lum[band]/area_abs_pc/
                                                    self.flux_0[ band ] )
        return maxdist
        
    def get_all_bbs_rise( self, tsample, obbands ):
        bnds = []
        for i,bandLabel in enumerate(obbands):
            this_bb = self.broadbands_rise[bandLabel]( tsample[i] )
            bnds.append(this_bb)
        return np.array(bnds)
    def get_all_bbs_decay( self, tsample, obbands ):
        bnds = []
        for i,bandLabel in enumerate(obbands):
            this_bb = self.broadbands_decay[bandLabel]( tsample[i] )
            bnds.append(this_bb)
        return np.array(bnds)
    def sample_all_LCs( self, obTimes, obBands, threshold ):
        # total time for observation in hours
        """temp:::::
        print obTimes, "obtimes"
        #obBands = [['U', 'B', 'V', 'R', 'I'] for i in range(10000/5)]
        #obTimes = [[i,i,i,i,i] for i in np.linspace(0,100. * 24. * 60. * 60., 10000 / 5)]
        obBands = [['U', 'B', 'V', 'R', 'I'] for i in range(10000/5)]
        obTimes = [[i,i,i,i,i] for i in np.linspace(0,1. * 24. * 60. * 60., 10000 / 5)]
        obBands = np.array(obBands).ravel()
        obTimes = np.array(obTimes).ravel()
        #obBands = np.append(obBands, ['g', 'r'])
        #obTimes = np.linspace(0,100. * 24. * 60. * 60., 1000 / 5)
        for mdwarf in self.transients:
            mdwarf.visibleframes = np.ones(len(obTimes), dtype=bool)
        """
        """
        """
        #print obTimes
        maxtime = 1.e-99
        mintime = 1.e99
        for frame in obTimes:
            maxtime = max(maxtime, max(frame))
            mintime = min(mintime, min(frame))
        tWindow = ( maxtime - mintime )/3600.0 
        for mdwarf in self.transients:
            visible = False
            active = False
            
            vis_frames = np.array(mdwarf.visibleframes)
            """
            obT = []
            obB = []
            for i,frame in enumerate(obTimes):
                if i in vis_frames:
                    obBB = obBands[i]
                    for j,f in enumerate(frame):
                        obT.append(f)
                        obB.append(obBB[j])
            obT = np.array(obT)
            obB = np.array(obB)
            """
            obT = np.concatenate(obTimes[vis_frames])
            obB = np.concatenate(obBands[vis_frames])
            #print "obT", obT, obB
            #print type(obT), np.shape(obT)
            
            # galactic height in parsecs
            hZ = cT.get_galactic_height( mdwarf.LC_RA, mdwarf.LC_DEC, 
                                         mdwarf.LC_Dist ) * 1.0e3 # parsecs
            P_A = self.aval * np.exp( -self.bval * abs(hZ) )  
            ztest = np.random.random()
            if ztest <= P_A: active = True


            #active, hZ = True, 1.

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
                eBins_m = self.dlogE_ina + 0.5 * (eBins_l[1] - eBins_l[0])
            #Determine Flare Frequency
            cum_ne = tWindow * np.power(10.0, alf) * np.power(10.0, bet*eBins_l)
            #print alf, bet, "cum_ne", cum_ne, "eBins_l", eBins_l
            #cum_ne is cumulative. We need to decumulate it:
            ne = [cum_ne[i] - cum_ne[i+1] for i in range(self.Ebins +1)]
            #ne = np.array(ne)
            #z = np.random.random(self.Ebins + 1)  #Round probabilities off to integers
            #cond = z < ne - ne.astype(int)
            #ne[cond] = 1
            #nFlares = ne.astype(int)
            #print nFlares
            nFlares = geometricP(ne)
            
            if max( nFlares ) < 1: 
                self.N_trans -=1
            else:
                lumlist, Qflux = self.get_all_lums( nFlares, eBins_m[:-1], obT, obB )
                maglist = np.zeros( len(obB) )
                Qmaglist= np.zeros( len(obB) )
                for j, lum in enumerate(lumlist):
                    bnd = obB[j]
                    area = 4.0 * np.pi * np.power( mdwarf.LC_Dist * 1e3 * pc_to_cm, 2.0 ) # cm^2
                    mag  = -2.5 * np.log10( lum/area/self.flux_0[ bnd ] )
                    Qmag = -2.5 * np.log10( Qflux[j]/area/self.flux_0[ bnd ] )
                    maglist[j]  = mag
                    Qmaglist[j] = Qmag
                    #the luminosity (in mags) of the quiescent star
                    if bnd in self.Lq.keys():
                        mag_non_outb = -2.5 * np.log10( self.Lq[bnd]/area/self.flux_0[ bnd ] )
                    else:
                        #The quiescent luminosity for band ", bnd, " is not included in the file with M dwarf data
                        mag_non_outb = 0
                    maglist[j]  += mdwarf.Extinction[ bnd ]
                    Qmaglist[j] += mdwarf.Extinction[ bnd ]
                    if mag < threshold[bnd] and mag < mag_non_outb - self.mag_resolution:
                        visible = True
                if visible:
                    """
                    print mdwarf.Extinction, mdwarf.LC_Dist, mdwarf.LC_RA, mdwarf.LC_DEC
                    import matplotlib.pyplot as plt	#plot a single Mdwarf lightcurve
                    BoolU = obBands == 'U'
                    BoolB = obBands == 'B'
                    BoolV = obBands == 'V'
                    BoolR = obBands == 'R'
                    BoolI = obBands == 'I'
                    tvals = np.array(obTimes)
                    plt.plot(tvals[BoolI] / 8.64e4 * 24., maglist[BoolI], color = '#c20000', linewidth=2, label='I')
                    plt.plot(tvals[BoolR] / 8.64e4 * 24., maglist[BoolR], color = '#ff0000', linewidth=2, label='R')
                    plt.plot(tvals[BoolV] / 8.64e4 * 24., maglist[BoolV], color = '#a6ff00', linewidth=2, label='V')
                    plt.plot(tvals[BoolB] / 8.64e4 * 24., maglist[BoolB], color = '#0028ff', linewidth=2, label='B')
                    plt.plot(tvals[BoolU] / 8.64e4 * 24., maglist[BoolU], color = '#610061', linewidth=2, label='U')
                    plt.xlabel('Time (days)')
                    plt.ylabel('Apparent magnitude')
                    plt.legend()
                    plt.gca().invert_yaxis()
                    plt.show()
                    """


                    self.radec_list.append([ mdwarf.LC_RA, mdwarf.LC_DEC ])
                    self.bandmags_list.append( maglist )   
                    self.Qmag_list.append(Qmaglist)
                    self.bandbands_list.append( obB ) 
                    self.obT_list.append(obT)
                else:
                    self.N_trans -=1
        self.bandmags_list = np.array(self.bandmags_list)
    def get_all_lums(self, nFlares, eBins, obT, obbands ):
        window  = obT[-1] - obT[0]
        lumlist = np.zeros( len(obT) )
        Qflux   = emptyval * np.ones( len(obT) ) #Quiescent flux
        eF, nF  = eBins[nFlares>0], nFlares[nFlares>0]
        for i, n in enumerate(nF):
            logEU = eF[i]
            trise  = np.power(10.0, np.random.normal(0.31*logEU - 7.6, 0.32))
            tdecay = np.power(10.0, np.random.normal(0.51*logEU - 13., 0.22))

            thisE = np.power(10.0, logEU)
            dts = np.random.random( n ) * window
            #generate outbursts that started after t0
            for j,dt in enumerate(dts):
                tSamp = obT-dt
                tr_scale = trise / self.LC_time['rise'] 
                td_scale = tdecay / self.LC_time['decay']
                Energy_scale = self.get_Energy_scale(thisE, tr_scale, td_scale)
                tSamp_rise  = tSamp / tr_scale
                tSamp_decay = (tSamp - trise) / td_scale
                new_flux_rise  = Energy_scale *\
                                  self.get_all_bbs_rise( tSamp_rise, obbands )
                new_flux_decay  = Energy_scale *\
                                  self.get_all_bbs_decay( tSamp_decay, obbands)
                lumlist += new_flux_rise + new_flux_decay
            dts = np.random.random( n ) * window
            dts = dts[ dts < (trise + tdecay) * 4. * 24.*60.*60.]
            for j,dt in enumerate(dts):
                tSamp = obT+dt
                tr_scale = trise / self.LC_time['rise'] 
                td_scale = tdecay / self.LC_time['decay']
                Energy_scale = self.get_Energy_scale(thisE, tr_scale, td_scale)
                tSamp_rise  = tSamp / tr_scale
                tSamp_decay = (tSamp - trise) / td_scale
                new_flux_rise  = Energy_scale *\
                                  self.get_all_bbs_rise( tSamp_rise, obbands )
                new_flux_decay  = Energy_scale *\
                                  self.get_all_bbs_decay( tSamp_decay, obbands)
                lumlist += new_flux_rise + new_flux_decay
        for i,bnd in enumerate(obbands):   #Add quiescence flux
            if bnd in self.Lq.keys():
                lumlist[i] += self.Lq[ bnd ]
                Qflux[i] = self.Lq[bnd]
        return lumlist, Qflux
    def get_Energy_scale(self, flare_E, tr_scale, td_scale):
        """
        Calculates the scale with which the flux of the template flare has to 
         be adjusted for this flare
        """
        Time_rise   = self.broadbands_rise['t']
        Time_decay  = self.broadbands_decay['t']
        
        #Find the energy change from the time adjustments alone
        E_new_rise  = np.trapz(self.U_rise, Time_rise * tr_scale)
        E_new_decay = np.trapz(self.U_decay, Time_decay * td_scale)
        E_new_inter = (E_new_rise + E_new_decay) * 640.4

        return flare_E / E_new_inter
                

class xgal_template( transientTemplate  ):
    def __init__( self, tag, trData, obRun ):
        transientTemplate.__init__( self, tag, trData, obRun )
        #self.stellarNrm = 2.47e-14	#for snIa
        self.NM         = float(trData[6])
        self.alpha      = float(trData[13])
        self.HGE_param  = trData[14]    #Host galaxy extinction parameter
        self.Kcormodel  = trData[5]
        self.colorScheme= obRun.colorScheme
        self.Kcorfuns   = cm.Kcor(self.Kcormodel, self.colorScheme, self.bands)
        self.Host_ext   = dust.Host_extinction(self.HGE_param , obRun)
        #The following are only applicable to SNIa
        self.scale      = 1		
        self.scale_0    = 1.09
        self.M_B_max    = -18.89
        
    def setUpLC( self, colorchoice ):
        print 'for %s' % self.tag
        lcdata_up = h5py.File( self.LCFile + '_UBVRI.hdf5', 'r' )
        lcdata_lo = h5py.File( self.LCFile + '_%s.hdf5' % colorchoice, 'r' )
        tms = lcdata_up['times'][:]    #times(UBVRI) = times(ugriz)
        self.broadbands['t'] = tms;
        self.peak_mags = {}   #Empty dict to be filled in later
        Kcorfile_up = h5py.File('LightCurveFiles/Kcorrections/%s_UBVRI.hdf5' % (self.Kcormodel),'r')
        Kcorfile_lo = h5py.File('LightCurveFiles/Kcorrections/%s_%s.hdf5' % (self.Kcormodel, self.colorScheme),'r')
        for band in self.bands:
            print "%s band" % band
            if band.islower():
                lcdata   = lcdata_lo
                Kcorfile = Kcorfile_lo
            else:
                lcdata   = lcdata_up
                Kcorfile = Kcorfile_up
            if not band in lcdata.keys():
                print "no data for %s band" % band
                bb = emptyval * np.ones( len(tms) )
            else:
                scale = self.peak_mag_R - min(lcdata_up['R'][:])
                bb = lcdata[band][:] + scale
            ff = interp1d( tms, bb, bounds_error = False, fill_value = emptyval )
            def fun( times, fcn = ff ):
                return fcn(times)
            self.broadbands[band] = fun
            #Now we need to determine the brightness maximum in each color band
            #This is used in get_blueprints to determine whether the transient is
            # visible or not
            #We need the maximum over ALL redshifts!
            minimum_band = emptyval
            for z in range(len(Kcorfile[band][0,:])):
                new_minimum = min(lcdata[band][:] + Kcorfile[band][:,z])
                minimum_band = min(minimum_band, new_minimum)
            self.peak_mags[band] = minimum_band + scale
        lcdata_lo.close()
        lcdata_up.close()
    def Phillips_scale( self, magDev ):
        """
        The Phillips peak-luminosity-broadness relation.
        Is only applicable to supernovae 1a
        """
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
        return np.array(bnds)
    def get_Ntransients( self, cell, Nr_of_years ):
        kpc3 = cell.vol
        N_per_y = kpc3 * cell.Grid.Cosmo.Nr_density_xtr( cell.z ) 
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
            coords = c.sampleCoords()
            if coords != []:
                tOutburst = -365.25 * Nr_of_years * np.random.random() 
                if tOutburst + self.deltaT_LC + dtObs > 0.0:
                    magDelta = np.random.normal( scale = self.std_mag_R )
                    peak_M = {band: self.peak_mags[band] + magDelta for band in self.bands}
                    newLC = LC_blueprint( coords, peak_M, magDelta, tOutburst )
                    if self.canSee_LC( newLC ):
                        _Coords = (newLC.LC_RA, newLC.LC_DEC, newLC.LC_Dist)
                        Host_ext = self.Host_ext.Sample_host_extinction()
                        MW_ext   = c.Grid.Xtr_dust.Sample_extinction(*_Coords)
                        newLC.Extinction = {band: MW_ext[band] + Host_ext[band] 
                                            for band in MW_ext.keys()}
                        if self.canSee_LC( newLC ):
                            newLC.Redshift = c.Grid.Cosmo.get_redshift(newLC.LC_Dist / 1.e3)	#convert to Mpc
                            self.transients.append(newLC)
                            self.N_trans += 1
    def sample_all_LCs( self, obTimes, obBands, threshold ):
        for lc in self.transients: 
            visible = False
            #print type(obTimes), obTimes
            vis_frames = np.array(lc.visibleframes)
            obT = np.concatenate(obTimes[vis_frames])
            #print "obT", type(obT), obT
            obB = np.concatenate(obBands[vis_frames])
            obT_f = max(obT)
            dts = obT_f - obT
            
            tSamp = -lc.tExplode * 8.64e4 - dts
            dval = lc.LC_Dist	#is in kpc
            dist_mod = 5.0 * np.log10(dval*1.0e3) - 5.0
            magDev = lc.magDelta
            if self.tag == 'SNIa':
                self.Phillips_scale(magDev)
            bandMags = self.get_all_bbs( tSamp, obB, lc.Redshift )
            for i,label in enumerate(obB):
                bandMags[i] += dist_mod + lc.Extinction[label] + magDev
                if bandMags[i] < threshold[label]: visible = True
            #for i,label, in enumerate(bandlabels):
            #    bandMags[i][bandMags[i] < emptyval] += dist_mod + lc.Extinction[label] + magDev
            #bbExt = np.array([np.min( band ) for band in bandMags])
            #if np.any( bbExt < threshold ) :
            #print "bandmags", self.tag, type(bandMags)
            if visible: 
                self.radec_list.append( [ lc.LC_RA, lc.LC_DEC ] )
                self.bandmags_list.append( bandMags )
                self.Qmag_list.append(emptyval * np.ones(len(bandMags)))
                self.bandbands_list.append( obB )
                self.obT_list.append(obT)
            else:
                self.N_trans -=1
            #print type(bandMags)
        self.bandmags_list = np.array(self.bandmags_list)
        
class xgal_no_Kcor_template( transientTemplate  ):
    """
    This template is meant for transients that either do not have any
     K-corrections or are so close by that (and fainter than supernovae) that
     they are only found at very low redshifts.
    The volumetric density therefore also does not need to scale with redshift.
    """
    def __init__( self, tag, trData, obRun ):
        transientTemplate.__init__( self, tag, trData, obRun )
        self.stellarNrm = float(trData[6])
        self.HGE_param  = trData[14]    #Host galaxy extinction parameter
        self.Host_ext   = dust.Host_extinction(self.HGE_param , obRun)
        self.colorScheme= obRun.colorScheme
        #The following are only applicable to SNIa
        self.scale      = 1		
        self.scale_0    = 1.09
        self.M_B_max    = -18.89
        
    def Phillips_scale( self, magDev ):
        """
        The Phillips peak-luminosity-broadness relation.
        Is only applicable to supernovae 1a
        """
        Mag_zero = self.M_B_max + magDev - 5 * np.log10(cm.H0) + 25
        scale = 1./(0.59 * Mag_zero + 3.0)	#The real scale factor
        self.scale = self.scale_0 / scale	#Adjusted for scale factor of template LC
    def get_Ntransients( self, cell, Nr_of_years ):
        kpc3 = cell.vol
        N_per_y = kpc3 * self.stellarNrm	 * Nr_of_years
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
            coords = c.sampleCoords()
            if coords != []:
                tOutburst = -365.25 * Nr_of_years * np.random.random() 
                if tOutburst + self.deltaT_LC + dtObs > 0.0:
                    magDelta = np.random.normal( scale = self.std_mag_R )
                    peak_M = {band: self.peak_mags[band] + magDelta for band in self.bands}
                    newLC = LC_blueprint( coords, peak_M, magDelta, tOutburst )
                    if self.canSee_LC( newLC ):
                        _Coords = (newLC.LC_RA, newLC.LC_DEC, newLC.LC_Dist)
                        Host_ext = self.Host_ext.Sample_host_extinction()
                        MW_ext   = c.Grid.Xtr_dust.Sample_extinction(*_Coords)
                        newLC.Extinction = {band: MW_ext[band] + Host_ext[band] 
                                            for band in MW_ext.keys()}
                        if self.canSee_LC( newLC ):
                            self.transients.append(newLC)
                            self.N_trans += 1
    def sample_all_LCs( self, obTimes, obBands, threshold ):
        for lc in self.transients: 
            visible = False
            
            vis_frames = np.array(lc.visibleframes)
            obT = np.concatenate(obTimes[vis_frames])
            obB = np.concatenate(obBands[vis_frames])
            obT_f = max(obT)
            dts = obT_f - obT
        
            tSamp = -lc.tExplode * 8.64e4 - dts
            dval = lc.LC_Dist	#is in kpc
            dist_mod = 5.0 * np.log10(dval*1.0e3) - 5.0
            magDev = lc.magDelta
            if self.tag == 'SNIa':
                self.Phillips_scale(magDev)
            bandMags = self.get_all_bbs( tSamp, obB )
            for i,label in enumerate(obB):
                bandMags[i] += dist_mod + lc.Extinction[label] + magDev
                if bandMags[i] < threshold[label]: visible = True
            if visible: 
                self.radec_list.append( [ lc.LC_RA, lc.LC_DEC ] )
                self.bandmags_list.append( bandMags )
                self.Qmag_list.append(emptyval * np.ones(len(bandMags)))
                self.bandbands_list.append( obB )
                self.obT_list.append(obT)
            else:
                self.N_trans -=1
        self.bandmags_list = np.array(self.bandmags_list)

class kilonovaTemplate( transientTemplate ):
    def __init__(self, tMerge, dist, mvRed, mvBlue, obRun, Framenr):
        grid            = obRun.SkyGrids[Framenr]
        self.tag        = 'kilonova'
        self.bands      = obRun.bands
        self.t_merge    = tMerge
        if dist != None:
            self.dist = dist
        else: self.dist = sampleLIGODistance( obRun.maxLIGODist)
        self.DEC        = sampleDEC( grid.DEC_lo, grid.DEC_hi )
        RA_lo           = grid.RA_c - (0.5 * obRun.aperture_RA /
                            np.cos(self.DEC * (np.pi/180.)) )
        RA_hi           = grid.RA_c - (0.5 * obRun.aperture_RA /
                            np.cos(self.DEC * (np.pi/180.)) )
        self.RA         = sampleRA( RA_lo, RA_hi )
        self.N_trans    = 1
        self.N_transnonsee = 1
        self.mvRed      = mvRed
        self.mvBlue     = mvBlue
        self.f_band     = {}
        self.netLC      = {}
        self.radec_list    = []
        self.bandmags_list = []
        self.Qmag_list     = []
        self.bandbands_list= []
        self.obT_list      = []
        self.colorchoice= obRun.colorScheme
        self.transients = [LC_blueprint([self.RA, self.DEC, self.dist], 
                                       0., 1., self.t_merge)] #magPk and MagDelta don't matter here.

    def sample_all_LCs( self, obTimes, obBands, threshold ):
        vis_frames = np.array(self.transients[0].visibleframes)

        obT = np.concatenate(obTimes[vis_frames])
        obB = np.concatenate(obBands[vis_frames])
        obT_f = max(obT)
        dts = obT_f - obT
        bandmags_list = np.ones(len(dts)) * emptyval
        
        if self.mvRed != None:
            mred, vred = get_best_match( 'red', *self.mvRed)
            rlcf_up = 'kilonovaLib/knRed_m%.4f_vk%.2f_Xlan1e-1.0_UBVRI.h5' % (mred, vred)
            rlcf_lo = 'kilonovaLib/knRed_m%.4f_vk%.2f_Xlan1e-1.0_%s.h5' % (mred, vred, self.colorchoice)
            redLC_up = h5py.File(rlcf_up, 'r' )
            redLC_lo = h5py.File(rlcf_lo, 'r' )
            if self.mvBlue == None: 
                self.netLC['times'] = redLC_up['times'][:]
                for band in self.bands:
                    if band.islower():
                        self.netLC[band] =  redLC_lo[band][:]
                    else:
                        self.netLC[band] =  redLC_up[band][:]
                redLC_up.close()
                redLC_lo.close()
        if self.mvBlue != None:
            mblu, vblu = get_best_match( 'blue', *self.mvBlue)
            blcf_up = 'kilonovaLib/knBlu_m%.4f_vk%.2f_Xlan1e-5.0_UBVRI.h5' % (mblu, vblu)
            blcf_lo = 'kilonovaLib/knBlu_m%.4f_vk%.2f_Xlan1e-5.0_%s.h5' % (mblu, vblu, self.colorchoice)
            bluLC_up = h5py.File(blcf_up, 'r')
            bluLC_lo = h5py.File(blcf_lo, 'r')
            if self.mvRed == None:
                self.netLC['times'] = bluLC_up['times'][:]
                for band in self.bands:
                    if band.islower():
                        self.netLC[band] =  bluLC_lo[band][:]
                    else:
                        self.netLC[band] =  bluLC_up[band][:]
                bluLC_up.close()
                bluLC_lo.close()
        if not self.netLC: # combine kilonova light curves
            if len(redLC_up['times'][:]) > len(bluLC_up['times'][:]):
                self.netLC['times'] = redLC_up['times'][:]
            else:
                self.netLC['times'] = bluLC_up['times'][:]
            for band in self.bands:
                if band.islower():
                    m1, m2 = redLC_lo[band][:], bluLC_lo[band][:]
                else:
                    m1, m2 = redLC_up[band][:], bluLC_up[band][:]
                if len(m1) > len(m2):
                    m2 = np.append(m2, np.zeros(len(m1) - len(m2)))
                if len(m2) > len(m1):
                    m1 = np.append(m1, np.zeros(len(m2) - len(m1)))
                m12 = -2.5*np.log10( np.power(10.0, -.4*m1) + np.power(10.0, -.4*m2) )
                self.netLC[band] = m12
            redLC_lo.close()
            redLC_up.close()
            bluLC_lo.close()
            bluLC_up.close()
        # now that we have the light curve, sample it at the appropriate times
        tSample = (obT-self.t_merge)
        dist_mod = 5.0 * np.log10( self.dist/pc_to_cm  ) - 5.0
        for j, band in enumerate(self.bands):
            bandLC = self.netLC[band]
            self.f_band[band] = interp1d( self.netLC['times'][ np.isinf(bandLC) == False],
                                   bandLC[np.isinf(bandLC) == False], 
                                   bounds_error=False, fill_value=np.inf )
        visible = False
        for j,band in enumerate(obB):
            #Host_ext = self.Host_ext.Sample_host_extinction()
            #MW_ext   = c.Grid.Xtr_dust.Sample_extinction(*coords)
            #Extinction = { band: MW_ext[band] + Host_ext[band] 
            #               for band in MW_ext.keys() }
            bandmags_list[j] = self.f_band[band]( tSample[j] ) + dist_mod #+ Extinction
            if bandmags_list[j] < threshold[band]:visible = True
        #The following three items should be put into an array like with all other transients
        if visible:
            self.radec_list.append([self.RA, self.DEC])
            self.bandmags_list.append(bandmags_list)
            self.bandbands_list.append(obB)
            self.Qmag_list.append(emptyval * np.ones(len(bandmags_list)))
            self.obT_list.append(obT)
        else:	
            self.N_trans = 0
        self.bandmags_list = np.array(self.bandmags_list)


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
        self.visibleframes = []

class DwarfNova:
    """
    A Dwarf nova class containing all of the relevant parameters.
    All of this is based on either the article or the database of Otulakowska-
     Hypka (2016)
    P_c is the recurrence time of normal outbursts
    P_sc is the recurrence time of super/longer outbursts
    P_ss is the recurrence time for standstills
    A_no is the amplitude of normal outbursts
    A_so is the amplutide of super/longer outbursts
    D_no is the duration of normal outbursts
    D_so is the duration of super/longer outbursts
    D_ss is the duration of the standstill minus rise and decay (which scale with D_no)
    no_per_so is the average number of normal outbursts per super/long outbursts
    """
    def __init__ (self, Type, P_c_s, P_c_loc, P_c_scale, q_mag, bands):
        self.type = Type
        
        self.P_c = lognorm.rvs(P_c_s, P_c_loc, P_c_scale) #days
        mu1       = 1.32 + 0.66 * np.log10( self.P_c )
        self.P_sc = 10**( np.random.normal(loc = mu1, scale=0.24) )   #days
        self.P_ss = (365.*5. - 10.) * np.random.random_sample() + 10.  #days
        
        self.A_no = lognorm.rvs(0.32, -0.89, 3.84)
        mu2       = 1.043 * self.A_no + 0.68
        self.A_so = np.random.normal(loc = mu2, scale=0.41)    #mags
        if self.A_so < 0 : self.A_so = 0
        
        if Type == 'SUUMa':
            self.D_no = lognorm.rvs(0.49, 0.36, 3.35)   #days
        elif Type == 'UGem':
            self.D_no = lognorm.rvs(0.74, 1.83, 6.60)   #days
        elif Type == 'ZCam':
            self.D_no = lognorm.rvs(0.49, 0.35, 8.31)   #days
        else:
            print "Warning: Class DwarfNova in transient.py did not receive a correct type. Using U Gem as default."
            self.D_no = lognorm.rvs(0.74, 1.83, 6.60)   #days
        mu3       = 2.09 * self.A_so + 9.7
        self.D_so = np.random.normal(loc = mu3, scale=5.7)
        if self.D_so < 0 : self.D_so = 0
        self.D_ss = 365.  * np.random.random_sample()  #days
        self.D_rise = 1   #to be changed later
        
        self.D_decay = 1  #to be changed later

        if self.P_sc >  self.P_c:
            self.nr_no_per_so = self.P_sc / self.P_c - 1
            self.nr_so_per_no = 1./self.nr_no_per_so
        else:
            self.nr_so_per_no = self.P_c / self.P_sc - 1
            self.nr_no_per_so = 1./self.nr_so_per_no
            
        qmag_mu = 6.33 + 1.64 * np.log10(min(self.P_c, self.P_sc))  
        self.quiescent_mag = np.random.normal(loc = qmag_mu, scale = 0.3) #in V-band
        self.Q_mag = {band: q_mag[band] + self.quiescent_mag 
                 for band in np.append('V',bands)}
        
        self.maxoutblen = max(self.D_no, self.D_so)
        self.maxrecurt  = max(self.P_c, self.P_sc)
        self.minrecurt  = min(self.P_c, self.P_sc)
    def Change_standstill( self ):
        """
        Standstills change in duration and frequency so much more than other
         outburst types, that we need to change them after every standstill.
        This is done in this function.
        """
        self.D_ss = (365. - 10.) * np.random.random_sample() + 10.  #days
        self.P_ss = (365.*5. - 10.) * np.random.random_sample() + 10.  #days
        """
        temp::::::
        self.P_ss = 100.
        self.D_ss = 100.
        """
        """
        temp
        """
        
    def Get_max_lookbacktime( self, D_rise = 0, D_decay = 0 ):
        """
        Calculates the maximum time we have to look back in to sample the first
         outburst time
        """
        if self.type != 'ZCam':
            return self.maxoutblen + self.maxrecurt
        else:
            if self.P_ss > self.maxrecurt:
                #print "P_ss, D_ss, D_rise, D_decay", self.P_ss, self.D_ss, D_rise, D_decay
                return self.P_ss + self.D_ss + D_rise + D_decay
            else:
                return self.maxoutblen + self.maxrecurt
            
        


def Maximum_observing_distance(mag_lim, pkmag, devmag):
    """
    Calculates up to what distance we should generate transients. 
    This is done by calculating the peak magnitude with a 3 sigma standard deviation.
    """
    Pkmag = np.array([ pkmag[ band ] for band in pkmag.keys()]).astype(float)
    mag_lim = np.array([ mag_lim[ band ] for band in pkmag.keys()]).astype(float)
    exp_term = 0.2*( max( mag_lim - ( Pkmag - 3.0*devmag ) ) ) + 1.0
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
    f_0 = fc.flux_0[obRun.colorScheme]
    f_0.update(fc.flux_0['UBVRI'])
    Qmag = { band : flux2mag(Qflux[band], f_0[band]) for band in colorbands}

    #The difference between the magnitude limit and quiescent magnitude of the M dwarf
    difmag = {band : Qmag[band] - obRun.mag_limit[band] for band in colorbands}
    mindifmag, mindifband = min( zip(difmag.values(), difmag.keys()) )
    if mindifmag < 0:	#The M dwarf is visible in quiescence
        logminE = 28
        #####################################
        #Here a function should come to investigate the minimum flare required for an M dwarf that is visible. I.e. work with the sensitivity of the telescope: obRun.mag_resolution
        ######################33
    else:		#The M dwarf is not visible in quiescence
        #Now calculate how large a flare would be in the U-band to increase 
        # the M dwarf's brightness by mindifmag magnitudes
        Qflux_U = Lq['U'] / Area
        Qmag_U = flux2mag( Qflux_U, f_0['U'] )
        Outburstflux_U = mag2flux( Qmag_U - mindifmag, f_0['U'])
        L_outburst_U = Outburstflux_U * Area
        #print Qmag_U, mindifmag, mag2flux(Outburstflux_U, f_0['U']), L_outburst_U
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
    There is a non-zero probability that the number of transients is higher than nt.
    nt: the average number of transients in this cell (float)
    returns: the actual number of transients in this cell (int)
    """
    return np.random.poisson(nt)

def PeakMagData( File, Tag, bands ):
    """
    Extracts the mean (+std) peak magnitudes for a transient from
     the file.
    File: The file with peak magnitudes
    Tag: The name of the transient
    bands: the color bands in use 
    returns: the peak magnitudes in the color bands in use + the standard deviation
    """
   # transientsData = getFileLine( File, Tag )
    PD = getFileLine( File, Tag )
    PeakMagData = {'U':PD[0],  'B':PD[1],  'V':PD[2],  'R':PD[3],  'I':PD[4], 
                   'J':PD[5],  'H':PD[6],  'K':PD[7],  'u':PD[8],  'g':PD[9], 
                   'r':PD[10], 'i':PD[11], 'z':PD[12], 'y':PD[13], 'q':PD[14], 
                   'std':PD[15]}
    PMData = {band: float(PeakMagData[band]) for band in bands}
    return PMData, float(PeakMagData['std'])

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
    #ra_conv = np.array([15.0, 0.25, 1.0/240.])
    #hmslow = np.array([float(val) for val in ra_lo.split(':')])
    #hmshi= np.array([float(val) for val in ra_hi.split(':')])
    #ralo = np.sum( hmslow*ra_conv )
    #rahi = np.sum( hmshi*ra_conv )
    #if ralo > rahi: ralo, rahi = rahi, ralo
    if ra_lo > ra_hi: ra_lo, ra_hi = ra_hi, ra_lo
    #RA_sampled = np.random.random() * (rahi - ralo ) + ralo
    RA_sampled = np.random.random() * (ra_hi - ra_lo) + ra_lo
    return RA_sampled

def sampleDEC( dec_lo, dec_hi ):
    #dec_conv = np.array([1.0, 1.0/60.0, 1.0/3600.0])
    #dmslow = np.array([float(val) for val in dec_lo.split(':')])
    #dmshi= np.array([float(val) for val in dec_hi.split(':')])
    #declo = np.sum( dmslow*dec_conv )
    #dechi = np.sum( dmshi*dec_conv )
    #if declo > dechi: declo, dechi = dechi, declo
    if dec_lo > dec_hi: dec_lo, dec_hi = dec_hi, dec_lo
    #DEC_sampled = np.random.random() * (dechi - declo) + declo
    DEC_sampled = np.random.random() * (dec_hi - dec_lo) + dec_lo
    return DEC_sampled

def sampleLIGODistance(  dmax ):
    while True:
        x = np.random.random()*dmax
        y = np.random.random()*dmax
        z = np.random.random()*dmax
        r = np.sqrt( x*x + y*y + z*z )
        if r <= dmax: break
    return r

