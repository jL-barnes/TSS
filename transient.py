import numpy as np
import h5py
import MilkyWay as MW
import MilkyWay2 as MW2
import cosmology as cm
import Filterchar as fc
import dust
import glob
from scipy.stats import lognorm
from scipy.interpolate import interp1d
from utils import mag2flux, flux2mag, geometricP, getFileLine, get_galactic_height


year_to_seconds = 8.64e4 * 365.25
day_to_seconds = 8.64e4
pc_to_cm = 3.08567758e18
emptyval = 1000.0



class TransientSet:
    """ A Transient set contains all transients for a single 
    observation frame
    ...

    Attributes
    ----------
    obRun : object
        Observation object instance
    Framenr : int
        The number of this frame
    """
    def __init__ (self, obRun, Framenr):
        self.Framenr = Framenr
        self.grid = obRun.SkyGrids[self.Framenr]
        self.transient_templates = []
        for tr in obRun.transientsList:
            new_template = makeTemplate(tr, obRun)
            self.transient_templates.append(new_template)
        self.deltaT_obs  = ( (obRun.TimeParams.t_end 
                            - obRun.TimeParams.t_start) / 8.64e4 )
        self.setUpLC_funs(obRun.colorScheme)

    def setUpLC_funs(self, colorchoice):
        """ For all transient types set up a functional broadband form 
             from which we can sample at the times of observation
        """
        print('setting up functional broadband forms...')
        for trtemp in self.transient_templates:
            trtemp.setUpLC(colorchoice)
        print('done!')

    def populate(self):
        """ Populates an observation grid with transient objects and 
             obtains so called blueprints that are used to generate 
             observations when take_data() is called
        """
        galIdx = self.grid.N_RA * self.grid.N_DEC * self.grid.N_DMW
        if self.grid.Gal_trans:
            self.populate_galactic(galIdx)
        if self.grid.Xgal_trans:
            self.populate_xgal(galIdx)

    def populate_galactic(self, galIdx):
        """ Populates a frame with Galactic transients

           Stops populating once a distance > D_gal_max is reached. 
            Past that distance the transients won't be visible anyway.
        """
        galTrans = [tr for tr in self.transient_templates if tr.galType != 3]
        cells_galactic = self.grid.cellGrid[0 : galIdx]
        for gtr in galTrans:
            D_gal_max = gtr.Max_obs_D * 1.e-3 #kpc
            for c in cells_galactic:
                if D_gal_max < c.DMid + c.hD: 
                    break
                gtr.get_blueprints(c, self.deltaT_obs)

    def populate_xgal(self, galIdx):
        """ Populates a frame with extragalactic transients
        """
        xgalTrans3 = [tr for tr in self.transient_templates 
                      if tr.galType == 3 ]
        xgalTrans4 = [tr for tr in self.transient_templates 
                      if tr.galType == 4 ]
        cells_xgal = self.grid.cellGrid[galIdx:]
        for xtr in xgalTrans3:
            self.grid.Cosmo.create_Nr_density_funct(xtr.NM, xtr.alpha) 
            for c in cells_xgal:
                xtr.get_blueprints(c, self.deltaT_obs)
        for xtr in xgalTrans4:
            for c in cells_xgal:
                xtr.get_blueprints(c, self.deltaT_obs)

    def inject_kilonova(self, tMerge, dist, mvRed, mvBlue, obRun, Framenr):
        """ Injects a single kilonova
        """
        self.transient_templates.append( kilonovaTemplate(
            tMerge, dist, mvRed, mvBlue, obRun, Framenr) )
        if obRun.Trlegend == {}:
            obRun.Trlegend['kilonova'] = 0
        else:
            max_val = max(obRun.Trlegend, key=obRun.Trlegend.get)
            obRun.Trlegend['kilonova'] = obRun.Trlegend[max_val] + 1
            

def makeTemplate( tag, obRun ):
    """ Create the correct template for each transient type

        -------
        returns: Transient template
    """
    transientsData = getFileLine(obRun.transientFile, tag)
    Type = int(transientsData[0])
    if Type == 0 :
        temp = galactic_nRec_template(tag, transientsData, obRun)
    elif Type == 1:
        if tag == 'UGem':
            temp = UGem_template(tag, transientsData, obRun)
        elif tag == 'SUUMa':
            temp = SUUMa_template(tag, transientsData, obRun)
        elif tag == 'ZCam':
            temp = ZCam_template(tag, transientsData, obRun)
        else:
            temp = galactic_recur_template(tag, transientsData, obRun)
    elif Type == 2:
        temp = Mdwarf_template(tag, transientsData, obRun)
    elif Type == 3:
        temp = xgal_template(tag, transientsData, obRun)
    elif Type == 4:
        temp = xgal_no_Kcor_template(tag, transientsData, obRun)
    return temp


class transientTemplate:
    """ A transient template. Each type of transient has its own 
         template. 

        This class contains the basic elements that a template must 
         have. It is basically the template for Galactic non-recurrent 
         transients.
        This template basically reads out dataTransients.dat, the file 
         with all transient parameters.
    """
    def __init__(self, tag, transientsData, obRun):
        self.tag         = tag
        self.galType     = int(transientsData[0])
        self.galkey      = [int(tD) for tD in transientsData[1:5]]
        self.bands       = obRun.bands
        self.mag_lim     = obRun.threshold
        self.broadbands  = {} 
        self.LCFile      = 'LightCurveFiles/' + transientsData[5]
        self.stellarNrm  = float(transientsData[6]) / MW.rho_stellar_sun
        self.scaleH      = float(transientsData[7])
        self.N_trans     = 0    # how many transients each template has
        self.transients  = []   # list of light curves blueprints here
        self.N_transnonsee = 0  # Nr of generated transients, regardless
                                #  of visibility
        self.peak_mag_R    = float(transientsData[10])
        self.std_mag_R     = float(transientsData[11])
        self.radec_list    = []
        self.bandmags_list = []
        self.bandbands_list= []
        self.Qmag_list     = []
        self.obT_list      = []
        if self.tag in ['UGem', 'SUUMa', 'ZCam']:
            self.deltaT_LC = transientsData[8].split(',')
        else:
            self.deltaT_LC = float(transientsData[8])

    def setUpLC(self, colorchoice):
        """ An interpolation function of the light curve data is made 
             for each filter that will be used in the observation.

            Later, this function can be sampled to obtain the transient 
             brightness at the time of observation.
            Also, the peak magnitude of the LC and the maximum 
             observation distance up to which a transient is still 
             visible are determined.
        """
        print('for %s' % self.tag)
        lcdata_up = h5py.File(self.LCFile + '_UBVRI.hdf5', 'r')
        lcdata_lo = h5py.File(self.LCFile + '_%s.hdf5' % colorchoice, 'r')
        if not np.all(lcdata_up['times'][:] == lcdata_lo['times'][:]):
            print('Warning: the times array in the LC data for UBVRI is', \
                   'not equal to %s' % (colorchoice))
        tms = lcdata_up['times'][:]    #times(UBVRI) = times(ugriz)
        self.broadbands['t'] = tms
        self.peak_mags = {}   #Empty dict to be filled in later
        for band in self.bands:
            print("%s band" % band)
            if band.islower():
                lcdata = lcdata_lo
            else:
                lcdata = lcdata_up
            if not band in lcdata.keys():
                print("no data for %s band" % band)
                bb = emptyval * np.ones(len(tms))
            else:
                scale = self.peak_mag_R - min(lcdata_up['R'][:])
                bb = lcdata[band][:] + scale
            ff = interp1d(tms, bb, bounds_error=False, fill_value=emptyval)
            def fun(times, fcn=ff):
                return fcn(times)
            self.broadbands[band] = fun
            #Now we need to determine the brightness maximum in each color band
            #This is used in get_blueprints to determine whether the transient 
            # is visible or not
            self.peak_mags[band] = min(lcdata[band][:])
        lcdata_lo.close()
        lcdata_up.close()
        self.Max_obs_D =  Maximum_observing_distance(self.mag_lim, 
            self.peak_mags, self.std_mag_R)

    def get_all_bbs(self, tsample, obbands):
        """At Time=tsample sample the light curve function for band = obbands

            -------
            returns: The magnitudes of this outburst
        """
        bnds = []
        for i,bandLabel in enumerate(obbands):
            this_bb = self.broadbands[bandLabel](tsample[i])
            bnds.append(this_bb)
        return np.array(bnds)

    def get_Ntransients(self, cell, Nr_of_years):
        """ Generate a certain number of transients in the given cell.

            This depends on the transient density, the location of the 
             cell and the number of years the simulation spans (many 
              densities are given in N per year)

            rho_ddh is the stellar density in this cell that is part of both
             disks and the halo
            rho_b is the stellar density in this cell that is part of the 
             bulge

            -------
            returns: The number of transients in this cell
        """
        norm = self.get_N_normalization()
        if abs(self.scaleH - MW.MW_Hthin) < 0.001:
            rho_tnd, rho_tkd, rho_b, rho_h = cell.rho * self.galkey
            rho_total = norm * (rho_tnd + rho_tkd + rho_h) + rho_b
            N_per_y = self.stellarNrm * cell.vol * rho_total
        else:
            this_rho = MW.get_MW_dens([cell.raMid, cell.decMid, cell.DMid],
                                      Hthin=self.scaleH)
            rho_tnd, rho_tkd, rho_b, rho_h = this_rho * self.galkey
            rho_total = norm * (rho_tnd + rho_tkd + rho_h) + rho_b
            N_per_y = self.stellarNrm * cell.vol * rho_total
        N_per_y = N_per_y * Nr_of_years
        return geometricP(N_per_y) 

    def get_N_normalization(self):
        """ Calculate the normalization factor for a transient density.

        Some Galactic transients require e.g. a single disk MW profile.
        With galkey one can switch on/off the four components of the 
         Milky Way. This function calculates the normalization factor
         that corresponds to the components that are switched on.
        Because the stellar density in the Solar neighborhood of the 
         bulge component is so low (e-101), the bulge normalization
         factor is independent on whether the other components are 
         switched on or off.
        This normalization factor needs not to be applied to the bulge.

        galkey = [thin disk, thick disk, bulge, halo]

        -------
        returns: The normalization factor norm
        """
        gal_factor = np.array([1.0, MW.MW_fThick, 0.0, MW.MW_fHalo])
        norm = 1. / np.sum( gal_factor[np.where(self.galkey)] )
        return norm
        
 
    def get_blueprints(self, c, dtObs):
        """ Create the blueprints for this type of transient

            For every transient that is generated, we need to generate 
             a time of outburst, draw a peak magnitude and test if we 
             can still see this transient given the position, dust 
             extinction and peak mag.
        """
        Nr_of_years = (self.deltaT_LC + dtObs) / 365.25   
        #At first observation the first LC should have just ended
        N_transients = self.get_Ntransients(c, Nr_of_years)
        self.N_transnonsee += N_transients
        for i in range(N_transients):
            coords = c.sampleCoords()
            if coords != []:
                tOutburst = -365.25 * Nr_of_years * np.random.random()
                if tOutburst + self.deltaT_LC + dtObs > 0.0:
                    magDelta = np.random.normal(scale=self.std_mag_R)
                    peak_M = {band: self.peak_mags[band] + magDelta 
                              for band in self.bands}
                    newLC = LC_blueprint(coords, peak_M, magDelta, tOutburst)
                    if self.canSee_LC(newLC):
                        _Coords = (newLC.LC_RA, newLC.LC_DEC, newLC.LC_Dist)
                        newLC.Extinction = c.Grid.Gal_dust.Sample_extinction(
                            *_Coords)
                        if self.canSee_LC(newLC):
                            self.transients.append(newLC)
                            self.N_trans += 1

    def sample_all_LCs(self, obTimes, obBands, threshold):
        """ For every transient, sample the light curve function at 
            obTimes in filters obBands.

           Only use the observations of the frames in which the 
            transient is visible.
           If the brightness at the observation times meets the 
            threshold, keep it, otherwise throw it away.
        """
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
            bandMags = self.get_all_bbs(tSamp, obB)
            for i,label, in enumerate(obB):
                bandMags[i] += dist_mod + lc.Extinction[label] + magDev
                if bandMags[i] < threshold[label]: visible = True
            if visible:
                self.radec_list.append([lc.LC_RA, lc.LC_DEC])
                self.bandmags_list.append(bandMags)
                self.bandbands_list.append(obB)
                self.Qmag_list.append(emptyval * np.ones(len(bandMags)))
                self.obT_list.append(obT)
            else:
                self.N_trans -=1
        self.bandmags_list = np.array(self.bandmags_list)

    def take_data(self, colorchoice, colorsys):
        """ Take the data from the transient template object.

            Converts Vega magnitudes to AB magnitudes or vice versa if 
             needed
            Returns the observational data to observ.py. This data
             contains an entry/row for every observation.

            -------
            returns: A list of RA/DEC coordinates, all magnitudes in 
                     all bands, quiescent magnitudes, the corresponding 
                     bands and observation times
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
        return (self.radec_list, self.bandmags_list, self.Qmag_list, 
                self.bandbands_list, self.obT_list)

    def canSee_LC(self, lightcurve):
        """ A crude way to estimate whether the light curve should still  
             be visible

            If the Extinction is not zero, it has been adjusted to a real 
             value.
            That lets us calculate the visibility of the transient WITH 
             extinction.

            -------
            returns: A bool whether the transient can be seen or not
        """
        seen = False
        d_10pc = 0.01
        dist_ratio = lightcurve.LC_Dist / d_10pc
        apparent_mag = {band: lightcurve.PeakMag[band]
                        + 5.0 * np.log10(dist_ratio) for band in self.bands}
        if lightcurve.Extinction == {}:
            for band in self.bands:
                if apparent_mag[band] < self.mag_lim[band]: 
                    seen = True
                    return seen
        else:
            for band in self.bands:
                if (apparent_mag[band] + lightcurve.Extinction[band] 
                        < self.mag_lim[band]): 
                    seen = True
                    return seen
        return seen


class galactic_nRec_template(transientTemplate):
    """ The template for non-recurrent transients.

        Inherits the standard transient template and adds nothing
    """
    def __init__(self, tag, trData, obRun):
        transientTemplate.__init__(self, tag, trData, obRun)


class galactic_recur_template(transientTemplate):
    """ The template for recurrent transients.

        Inherits the standard transient template and adds a few extra
         parameters. In this template P_c_s/loc/scale are parameters
         for scipy lognormal fits for the recurrence time of outbursts.
    """
    def __init__(self, tag, trData, obRun):
        transientTemplate.__init__(self, tag, trData, obRun)
        freq_normal = trData[9].split(',')
        self.P_c_s    = float(freq_normal[0].replace("(", ""))
        self.P_c_loc  = float(freq_normal[1])
        self.P_c_scale= float(freq_normal[2].replace(")", ""))

    def get_blueprints(self, c, dtObs):
        """ Create the blueprints for this type of transient

            For every transient that is generated, we need to generate 
             a time of outburst, draw a peak magnitude and test if we 
             can still see this transient given the position, dust 
             extinction and peak mag.
        """
        Nr_of_years = 1.	#Recurrent transient => Nr. of transients does 
                         # not depend on observation time
        N_transients = self.get_Ntransients(c, Nr_of_years)
        self.N_transnonsee += N_transients
        for i in range(N_transients):
            coords = c.sampleCoords()
            if coords != []:
                tOutburst = -365.25 * Nr_of_years * np.random.random()
                magDelta = self.std_mag_R * 3.0
                peak_M = {band: self.peak_mags[band] + magDelta 
                          for band in self.bands}
                newLC = LC_blueprint(coords, peak_M, magDelta, tOutburst)
                if self.canSee_LC(newLC):
                    _Coords = (newLC.LC_RA, newLC.LC_DEC, newLC.LC_Dist)
                    newLC.Extinction = c.Grid.Gal_dust.Sample_extinction(
                        *_Coords)
                    if self.canSee_LC(newLC):
                        self.transients.append(newLC)
                        self.N_trans += 1

    def sample_all_LCs(self, obTimes, obBands, threshold):
        """ For every transient, sample the light curve function at 
            obTimes in filters obBands.

           Only use the observations of the frames in which the 
            transient is visible.
           If the brightness at the observation times meets the 
            threshold, keep it, otherwise throw it away.
        """
        for lc in self.transients:
            vis_frames = np.array(lc.visibleframes)
            obT = np.concatenate(obTimes[vis_frames])
            obB = np.concatenate(obBands[vis_frames])
            obT_f = max(obT)
            dts = obT_f - obT
            radec, mags, Qmags = self.sample_single_LC(lc, dts, obB, 
                                                       threshold)
            if len(radec) > 0:
                self.radec_list.append(radec)
                self.bandmags_list.append(mags)
                self.Qmag_list.append(Qmags)
                self.bandbands_list.append(obB)
                self.obT_list.append(obT)
        self.bandmags_list = np.array(self.bandmags_list)

    def sample_single_LC(self, lc, tvals, obBands, threshold ):
        """ Sample the light curve of a single transient

            For these recurrent transients we first sample a normal 
             outburst frequency from a scipy lognormal distribution.
            The initial outburst explosion time is noted and we add
             outbursts after and before this initial outburst.
            We construct an aggregate light curve. If outbursts overlap
             their fluxes are added to each other.
            At last the distance modulus is applied plus dust extinction
             and the transient is tested for visibility.

            -------
            returns: A list of RA,DEC coordinates, magnitudes and 
                     quiescent magnitudes
        """
        #Sample outburst frequency
        self.P_c = lognorm.rvs(self.P_c_s, self.P_c_loc, self.P_c_scale)
        # figure out the times for each outburst
        timeOuts = []
        t_OutBurst = lc.tExplode
        # step forward in time
        while t_OutBurst < self.deltaT_LC + (max(tvals) / 8.64e4):
            timeOuts.append(t_OutBurst)
            t_OutBurst += self.P_c
        # step backward in time
        t_OutBurst = lc.tExplode - self.P_c
        while t_OutBurst + self.deltaT_LC > (-min(tvals) / 8.64e4): 
            timeOuts.insert(0, t_OutBurst)
            t_OutBurst -= self.P_c
        # construct the aggregate light curve
        if len( timeOuts ) > 0:
            timeOuts = timeOuts[::-1]
            tO = timeOuts[0]
            tSamp = -tvals - tO * 8.64e4
            # deviation for this occurrence of the transient
            magDev = np.random.normal(scale=self.std_mag_R)
            bandMags = self.get_all_bbs(tSamp, obBands)
            for band in bandMags:
                band[band < emptyval] +=  magDev
            for tO in timeOuts[1:]:
                if tO > (-tvals[0] / 8.64e4) - self.deltaT_LC:
                    tSamp = -tvals - tO * 8.64e4
                    magDev = np.random.normal(scale=self.std_mag_R)
                    these_bandMags = self.get_all_bbs(tSamp, obBands)
                    #Check for overlap with previous outbursts
                    for bi, band in enumerate(these_bandMags):    
                        band[band < emptyval] += magDev
                        oband = bandMags[bi]
                        i3 = (band != emptyval) & (oband == emptyval)
                        i4 = (band != emptyval) & (oband != emptyval)
                        bandMags[bi][i3] = band[i3]
                        bandMags[bi][i4] = -2.5 * np.log10(
                            np.power(10, -.4*band[i4]) 
                            + np.power(10, -.4*oband[i4]))
        # apply distance modulus
        dval = lc.LC_Dist
        dist_mod = 5.0 * np.log10(dval*1.0e3) - 5.0
        visible = False
        Qmags = emptyval * np.ones(len(bandMags))
        for i,label, in enumerate(obBands):
            bandMags[i] += dist_mod + lc.Extinction[label]
            if bandMags[i] < threshold[label]: visible = True
        if visible:
            return [lc.LC_RA, lc.LC_DEC], bandMags, Qmags
        else:
            self.N_trans -= 1
            return [],[],[]


class SUUMa_template(galactic_recur_template):
    """ The template for SU UMa dwarf nova transients.

        Inherits the standard recurrent transient template.
    """
    def __init__(self, tag, trData, obRun):
        galactic_recur_template.__init__(self, tag, trData, obRun)
        self.type              = 'SUUMa'
        self.LCFile            = 'LightCurveFiles/Dwarfnovae/' + trData[5]
        self.broadbands_normal = {}
        self.broadbands_super  = {}
        self.quiescentmag      = {}
        self.BP_Q_mag          = {}

    def setUpLC(self, colorchoice):
        """ An interpolation function of the light curve data is made 
             for each filter that will be used in the observation.

            This is done for normal as well as super outbursts.
            Later, this function can be sampled to obtain the transient 
             brightness at the time of observation.
            Also, the peak magnitude of the LC and the maximum 
             observation distance up to which a transient is still 
             visible are determined.
        """
        print('for %s' % self.tag)
        self.BP_outb_amp = {}
        self.BP_D_no = float(self.deltaT_LC[0].replace("(", ""))
        self.BP_D_so = float(self.deltaT_LC[1].replace(")", ""))
        self.BP_t_s  = {}
        self.BP_t_e  = {}
        for TYPE in ['normal', 'super']:
            lcdata_up = h5py.File(self.LCFile 
                                  + '_UBVRI_%s.hdf5' % (TYPE), 'r')
            lcdata_lo = h5py.File(self.LCFile 
                                  + '_%s_%s.hdf5' % (colorchoice, TYPE), 'r')
            self.BP_outb_amp[TYPE], self.BP_Q_mag[TYPE] = self.Get_outb_amp(
                lcdata_up, lcdata_lo, np.append('V',self.bands))
            if not np.all(lcdata_up['times'][:] == lcdata_lo['times'][:]):
                print('Warning: the times array in the LC data for UBVRI',\
                      'is not equal to %s' % (colorchoice))
            tms = lcdata_up['times'][:]    #times(UBVRI) = times(ugriz)
            self.BP_t_s[TYPE] = tms[0]
            self.BP_t_e[TYPE] = tms[-1]
            if TYPE == 'normal':
                self.broadbands_normal['t'] = tms
            if TYPE == 'super':
                self.broadbands_super['t'] = tms
            self.peak_mags = {}   #Empty dict to be filled in later
            for band in np.append('V',self.bands):
                print("%s band" % band)
                if band.islower():
                    lcdata = lcdata_lo
                else:
                    lcdata = lcdata_up
                if not band in lcdata.keys():
                    print("no data for %s band" % band)
                    bb = emptyval * np.ones(len(tms))
                else:
                    #Only record relative change
                    bb = lcdata[band][:] - max(lcdata[band][:]) 
                ff = interp1d(tms, bb, bounds_error=False, fill_value=emptyval)
                this_Q_mag = max(lcdata[band][:]) - self.BP_Q_mag[TYPE]
                if band not in self.quiescentmag.keys():
                    self.quiescentmag[band] = this_Q_mag
                else:
                    self.quiescentmag[band] = max(self.quiescentmag[band], 
                                                  this_Q_mag)
                def fun(times, fcn=ff):
                    return fcn(times)
                if TYPE == 'normal':
                    self.broadbands_normal[band] = fun
                elif TYPE == 'super':
                    self.broadbands_super[band] = fun
                #Now we need to determine the brightness maximum in each band
                #This is used in get_blueprints to determine whether the 
                # transient is visible or not
                self.peak_mags[band] = min(lcdata[band][:])
                
            lcdata_lo.close()
            lcdata_up.close()
        self.Max_obs_D = Maximum_observing_distance(self.mag_lim, 
            self.peak_mags, self.std_mag_R)

    def get_all_normal_bbs(self, tsample, obbands, DN):
        """ At Time=tsample sample the normal outburst light curve 
             function for band = obbands.

            In other color filters than 'V' the outburst is scaled such
             that the color at peak magnitude is maintained the same as
             for the 'standard' dwarf nova.
            What is added to the quiescent magnitude (i.e. the outburst)
             may never be negative (=> scale must be positive)
   
            -------         
            returns: a numpy array with normal outburst magnitudes at t=tsample
        """
        bnds = []
        for i,bandLabel in enumerate(obbands):
            outb_amp_old = self.BP_outb_amp['normal']
            scale = ((DN.A_no + outb_amp_old[bandLabel] - outb_amp_old['V'])
                     / (outb_amp_old[bandLabel]))
            if scale < 0.: scale = 0
            this_bb = (self.broadbands_normal[bandLabel](tsample[i]) * scale
                       + DN.Q_mag[bandLabel])
            if this_bb > DN.Q_mag[bandLabel]:
                this_bb = DN.Q_mag[bandLabel]
            bnds.append(this_bb)
        return np.array(bnds)

    def get_all_super_bbs(self, tsample, obbands, DN):
        """ At Time=tsample sample the super outburst light curve 
             function for band = obbands.

            In other color filters than 'V' the outburst is scaled such
             that the color at peak magnitude is maintained the same as
             for the 'standard' dwarf nova.
            What is added to the quiescent magnitude (i.e. the outburst)
             may never be negative (=> scale must be positive)
            
            -------
            returns: a numpy array with super outburst magnitudes at t=tsample
        """
        bnds = []
        for i,bandLabel in enumerate(obbands):
            outb_amp_old = self.BP_outb_amp['super']
            scale = ((DN.A_no + outb_amp_old[bandLabel] - outb_amp_old['V']) 
                     / (outb_amp_old[bandLabel]))
            if scale < 0.: scale = 0
            this_bb = (self.broadbands_super[bandLabel](tsample[i]) * scale
                       + DN.Q_mag[bandLabel])
            if this_bb > DN.Q_mag[bandLabel]:
                this_bb = DN.Q_mag[bandLabel]
            bnds.append(this_bb)
        return np.array(bnds)

    def sample_single_LC(self, lc, tvals, obBands, threshold):
        """ Sample the light curve of a single dwarf nova

            We first create a Dwarf Nova object.
            The initial outburst explosion time is noted and we add
             outbursts after and before this initial outburst.
            This is done in a specific manner: taking the difference
             between the normal outburst and superoutburst interval into
             account, we progressively add outburst times. After a 
             certain number of normal outbursts, a superoutburst occurs.
            The duration of outbursts are scaled by scaling the time
             array that samples the outburst light curves.
            We construct an aggregate light curve. If outbursts overlap
             the brightest 'measurement' is taken.
            At last the distance modulus is applied plus dust extinction
             and the transient is tested for visibility.

            -------
            returns: A list of RA,DEC coordinates, magnitudes and 
                     quiescent magnitudes
        """
        DN = DwarfNova(self.type, self.P_c_s, self.P_c_loc, self.P_c_scale,
                       self.quiescentmag, self.bands)
        max2 = DN.Get_max_lookbacktime() #max time before t=0 to be simulated
        #Generate the time of the first outburst
        t_OutBurst = -max2 * np.random.random_sample() - max2   
        if DN.nr_no_per_so > 1.:
            init_outbtype = np.random.choice([0,1], p=[1 - 1./DN.nr_no_per_so, 
                                                       1./DN.nr_no_per_so])
        else:
            init_outbtype = np.random.choice([0,1], p=[1./DN.nr_so_per_no, 
                                                       1 - 1./DN.nr_so_per_no])
        # 0 = normal outburst
        # 1 = super outburst
        outb_interval = {0 : self.Get_outbinterval(DN, DN.minrecurt, 0), 
                         1 : self.Get_outbinterval(DN, DN.minrecurt, 1)}
        
        
        t_OutBurst += outb_interval[init_outbtype]
        timeOuts = [t_OutBurst]
        typeOuts = [init_outbtype]
        

        # Now add outburst times after the initial one.
        # There are different cases for if the super outburst recurrence time
        #  is longer than twice the normal outburst recurrence time or shorter      
        if DN.P_sc > 2 * DN.P_c:
            nr_no_outb = 0
            if init_outbtype == 0:
                nr_no_outb = 1
             # step forward in time
            while t_OutBurst < 1.1 * max(tvals) / 8.64e4:
                outbtype = self.Get_outbtype(nr_no_outb, DN.nr_no_per_so)
                if outbtype == 0: 
                    nr_no_outb += 1
                elif outbtype == 1: 
                    nr_no_outb = 0
                t_OutBurst += outb_interval[outbtype]
                typeOuts.append(outbtype)
                timeOuts.append(t_OutBurst)
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
                t_OutBurst += outb_interval[outbtype]
                typeOuts.append(outbtype)
                timeOuts.append(t_OutBurst)
        else:   # The type of outburst switches every time
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
                timeOuts.append(t_OutBurst)
                        
                        
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
                Bool1 = tSamp < self.BP_t_e['normal']
                Bool2 = tSamp > self.BP_t_s['normal']
            if typeOuts[i] == 1:
                tscale = DN.D_so / self.BP_D_so
                tSamp = np.array(-tvals - tO*8.64e4) / tscale #Scale duration
                Bool1 = tSamp < self.BP_t_e['super']
                Bool2 = tSamp > self.BP_t_s['super']
            if any(Bool1 & Bool2): #i.e. the outburst is *in* the observations.
                if typeOuts[i] == 0:
                    these_bandMags = np.array(self.get_all_normal_bbs(tSamp, 
                        obBands, DN))
                if typeOuts[i] == 1:
                    these_bandMags = np.array(self.get_all_super_bbs(tSamp, 
                        obBands, DN))
                    
                i5 = (these_bandMags != Q_bandMags)
                bandMags[i5] = np.min([these_bandMags[i5],bandMags[i5]], 
                                      axis=0)

        # apply distance modulus
        dval = lc.LC_Dist
        dist_mod = 5.0 * np.log10(dval*1.0e3) - 5.0
        visible = False
        for i,label, in enumerate(obBands):
            bandMags[i] += dist_mod + lc.Extinction[label]
            Q_bandMags[i] += dist_mod + lc.Extinction[label]
            if bandMags[i] < threshold[label]: visible = True

        if visible :
            return [lc.LC_RA, lc.LC_DEC], bandMags, Q_bandMags
        else:
            self.N_trans -= 1
            return [],[], []
                        
    def Get_outbtype(self, nr_no_outb, nr_no_per_so):
        """ Function that returns the outburst type of the next outburst

            What type of outburst do we expect? if e.g. nr_no_per_so=3.2
             then we expect to have 3 normal outbursts and then a super
             outburst with a probability of 0.2. If there have already 
             been 4 normal outbursts in a row, we get a super outburst.
           
            0 = normal outburst
            1 = super outburst

            -------
            returns: The type of the next outburst
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

    def Get_outbinterval(self, DN, outb_interval, outbtype, D_ss = 0):
        """Calculates the interval in between two outbursts.

            The interval between the start of two outbursts differs for 
             a normal and a superoutburst. After all, the time between 
             outbursts is the same, but their duration is different.

            -------
            returns: the interval between two outbursts.
        """
        if DN.P_c < DN.P_sc:
            if outbtype == 1:
                new_outb_interval = outb_interval + (DN.D_so - DN.D_no)  
                                                  #Add duration of superoutburst
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
                new_outb_interval = outb_interval + (DN.D_no - DN.D_so)  
                                                 #Add duration of superoutburst
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

    def Get_outb_amp(self, lcdata_up, lcdata_lo, bands):
        """ A little function that finds the outburst amplitude for 
             every band in the blueprint light curve.

            It also calculates the quiescence magnitude in the V-band 
             of the light curve data.

            --------
            returns: the outburst amplitude in every band
                  and the quiescence mag in V
        """
        if 'V' not in lcdata_up.keys(): 
            raise TypeError("There is no V-band data in your UBVRI light"  
                            + " curve file for %s" % (self.tag))
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
        
        
class UGem_template(SUUMa_template):
    """
        The U Gem template is the same as for SU UMa.
        Just like with SU UMa there are two different types of outbursts:
          - The normal one
          - A longer one, but that isn't denoted as "super" in the outburst 
          because it is not brighter than its 'normal' counterpart.
    """
    def __init__(self, tag, trData, obRun):
        SUUMa_template.__init__(self, tag, trData, obRun)
        self.type = 'UGem'
        

class ZCam_template( SUUMa_template ):
    """ The template for Z Cam dwarf nova transients.

        Inherits the SU UMa transient template.
        The only addition is that of the standstills.
        Standstills are modelled with a rise into standstill, a decay 
         from standstill and a time of constant brightness in between.
    """
    def __init__(self, tag, trData, obRun):
        SUUMa_template.__init__(self, tag, trData, obRun)
        self.type = 'ZCam'
        self.broadbands_rise   = {}   #for stillstand
        self.broadbands_decay  = {}   #for stillstand

    def setUpLC(self, colorchoice):
        """ An interpolation function of the light curve data is made 
             for each filter that will be used in the observation.

            This is done for normal, super outbursts and standstills.
            The standstills are split into a rise into standstill and a
             decay from standstill.
            Later, this function can be sampled to obtain the transient 
             brightness at the time of observation.
            Also, the peak magnitude of the LC and the maximum 
             observation distance up to which a transient is still 
             visible are determined.
        """
        print('for %s' % self.tag)
        self.BP_outb_amp = {}
        self.BP_D_no = float(self.deltaT_LC[0].replace("(", ""))
        self.BP_D_so = float(self.deltaT_LC[1].replace(")", ""))
        self.BP_t_s  = {}
        self.BP_t_e  = {}
        for TYPE in ['normal', 'super', 'rise', 'decay']:
            print("file", self.LCFile + '_UBVRI_%s.hdf5' % (TYPE), 'r')
            lcdata_up = h5py.File(self.LCFile + '_UBVRI_%s.hdf5' % (TYPE), 'r')
            lcdata_lo = h5py.File(self.LCFile + '_%s_%s.hdf5' 
                                  % (colorchoice, TYPE), 'r')
            self.BP_outb_amp[TYPE], self.BP_Q_mag[TYPE] = self.Get_outb_amp(
                lcdata_up, lcdata_lo, np.append('V',self.bands))
            if not np.all(lcdata_up['times'][:] == lcdata_lo['times'][:]):
                print('Warning: the times array in the LC data for UBVRI is',\
                      ' not equal to %s' % (colorchoice))
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
                print("%s band" % band)
                if band.islower():
                    lcdata = lcdata_lo
                else:
                    lcdata = lcdata_up
                if not band in lcdata.keys():
                    print("no data for %s band" % band)
                    bb = emptyval * np.ones(len(tms))
                else:
                    #Only record relative change
                    bb = lcdata[band][:] - max(lcdata[band][:]) 
                ff = interp1d(tms, bb, bounds_error=False, fill_value=emptyval)
                this_Q_mag = max(lcdata[band][:]) - self.BP_Q_mag[TYPE]
                if band not in self.quiescentmag.keys():
                    self.quiescentmag[band] = this_Q_mag
                else:
                    self.quiescentmag[band] = max(self.quiescentmag[band], 
                                                  this_Q_mag)
                def fun(times, fcn=ff):
                    return fcn(times)
                if TYPE == 'normal':
                    self.broadbands_normal[band] = fun
                elif TYPE == 'super':
                    self.broadbands_super[band] = fun
                elif TYPE == 'rise':
                    self.broadbands_rise[band] = fun
                    self.BP_A_ss[band] = min(self.BP_A_ss[band], 
                        lcdata[band][-1] - max(lcdata[band][:]))
                elif TYPE == 'decay':
                    self.broadbands_decay[band] = fun
                    self.BP_A_ss[band] = min(self.BP_A_ss[band], 
                        lcdata[band][0] - max(lcdata[band][:]))
                #Now we need to determine the brightness maximum in each band
                #This is used in get_blueprints to determine whether the 
                # transient is visible or not
                self.peak_mags[band] = min(lcdata[band][:])
                
            lcdata_lo.close()
            lcdata_up.close()
        self.BP_d_rise = self.BP_t_e['rise'] - self.BP_t_s['rise'] #in s
        self.BP_d_decay = self.BP_t_e['decay'] - self.BP_t_s['decay'] #in s
        self.Max_obs_D = Maximum_observing_distance(self.mag_lim, 
            self.peak_mags, self.std_mag_R)

    def get_all_standstill_bbs(self, tsample, obbands, DN, D_ss, tO):
        """ At Time=tsample sample the standstill light curve function
             for band = obbands.

            A standstill comprises a rise, a plateau phase and a decay.
            The luminosity at the plateau phase is scaled similarly
             to the peak brightness of a normal outburst.
            In other color filters than 'V' the standstill is scaled such
             that the color at peak magnitude is maintained the same as
             for the 'standard' dwarf nova.
            What is added to the quiescent magnitude (i.e. the outburst)
             may never be negative (=> scale must be positive)
            
            returns: a numpy array with standstill magnitudes at t=tsample
        """
        bnds = []
        rise =  ( self.broadbands_rise['V'](tsample) * 
                       DN.A_no / self.BP_outb_amp['normal']['V']  
                       + DN.Q_mag['V'])
        lowvalue = (emptyval * DN.A_no / self.BP_outb_amp['normal']['V']  
                       + DN.Q_mag['V'])
        # Identify the end of the rise into standstill
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
        # d_till_decay is the duration of the standstill up till the decay
        #  phase in the blueprint ``time-frame''
        d_till_decay = self.BP_D_no / DN.D_no * D_ss + (endrise / 8.64e4)
        for i,bandLabel in enumerate(obbands):
            outb_amp_old = self.BP_outb_amp['normal']
            scale = ((DN.A_no + outb_amp_old[bandLabel] - outb_amp_old['V']) /
                     (outb_amp_old[bandLabel]))
            if scale < 0. : scale = 0
            this_bb = (self.broadbands_rise[bandLabel](tsample[i]) * scale
                       + DN.Q_mag[bandLabel])
            # Add the plateau phase luminosities.
            if tsample[i] > endrise and tsample[i] < d_till_decay * 8.64e4:
                this_bb = self.BP_A_ss[bandLabel] * scale + DN.Q_mag[bandLabel]
            if tsample[i] > d_till_decay * 8.64e4 :
                time = tsample[i] - d_till_decay * 8.64e4
                this_bb = (self.broadbands_decay[bandLabel](time) * scale
                           + DN.Q_mag[bandLabel])
            if this_bb > DN.Q_mag[bandLabel]:
                this_bb = DN.Q_mag[bandLabel]
            bnds.append(this_bb)
        return np.array(bnds)

    def sample_single_LC(self, lc, tvals, obBands, threshold):
        """ Sample the light curve of a single dwarf nova

            We first create a Dwarf Nova object.
            The initial outburst explosion time is noted and we add
             outbursts after and before this initial outburst.
            This is done in a specific manner: taking the difference
             between the normal outburst and superoutburst interval into
             account, we progressively add outburst times. After a 
             certain number of normal outbursts, a superoutburst occurs.
            However, when it is time for a standstill (the standstill
             recurrence time has passed), a standstill is generated in 
             its place.
            When a standstill has occured, a new standstill recurrence
             time and duration are generated.
            The duration of outbursts are scaled by scaling the time
             array that samples the outburst light curves.
            We construct an aggregate light curve. If outbursts overlap
             the brightest 'measurement' is taken.
            At last the distance modulus is applied plus dust extinction
             and the transient is tested for visibility.

            -------
            returns: A list of RA,DEC coordinates, magnitudes and 
                     quiescent magnitudes
        """
        DN = DwarfNova(self.type, self.P_c_s, self.P_c_loc, self.P_c_scale,
                       self.quiescentmag, self.bands)
        max2 = DN.Get_max_lookbacktime(D_rise = self.BP_d_rise / 8.64e4, 
                                       D_decay = self.BP_d_decay / 8.64e4)
        #Generate the time of the first outburst
        t_OutBurst = -max2 * np.random.random_sample() - max2   
        init_outbtype = 2
        #scale rise and decay in length
        DN.D_rise  = (self.BP_d_rise / 8.64e4) * DN.D_no / self.BP_D_no    
        DN.D_decay = (self.BP_d_decay / 8.64e4) * DN.D_no / self.BP_D_no  
        outblen_ss = DN.D_rise + DN.D_decay + DN.D_ss #Standstill duration

        # 0 = normal outburst
        # 1 = super outburst
        # 2 = standstill
        outb_interval = {0 : self.Get_outbinterval(DN, DN.minrecurt, 0), 
                         1 : self.Get_outbinterval(DN, DN.minrecurt, 1),
                         2 : self.Get_outbinterval(DN, DN.minrecurt, 2, 
                                                   D_ss = outblen_ss)}
        
        D_ss     = [DN.D_ss]
        DN.Change_standstill() #Generate a ss recurrence time and duration.
        next_ss_outb = t_OutBurst + DN.P_ss
        
        t_OutBurst += outb_interval[init_outbtype]
        timeOuts = [t_OutBurst]
        typeOuts = [init_outbtype]
        t_since_so = 0
        
        # Now add outburst times after the initial one.
        # There are different cases for if the super outburst recurrence time
        #  is longer than twice the normal outburst recurrence time or shorter 
        if DN.P_sc > 2 * DN.P_c:
            nr_no_outb = 0
            if init_outbtype == 0:
                nr_no_outb = 1
            # step forward in time
            while t_OutBurst < 1.1 * max(tvals) / 8.64e4:
                outbtype = self.Get_outbtype(nr_no_outb, DN.nr_no_per_so)
                if t_since_so > DN.P_sc: outbtype = 1
                if t_OutBurst + outb_interval[outbtype] > next_ss_outb:
                    outbtype = 2
                if outbtype in [1,2]:
                    t_since_so = 0
                if outbtype == 0:  
                    nr_no_outb += 1
                    t_since_so += outb_interval[outbtype]
                elif outbtype == 1: nr_no_outb = 0
                t_OutBurst += outb_interval[outbtype]
                D_ss.append(DN.D_ss)
                if outbtype == 2:
                    DN.Change_standstill()
                    outblen_ss = DN.D_rise + DN.D_decay + DN.D_ss
                    outb_interval[2] = self.Get_outbinterval(DN, DN.minrecurt, 
                        2, D_ss=outblen_ss)
                    next_ss_outb = t_OutBurst + DN.P_ss
                typeOuts.append(outbtype)
                timeOuts.append(t_OutBurst)
        elif DN.P_c < 2 * DN.P_sc:
            nr_so_outb = 0
            if init_outbtype == 1:
                nr_so_outb = 1
             # step forward in time
            while t_OutBurst < 1.1 * max(tvals) / 8.64e4:
                outbtype_so = self.Get_outbtype(nr_so_outb, DN.nr_so_per_no)
                if t_since_so > DN.P_sc: outbtype_so = 0
                if t_OutBurst + outb_interval[outbtype_so] > next_ss_outb:
                    outbtype_so = 2
                if outbtype_so in [0,2]:
                    t_since_so = 0
                if outbtype_so == 1: 
                    outbtype = 0
                elif outbtype_so == 0: 
                    outbtype = 1
                    nr_so_outb = 0
                    t_since_so += outb_interval[outbtype]
                elif outbtype_so == 2:
                    outbtype = 2
                    nr_so_outb = 0
                t_OutBurst += outb_interval[outbtype]
                D_ss.append(DN.D_ss)
                if outbtype == 2:
                    DN.Change_standstill()
                    outblen_ss = DN.D_rise + DN.D_decay + DN.D_ss
                    outb_interval[2] = self.Get_outbinterval(DN, DN.minrecurt, 
                        2, D_ss=outblen_ss)
                    next_ss_outb = t_OutBurst + DN.P_ss
                typeOuts.append(outbtype)
                timeOuts.append(t_OutBurst)
        else:   # The type of outburst switches every time
            while t_OutBurst < 1.1 * max(tvals) / 8.64e4:
                if typeOuts[-1] == 0:
                    outbtype = 1
                elif typeOuts[-1] == 1:
                    outbtype = 0
                else:
                    p = np.random.random_sample()
                    if p < 0.5: outbtype = 0
                    else: outbtype = 1
                if t_OutBurst + outb_interval[outbtype] > next_ss_outb:
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
                    outb_interval[2] = self.Get_outbinterval(DN, DN.minrecurt, 
                        2, D_ss=outblen_ss)
                    next_ss_outb = t_OutBurst + DN.P_ss
                typeOuts.append(outbtype)
                timeOuts.append(t_OutBurst)
        
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
                Bool1 = tSamp < self.BP_t_e['normal']
                Bool2 = tSamp > self.BP_t_s['normal']
            if typeOuts[i] == 1:
                tscale =  DN.D_so / self.BP_D_so
                tSamp = np.array(-tvals - tO*8.64e4) / tscale #Scale duration
                Bool1 = tSamp < self.BP_t_e['super']
                Bool2 = tSamp > self.BP_t_s['super']
            if typeOuts[i] == 2:
                tscale =  DN.D_no / self.BP_D_no
                tSamp = np.array(-tvals - tO * 8.64e4) / tscale #Scale duration
                #always sample standstills
                Bool1 = np.ones(len(tSamp), dtype=bool)   
                Bool2 = np.ones(len(tSamp), dtype=bool) 
            if any(Bool1 & Bool2): #i.e. the outburst is *in* the observations.
                if typeOuts[i] == 0:
                    these_bandMags = np.array(self.get_all_normal_bbs(tSamp, 
                        obBands, DN))
                if typeOuts[i] == 1:
                    these_bandMags = np.array(self.get_all_super_bbs(tSamp, 
                        obBands, DN))
                if typeOuts[i] == 2:
                    these_bandMags = np.array(self.get_all_standstill_bbs(
                        tSamp, obBands, DN, D_ss[i], tO))

                i5 = (these_bandMags != Q_bandMags)
                bandMags[i5] = np.min([these_bandMags[i5],bandMags[i5]], 
                                      axis=0)
        
        # apply distance modulus
        dval = lc.LC_Dist
        dist_mod = 5.0 * np.log10(dval*1.0e3) - 5.0
        visible = False
        for i,label, in enumerate(obBands):
            bandMags[i] += dist_mod + lc.Extinction[label]
            Q_bandMags[i] += dist_mod + lc.Extinction[label]
            if bandMags[i] < threshold[label]: visible = True
        
        if visible :
            return [lc.LC_RA, lc.LC_DEC], bandMags, Q_bandMags
        else:
            self.N_trans -= 1
            return [], [], []


class Mdwarf_template(transientTemplate):
    """ The template for flaring M-dwarfs

        Inherits the standard transient template and adds a few extra
         parameters. These are:
         -aval and bval: parameters for determination of activeness
         -alf_ac and bet_ac: Flare Frequency Distribution parameters 
          for active M-dwarfs 
         -alf_in and bet_in: Flare Frequency Distribution parameters 
          for inactive M-dwarfs  
         -da_dz: Flare Frequency Distribution slope parameter
         -logE_min_ac and logE_max_ac: boundaries of flare energy (in
          U-band) for active M-dwarfs
         -logE_min_in and logE_max_in: boundaries of flare energy (in
          U-band) for inactive M-dwarfs
         -Lq: The quiescent luminosity in each colorfilter
    """
    def __init__(self, tag, trData, obRun):
        transientTemplate.__init__(self, tag, trData, obRun)
        vals = np.genfromtxt(obRun.MdwarfFile, names=True ).T[tag]
        self.aval, self.bval, self.alf_ac, self.alf_in,\
            self.bet_ac, self.bet_in, self.da_dz = vals[0:7]
        self.Ebins         = 100	#Nr of energy bins
        
        logE_min_ac, logE_max_ac = vals[7:9]
        logE_min_in, logE_max_in = vals[9:11]
        qs                 = vals[11:]
        last_logE_ac       = (logE_max_ac - (logE_max_ac - logE_min_ac) 
                              / self.Ebins)
                             #last_logE is the the start of the last logE-bin
        last_logE_in       = (logE_max_in - (logE_max_in - logE_min_in) 
                              / self.Ebins)        
        self.dlogE_act     = np.linspace(logE_min_ac, last_logE_ac, 
                                         self.Ebins + 2)
        self.dlogE_ina     = np.linspace(logE_min_in, last_logE_in, 
                                         self.Ebins + 2)
        self.Lq            = {'U':qs[0], 'B':qs[1], 'V':qs[2], 'R':qs[3],
                              'I':qs[4], 'J':qs[5], 'H':qs[6], 'K':qs[7],
                              'u':qs[8], 'g':qs[9], 'r':qs[10], 'i':qs[11],
                              'z':qs[12], 'q':qs[13], 'y':qs[14]}
        self.Epk           = {}
        self.mag_resolution = obRun.mag_resolution
        self.broadbands_rise  = {}
        self.broadbands_decay = {}
        self.LC_time = {'decay': self.deltaT_LC * 24. * 3600.}

    def get_blueprints(self, c, dtObs):
        """ Create the blueprints for this type of transient

            For every transient that is generated, we need to generate 
             a time of outburst, draw a peak magnitude and test if we 
             can still see this transient given the position, dust 
             extinction and peak mag.
        """
        Nr_of_years = 1.	#Recurrent transient => Nr. of transients does not 
                         # depend on observation time
        N_transients = self.get_Ntransients(c, Nr_of_years)
        self.N_transnonsee += N_transients
        for i in range(N_transients):
            coords = c.sampleCoords()
            if coords != []:
                t_dummy = -365.25 * Nr_of_years * np.random.random()
                newLC = LC_blueprint(coords, self.peak_mags, 0.0, t_dummy)
                if self.canSee_LC(newLC):
                    _Coords = (newLC.LC_RA, newLC.LC_DEC, newLC.LC_Dist)
                    newLC.Extinction = c.Grid.Gal_dust.Sample_extinction(
                        *_Coords)
                    if self.canSee_LC(newLC):
                        self.transients.append(newLC)
                        self.N_trans += 1

    def setUpLC(self, colorchoice):
        """ An interpolation function of the light curve data is made 
             for each filter that will be used in the observation.

            An outburst consists of a rise and a decay.
            Later, this function can be sampled to obtain the transient 
             brightness at the time of observation.
            Also, the peak luminosity of the LC and the maximum 
             observation distance up to which a transient is still 
             visible are determined.
        """
        print('for %s' % self.tag)
        self.flux_0 = fc.flux_0[colorchoice]
        self.flux_0.update(fc.flux_0['UBVRI'])
        self.Pk_rise = {band:0 for band in np.append('U',self.bands)}
        self.Pk_decay = {band:0 for band in np.append('U',self.bands)}
        for TYPE in ['rise', 'decay']:
            lcdata_lo = h5py.File(self.LCFile + '_%s_%s.hdf5' 
                                  % (colorchoice, TYPE), 'r')
            lcdata_up= h5py.File(self.LCFile + '_UBVRI_%s.hdf5' % (TYPE), 'r')
            if not np.all(lcdata_up['times'][:] == lcdata_lo['times'][:]):
                print ('Warning: the times array in the LC data for UBVRI is',\
                       ' not equal to %s' % (colorchoice))
            tms = lcdata_up['times'][:]
            if TYPE == 'rise':self.broadbands_rise['t'] = tms;
            else: self.broadbands_decay['t'] = tms
            for band in np.append('U',self.bands):
                print("%s band" % band)
                if band.islower():
                    lcdata = lcdata_lo
                else:
                    lcdata = lcdata_up                
                if not band in lcdata.keys():
                    print("no data for %s band" % band)
                    bb = np.zeros(len(tms))
                else:
                    bb = lcdata[band][:]
                ff = interp1d(tms, bb, bounds_error=False, fill_value=0.0)
                def fun(times, fcn=ff):
                    return fcn(times)
                if TYPE == 'rise':
                    self.broadbands_rise[band] = fun
                    self.Pk_rise[band] = np.max(bb)
                    self.LC_time['rise'] = tms[-1]
                elif TYPE == 'decay':
                    self.broadbands_decay[band] = fun
                    self.Pk_decay[band] = np.max(bb)
                self.Epk[band] = max(self.Pk_rise[band], self.Pk_decay[band])
            if TYPE == 'rise':
                E_rise = np.trapz(lcdata_up['U'][:], lcdata_up['times'][:])
                self.U_rise = lcdata_up['U'][:]
            if TYPE == 'decay':
                E_decay = np.trapz(lcdata_up['U'][:], lcdata_up['times'][:])
                self.U_decay = lcdata_up['U'][:]
        self.E_LC = (E_rise + E_decay) * 640.4
        lcdata.close()
        self.Max_obs_D = self.Max_Mdwarf_distance()

    def Max_Mdwarf_distance(self):
        """ Calculates the maximum distance up to which an Mdwarf can be seen.

        This is done by looking at the maximum energy an Mdwarf can have in
         each filter band and taking the maximum over them.

            -------
            returns: The maximum distance (in pc)
        """
        max_energy = max(self.dlogE_act[-1], self.dlogE_ina[-1])
        scale      = 10**(max_energy) / self.E_LC
        max_lum = {band: (self.Epk[band] * scale + self.Lq[band]) 
                   for band in self.bands}
        maxdist = 0
        self.peak_mags = {}    #Empty dict to be filled in later
        #area in cm^2 of sphere of 10 pc (to calculate absolute magnitudes)
        area_abs_pc = 4 * np.pi * (10 * pc_to_cm) ** 2.    
        for band in self.bands:
            #Calculate the area of the sphere with the distance up to which we
            # can see the M dwarf in case of the largest outburst.            
            area = (max_lum[band] / self.flux_0[band]
                    * 10**(0.4 * self.mag_lim[band]))    #in cm^2
            dist = np.sqrt(area / (4 * np.pi)) / (pc_to_cm)
            if dist > maxdist: maxdist = dist     #in pc
            self.peak_mags[band] = -2.5 * np.log10(max_lum[band]/ area_abs_pc
                                                   / self.flux_0[band])
        return maxdist

    def get_all_bbs_rise(self, tsample, obbands):
        """ At Time=tsample sample the rise light curve function for
             band = obbands

            -------
            returns: a numpy array with super outburst magnitudes at t=tsample
        """
        bnds = []
        for i,bandLabel in enumerate(obbands):
            this_bb = self.broadbands_rise[bandLabel](tsample[i])
            bnds.append(this_bb)
        return np.array(bnds)

    def get_all_bbs_decay(self, tsample, obbands):
        """ At Time=tsample sample the decay light curve function for
             band = obbands

            -------
            returns: a numpy array with super outburst magnitudes at t=tsample
        """
        bnds = []
        for i,bandLabel in enumerate(obbands):
            this_bb = self.broadbands_decay[bandLabel](tsample[i])
            bnds.append(this_bb)
        return np.array(bnds)

    def sample_all_LCs(self, obTimes, obBands, threshold):
        """ For every transient, sample the light curve function at 
            obTimes in filters obBands.

           Only use the observations of the frames in which the 
            transient is visible.
           It is first determined whether the M-dwarf is active or not
           Then a Flare Frequency Distribution is determined.
           If there are more than 0 flares, their luminosity is sampled
            from the light curves and added on top of the quiescent lum
           If the brightness at the observation times meets the 
            threshold, keep it, otherwise throw it away.
        """
        maxtime = 1.e-99
        mintime = 1.e99
        for frame in obTimes:
            maxtime = max(maxtime, max(frame))
            mintime = min(mintime, min(frame))
        tWindow = (maxtime - mintime)/3600.0 
        for mdwarf in self.transients:
            visible = False
            active = False
            
            vis_frames = np.array(mdwarf.visibleframes)
            obT = np.concatenate(obTimes[vis_frames])
            obB = np.concatenate(obBands[vis_frames])
            
            # galactic height in parsecs
            hZ = get_galactic_height(mdwarf.LC_RA, mdwarf.LC_DEC, 
                                     mdwarf.LC_Dist) * 1.0e3 # parsecs
            P_A = self.aval * np.exp(-self.bval * abs(hZ))  
            ztest = np.random.random()
            if ztest <= P_A: 
                active = True
            if active: 
                alf = self.alf_ac + self.da_dz * abs(hZ) / 1.e3	
                                                         #da_dz is in kpc^-1
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
            cum_ne = (tWindow * np.power(10.0, alf) 
                      * np.power(10.0, bet * eBins_l))
            #cum_ne is cumulative. We need to decumulate it:
            ne = [cum_ne[i] - cum_ne[i+1] for i in range(self.Ebins +1)]
            nFlares = geometricP(ne)
            
            if max( nFlares ) < 1: 
                self.N_trans -=1
            else:
                lumlist, Qflux = self.get_all_lums(nFlares, eBins_m[:-1], 
                                                   obT, obB)
                maglist = np.zeros(len(obB))
                Qmaglist= np.zeros(len(obB))
                for j, lum in enumerate(lumlist):
                    bnd = obB[j]
                    area = 4.0 * np.pi * np.power(mdwarf.LC_Dist * 1e3 
                                                  * pc_to_cm, 2.0)  # cm^2
                    mag  = -2.5 * np.log10(lum / area / self.flux_0[bnd])
                    Qmag = -2.5 * np.log10(Qflux[j] / area / self.flux_0[bnd])
                    maglist[j]  = mag
                    Qmaglist[j] = Qmag
                    #the luminosity (in mags) of the quiescent star
                    if bnd in self.Lq.keys():
                        mag_non_outb = -2.5 * np.log10(self.Lq[bnd] / area 
                                                       / self.flux_0[bnd])
                    else:
                        print("Warning: the quiescent luminosity for band ", 
                              bnd,
                              " is not included in the file with M dwarf data")
                        mag_non_outb = 0
                    maglist[j]  += mdwarf.Extinction[bnd]
                    Qmaglist[j] += mdwarf.Extinction[bnd]
                    if (mag < threshold[bnd] and mag < mag_non_outb 
                            - self.mag_resolution):
                        visible = True
                if visible:
                    self.radec_list.append([mdwarf.LC_RA, mdwarf.LC_DEC])
                    self.bandmags_list.append(maglist)   
                    self.Qmag_list.append(Qmaglist)
                    self.bandbands_list.append(obB) 
                    self.obT_list.append(obT)
                else:
                    self.N_trans -= 1
        self.bandmags_list = np.array(self.bandmags_list)

    def get_all_lums(self, nFlares, eBins, obT, obbands):
        """ Calculate the luminosities of all M-dwarf flares for this
            single M-dwarf.

            For every flare (in nF) draws a rise time and decay time.
            Outbursts are generated by scaling the flare light curve
             and energy.
            If multiple flares occur simultaneously their fluxes are
             added to each other. The quiescent flux is added on top of 
             this.

            -------
            returns: a list of observations of the luminosity of this
                     M-dwarf plus a list of quiescent luminosities of 
                     equal length.
        """
        window  = obT[-1] - obT[0]
        lumlist = np.zeros(len(obT))
        Qflux   = emptyval * np.ones(len(obT)) #Quiescent flux
        eF, nF  = eBins[nFlares > 0], nFlares[nFlares > 0]
        for i, n in enumerate(nF):
            logEU = eF[i]
            trise  = np.power(10.0, np.random.normal(0.31 * logEU - 7.6, 0.32))
            tdecay = np.power(10.0, np.random.normal(0.51 * logEU - 13., 0.22))

            thisE = np.power(10.0, logEU)
            dts = np.random.random(n) * window
            #generate outbursts that started after t0
            for j,dt in enumerate(dts):
                tSamp = obT - dt
                tr_scale = trise / self.LC_time['rise'] 
                td_scale = tdecay / self.LC_time['decay']
                Energy_scale = self.get_Energy_scale(thisE, tr_scale, td_scale)
                tSamp_rise  = tSamp / tr_scale
                tSamp_decay = (tSamp - trise) / td_scale
                new_flux_rise  = (Energy_scale
                                  * self.get_all_bbs_rise(tSamp_rise, obbands))
                new_flux_decay = (Energy_scale 
                                  * self.get_all_bbs_decay(tSamp_decay, 
                                                            obbands))
                lumlist += new_flux_rise + new_flux_decay
            dts = np.random.random(n) * window
            dts = dts[dts < (trise + tdecay) * 4. * 24. * 60. * 60.]
            for j,dt in enumerate(dts):
                tSamp = obT + dt
                tr_scale = trise / self.LC_time['rise'] 
                td_scale = tdecay / self.LC_time['decay']
                Energy_scale = self.get_Energy_scale(thisE, tr_scale, td_scale)
                tSamp_rise  = tSamp / tr_scale
                tSamp_decay = (tSamp - trise) / td_scale
                new_flux_rise  = (Energy_scale
                                  * self.get_all_bbs_rise(tSamp_rise, obbands))
                new_flux_decay = (Energy_scale 
                                  * self.get_all_bbs_decay(tSamp_decay, 
                                                            obbands))
                lumlist += new_flux_rise + new_flux_decay
        for i,bnd in enumerate(obbands):   #Add quiescence flux
            if bnd in self.Lq.keys():
                lumlist[i] += self.Lq[bnd]
                Qflux[i] = self.Lq[bnd]
        return lumlist, Qflux

    def get_Energy_scale(self, flare_E, tr_scale, td_scale):
        """ Calculates the scale with which the flux of the template
            flare has to be adjusted for this flare

            -------
            returns: The scale with which the light curve has to be 
                     scaled.
        """
        Time_rise   = self.broadbands_rise['t']
        Time_decay  = self.broadbands_decay['t']
        
        #Find the energy change from the time adjustments alone
        E_new_rise  = np.trapz(self.U_rise, Time_rise * tr_scale)
        E_new_decay = np.trapz(self.U_decay, Time_decay * td_scale)
        E_new_inter = (E_new_rise + E_new_decay) * 640.4

        return flare_E / E_new_inter
                

class xgal_template(transientTemplate ):
    """ The template for extragalactic transients

        Inherits the standard transient template and adds a few extra
         parameters. They are:
         - NM: The number of transients per formed solar mass
         - alpha: The Delay Time Distribution parameter.
         - HGE_param: The host galaxy extinction parameter
         - Kcormodel: The prefix of Kcorrection file names
         - Host_ext: A host extinction object
         - scale: The Phillips relation scale (preliminarily set to 1)
         - scale_0: The Phillips relation scale of the blueprint LC
         - M_B_max: Maximum brightness in B-band of the blueprint LC
    """
    def __init__(self, tag, trData, obRun):
        transientTemplate.__init__(self, tag, trData, obRun)
        self.NM         = float(trData[6])
        self.alpha      = float(trData[13])
        self.HGE_param  = trData[14]
        self.Kcormodel  = trData[5]
        self.colorScheme= obRun.colorScheme
        self.Kcorfuns   = cm.Kcor(self.Kcormodel, self.colorScheme, self.bands)
        self.Host_ext   = dust.Host_extinction(self.HGE_param , obRun)
        #The following are only applicable to SNIa
        self.scale      = 1		
        self.scale_0    = 1.09
        self.M_B_max    = -18.89
        
    def setUpLC(self, colorchoice):
        """ An interpolation function of the light curve data is made 
             for each filter that will be used in the observation.

            Later, this function can be sampled to obtain the transient 
             brightness at the time of observation.
            Also, the peak magnitude of the LC and the maximum 
             observation distance up to which a transient is still 
             visible are determined.
        """
        print('for %s' % self.tag)
        lcdata_up = h5py.File(self.LCFile + '_UBVRI.hdf5', 'r')
        lcdata_lo = h5py.File(self.LCFile + '_%s.hdf5' % colorchoice, 'r')
        tms = lcdata_up['times'][:]    #times(UBVRI) = times(ugriz)
        self.broadbands['t'] = tms;
        self.peak_mags = {}   #Empty dict to be filled in later
        Kcorfile_up = h5py.File('LightCurveFiles/Kcorrections/%s_UBVRI.hdf5' 
                                % (self.Kcormodel),'r')
        Kcorfile_lo = h5py.File('LightCurveFiles/Kcorrections/%s_%s.hdf5' 
                                % (self.Kcormodel, self.colorScheme),'r')
        for band in self.bands:
            print("%s band" % band)
            if band.islower():
                lcdata   = lcdata_lo
                Kcorfile = Kcorfile_lo
            else:
                lcdata   = lcdata_up
                Kcorfile = Kcorfile_up
            if not band in lcdata.keys():
                print("no data for %s band" % band)
                bb = emptyval * np.ones(len(tms))
            else:
                scale = self.peak_mag_R - min(lcdata_up['R'][:])
                bb = lcdata[band][:] + scale
            ff = interp1d(tms, bb, bounds_error=False, fill_value=emptyval)
            def fun(times, fcn=ff):
                return fcn(times)
            self.broadbands[band] = fun
            #Now we need to determine the brightness maximum in each color band
            #This is used in get_blueprints to determine whether the transient
            # is visible or not
            #We need the maximum over ALL redshifts!
            minimum_band = emptyval
            for z in range(len(Kcorfile[band][0,:])):
                new_minimum = min(lcdata[band][:] + Kcorfile[band][:,z])
                minimum_band = min(minimum_band, new_minimum)
            self.peak_mags[band] = minimum_band + scale
        lcdata_lo.close()
        lcdata_up.close()

    def Phillips_scale(self, magDev):
        """ The Phillips peak-luminosity-broadness relation.

            It is only applicable to supernovae 1a
        """
        Mag_zero = self.M_B_max + magDev - 5 * np.log10(cm.H0) + 25
        scale = 1. / (0.59 * Mag_zero + 3.0)	#The real scale factor
        #Adjust for scale factor of template LC
        self.scale = self.scale_0 / scale 

    def get_all_bbs(self, tsample , obbands, redshift):
        """ At Time=tsample sample the light curve function for band = obbands

            Also takes K-corrections into account.

            -------
            returns: The magnitudes of this outburst
        """
        bnds = []
        for i, bandLabel in enumerate(obbands):
            Kcor = self.Kcorfuns.Sample_Kcor(bandLabel, 
                                             tsample[i] * self.scale, redshift)
            this_bb = (self.broadbands[bandLabel](tsample[i] * self.scale) 
                       + float(Kcor))
            bnds.append(float(this_bb))
        return np.array(bnds)

    def get_Ntransients(self, cell, Nr_of_years):
        """ Generate a certain number of transients in the given cell.

            This depends on the transient density, the location of the 
             cell and the number of years the simulation spans (many 
              densities are given in N per year)

            -------
            returns: The number of transients in this cell
        """
        kpc3 = cell.vol
        N_per_y = kpc3 * cell.Grid.Cosmo.Nr_density_xtr(cell.z) 
        N_per_y = N_per_y * Nr_of_years
        if N_per_y > 2.0:
            NTr = int(N_per_y)
        else: NTr = geometricP(N_per_y)
        return NTr

    def get_blueprints(self, c, dtObs):
        """ Create the blueprints for this type of transient

            For every transient that is generated, we need to generate 
             a time of outburst, redshift, draw a peak magnitude and 
             test if we can still see this transient given the position, 
             dust extinction and peak mag.
        """
        #At first observation the first LC should have just ended
        Nr_of_years = (self.deltaT_LC + dtObs) / 365.25   
        N_transients = self.get_Ntransients( c, Nr_of_years )
        self.N_transnonsee += N_transients
        for i in range( N_transients ):
            coords = c.sampleCoords()
            if coords != []:
                tOutburst = -365.25 * Nr_of_years * np.random.random() 
                if tOutburst + self.deltaT_LC + dtObs > 0.0:
                    magDelta = np.random.normal(scale=self.std_mag_R)
                    peak_M = {band: self.peak_mags[band] + magDelta 
                              for band in self.bands}
                    newLC = LC_blueprint(coords, peak_M, magDelta, tOutburst)
                    if self.canSee_LC(newLC):
                        _Coords = (newLC.LC_RA, newLC.LC_DEC, newLC.LC_Dist)
                        Host_ext = self.Host_ext.Sample_host_extinction()
                        MW_ext   = c.Grid.Xtr_dust.Sample_extinction(*_Coords)
                        newLC.Extinction = {band: MW_ext[band] + Host_ext[band] 
                                            for band in MW_ext.keys()}
                        if self.canSee_LC(newLC):
                            newLC.Redshift = c.Grid.Cosmo.get_redshift(
                                newLC.LC_Dist / 1.e3)	#convert to Mpc
                            self.transients.append(newLC)
                            self.N_trans += 1

    def sample_all_LCs(self, obTimes, obBands, threshold):
        """ For every transient, sample the light curve function at 
            obTimes in filters obBands.

           Only use the observations of the frames in which the 
            transient is visible.
           If the transients are supernovae Ia: apply the Phillips
            relations.
           If the brightness at the observation times meets the 
            threshold, keep it, otherwise throw it away.
        """
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
            bandMags = self.get_all_bbs(tSamp, obB, lc.Redshift)
            for i,label in enumerate(obB):
                bandMags[i] += dist_mod + lc.Extinction[label] + magDev
                if bandMags[i] < threshold[label]: visible = True
            if visible: 
                self.radec_list.append([lc.LC_RA, lc.LC_DEC])
                self.bandmags_list.append(bandMags)
                self.Qmag_list.append(emptyval * np.ones(len(bandMags)))
                self.bandbands_list.append(obB)
                self.obT_list.append(obT)
            else:
                self.N_trans -=1
        self.bandmags_list = np.array(self.bandmags_list)
      
  
class xgal_no_Kcor_template(transientTemplate):
    """ The template for extragalactic transients without K-corrections

        This template is meant for transients that either do not have 
         any K-corrections or are so close by (and fainter than 
         supernovae) that they are only found at very low redshifts.
        The volumetric density therefore also does not need to scale 
         with redshift.
    """
    def __init__(self, tag, trData, obRun):
        transientTemplate.__init__(self, tag, trData, obRun)
        self.stellarNrm = float(trData[6])
        self.HGE_param  = trData[14]    #Host galaxy extinction parameter
        self.Host_ext   = dust.Host_extinction(self.HGE_param , obRun)
        self.colorScheme= obRun.colorScheme
        #The following are only applicable to SNIa
        self.scale      = 1		
        self.scale_0    = 1.09
        self.M_B_max    = -18.89
        
    def Phillips_scale(self, magDev):
        """ The Phillips peak-luminosity-broadness relation.

            Is only applicable to supernovae 1a
        """
        Mag_zero = self.M_B_max + magDev - 5 * np.log10(cm.H0) + 25
        scale = 1./(0.59 * Mag_zero + 3.0)	#The real scale factor
        self.scale = self.scale_0 / scale	#Adjusted for scale factor of template LC

    def get_Ntransients(self, cell, Nr_of_years):
        """ Generate a certain number of transients in the given cell.

            This depends on the transient density, the location of the 
             cell and the number of years the simulation spans (many 
             densities are given in N per year)

            -------
            returns: The number of transients in this cell
        """
        kpc3 = cell.vol
        N_per_y = kpc3 * self.stellarNrm * Nr_of_years
        if N_per_y > 2.0:
            NTr = int(N_per_y)
        else: NTr = geometricP(N_per_y)
        return NTr

    def get_blueprints(self, c, dtObs):
        """ Create the blueprints for this type of transient

            For every transient that is generated, we need to generate 
             a time of outburst, draw a peak magnitude and test if we 
             can still see this transient given the position, (host
             galaxy) dust extinction and peak mag.
        """
        #At first observation the first LC should have just ended
        Nr_of_years = (self.deltaT_LC + dtObs) / 365.25   
        N_transients = self.get_Ntransients(c, Nr_of_years)
        self.N_transnonsee += N_transients
        for i in range(N_transients):
            coords = c.sampleCoords()
            if coords != []:
                tOutburst = -365.25 * Nr_of_years * np.random.random() 
                if tOutburst + self.deltaT_LC + dtObs > 0.0:
                    magDelta = np.random.normal(scale=self.std_mag_R)
                    peak_M = {band: self.peak_mags[band] + magDelta 
                              for band in self.bands}
                    newLC = LC_blueprint( coords, peak_M, magDelta, tOutburst)
                    if self.canSee_LC(newLC):
                        _Coords = (newLC.LC_RA, newLC.LC_DEC, newLC.LC_Dist)
                        Host_ext = self.Host_ext.Sample_host_extinction()
                        MW_ext   = c.Grid.Xtr_dust.Sample_extinction(*_Coords)
                        newLC.Extinction = {band: MW_ext[band] + Host_ext[band] 
                                            for band in MW_ext.keys()}
                        if self.canSee_LC(newLC):
                            self.transients.append(newLC)
                            self.N_trans += 1

    def sample_all_LCs(self, obTimes, obBands, threshold):
        """ For every transient, sample the light curve function at 
            obTimes in filters obBands.

           Only use the observations of the frames in which the 
            transient is visible.
           If the brightness at the observation times meets the 
            threshold, keep it, otherwise throw it away.
        """
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
            bandMags = self.get_all_bbs(tSamp, obB)
            for i,label in enumerate(obB):
                bandMags[i] += dist_mod + lc.Extinction[label] + magDev
                if bandMags[i] < threshold[label]: 
                    visible = True
            if visible: 
                self.radec_list.append([lc.LC_RA, lc.LC_DEC])
                self.bandmags_list.append(bandMags)
                self.Qmag_list.append(emptyval * np.ones(len(bandMags)))
                self.bandbands_list.append(obB)
                self.obT_list.append(obT)
            else:
                self.N_trans -=1
        self.bandmags_list = np.array(self.bandmags_list)


class kilonovaTemplate(transientTemplate):
    """ The template for kilonovatransients.
        
        This template basically reads the parameters from params.py and
         the parameters that belong the observation frame.
    """
    def __init__(self, tMerge, dist, mvRed, mvBlue, obRun, Framenr):
        grid            = obRun.SkyGrids[Framenr]
        self.tag        = 'kilonova'
        self.bands      = obRun.bands
        self.t_merge    = tMerge
        if dist != None:   #Set distance of kilonova
            self.dist = dist
        else: self.dist = sampleLIGODistance(obRun.maxLIGODist)
        # Generate RA/DEC coordinates within the frame
        self.DEC        = sampleDEC(grid.DEC_lo, grid.DEC_hi)
        RA_lo           = (grid.RA_c - (0.5 * obRun.aperture_RA 
                           / np.cos(self.DEC * (np.pi/180.))) )
        RA_hi           = (grid.RA_c - (0.5 * obRun.aperture_RA 
                           / np.cos(self.DEC * (np.pi/180.))) )
        self.RA         = sampleRA(RA_lo, RA_hi)
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
                                         0., 1., self.t_merge)] 
                                         #magPk and MagDelta don't matter here.

    def sample_all_LCs(self, obTimes, obBands, threshold):
        """ For the kilonova sample the light curve function at obTimes
             in filters obBands.

           Only use the observations of the frames in which the 
            transient is visible.
           If only the red component is enabled we sample the red light 
            curve and if only the blue component is enabled we sample 
            the blue light curve. If both are enabled we combine their
            light curves by adding fluxes and then sample it.
        """
        vis_frames = np.array(self.transients[0].visibleframes)

        obT = np.concatenate(obTimes[vis_frames])
        obB = np.concatenate(obBands[vis_frames])
        obT_f = max(obT)
        dts = obT_f - obT
        bandmags_list = np.ones(len(dts)) * emptyval
        
        if self.mvRed != None:
            mred, vred = get_best_match('red', *self.mvRed)
            rlcf_up = ('kilonovaLib/knRed_m%.4f_vk%.2f_Xlan1e-1.0_UBVRI.h5' 
                       % (mred, vred))
            rlcf_lo = ('kilonovaLib/knRed_m%.4f_vk%.2f_Xlan1e-1.0_%s.h5' 
                       % (mred, vred, self.colorchoice))
            redLC_up = h5py.File(rlcf_up, 'r' )
            redLC_lo = h5py.File(rlcf_lo, 'r' )
            if self.mvBlue == None: 
                self.netLC['times'] = redLC_up['times'][:]
                for band in self.bands:
                    if band.islower():
                        self.netLC[band] = redLC_lo[band][:]
                    else:
                        self.netLC[band] = redLC_up[band][:]
                redLC_up.close()
                redLC_lo.close()
        if self.mvBlue != None:
            mblu, vblu = get_best_match('blue', *self.mvBlue)
            blcf_up = ('kilonovaLib/knBlu_m%.4f_vk%.2f_Xlan1e-5.0_UBVRI.h5' 
                       % (mblu, vblu))
            blcf_lo = ('kilonovaLib/knBlu_m%.4f_vk%.2f_Xlan1e-5.0_%s.h5' 
                       % (mblu, vblu, self.colorchoice))
            bluLC_up = h5py.File(blcf_up, 'r')
            bluLC_lo = h5py.File(blcf_lo, 'r')
            if self.mvRed == None:
                self.netLC['times'] = bluLC_up['times'][:]
                for band in self.bands:
                    if band.islower():
                        self.netLC[band] = bluLC_lo[band][:]
                    else:
                        self.netLC[band] = bluLC_up[band][:]
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
                m12 = -2.5 * np.log10(np.power(10.0, -.4 * m1) 
                                      + np.power(10.0, -.4 * m2))
                self.netLC[band] = m12
            redLC_lo.close()
            redLC_up.close()
            bluLC_lo.close()
            bluLC_up.close()
        # Now that we have the light curve, we first interpolate it
        tSample = (obT - self.t_merge)
        dist_mod = 5.0 * np.log10(self.dist / pc_to_cm) - 5.0
        for j, band in enumerate(self.bands):
            bandLC = self.netLC[band]
            self.f_band[band] = interp1d(
                self.netLC['times'][np.isinf(bandLC) == False],
                bandLC[np.isinf(bandLC) == False], 
                bounds_error=False, 
                fill_value=np.inf 
                )
        # And then sample it at the appropriate times + add distance modulus
        visible = False
        for j,band in enumerate(obB):
            bandmags_list[j] = self.f_band[band](tSample[j]) + dist_mod
            if bandmags_list[j] < threshold[band]:
                visible = True
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
    """ An instance of this class is a blueprint transient

        tExplode is the (random) time of (first) explosion
        LC_RA is the RA coordinate of the transient
        LC_DEC is the DEC coordinate of the transient
        LC_Dist is the distance coordinate of the transient
        PeakMag is the peak magnitude of this transient
        magDelta is the deviation of the peak magnitude from the 
         mean of this transient type
        Extinction is the dust extinction in every band
        Redshift is the redshift of the transient
        visibleframes is a list of framenumbers in which the transient
         is visible
        
    """
    def __init__ (self, coords, magPk, magDelta, tBurst):
        self.tExplode   = tBurst
        self.LC_RA      = coords[0]
        self.LC_DEC     = coords[1]
        self.LC_Dist    = coords[2]
        self.PeakMag    = magPk
        self.magDelta   = magDelta
        self.Extinction = {}
        self.Redshift   = 0	#preliminarily
        self.visibleframes = []


class DwarfNova:
    """ A Dwarf nova class containing all of the relevant parameters.

        All of this is based on either the article or the database of 
         Otulakowska-Hypka (2016)
        P_c is the recurrence time of normal outbursts
        P_sc is the recurrence time of super/longer outbursts
        P_ss is the recurrence time for standstills
        A_no is the amplitude of normal outbursts
        A_so is the amplutide of super/longer outbursts
        D_no is the duration of normal outbursts
        D_so is the duration of super/longer outbursts
        D_ss is the duration of the standstill minus rise and decay 
         (which scale with D_no)
        no_per_so is the average number of normal outbursts per 
         super/long outbursts
    """
    def __init__ (self, Type, P_c_s, P_c_loc, P_c_scale, q_mag, bands):
        self.type = Type
        self.P_c = lognorm.rvs(P_c_s, P_c_loc, P_c_scale) #days
        mu1       = 1.32 + 0.66 * np.log10(self.P_c)
        self.P_sc = 10**(np.random.normal(loc=mu1, scale=0.24))   #days
        self.P_ss = (365.*5. - 10.) * np.random.random_sample() + 10.  #days
        self.A_no = lognorm.rvs(0.32, -0.89, 3.84)
        mu2       = 1.043 * self.A_no + 0.68
        self.A_so = np.random.normal(loc=mu2, scale=0.41)    #mags
        if self.A_so < 0 : self.A_so = 0
        if Type == 'SUUMa':
            self.D_no = lognorm.rvs(0.49, 0.36, 3.35)   #days
        elif Type == 'UGem':
            self.D_no = lognorm.rvs(0.74, 1.83, 6.60)   #days
        elif Type == 'ZCam':
            self.D_no = lognorm.rvs(0.49, 0.35, 8.31)   #days
        else:
            print("Warning: Class DwarfNova in transient.py did not receive",\
                  " a correct type. Using U Gem as default.")
            self.D_no = lognorm.rvs(0.74, 1.83, 6.60)   #days
        mu3       = 2.09 * self.A_so + 9.7
        self.D_so = np.random.normal(loc=mu3, scale=5.7)
        if self.D_so < 0 : self.D_so = 0
        self.D_ss = 365. * np.random.random_sample()  #days
        self.D_rise = 1   #to be changed later
        self.D_decay = 1  #to be changed later

        if self.P_sc > self.P_c:
            self.nr_no_per_so = self.P_sc / self.P_c - 1
            self.nr_so_per_no = 1. / self.nr_no_per_so
        else:
            self.nr_so_per_no = self.P_c / self.P_sc - 1
            self.nr_no_per_so = 1. / self.nr_so_per_no
            
        qmag_mu = 6.33 + 1.64 * np.log10(min(self.P_c, self.P_sc))  
        self.quiescent_mag = np.random.normal(loc=qmag_mu, scale=0.3)#in V-band
        self.Q_mag = {band: q_mag[band] + self.quiescent_mag 
                      for band in np.append('V',bands)}
        
        self.maxoutblen = max(self.D_no, self.D_so)
        self.maxrecurt  = max(self.P_c, self.P_sc)
        self.minrecurt  = min(self.P_c, self.P_sc)

    def Change_standstill( self ):
        """ Standstills change in duration and frequency so much more 
             than other outburst types, that we need to change them 
             after every standstill. This is done in this function.
        """
        self.D_ss = (365. - 10.) * np.random.random_sample() + 10.  #days
        self.P_ss = (365. * 5. - 10.) * np.random.random_sample() + 10.  #days
        
    def Get_max_lookbacktime( self, D_rise = 0, D_decay = 0 ):
        """ Calculates the maximum time we have to look back in to 
             sample the first outburst time
        """
        if self.type != 'ZCam':
            return self.maxoutblen + self.maxrecurt
        else:
            if self.P_ss > self.maxrecurt:
                return self.P_ss + self.D_ss + D_rise + D_decay
            else:
                return self.maxoutblen + self.maxrecurt
            
       
def Maximum_observing_distance(mag_lim, pkmag, devmag):
    """ Calculates up to what distance we should generate transients. 

        This is done by calculating the distance of an object with the
         peak magnitude plus a 3 sigma standard deviation.
    """
    Pkmag = np.array([pkmag[band] for band in pkmag.keys()]).astype(float)
    mag_lim = np.array([mag_lim[band] for band in pkmag.keys()]).astype(float)
    exp_term = 0.2 * (max(mag_lim - (Pkmag - 3.0 * devmag))) + 1.0
    Max_D = np.power(10.0, exp_term) # parsecs
    return Max_D

def get_best_match(br, m_msun, v_c):
    """ Get the best match kilonova light curve files for:

        br: blue or red: the color of the component 
        m_msun: the ejecta mass in units of M_Sun
        v_c: the ejecta velocity in units of light speed
        This function picks the best matching (m,v) out of the
         possibilities given by all the files present in the
         kilonovaLib folder.

        -------
        returns: the best matching mass and velocity
    """
    ms, vs = order_mv()
    m_match = ms[np.argmin(np.abs(ms - m_msun))]
    v_match = vs[np.argmin(np.abs(vs - v_c))]
    print('Best match for %s kilonova part: M = %.2e Msun, v = %.2f c.' 
           % (br, m_match, v_match) )
    return m_match, v_match
    
def order_mv():
    """ Finds all the possible masses and velocities of the kilonova 
         ejecta as posed by the filenames in the kilonovaLib folder.

        -------
        returns: All possible masses and velocities (in solar masses
                 and light speed)
    """
    ms, vs = [], []
    for f in (glob.glob('kilonovaLib/knBlu*.h5') 
              + glob.glob('kilonovaLib/knRed*.h5') ):
        this_m = f[f.find('_m')+2 : f.find('_v')]
        this_v = f[f.find('_v')+3 : f.find('_X')]
        if this_m not in ms: 
            ms.append(this_m)
        if this_v not in vs: 
            vs.append(this_v)
    ms = np.sort(np.array([float(m) for m in ms]))
    vs = np.sort(np.array([float(v) for v in vs]))
    return ms, vs
    
def sampleRA(ra_lo, ra_hi):
    """ Sample an RA coordinate within the frame's boundaries

        -------
        returns: RA coordinate in degrees
    """
    if ra_lo > ra_hi: ra_lo, ra_hi = ra_hi, ra_lo
    RA_sampled = np.random.random() * (ra_hi - ra_lo) + ra_lo
    return RA_sampled

def sampleDEC(dec_lo, dec_hi):
    """ Sample an DEC coordinate within the frame's boundaries


        -------
        returns: DEC coordinate in degrees
    """
    if dec_lo > dec_hi: dec_lo, dec_hi = dec_hi, dec_lo
    DEC_sampled = np.random.random() * (dec_hi - dec_lo) + dec_lo
    return DEC_sampled

def sampleLIGODistance(  dmax ):
    """ Sample a distance within the LIGO boundaries

        -------
        returns: a distance in centimeters
    """
    while True:
        x = np.random.random()*dmax
        y = np.random.random()*dmax
        z = np.random.random()*dmax
        r = np.sqrt( x*x + y*y + z*z )
        if r <= dmax: break
    return r

