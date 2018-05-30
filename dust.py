import numpy as np
import params as pm
import json, requests
import h5py
import math
from astLib import astCoords
from astropy import units as u
from itertools import chain
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import RegularGridInterpolator

#https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/LCYHJG
#source of dustmaps-package fro Green et al.
#http://dustmaps.readthedocs.io/en/latest/_modules/dustmaps/bayestar.html#BayestarQuery


RV_BV = {'U':4.334, 'B':3.626, 'V':2.742, 'R':2.169, 'I':1.505, 'J':0.764, 'H':0.449, 'K':0.302, 'u':4.239, 'g':3.303, 'r':2.285, 'i':1.698, 'z':1.263}
#3.1 at http://iopscience.iop.org/article/10.1088/0004-637X/737/2/103#apj398709t6
RV_JK = { key: RV_BV[key] * 1.748 for key in RV_BV.keys() }
emptyval = 1000.0

class Green_Extinction:
    def __init__(self, Xtr_dust, RA_lo, RA_hi, DEC_lo, DEC_hi):
        self.Name        = "Green et al. (2017)"
        self.Ang_Res     = 0.05667	#deg (=3.4 arcmin)
        self.bands       = list( pm.showbands )
        self.RA_lo       = RA_lo
        self.RA_hi       = RA_hi
        self.DEC_lo      = DEC_lo
        self.DEC_hi      = DEC_hi
        self.queried     = False
        self.f           = self.Setup_dust_grid()
        self.Xtr_dust    = Xtr_dust

    def Sample_extinction(self, ra, dec, D):
        """
        Sample the 3D dust grid to obtain the EBV at
         the location of the transient
        This is then converted to an extinction
        The result might be either emptyval or 'nan'. This means that
         the result is out of the interpolation/query range. 
         In that case we switch to the extragalactic EBV
        Input:
        ra: RA-coordinate in degrees
        dec: DEC-coordinate in degrees
        D: distance to transient in kpc
        """
        dmod = 5. * np.log10(D * 1000.) - 5.
        A = {}
        EBV = self.f( [ra, dec, dmod] )
        if EBV == emptyval or math.isnan(EBV):	#Out of intpol range
            A = self.Xtr_dust.Sample_extinction(ra, dec)
        else:
            for color in self.bands:
                A[color] = float(EBV) * RV_BV[color] 
        return A   

    def Setup_dust_grid(self):
        """
        Here we set up the 3D dust grid. We use the 
         RegularGridInterpolator because it's the fasted
        Our grid is regular (and rectangular).
        To also have data just across the border, we add 1 RA and DEC
         and add grid points 0.5*Ang_Res across this border
        We add 1 extra RA-coord to have an angular resolution slightly
         smaller than angres
        We then query the Green dust map with the 'full' mode (which
         is actually faster than the 'lite' mode)
        Subsequently we convert the 2D grid of Bestfit into a 3D grid
        We add data points for a distance modulus of zero (we assume 
         the EBV at d=0 to be zero)
        And finally we Interpolate this grid with linear interpolation
        If one samples a coordinate outside the interpolation 
         boundaries, it will return the value -100 
        """

        Nr_RA = np.ceil((self.RA_hi - self.RA_lo) / self.Ang_Res) + 2
        Nr_DEC = np.ceil((self.DEC_hi - self.DEC_lo) / self.Ang_Res) + 2

        ra  = np.linspace(self.RA_lo  - 0.5 * self.Ang_Res, 
                          self.RA_hi  + 0.5 * self.Ang_Res, Nr_RA)
        dec = np.linspace(self.DEC_lo - 0.5 * self.Ang_Res, 
                          self.DEC_hi + 0.5 * self.Ang_Res, Nr_DEC)
        RA, DEC = np.meshgrid(ra,dec)
        RA = list(chain.from_iterable(RA))
        DEC = list(chain.from_iterable(DEC))

        DQ = Greendustquery(RA,DEC, mode='full')	#Query Argonaut

        D = [0]		#Distance modulus zero isn't included in query
        D.extend(DQ['distmod'])
        Dlen = len(D)
        Bestfit = DQ['best']

        EBV = np.zeros((len(ra), len(dec), Dlen))
        for i in range( len(dec) ):
            for j in range( len(ra) ):
                Onecoord = [0]			#The EBV at d=0 at one coordinate
                Onecoord.extend(Bestfit[i+j])	#extend with EBV at other ds
                #print Onecoord
                EBV[j,i,:] = Onecoord

        self.queried = True

        return RegularGridInterpolator((ra,dec,D), EBV,
                                       method='linear',  
                                       bounds_error = False, 
                                       fill_value = emptyval)




class Schlegel_Extinction:
    def __init__(self, RA_lo, RA_hi, DEC_lo, DEC_hi):
        self.Name        = "Schlegel (1998) with the query API from Green et al. (2015)"
        self.Ang_Res     = 0.1017	#deg (=6.1 arcmin)
        self.bands       = list( pm.showbands ) 
        self.RA_lo       = RA_lo
        self.RA_hi       = RA_hi
        self.DEC_lo      = DEC_lo
        self.DEC_hi      = DEC_hi
        self.queried     = False
        self.f           = self.Setup_dust_grid()

    def Sample_extinction(self, ra, dec):
        """
        Sample the 2D dust grid to obtain the EBV at
         the location of the transient
        This is then converted to an extinction
        """
        EBV = self.f(ra,dec)
        A = {}
        for color in self.bands:
            A[color] = float(EBV) * RV_BV[color] 
        return A   

    def Setup_dust_grid(self):
        """
        Here we set up the 2D dust grid.
        Our grid is regular (and rectangular), which allows us to use
         RectBivariateSpline, which is faster than interp2d
        To also have data just across the border, we add 1 RA and DEC
         and add grid points 0.5*Ang_Res across this border
        We add 1 extra RA-coord to have an angular resolution slightly
         smaller than angres
        We then query the Schlegel dust map with the Argonaut query
        Subsequently we convert the 1D grid of SFD into a 2D grid
        And finally we Interpolate this grid
        """
        Nr_RA = np.ceil((self.RA_hi - self.RA_lo) / self.Ang_Res) + 2
        Nr_DEC = np.ceil((self.DEC_hi - self.DEC_lo) / self.Ang_Res) + 2
        ra  = np.linspace(self.RA_lo  - 0.5 * self.Ang_Res, 
                          self.RA_hi  + 0.5 * self.Ang_Res, Nr_RA)
        dec = np.linspace(self.DEC_lo - 0.5 * self.Ang_Res, 
                          self.DEC_hi + 0.5 * self.Ang_Res, Nr_DEC)
        RA, DEC = np.meshgrid(ra,dec)
        RA = list(chain.from_iterable(RA))
        DEC = list(chain.from_iterable(DEC))

        DQ = Greendustquery(RA, DEC, coordsys='equ', mode='sfd')
        SFD = DQ['EBV_SFD']
        EBV = np.zeros( (len(ra), len(dec)) )
        ra_len = len(ra)
        for i in range( len(dec) ):
            EBV[:,i] = SFD[i*ra_len: (i+1)*ra_len]

        self.queried = True        
        return RectBivariateSpline(ra, dec, EBV)
 



class Schultheis_Extinction:
    def __init__(self, Xtr_dust, RA_lo, RA_hi, DEC_lo, DEC_hi):
        self.Name        = "Schultheis et al. (2014)"
        self.Ang_Res     = 0.1	#deg (=6 arcmin)
        self.bands       = list( pm.showbands )
        self.RA_lo       = RA_lo
        self.RA_hi       = RA_hi
        self.DEC_lo      = DEC_lo
        self.DEC_hi      = DEC_hi
        self.Bulge_file  = h5py.File('Dustmaps/dustbulge.hdf5','r')
        self.queried     = False
        self.f           = self.Setup_dust_grid()
        self.Xtr_dust    = Xtr_dust

    def Sample_extinction(self, ra, dec, D):
        """
        Sample the 3D dust grid to obtain the EJK at
         the location of the transient
        This is then converted to an extinction
        The result might be either '-100' or 'nan'. Don't worry
         this will be fixed later
        Input:
        ra: RA-coordinate in degrees
        dec: DEC-coordinate in degrees
        D: distance to transient in kpc
        """
        D = D 
        #lon, lat = 0,0
        lon, lat = astCoords.convertCoords( "J2000", "GALACTIC", ra, dec, 2000 )
        if lon > 180: lon -= 360	#Here lon runs from -180 to 180
        A = {}
        EJK = self.f( [lon, lat, D] )
        if EJK == emptyval or math.isnan(EJK):	#Out of intpol range
            A = self.Xtr_dust.Sample_extinction(ra, dec)
        else:
            for color in self.bands:
                A[color] = float(EJK) * RV_JK[color] 
        return A 

    def Setup_dust_grid(self):
        """
        Loads the file with the dust data
        This file is a preprocessed version of the file provided
         by Schultheis et al. (2014)
        We interpolate this grid using the RegularGrindInterpolator
        This is done with the nearest neighbor method, because
         Schultheis et al. calculated extinction in bins
        If one samples a coordinate outside the interpolation 
         boundaries, it will return the value -100 
        """
        EJK = np.array(self.Bulge_file['EJK'][:])
        lon = np.array(self.Bulge_file['LON'][:])
        lat = np.array(self.Bulge_file['LAT'][:])
        D   = np.array(self.Bulge_file['DIST'][:]) /1000.	#to kpc
        self.queried = True
        return RegularGridInterpolator((lon, lat, D), EJK,
                                       method='nearest', 
                                       bounds_error = False, 
                                       fill_value = emptyval)




class No_dust:
    def __init__(self):
        self.Name        = "No dust"
        self.queried     = False
    def Sample_extinction(self, ra, dec, d=0):
        """
        We ignore dust, so, A=0
        """
        self.bands       = list( pm.showbands )
        A = {}
        for color in self.bands:
            A[color] = 0
        return A




def query_works ():
    """
    Test if we can connect to the Argonaut server by pinging it
    """
    try:
        requests.get('http://argonaut.skymaps.info/gal-lb-query-light')
    except requests.exceptions.ConnectionError:
        return False
    return True


def InGreenBoundary(RA_lo, RA_hi, DEC_lo, DEC_hi):
    """
    Test if the RA and DEC are within the boundaries
    This is tested by looking at the converged entry for this
     sightline. If any part of the observation window is
     within the boundary, we'll get green light to use this
     3D map.
    """
    RAs  = [RA_lo, RA_lo, RA_hi, RA_hi]
    DECs = [DEC_lo, DEC_hi, DEC_lo, DEC_hi]
    DQ   = dustquery(RAs, DECs, coordsys='equ', mode='full')
    return np.any(DQ['converged'])

def InSchultheisBoundary(RA_lo, RA_hi, DEC_lo, DEC_hi):
    """
    Test if any of the outer coordinates of the field of view
     are within the Schultheis+ dust map boundary
    If any of them is inside, return True, else return False
    """
    RAs = [RA_lo, RA_hi, RA_lo, RA_hi]
    DECs = [DEC_lo, DEC_lo, DEC_hi, DEC_hi]
    inSchultheis = False
    for i, RA in enumerate(RAs):
        lon, lat = astCoords.convertCoords("J2000", "GALACTIC", 
                                           RA, DECs[i], 2000)
        lonbool = 0. < lon < 10. or 350. < lon < 360
        latbool = -10. < lat < 5.
        if lonbool and latbool:
            inSchultheis = True
            break
    return inSchultheis

def Greendustquery(lon, lat, coordsys='equ', mode='full'):
    """
    A wrapper to make sure that the data is in the right format for
     the Argonaut server.
    The Argonaut server only takes arrays of length <5000. 
    This means that the coordinate arrays have to be cut up if they're
     larger than 5000.
    It also only extracts the important entries of the dustquery.
     Likewise, 'distmod' is only needed once in Dust_Ext
    """
    Total_len = len(lon)
    if Total_len > 5000:
        if mode == 'lite' or mode == 'full':
            Dust_Ext = {'best': [], 'distmod' : []}
        else:
            Dust_Ext = {'EBV_SFD': []}
        while len(lon) > 0:
            lon_part_len = min(5000, len(lon) )
            lon_part = lon[ 0:lon_part_len ]
            lon = lon[ lon_part_len: ]
            lat_part = lat[ 0:lon_part_len ]
            lat = lat[ lon_part_len: ]
            print ("Querying the Argonaut server for dust data in %d"
                   " out of %d data points" % (lon_part_len, Total_len))
            Dust = dustquery(lon_part, lat_part, coordsys, mode)
            if mode == 'lite' or mode == 'full':
                best = Dust_Ext['best']
                Dust_Ext['best'] += Dust['best']
                Dust_Ext['distmod'] = Dust['distmod']
            elif mode == 'sfd':
                Dust_Ext['EBV_SFD'] += Dust['EBV_SFD'] 
    else:
        print "Querying the Argonaut server for dust data..."
        Dust_Ext = dustquery(lon, lat, coordsys, mode)

    return Dust_Ext

def dustquery(lon, lat, coordsys='equ', mode='full'):
    '''
    Send a line-of-sight reddening query to the Argonaut web server.

    Inputs:
      lon, lat: longitude and latitude, in degrees.
      coordsys: 'gal' for Galactic, 'equ' for Equatorial (J2000).
      mode: 'full', 'lite' or 'sfd'
    
    In 'full' mode, outputs a dictionary containing, among other things:
      'distmod':    The distance moduli that define the distance bins.
      'best':       The best-fit (maximum proability density)
                    line-of-sight reddening, in units of SFD-equivalent
                    E(B-V), to each distance modulus in 'distmod.' See
                    Schlafly & Finkbeiner (2011) for a definition of the
                    reddening vector (use R_V = 3.1).
      'samples':    Samples of the line-of-sight reddening, drawn from
                    the probability density on reddening profiles.
      'success':    1 if the query succeeded, and 0 otherwise.
      'converged':  1 if the line-of-sight reddening fit converged, and
                    0 otherwise.
      'n_stars':    # of stars used to fit the line-of-sight reddening.
      'DM_reliable_min':  Minimum reliable distance modulus in pixel.
      'DM_reliable_max':  Maximum reliable distance modulus in pixel.
  
    Less information is returned in 'lite' mode, while in 'sfd' mode,
    the Schlegel, Finkbeiner & Davis (1998) E(B-V) is returned.
    '''
 
    url = 'http://argonaut.skymaps.info/gal-lb-query-light'
   
    payload = {'mode': mode}
    
    if coordsys.lower() in ['gal', 'g']:
        payload['l'] = lon
        payload['b'] = lat
    elif coordsys.lower() in ['equ', 'e']:
        payload['ra'] = lon
        payload['dec'] = lat
    else:
        raise ValueError("coordsys '{0}' not understood.".format(coordsys))
  
    headers = {'content-type': 'application/json'}
  
    r = requests.post(url, data=json.dumps(payload), headers=headers)
  
    try:
        r.raise_for_status()
    except requests.exceptions.HTTPError as e:
        print('Response received from Argonaut:')
        print(r.text)
        raise e
   
    return json.loads(r.text)


