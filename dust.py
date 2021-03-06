import h5py
import math
import numpy as np
import Filterchar as fc
import astropy.units as u
from dustmaps import sfd
from dustmaps import bayestar
from astLib import astCoords
from itertools import chain
from astropy.coordinates import SkyCoord
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import RegularGridInterpolator



emptyval = 1000.0

class Green_Extinction:
    """ The Green extinction class handles dust extinction for the 
     Green et al. dust extinction.
    This is a 3D dust map for the Milky Way in the Northern Hemisphere

    Parameters
    ----------
    Xtr_dust : Schlegel_Extinction object
       An object of an extragalactic dust extinction class
    bands : numpy 1D-array
       An array of all color filters in use
    colorscheme : str
        The colorscheme in use as secondary (what ugriz color system)
    offline : bool
        Do computations offline?
    RA_lo : float
        Lowest RA coordinate of the frame
    RA_hi : float
        Highest RA coordinate of the frame
    DEC_lo : float
        Lowest DEC coordinate of the frame
    DEC_hi : float
        Highest DEC coordinate of the frame

    Attributes
    ----------
    Ang_res : float
        The angular resolution of this dust map in degrees
    queried : bool
        Whether the map has already been queried
    f : RegularGridInterpolator object
        This will contain a function which you can interpolate to get
         a dust extinction for a single coordinate
    RV_BV : dict
        Dictionary with the R_V(B-V) conversion factor
    """
    def __init__(self, Xtr_dust, bands, colorscheme, offline, RA_lo, RA_hi, DEC_lo, DEC_hi):
        self.Name        = "Green et al. (2019)"
        self.Ang_Res     = 0.05667	#deg (=3.4 arcmin) = minimum angular res.
        self.bands       = bands
        self.RA_lo       = RA_lo
        self.RA_hi       = RA_hi
        self.DEC_lo      = DEC_lo
        self.DEC_hi      = DEC_hi
        self.queried     = False
        self.offline     = offline
        self.f           = self.Setup_dust_grid()
        self.Xtr_dust    = Xtr_dust
        self.RV_BV       = fc.RV_BV[colorscheme]
        self.RV_BV.update(fc.RV_BV['UBVRI'])

    def Sample_extinction(self, ra, dec, D):
        """ Sample the 3D dust grid to obtain the EBV (dust extinction 
         coefficients) at the location of the transient.

        This is then converted to an extinction in each color band
        The result might be emptyval or 'nan'. This means that the
         result is out of the interpolation/query range. 
        In that case we switch to the extragalactic EBV

        Parameters
        ----------
        ra : float
            RA-coordinate in degrees
        dec : float
            DEC-coordinate in degrees
        D : float
            distance to transient in kpc

        Returns
        -------
        A : dict
            The dust extinction in units of magnitude
        """
        dmod = 5. * np.log10(D * 1000.) - 5.
        A = {}
        EBV = self.f( [ra, dec, dmod] )
        if EBV == emptyval or math.isnan(EBV):	#Out of intpol range
            A = self.Xtr_dust.Sample_extinction(ra, dec, D)
        else:
            for color in self.bands:
                A[color] = float(EBV) * self.RV_BV[color] 
        return A   

    def Setup_dust_grid(self):
        """ Here we set up the 3D dust grid. We use the 
         RegularGridInterpolator because the grid is regular and
         rectangular

        To also have data just across the border, we add 1 RA and DEC
         coordinate and add grid points 0.5*Ang_Res across this border
        We add 1 extra RA-coord to have an angular resolution slightly
         smaller than angres
        We then query the Green dust map
        Subsequently we convert the 2D object of DQ into a 3D grid
        We add data points for a distance modulus of zero (we assume 
         the EBV at d=0 to be zero)
        And finally we Interpolate this grid with linear interpolation
        If one samples a coordinate outside the interpolation 
         boundaries, it will return the value -100 

        Returns
        -------
        A function which you can interpolate to get a dust extinction
         for a single coordinate
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
        
        #Query Argonaut
        DQ = Greendustquery(RA,DEC, self.offline, Mode = 'galactic')	

        D = [0]		#Distance modulus zero isn't included in query
        D.extend(np.linspace(4.,19.,120))    #The distance moduli as queried
        Dlen = len(D)

        EBV = np.zeros((len(ra), len(dec), Dlen))
        for i in range( len(dec) ):
            for j in range( len(ra) ):
                Onecoord = [0]			#The EBV at d=0 at one coordinate
                Onecoord.extend(DQ[i+j])	#extend with EBV at other ds
                EBV[j,i,:] = np.array(Onecoord)

        self.queried = True
        
        return RegularGridInterpolator((ra,dec,D), EBV,
                                       method='linear',  
                                       bounds_error = False, 
                                       fill_value = emptyval)


class Schlegel_Extinction:
    """ The Schlegel extinction class handles dust extinction for the 
     Schlegel et al. dust extinction.
    This is an all-sky 2D dust map for the cumulative dust extinction.

    Parameters
    ----------
    bands : numpy 1D-array
       An array of all color filters in use
    colorscheme : str
        The colorscheme in use as secondary (what ugriz color system)
    offline : bool
        Do computations offline?
    RA_lo : float
        Lowest RA coordinate of the frame
    RA_hi : float
        Highest RA coordinate of the frame
    DEC_lo : float
        Lowest DEC coordinate of the frame
    DEC_hi : float
        Highest DEC coordinate of the frame

    Attributes
    ----------
    Ang_res : float
        The angular resolution of this dust map in degrees
    queried : bool
        Whether the map has already been queried
    f : RegularGridInterpolator object
        This will contain a function which you can interpolate to get
         a dust extinction for a single coordinate
    RV_BV : dict
        Dictionary with the R_V(B-V) conversion factor
    """
    def __init__(self, bands, colorscheme, offline, RA_lo, RA_hi, DEC_lo, DEC_hi):
        self.Name        = "Schlegel (1998)"
        self.Ang_Res     = 0.1017	#deg (=6.1 arcmin)
        self.bands       = bands
        self.RA_lo       = RA_lo
        self.RA_hi       = RA_hi
        self.DEC_lo      = DEC_lo
        self.DEC_hi      = DEC_hi
        self.queried     = False
        self.offline     = offline
        self.f           = self.Setup_dust_grid()
        self.RV_BV       = fc.RV_BV[colorscheme]
        self.RV_BV.update(fc.RV_BV['UBVRI'])

    def Sample_extinction(self, ra, dec, D):
        """ Sample the 2D dust grid to obtain the EBV at the location 
         of the transient

        Parameters
        ----------
        ra : float
            RA-coordinate in degrees
        dec : float
            DEC-coordinate in degrees
        D : float
            distance to transient in kpc

        Returns
        -------
        A : dict
            The dust extinction in units of magnitude
        """
        EBV = self.f(ra,dec)
        A = {}
        for color in self.bands:
            A[color] = float(EBV) * self.RV_BV[color] 
        return A   

    def Setup_dust_grid(self):
        """Here we set up the 2D dust grid.

        Our grid is regular (and rectangular), which allows us to use
         RectBivariateSpline, which is faster than interp2d
        To also have data just across the border, we add 1 RA and DEC
         and add grid points 0.5*Ang_Res across this border
        We add 1 extra RA-coord to have an angular resolution slightly
         smaller than angres
        We then query the Schlegel dust map with the Argonaut query
        Subsequently we convert the 1D grid of SFD into a 2D grid
        And finally we Interpolate this grid

        Returns
        -------
        A function which you can interpolate to get a dust extinction
         for a single coordinate
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

        DQ = Greendustquery(RA, DEC, self.offline, Mode='extragalactic')
        EBV = np.zeros( (len(ra), len(dec)) )
        ra_len = len(ra)
        for i in range( len(dec) ):
            EBV[:,i] = DQ[i*ra_len: (i+1)*ra_len]

        self.queried = True        
        return RectBivariateSpline(ra, dec, EBV)
 

class Schultheis_Extinction:
    """ The Schultheis extinction class handles dust extinction for the 
     Schultheis et al. dust extinction.
    This is a 3D dust map for the Milky Way bulge

    Parameters
    ----------
    Xtr_dust : Schlegel_Extinction object
       An object of an extragalactic dust extinction class
    bands : numpy 1D-array
       An array of all color filters in use
    colorscheme : str
        The colorscheme in use as secondary (what ugriz color system)
    offline : bool
        Do computations offline?
    RA_lo : float
        Lowest RA coordinate of the frame
    RA_hi : float
        Highest RA coordinate of the frame
    DEC_lo : float
        Lowest DEC coordinate of the frame
    DEC_hi : float
        Highest DEC coordinate of the frame

    Attributes
    ----------
    Ang_res : float
        The angular resolution of this dust map in degrees
    Bulge_file : h5py file object
        The file with the bulge dust data.
    queried : bool
        Whether the map has already been queried
    f : RegularGridInterpolator object
        This will contain a function which you can interpolate to get
         a dust extinction for a single coordinate
    RV_BV : dict
        Dictionary with the R_V(B-V) conversion factor
    useGreen : bool
        Whether we can use the Green et al. extinction if a coordinate
         is outside the Schultheis boundary
    """
    def __init__(self, Xtr_dust, bands, colorscheme, offline, RA_lo, RA_hi, DEC_lo, DEC_hi):
        self.Name        = "Schultheis et al. (2014)"
        self.Ang_Res     = 0.1	#deg (=6 arcmin)
        self.bands       = bands
        self.RA_lo       = RA_lo
        self.RA_hi       = RA_hi
        self.DEC_lo      = DEC_lo
        self.DEC_hi      = DEC_hi
        self.Bulge_file  = h5py.File('Dustmaps/dustbulge.hdf5','r')
        self.queried     = False
        self.f           = self.Setup_dust_grid()
        self.Xtr_dust    = Xtr_dust
        self.RV_JK       = fc.RV_JK[colorscheme]
        self.RV_JK.update(fc.RV_JK['UBVRI'])
        self.useGreen    = False
        if InGreenBoundary(RA_lo, RA_hi, DEC_lo, DEC_hi):
            self.useGreen = True
            self.Green = Green_Extinction(Xtr_dust, bands, colorscheme, offline, RA_lo, RA_hi, DEC_lo, DEC_hi)

    def Sample_extinction(self, ra, dec, D):
        """ Sample the 3D dust grid to obtain the EJK at
         the location of the transient

        This is then converted to an extinction
        The result might be either '-100' or 'nan'. In that case we sample
         the Green extinction map if possible

        Parameters
        ----------
        ra : float
            RA-coordinate in degrees
        dec : float
            DEC-coordinate in degrees
        D : float
            distance to transient in kpc

        Returns
        -------
        A : dict
            The dust extinction in units of magnitude
        """
        lon, lat = astCoords.convertCoords( "J2000", "GALACTIC", ra, dec, 2000 )
        if lon > 180: lon -= 360	#Here lon runs from -180 to 180
        A = {}
        EJK = self.f( [lon, lat, D] )
        if EJK == emptyval or math.isnan(EJK):	#Out of intpol range
            if self.useGreen:
                A = self.Green.Sample_extinction(ra, dec, D)
            else:
                A = self.Xtr_dust.Sample_extinction(ra, dec, D)
        else:
            for color in self.bands:
                A[color] = float(EJK) * self.RV_JK[color] 
        return A 

    def Setup_dust_grid(self):
        """Here we set up the 3D dust grid.


        Loads the file with the dust data
        This file is a preprocessed version of the file provided
         by Schultheis et al. (2014)
        We interpolate this grid using the RegularGrindInterpolator
        This is done with the nearest neighbor method, because
         Schultheis et al. calculated extinction in bins
        If one samples a coordinate outside the interpolation 
         boundaries, it will return the value -100 

        Returns
        -------
        A function which you can interpolate to get a dust extinction
         for a single coordinate
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
    """ A class that can serve as a replacement of the classes above
     when there is no dust extinction

    Parameters
    ----------
    bands : numpy 1D-array
       An array of all color filters in use


    Attributes
    ----------
    queried : bool
        Whether the map has already been queried
    """
    def __init__(self, bands):
        self.Name        = "No dust"
        self.queried     = False
        self.bands       = bands

    def Sample_extinction(self, ra, dec, d):
        """ We ignore dust, so, A=0

        Returns
        -------
        A : dict
            The dust extinction in units of magnitude
        """
        A = {}
        for color in self.bands:
            A[color] = 0
        return A

#%%
class Host_extinction:
    """ A class for host galaxy extinction.

    Parameters
    ----------
    In : str
        The type of host galaxy extinction. Can be (no, G%f, F%f) for 
         (no dust extinction, a Gaussian with sigma=%f, an exponential
         with sigma=\f)
    obRun : Observation instance
        The parent observation instance of which the grid instance is a
         child

    Attributes
    ----------
    bands : numpy 1D-array
       An array of all color filters in use
    RV_BV : dict
        Dictionary with the R_V(B-V) conversion factor
    """
    def __init__(self, In, obRun):
        self.no_host = False
        self.obRun   = obRun
        self.bands   = self.obRun.bands
        self.RV_BV   = fc.RV_BV[obRun.colorScheme]
        self.RV_BV.update(fc.RV_BV['UBVRI'])
        if In == 'no':
            self.no_host = True
        else:
            self.Exp_Gauss = In[0]
            self.sigma     = float(In[1:])
            if self.Exp_Gauss not in ['E', 'G', 'e', 'g']:
                print("Warning: Invalid distribution type entered.")
                print("Please enter E... or G... in the host extinction column in", self.obRun.transientFile)
                print("No host galaxy extinction is assumed")
                self.no_host = True

    def Sample_host_extinction(self):
        """ Sample the host galaxy extinction distribution.

        If no host galaxy extinction is allowed for this transient, 
         A=0 is returned.
        Otherwise the dust extinction is sampled from either:
            - An exponential decay
            - A one-sided gaussian centered at A=0
            Both with self.sigma as scale parameter
        We assume a R_V=3.1 extinction law similar to the Milky Way.

        Returns
        -------
        A : dict
            The dust extinction in units of magnitude
        """
        A = {}
        if self.obRun.nodust or self.no_host:
            for color in self.bands:
                A[color] = 0
            return A
        if self.Exp_Gauss in ['E', 'e']:
            AV = np.random.exponential(scale=self.sigma)
        elif self.Exp_Gauss in ['G', 'g']:
            AV = abs(np.random.normal(loc=0.0, scale=self.sigma))
        EBV = AV /  self.RV_BV['V']
        for color in self.bands:
            A[color] = float(EBV) * self.RV_BV[color] 
        return A
    

def InGreenBoundary(RA_lo, RA_hi, DEC_lo, DEC_hi):
    """ Test if the RA and DEC are within the boundaries

    This is tested by looking at the converged entry for this
     sightline. If any part of the observation window is
     within the boundary, we'll get green light to use this
     3D map.

    Parameters
    ----------
    RA_lo : float
        Lowest RA coordinate of the frame
    RA_hi : float
        Highest RA coordinate of the frame
    DEC_lo : float
        Lowest DEC coordinate of the frame
    DEC_hi : float
        Highest DEC coordinate of the frame

    Returns
    -------
    A boolean on whether any coordinate is within the boundaries
    """
    RAs  = np.array([RA_lo, RA_lo, RA_hi, RA_hi]) * u.deg
    DECs = np.array([DEC_lo, DEC_hi, DEC_lo, DEC_hi]) * u.deg
    Coords = SkyCoord(RAs, DECs, frame='icrs')
    Bayestar = bayestar.BayestarWebQuery(version='bayestar2019')
    DQ   = Bayestar(Coords, mode='best')
    return np.any(DQ > 0.)

def InSchultheisBoundary(RA_lo, RA_hi, DEC_lo, DEC_hi):
    """ Test if any of the outer coordinates of the field of view
     are within the Schultheis+ dust map boundary

    Parameters
    ----------
    RA_lo : float
        Lowest RA coordinate of the frame
    RA_hi : float
        Highest RA coordinate of the frame
    DEC_lo : float
        Lowest DEC coordinate of the frame
    DEC_hi : float
        Highest DEC coordinate of the frame

    Returns
    -------
    A boolean. If any of the coordinates is inside the boundary, return
     True, else return False
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

def Greendustquery(ra, dec, offline, Mode='galactic'):
    """ A wrapper to make sure that the data is in the right format for
     the Argonaut server.

    The Argonaut server only takes arrays of length <5000. 
    This means that the coordinate arrays have to be cut up if they're
     larger than 5000.
    It also only extracts the important entries of the dustquery.
     Likewise, 'distmod' is only needed once in Dust_Ext

    Parameters
    ----------
    ra : list
        List of RA coordinates of the points to be queried
    dec : list
        List of DEC coordinates of the points to be queried
    offline : bool
        Do computations offline?
    Mode : str, optional
        Either 'Galactic' or 'extragalactic'. If galactic, the 
         bayestar2019 (Green et al.) map will be queried, otherwise
         the SFD (Schlegel) map is queried.

    Returns
    -------
    Dust_Ext : numpy 1D-array
        An array of Dust extinction coefficients in units of E(B-V)_SFD
         for every (ra[i],dec[i]) coordinate.
    """
    Total_len = len(ra)
    if offline:
        if Mode == 'galactic':
            Query = bayestar.BayestarQuery(version='bayestar2019')
        elif Mode == 'extragalactic':
            Query = sfd.SFDQuery()    
    else:
        if Mode == 'galactic':
            Query = bayestar.BayestarWebQuery(version='bayestar2019')
        elif Mode == 'extragalactic':
            Query = sfd.SFDWebQuery()
    if Total_len > 5000:
        Dust_Ext = []
        while len(ra) > 0:
            ra_part_len = min(5000, len(ra) )
            ra_part = np.array(ra[ 0:ra_part_len ]) * u.deg
            ra = ra[ ra_part_len: ]
            dec_part = np.array(dec[ 0:ra_part_len ]) * u.deg
            dec = dec[ ra_part_len: ]
            print ("Querying remote server for Bayestar dust data in %d"
                   " out of %d data points" % (ra_part_len, Total_len))
            Coords = SkyCoord(ra_part, dec_part, frame='icrs')
            if Mode == 'galactic':
                Dust = Query(Coords, mode='best')
            elif Mode == 'extragalactic':
                Dust = Query(Coords)
            else: raise ValueError("Incorrect mode entered in the Dust query")
            Dust_Ext.extend(Dust)
        Dust_Ext = np.array(Dust_Ext)
    else:
        ra = np.array(ra) * u.deg
        dec = np.array(dec) * u.deg
        Coords = SkyCoord(ra, dec, frame='icrs')
        if Mode == 'galactic':
            print("Querying remote server for Galactic Bayestar dust data...")
            Dust_Ext = Query(Coords, mode='best')
        elif Mode == 'extragalactic':
            print("Querying remote server for extragalactic Bayestar dust data...")
            Dust_Ext = Query(Coords)

    return Dust_Ext
