from scipy.integrate import quad
from scipy import interpolate
from scipy.interpolate import RectBivariateSpline
import h5py
import numpy as np
import localcosmolopy as cd

Mpc = 3.0857e22		#meter
H0 = 67.3

class Cosmology:
    """ This class defines the cosmology in TSS

    Because we use a flat Lambda-CDM cosmology, Omega_k=0 and we have
     no curvature
    We make this a class instead of loose functions in order to speed
     up the functions: with loose functions one would have to 
     calculate a redshift-function every iteration.

    Parameters
    ----------
    MaxDist : float
        The maximum luminosity distance of the grid in kpc

    Attributes
    ----------
    s2yr : float
        Conversion factor between seconds and years
    H0 : float
        The Hubble-Lemaitre constant
    Omega_M : float
        The Omega_M that belongs to our Lambda-CDM cosmology
    Omega_L : float
        The Omega_Lambda that belongs to our Lambda-CDM cosmology
    Omega_K : float
        The Omega_K that belongs to our Lambda-CDM cosmology
    Zarray : numpy 1D-array
        All possible redshifts
    Darray : numpy 1D-array
        The luminosity distances that correspond to the redshifts of
         Zarray
    fn : interp1d function
        An interpolatable function to convert a luminosity distance 
         into a redshift
    factor : float
        The Delay Time Distribution normalization factor
    """
    def __init__(self, MaxDist):
        self.s2yr = 1/(60*60.*24*365.25)
        self.H0 = H0 * 1e3 / Mpc	#Planck2015
        self.Omega_M = 0.315
        self.Omega_L = 0.685
        self.Omega_K = 0
        self.MaxDist = MaxDist		#Maximum Lum. Dist. of grid in kpc
        self.Zarray = np.logspace(-5,2,1e3)	#Range chosen to have all cosmological possibilities
        self.Darray = cd.luminosity_distance(self.Zarray, 
                                             omega_M_0 = self.Omega_M, 
                                             omega_lambda_0 = self.Omega_L, 
                                             omega_k_0 = self.Omega_K, 
                                             h= self.H0/1e5 * Mpc)
        self.fn = interpolate.interp1d(self.Darray, self.Zarray)
        self.factor = 0.001759927	#Default is the factor for SN Ia
        self.initialize_redshift_age()

    def get_redshift(self, DL):
        """ Computes the redshift at luminosity distance DL by 
         interpolating fn

        Parameters
        ----------
        DL : float
            Luminosity distance in Mpc

        Returns
        -------
        Redshift
        """

        return self.fn(DL)

    def initialize_redshift_age (self):
        """ Initialize a relation between redshift and age that can 
         later be interpolated

        Attributes (created)
        ----------
        z_age : interp1d function
            An interpolatable function in which you can put an array of
             universe ages (in years) and return the corresponding 
             redshift.
        """
        Times = cd.age(self.Zarray, use_flat=False,
                            omega_M_0 = self.Omega_M, 
                            omega_lambda_0 = self.Omega_L, 
                            omega_k_0 = self.Omega_K, 
                            h= self.H0/1e5 * Mpc) * self.s2yr #Current age of universe


        self.z_age = interpolate.interp1d(Times, self.Zarray, fill_value="extrapolate")

    def redshift_from_age (self, t):
        """ Compute the redshift that corresponds to an age of the 
         universe

        Parameters
        ----------
        t : numpy 1D-array or float
            The age of the universe at that point in years

        Returns
        -------
        z : numpy 1D-array or float
           Redshift at universe age t
        """
        return self.z_age(t)

    def SFH (self, T):
        """ Returns the Star Forming History at universe age T
         (amount of stars formed per year per volume)

        Parameters
        ----------
        T : numpy 1D-array or float
            The age of the universe at that point in years

        Returns
        ----------
        A numpy 1-D array or float with the Star Forming History in 
         units of M_odot yr^-1 Mpc-3
        """
        return (0.01 * ( ( 1 + self.redshift_from_age(T) ) **2.6) / 
               (1 + ( (1 + self.redshift_from_age(T)) / 3.2 ) **6.2) )

    def PSI(self, t):
        """ Returns the Delay Time Distribution (DTD) 

        It is normalized with Normalize_PSI

        Parameters
        ----------
        t : numpy 1D-array or float
            Time after stellar birth in years

        Returns
        -------
        A numpy 1D-array or float that is the DTD
        """
        return ( t**self.alpha ) * self.factor 

    def Normalize_PSI(self):
        """ Normalizes PSI such that the integral is equal to N/M

        The first white dwarfs/stellar explosions happen after 40 Myr
        """
        times = np.logspace( np.log10(0.04e9), 
                             np.log10(self.tmax), 1000)
        psidt = 0
        for i,T in enumerate(times[:-2]):
            psidt += T**self.alpha * (times[i+1] - times[i]) 
        self.factor = self.NM / psidt

    def create_Nr_density_funct(self, NM, alpha):
        """ Creates the number density function so that one doesn't
         have to run it each time one wants to sample the function

        Parameters
        ----------
        NM : float
            Number Transients per M_\odot
        alpha : float
            Exponential index of the Delay Time Distribution for this 
             type of transient.
        """
        self.Nfn = self.set_Nr_density_funct(NM, alpha)

    def integrand(self, tau, t):
        """ The function that needs to be integrated to get a 
         convolution of the star forming history with the delay-time
         distribution.

        Parameters
        ----------
        tau : numpy 1D-array or float
            The time (in years) after stellar birth
        t : numpy 1D-array or float
            The age of the universe at that point in years

        Returns
        -------
        An array or a float with the result of that convolution
        """
        return self.SFH(t - tau) * self.PSI(tau)

    def set_Nr_density_funct(self, NM, alpha): 
        """ Sets the number density function for extragalctic 
         transients.

        This needs to be done for every type of transients seperately
        The first white dwarfs are only formed after 40Myr
        t is an array with ages of the universe at certain redshifts
        Because we want to create an interpolation at the end, we
         need to take the range of universe ages a bit more
         spaciously, hence the *1.1 in the redshift range.
        Vol_rate is the total (comoving) volumetric transient rate
         It is in integration over the SFH * DTD (in number yr^-1 Mpc^-3)


        Attributes (created)
        ----------
        NM : float
            Number Transients per M_\odot
        alpha : float
            Exponential index of the Delay Time Distribution for this 
             type of transient.
        tmax : float
            Current age of universe

        Returns:
            An interp1d function which you can interpolate with a
             redshift to get a volumetric density rate        
        """
        self.NM = NM		
        self.alpha = alpha		#exp. index of DTD
        zmin = self.get_redshift(self.MaxDist / 1.e3)
        Z = np.linspace(1.1*zmin,0, 1000)
        t = cd.age(Z, use_flat=False,
                            omega_M_0 = self.Omega_M, 
                            omega_lambda_0 = self.Omega_L, 
                            omega_k_0 = self.Omega_K, 
                            h= self.H0/1e5 * Mpc) * self.s2yr 
        self.tmax = max(t)

        if self.alpha == 0.0:	#No DTD
            Vol_rate = self.NM * self.SFH(t)
        else:
            self.Normalize_PSI()
            Integral = np.zeros(len(t))
            for i in range(len(t) - 1):
                Integral[i] = quad(self.integrand, 0.04e9, t[i], args=(t[i]))[0]
            Vol_rate = Integral
        #The above will result in Vol_rate[z=0]=0. This is of course not true
        # So we set Vol_rate[z=0] to be the previous Vol_rate
        Vol_rate[-1] = Vol_rate[-2]
        return interpolate.interp1d(Z, Vol_rate)

    def Nr_density_xtr(self, z):
        """ Get the comoving number density of an extragalactic 
         transient

        Parameters
        ----------
        z : numpy 1D-array or float
         redshift
        
        Returns
        -------
        numpy 1D-array or float
         with the number density in per year per kpc^-3 for redshift z
        """ 
        return self.Nfn(z) * 1e-9	#Convert to kpc^-3


class Kcor:
    def __init__(self, Kcormodel, colorsys, obBands):
        """ The setup of the K-corrections function. This includes
         time dilation

        We take the file with the K-corrections (plus time dilation)
         and make an interpolation for every color band.
        The K-correction files have to be constructed with the jupyter 
         notebook called Auxiliary/Kcor.ipynb. 

        Parameters
        ----------
        Kcormodel : str
            The name of the transient model as used in SNCosmo
        colorsys : str
            The color system in use: UBVRI, sdss, blackgem, ztf or lsst
        obBands : numpy 1D-array
            The color bands in use

        Attributes
        ----------
        Kcorfuns : dict
            A dictionary with interpolatable Kcorrection functions for
             each colorbands. When you input a time plus redshift, you
             receive a magnitude
        """
        Kcorfile_up = h5py.File('LightCurveFiles/Kcorrections/%s_UBVRI.hdf5' % (Kcormodel),'r')
        Kcorfile_lo = h5py.File('LightCurveFiles/Kcorrections/%s_%s.hdf5' % (Kcormodel, colorsys),'r')
        self.Kcorfuns = {}
        for band in obBands:
            if band.islower():
                Kcorfile = Kcorfile_lo
            elif band.isupper():
                Kcorfile = Kcorfile_up
            redshifts = Kcorfile['Z']
            times = Kcorfile['times']
            K_band = Kcorfile[band][:]
            self.Kcorfuns[band] = RectBivariateSpline(times, redshifts, K_band) 

    def Sample_Kcor(self, band, time, redshift):
        """ Sample the K-correction interpolation functions that were 
         constructed in the function above.

        The redshift cannot be larger than 2 due to limitations of the K-correction data

        Parameters
        ----------
        band : str
            Name of the colorband of your observation
        time : float
            The time of the observation relative to the transient
        redshift : float
            The redshift of the transient of which you want a 
             K-correction

        Returns
        -------
        float
        An interpolation of Kcorfuns in band band for time time and 
         redshift redshift
        """
        if redshift > 2: redshift = 2
        return self.Kcorfuns[band](time, redshift)
    


