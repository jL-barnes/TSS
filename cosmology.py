from scipy.integrate import quad
from scipy import interpolate
from scipy.interpolate import RectBivariateSpline
import h5py
import numpy as np
import localcosmolopy as lcm

Mpc = 3.0857e22		#meter
H0 = 67.3

class Cosmology:
    """
    This class defines the cosmology
    Because we use a flat Lambda-CDM cosmology, Omega_k=0
     and we have no curvature
    We make this a class instead of loose functions in order
     speed up the functions: with loose functions one would have
     to calculate a redshift-function every iteration.
    """
    def __init__(self, MaxDist):
        self.s2yr = 1/(60*60.*24*365.25)
        self.H0 = H0 * 1e3 / Mpc	#Planck2015
        self.Omega_M = 0.315
        self.Omega_L = 0.685
        self.Omega_K = 0
        self.MaxDist = MaxDist		#Maximum Lum. Dist. of grid in kpc
        self.Zarray = np.logspace(-5,2,1e3)	#Range chosen to have all cosmological possibilities
        self.Darray = lcm.luminosity_distance(self.Zarray, 
                                              omega_M_0 = self.Omega_M, 
                                              omega_lambda_0 = self.Omega_L, 
                                              omega_k_0 = self.Omega_K, 
                                              h= self.H0/1e5 * Mpc)
        self.fn = interpolate.interp1d(self.Darray, self.Zarray)#,
                                       #fill_value = "extrapolate")
        self.factor = 0.001759927	#The normalization factor for PSI
                                    #Default is the factor for SN Ia

    def get_redshift(self, DL):
        """
        Computes the redshift at luminosity distance DL
        DL: Luminosity distance in Mpc
        returns: redshift
        """
        #print "maxdist", self.MaxDist
        #if self.fn(DL) > 0.1:
            #print DL, self.fn(DL)
        return self.fn(DL)


    def redshift_from_age (self, t):
        """
        With a flat Lambda-CDM cosmology the function of age to redshift
         can be simplified with a single analytial function.
        t: the age of the universe at that point
        returns: z: redshift at universe age t
        """
        return (np.sqrt( self.Omega_M / (1 - self.Omega_M) ) *
                np.sinh( 1.5 * self.H0 / self.s2yr * t *
                        np.sqrt(1 - self.Omega_M) ) ) **(-2./3.) -1

    def age_from_redshift (self, z):
        """
        With a flat Lambda-CDM cosmology the function of redshift to age
         can be simplified with a single analytial function.
        z: a redshift
        returns: t: age of universe at redshift z
        """     
        return ( 2 / (3 * self.H0 / self.s2yr * np.sqrt(1 -  self.Omega_M) ) *
               np.arcsinh( np.sqrt( (1 - self.Omega_M) / self.Omega_M) *
                          (1 + z)**(-3./2.) )  )

    def SFH (self, T):
        """
        Returns the Star Forming History at universe age T
        returns in units of M_odot yr^-1 Mpc-3 
        """
        return (0.01 * ( ( 1 + self.redshift_from_age(T) ) **2.6) / 
               (1 + ( (1 + self.redshift_from_age(T)) / 3.2 ) **6.2) )

    def PSI(self, t):
        """
        Returns the Delay Time Distribution (DTD) at time t
        It is normalized with Normalize_PSI
        """
        return ( t**self.alpha ) * self.factor 

    def Normalize_PSI(self):
        """
        Normalizes PSI such that the integral is equal to N/M
        The first white dwarfs/stellar explosions happen after 40 Myr
        """
        times = np.logspace( np.log10(0.04e9), 
                             np.log10(self.tmax), 1000)
        psidt = 0
        for i,T in enumerate(times[:-2]):
            psidt += T**self.alpha * (times[i+1] - times[i]) 
        self.factor = self.NM / psidt

    def create_Nr_density_funct(self, NM, alpha):
        """
        Creates the number density function so that one doesn't
         have to run it each time one wants to sample the function
        """
        self.Nfn = self.set_Nr_density_funct(NM, alpha)

    def integrand(self, tau, t):
        return self.SFH(t - tau) * self.PSI(tau)

    def set_Nr_density_funct(self, NM, alpha): 
        """
        Sets the number density function for extragalctic transients
        The first white dwarfs are only formed after 40Myr
        t is an array with ages of the universe at certain redshifts
        Because we want to create an interpolation at the end, we
         need to take the range of universe ages a bit more
         spaciously, hence the /1.1
        Vol_rate is the total (comoving) volumetric transient rate
         It is in integration over the SFH * DTD (in number yr^-1 Mpc^-3)
        """
        self.NM = NM			#Nr. Transients per M_\odot
        self.alpha = alpha		#exp. index of DTD
        self.tmax = self.age_from_redshift (0) #Current age of universe
        zmin = self.get_redshift(self.MaxDist / 1.e3)
        print "Calculating cosmology up to z = ", zmin
        tmin = lcm.age(zmin, omega_M_0 = self.Omega_M, 
                             omega_lambda_0 = self.Omega_L, 
                             omega_k_0 = self.Omega_K, 
                             h= self.H0/1e5 * Mpc) * self.s2yr
        t = np.linspace(tmin / 1.1, self.tmax,1000)
        Z = self.redshift_from_age(t) 
        Z[-1] = 0.0	#Due to rounding errors Python may calculate 
                        #the last Z as ~ -1e-16, which creates a problem

        if self.alpha == 0.0:	#No DTD
            Vol_rate = self.NM * self.SFH(t)
        else:
            self.Normalize_PSI()
            Integral = np.zeros(len(t))
            for i in range(len(t) - 1):
                Integral[i] = quad(self.integrand, 0.04e9, t[i], args=(t[i]))[0]
            Vol_rate = Integral

        return interpolate.interp1d(Z, Vol_rate)


    def Nr_density_xtr(self, z):
        """
        Returns the comoving number density of an extragalactic transient
         per year for redshift z
        """ 
        return self.Nfn(z) * 1e-9	#Convert to kpc^-3

class Kcor:
    def __init__(self, Kcormodel, colorsys, obBands):
        """
        The setup of the K-corrections function. This includes time dilation
        We take the file with the K-corrections (plus time dilation) and
         make an interpolation for every color band.
        The K-correction files have to be constructed with the jupyter notebook
         called Auxiliary/Kcor.ipynb. 
        Kcormodel: the name of the transient model as used in SNCosmo
        colorsystem: The color system in use: UBVRI, sdss, blackgem or lsst
        """
        Kcorfile = h5py.File('LightCurveFiles/Kcorrections/%s_%s.hdf5' % (Kcormodel, colorsys),'r')
        self.Kcorfuns = {}
        bands = obBands
        redshifts = Kcorfile['Z']
        times = Kcorfile['times']
        for band in bands:
            K_band = Kcorfile[band][:]
            self.Kcorfuns[band] = RectBivariateSpline(times, redshifts, K_band) 

    def Sample_Kcor(self, band, time, redshift):
        """
        Sample the K-correction interpolation functions that were constructed
         in the function above.
        Kcorfuns: the K-corrections function as defined in Setup_Kcor
        band: Color band in which to interpolate
        time: The time in days after explosion that has to be sampled
        redshift: redshift of the transient
        The redshift cannot be larger than 2 due to limitations of the K-correction data
        """
        if redshift > 2: redshift = 2
        return self.Kcorfuns[band](time, redshift)
    


