from scipy import interpolate
import localcosmolopy as lcm
import numpy as np


class Cosmology:
    def __init__(self):
        self.H0 = 67.8 #km/s/Mpc	#Planck2015
        self.Omega_M = 0.308
        self.Omega_L = 0.692
        self.Omega_K = 0
        #self.cosm = {'omega_M_0' : self.Omega_M, 'omega_lambda_0' : self.Omega_L, 
        #             'omega_k_0': self.Omega_K, 'h': self.H0/100.}
        self.Zarray = np.logspace(-5,2,1e3)	#Range chosen to have all cosmological possibilities
        self.Darray = lcm.luminosity_distance(self.Zarray, 
                                              omega_M_0 = self.Omega_M, 
                                              omega_lambda_0 = self.Omega_L, 
                                              omega_k_0 = self.Omega_K, 
                                              h= self.H0/100.)
        self.fn = interpolate.interp1d(self.Darray, self.Zarray)

    def get_redshift(self, DL):
        """
        Computes the redshift at luminosity distance DL
        DL: Luminosity distance in Mpc
        returns: redshift
        """
        return self.fn(DL)



class NoCosmology:
    def __init__(self):
        self.nothing = 0

    def get_redshift_at_D( DL ):
        """
        For intragalactic transients: there is no redshift
        """
        f = interpolate.interp1d(D2, Z2)
        return 0
