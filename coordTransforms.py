import numpy as np
import MilkyWay as MW
from astLib import astCoords



def get_galactic_height( alpha, delta, dist_kpc ):
    l, b = astCoords.convertCoords( "J2000", "GALACTIC", alpha, delta, 2000 )
    r, z = galactic_to_cylindrical( b, l, dist_kpc )
    return z


def galactic_to_cylindrical( B, L, D ):
	b = np.deg2rad(B)
	l = np.deg2rad(L)

	Z = D*np.sin(b)
        rt1 = MW.MW_Rsun - D * np.cos(l) * np.cos(b)
        rt2 = -1.0 * D * np.sin(l)*np.cos(b)
	R = np.sqrt( rt1*rt1 + rt2*rt2 )
	return R,Z

def gal_to_cart( B, L, D ):

	b = np.deg2rad(B)
	l = np.deg2rad(L)

	
	X = MW.MW_Rsun - D * np.cos(l) * np.cos(b)
	Y = -1.0 * D * np.sin(l) * np.cos(b)
	Z = D * np.sin(b)
	return X,Y,Z
