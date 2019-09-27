""" A few standard utility functions used throughout TSS
"""
import numpy as np
import MilkyWay as MW
from astLib import astCoords



def mag2flux(mag,f_0):
    """ Converts magnitude to flux in Jy
    """
    flux = f_0 * 10**(-0.4*mag)
    return flux

def flux2mag(flux,f_0):
    """ Converts flux (Jy) to magnitudes
    """
    mag = -2.5*np.log10(flux/f_0)
    return mag


def geometricP( nt ):
    """ Use Poisson statistics to obtain a number of transients inside
         a cell.

        nt: the average number of transients in this cell (float)

        -------
        returns: the actual number of transients in this cell (int)
    """
    return np.random.poisson(nt)

def getFileLine( file, tag ):
    """ Read a line from the dataTransients.dat file

        -------
        returns: The data belonging to transient type 'tag' 
    """
    f = open( file, 'r' )
    while f.readline() != '\n':
        pass
    data = f.readline()
    while data.split()[0] != tag:
        data = f.readline()
    restOfData = data.split()[1:]
    f.close()
    return restOfData

def get_galactic_height( alpha, delta, dist_kpc ):
    """ Get the galactic height that corresponds to equatorial
         coordinates (RA,DEC) = (alpha, delta) at distance dist_kpc

        -------
        returns: galactic height in kpc
    """
    l, b = astCoords.convertCoords( "J2000", "GALACTIC", alpha, delta, 2000 )
    r, z = galactic_to_cylindrical( b, l, dist_kpc )
    return z


def galactic_to_cylindrical( B, L, D ):
    """ Convert Galactic coordinates to cylindrical coordinates

        This is an approximation where that would be correct if Z_sun=0.
        Using Astropy's coordinate transforms takes 1000 times longer than this

        B and L: the galactic latitude and longitude
        D: Distance from the Sun

        -------
        returns: the cylindrical coordinates radius and height R and Z
    """
    b = np.deg2rad(B)
    l = np.deg2rad(L)

    Z = D * np.sin(b)
    rt1 = MW.MW_Rsun - D * np.cos(l) * np.cos(b)
    rt2 = -1.0 * D * np.sin(l)*np.cos(b)
    R = np.sqrt( rt1*rt1 + rt2*rt2 )
    return R,Z

def gal_to_cart( B, L, D ):
    """ Convert Galactic coordinates to carthesian coordinates with the
         center of the Milky Way as the origin of the carthesian
         coordinates.

        B and L: the galactic latitude and longitude
        D: Distance from the Sun

        -------
        returns: the carthesian coordinates X,Y,Z
    """
    b = np.deg2rad(B)
    l = np.deg2rad(L)

	
    X = MW.MW_Rsun - D * np.cos(l) * np.cos(b)
    Y = -1.0 * D * np.sin(l) * np.cos(b)
    Z = D * np.sin(b)
    return X,Y,Z
