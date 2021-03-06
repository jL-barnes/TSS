import numpy as np
from astLib import astCoords
from utils import get_galactic_height, galactic_to_cylindrical, gal_to_cart


# The Milky Way model (Juric et al. 2008)

MW_fThick            = 0.12       # thick disk ratio
MW_fHalo             = 0.0051     # halo ratio
MW_qHalo             = 0.7        # halo elipsidity
MW_pHalo             = 2.77       # halo power

MW_Rsun              = 8.0        # solar cylindrical R (kpc) 
MW_Zsun              = 0.025      # solar cylindrical Z (kpc)
MW_Lthin             = 2.6        # scale lengths and heights (all in kpc)
MW_Lthick            = 3.6
MW_Hthick            = 0.9        
MW_Hthin             = 0.3
MW_Lbulge            = 0.5        # source: Nelemans et al. (2004)
MW_rho_bulge         = 1.64e10
MW_r_break           = 18.0
MW_alpha_in          = 2.1
MW_alpha_out         = 3.9
MW_endHalo           = 70.0       # end of halo
MW_Rbulge            = 1.0        # end of bulge in R (kpc) 


DMax_MW_kpc = 200        # cut off for transients in the Milky Way, in kpc

# stellar density in the solar neighborhood
rho_stellar_sun = 8.5e7
rho_sun_thin = rho_stellar_sun/( 1.0 + MW_fThick + MW_fHalo )

def get_MW_dens( raDec_Coord, Hthin = 0.3 ):
    """ Calculate the Milky Way density at coordinates raDEC_Coord

    Parameters
    ----------
    reDEC_Coord : list of 3 floats
        The coordinates of the cell for which to calculate the density. It is
         in the form of [RA, DEC, Distance]
    Hthin : float
        Scale height of the thin disk in kpc
  
    Returns
    -------
    dens_tot : 1D-numpy array of 4 floats
        The Milky Way density as divided into its three components:
         [The thin disk contribution, thick disk, bulge, halo contribution]
    """
    this_RA, this_DEC, this_Dist = raDec_Coord
    lval, bval = astCoords.convertCoords("J2000", "GALACTIC", 
                                         this_RA, this_DEC, 2000)
    R, Z = galactic_to_cylindrical( bval, lval, this_Dist )
    rho_thin  = rho_stellar_sun * MW_dens_thin( R, Z, Hthin )
    rho_thick = rho_stellar_sun * MW_fThick * MW_dens_thick( R, Z )
    rho_bulge = MW_dens_bulge( R, Z)
    rho_halo  = rho_stellar_sun * MW_fHalo * MW_dens_halo( R, Z )
    dens_tot = [rho_thin, rho_thick, rho_bulge, rho_halo]
    return np.array(dens_tot)

def MW_dens_thin( R, Z, Hthin ):
    """ Calculate the transient density that is part of the thin disk

    Parameters
    ----------
    R : float
        Cylindrical radius coordinate in kpc
    Z : float
        Cylindrical height coordinate in kpc
    Hthin : float
        Scale height of the thin disk in kpc

    Returns
    -------
    e^expTerm : float
        The transient density at coordinates (R,Z) as part of the thin disk
    """
    expTerm = (MW_Rsun - R) / MW_Lthin - abs(Z) / Hthin
    return np.exp(expTerm)

def MW_dens_thick( R, Z ):
    """ Calculate the transient density that is part of the thick disk

    Parameters
    ----------
    R : float
        Cylindrical radius coordinate in kpc
    Z : float
        Cylindrical height coordinate in kpc

    Returns
    -------
    e^expTerm : float
        The transient density at coordinates (R,Z) as part of the thick disk
    """
    expTerm = (MW_Rsun - R) / MW_Lthick - abs(Z) / MW_Hthick
    return np.exp(expTerm)

def MW_dens_halo( R, Z ):
    """ Calculate the transient density that is part of the halo

    The halo is modeled with a broken power law with a break at MW_r_break
    The flex-point of the power law is at "cond"

    Parameters
    ----------
    R : float
        Cylindrical radius coordinate in kpc
    Z : float
        Cylindrical height coordinate in kpc

    Returns
    -------
    dens_halo : float
        The transient density at coordinates (R,Z) as part of the halo
    """
    cond = np.sqrt(R*R + (Z/MW_qHalo)*(Z/MW_qHalo))
    if cond > MW_Rbulge and cond < MW_r_break:
        t1 = MW_Rsun/np.sqrt( R*R + Z*Z/MW_qHalo/MW_qHalo )
        dens_halo = np.power( t1, MW_alpha_in )
    elif cond > MW_r_break and cond < MW_endHalo:
        t1 = np.power( MW_Rsun, MW_alpha_in )
        t2 = np.power( MW_r_break, MW_alpha_out - MW_alpha_in )
        t3 = 1.0/np.sqrt( R*R + Z*Z/MW_qHalo/MW_qHalo )
        t3 = np.power( t3, MW_alpha_out ) 
        dens_halo = t1 * t2 * t3
    else:
        dens_halo = 1.0e-99
    return dens_halo


def MW_dens_bulge( R, Z ):
    """ Calculate the transient density that is part of the bulge

    For this we do the same as Nelemans et al. (2004)

    Parameters
    ----------
    R : float
        Cylindrical radius coordinate in kpc
    Z : float
        Cylindrical height coordinate in kpc

    Returns
    -------
    e^expTerm : float
        The transient density at coordinates (R,Z) as part of the bulge
    """
    expTerm = -1.0*( R*R + Z*Z )/( MW_Lbulge*MW_Lbulge )
    return MW_rho_bulge * np.exp( expTerm )


