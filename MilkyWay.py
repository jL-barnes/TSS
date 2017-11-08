import numpy as np
from astLib import astCoords

# Milky Way part keys

THIN  = 0
THICK = 1 
BULGE = 2
HALO  = 3


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
MW_rho_bulge         = 3.735e10
MW_r_break           = 18.0
MW_alpha_in          = 2.1
MW_alpha_out         = 3.9
MW_endHalo           = 70.0       # end of halo

DMax_MW_kpc = 200        # cut off for transients in the Milky Way, in kpc

# stellar density in the solar neighborhood
rho_stellar_sun = 8.5e7
rho_sun_thin = rho_stellar_sun/( 1.0 + MW_fThick + MW_fHalo )


def get_MW_dens( raDec_Coord, piecewise, Hthin = 0.3 ):
    this_RA, this_DEC, this_Dist = raDec_Coord
    lval, bval = astCoords.convertCoords( "J2000", "GALACTIC", this_RA, this_DEC, 2000 )
    R, Z = galactic_to_cylindrical( bval, lval, this_Dist )
    rho = MW_dens( R, Z, piecewise, Hthin )
    return rho
    
# set piecewise = 1 to get a list of densities: [ rho_thin, rho_thick, ... ]
# set piecewise == 0 to get a sum over all densities
def MW_dens( R, Z, piecewise, Hthin ):
    rho_thin  = MW_dens_thin( R, Z, Hthin )
    rho_thick = MW_fThick * MW_dens_thick( R, Z )
    rho_bulge = MW_dens_bulge( R, Z )
    rho_halo  = MW_fHalo * MW_dens_halo( R, Z )
    if piecewise == True:
        dens_tot = [rho_thin, rho_thick, rho_bulge, rho_halo] 
    else:
        dens_tot = rho_thin + rho_thick + rho_bulge + rho_halo
    return dens_tot

def MW_dens_thin( R, Z, Hthin ):
    expTerm = (MW_Rsun-R)/MW_Lthin - (abs(Z) + MW_Zsun)/Hthin
    return rho_sun_thin * np.exp( expTerm )

def MW_dens_thick( R, Z ):
    expTerm = (MW_Rsun-R)/MW_Lthick - (abs(Z) + MW_Zsun)/MW_Hthick
    return rho_sun_thin * np.exp( expTerm )

def MW_dens_bulge( R, Z ):
    expTerm = -1.0*( R*R + Z*Z )/( MW_Lbulge*MW_Lbulge )
    return MW_rho_bulge * np.exp( expTerm )

def MW_dens_halo( R, Z ):
    cond = np.sqrt(R*R + (Z/MW_qHalo)*(Z/MW_qHalo))
    if cond < MW_r_break:
        t1 = MW_Rsun/np.sqrt( R*R + Z*Z/MW_qHalo/MW_qHalo )
        dens_halo = rho_sun_thin * np.power( t1, MW_alpha_in )
    elif cond > MW_r_break and cond < MW_endHalo:
        t1 = np.power( MW_Rsun, MW_alpha_in )
        t2 = np.power( MW_r_break, MW_alpha_out - MW_alpha_in )
        t3 = 1.0/np.sqrt( R*R + Z*Z/MW_qHalo/MW_qHalo )
        t3 = np.power( t3, MW_alpha_out ) 
        dens_halo = rho_sun_thin * t1 * t2 * t3
    else:
        dens_halo = 1.0e-99
    return dens_halo

'''
Coordinate transformations
Sources:
http://rhodesmill.org/pyephem/coordinates.html
matrix http://arxiv.org/pdf/1010.3773.pdf
'''

def eq_to_gal( alpha, delta ):
	# Convert to radians
	a = np.deg2rad(alpha)
	d = np.deg2rad(delta)
        
        # Liu, Zhu, Zhang 2010 J2000
	a_G = 3.367
	d_G = .47304
        
        l0 = 5.2883

	atn = np.arctan( (np.sin(a_G-a))/(np.cos(a_G-a)*np.sin(d_G)-np.tan(d)*np.cos(d_G)))
	l = l0 - atn
        b = np.arcsin( np.sin(d)*np.sin(d_G) + np.cos(d)*np.cos(d_G)*np.cos(a_G-a) )
	while b > 2.0*np.pi:
		b -= 2*pi
	return b,l

def get_galactic_height( alpha, delta, dist_kpc ):
    b, l = eq_to_gal( alpha, delta )
    r, z = galactic_to_cylindrical( b, l, dist_kpc )
    return z

def galactic_to_cylindrical( B, L, D ):
	b = np.deg2rad( B )
	l = np.deg2rad( L )

	Z = D*np.sin(b)
        rt1 = MW_Rsun - D * np.cos(l) * np.cos(b)
        rt2 = -1.0 * D * np.sin(l)*np.cos(b)
	R = np.sqrt( rt1*rt1 + rt2*rt2 )
	return R,Z

def gal_to_cart( B, L, D ):

	b = np.deg2rad(B)
	l = np.deg2rad(L)
	
	X = MW_Rsun - D * np.cos(l) * np.cos(b)
	Y = -1.0 * D * np.sin(l) * np.cos(b)
	Z = D * np.sin(b)
	return X,Y,Z

