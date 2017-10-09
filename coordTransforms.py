import numpy as np
from astLib import astCoords
from MilkyWay import MW_Rsun


'''
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
        b = np.arcsin( np.sin(d)*sin(d_G) + np.cos(d)*np.cos(d_G)*np.cos(a_G-a) )
	while b > 2.0*np.pi:
		b -= 2*pi
	return b,l


def galactic_to_cylindrical( B, L, D ):
	b = np.deg2rad(B)
	l = np.deg2rad(L)

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
