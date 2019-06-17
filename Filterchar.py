"""
Some telescope filter characteristics
"""


"""
RV_BV is the way a filter reacts to dust extinction. 
It is the conversion from E(B-V)_SFD to an extinction A_b
It's taken from table 6 (with R_V=3.1) from Schlafly & Finkbeiner (2011)
http://iopscience.iop.org/article/10.1088/0004-637X/737/2/103#apj398709t6
RV_JK is the conversion from E(J-K) to E(B-V)_SFD
"""

RV_BV = {'UBVRI':   {'U':4.334, 'B':3.626, 'V':2.742, 'R':2.169, 'I':1.505},
         'sdss':    {'u':4.239, 'g':3.303, 'r':2.285, 'i':1.698, 'z':1.263},
         'blackgem':{'u':4.091, 'g':3.277, 'r':2.286, 'i':1.672, 'z':1.230, 'q':2.625},
         'lsst':    {'u':4.145, 'g':3.237, 'r':2.273, 'i':1.684, 'z':1.323, 'y':1.088},
         'ztf':     {'g':3.265, 'r':2.203, 'i':1.590}
        } 
RV_JK = { key: { k: RV_BV[key][k] * 1.748 * 0.884 for k in RV_BV[key].keys() } 
          for key in RV_BV.keys()}

"""
The flux zero points of each filter in in erg cm^-2 s^-1 A^-1
UBVRI is in Vega magnitudes, the others in AB magnitudes, calculated with sncosmo
"""
flux_0 = {'UBVRI': {'U':417.5e-11, 'B':632.0e-11, 'V':363.1e-11, 'R':217.7e-11, 'I':112.6e-11},            
          'sdss':  {'u':847.e-11, 'g':490.e-11, 'r':287.e-11, 'i':195.e-11, 'z':137.e-11},
          'blackgem': {'u':754.e-11, 'g':472.e-11, 'r':279.e-11, 'i':186.e-11, 'z':130.e-11, 'q':324.e-11},
          'lsst': {'u':799.e-11, 'g':473.e-11, 'r':284.e-11, 'i':193.e-11, 'z':145.e-11, 'y':115.e-11},     
          'ztf':  {'g':470.e-11, 'r':264.e-11, 'i':175.e-11}
         }

"""
Conversion of Vega magnitudes to AB
Computed with the bandpasses in sncosmo
"""
Vega_to_AB = {'UBVRI': {'U': -0.81, 'B':0.10, 'V':-0.01, 'R':-0.19, 'I':-0.44}}

"""
Conversion of AB magnitudes to Vega
Computed with the bandpasses in sncosmo
"""

AB_to_Vega = {'sdss': {'u':0.89, 'g':-0.10, 'r':0.15, 'i':0.36, 'z':0.52},
              'blackgem': {'u':0.42, 'g':-0.10, 'r':0.16, 'i':0.38, 'z':0.52, 'q':0.05},
              'lsst': {'u':0.66, 'g':-0.09, 'r':0.15, 'i':0.37, 'z':0.51, 'y':0.55},
              'ztf': {'g':-0.10, 'r':0.18, 'i':0.42}
             }
