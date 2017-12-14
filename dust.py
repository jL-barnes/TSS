import numpy as np
import testparams as pm
#https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/LCYHJG
#source of dustmaps-package fro Green et al.
#http://dustmaps.readthedocs.io/en/latest/_modules/dustmaps/bayestar.html#BayestarQuery


RV = {'U':4.334, 'B':3.626, 'V':2.742, 'R':2.169, 'I':1.505, 'J':0.764, 'H':0.449, 'K':0.302, 'u':4.239, 'g':3.303, 'r':2.285, 'i':1.698, 'z':1.263}
#3.1 at http://iopscience.iop.org/article/10.1088/0004-637X/737/2/103#apj398709t6

class Green_Extinction:
    def __init__(self):
        self.Name        = "Green et al. (2017)"
        self.Ang_Res     = 9.8902e-4	#rad (=3.4 arcmin)
        self.query_works = 
        self.bands       = list( pm.showbands )

    def dustquery(lon, lat, coordsys='equ', mode='full'):
        '''
        Send a line-of-sight reddening query to the Argonaut web server.
    
        Inputs:
          lon, lat: longitude and latitude, in degrees.
          coordsys: 'gal' for Galactic, 'equ' for Equatorial (J2000).
          mode: 'full', 'lite' or 'sfd'
    
        In 'full' mode, outputs a dictionary containing, among other things:
          'distmod':    The distance moduli that define the distance bins.
          'best':       The best-fit (maximum proability density)
                        line-of-sight reddening, in units of SFD-equivalent
                        E(B-V), to each distance modulus in 'distmod.' See
                        Schlafly & Finkbeiner (2011) for a definition of the
                        reddening vector (use R_V = 3.1).
          'samples':    Samples of the line-of-sight reddening, drawn from
                        the probability density on reddening profiles.
          'success':    1 if the query succeeded, and 0 otherwise.
          'converged':  1 if the line-of-sight reddening fit converged, and
                        0 otherwise.
          'n_stars':    # of stars used to fit the line-of-sight reddening.
          'DM_reliable_min':  Minimum reliable distance modulus in pixel.
          'DM_reliable_max':  Maximum reliable distance modulus in pixel.
    
        Less information is returned in 'lite' mode, while in 'sfd' mode,
        the Schlegel, Finkbeiner & Davis (1998) E(B-V) is returned.
        '''
    
        url = 'http://argonaut.skymaps.info/gal-lb-query-light'
    
        payload = {'mode': mode}
    
        if coordsys.lower() in ['gal', 'g']:
            payload['l'] = lon
            payload['b'] = lat
        elif coordsys.lower() in ['equ', 'e']:
            payload['ra'] = lon
            payload['dec'] = lat
        else:
            raise ValueError("coordsys '{0}' not understood.".format(coordsys))
    
        headers = {'content-type': 'application/json'}
    
        r = requests.post(url, data=json.dumps(payload), headers=headers)
    
        try:
            r.raise_for_status()
        except requests.exceptions.HTTPError as e:
            print('Response received from Argonaut:')
            print(r.text)
            raise e
    
        return json.loads(r.text)

    def InBoundary(self):
        #Simply query, then search in dict for 'success'

    def Extinction(self, RA, DEC):
        #Still have to work with max/min-dist + ang res

        self.A = {}
        DQ = self.dustquery(RA, DEC, coordsys='equ', mode='lite')
        EBV = np.array( DQ['best'] )     ########Use the samples instead of best!!!!!!!
#Samples should give a better probabilistic representation of the data
        distmod = DQ['distmod']

        for color in self.bands:
            self.A[color] = EBV * RV[color] 

        return self.A   



class Schlegel_Extinction:
    def __init__(self, GreenExt):
        self.Name        = "Schlegel (1998) with the query API from Green et al. (2015)"
        self.Ang_Res     = 1.7744e-3	#rad (=6.1 arcmin)
        self.bands       = list( pm.showbands ) 

    def Extinction(self, RA, DEC):

        self.A = {}
        DQ = self.dustquery(RA, DEC, coordsys='equ', mode='lite')
        EBV = np.array( DQ['best'] )

        for color in self.bands:
            self.A[color] = EBV * RV[color] 

        return self.A   

class Schultheis_Extinction:
    def __init__(self):
        self.Name        = "Schultheis et al. (2014)"
        self.Ang_Res     = 1.7453e-3	#rad (=6 arcmin)
        self.bands       = list( pm.showbands )

    def InBoundary(self):




def Get_Extinction():

    successful = True

    GreenExt = Green_Extinction()
    SchlegelExt = Schlegel_Extinction()
    SchultheisExt = Schultheis_Extinction()

    if extraglactic():
         Extinction = SchlegelExt.extinction()
    elif GreenExt.InBoundary():
         Extinction = np.append( GreenExt.extinction(), SchlegelExt.extinction() ) #+ some interpolation using min./max. reliable distance
    elif SchultheisExt.InBoundary():
         Extinction = np.append( Schultheis.extinction(), SchlegelExt.extinction() ) #+ some interpolation         

#Make this a try, except
    if not Extinction.query_works:
         print "The extinction map from ", self.Name, " isn't working. "
         print "We'll exclude dust extinction in the simulation. "
         successful = False
         A = 0






