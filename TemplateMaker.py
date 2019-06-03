
#CODE TO PRODUCE A LIBRARY OF GALAXY TEMPLATES FOR USE LATER

# coding: utf-8

# In[1]:

from multiprocessing import Pool
import numpy as np
import os
import scipy
#from pylab import *
from astropy.io import fits
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u
import time
from scipy.interpolate import interp2d
from scipy.ndimage.interpolation import zoom
from scipy import optimize
from scipy import integrate
import tqdm


#os.system("taskset -p 0xff %d" % os.getpid())
# Define a series of functions to perfrom modelling of PSF

# In[2]:


# Define a series of functions used to help generate non-centered galaxies

# In[3]:


def gamma(s):
    """
    Define and return the value of the gamma function.
    """
    def integrand(t, s): return t**(s-1) * np.exp(-1*t)
    gi = scipy.integrate.quad(integrand, 0, np.inf, args=s)[0]
    return gi


def get_bn(n0):
    """
    Calculates the parameter bn from the Sersic profile. Original from GLACiAR code.

    Args:
        n0 (int) = Sersic index.

    Returns:
        bn (float) = Value of bn.
    """
    def errfunc(bn, n):
        return abs(scipy.special.gamma(2*n0) -
                   2*scipy.special.gammainc(2*n0, bn) *
                   scipy.special.gamma(2*n0))
    bn = scipy.optimize.fmin(errfunc, 1., args=(n0,), disp=False)[0]
    return bn

def makeSersic(n0, bn0, re, ell, inc_angle, size_galaxy, disy, disx):
    """
    Calculates the flux for each pixel following a Sersic profile.

    Args:
        n0 (int) = Sersic index.
        bn0 (float) = bn parameter.
        re (float) = Effective radius in pixels.
        ell (float) = Eccentricity. Varies between 0 and 1.
        inc_angle (float) = Inclination angle in radians. Varies
                            between 0 and Pi/2.
        size_galaxy (int) = Diameter of the galaxy stamp in pixels.
    Returns:
        fl (float) = Flux for each pixel.
	disx/y = displacement of galaxy from center in pixels (between -0.5 - +0.5)
    """
    stamp = np.zeros((size_galaxy,size_galaxy))
    s2 = size_galaxy / 2
    major_axis = re
    minor_axis = re * (1.-ell)
    I_e = ((bn0)**(2*n0)) / (2*np.pi*n0*major_axis*minor_axis*gamma(2*n0))
    def f(x, y):
        x_aux = (x-s2)*np.cos(inc_angle) + (y-s2)*np.sin(inc_angle)
        y_aux = -(x-s2)*np.sin(inc_angle) + (y-s2)*np.cos(inc_angle)
        radius = np.sqrt((x_aux/major_axis)**2 + (y_aux/minor_axis)**2)
        return I_e * np.exp(-bn0*((radius)**(1./n0)))

    for i in range(size_galaxy):
        def g(x):
            return i - 1./2. + disy
        def h(x):
            return i + 1./2. + disy
        for j in range(size_galaxy):
            fl = scipy.integrate.dblquad(f, j-(1./2.)+disx, j+(1./2.)+disx, g, h, 
                                 epsabs=1.49e-08, epsrel=1.49e-08)[0]
            stamp[i,j] = fl
    return stamp

def Templates(i):
	print("galaxy", i+1, "/", len(comb), comb[i])
	###GALAXY GENERATION
	t0=time.time()

	ind = comb[i][0] #Sersic Index n
	b = get_bn(ind) #Normalisation Based on n
	dx = comb[i][1]
	dy = comb[i][2]
	ell = comb[i][3]
	a = comb[i][4]
	Re= comb[i][5] * u.pc #size of galaxy half light
	Re= Re.to(u.Mpc)
	ReAng = ((Re* (1+Z)**2) / cosmo.luminosity_distance(Z)) * 57.2958 * 3600#size of galaxy in arcsec (converting from rad) via angular distance
	#print(cosmo.luminosity_distance(Z))
	galhalflight = ReAng/res
	#print("arcsec size =", ReAng)
	#print("pixel rad =", galhalflight)
	headers = '%s' % comb[i]
	stamps = makeSersic(ind, b, galhalflight, ell, a, stampsize, -dy, -dx)
	#CountConv = CountList[i]/stamps[1]
	Galaxy = stamps 
	np.savetxt('GalaxyTemplates/%s_%s.txt' % (name, i+1), Galaxy, header=headers)
	#tot+=1
	t1=time.time()
	t01 = time.time()
	print("time taken", t1-t0)
	tottime = t01-t00
	print(i, "galaxies in", tottime, "total time") 

# Begin MAIN component of the code

# In[4]:


"""
MAIN
"""
if __name__ == '__main__':
	t00 = time.time()
	### SET UP BASIC PARAMTERS AND PRIORS
	indlist = [1] #sersic index
	dispxlist = [-0.5, -0.25,  0, 0.25 ]
	dispylist = [-0.5, -0.25,  0,  0.25] #sub-pixel displacement values of galaxy center
	e= [0, 0.2, 0.4, 0.6, 0.8] #ellipticity
	ang = [0, (np.pi * 1./5.),  (np.pi * 1./5.), (np.pi * 2./5.), (np.pi * 3./5.), (np.pi * 4./5.)] #rotation angle from observer perspective
	GalSize = np.linspace(800, 12000, 5)#pc half light size

	Parameters = [indlist, dispxlist, dispylist, e, ang, GalSize]
	comb = [[]]
	for x in Parameters: #create list of all possible combinations of above parameters
		t=[]
		for y in x:
			for i in comb:
				t.append(i+[y])
		comb = t

	print(len(comb))
	total=len(comb) 
	name = 'series1' #reference name for other scripts to call upon to open these templates.

	np.savetxt('GalaxyTemplates/%s_params.txt' % name, comb)
	stampsize=30 #size of side of square stamp generated in pixels.


	Z = 4 #Redshift for these galaxies (Used for galaxy radii, recommend use center of range you want to consider as size can be rescaled later)

	res = 0.15 # resolution of image arcsec/pixel


	with Pool(processes=2) as pool: #number of cores
		list(tqdm.tqdm(pool.imap(Templates, range(0,total))))

