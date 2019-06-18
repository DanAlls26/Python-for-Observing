#Importing packages
import numpy as np 
from astropy import wcs
from  astropy.io import fits
import matplotlib.pyplot as plt 
import glob 
import argparse
import sys
from scipy.stats import scoreatpercentile 
from astropy.utils.data import get_pkg_data_filename
from astroplan.plots import plot_finder_image
from astroplan.plots import plot_sky
from astroplan import FixedTarget
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.patches as patches 
# note to self debating usefulness of argparse, sys, and glob(may not be needed)
# scipy may also not be necessary 
'''
#################################################################################
--------------------------------------------------------------------------------
################################################################################
Generating a Finding chart for the Breyo Observatory:

    - Current version for the SBIG STL-11000M ccd ONLY!
    - required packages:
        - Numpy
        - Astropy
        - Matplotlib
        - scipy
        - astroplan  (requires c++ build tools)
        - astroquery (requires c++ build tools)
        - argparse
        - glob 
    - Functionality:
        - generates a finding chart 

    - To do:
        - offset
        - rotation
        - loop 

################################################################################
--------------------------------------------------------------------------------
################################################################################
'''

# angular size/initial conditions of the rectangles 

#main imaging ccd
imcordx = np.array([0.20,0.20,0.6,0.6,0.4])
imcordy = np.array([0.17,0.62,0.62,0.17,0.4])

#guiding ccd
gcordx = np.array([0.12,0.12,0.18,0.18,0.16])
gcordy = np.array([0.38, 0.42, 0.42, 0.38,0.4])

#field of view of targeted space map
fov = 0.4
degpix = fov/150

#tra = float(input( "Target's Right Ascension in degrees:"))
#tdec = float(input("Target's Declination in degrees:"))

# for debugging -------------------------(remove when done)--------
tra = 267.9920073160415
tdec = 37.55208262984203
# -------------------------------------------------------------------

center_coord = SkyCoord(ra = tra*u.deg, dec = tdec*u.deg)
center = FixedTarget(coord=center_coord, name="Center")


imcordx = imcordx/degpix - 0.5
imcordy = imcordy/degpix - 0.5
gimcordx = gcordx/degpix - 0.5
gimcordy = gcordy/degpix - 0.5

# Ploting 
plt.figure()
ax,hdu = plot_finder_image(center,fov_radius = fov*u.deg)

#box for imaging ccd:
line1 = patches.FancyArrow(imcordx[0],imcordy[0],imcordx[1]-imcordx[0],imcordy[1]-imcordy[0],color='blue')
line2 = patches.FancyArrow(imcordx[1],imcordy[1],imcordx[2]-imcordx[1],imcordy[2]-imcordy[1],color='blue')
line3 = patches.FancyArrow(imcordx[2],imcordy[2],imcordx[3]-imcordx[2],imcordy[3]-imcordy[2],color='blue')
line4 = patches.FancyArrow(imcordx[3],imcordy[3],imcordx[0]-imcordx[3],imcordy[0]-imcordy[3],color='blue')

ax.add_patch(line1)
ax.add_patch(line2)
ax.add_patch(line3)
ax.add_patch(line4)

#box for guiding ccd: 
gline1 = patches.FancyArrow(gimcordx[0],gimcordy[0],gimcordx[1]-gimcordx[0],gimcordy[1]-gimcordy[0],color='red')
gline2 = patches.FancyArrow(gimcordx[1],gimcordy[1],gimcordx[2]-gimcordx[1],gimcordy[2]-gimcordy[1],color='red')
gline3 = patches.FancyArrow(gimcordx[2],gimcordy[2],gimcordx[3]-gimcordx[2],gimcordy[3]-gimcordy[2],color='red')
gline4 = patches.FancyArrow(gimcordx[3],gimcordy[3],gimcordx[0]-gimcordx[3],gimcordy[0]-gimcordy[3],color='red')

ax.add_patch(gline1)
ax.add_patch(gline2)
ax.add_patch(gline3)
ax.add_patch(gline4)

plt.show()

