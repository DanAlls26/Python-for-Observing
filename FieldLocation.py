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


'''
#################################################################################
--------------------------------------------------------------------------------
################################################################################
Generating a Finding chart for the Breyo Observatory:

    - Current version for the SBIG STL-11000M ccd ONLY!


################################################################################
--------------------------------------------------------------------------------
################################################################################
'''

# functions for use later in code

def cornerCords(image,length,height,binning):

    length = length/binning
    height = height/binning
    # Load the FITS hdulist using astropy.io.fits
    hdulist = fits.open(image)

    # Parse the WCS keywords in the primary HDU
    w = wcs.WCS(hdulist[0].header)

    # finding position of the corners of the image and center
    pixelList = np.array([[0,0],[length,0],[length,height],[0,height],
    [length/2,height/2]])
    world = w.wcs_pix2world(pixelList,1)

    return(world) 
 # 0 deg is north
 
#main finding chart generation:
mlength = 4008
mheight = 2745
binning = 2

#path = input (" Enter path:")
path = 'C:/Users/Daniel/Desktop/2019-06-03/'
print(path)
#mimage = input(" Enter the name of the target image(include .fit/.fits extension):")
mimage = 'fdzTrES-3bcal.fit'
print(mimage)
#gimage = input(" Enter the name of the guider image(include .fit/.fits extension):")
gimage = 'Autoguider-Image2.fit'
print(gimage)

gim = fits.getdata(path+gimage)
a = np.shape(gim)
glength = a[0]*2
gheight = a[1]*2


mainW = cornerCords(path+mimage,mlength,mheight,binning)
guideW = cornerCords(path+gimage,glength,gheight,binning)
center_coord = SkyCoord(ra = mainW[4,0]*u.deg, dec = mainW[4,1]*u.deg)
center = FixedTarget(coord=center_coord, name="Center")

gcenter_coord = SkyCoord(ra = guideW[4,0]*u.deg, dec = guideW[4,1]*u.deg)
gcenter = FixedTarget(coord = gcenter_coord, name="GuideCenter")


raW = mainW[:,0]
decW = mainW[:,1]
fov = 1
leftraW = (mainW[4,0]-fov)
leftdecW = (mainW[4,1]-fov)
degpix = fov/150
imcordx = (raW-leftraW)/degpix - 0.5
imcordy = (decW-leftdecW)/degpix - 0.5


graW = guideW[:,0]
gdecW = guideW[:,1]
gimcordx = (graW-leftraW)/degpix - 0.5
gimcordy = (gdecW-leftdecW)/degpix - 0.5
# print(gimcordx)

# print(graW-leftraW)
# print(gdecW-leftdecW)

print(mainW[4,0])
print(mainW[4,1])
'''
plt.figure()
ax,hdu = plot_finder_image(center,fov_radius = fov*u.deg)
point = patches.Circle( (imcordx[0],imcordy[0]), radius = 3,color='red')
point1 = patches.Circle( (imcordx[1],imcordy[1]), radius = 3,color='blue')
point2 = patches.Circle( (imcordx[2],imcordy[2]), radius = 3,color='black')
point3 = patches.Circle( (imcordx[3],imcordy[3]), radius = 3,color='green')

gpoint = patches.Circle( (gimcordx[0],gimcordy[0]), radius = 3)
gpoint1 = patches.Circle( (gimcordx[1],gimcordy[1]), radius = 3)
gpoint2 = patches.Circle( (gimcordx[2],gimcordy[2]), radius = 3)
gpoint3 = patches.Circle( (gimcordx[3],gimcordy[3]), radius = 3)

ax.add_patch(point)
ax.add_patch(point1)
ax.add_patch(point2)
ax.add_patch(point3)

ax.add_patch(gpoint)
ax.add_patch(gpoint1)
ax.add_patch(gpoint2)
ax.add_patch(gpoint3)

plt.show()


print('Finding Chart Completed')

'''
