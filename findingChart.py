# To Do List ----------------------------------------------------------------------------

'''
- Turn rotation and footprint masks into loops 

- Think of a better way to set up rotation 


'''


# Import Libraries ----------------------------------------------------------------------

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches

from astroplan.plots import plot_finder_image
from astroplan.plots import plot_sky
from astroplan import FixedTarget

from astropy.coordinates import SkyCoord
import astropy.units as u

import argparse

<<<<<<< HEAD
# Setting Up Command Line Arguments -----------------------------------------------------

parser = argparse.ArgumentParser()

# Positional Argument

parser.add_argument( "tpos", type = list,
    help =" Target coordinates of the finding image in degrees: [ra,dec]")

# Optional Arguments

parser.add_argument( "-fov", "--fieldOfView", type = int, default = 0.4,
    help = " Field of view in degrees. Defaults to 0.4")

parser.add_argument( "-off", "--offset", type = list, default = [0,0],
    help = " Vector distance from the center of the finding image to the target position in arcminutes. defaults to [ra=0,dec=0]")

parser.add_argument( "-r", "--rotate", type =float, default = 0,
    help = " Angle to rotate the image from north. Defaults to north at 0 degrees.")

parser.add_argument( "-tn","--targetName", type = str, default = 'Target',
    help = "Name of the target object")

parser.add_argument( "-v", "--verbose", action = 'store_true', default = False,
    help = " Program will talk while running")

# Description in help

args = parser.parse_args()


#description = '''Generates a finding chart for an target celestial object
#   for the Breyo Observatory's SBIG stl-11000m CCD.'''

# Unloading Arguments -------------------------------------------------------------------

target = args.tpos

fov = args.fieldOfView

off = args.offset

rotate = args.rotate

targetName = args.targetName

verbose = args.verbose

# Unpacking

tra = float(target[0])
tdec = float(target[1])

offra = float(off[0])
offdec = float(off[1])

# Defining CCD Shapes -------------------------------------------------------------------

#main ccd

if verbose == True:
    print("Defining the shape of the main ccd's sky footprint.")
    
=======
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
        - astroplan  (requires visual C++ build tools)
    - Functionality:
        - generates a finding chart
    - useage:
        - open an python environment terminal (ie anaconda prompt)
        - move to directory with this folder in it
        - type %run findingChart.py
        - follow the directions appearing in the terminal 
    - To do:
        - improve fact list
    - Future programs: 
        - full sky map 
        - airmass chart
        - utilize argparse (write help) 
        - overarching program that runs everything at once

################################################################################
--------------------------------------------------------------------------------
################################################################################
'''
print('---------------------------------------------')
targetname = input('Please input the name of your target:')
print('---------------------------------------------')
print('')
print('')
# angular size/initial conditions of the rectangles 

#main imaging ccd
>>>>>>> 39c9b0a7315465f4ea6f297b718d1d92138b0ec5
imcordx = np.array([0.25,0.25,0.55,0.55,0.4])
imcordy = np.array([0.17833335000000003,0.62166665,0.62166665,0.17833335000000003,0.4])

#guiding ccd

if verbose == True:
    print("Defining the shape of the guide ccd's sky footprint.")
   
gcordx = np.array([0.12,0.12,0.18,0.18,0.16])
gcordy = np.array([0.38, 0.42, 0.42, 0.38,0.4])

# degrees per pixel

if verbose == True:
    print("Determining image scale.")
    
degpix = fov/150

imcordx = ((imcordx - 0.4)/degpix)
imcordy = ((imcordy - 0.4)/degpix)
gimcordx = ((gcordx - 0.4)/degpix)
gimcordy = ((gcordy - 0.4)/degpix)

# Offset --------------------------------------------------------------------------------

if verbose == True:
    print('Offsetting image coordinates form target coordinates.')
    
    
offra = (offra*u.arcminute).to(u.deg)/u.deg
offdec = (offdec*u.arcminute).to(u.deg)/u.deg

cra = tra + offra
cdec = tdec +offdec

# Rotate Field --------------------------------------------------------------------------
if verbose == True and not rotate ==0 :
    print("Rotating orientation of ccd footprints on the sky.")
    
    
imangle = (np.arctan(imcordy[0:4]/imcordx[0:4])*u.rad).to(u.deg)
gimangle = (np.arctan(gimcordy[0:4]/gimcordx[0:4])*u.rad).to(u.deg)
imrad = imcordy[0:4]/np.sin(imangle.to(u.rad))
gimrad = gimcordy[0:4]/np.sin(gimangle.to(u.rad))

imangle = (imangle/u.deg + rotate)*u.deg
gimangle = (gimangle/u.deg + rotate)*u.deg
imangle = imangle.to(u.rad)
gimangle = gimangle.to(u.rad) 

imcordx = imrad*np.cos(imangle)
imcordy = imrad*np.sin(imangle)
gimcordx = gimrad*np.cos(gimangle)
gimcordy = gimrad*np.sin(gimangle)

# Pulling in Sky Coords -----------------------------------------------------------------

if verbose == True :
    print("Rotating orientation of ccd footprints on the sky.")
       
center_coord = SkyCoord(ra = cra*u.deg, dec = cdec*u.deg)
center = FixedTarget(coord=center_coord, name=targetName+' Finding Chart')

# Translation ---------------------------------------------------------------------------

imcordx = imcordx + 149.5
imcordy = imcordy + 149.5
gimcordx = gimcordx + 149.5
gimcordy = gimcordy + 149.5

# Plotting ------------------------------------------------------------------------------
if verbose == True :
    print("Pulling star field for finding chart from DSS.")


plt.figure()
ax, hdu = plot_finder_image(center,fov_radius = fov*u.deg)
    

if verbose == True :
    print("Creating footprint masks.")
#box for imaging ccd:

line1 = patches.FancyArrow(imcordx[0],imcordy[0],imcordx[1]-imcordx[0],imcordy[1]-imcordy[0],color='blue')
line2 = patches.FancyArrow(imcordx[1],imcordy[1],imcordx[2]-imcordx[1],imcordy[2]-imcordy[1],color='blue')
line3 = patches.FancyArrow(imcordx[2],imcordy[2],imcordx[3]-imcordx[2],imcordy[3]-imcordy[2],color='blue')
line4 = patches.FancyArrow(imcordx[3],imcordy[3],imcordx[0]-imcordx[3],imcordy[0]-imcordy[3],color='blue',label = 'Imaging CCD')

ax.add_patch(line1)
ax.add_patch(line2)
ax.add_patch(line3)
ax.add_patch(line4)

#box for guiding ccd:

gline1 = patches.FancyArrow(gimcordx[0],gimcordy[0],gimcordx[1]-gimcordx[0],gimcordy[1]-gimcordy[0],color='red')
gline2 = patches.FancyArrow(gimcordx[1],gimcordy[1],gimcordx[2]-gimcordx[1],gimcordy[2]-gimcordy[1],color='red')
gline3 = patches.FancyArrow(gimcordx[2],gimcordy[2],gimcordx[3]-gimcordx[2],gimcordy[3]-gimcordy[2],color='red')
gline4 = patches.FancyArrow(gimcordx[3],gimcordy[3],gimcordx[0]-gimcordx[3],gimcordy[0]-gimcordy[3],color='red',label='Guide CCD')

ax.add_patch(gline1)
ax.add_patch(gline2)
ax.add_patch(gline3)
ax.add_patch(gline4)

# circle target

circle = patches.Circle((149.5+(offra/degpix),149.5+(offdec/degpix)),
        radius=4,fill = False, color ='green',label=targetName)

ax.add_patch(circle)

# Legend

plt.legend(handles=[circle,line4,gline4])

# Generate plot

plt.show()

# finnished -----------------------------------------------------------------------------

