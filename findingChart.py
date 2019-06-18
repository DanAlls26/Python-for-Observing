#Importing packages
import numpy as np 
import matplotlib.pyplot as plt 
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
    - required packages:
        - Numpy
        - Astropy
        - Matplotlib
        - scipy
        - astroplan  (requires visual C++ build tools)
        - astroquery (requires visual C++ build tools)
        - argparse
        - glob 
    - Functionality:
        - generates a finding chart
    - useage:
        - open an python environment terminal (ie anaconda prompt)
        - move to directory with this folder in it
        - type %run findingChart.py
        - follow the directions appearing in the terminal 
    - To do:
        - improve fact list 

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
imcordx = np.array([0.2460,0.2460,0.5541,0.5541,0.4])
imcordy = np.array([0.17,0.62,0.62,0.17,0.4])

#guiding ccd
gcordx = np.array([0.12,0.12,0.18,0.18,0.16])
gcordy = np.array([0.38, 0.42, 0.42, 0.38,0.4])

#field of view of targeted space map
fov = 0.4
degpix = fov/150
print('---------------------------------------------')
print('Please put in the target Coordinates below.')
print('')
tra = float(input( "Target's Right Ascension in degrees:"))
tdec = float(input("Target's Declination in degrees:"))
print('---------------------------------------------')
print('')
print('')
# for debugging -------------------------(comment out when done)--------
#tra = 267.9920073160415
#tdec = 37.55208262984203
# -------------------------------------------------------------------

print('---------------------------------------------')
print('For the offset settings if you desire no offset at all input a 0')
print('Leaving these blank will BREAK the program!!!')
print('For moving the image to the left please make the offset RA positive.')
print('For moving the image to the right please make the offset RA negative.')
print('For moving the image to up please make the offset DEC positive.')
print('For moving the image to down please make the offset DEC negative.')
print('---------------------------------------------')
print('')
offra = float(input( "Right Ascension offset in arcminutes:"))
offdec = float(input("Declination offset in arcminutes:"))
print('---------------------------------------------')
print('')
print('')
#for debugging -------------------------(comment out when done)--------
#offra = 3
#offdec = 0.3
# -------------------------------------------------------------------

offra = (offra*u.arcminute).to(u.deg)/u.deg
offdec = (offdec*u.arcminute).to(u.deg)/u.deg

tra = tra + offra
tdec = tdec +offdec

#rotate field
print('---------------------------------------------')
print('Please input a rotation angle. If you do not want a rotation angle')
print('please enter a 0 below')
print('---------------------------------------------')
print('')
rotate = float(input('Rotation angle in degrees:'))
print('---------------------------------------------')
print('')
print('')

    
imcordx = ((imcordx - 0.4)/degpix)
imcordy = ((imcordy - 0.4)/degpix)
gimcordx = ((gcordx - 0.4)/degpix)
gimcordy = ((gcordy - 0.4)/degpix)

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

# setting up for ploting 
center_coord = SkyCoord(ra = tra*u.deg, dec = tdec*u.deg)
center = FixedTarget(coord=center_coord, name=targetname)

imcordx = imcordx + 149.5
imcordy = imcordy + 149.5
gimcordx = gimcordx + 149.5
gimcordy = gimcordy + 149.5
  
plt.figure()
ax,hdu = plot_finder_image(center,fov_radius = fov*u.deg)

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
circle = patches.Circle((149.5+(offra/degpix),149.5+(offdec/degpix)),radius=4,fill = False, color ='green',label='target')
ax.add_patch(circle)

# Legend
plt.legend(handles=[circle,line4,gline4])
plt.show()

# just one more thing
#fact list:
print('---------------------------------------------')
print('Fact List: (WIP)')
print('---------------------------------------------')
print('The center RA of this finding image is %(ra)f degrees' %{'ra': tra})
print('The center DEC of this finding image is %(dec)f degrees' %{'dec': tdec})
print('---------------------------------------------')

  




