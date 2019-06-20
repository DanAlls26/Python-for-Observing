 #Importing packages
import numpy as np 
import matplotlib.pyplot as plt 
from astroplan.plots import plot_airmass
from astroplan.plots import plot_sky
from astroplan import FixedTarget
from astroplan import Observer
from astropy.coordinates import SkyCoord
from astropy.coordinates import EarthLocation
import astropy.units as u

'''
#################################################################################
--------------------------------------------------------------------------------
################################################################################
Generating an airmass plot for the Breyo Observatory: WIP!!!!!!

    - Current version for the Breyo Observatory in Loudonville, NY!
    - required packages:
        - Numpy
        - Astropy
        - Matplotlib
        - astroplan  (requires visual C++ build tools)
        - astroquery (requires visual C++ build tools)
    - Functionality:
        - generates an airmass chart
    - useage:
        - open an python environment terminal (ie anaconda prompt)
        - move to directory with this folder in it
        - type %run airMass.py
        - follow the directions appearing in the terminal 
    - To do:
        - make graph times EST instead of UTC 
        - impliment observing window
        - make sure DST is taken into account 
    - Future programs: 
        - full sky map 
        - utilize argparse (write help) 
        - overarching program that runs everything at once

################################################################################
--------------------------------------------------------------------------------
################################################################################
'''
print('Please note this code is very much a work in progress') 
# Defining location of observatory
loc = EarthLocation.from_geodetic(-73.751433*u.deg, 42.719546*u.deg, 106*u.m)
Breyo = Observer(location = loc, name = "Breyo", timezone ="US/Eastern")

# Defining time of observation
obsDate = input('Date of Observation run [yyyy-mm-dd]:')
obsTime = input('Time of Observation in military [hh:mm:ss]:')
obsDateTime = obsDate + ' ' +obsTime
print(obsDateTime)

# Target info
targetname = input('Name of the target:')
tra = float(input('Target RA in degrees:'))
tdec = float(input('Target DEC in degrees:'))

center_coord = SkyCoord(ra = tra*u.deg, dec = tdec*u.deg)
center = FixedTarget(coord=center_coord, name=targetname)

plot_airmass(center, Breyo, obsDateTime, brightness_shading=True, altitude_yaxis=True)
plt.show()
