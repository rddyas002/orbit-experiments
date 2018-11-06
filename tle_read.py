from sgp4.earth_gravity import wgs84
from sgp4.io import twoline2rv
from datetime import datetime, timedelta
import julian

import matplotlib.pyplot as plt
import numpy as np
import csv
import math

from skyfield.api import EarthSatellite, Topos, load

import pymap3d as pm

ts = load.timescale()

filepath = 'SUMBANDILA.TLE'
count = 1
with open(filepath) as fp:  
    line1 = fp.readline()
    line2 = fp.readline()
    while line2:
        sat = EarthSatellite(line1, line2, name='Sumbandila-1')
        print("Epoch ", count, ": ", sat.epoch.utc_datetime())
        line1 = fp.readline()
        line2 = fp.readline()
        count = count + 1

fp.close()      
