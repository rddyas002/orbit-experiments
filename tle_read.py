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
        state = sat.at(sat.epoch)
        print("Epoch ", count, ": ", sat.epoch.utc_datetime(), " r,v: ", state.position.km, state.velocity.km_per_s)
        line1 = fp.readline()
        line2 = fp.readline()
        count = count + 1

fp.close()

line1 = ('1 12345U 12345A   19020.37500000  .00000000  00000-0  16883-3 0     7')
line2 = ('2 12345  51.6000   8.6001 0003960  33.0698 346.9303 15.58732945    08')
psuedosat = EarthSatellite(line1, line2, name='pseudosat-1')
state = psuedosat.at(psuedosat.epoch)
print("Epoch ", count, ": ", psuedosat.epoch.utc_datetime(), " r,v: ", state.position.km, state.velocity.km_per_s)
