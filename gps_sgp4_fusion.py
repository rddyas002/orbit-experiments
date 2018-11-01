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

line1 = ('1 12345U 12345A   16111.34375000 0.00000000  00000-0  10000-4 0     0')
line2 = ('2 12345  98.1900 185.8460 0000738  94.0712 285.9289 14.57010811    08')
sat = EarthSatellite(line1, line2, name='EOSat1')

dt = datetime(2016, 4, 20, 8, 15, 0)
t_simulation = 60000
delta_t = 1.0
sgp4_time = []
sgp4_time_datetime = []
position = []
velocity = []
geodetic_sgp4 = []
# initialise sgp4
satellite = twoline2rv(line1, line2, wgs84)
counter = 0
# run SGP$ first
while counter < t_simulation:
    # get gregorian time
    ps, vl = satellite.propagate(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
    position.append([ps[0],ps[1],ps[2]])
    velocity.append([vl[0],vl[1],vl[2]])
    ecef = pm.eci2ecef((ps[0]*1e3, ps[1]*1e3, ps[2]*1e3) ,dt)
    geodetic_sgp4.append(pm.ecef2geodetic(ecef[0],ecef[1],ecef[2]))
    sgp4_time.append(float(counter))
    sgp4_time_datetime.append(dt)
    #print('SGP4: ',geodetic_sgp4[counter])
    counter = counter + 1
    dt = dt + timedelta(seconds=delta_t)
    if (counter % 1000) == 0:
        print(counter, ' seconds running')

np_SGP4_position = np.array(position)
np_SGP4_velocity = np.array(velocity)

ElapsedSecs = []
GPS_position = []
GPS_velocity = []
GPS_altitude = []
X = []
Y = []
Z = []
Vx = []
Vy = []
Vz = []
counter = 0
filename = 'data\ReportFile1.txt'
with open(filename, 'r') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        if counter != 0:
            ElapsedSecs.append(float(row[0]))
            GPS_altitude.append([float(row[1]), float(row[2]), float(row[3])/1e3])
            GPS_position.append([float(row[4]), float(row[5]), float(row[6])])
            GPS_velocity.append([float(row[7]), float(row[8]), float(row[9])])

        counter = counter + 1
        if counter > t_simulation:
            break

np_GPS_position = np.array(GPS_position)
np_GPS_velocity = np.array(GPS_velocity)
np_Pos_error = (np_GPS_position - np_SGP4_position)*1e3
np_Vel_error = np_GPS_velocity - np_SGP4_velocity

tm = ts.utc(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
geo = sat.at(tm)

# plot sgp4 position
plt.figure(1,[5,5])
#plt.subplot(2,1,1)
plt.plot(np_Pos_error)
plt.show()

    
#plt.figure(1,[5,5])
#plt.plot(ElapsedSecs,Latitude,'r')
#plt.xlabel('Time (seconds)')
#plt.ylabel('LLA (degrees)')
#plt.show()
#plt.savefig(final_directory+'/'+filename + '_plot_' + test + '_.png', bbox_inches='tight')

#print(satellite.error)    # nonzero on error
#print(satellite.error_message)
#print(position)
#print(velocity)
