import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import Bounds

from sgp4.earth_gravity import wgs84
from sgp4.io import twoline2rv
from sgp4.ext import days2mdhms, jday
from sgp4.propagation import sgp4init, sgp4
from sgp4.model import Satellite

import csv
from math import pi, pow

#def rosen(x):
    #"""The Rosenbrock function"""
    #return sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0)
#
#x0 = np.array([1.3, 0.7, 0.8, 1.9, 1.2])
#bounds = Bounds([0, 0, 0, 0, 0], [2.0, 2.0, 2.0, 2.0, 2.0])
#res = minimize(rosen, x0, method='nelder-mead', options={'xtol': 1e-8, 'disp': True},bounds=bounds)

def sgp4_compute(jdsatepoch,ecco,inclo,nodeo,argpo,mo,no,bstar,delta_min):
    deg2rad = pi/180.0;         
    xpdotp = 1440.0/(2.0*pi);

    satrec = Satellite();
    satrec.error = 0;
    whichconst = wgs84;
    satrec.whichconst = whichconst;
    tumin = whichconst.tumin
    satrec.satnum = 8195;
    satrec.jdsatepoch = jdsatepoch
    satrec.no = no;
    satrec.ecco = ecco;
    satrec.inclo = inclo;
    satrec.nodeo = nodeo;
    satrec.argpo = argpo;
    satrec.mo = mo;
    satrec.nddot = 0.00000e0;
    satrec.bstar = bstar;
    satrec.ndot = 0.000000;

    satrec.a = pow(satrec.no*tumin, (-2.0/3.0));
    satrec.alta = satrec.a*(1.0 + satrec.ecco) - 1.0;
    satrec.altp = satrec.a*(1.0 - satrec.ecco) - 1.0;
    
    satrec.elnum = 813;
    satrec.revnum = 22565;
    satrec.classification = 'U';
    satrec.intldesg = '75081A';
    satrec.ephtype = 0;

    # convert units and initialize
    satrec.no = satrec.no / xpdotp; 
    satrec.ndot  = satrec.ndot  / (xpdotp*1440.0);  
    satrec.nddot = satrec.nddot / (xpdotp*1440.0*1440);
    satrec.inclo = satrec.inclo  * deg2rad;
    satrec.nodeo = satrec.nodeo  * deg2rad;
    satrec.argpo = satrec.argpo  * deg2rad;
    satrec.mo    = satrec.mo     * deg2rad;

    afspc_mode = False

    sgp4init(whichconst, afspc_mode, satrec.satnum, satrec.jdsatepoch-2433281.5, satrec.bstar, satrec.ecco, satrec.argpo, satrec.inclo, satrec.mo, satrec.no, satrec.nodeo, satrec)
    return sgp4(satrec,delta_min)

def read_csv_data(filepath):
    t_gps = []
    x_gps = []
    v_gps = []

    # read in csv file
    with open(filepath) as csvfile:
        csv_rows = csv.reader(csvfile, delimiter=',')
        for row in csv_rows:
            t_gps.append(float(row[0]))
            x_gps.append([float(row[1]),float(row[2]),float(row[3])])
            v_gps.append([float(row[4]),float(row[5]),float(row[6])])
    
    np_t = np.array(t_gps)
    np_x_gps = np.array(x_gps)
    np_v_gps = np.array(v_gps)    
    return np_t, np_x_gps, np_v_gps

def get_sgp4_propagation(np_t,jdsatepoch,ecco,inclo,nodeo,argpo,mo,no,bstar,delta_min):
    x_sgp4 = []
    v_sgp4 = []

    for t in np_t:
        delta_min = t/60
        [r,v] = sgp4_compute(jdsatepoch,ecco,inclo,nodeo,argpo,mo,no,bstar,delta_min)
        x_sgp4.append(r)
        v_sgp4.append(v)

    np_x_sgp4 = np.array(x_sgp4)
    np_v_sgp4 = np.array(v_sgp4)
    return np_t, np_x_sgp4, np_v_sgp4

year = 2016
[mon,day,hr,minute,sec] = days2mdhms(year,111.34375000)
jdsatepoch = jday(year,mon,day,hr,minute,sec)
ecco = 85.1260441997737e-006;
inclo = 97.5200151325645e+000;
nodeo = 185.845569746093e+000;
argpo = 106.542369211203e+000;
mo = 273.306829557493e+000;
no = 15.1105221216793e+000;
bstar = 1.06508772985060e-003;
delta_min = 0

bounds = Bounds([0, 96, 180, 0, 0, 15, 1e-6], [0.1, 99, 190, 360, 360, 16, 1e-2])

filepath = './/data//sat_gps.csv'

[np_t, np_x_gps, np_v_gps] = read_csv_data(filepath)
[np_t, np_x_sgp4, np_v_sgp4] = get_sgp4_propagation(np_t,jdsatepoch,ecco,inclo,nodeo,argpo,mo,no,bstar,delta_min)

np_error = np.absolute(np_x_sgp4 - np_x_gps)

plt.title('GPS vs. SGP4 error')
plt.plot(np_t/84600,np_error)
#plt.plot(np_t/84600,np_x_sgp4,'--')
plt.grid()
plt.show()
