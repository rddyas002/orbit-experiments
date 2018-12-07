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
from math import pi, pow, sqrt

import time

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

def get_sgp4_propagation(np_t,jdsatepoch,ecco,inclo,nodeo,argpo,mo,no,bstar):
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
    
    x_sgp4 = []
    v_sgp4 = []

    for t in np_t:
        delta_min = t/60
        #[r,v] = sgp4_compute(jdsatepoch,ecco,inclo,nodeo,argpo,mo,no,bstar,delta_min)
        [r,v] = sgp4(satrec,delta_min)
        x_sgp4.append(r)
        v_sgp4.append(v)

    np_x_sgp4 = np.array(x_sgp4)
    np_v_sgp4 = np.array(v_sgp4)
    return np_t, np_x_sgp4, np_v_sgp4

def sgp4_error(x,jdsatepoch,np_t,r_gps,v_gps):
    e_sq = 0
    weighting = np.arange(len(np_t),0,-1)
    [np_t, np_x_sgp4, np_v_sgp4] = get_sgp4_propagation(np_t,jdsatepoch,x[0],x[1],x[2],x[3],x[4],x[5],x[6])

    pos_error = np.linalg.norm(np.transpose(r_gps - np_x_sgp4), axis=0)
    vel_error = np.linalg.norm(np.transpose(v_gps - np_v_sgp4), axis=0)
    pos_err_weighted = np.multiply(pos_error,weighting)
    vel_err_weighted = np.multiply(pos_error,weighting)
    pos_err_sq = np.square(pos_err_weighted)
    vel_err_sq = np.square(vel_err_weighted)
    err = sqrt(np.sum(pos_err_sq)) + sqrt(np.sum(vel_err_sq))

    return err

def sgp4_error_print(x,jdsatepoch,np_t,r_gps,v_gps):
    e_sq = 0
    weighting = np.arange(len(np_t),0,-1)
    [np_t, np_x_sgp4, np_v_sgp4] = get_sgp4_propagation(np_t,jdsatepoch,x[0],x[1],x[2],x[3],x[4],x[5],x[6])

    pos_error = np.linalg.norm(np.transpose(r_gps - np_x_sgp4), axis=0)
    vel_error = np.linalg.norm(np.transpose(v_gps - np_v_sgp4), axis=0)
    plt.plot(np_t/84600,pos_error)
    plt.show()
    pos_err_weighted = np.multiply(pos_error,weighting)
    pos_err_sq = np.square(pos_err_weighted)
    err = sqrt(np.sum(pos_err_sq))

    return err
        

year = 2016
[mon,day,hr,minute,sec] = days2mdhms(year,111.34375000)
jdsatepoch = jday(year,mon,day,hr,minute,sec)
ecco = 0.001
inclo = 97.5
nodeo = 185.0
argpo = 106.0
mo = 273.0
no = 15.11
bstar = 1.06508772985060e-004

ecco = 8.373941923917806e-05
inclo = 97.51977154369243
nodeo = 185.84556875835602
argpo = 136.12512076314346
mo = 243.7244849621949
no = 15.110521668513265
bstar = 0.0010642643733820103

print('Initial parameters:')
print('ecc = {0}'.format(ecco))
print('inclo = {0}'.format(inclo))
print('nodeo = {0}'.format(nodeo))
print('argpo = {0}'.format(argpo))
print('mo = {0}'.format(mo))
print('no = {0}'.format(no))
print('bstar = {0}'.format(bstar))

x0 = np.array([ecco,inclo,nodeo,argpo,mo,no,bstar])

bounds = Bounds([0, 96, 180, 0, 0, 15, 1e-6], [0.1, 99, 190, 360, 360, 16, 1e-1])

filepath = './/data//sat_gps.csv'

[np_t, np_x_gps, np_v_gps] = read_csv_data(filepath)
#err = sgp4_error_print(x0,jdsatepoch,np_t,np_x_gps,np_v_gps)

start = time.time()
print('Running optimisation...')
res = minimize(sgp4_error, x0, args=(jdsatepoch,np_t,np_x_gps,np_v_gps), method='nelder-mead', options={'xtol': 1e-12, 'disp': True})
end = time.time()
print('Time taken for minimization: {0:2f} seconds'.format(end-start))
#res = minimize(sgp4_error, x0, args=(jdsatepoch,np_t,np_x_gps,np_v_gps), method='L-BFGS-B', options={'disp': True},bounds=bounds)
if 'res' in locals():
    print('minimize operation succeed: {0}'.format(res.success))
    print('ecco = {0}'.format(res.x[0]))
    print('inclo = {0}'.format(res.x[1]))
    print('nodeo = {0}'.format(res.x[2]))
    print('argpo = {0}'.format(res.x[3]))
    print('mo = {0}'.format(res.x[4]))
    print('no = {0}'.format(res.x[5]))
    print('bstar = {0}'.format(res.x[6]))
#[np_t, np_x_sgp4, np_v_sgp4] = get_sgp4_propagation(np_t,jdsatepoch,ecco,inclo,nodeo,argpo,mo,no,bstar)

#np_error = np.absolute(np_x_sgp4 - np_x_gps)

#plt.title('GPS vs. SGP4 error')
#plt.plot(np_t/84600,np_error)
#plt.plot(np_t/84600,np_x_sgp4,'--')
#plt.grid()
#plt.show()
