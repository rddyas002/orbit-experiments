clear all;
% convert ASCII tle to bin
run('C:\Users\Yashren\work\orbit-experiments\lte2bin.m')

% import csv of GMAT orbit
GMAT = csvread('data\ReportFile1.txt');
% format:
% Secs,lat,long,alt,X,Y,Z,Vx,Vy,Vz
GMAT_simdata.time = GMAT(:,1);
GMAT_simdata.signals.values = GMAT(:,5:10);
orbit_start_UTC = [2019 01 20 09 00 00];

Orbit_period = 24*60*60/KEP_TLE(7);     % seconds
GPS_on_percent = 0.05;
GPS_time_sequence = [0 Orbit_period*GPS_on_percent Orbit_period*GPS_on_percent+1 Orbit_period/2-1 Orbit_period/2 Orbit_period/2+Orbit_period*GPS_on_percent Orbit_period/2+Orbit_period*GPS_on_percent+1 Orbit_period];
GPS_value_sequence = [1 1 0 0 1 1 0 0];
%GPS_time_sequence = [0 Orbit_period*GPS_on_percent Orbit_period*GPS_on_percent+1 Orbit_period];
%GPS_value_sequence = [1 1 0 0];
sim_stop = Orbit_period*5*KEP_TLE(7);
sim_stop_max = max(GMAT_simdata.time);
mean_motion = KEP_TLE(7)*2*pi/(24*60*60);
mu = 398600.4415;
semi_major = (mu/mean_motion^2)^(1/3);