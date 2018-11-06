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
sim_stop = KEP_TLE(7)*Orbit_period*5;
sim_stop_max = max(GMAT_simdata.time);