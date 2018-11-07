% get TLE data
tle_loc = 'REFERENCE.TLE';
fd = fopen(tle_loc, 'r');
longstr1 = fgets(fd, 130);
longstr2 = fgets(fd, 130);
fclose(fd);
KEP_TLE(1) = str2num(longstr1(19:32));      % epoch
KEP_TLE(2) = str2num(longstr2(9:16));       % inclination
KEP_TLE(3) = str2num(longstr2(18:25));      % RAAN
KEP_TLE(4) = str2num(longstr2(27:33))/1e7;  % Eccentricity
KEP_TLE(5) = str2num(longstr2(35:42));      % Argument of perigee
KEP_TLE(6) = str2num(longstr2(44:51));      % Mean anomaly
KEP_TLE(7) = str2num(longstr2(53:63));      % Mean motion
Bstar_exp  = str2num(longstr1(60:61));
KEP_TLE(8) = str2num(longstr1(54:59))*10^(Bstar_exp)/1e5;  
% Bstar
fd = fopen('AdcsTle.bin','w');
fwrite(fd,KEP_TLE,'double');
fclose(fd);