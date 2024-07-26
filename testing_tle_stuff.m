% %% PREAMBEL
clear;
clc;
close all;


% % READ IN TLE
tle = read_in_tle("eive.tle");



% transformiere TLE Datenstruktur in sgp4 lesbare Datenstruktur
satdata = tle_to_satdata(tle);

% Abschätzung Algorithmus
% Berechne Orbitdauer
mean_orbit_duration = 2*pi/satdata.xno
if mean_orbit_duration<=225 %[min]
    disp("SGP4")
else
    disp("SDP4")
end
% ====================== CONFIGURE PROPAGATION ============================                                            
propagate_min = (24*60); % [minuten] Zeit, die man den Satellitenvektor propagieren möchte 

% SET FAVORITE LOCATION
latitude  = 48.7495 % [°]
longitude =  9.1038 % [°]


[rteme, vteme] = sgp4(propagate_min, satdata);


% % TRAFO


% 
% % SGP4
% [r_ecef,v_ecef] = teme2ecef(r_teme,v_teme,ttt,jdut1,lod,xp,yp,eqeterms)
% %  this function trsnforms a vector from the true equator mean equniox frame
% %    (teme), to an earth fixed (ITRF) frame.  the results take into account
% %    the effects of sidereal time, and polar motion.



% FROM CELESTIAL 

% from true equinox mean equator to Earth Centered Earth Fixed
% ECEF zu X


% conventional inertial i0 to mean inertial i_  = ADD PRECISION
% mean inertial i_ to instantaneous inertial T  = ADD NUTATION
% instantaneous inertial T to instantaneous terrestrial e = ADD GAST
% instantaneous terrestrial e to



% TO TOPOCENTRIC


% calculate checksum