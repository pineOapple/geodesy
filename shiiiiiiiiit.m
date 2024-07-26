clc
close all
clear

%% COORDINATES AND SHIT
phi = -70.6734; % latitude in degrees
lambda = -8.2741; % longitude in degrees

%% DATE AND SHIT
yr = 2024;
mt = 9;
dy = 21;
ut1 = 12;
startTime = datetime(yr, mt, dy, ut1, 0, 0);

[R0, P0, N0, G0, R1, P1, N1, G1, t, mjd, jd, gast, gmst, eps0, deleps, delpsi, za, thetaa, zetaa] = GDS_INT_TO_LCL(lambda, phi, yr, mt, dy, ut1);

tle_path = 'data/eive.tle'
tleStruct = tleread(tle_path)

%% Getting satellite radius
stopTime = startTime + days(2);
sampleTime = 60;
sc = satelliteScenario(startTime,stopTime,sampleTime);
[r,v] = propagateOrbit(startTime,tleStruct);

EIVE = satellite(sc,tle_path, ...
    "Name","EIVE", ...
    "OrbitPropagator","two-body-keplerian");

% r_e = R0*r-[6378135 6378135 6378135]';



v = satelliteScenarioViewer(sc);
plot3(r(1),r(2),r(3));
satSGP4.MarkerColor = [0 1 0];
satSGP4.Orbit.LineColor = [0 1 0];
satSGP4.LabelFontColor = [0 1 0];
satSDP4.MarkerColor = [1 0 1];
satSDP4.Orbit.LineColor = [1 0 1];
satSDP4.LabelFontColor = [1 0 1];