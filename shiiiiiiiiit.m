clc
close all
clear

%% shiasdf
tle_path = 'data/eive.tle'
tleStruct = tleread(tle_path)
startTime = datetime(2024, 7, 25, 12, 0, 0);
stopTime = startTime + days(2);
sampleTime = 60;
sc = satelliteScenario(startTime,stopTime,sampleTime)
[r,v] = propagateOrbit(startTime,tleStruct);
EIVE = satellite(sc,tle_path, ...
    "Name","EIVE", ...
    "OrbitPropagator","two-body-keplerian");

v = satelliteScenarioViewer(sc);
satSGP4.MarkerColor = [0 1 0];
satSGP4.Orbit.LineColor = [0 1 0];
satSGP4.LabelFontColor = [0 1 0];
satSDP4.MarkerColor = [1 0 1];
satSDP4.Orbit.LineColor = [1 0 1];
satSDP4.LabelFontColor = [1 0 1];