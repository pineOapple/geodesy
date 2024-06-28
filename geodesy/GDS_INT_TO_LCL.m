function [R, G] = GDS_INT_TO_LCL(lambda, phi, yr, m, d, ut1)
% GAST, NUTATION, PRECESSION, JULICANCENTURIES
[t, ~, ~] = GDS_JULIANC(yr, m, d, ut1);
[gast, ~] = GDS_JUL_TO_GAST(ut1, t);
[eps0, deleps, delpsi] = GDS_NUTATION(t);
[za, thetaa, zetaa] = GDS_PRECESSION(t);

% Systemtransformation, wir drehen das System und arbeiten mit
% Objektrotationsmatrizen, daher überall inverse Vorzeichen.
N = rotx(eps0+deleps)*rotz(delpsi)*rotx(-eps0);
P = rotz(za)*roty(-thetaa)*rotz(zetaa);
R = roty(-90+phi)*rotz(-lambda)*rotz(-(360*gast)/23.93447192)*N*P;

% Anderes Vorzeichen, da wir unseren GAST im Inertialsystem und nicht im
% Lokalsystem sehen wollen!
G = rotz(-(360*gast)/23.93447192);