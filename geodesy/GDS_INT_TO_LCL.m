function [R0,P0, N0, G0,R1,P1, N1, G1, t, mjd, jd, gast, gmst, eps0, deleps, delpsi, za, thetaa, zetaa] = GDS_INT_TO_LCL(lambda, phi, yr, m, d, ut1)
% GAST, NUTATION, PRECESSION, JULICANCENTURIES
[t, mjd, jd] = GDS_JULIANC(yr, m, d, ut1);
[gast, gmst] = GDS_JUL_TO_GAST(ut1, t);
[eps0, deleps, delpsi] = GDS_NUTATION(t);
[za, thetaa, zetaa] = GDS_PRECESSION(t);

% Systemtransformation, wir drehen das System und arbeiten mit
% Objektrotationsmatrizen, daher überall inverse Vorzeichen.
N0 = rotx(eps0+deleps)*rotz(delpsi)*rotx(-eps0);
P0 = rotz(za)*roty(-thetaa)*rotz(zetaa);
G0 = rotz(-(360*gast)/24);
R0 = roty(-90+phi)*rotz(-lambda)*G0*N0*P0;

N1 = rotx(-eps0-deleps)*rotz(-delpsi)*rotx(eps0);
P1 = rotz(-za)*roty(thetaa)*rotz(-zetaa);
G1 = rotz((360*gast)/24);
R1 = roty(90-phi)*rotz(lambda)*G1*N1*P1;

fak = 1;
% Systemtransformation, wir drehen das System und arbeiten mit
% Objektrotationsmatrizen, daher überall inverse Vorzeichen.
N0 = rotx(fak*(eps0+deleps))*rotz(fak*delpsi)*rotx(-fak*eps0);
P0 = rotz(fak*za)*roty(-fak*thetaa)*rotz(fak*zetaa);
G0 = rotz(-(360*gast)/24);
R0 = roty(-90+phi)*rotz(-lambda)*G0*N0*P0;

N1 = rotx(-fak*(eps0-deleps))*rotz(-fak*delpsi)*rotx(fak*eps0);
P1 = rotz(-fak*za)*roty(fak*thetaa)*rotz(-fak*zetaa);
G1 = rotz((360*gast)/24);
R1 = roty(90-phi)*rotz(lambda)*G1*N1*P1;