function r_g = GDS_INT_TO_LCL(lambda, phi, y, m, d, ut1)
% GAST, NUTATION, PRECESSION, JULICANCENTURIES
[t, ~, ~] = gds.GDS_JULIANC(cyear, cmonth, cday, ut1);
[gast, ~] = gds.GDS_JUL_TO_GAST(ut1, t);
[eps0, deleps, delpsi] = gds.GDS_NUTATION(t);
[za, thetaa, zetaa] = gds.GDS_PRECESSION(t);
gds.
% POSITION OF OBJECT
[r_i(1),r_i(2),r_i(3)] = sph2cart(deg2rad(alpha), deg2rad(delta), radius);

% TRANSFORMATION
N = rotx(eps0+deleps)*rotz(delpsi)*rotx(-eps0);
P = rotz(za)*roty(-thetaa)*rotz(zetaa);
r_g = roty(-90+phi)*rotz(-lambda)*rotz(-(360*gast)/24)*N*P*r_i';