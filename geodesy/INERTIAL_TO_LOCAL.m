function r_g = INERTIAL_TO_LOCAL(lambda, phi, cyear, cmonth, cday, ut1, alpha, delta, radius)
% GAST, NUTATION, PRECESSION, JULICANCENTURIES
[t, ~, ~] = GDS_JULIANC(cyear, cmonth, cday, ut1);
[gast, ~] = GDS_JUL_TO_GAST(ut1, t);
[eps0, deleps, delpsi] = GDS_NUTATION(t);
[za, thetaa, zetaa] = GDS_PRECESSION(t);

% POSITION OF OBJECT
[r_i(1),r_i(2),r_i(3)] = sph2cart(deg2rad(alpha), deg2rad(delta), radius);

% TRANSFORMATION
N = rotx(eps0+deleps)*rotz(delpsi)*rotx(-eps0);
P = rotz(za)*roty(-thetaa)*rotz(zetaa);
r_g = roty(-90+phi)*rotz(-lambda)*rotz(-(360*gast)/24)*N*P*r_i';