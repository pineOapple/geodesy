[R, G, t, mjd, jd, gast, gmst, eps0, deleps, delpsi, za, thetaa, zetaa] = GDS_INT_TO_LCL(lambda, phi, yr, m, d, ut1);

%% Position of the Sun
[aapp,dapp,r] = GDS_SOLARPOS(t);