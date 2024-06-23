function [aapp,dapp,r,lam,beta] = lunarpos(t)

% LUNARPOS calculates the apparent position of the Moon 
%          referring to the mean equinox of date
%
% HOW [aapp,dapp,r,lam,beta] = lunarpos(t)
%
% IN  t    - time relative to J2000.0 [julian centuries (of 36525 days)]
% OUT aapp - right ascension [deg] 
%     dapp - declination [deg]
%     r    - distance  [m]
%     lam  - ecliptical longitude [deg]
%     beta - ecliptical latitude [deg]
%
% NB  accuracy: 10" in longitude, 4" in latitude
%
% Nico Sneeuw			Munich				01/98
% Raymond Tsoi          Calgary             05/03

% uses NUTATION, LUNARPOSSUB
%
% revision history:
%  - 980120NS: thorough revision, vectorization and translation of original
%              MOONPOS.M by Roland Schmidt 07/97
%  - 981204NS: output wrong for vectorial input t. correct for scalar input.
%              A for-loop wrap is introduced as simple bug-fix. Ugly, though.
%  - 2001.11.NS: nutwink -> nutation
%  - 2003.05.RT: subroutine, lunarpossub.dll, has been implemented to improve speed performance
%                tested. Differences below mas level.


% Modified version
[eps0,deps,dpsi] = nutation(t);
[aapp,dapp,r,lam,beta] = lunarpossub(t,eps0,deps,dpsi);


