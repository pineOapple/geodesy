function [aapp,dapp,r] = solarpos(t)

% SOLARPOS calculates the apparent position of the Sun 
%          referring to the true equinox of date
%
% HOW [aapp,dapp,r] = solarpos(t)
%
% IN  t    - time relative to J2000.0 [julian centuries (of 36525 days)]
% OUT aapp - apparent right ascension [deg]
%     dapp - apparent declination [deg]
%     r    - distance  [AU]
%
% NB relatively low accuracy: 0.01 [degree] 
%
% Nico Sneeuw			Munich				01/98

% uses NUTATION
%
% revision history:
%  - 980120NS: thorough revision, vectorization and translation of original
%              HELPOS by Roland Schmidt 07/97
%  - 2001.11.NS: nutwink -> nutation


% mean longitude relative to mean equinox of date
Lo   = 280.46645 + 36000.76983*t + 0.0003032*t.^2;	  % [deg]

% mean anomaly of Sun
M    = 357.52910 + 35999.05030*t - 0.0001559*t.^2;	  % [deg]
M    = M*pi/180;                                      % [rad]

% eccentricity of Earth orbit
e    = 0.016708617 - 0.000042037*t - 0.0000001236*t.^2;

% Sun's equation of the center
C    = (1.9146-0.004817*t-0.000014*t.^2).*sin(M) + ...
       (0.019993-0.000101*t).*sin(2*M) + ...
        0.00029*sin(3*M);                             % [deg]

% true longitude and true anomaly of Sun 
Ltru = Lo + C;                                        % [deg]
v    = M + C*pi/180;	                              % [rad]

% solar radiusvector
r    = 1.000001018*(1-e.^2)./(1+e.*cos(v));           % [AU]

% apparent longitude, relative to true equinox of date
om   = 125.04 - 1934.136*t;                           % [deg]
om   = om * pi/180;                                   % [rad]
Lapp = Ltru - 0.00569 - 0.00478*sin(om);              % [deg]
Lapp = Lapp * pi/180;                                 % [rad]

% obliquity of ecliptic
[eps0,deleps,delpsi] = nutation(t);                   % [deg]
eps  = ( eps0 + 0.00256*cos(om) ) *pi/180;            % [rad]

% apparent right ascension and declination
aapp = atan2((cos(eps).*sin(Lapp)),cos(Lapp));        % [rad]
aapp = rem(aapp*180/pi+360,360);	                  % [deg]
dapp = asin(sin(eps).*sin(Lapp))*180/pi;              % [deg]
