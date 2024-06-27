function [dalpha,ddelta] = parallax(alpha,delta,y,m,d,ut1,phi,lambda)

% PARALLAX computes corrections to geocentric lunar coordinates, due to parallax,
%          i.e. due to the difference geocentric <--> topocentric
%          (depending on observer position and time)
%
% HOW [dalpha,ddelta] = parallax(alpha,delta,y,m,d,ut1,phi,lambda)
% IN  alpha  - right ascension of Moon        [deg,vector]
%     delta  - declination     "   "          [deg,vector]
%     y      - year (4-digit!)
%     m      - month
%     d      - day
%     ut1    - universal time Greenwich       [decimal hours]
%     phi    - latitude of observer           [deg,skalar]
%     lambda - longitude "    "               [deg,skalar]
% OUT dalph  - correction to right ascension  [deg,vektor]
%     ddelta - correction to declination      [deg,vektor]


%----------------------------------------------------------------------------
% Sneeuw/Zebhauser 			            Muenchen			            12/98 
%----------------------------------------------------------------------------
% uses JULIAN2000, JUL2GAST
%
%----------------------------------------------------------------------------
% revision history
% .Oct. 2001: translation
%----------------------------------------------------------------------------
% remarks
% .Formulae from Meeus, chapter 39
% .Might be adapted to other celestial bodies by changing the variable rm.
% .Caveats:
%  - rm could be more precise
%  - Earth flattening not considered
%----------------------------------------------------------------------------

% preliminaries
error(nargchk(8,8,nargin))
rm     = 0.00258;							% distance to Moon [AU]
sinpi  = sin(8.794/3600*pi/180)/rm;	% sin(maximum parallax), about 1 deg.

% time computations
T      = julian2000(y,m,d,ut1);     % time since epoch J2000 [jul. centuries] 
gast   = jul2gast(ut1,T);           % GAST [dec. hours]
tm     = 15*gast + lambda - alpha;  % hour angle of the Moon [deg]

% deg -> rad
phi    = phi/180*pi;
lambda = lambda/180*pi;
tm     = tm/180*pi;
delta  = delta/180*pi;

% here we go
denom  = cos(delta)-cos(phi)*sinpi*cos(tm);
dalpha = atan(-cos(phi)*sinpi*sin(tm)./denom);	% correction alpha [rad]
cda    = cos(dalpha);
deltap = atan2((sin(delta)-sin(phi)*sinpi).*cda,denom);	% delta prime [rad]

% Fazit [deg]
dalpha = dalpha*180/pi;
ddelta = (deltap-delta)*180/pi;

