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

% uses NUTATION, slsr.mat, sb.mat
%
% revision history:
%  - 980120NS: thorough revision, vectorization and translation of original
%              MOONPOS.M by Roland Schmidt 07/97
%  - 981204NS: output wrong for vectorial input t. correct for scalar input.
%              A for-loop wrap is introduced as simple bug-fix. Ugly, though.
%  - 2001.11.NS: nutwink -> nutation

% preliminaries
[tr,tc] = size(t);
%t  = t(:);

% bug-fix 981204NS: initialize output-variables and wrap the original m-code 
%                   into for-loop
aapptmp = zeros(size(t));
dapptmp = zeros(size(t));
rtmp    = zeros(size(t));
lamtmp  = zeros(size(t));
betatmp = zeros(size(t));
ttmp    = t;
for i = 1:length(ttmp)
   t = ttmp(i);
%---------------------------------------------------------------------------   
   
   
%-------------------------------------------------------
% the angles Ls,D,M,Ms,F from Meeus S.336/337, all in [rad]
%-------------------------------------------------------
tt = [t.^4 t.^3 t.^2 t ones(size(t))] *pi/180;		% sort of Vandermonde
Ls = tt*[-1/65194000   1/538841   -0.0013268 481267.88134236 218.3164591]';
D  = tt*[-1/113065000  1/545868   -0.0016300 445267.1115168  297.8502042]';
M  = tt*[ 0            1/24490000 -0.0001536  35999.0502909  357.5291092]';
Ms = tt*[-1/14712000   1/69699     0.0089970 477198.8676313  134.9634114]';
F  = tt*[ 1/863310000 -1/3526000  -0.0034029 483202.0175273   93.2720993]';

A1 = (119.75 +    131.849*t) * pi/180;
A2 = ( 53.09 + 479264.290*t) * pi/180;
A3 = (313.45 + 481266.484*t) * pi/180;


%----------------------------------------------------
% computation of the sums sigl, sigr, sigb
% slsr.mat and sb.mat contain tables from Meeus
%
%----------------------------------------------------
% periodic terms for long. and dist. to Moon
load slsr
[rr,cc]=size(slsr);

E  = 1 - 0.002516*t - 0.0000074*t.^2;                 % eccentricity corr.
Eq = E.^2;
i1 = find(abs(slsr(:,2)) == 1);
i2 = find(abs(slsr(:,2)) == 2);

scarg = slsr(:,1:4)*[D M Ms F]';
scoef = slsr(:,5) * ones(1,length(t));                % sine coefficient
ccoef = slsr(:,6) * ones(1,length(t));	               % cosine coefficient
scoef(i1,:) = scoef(i1,:).*(ones(length(i1),1)*E');	% correction for
ccoef(i1,:) = ccoef(i1,:).*(ones(length(i1),1)*E');	%   changing ecc.
scoef(i2,:) = scoef(i2,:).*(ones(length(i2),1)*Eq');	%   of earth orbit.
ccoef(i2,:) = ccoef(i2,:).*(ones(length(i2),1)*Eq');

sigl = sum(scoef.*sin(scarg))';
sigl = sum(sigl) + 3958*sin(A1) + 1962*sin(Ls-F) + 318*sin(A2);
sigl = sigl/1000000;                                  % [deg]

sigr = sum(ccoef.*cos(scarg))';
sigr = sigr/1000;	                                    % [km]

% periodic terms for latitude of Moon
load sb
[rr,cc]=size(sb);

scarg = sb(:,1:4)*[D M Ms F]';
scoef = sb(:,5) * ones(1,length(t));                  % sine coefficient
scoef(i1,:) = scoef(i1,:).*(ones(length(i1),1)*E');	% correction for ecc.
scoef(i2,:) = scoef(i2,:).*(ones(length(i2),1)*Eq');	%   of earth orbit

sigb = sum(scoef.*sin(scarg))';
sigb = sum(sigb) - 2235*sin(Ls) + 382*sin(A3) + 175*sin(A1-F) + ...
       175*sin(A1+F) + 127*sin(Ls-Ms) - 115*sin(Ls+Ms);
sigb = sigb/1000000;	                                 % [deg]


%-------------------------
% coordinates lam,beta,r
%-------------------------
lam  = Ls*180/pi+sigl;                                % [deg]
lam  = rem(rem(lam,360)+360,360);                     % (0-360)
beta = sigb;                                          % [deg]
r    = (385000.56 + sigr)*1000;                       % [m]


%--------------------------------------------
% apparent right ascension and declination
%--------------------------------------------

% nutation in longitude (accuracy ~ 1")
% in Meeus more exact formulae (e.g. p.339)

[eps0,deleps,delpsi] = nutation(t);
lams = (lam+delpsi)*pi/180;                           % apparent longitude
eps  = (eps0+deleps)*pi/180;                          % true obl. of ecliptic
b2   = beta*pi/180;                                   % [rad]
aapp = atan2((sin(lams).*cos(eps)-tan(b2).*sin(eps)),cos(lams));
dapp = asin(sin(b2).*cos(eps)+cos(b2).*sin(eps).*sin(lams));
aapp = rem(aapp*180/pi+360,360);
dapp = dapp*180/pi;


%--------------------------------------------
% reshape output according to [tr,tc]
%--------------------------------------------
%lam  = reshape(lam,tr,tc);
%beta = reshape(beta,tr,tc);
%r    = reshape(r,tr,tc);
%aapp = reshape(aapp,tr,tc);
%dapp = reshape(dapp,tr,tc);


%---------------------------------------------------------------------------   
% here we continue with bug-fix 981204NS
aapptmp(i)=aapp;
dapptmp(i)=dapp;
rtmp(i)   =r;
lamtmp(i) =lam;
betatmp(i)=beta;
end
aapp = aapptmp;
dapp = dapptmp;
r    = rtmp;
lam  = lamtmp;
beta = betatmp;