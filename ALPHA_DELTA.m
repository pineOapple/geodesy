function [delta, alpha] = ALPHA_DELTA(year)
    % Function to calculate the declination (delta) and right ascension (alpha)
    % of the Sun for each day of a given year.
    
    % Number of days in the year
    if mod(year, 4) == 0 && (mod(year, 100) ~= 0 || mod(year, 400) == 0)
        daysInYear = 366; % Leap year
    else
        daysInYear = 365; % Non-leap year
    end
    
    % Allocate arrays for declination and right ascension
    delta = zeros(1, daysInYear);
    alpha = zeros(1, daysInYear);
    
    % Loop over each day of the year
    for day = 1:daysInYear
        % Julian date calculation
        JD = 367 * year - floor(7 * (year + floor((day + 9) / 12)) / 4) + floor(275 * day / 9) + 1721013.5;
        
        % Number of centuries since J2000.0
        T = (JD - 2451545.0) / 36525;
        
        % Mean longitude of the Sun
        L0 = mod(280.46646 + 36000.76983 * T + 0.0003032 * T^2, 360);
        
        % Mean anomaly of the Sun
        M = 357.52911 + 35999.05029 * T - 0.0001537 * T^2;
        
        % Sun's equation of center
        C = (1.914602 - 0.004817 * T - 0.000014 * T^2) * sind(M) + (0.019993 - 0.000101 * T) * sind(2 * M) + 0.000289 * sind(3 * M);
        
        % True longitude of the Sun
        lambda = L0 + C;
        
        % Obliquity of the ecliptic
        epsilon = 23.439292 - 0.000013 * T;
        
        % Declination of the Sun
        delta(day) = asind(sind(epsilon) * sind(lambda));
        
        % Right ascension of the Sun
        alpha(day) = mod(atan2d(cosd(epsilon) * sind(lambda), cosd(lambda)), 360);
    end
end
