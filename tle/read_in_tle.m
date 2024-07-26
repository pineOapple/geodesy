function [tle] = read_in_tle(path_to_tle)
% STRONGLY INSPIRED BY 
    %READ_IN_TLE Summary of this function goes here
    % open file
    tle_file = fopen(path_to_tle, 'r');

    % read first line
    tline   = fgetl(tle_file);
    tle.norad   = tline(3:7);      			          % Satelliten ID (NORAD Katalog Nr.)
    tle.clas    = tline(8);					          % Klassifikation (U = öffentlich)
    tle.id      = tline(10:17);			              % Internationale Bezeichnung
    tle.epoch   = str2num(tline(19:32));              % Epoche 
    tle.year    = str2num(tline(19:20));              % Jahr
    tle.doy     = str2num(tline(21:32));              % Tag des Jahres (Float => Nachkommastellen in h, min, s umrechenbar)
    tle.bal_co  = str2num(tline(34:43));              % TD1 Ballistischer Koeefizient, erste (zeitliche) Ableitung der mittleren Bewegung n
    tle.TD2     = str2num(tline(45:50));              % 2. (zeitliche) Ableitung mittlere Bewegung n 
    tle.ExTD2   = str2num(tline(51:52));              % Exponent of 2nd Time Derivative
    tle.bstar   = str2num(tline(54:59));              % vereinfachter Strahlungsdruckterm B*
    tle.exp_bstar= str2num(tline(60:61));             % B*Exponent 
    tle.eph     = tline(63);                          % Ephemeris Type (0 = SGP4)
    tle.Enum    = str2num(tline(65:end-1))            % Element Number#
    tle.cs_1    = str2num(tline(end))                 % checksumme line1 modulo 10

    % read second line
    tline       = fgetl(tle_file);
    tle.i       = str2num(tline(9:16));                 % Inklination [°] 
    tle.raan    = str2num(tline(18:25));                % Rektazension des aufsteigenden Knotens (Right Ascension of Ascending Node) [°]
    tle.e       = str2num(strcat('0.',tline(27:33)));   % (numerische) Exzentrizität
    tle.omega   = str2num(tline(35:42));                % Argument des Perigäums [°]
    tle.M       = str2num(tline(44:51));                % Mittlere Anomalie [°]
    tle.n      = str2num(tline(53:63));                % Mittlere Bewegung = Anzahl an Orbits/Tag
    tle.rNo     = str2num(tline(65:end-1));             % Anzahl an Orbits seit Start (kann zu Overflows kommen)
    tle.cs_2    = str2num(tline(end))                   % checksumme line2 modulo 10

    fclose(tle_file);
    %print TLE
    tle
end

