%% Applied Orbital Homework #5
% Satellite: GRACE Follow-On (google it, if you want more information)
% 
% Epoch: Specified in UTC Time System
% 	EPOCH1:   Year     Month      Day
% 	EPOCH2:   Hour     Minute     Seconds
% 
% POS (Position, ECI, in meters)          X       Y         Z
% VEL (Velocity, ECI, in meters/second)  Xdot     Ydot      Zdot

EPOCH1 = [2018.0, 11.0, 30.0];
EPOCH2 = [23.0, 59.0, 42.0];

POS = [5471639.55639308, -4009260.88949393, -1113125.19797190];
VEL = [-1210.66329250225, 440.63844202350, -7515.04479126903];
rhat = POS/norm(POS);
%constants
mu = 3.986004415*10^14;
ae = 6378136.3;
we = 7.292115*10^-5;
grav=9.81;
j2=1.082*10^-3;

%% Problem 1. (10 points)
% converting earth rotation 1c to degree


%% Problem 2. (50 points)
% Epoch Noon Jan 1, 2000
P2epoch = [2000 1 1 12 0 0];
P3epoch = [2021 2 15 21 0 0];
P4epoch = [2018.0, 11.0, 30.0, 23.0, 59.0, 42.0];
%A. (10) Calculate the mean longitude of the sun (L)
JD = juliandate(P4epoch);
D = JD - 2451545.0; %in units of days
%Mean Longitude of Sun:
L = (pi/180)*(280.460 + 0.9856474*D) % rad
Ldeg = rad2deg(L);
%B. Calculate the RA and Dec of the Sun
%Mean Anomaly:
g = pi/180*(357.528 + 0.9856003*D); % rad
%ecliptic longitude:
lambda = L + deg2rad(1.915)*sin(g) + deg2rad(0.020)*sin(2*g);
%Ecliptic Obliquity:
eps = (pi/180)*(23.439 - (4e-7)*D);
%Right Ascension: 
RA = atan(cos(eps)*tan(lambda))
RAdeg = rad2deg(RA)
Dec = asin(sin(eps)*sin(lambda))
Decdeg = rad2deg(Dec)
%C. Calculate GMST
GMST = 15*(pi/180)*(18.697374458 + 24.06570982*D); %in units of rad


%% Problem 4.
% Calculate the (classical) Keplerian orbital elements at epoch (present angles in degrees).
rvGraceFO = [5471639.55639308   -4009260.88949393   -1113125.19797190, -1210.66329250225     440.63844202350   -7515.04479126903];
oeGFO = hw6rv2oe(rvGraceFO,mu);
oeDeg = [oeGFO(1:2)' oeGFO(3:6)'.*180./pi]

sunECIhat = [sin(111.7549)*cos(67.0102), sin(111.7549)*sin(67.0102), cos(111.7549)];
R3CapOmega = [cos(oeDeg(5)) sin(oeDeg(5)) 0; -sin(oeDeg(5)) cos(oeDeg(5)) 0; 0 0 1];
R1Inclination = [1 0 0; 0 cos(oeDeg(3)) sin(oeDeg(3)); 0 -sin(oeDeg(3)) cos(oeDeg(3))];

sunOrbitHat = R3CapOmega*R1Inclination*sunECIhat'
