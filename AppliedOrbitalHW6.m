%% Applied Orbital HW#6
%constants
mu = 3.986004415*10^14;
ae = 6378136.3;
we = 7.292115*10^-5;
g=9.81;
j2=1.082*10^-3;

%%Preallocations
n = zeros(1,4);
bigWbar = zeros(1,4);
wbar = zeros(1,4);
Mbar = zeros(1,4);
T_kepler = zeros(1,4);
T_anomalistic = zeros(1,4);
u = zeros(1,4);
T_nodal = zeros(1,4);
bigW = zeros(4,361);
beta = zeros(4,361);
%% Problem 1
% a. Calculate Precession rates in units of degrees/day
%       1       2       3       4
%   [Topex   Grace   ERS-1   Lageos]
a = [7705    6820    7156    12271]*1000; % m
e = [0.001   0.0016  0.001   0.004]; % 
i = [65.99   89.02   98.6    109.83]*pi/180; % rad

for j=1:4
n(j) = sqrt(mu/a(j)^3);
bigWbar(j) = -1.5*n(j)*(ae/a(j))^2*j2/sqrt(1-e(j)^2)*cos(i(j));
wbar(j) = -0.75*n(j)*(ae/a(j))^2*j2/(1-e(j)^2)^2*(1-5*(cos(i(j)))^2);
Mbar(j) = n(j)*(1-(0.75*(ae/a(j))^2*j2/(1-e(j)^2)^(3/2)*(1-3*(cos(i(j)))^2)));
end
%convert from rad/s to deg/day
Node_Precession_Rate = (bigWbar)*180/pi*60*60*24
Perigree_Precession_Rate = (wbar)*180/pi*60*60*24
Mean_Anomaly_Precession_Rate = (Mbar)*180/pi*60*60*24

% b. Calculate the Kepler, Anomalistic, and Draconitic/Nodal Periods
for j=1:4
    T_kepler(j) = 2*pi/n(j);
    T_anomalistic(j) = 2*pi/Mbar(j);
    u(j) = wbar(j) + Mbar(j);
    T_nodal(j) = 2*pi/u(j);
end
Keplerian_Period = T_kepler/60
Anomalistic_Period = T_anomalistic/60
Nodal_Period = T_nodal/60

% c.Calculate the Sun-Cycle Duration in units of days
bigWs = (360/365.2426); %deg/day
for j=1:4
Cs(j) = 360/(Node_Precession_Rate(j) - bigWs);
end
SunCycle_Duration = Cs

% d. Identify if Satellite is in sun-synchronous orbit
SunSyncTest = Node_Precession_Rate - bigWs
%Only the satellite ERS-1 (Satellite #3) is nearly sun-synchronous because
%its mean sun and nodal precession rate difference is near zero.

%% Problem 2.
% Calculate the B' angle at 1-day intervals and plot its variaton within
% one complete cycle of 360 degrees.
% Epoch Noon Jan 1, 2000
 P2epoch = [2000 3 20 0 0 0];
% %A. (10) Calculate the mean longitude of the sun (L)
 JD = juliandate(P2epoch);
 D = JD - 2451545.0; %in units of days
%Mean Longitude of Sun:

  
L = (0.9856474*D)*pi/180; % rad
Ldeg = rad2deg(L);
%B. Calculate the RA and Dec of the Sun
%Mean Anomaly:
g = pi/180*(357.528 + 0.9856003*D); % rad
%ecliptic longitude:
lambda = L + deg2rad(1.915)*sin(g) + deg2rad(0.020)*sin(2*g);
%Ecliptic Obliquity:
eps = (pi/180)*(23.439 - (4e-7)*D);
%Right Ascension: 
RA = atan(cos(eps)*tan(lambda));
RAdeg = rad2deg(RA)
Dec = asin(sin(eps)*sin(lambda));
Decdeg = rad2deg(Dec)
%C. Calculate GMST
GMST = 15*(pi/180)*(18.697374458 + 24.06570982*D); %in units of rad

% Use RA and Dec to calculate e_sun (remains constant)
e_sun = [sin(pi/2 - Dec)*cos(RA),
        sin(pi/2 - Dec)*cos(RA),
        cos(pi/2 - Dec)];
% e_sun1 = [cos(lambda),
%     cos(eps)*sin(lambda),
%     sin(eps)*sin(lambda)];

% To calculate e_h, need the change in RAAN
%   Change in time = sun cycle duration
bigWo = 0;
%Topex B'
dt1 = 0:(Cs(1)/360):Cs(1);
bigW(1,:) = bigWo + Node_Precession_Rate(1)*(dt1); % degrees
for j=1:length(dt1)
eh_Topex = [sind(bigW(1,j))*sin(i(1)),
       -1*cosd(bigW(1,j))*sin(i(1)),
       cos(i(1))];
beta(1,j) = 90 - acosd(dot(e_sun1,eh_Topex)); % in degrees
end
%GRACE B'
dt2 = 0:(Cs(2)/360):Cs(2);
bigW(2,:) = bigWo + Node_Precession_Rate(2)*(dt2); % degrees
for j=1:361
eh_Grace = [sind(bigW(2,j))*sin(i(2)),
       -1*cosd(bigW(2,j))*sin(i(2)),
       cos(i(2))];
beta(2,j) = 90 - acosd(dot(e_sun1,eh_Grace)); % in degrees
end
%Lageos B'
dt4 = 0:(Cs(4)/360):Cs(4);
bigW(4,:) = bigWo + Node_Precession_Rate(4)*(dt4); % degrees
for j=1:361
eh_Lageos = [sind(bigW(4,j))*sin(i(4)),
       -1*cosd(bigW(4,j))*sin(i(4)),
       cos(i(4))];
beta(4,j) = 90 - acosd(dot(e_sun1,eh_Lageos)); % in degrees
end
Cycleof360 = 1:361;
figure(1)
plot(Cycleof360,beta(1,:),Cycleof360,beta(2,:),Cycleof360,beta(4,:))
title("B' of satellites vs 360 degree cycle")
legend("B'-Topex","B'-GRACE","B'-Lageos")
xlabel("360 Degree Cycle (deg)")
ylabel("B' (deg)")
grid on
figure(2)
plot(dt1,beta(1,:),dt2,beta(2,:),dt4,beta(4,:))
title("B' of satellites vs Sun-cycle Duration")
legend("B'-Topex","B'-GRACE","B'-Lageos")
xlabel("Sun-cycle Duration (days)")
ylabel("B' (deg)")
grid on