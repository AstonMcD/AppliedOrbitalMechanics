%% Applied Orbital Mechanics HW#1

%constants
mu = 3.986004415*10^14;
ae = 6378136.3;
we = 7.292115*10^-5;

% Satellite: GRACE Follow-On (google it, if you want more information)
% 
% Epoch: Specified in UTC Time System
% 	EPOCH1:   Year     Month      Day
% 	EPOCH2:   Hour     Minute     Seconds
% 
% POS (Position, ECI, in meters)          X       Y         Z
% VEL (Velocity, ECI, in meters/second)  Xdot     Ydot      Zdot
% 
%           EPOCH1                  2018.0                11.0                30.0
%           EPOCH2                    23.0                59.0                42.0
%           POS           5471639.55639308   -4009260.88949393   -1113125.19797190
%           VEL          -1210.66329250225     440.63844202350   -7515.04479126903

%% PROBLEM 1
% Calculate the (classical) Keplerian orbital elements at epoch (present angles in degrees).
rvGraceFO = [5471639.55639308   -4009260.88949393   -1113125.19797190, -1210.66329250225     440.63844202350   -7515.04479126903];
oeGFO = hw6rv2oe(rvGraceFO,mu);
oeGraceFOinDegrees = [oeGFO(1:2)' oeGFO(3:6)'.*180./pi]

%% PROBLEM 2

%% PROBLEM 3
instantAngVel=norm(rvGraceFO(4:6))/norm(rvGraceFO(1:3))
meanAngVel=sqrt(mu/oeGFO(1)^3)
    %instantaneous velocity will have minute variations on different instances due to the
    %eccentricity making it near circular but not perfect.

%% PROBLEM 4
rvRecalc = hw6oe2rv(oeGFO,mu);
difference = rvGraceFO' - rvRecalc
magRdiff = norm(difference(1:3))
magVdiff = norm(difference(4:6))

%% PROBLEM 5
specificEnergy = -mu/(2*oeGFO(1))
specificAngMom = cross(rvGraceFO(1:3), rvGraceFO(4:6))
angMomMag = norm(specificAngMom)

%% PROBLEM 6
%Propagation
T = 2*pi*sqrt(oeGFO(1)^3/mu);
Tquarter = T/4;

E0=2*atan(sqrt((1-oeGFO(2))/(1+oeGFO(2)))*tan(oeGFO(6)));
M=sqrt(mu/(oeGFO(1)^3))*Tquarter-2*pi+E0-oeGFO(2)*sin(E0);
Ef=hw7kepler(M,oeGFO(2));
nuF=2*atan(sqrt((1+oeGFO(2))-(1-oeGFO(2)))*tan(Ef/2));
nuSpaced=linspace(oeGFO(6),nuF,109)';

oe=zeros(109,6);
oe(:,1)=oeGFO(1);
oe(:,2)=oeGFO(2);
oe(:,3)=oeGFO(3);
oe(:,4)=oeGFO(4);
oe(:,5)=oeGFO(5);
oe(:,6)=nuSpaced;

rv = zeros(109,6);
for i=1:109
    rv(i,:)=hw6oe2rv(oe(i,:),mu);
end

options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
[t,rv6]=ode45(@(t,r) body2(t,r,mu),[0 Tquarter], rvGraceFO, options);
figure(1)
plot3(rv6(:,1),rv6(:,2),rv6(:,3),"r",rv(:,1),rv(:,2),rv(:,3),"b",0,0,0,"g")
grid
title("Quarter Period of Orbit Plot")
xlabel("x-position (m)")
ylabel("y-position (m)")
zlabel("z-position (m)")
legend("x","y","z")

%Energy
G=6.67408e-11;
mE=mu/G;
mg=600;
[k,u]=energy(rv,mE,mg,G)
E=k-u;
figure(2)
plot(t,E,"k")
title("Specific Energy Plot")
xlabel("Time (s)")
ylabel("Specific Energy (m^2/s^2)")

%Angular Momentum
h=angmom(rv6,mg);
figure(3)
plot(t,h(:,1),t,h(:,2),t,h(:,3))
title("Specific Angular Momentum Plot")
xlabel("Time (s)")
ylabel("Specific Angular Momentum (m^2/s^2)")
legend("x","y","z")

%% PROBLEM 7
r7=[5471639.55639308   -4009260.88949393   -1113125.19797190];
v7=[-1210.66329250225     440.63844202350   -7515.04479126903];

r = R2(23.5)*R3(144)*r7';
v = R2(23.5)*R3(144)*v7';

rvRot = [r' v'];
oeRot = hw6rv2oe(rvRot,mu);
oeRotinDegrees = [oeRot(1:2)' oeRot(3:6)'*180/pi]

function R = R3(theta)
    R = [cosd(theta),sind(theta),0; -sind(theta), cosd(theta), 0; 0, 0, 1];
end
function R = R2(theta)
    R = [cosd(theta), 0, -sind(theta);0, 1, 0; sind(theta), 0, cosd(theta)];
end