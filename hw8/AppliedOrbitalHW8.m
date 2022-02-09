clear;clc;
%Constants
mu = 3.986e14;
ae = 6378136.3; %m
J2 = 1.082e-3;

a = 1275e3+ae; %m
e = 0;
I = 81; %deg

oe0 = [a e I 0 0 0];

n = sqrt(mu/oe0(1)^3);
%Ground Station Geogrpahical Coordinates
gsLon1Lat2 = [30.2912, -97.7377;
            78.9282 11.8787;
            5.16866 -52.7750];
%Timeset
t = linspace(0,3600*24, 20000);

%Preallocation
RiseSet = zeros(1,6);
table = [0 0 0];
InSky = true;
MaxElevation = zeros(1,3);

for i = 1:(length(t)-1)
    %Define Satellite Orbital Elements
    M = n*t(i) + M_dot(a,e,I)*t(i);
    nu = E2nu(kepler(M,e),e); %in rad
    oe = [a e I/180*pi w_dot(a,e,I)*t(i) oe0(5)/180*pi-bigW_dot(a,e,I)*t(i) nu];
    rv = hw6oe2rv(oe,mu);
    rECI = rv(1:3)';
    
    %Define OE one step forward
    M = n*t(i+1) + M_dot(a,e,I)*t(i+1);
    nuNext = E2nu(kepler(M,e),e); %in rad
    oeNext = [a e I/180*pi w_dot(a,e,I)*t(i+1) oe0(5)/180*pi-bigW_dot(a,e,I)*t(i+1) nuNext];
    rv = hw6oe2rv(oeNext,mu);
    rECI_next = rv(1:3)';
    
    for gs = 1:3
        %Define slant range now and one step forward
        pSEZ_now = getSEZ(t(i), rECI, gsLon1Lat2(gs,1), gsLon1Lat2(gs,2));
        pSEZ_next = getSEZ(t(i+1), rECI_next, gsLon1Lat2(gs,1), gsLon1Lat2(gs,2));
        
        %Z-coord
        Z = pSEZ_now(3);
        Z_next = pSEZ_next(3);
        
        %Test to see if satellite crosses horizon:
        %Rise (Negative to Positive)
        if Z < 0 && Z_next > 0
            RiseSet(table(gs)+1, 2*gs-1) = t(i+1);
            MaxElevation(table(gs)+1,gs) = asind(Z/norm(pSEZ_now));
            InSky = true;
        %Set (Positive to Negative)
        elseif Z > 0 && Z_next < 0
            RiseSet(table(gs)+1, 2*gs) = t(i);
            table(gs) = table(gs)+1;
            MaxElevation(table(gs)+1,1) = 0;
        %InSky marked as no longer visible
            InSky = false;
        end
        %Determine Max Elevation
        if InSky && asind(Z/norm(pSEZ_now))>MaxElevation(table(gs)+1,gs)
            MaxElevation(table(gs)+1,gs) = asind(Z/norm(pSEZ_now));
        end
    end
end
%% Functions:
% SEZ Frame
function SEZ = getSEZ(t, r_ECI, gs_lat, gs_lon)
    ae = 6378136.3; %m
    we = 7.292115e-5; %rad/s
    rECEF = R3(we*t*180/pi)*r_ECI';
    rSEZ = R2(90-gs_lat)*R3(gs_lon)*rECEF;
    SEZ = rSEZ - ae*[cosd(gs_lon)*cosd(gs_lat);sind(gs_lon)*cosd(gs_lat);sind(gs_lat)];
end

% Rotation Matrices
function R1 = R1(theta)
    R1  = [1,       0,         0;
                0, cosd(theta),  sind(theta);
                0, -sind(theta), cosd(theta)];
end
function R2 = R2(theta)
    R2  = [cosd(theta),  0 , -sind(theta);
                    0,      1,      0 ;
               sind(theta),  0, cosd(theta)];
end
function R3 = R3(theta)
    R3  = [cosd(theta),   sind(theta) , 0;
                -sind(theta), cosd(theta),  0;
                    0,          0 ,       1];
end

%Precession Rates
function rate = bigW_dot(a,e,i) %i in degrees
    mu = 3.986e14;
    ae = 6378136.3; %m
    j2 = 1.082e-3;
    rate = -1.5*sqrt(mu/(a^3))*(ae/a)^2*j2*(1/sqrt(1-e^2))*cos((i*pi/180));
end
function rate = w_dot(a,e,i)
    mu = 3.986e14;
    ae = 6378136.3; %m
    j2 = 1.082e-3;
    rate = -.75*sqrt(mu/(a^3))*(ae/a)^2*j2*(1/(1-e^2)^2)*(1-5*(cos((i*pi/180))^2));
end
function rate = M_dot(a,e,i)
    mu = 3.986e14;
    ae = 6378136.3; %m
    j2 = 1.082e-3;
    rate = sqrt(mu/(a^3))*(1-.75*(ae/a)^2*j2*(1/sqrt((1-e^2)^3))*(1-3*(cos((i*pi/180))^2)));
end