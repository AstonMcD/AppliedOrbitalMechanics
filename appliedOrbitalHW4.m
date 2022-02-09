%% Applied Orbital Mechanics HW#3

%constants
mu = 3.986004415*10^14;
ae = 6378136.3;
we = 7.292115*10^-5;
g=9.81;
j2=1.082*10^-3;
% Preallocating

%% Intro Functions for syntax and copying:
%fd = 0.5*rho*vr^2*Cd*A/m;
% v = vr + cross(we, r);
% vwind = vr + cross(we, r) + vw;
% dadt = -2*vr*a^2*fd/mu;
% dedt = -2*(e+cos(f))*fd/vr;

%% 1. Orders of Magnitude:
%A
i = deg2rad(66);
Cd = 2.0;
A = 1.3; %(m^2)
m = 350; %(kg)
rho = 1.6*10^-12; %(kg/m^3)
alt = 450*1000; %(m)
a = ae + alt;
e = 0; %circular orbit
v_1a = sqrt(mu/a);
fd_1a = 0.5*rho*Cd*A*(v_1a^2)/m;
dadt_1a = -2*v_1a*a^2*fd_1a/mu;
dadt_1a_mday = dadt_1a*86400;
%B
f = 0:3:360;
for j = 1:length(f)
   
oe = [a e i 0 0 f(j)/180*pi];
rv = hw6oe2rv(oe, mu)';
v_1b = norm(rv(4:6) - cross([0 0 we],rv(1:3)));
fd_1b = 0.5*rho*Cd*A*(v_1b^2)/m;
dadt_1b_mday(j) = 86400*-2*v_1b*(a^2)*fd_1b/mu;

end
figure(1)
plot(f,dadt_1b_mday)
title("Semimajor Axis vs. True Anomaly")
xlabel("True Anomaly (degrees)")
ylabel("Semimajor Axis Decay Rate (m/day)")

%C
for j = 1:length(f)
   
oe = [a e i 0 0 f(j)/180*pi];
rv = hw6oe2rv(oe, mu)';
v_1b = norm(rv(4:6) - cross([0 0 we],rv(1:3)));
v_1c = v_1b - 500*(v_1b/norm(v_1b));
fd_1c = 0.5*rho*Cd*A*(v_1c^2)/m;
dadt_1b_mday(j) = 86400*-2*v_1c*(a^2)*fd_1c/mu;

end
figure(2)
plot(f,dadt_1b_mday)
title("Semimajor Axis vs. True Anomaly with Wind Gust")
xlabel("True Anomaly (degrees)")
ylabel("Semimajor Axis Decay Rate (m/day)")

%% 2. Eccentric Orbits
e2 = 0.032;
a2i = 6928*1000;
Cd2 = 2.0;
A2 = 0.8;
m2 = 500;

e2 =zeros(1,121);
a2 = zeros(1,121);
e2(1) = 0.032;
a2(1) = a2i;
for k = 1:length(f)
    p = a2(k)*(1-e2(k)^2);
    v = mu*(1+e2(k)+2*e2(k)*cos(f(k)/180*pi))/p;
    fd_2a = 0.5*rho*Cd*A*(v^2)/m;
dadt_2a = -2*v*a^2*fd_2a/mu;
    dedt = -2*(e2(k) + cos(f(k)/180*pi))*fd_2a;
    if k > 1
    a2(k) = a2(k-1) + dadt_2a;
    e2(k) = e2(k-1) +dedt;
    end
end
e2(1) = 0.032;
a2(1) = a2i;


figure(3)
plot(f,a2,f,e2)

% Cd = 2;
% A = .8; %m^2
% m = 500; %kg
% rho = 1.6e-12; %kg/m^3
% mu = 3.986004415e14; %m^3/s^2
% r_e = 6378136.3; %m
% a = 6928000; %m
% e = .032;
% 
% 
% resolution = 1; %deg
% nu = 0:resolution:360-resolution; %deg
% a_plot = [a];
% e_plot = [e];
% alt_plot = [a-r_e];
% t_prev = 0;
% 
% for i = 2:length(nu)
%     p = a_plot(i-1)*(1-e_plot(i-1)^2);
%     n = sqrt(mu/a_plot(i-1)^3);
%     v = sqrt(mu/p*(1+e_plot(i-1)^2+2*e_plot(i-1)*cosd(nu(i))));
%     dadt = -2*v*a_plot(i-1)^2/mu*(.5*rho*v^2*A*Cd/m);
%     dedt = -2/v*(.5*rho*v^2*A*Cd/m)*(e_plot(i-1)+cosd(nu(i)));
%     E = nu2E(nu(i)/180*pi,e_plot(i-1));
%     t = (E-e_plot(i-1)*sin(E))/n;
%     dt = t-t_prev;
%     t_prev = t;
%     a_plot(i) = a_plot(i-1) + dadt*dt;
%     e_plot(i) = e_plot(i-1) + dedt*dt;
%     r = (a_plot(i)*(1-e_plot(i)^2))/(1+e_plot(i)*cosd(nu(i)));
%     alt_plot(i) = r-r_e;
% end
% figure(3)
% plot(nu,alt_plot);
% xlabel("True Anomaly (deg)");
% ylabel("Altitude (m)");
% figure(4);
% plot(nu,e_plot);    
% xlabel("True Anomaly (deg)");
% ylabel("Eccentricity");

    