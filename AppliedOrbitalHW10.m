clear all, clc
%% Applied Orbital HW#10
%constants
mu = 3.986004415*10^14;
ae = 6378136.3;
we = 7.292115*10^-5;
g=9.81;
j2=1.082*10^-3;
d2r = pi/180;
r2d2 = 180/pi; %r2d2 cause why not lol (rad2deg)
%% Target Satellite "True" Orbit:
%Set up
t = 0:(60*60*24); %duration of one day
rtgt = 6778*1000;
w = sqrt(mu/rtgt^3);
e = 0; %circular orbit
i = 90*d2r; %polar orbit
% 2 body propagation
OEtgt=[rtgt e i 0 0 0];
RVtgt= hw6oe2rv(OEtgt, mu);
options = odeset('RelTol', 1e-14, 'AbsTol', 1e-14);
[t,RVtgt_prop]=ode45(@(t,r) body2(t,r,mu),[0:86400], RVtgt, options);
%% CASE 1: "True" Orbit
rint = rtgt-100;
eint = 0;
iint = 90;
% 2 body propagation
OEint1 =[rint e i 0 0 0];
RVint1 = hw6oe2rv(OEint1, mu);
options = odeset('RelTol', 1e-14, 'AbsTol', 1e-14);
[t,RVint1_prop]=ode45(@(t,r) body2(t,r,mu),[0:86400], RVint1, options);
%% "True" Orbit difference in ECI frame
RVdiff1 = RVint1_prop - RVtgt_prop;
figure(1)
plot3(RVdiff1(:,1),RVdiff1(:,2),RVdiff1(:,3))
title('CASE 1: 3D difference between "True" Orbits [ECI frame]')

posMagListTrue1 = zeros(1,length(t));
posVecVecTrue1 = [RVdiff1(:,1),RVdiff1(:,2),RVdiff1(:,3)];
for j = 1:length(t)
    posMagListTrue1(j) = norm(posVecVecTrue1(j,:));
end
figure(2)
plot(t,posMagListTrue1);
title("CASE 1: Distance Magnitude difference for True Orbits [ECI frame]")
xlabel("Time (sec)")
ylabel("Distance (m)")
%% RTN Conversion
R_rel1 = zeros(3,length(t));
Rt = RVtgt_prop(:,1:3)';
Vt = RVtgt_prop(:,4:6)';
Rc1 = RVint1_prop(:,1:3)';
for k= 1:length(t)
    rh = Rt(:,k) / norm(Rt(:,k));
    H = cross(Rt(:,k),Vt(:,k));
    Nh = H/norm(H);
    Th = cross(Nh,rh);
    RTNtoECI = [rh, Th, Nh];
    R_rel1(:,k) = transpose(RTNtoECI) * (Rc1(:,k) - Rt(:,k));
    %relative position of chase vehicle wrt to tgt vehicle
end
figure(3)
plot3(R_rel1(1,:),R_rel1(2,:),R_rel1(3,:))
title('CASE 1: 3D Relative Position [RTN frame]')

posMagListRel1 = zeros(1,length(t));
posVecVecRel1 = [R_rel1(1,:),R_rel1(2,:),R_rel1(3,:)];
for j = 1:length(t)
    posMagListRel1(j) = norm(posVecVecRel1(:,j));
end
figure(4)
plot(t,posMagListRel1);
title("CASE 1: Relative distance in [RTN frame]")
xlabel("Time (sec)")
ylabel("Distance (m)")   
%% CASE 1: CW Equations
%vrel = (rint*sqrt(mu/rint^3))-(rtgt*w);
xo = -100;
yo = 0;
zo =0;
vxo = 0;
vyo = 0;
vzo = 0;
xt = (vxo/w)*sin(w*t)-(3*xo+(2*vyo/w))*cos(w*t)+(4*xo+(2*vyo/w));
yt = (6*xo+(4*vyo/w))*sin(w*t)+(2*vxo/w)*cos(w*t)-(6*w*xo+3*vyo)*t+(yo-(2*vxo/w));
zt = -zo*w*sin(w*t)+vzo*cos(w*t);
vxt = vxo*cos(w*t)+(3*w*xo+2*vyo)*sin(w*t);
vyt = (6*w*xo+4*vyo)*cos(w*t)-2*vxo*sin(w*t)-(6*w*xo+3*vyo);
vzt = -zo*w*sin(w*t)+vzo*cos(w*t);
% Grows infinitely from secular growth.
figure(5)
plot3(xt,yt,zt)
title("CASE 1: 3D position from CW Equations for one day duration")

posMagListCW1 = zeros(1,length(t));
posVecVecCW1 = [xt,yt,zt];
for j = 1:length(t)
    posMagListCW1(j) = norm(posVecVecCW1(j,:));
end
figure(6)
plot(t,posMagListCW1);
title("CASE 1: Distance Magnitude for CW Equations")
xlabel("Time (sec)")
ylabel("Distance (m)")
%% CASE 1: Comparisons
% Rcomp1=abs(posMagListRel1-posMagListCW1);
% figure(7)
% plot(t,Rcomp1);
% Rcomp1_3D = R_rel1 - posVecVecCW1';
% figure(8)
% plot3(Rcomp1_3D(:,1),Rcomp1_3D(:,2),Rcomp1_3D(:,3))

%The "true" orbit 3D difference demonstrates a beautiful spiral which
%aligns with the CW equation found to have oscillation and constant secular
%growth in the y direction. Their magnitudes increase continously because
%of this secular growth with the "True" orbits growth entirely linear. The
%CW equation oscillates but obviosuly trends upwards with linear growth. It
%has about 45000 meters greater magnitude over the duration of a day.

%% CASE 2: "True" Orbit
M2 = (360-0.00085)*d2r;
nu2 = E2nu(kepler(M2,e),e); %in rad
OEint2 = [rtgt e i 0 0 nu2];
RVint2 = hw6oe2rv(OEint2, mu);
options = odeset('RelTol', 1e-14, 'AbsTol', 1e-14);
[t,RVint2_prop]=ode45(@(t,r) body2(t,r,mu),[0:86400], RVint2, options);
R_rel2 = zeros(3,length(t));
Rc2 = RVint2_prop(:,1:3)';
for k= 1:length(t)
    rh = Rt(:,k) / norm(Rt(:,k));
    H = cross(Rt(:,k),Vt(:,k));
    Nh = H/norm(H);
    Th = cross(Nh,rh);
    RTNtoECI = [rh, Th, Nh];
    R_rel2(:,k) = transpose(RTNtoECI) * (Rc2(:,k) - Rt(:,k));
    %relative position of chase vehicle wrt to tgt vehicle
end
figure(7)
plot3(R_rel2(1,:),R_rel2(2,:),R_rel2(3,:))
title('CASE 2: 3D Relative Position [RTN frame]')

posMagListRel2 = zeros(1,length(t));
posVecVecRel2 = [R_rel2(1,:),R_rel2(2,:),R_rel2(3,:)];
for j = 1:length(t)
    posMagListRel2(j) = norm(posVecVecRel2(:,j));
end
figure(8)
plot(t,posMagListRel2);
title("CASE 2: Relative distance in [RTN frame]")
xlabel("Time (sec)")
ylabel("Distance (m)")   
%% CASE 2: CW Equations
xo2 = 0;
yo2 = -0.00085*d2r*rtgt;
zo2 = 0;
vxo2 = 0;
vyo2 = 0;
vzo2 = 0;
xt2 = (vxo2/w)*sin(w*t)-(3*xo2+(2*vyo2/w))*cos(w*t)+(4*xo2+(2*vyo2/w));
yt2 = (6*xo2+(4*vyo2/w))*sin(w*t)+(2*vxo2/w)*cos(w*t)-(6*w*xo2+3*vyo2)*t+(yo2-(2*vxo2/w));
zt2 = -zo2*w*sin(w*t)+vzo2*cos(w*t);
vxt2 = vxo2*cos(w*t)+(3*w*xo2+2*vyo2)*sin(w*t);
vyt2 = (6*w*xo2+4*vyo2)*cos(w*t)-2*vxo2*sin(w*t)-(6*w*xo2+3*vyo2);
vzt2 = -zo2*w*sin(w*t)+vzo2*cos(w*t);
% Grows infinitely from secular growth.
figure(9)
plot3(xt2,yt2,zt2)
title("CASE 2: 3D position from CW Equations for one day duration")

posMagListCW2 = zeros(1,length(t));
posVecVecCW2 = [xt2,yt2,zt2];
for j = 1:length(t)
    posMagListCW2(j) = norm(posVecVecCW2(j,:));
end
figure(10)
plot(t,posMagListCW2);
title("CASE 2: Distance Magnitude for CW Equations")
xlabel("Time (sec)")
ylabel("Distance (m)")
%% CASE 3: "True" Orbit
OEint3 = [rtgt e i 0 0 0];
RVint3 = hw6oe2rv(OEint3, mu);
RVint3 = [RVint3(1) RVint3(2)+100 RVint3(3) RVint3(4) RVint3(5) RVint3(6)]; 
%approx seperated by 100m along positive orbit normal?
[t,RVint3_prop]=ode45(@(t,r) body2(t,r,mu),[0:86400], RVint3, options);
R_rel3 = zeros(3,length(t));
Rc3 = RVint1_prop(:,1:3)';
for k= 1:length(t)
    rh = Rt(:,k) / norm(Rt(:,k));
    H = cross(Rt(:,k),Vt(:,k));
    Nh = H/norm(H);
    Th = cross(Nh,rh);
    RTNtoECI = [rh, Th, Nh];
    R_rel3(:,k) = transpose(RTNtoECI) * (Rc3(:,k) - Rt(:,k));
    %relative position of chase vehicle wrt to tgt vehicle
end
figure(11)
plot3(R_rel3(1,:),R_rel3(2,:),R_rel3(3,:))
title('CASE 3: 3D Relative Position [RTN frame]')

posMagListRel3 = zeros(1,length(t));
posVecVecRel3 = [R_rel3(1,:),R_rel3(2,:),R_rel3(3,:)];
for j = 1:length(t)
    posMagListRel3(j) = norm(posVecVecRel3(:,j));
end
figure(12)
plot(t,posMagListRel3);
title("CASE 3: Relative distance in [RTN frame]")
xlabel("Time (sec)")
ylabel("Distance (m)")   

%% CASE 3: CW Equations
xo3 = 0;
yo3 = 0;
zo3 = 100;
vxo3 = 0;
vyo3 = 0;
vzo3 = 0;
xt3 = (vxo3/w)*sin(w*t)-(3*xo3+(2*vyo3/w))*cos(w*t)+(4*xo3+(2*vyo3/w));
yt3 = (6*xo3+(4*vyo3/w))*sin(w*t)+(2*vxo3/w)*cos(w*t)-(6*w*xo3+3*vyo3)*t+(yo3-(2*vxo3/w));
zt3 = -zo3*w*sin(w*t)+vzo3*cos(w*t);
vxt3 = vxo3*cos(w*t)+(3*w*xo3+2*vyo3)*sin(w*t);
vyt3 = (6*w*xo3+4*vyo3)*cos(w*t)-2*vxo3*sin(w*t)-(6*w*xo3+3*vyo3);
vzt3 = -zo3*w*sin(w*t)+vzo3*cos(w*t);
% Grows infinitely from secular growth.
figure(13)
plot3(xt3,yt3,zt3)
title("CASE 3: 3D position from CW Equations for one day duration")

posMagListCW3 = zeros(1,length(t));
posVecVecCW3 = [xt3,yt3,zt3];
for j = 1:length(t)
    posMagListCW3(j) = norm(posVecVecCW3(j,:));
end
figure(14)
plot(t,posMagListCW3);
title("CASE 3: Distance Magnitude for CW Equations")
xlabel("Time (sec)")
ylabel("Distance (m)")
%% All "True" Orbits
figure(15)
plot3(RVtgt_prop(:,1),RVtgt_prop(:,2),RVtgt_prop(:,3),"r")
hold on
%figure(10)
plot3(RVint1_prop(:,1),RVint1_prop(:,2),RVint1_prop(:,3))
%figure(11)
plot3(RVint2_prop(:,1),RVint2_prop(:,2),RVint2_prop(:,3))
% figure(12)
plot3(RVint3_prop(:,1),RVint3_prop(:,2),RVint3_prop(:,3))
hold off
legend("Tgt","Chase-1","Chase-2","Chase-3")
%% Functions
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
%Keplerian Motion
function drdt=body2(t,r,mu)
drdt(1)=r(4);
drdt(2)=r(5);
drdt(3)=r(6);
drdt(4)=-1*mu*r(1)/(r(1)^2+r(2)^2+r(3)^2)^1.5;
drdt(5)=-1*mu*r(2)/(r(1)^2+r(2)^2+r(3)^2)^1.5;
drdt(6)=-1*mu*r(3)/(r(1)^2+r(2)^2+r(3)^2)^1.5;
drdt=drdt';
end
function E = kepler(M,e)
    E0=M;
    deltaE = 1;
    tol = 1e-4;
    n=0;
    while abs(deltaE) > tol
        deltaE = -(M - E0 + e*sin(E0))/(-1 + e*cos(E0));
        E0 = E0 + deltaE;
        n = n + 1; 
    end
    E = E0;
end
function nu = E2nu(E,e)
    %using rad now
    nu = 2*atan2(tan(E/2)*sqrt((1+e)/(1-e)),1);
    if nu < 0
        nu = 2*pi-nu;
    end
end