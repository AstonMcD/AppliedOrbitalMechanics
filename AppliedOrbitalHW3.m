%% Applied Orbital Mechanics HW#3

%Mo & to specified by initial conditions

%constants
mu = 3.986004415*10^14;
ae = 6378136.3;
we = 7.292115*10^-5;
g=9.81;
j2=1.082*10^-3;
%

%% PART B
%t = 0.0 sec
rv0=[0.447927332846917398E+07  0.245170029746151669E+07 -0.461804936123361625E+07 -0.455055498584698034E+04 -0.227543891441463575E+04 -0.564696384674148430E+04];
oe0 = hw6rv2oe(rv0,mu);


options = odeset('RelTol', 1e-14, 'AbsTol', 1e-14);
[t,rvf_prop]=ode45(@(t,r) body2(t,r,mu),[0:180:259200], rv0, options);

rvf_propagated = [rvf_prop(end,1), rvf_prop(end,2), rvf_prop(end,3), rvf_prop(end,4), rvf_prop(end,5), rvf_prop(end,6)];
oef_prop = hw6rv2oe(rvf_propagated,mu);

%t = 259200.0 sec
rvf =[0.187794820127492165E+07  0.835255100187982782E+06  0.653630524325667322E+07  0.645165445904742865E+04  0.338219169479819857E+04 -0.227861972090593599E+04];
oef = hw6rv2oe(rvf,mu);
%[a e i w bigOmega v(nu)]

rvDifference = rvf_propagated - rvf;
oeDifference = oef_prop - oef;

%% PART C
file = load('hw2-p01.txt');
timescale = (file(:,1));
rvTE = (file(:,2:7));


for b = 1:size(rvTE,1)
    range_prop = norm(rvf_prop(b,1:3));
    range_TE = norm(rvTE(b,1:3));
    rvSPEdiff_plot(b) = range_TE - range_prop;
end

figure(1)
plot(timescale,rvSPEdiff_plot)
xlabel('time (s)')
ylabel('Difference')
title('Difference between TE and IPE')
%The graph diverges, IPE value becomes less and less accurate as the orbit
%is propagated.

%% PART D
a0=oe0(1);
e0=oe0(2);
i0=oe0(3);
w0=oe0(4);
bigW0=oe0(5);

E0=2*atan(sqrt((1-oe0(2))/(1+oe0(2)))*tan(oe0(6)));
M0=E0-e0*sin(E0);
nu0= oe0(6); 

oe0_M=oe0;
oe0_M(6) = M0;
n = sqrt(mu/a0^3);

% T = 2*pi*sqrt(oe0(1)^3/mu);
% Tquarter = T/4;
% 
% E0=2*atan(sqrt((1-oe0(2))/(1+oe0(2)))*tan(oe0(6)));
% M0=sqrt(mu/(oe0(1)^3))*T-2*pi+E0-oe0(2)*sin(E0);

bigWbar = -1.5*n*(ae/a0)^2*j2/sqrt(1-e0^2)*cos(i0);
wbar = -0.75*n*(ae/a0)^2*j2/(1-e0^2)^2*(1-5*(cos(i0))^2);
Mbar = n*(1-(0.75*(ae/a0)^2*j2/(1-e0^2)^(3/2)*(1-3*(cos(i0))^2)));

bigWspe = bigW0 + bigWbar*(timescale);
wspe = w0 + wbar*(timescale);
Mspe = M0 + Mbar*(timescale);

spe=zeros(1441,6);
spe(:,1)=a0;
spe(:,2)=e0;
spe(:,3)=i0;
spe(:,4)=wspe;
spe(:,5)=bigWspe;
spe(:,6)=Mspe;

% Ef_M=zeros(1441,1)
% for iter=1:1441
%     Ef_M(i)=hw7kepler(Mspe(iter),e0);
% end
% 
% nu_M = 2*atan(sqrt((1+e0)/(1-e0))*tan(Ef_M/2));
% speNu=spe
oeTE=zeros(1441,6);
for i = 1:1441
    oeTE(i,:)= hw6rv2oe(rvTE(i,:), mu);
end
Ee=2*atan(sqrt((1-e0)/(1+e0)))*tan(oeTE(:,6)/2);
Me=Ee-e0*sin(Ee);
oeTE_M(:,6)=Me;

speFinal = [a0 e0 i0 wspe(end) bigWspe(end) Mspe(end)];
SPEdifference = speFinal - oe0_M'

diffwe=spe(:,4)-oeTE(:,4);
diffWe=spe(:,5)-oeTE(:,5);
diffMe=spe(:,6)-oeTE(:,6);
diffMe=diffMe+2*pi*(diffMe<-5);
% speInitial=zeros(1440,6);
% speInitial(:,1)=a0;
% speInitial(:,2)=e0;
% speInitial(:,3)=i0;
% speInitial(:,4)=wspe;
% speInitial(:,5)=bigWspe;
% speInitial(:,6)=Mspe;
% spe=[oe0,speInitial];

%% PROBLEM E

figure(2)
plot(timescale,diffwe)
title('Difference of Argument of Periapse for SPE and TE')
xlabel('time (s)')
ylabel('Difference (radians)')

figure(3)
plot(timescale,diffWe)
title('Difference of LAN for SPE and TE')
xlabel('time (s)')
ylabel('Difference (radians)')

figure(4)
plot(timescale,diffMe)
title('Difference of Mean Anomaly for SPE and TE')
xlabel('time (s)')
ylabel('Difference (radians)')
% rv_spe(1,:) = rv0;
% 
% w_spe = zeros(1441, 1);
% M_spe = zeros(1441, 1);
% W_spe = zeros(1441, 1);
% E_spe = zeros(1441, 1);
% nu_spe = zeros(1441, 1);
% for k = 2:size(timescale,1)
%     w_spe(k) = w0 + wbar*timescale(k);
%     M_spe(k) = M0 + Mbar*timescale(k);
%     W_spe(k) = bigW0 + bigWbar*timescale(k);
%     E_spe(k) = keplerEqn(M_spe(k),e0);
%     nu_spe(k) = 2*atan(sqrt((1+e0)/(1-e0))*tan(E_spe(k)/2));
%     oe_spe = [a0 e0 i0 w_spe(k) W_spe(k) nu_spe(k)];
%     rv_spe(k,:) = hw6oe2rv(oe_spe,mu)';
% end
% 
% for k = 1:size(rvTE,1)
%     range_spe = norm(rv_spe(k,1:3));
%     range_TE = norm(rvTE(k,1:3));
%     rvSPEdiff_plot(b) = range_TE - range_spe;
% end
% 
% figure(2)
% plot(timescale,rvSPEdiff_plot)
% xlabel('time (s)')
% ylabel('Difference')
% title('Difference between TE and SPE')

%% PROBLEM 2
muMars = 4.282837e13;
aMars = 3396190;
J2Mars = 1.932;
am2 = aMars + 100000;
I2 = linspace(0,pi,100)';
e2 = [0 .09 .148 .283 .331 .5];
 

nM = sqrt(muMars/(am2^3));
 

O2 = zeros(100,6);
w2 = zeros(100,6);
 

for p = 1:6
    O2(:,p) = (-3/2)*nM*((aMars/am2)^2)*J2Mars/sqrt(1-(e2(p)^2))*cos(I2);
    w2(:,p) = (-3/4)*nM*((aMars/am2)^2)*J2Mars/((1-(e2(p)^2))^2)*(1-5*(cos(I2).^2));
end
 

I2plot = I2*180/pi;
Omega2plot = O2*180/pi;
w2plot = w2*180/pi;
 

figure(5)
plot(I2plot, Omega2plot)
grid on
title('Daily Nodal Regression (Degrees/Day)')
xlabel('Inclination in degrees')
ylabel('Capital Omega dot in degrees/day')
 

figure(6)
plot(I2plot, w2plot)
grid on
title('Daily Apsidal Regression (Degrees/Day)')
xlabel('Inclination in degrees')
ylabel('Omega dot in degrees/day')