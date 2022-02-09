function[rv] = hw6oe2rv(oe, mu)

rMag = oe(1)*(1-oe(2)^2)/(1+oe(2)*cos(oe(6)));
pHat = [1 0 0];
qHat = [0 1 0];

rPQW = rMag*cos(oe(6))*pHat+rMag*sin(oe(6))*qHat;
RbigW = [cos(oe(5)) -1*sin(oe(5)) 0; sin(oe(5)) cos(oe(5)) 0; 0 0 1];
Ri = [1 0 0;0 cos(oe(3)) sin(-1*oe(3));0 sin(oe(3)) cos(oe(3))];
Rw = [cos(oe(4)) -1*sin(oe(4)) 0; sin(oe(4)) cos(oe(4)) 0; 0 0 1];

rIJK = RbigW*Ri*Rw*rPQW';
rv(1,1)=rIJK(1);
rv(2,1)=rIJK(2);
rv(3,1)=rIJK(3);

p = oe(1)*(1-oe(2)^2);
vPQW = -1*sin(oe(6))*sqrt(mu/p)*pHat+sqrt(mu/p)*(oe(2)+cos(oe(6)))*qHat;
vIJK=RbigW*Ri*Rw*vPQW';
rv(4,1)=vIJK(1);
rv(5,1)=vIJK(2);
rv(6,1)=vIJK(3);

end