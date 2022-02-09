function[oe] = hw6rv2oe(rv,mu)
r = [rv(1) rv(2) rv(3)];
rMag = norm(r);

v = [rv(4) rv(5) rv(6)];
vMag = norm(v);

oe(1,1) = -1*mu / (2*(vMag^2/2-mu/rMag));

h = cross(r,v);
e = cross(v,h)/mu-(r/rMag);
oe(2,1) = norm(e);

zHat = [0 0 1];

oe(3,1) = acos(dot(zHat,h/norm(h)));
nHat = cross(zHat,h) / norm(cross(zHat,h));
oe(4,1) = acos(dot(nHat,e)/norm(e));
if e(3)<0
    oe(4,1) = -1*oe(4,1);
end

oe(5,1) = atan2(nHat(2), nHat(1));

oe(6,1) = acos(dot(r,e)/(rMag*norm(e)));
if dot(r,v)<0
    oe(6,1) = -1*oe(6,1);
end
end

