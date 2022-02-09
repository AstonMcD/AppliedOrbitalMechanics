function drdt=body2(t,r,mu)
drdt(1)=r(4);
drdt(2)=r(5);
drdt(3)=r(6);
drdt(4)=-1*mu*r(1)/(r(1)^2+r(2)^2+r(3)^2)^1.5;
drdt(5)=-1*mu*r(2)/(r(1)^2+r(2)^2+r(3)^2)^1.5;
drdt(6)=-1*mu*r(3)/(r(1)^2+r(2)^2+r(3)^2)^1.5;
drdt=drdt';
end