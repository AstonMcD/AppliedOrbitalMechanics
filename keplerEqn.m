function E = keplerEqn(M,e)
n=1;
eps = 1e-14;

if -pi < M < 0 | M > pi
      E(n) = M - e;
else
        E(n) = M + e;
end

E(n+1) = E(n) + ((M - E(n) +e*sin(E(n))) / (1-e*cos(E(n))));

while abs(E(n+1) - E(n)) > eps
    n = n+1;
    E(n+1) = E(n) + ((M - E(n) + e*sin(E(n)) / (1 - e*cos(E(n)))));
end

E= E(end);
end