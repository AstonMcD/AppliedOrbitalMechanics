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