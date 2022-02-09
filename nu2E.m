function E = nu2E(nu,e)
%using rad now
    E = 2*atan2(tan(nu/2)*sqrt((1-e)/(1+e)),1);
    if E < 0
        E = 2*pi-E;
    end
end