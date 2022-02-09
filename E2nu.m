function nu = E2nu(E,e)
    %using rad now
    nu = 2*atan2(tan(E/2)*sqrt((1+e)/(1-e)),1);
    if nu < 0
        nu = 2*pi-nu;
    end
end
    