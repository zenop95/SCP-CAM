function  theta = M2theta(M,e)
E = M;
if e<1
    while abs((M-E + e*sin(E)))>1e-13;
        ddf = (e*cos(E)-1);
        E  = E-(M-E + e*sin(E))/ddf;
    end
    theta = 2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));    
else
    for i=1:20;
        ddf = (1 - e*cosh(E));
        E  = E-(M + E - e*sinh(E))/ddf;
        theta  = 2*atan(sqrt((1+e)/(e-1))*tanh(E/2));
    end
end