function  E = M2E(M,e)
E = M;
if e<1
for i=1:8;
    ddf = (e*cos(E)-1);
    E  = E-(M-E + e*sin(E))/ddf;
end
else 
for i=1:8;
    ddf = (1 - e*cosh(E));
    E  = E-(M + E - e*sinh(E))/ddf;
end
end