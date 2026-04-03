function [z1, z2, z3] = cubic_solver(a0, a1, a2)
    q = a1/3 - a2^2/9;
    r = (a1*a2 - 3*a0)/6 - a2^3/27;
    
    % case 1:
    if (r^2 + q^3 > 0)
        A = (abs(r) + sqrt(r^2 + q^3))^(1/3);
        if (r >= 0)
            t1 = A - q/A;
        else
            t1 = q/A - A;
        end
        z1 = t1 - a2/3;
        x2 = -t1/2 - a2/3;
        y2 = sqrt(3)/2 * (A + q/A);
        z2 = x2 + 1i*y2;
        z3 = x2 - 1i*y2;
    else
        if (q == 0)
            teta = 0;
        else
            teta = acos(r/sqrt(-q^3));
        end

        phi1 = teta/3;
        phi2 = phi1 - 2*pi/3;
        phi3 = phi1 + 2*pi/3;

        z1 = 2*sqrt(-q)*cos(phi1) - a2/3;
        z2 = 2*sqrt(-q)*cos(phi2) - a2/3;
        z3 = 2*sqrt(-q)*cos(phi3) - a2/3;
    end
end