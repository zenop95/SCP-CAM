function [Z1, Z2, Z3, Z4] = quartic_solver(A0, A1, A2, A3)

    C  = A3/4;
    b2 = A2 - 6*C^2;
    b1 = A1 - 2*A2*C + 8*C^3;
    b0 = A0 - A1*C + A2*C^2 - 3*C^4;

    a2 = b2;
    a1 = (b2^2/4 - b0);
    a0 = -b1^2/8;

    m0 = cubic_solver(a0, a1, a2);

    SIGMA = -1 + 2*(b1 > 0);
    R     = SIGMA * sqrt(m0^2 + b2*m0 + b2^2/4 - b0);
    c11   = sqrt(m0/2) - C;
    c12   = sqrt(-m0/2 - b2/2 - R);
    Z1    = c11 + c12;
    Z2    = c11 - c12;
    c21   = -sqrt(m0/2) - C;
    c22   = sqrt(-m0/2 - b2/2 + R);
    Z3    = c21 + c22;
    Z4    = c21 - c22;
end