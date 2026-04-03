
from cmath import *

def quartic_solver(A0,A1,A2,A3):

    def cubic_solver(a0,a1,a2):
        q=a1/3-a2**2/9
        r=(a1*a2-3*a0)/6-a2**3/27
        # case 1:
        if (r**2+q**3>0):
            A=(abs(r)+sqrt(r**2+q**3))**(1/3)
            if (r>=0):
                t1=A-q/A
            else:
                t1=q/A-A
            z1=t1-a2/3
            x2=-t1/2-a2/3

            y2=sqrt(3)/2*(A+q/A)

            z2=x2+sqrt(-1)*y2
            z3=x2-sqrt(-1)*y2
        else:
            if (q==0):
                teta=0
            else:
                teta=acos(r/sqrt(-q**3))
    
            phi1=teta/3
            phi2=phi1-2*pi/3
            phi3=phi1+2*pi/3

            z1=2*sqrt(-q)*cos(phi1)-a2/3

            z2=2*sqrt(-q)*cos(phi2)-a2/3

            z3=2*sqrt(-q)*cos(phi3)-a2/3
        return z1, z2, z3


    C=A3/4
    b2=A2-6*C**2
    b1=A1-2*A2*C+8*C**3
    b0=A0-A1*C+A2*C**2-3*C**4

    a2=b2
    a1=(b2**2/4-b0)
    a0=-b1**2/8

    m0, m1, m2=cubic_solver(a0,a1,a2);

    SIGMA=-1+2*(b1>0);
    R=SIGMA*sqrt(m0**2+b2*m0+b2**2/4-b0);
    c11=sqrt(m0/2)-C;
    c12=sqrt(-m0/2-b2/2-R);
    Z1=c11+c12;
    Z2=c11-c12;
    c21=-sqrt(m0/2)-C;
    c22=sqrt(-m0/2-b2/2+R);
    Z3=c21+c22;
    Z4=c21-c22;

    return Z1, Z2, Z3, Z4
