#define _USE_MATH_DEFINES
#define WITH_ALGEBRAICMATRIX

#include <dace/dace.h>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include "astroRoutines.h"

using namespace std;
using namespace DACE;

int main(void)
{
    int i, j, order, nvar, maxImpNum;
    const double mu = 3.986004418e14; //{[m^3/s^2]}
    const double rE = 6378137; //{m}

    // discretization, Hard Body Radius and time discretization
    double HBR, dt, t_integration, a, e, th0;

   // dummy variables for reading
    AlgebraicVector<double> xdum(6); 
    AlgebraicVector<double> Pdum(21); 
    AlgebraicMatrix<double> P0(6,6); 
    AlgebraicMatrix<double> P(6,6); 
    AlgebraicMatrix<double> Ppos(3,3); 
    AlgebraicMatrix<double> STM(6,6);

    // Relative position random variable and hard body radius
    ifstream infile;
	infile.open("initial_state.dat");
    if (infile.is_open())
	{
		for (int j = 0; j < 6; j++){
            infile >> xdum[j];}
		for (int j = 0; j < 21; j++){
            infile >> Pdum[j];}
        infile >> a; // semi-major axis
        infile >> e; // eccentricity
        infile >> th0; // initial true anomaly
        infile >> t_integration;
        infile >> dt;
        infile >> HBR; //Hard Body Radius
        infile >> order; //DA order
}
	else
	{
		cout << "initial_state.dat file not found!" << endl;
		return 1;
	}
	infile.close();

    
    P0.at(0,0) = Pdum[0];  P0.at(0,1) = Pdum[6];  P0.at(0,2) = Pdum[7];  P0.at(0,3) = Pdum[8];  P0.at(0,4) = Pdum[9];  P0.at(0,5) = Pdum[10];
    P0.at(1,0) = Pdum[6];  P0.at(1,1) = Pdum[1];  P0.at(1,2) = Pdum[11]; P0.at(1,3) = Pdum[12]; P0.at(1,4) = Pdum[13]; P0.at(1,5) = Pdum[14];
    P0.at(2,0) = Pdum[7];  P0.at(2,1) = Pdum[11]; P0.at(2,2) = Pdum[2];  P0.at(2,3) = Pdum[15]; P0.at(2,4) = Pdum[16]; P0.at(2,5) = Pdum[17];    
    P0.at(3,0) = Pdum[8];  P0.at(3,1) = Pdum[12]; P0.at(3,2) = Pdum[15]; P0.at(3,3) = Pdum[3];  P0.at(3,4) = Pdum[18]; P0.at(3,5) = Pdum[19];    
    P0.at(4,0) = Pdum[9];  P0.at(4,1) = Pdum[13]; P0.at(4,2) = Pdum[16]; P0.at(4,3) = Pdum[18]; P0.at(4,4) = Pdum[4];  P0.at(4,5) = Pdum[20];    
    P0.at(5,0) = Pdum[10]; P0.at(5,1) = Pdum[14]; P0.at(5,2) = Pdum[17]; P0.at(5,3) = Pdum[19]; P0.at(5,4) = Pdum[20]; P0.at(5,5) = Pdum[5];  

    Ppos.at(0,0) = Pdum[0];  Ppos.at(0,1) = Pdum[6];  Ppos.at(0,2) = Pdum[7];
    Ppos.at(1,0) = Pdum[6];  Ppos.at(1,1) = Pdum[1];  Ppos.at(1,2) = Pdum[11];
    Ppos.at(2,0) = Pdum[7];  Ppos.at(2,1) = Pdum[11]; Ppos.at(2,2) = Pdum[2];

    // DA initialisation
    nvar = 6;
    DA::init(order, nvar);
    DA::setEps(1e-30);
    
    int N  = floor(t_integration/cons(dt));

    // DA vector initialisation
    AlgebraicVector<DA> x0(6), x(6), xPos(3);

    for (int i = 0; i < 6; i++){
        x0[i] = xdum[i] + DA(i+1);}

    double t0 = 0;
    // write the output
    ofstream outfile;
    outfile.open("sqrMahaPoly.dat");
    outfile << setprecision(16);
    // N propagation
    for (int i = 0; i < N+1; i++)
    {
        if (e < 0.001){
            STM = CWStateTransition(a, i*dt, t0); // propagate covariance using state transition matrix
        }
        else{
            STM = YAStateTransition(a,e,th0,i*dt,t0);
        }
        x = STM*x0;
        P = STM*(P0*STM.transpose());

        Ppos.at(0,0) = P.at(0,0);  Ppos.at(0,1) = P.at(0,1);  Ppos.at(0,2) = P.at(0,2);
        Ppos.at(1,0) = P.at(1,0);  Ppos.at(1,1) = P.at(1,1);  Ppos.at(1,2) = P.at(1,2);
        Ppos.at(2,0) = P.at(2,0);  Ppos.at(2,1) = P.at(2,1);  Ppos.at(2,2) = P.at(2,2);
            
        xPos[0] = x[0];
        xPos[1] = x[1];
        xPos[2] = x[2];

        DA sqrMaha = SqrMaha(Ppos,xPos);
        outfile << sqrMaha << endl;
        }
    outfile.close();
}