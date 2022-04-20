// Compute Three Degree of Freedom Trajectory Model for a Debris Object within Earth's Spherical Gravity
// Colin Tierney
//
// This function solves the 3DoF equations of flight for a debris object
// in rectangular coordinates(x, y, z) with z going through the North Pole
//
// The ODEs
// x' = vx;
// y' = vy;
// z' = vz;
//
// x'' = vx' = -(G*M* z/r^3)
// y'' = vy' = -(G*M* z/r^3)
// z'' = vz' = -(G*M* z/r^3)
//
// r = sqrt(x ^ 2 + y ^ 2 + z ^ 2), rsq = r ^ 2


//libaries
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <time.h>

//header files
#include "cauchy_power.h"
#include "cauchy_prod.h"
#include "cauchy_square.h"
#include "horners.h"
#include "compute_spherical.h"
#include "compute_PSM_Method.h"
#include "outputData.h"

using namespace std;

//Function Declarations
void PSM_3DoF(double ipx, double ivx, double ipy, double ivy, double ipz, double ivz, double iax, double iay,
    double iaz, int mDeg, double dt, int endTime, double& gx, double& gy, double& gz, int k);

int main()
{
    //to compute runtime
    clock_t startTime = clock();

    //max degree to use for PSM Order; "PSM4" ==> mDeg = 4
    int mDeg = 5;

    //to zero out poly vectors after a certain point
    int pDeg = mDeg;

    //end time
    int endTime = 510;

    //time step
    double dt = 1;

    //declare gravity variables
    double gx, gy, gz;

    //initial conditions; these are from a test using TGx
    double ipx = 356403.951;        //x position
    double ipy = 4076924.673;       //y position
    double ipz = 5034509.72;        //z position
    double ivx = 3068.498667;       //x velocity
    double ivy = 1991.659005;       //y velocity
    double ivz = 582.244318;        //z velocity
    double iax = 1;                 //x acceleration
    double iay = 1;                 //y acceleration
    double iaz = 1;                 //z acceleration

    //cauchy degree
    int k = mDeg + 1;

    //power for cauchypower function
    int p = 2;

    //function calls
    PSM_3DoF(ipx, ivx, ipy, ivy, ipz, ivz, iax, iay, iaz, mDeg, dt, endTime, gx, gy, gz, k);

    //output runtime in the console window
    cout << endl << double(clock() - startTime) / (double)CLOCKS_PER_SEC << " seconds." << endl;
}

void PSM_3DoF(double ipx, double ivx, double ipy, double ivy, double ipz, double ivz, double iax,
    double iay, double iaz, int mDeg, double dt, int endTime, double& gx,
    double& gy, double& gz, int k)
{
    //variable declarations
    //psm vectors
    vector <double> xc(mDeg + 2);       //for x coordinate
    vector <double> yc(mDeg + 2);       //for y coordinate
    vector <double> zc(mDeg + 2);       //for z coordinate
    vector <double> vxc(mDeg + 2);      //for vx coordinate
    vector <double> vyc(mDeg + 2);      //for vy coordinate
    vector <double> vzc(mDeg + 2);      //for vz coordinate
    vector <double> rsq(mDeg + 2);      //for rsq = x^2 + y^2 + z^2 coordinate

    //initialize PSM Vectors
    vector <double> xPSM(endTime * 10);
    vector <double> yPSM(endTime * 10);
    vector <double> zPSM(endTime * 10);
    vector <double> vxPSM(endTime * 10);
    vector <double> vyPSM(endTime * 10);
    vector <double> vzPSM(endTime * 10);

    //initialize variables
    //PSM vectors will contain the values of the MacLauren Polynomials
    xPSM[1] = ipx;
    yPSM[1] = ipy;
    zPSM[1] = ipz;
    vxPSM[1] = ivx;
    vyPSM[1] = ivy;
    vzPSM[1] = ivz;

    //initialize the PSM variables for the first time step
    xc[1] = ipx;
    vxc[1] = ivx;
    yc[1] = ipy;
    vyc[1] = ivy;
    zc[1] = ipz;
    vzc[1] = ivz;
    rsq[1] = pow(ipx, 2) + pow(ipy, 2) + pow(ipz, 2);

    //csv file formatting
    ofstream file;
    file.open("output.csv");

    //define variables for computeMethodPSM
    int n = 1;                              //counter variable 
    int method = 1;                         //FSFO = 1; ASFO = 2; ASAO = 3
    int choice = 1;                         //corresponds to outputData.h function
    double tol = 1;                         //tolerance; the amount of precision to solve for
    int maxStep = 10;                        //the max timestep to go up to in the solver
    double& h = dt;                         //new variable for the timestep
    vector <int> degM(endTime * 10);        //max degree vector for storing the counter variable "ns"
    vector <double> timeVec(endTime * 10);  //time vector

    //intermediate variables for storing timestep value and psm order used at each step
    vector <double> timeStep(endTime * 10);
    vector <double> psmOrder(endTime * 10);

    while (timeVec[n] < endTime)
    {
        //PSM Order loop
        for (int k = 2; k <= mDeg + 1; k++)
        {
            //call gravity function
            computeSpherical(k, xc, yc, zc, rsq, mDeg, mDeg, gx, gy, gz);

            //main ODE's
            xc[k] = vxc[k - 1] / (k - 1);
            yc[k] = vyc[k - 1] / (k - 1);
            zc[k] = vzc[k - 1] / (k - 1);
            vxc[k] = gx / (k - 1);
            vyc[k] = gy / (k - 1);
            vzc[k] = gz / (k - 1);

            //auxilliary variables
            rsq[k] = cauchysquare(xc, k - 1) + cauchysquare(yc, k - 1) + cauchysquare(zc, k - 1);

            //if using 'ASAO', call funtion inside the PSM Order for loop
            if (method == 3)
            {
                //call computeMethodPSM
                computeMethodPSM(method, xc, yc, zc, vxc, vyc, vzc, k, mDeg, tol, maxStep, h, n, degM);

                //store each timestep
                timeStep[n] = dt;
                //intermediate variable to store the change in PSM Order
                psmOrder[n] = degM[n];
            }
        }

        //call method function if 'FSFO' or 'ASFO'
        if (method == 1 || method == 2)
        {
            //call computeMethodPSM
            computeMethodPSM(method, xc, yc, zc, vxc, vyc, vzc, k, mDeg, tol, maxStep, h, n, degM);
        }

        //output data
        if (choice == 1)
            outputData(file, method, xPSM, yPSM, zPSM, vxPSM, vyPSM, vzPSM, xc, yc, zc, vxc, vyc, vzc,
                psmOrder, timeStep, n, endTime, timeVec, mDeg, degM, h, dt, choice);

        //x coordinate coefficients
        xPSM[n + 1] = horners(h, xc, degM[n]);
        //update xPSM vector
        xc[1] = xPSM[n + 1];
        //y coordinate coefficients
        yPSM[n + 1] = horners(h, yc, degM[n]);
        //update yPSM vector
        yc[1] = yPSM[n + 1];
        //z coordinate coefficients
        zPSM[n + 1] = horners(h, zc, degM[n]);
        //update zPSM vector
        zc[1] = zPSM[n + 1];

        //output data
        if (choice == 2)
            outputData(file, method, xPSM, yPSM, zPSM, vxPSM, vyPSM, vzPSM, xc, yc, zc, vxc, vyc, vzc,
                psmOrder, timeStep, n, endTime, timeVec, mDeg, degM, h, dt, choice);

        //velocity coordinates; same method for solving as the xc,yc,zc from above
        vxPSM[n + 1] = horners(h, vxc, degM[n]);
        vxc[1] = vxPSM[n + 1];
        vyPSM[n + 1] = horners(h, vyc, degM[n]);
        vyc[1] = vyPSM[n + 1];
        vzPSM[n + 1] = horners(h, vzc, degM[n]);
        vzc[1] = vzPSM[n + 1];

        //auxilliary variable
        rsq[1] = pow(xc[1], 2) + pow(yc[1], 2) + pow(zc[1], 2);

        //update time vector
        timeVec[n + 1] = timeVec[n] + h;

        n += 1;
    }
    file.close();
}