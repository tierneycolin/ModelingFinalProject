// Compute Spherical Function
// Triton Systems inc. By Colin Tierney
//
// The function "void computeSpherical()" solves for the spherical gravity component with the below given parameters
//
// The numbers used below in the main function are for testing the algorithm

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

void computeSpherical(int k, vector<double>& xc, vector<double>& yc, vector<double>& zc, vector <double>& rsq, 
    int mDeg, int pDeg, double& gx, double& gy, double& gz)
{
    static vector<double> b(mDeg + 1);
    double const G = .3986005e15; // Gravitational constant m ^ 3 kg ^ -1 s ^ -2
    double const M = 1; 

    // setting initial condition
    if (k == 2)
    {
        b[1] = pow(rsq[1], -1.5);

        // x component
        gx = -(G * M * cauchyprod(xc, b, k - 1));

        // y component
        gy = -(G * M * cauchyprod(yc, b, k - 1));

        // z component
        gz = -(G * M * cauchyprod(zc, b, k - 1));
    }
    else
    {
        //Roger's idea to zero out higher order terms and look at accuracy
        for (int i = pDeg + 1; i <= mDeg; i++)
        {
            xc[i] = 0;
            yc[i] = 0;
            zc[i] = 0;
            rsq[i] = 0;
        }

        b[k - 1] = cauchypower(rsq, b, k - 2 - 1, -1.5);

        // x component
        gx = -(G * M * cauchyprod(xc, b, k - 1));

        // y component
        gy = -(G * M * cauchyprod(yc, b, k - 1));

        // z component
        gz = -(G * M * cauchyprod(zc, b, k - 1));
    }
}