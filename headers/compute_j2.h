// Compute J2 Function
// Triton Systems inc. By Colin Tierney
//
// The function "void computeJ2()" solves for the J2 component with the given parameters below
//
// The numbers used below in the main function are for testing the algorithm

#include <iostream>
#include <vector>

using namespace std;

void computeJ2(int k, vector<double>& xc, vector<double>& yc, vector<double>& zc, vector <double>& rsq,
    int mDeg, int pDeg, double& gx, double& gy, double& gz)
{
    //declare and intialize vectors
    vector<double> b(mDeg), d(mDeg), e(mDeg), zsq(mDeg), zsqe(mDeg);

    //declare g component variables
    double gx, gy, gz;

    //declare constants
    double const MU = .3986005e15; // Gravitational constant m ^ 3 kg ^ -1 s ^ -2
    double const M = 1; // Mass of earth in kg
    double const J2 = 1.755513178458136e+25; //J2 component

    //setting initial condition
    if (k == 2)
    {
        b[0] = pow(rsq[0], -1.5);

        //added for J2
        d[0] = pow(rsq[0], -2.5);
        e[0] = pow(rsq[0], -3.5);
        zsq[0] = pow(zc[0], 2);
        zsqe[0] = zsq[0] * e[0];

        //x component
        gx = -(MU * M * cauchyprod(xc, b, k - 1) + J2 * ((3 / 2) * cauchyprod(xc, d, k - 1) - (3 / 2) * 5 * cauchyprod(xc, zsqe, k - 1)));

        //y component
        gy = -(MU * M * cauchyprod(yc, b, k - 1) + J2 * ((3 / 2) * cauchyprod(yc, d, k - 1) - (3 / 2) * 5 * cauchyprod(yc, zsqe, k - 1)));

        //z component
        gz = -(MU * M * cauchyprod(zc, b, k - 1) + J2 * ((3 / 2) * 3 * cauchyprod(zc, d, k - 1) - (3 / 2) * 5 * cauchyprod(zc, zsqe, k - 1)));
    }
    else
    {
        //Roger's idea to zero out higher order terms and look at accuracy
        for (int i = pDeg; i <= mDeg; i++)
        {
            xc[i] = 0;
            yc[i] = 0;
            zc[i] = 0;
            rsq[i] = 0;
        }

        b[k - 1] = cauchypower(rsq, b, k - 2 - 1, -1.5);

        //added for J2
        d[k - 1] = cauchypower(rsq, d, k - 2 - 1, -2.5);
        e[k - 1] = cauchypower(rsq, e, k - 2 - 1, -3.5);
        zsq[k - 1] = cauchysquare(zc, k - 1 - 1);
        zsqe[k - 1] = cauchyprod(zsq, e, k - 1);

        //x component
        gx = -(MU * M * cauchyprod(xc, b, k - 1) + J2 * ((3 / 2) * cauchyprod(xc, d, k - 1) - (3 / 2) * 5 * cauchyprod(xc, zsqe, k - 1)));

        //y component
        gy = -(MU * M * cauchyprod(yc, b, k - 1) + J2 * ((3 / 2) * cauchyprod(yc, d, k - 1) - (3 / 2) * 5 * cauchyprod(yc, zsqe, k - 1)));

        //z component
        gz = -(MU * M * cauchyprod(zc, b, k - 1) + J2 * ((3 / 2) * 3 * cauchyprod(zc, d, k - 1) - (3 / 2) * 5 * cauchyprod(zc, zsqe, k - 1)));
    }

    //output values
    cout << "J2: \n" << endl;
    cout << "GX = " << gx << endl << "GY = " << gy << endl << "GZ = " << gz << endl;
}
