// Compute J4 Function
// Triton Systems inc. By Colin Tierney
//
// The function "void computeJ4()" solves for the J4 component with the given parameters below
//
// The numbers used below in the main function are for testing the algorithm

#include <iostream>
#include <vector>
#include <cmath>

//header files
#include "cauchy_prod.h"
#include "cauchy_power.h"
#include "cauchy_square.h"

using namespace std;

void computeJ4(int k, vector<double>& xc, vector<double>& yc, vector<double>& zc, vector <double>& rsq, int mDeg);

int main()
{
    //declare/intitialize test variables
    //use up to c_8 coefficients
    int mDeg = 8;

    //coefficient vector for x
    vector <double> xc(mDeg + 1);
    xc[1] = 6.4e6;

    //coefficient vector for y
    vector <double> yc(mDeg + 1);
    yc[1] = 1000;

    //coefficient vector for z
    vector <double> zc(mDeg + 1);
    zc[1] = 2000;

    //coefficient vector for radius^2
    vector <double> rsq(mDeg + 1);
    rsq[1] = pow(xc[1], 2) + pow(yc[1], 2) + pow(zc[1], 2);

    //degree
    int k = 2;

    //call computeJ2 funciton which contains other functions inside of it
    computeJ4(k, xc, yc, zc, rsq, mDeg);
}

void computeJ4(int k, vector<double>& xc, vector<double>& yc, vector<double>& zc, vector <double>& rsq, int mDeg)
{
    //declare constants
    const double G = 6.67408 * pow(10, -11);         // Gravitational constant m ^ 3 kg ^ -1 s ^ -2
    const double M = 5.972 * pow(10, 24);            // Mass of earth in kg
    const double J2 = 1.755513178458136e+25;
    const double J3 = -2.619361570252708e+29;
    const double J4 = -1.068381886873931e+36;

    //gravity coefficients
    double gx, gy, gz;

    //initiliaze vectors
    vector <double> b(mDeg);
    vector <double> d(mDeg);
    vector <double> e(mDeg);
    vector <double> zsq(mDeg);
    vector <double> zsqe(mDeg);
    vector <double> f(mDeg);
    vector <double> ze(mDeg);
    vector <double> zcub(mDeg);
    vector <double> zcubf(mDeg);
    vector <double> z4(mDeg);
    vector <double> z4f(mDeg);
    vector <double> g(mDeg);
    vector <double> z4g(mDeg);
    vector <double> zsqf(mDeg);
    vector <double> z5(mDeg);
    vector <double> z5g(mDeg);

    //setting initial condition
    if (k == 2)
    {
        b[1] = pow(rsq[1], (-1.5));

        // added for J2
        d[1] = pow(rsq[1], (-2.5));
        e[1] = pow(rsq[1], (-3.5));
        zsq[1] = pow(zc[1], 2);
        zsqe[1] = zsq[1] * e[1];

        // added for J3
        f[1] = pow(rsq[1], -4.5);
        ze[1] = zc[1] * pow(rsq[1], -3.5);
        zcub[1] = pow(zc[1], 3);
        zcubf[1] = zcub[1] * f[1];
        z4[1] = pow(zc[1], 4);
        z4f[1] = z4[1] * f[1];

        // added for J4
        g[1] = pow(rsq[1], (-5.5));
        z4g[1] = z4[1] * g[1];
        zsqf[1] = zsq[1] * f[1];
        z5[1] = pow(zc[1], 5);
        z5g[1] = z5[1] * g[1];

        // x component
        gx = -(G * M * cauchyprod(xc, b, k - 1) +  
        J2 * ((3 / 2) * cauchyprod(xc, d, k - 1) -        // J2
        (3 / 2) * 5 * cauchyprod(xc, zsqe, k - 1)) +      // J2
        J3 * ((35 / 2) * cauchyprod(xc, zcubf, k - 1) +        // J3
        - (15 / 2) * cauchyprod(xc, ze, k - 1)) +              // J3
        J4 * ((315 / 8) * cauchyprod(xc, z4g, k - 1) +    // J4
        - (105 / 4) * cauchyprod(xc, zsqf, k - 1) +       // J4
        (15 / 8) * cauchyprod(xc, e, k - 1)));            // J4

        // y component
        gy = -(G * M * cauchyprod(yc, b, k - 1) +  
        J2 * ((3 / 2) * cauchyprod(yc, d, k - 1) -        // J2
        (3 / 2) * 5 * cauchyprod(yc, zsqe, k - 1)) +      // J2
        J3 * ((35 / 2) * cauchyprod(yc, zcubf, k - 1) +         // J3
        - (15 / 2) * cauchyprod(yc, ze, k - 1)) +               // J3
        J4 * ((315 / 8) * cauchyprod(yc, z4g, k - 1) +    // J4
        - (105 / 4) * cauchyprod(yc, zsqf, k - 1) +       // J4
        (15 / 8) * cauchyprod(yc, e, k - 1)));            // J4

        // z component
        gz = -(G * M * cauchyprod(zc, b, k - 1) +
        J2 * ((3 / 2) * 3 * cauchyprod(zc, d, k - 1) -    // J2
        (3 / 2) * 5 * cauchyprod(zc, zsqe, k - 1))        // J2
        - J3 * ((35 / 2) * z4f[k - 1] +                         // J3
        15 * zsqe[k - 1] - (3 / 2) * d[k - 1]) +                // J3
        -J4 * ((315 / 8) * z5g[k - 1] +                   // J4
        (175 / 4) * zcubf[k - 1] +                        // J4
        -(75 / 8) * ze[k - 1]));                          // J4

        //test
        std::cout << endl << sqrt(pow(gx, 2) + pow(gy, 2) + pow(gz, 2)) << endl << gx << endl << gy << endl << gz << endl;
    }
    
    
    else
    {
        b[k - 1] = cauchypower(rsq, b, k - 2 - 1, -1.5);

        // added for J2
        d[k - 1] = cauchypower(rsq, d, k - 2 - 1, -2.5);
        e[k - 1] = cauchypower(rsq, e, k - 2 - 1, -3.5);
        zsq[k - 1] = cauchysquare(zc, k - 1 - 1);
        zsqe[k - 1] = cauchyprod(zsq, e, k - 1);

        // added for J3
        f[k - 1] = cauchypower(rsq, f, k - 2 - 1, -4.5);
        ze[k - 1] = cauchyprod(zc, e, k - 1);
        zcub[k - 1] = cauchypower(zc, zcub, k - 2 - 1, 3);
        zcubf[k - 1] = cauchyprod(zcub, f, k - 1);
        z4[k - 1] = cauchypower(zc, z4, k - 2 - 1, 4);
        z4f[k - 1] = cauchyprod(z4, f, k - 1);

        // added for  J4
        g[k - 1] = cauchypower(rsq, g, k - 2 - 1, -5.5);
        z4g[k - 1] = cauchyprod(z4, g, k - 1);
        zsqf[k - 1] = cauchyprod(zsq, f, k - 1);
        z5[k - 1] = cauchypower(zc, z5, k - 2 - 1, 5);
        z5g[k - 1] = cauchyprod(z5, g, k - 1);

        // x component
        gx = -(G * M * cauchyprod(xc, b, k - 1) +  
        J2 * ((3 / 2) * cauchyprod(xc, d, k - 1) -         // J2
        (3 / 2) * 5 * cauchyprod(xc, zsqe, k - 1)) +       // J2
        J3 * ((35 / 2) * cauchyprod(xc, zcubf, k - 1) +         // J3
        - (15 / 2) * cauchyprod(xc, ze, k - 1)) +               // J3
        J4 * ((315 / 8) * cauchyprod(xc, z4g, k - 1) +     // J4
        - (105 / 4) * cauchyprod(xc, zsqf, k - 1) +        // J4
        (15 / 8) * cauchyprod(xc, e, k - 1)));             // J4

        // y component
        gy = -(G * M * cauchyprod(yc, b, k - 1) +  
        J2 * ((3 / 2) * cauchyprod(yc, d, k - 1) -         // J2
        (3 / 2) * 5 * cauchyprod(yc, zsqe, k - 1)) +       // J2
        J3 * ((35 / 2) * cauchyprod(yc, zcubf, k - 1) +         // J3
        - (15 / 2) * cauchyprod(yc, ze, k - 1)) +               // J3
        J4 * ((315 / 8) * cauchyprod(yc, z4g, k - 1) +     // J4
        - (105 / 4) * cauchyprod(yc, zsqf, k - 1) +        // J4
        (15 / 8) * cauchyprod(yc, e, k - 1)));             // J4

        // z component
        gz = -(G * M * cauchyprod(zc, b, k - 1) + 
        J2 * ((3 / 2) * 3 * cauchyprod(zc, d, k - 1) -     // J2
        (3 / 2) * 5 * cauchyprod(zc, zsqe, k - 1))         // J2
        - J3 * ((35 / 2) * z4f[k - 1] +                         // J3
        15 * zsqe[k - 1] - (3 / 2) * d[k - 1]) +                // J3
        - J4 * ((315 / 8) * z5g[k - 1] +                   // J4
        (175 / 4) * zcubf[k - 1] +                         // J4
        - (75 / 8) * ze[k - 1]));                          //J4
    }
    std::cout << endl << sqrt(pow(gx, 2) + pow (gy, 2) + pow(gz, 2)) << endl << gx << endl << gy << endl << gz << endl;

}