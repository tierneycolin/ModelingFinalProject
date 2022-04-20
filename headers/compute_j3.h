// Compute J3 Function
// Triton Systems inc. By Colin Tierney
//
// The function "void computeJ3()" solves for the J3 component with the given parameters below
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

void computeJ3(int k, vector<double>& xc, vector<double>& yc, vector<double>& zc, vector <double>& rsq, int mDeg, int pDeg);

int main()
{
    //declare/intitialize test variables
    //use up to c_8 coefficients
    int mDeg = 8;

    //pDeg to zero out vectors after a certain point
    int pDeg = mDeg;

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
    computeJ3(k, xc, yc, zc, rsq, mDeg, pDeg);
}

void computeJ3(int k, vector<double>& xc, vector<double>& yc, vector<double>& zc, vector <double>& rsq, int mDeg, int pDeg)
{
    //declare variables
    double G = 6.67408 * pow(10,-11);               // Gravitational constant m ^ 3 kg ^ -1 s ^ -2
    double M = 5.972 * pow(10, 24);                  // Mass of earth in kg
    double J2 = 1.755513178458136e25;
    double C12 = 3.914137379344170e18;
    double S12 = -2.502032058856706e19;
    double C22 = -2.55288078060225e22;
    double S22 = 1.46545661274048e22;
    double J3 = -2.619361570252708e29;

    //declare g component variables
    double gx, gy, gz;

    // initialize vectors
    static vector <double> b(mDeg);        //for rsq^ { -1.5 } coordinate

    //added for J2
    static vector <double> d(mDeg);         // for rsq^ { -2.5 } coordinate
    static vector <double> e(mDeg);         // for rsq^ { -3.5 } coordinate
    static vector <double> zsq(mDeg);       // for z ^ 2
    static vector <double> zsqe(mDeg);      // for z ^ 2 * e

    // added for 12 term
    static vector <double> xsq(mDeg);
    static vector <double> ysq(mDeg);
    static vector <double> xy(mDeg);
    static vector <double> xyz(mDeg);
    static vector <double> xsqz(mDeg);
    static vector <double> ysqz(mDeg);

    // added for 22 term
    static vector <double> xysq(mDeg);
    static vector <double> xsqy(mDeg);
    static vector <double> ycub(mDeg);
    static vector <double> xcub(mDeg);

    // added for J3
    static vector <double> f(mDeg);            // for rsq^ { -4.5 }
    static vector <double> ze(mDeg);           // for z* rsq^ { -3.5 }
    static vector <double> zcub(mDeg);         // for z ^ 3
    static vector <double> zcubf(mDeg);        // for z ^ 3 * rsq^ { -4.5 }
    static vector <double> z4(mDeg);           // for z ^ 4;
    static vector <double> z4f(mDeg);          // for z ^ 4 * rsq^ { -4.5 }

    //setting intial condition
    if (k == 2)
    {
        b[1] = pow(rsq[1], -1.5);

        // added for J2
        d[1] = pow(rsq[1], -2.5);
        e[1] = pow(rsq[1], -3.5);
        zsq[1] = pow(zc[1], 2);
        zsqe[1] = zsq[1] * e[1];

        // added for 12 term
        xsq[1] = pow(xc[1], 2);
        ysq[1] = pow(yc[1], 2);
        xy[1] = xc[1] * yc[1];
        xyz[1] = xy[1] * zc[1];
        xsqz[1] = xsq[1] * zc[1];
        ysqz[1] = ysq[1] * zc[1];

        // added for 22 term
        xysq[1] = xy[1] * yc[1];
        xsqy[1] = xy[1] * xc[1];
        ycub[1] = ysq[1] * yc[1];
        xcub[1] = xsq[1] * xc[1];

        // added for J3
        f[1] = pow(rsq[1], -4.5);
        ze[1] = zc[1] * pow(rsq[1], -3.5);
        zcub[1] = pow(zc[1], 3);
        zcubf[1] = zcub[1] * f[1];
        z4[1] = pow(zc[1], 4);
        z4f[1] = z4[1] * f[1];
        
        // x component
        gx = -(G * M * cauchyprod(xc, b, k - 1) +
        J2 * ((3 / 2) * cauchyprod(xc, d, k - 1) -                  // J2
        (3 / 2) * 5 * cauchyprod(xc, zsqe, k - 1)) +                // J2
        +3 * C12 * (cauchyprod(zc, d, k - 1) - 5 * cauchyprod(xsqz, e, k - 1)) -        // 21
        15 * S12 * cauchyprod(xyz, e, k - 1)                                            // 21
        + 3 * C22 * (2 * cauchyprod(xc, d, k - 1) - 5 * cauchyprod(xcub, e, k - 1) +    // 22
        5 * cauchyprod(xysq, e, k - 1)) + 6 * S22 * (cauchyprod(yc, d, k - 1)           // 22
        - 5 * cauchyprod(xsqy, e, k - 1)) +                 // 22
        J3 * ((35 / 2) * cauchyprod(xc, zcubf, k - 1) +     // J3
        -(15 / 2) * cauchyprod(xc, ze, k - 1)));            // J3

        // y component
        gy = -(G * M * cauchyprod(yc, b, k - 1) +
        J2 * ((3 / 2) * cauchyprod(yc, d, k - 1) -               // J2
        (3 / 2) * 5 * cauchyprod(yc, zsqe, k - 1)) +             // J2
        3 * S12 * (cauchyprod(zc, d, k - 1) - 5 * cauchyprod(ysqz, e, k - 1)) -         // 21
        15 * C12 * cauchyprod(xyz, e, k - 1)                                            // 21
        + 3 * C22 * (-2 * cauchyprod(yc, d, k - 1) + 5 * cauchyprod(ycub, e, k - 1) -   // 22
        5 * cauchyprod(xsqy, e, k - 1)) + 6 * S22 * (cauchyprod(xc, d, k - 1)           // 22
        - 5 * cauchyprod(xysq, e, k - 1)) +               // 22
        J3 * ((35 / 2) * cauchyprod(yc, zcubf, k - 1) +   // J3
        -(15 / 2) * cauchyprod(yc, ze, k - 1)));          // J3


        // z component
        gz = -(G * M * cauchyprod(zc, b, k - 1) +
        J2 * ((3 / 2) * 3 * cauchyprod(zc, d, k - 1) -       // J2
        (3 / 2) * 5 * cauchyprod(zc, zsqe, k - 1))           // J2
        + 3 * C12 * (cauchyprod(xc, d, k - 1) - 5 * cauchyprod(xc, zsqe, k - 1))     // 21
        + 3 * S12 * (cauchyprod(yc, d, k - 1) - 5 * cauchyprod(yc, zsqe, k - 1))     // 21
        + 15 * C22 * (cauchyprod(xsqz, e, k - 1) - cauchyprod(ysqz, e, k - 1))       // 22
        + 30 * S22 * cauchyprod(xyz, e, k - 1)            // 22
        + J3 * ((-35 / 2) * z4f[k - 1] +                  // J3
        15 * zsqe[k - 1] - (3 / 2) * d[k - 1]));          // J3
                
        //output gravity coefficients
        cout << gx << endl << gy << endl << gz << endl;
        cout << sqrt(pow(gx, 2) + pow(gy, 2) + pow(gz, 2));
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

        // added for J2
        d[k - 1] = cauchypower(rsq, d, k - 2 - 1, -2.5);
        e[k - 1] = cauchypower(rsq, e, k - 2 - 1, -3.5);
        zsq[k - 1] = cauchysquare(zc, k - 1 - 1);
        zsqe[k - 1] = cauchyprod(zsq, e, k - 1);

        // added for C12
        xsq[k - 1] = cauchysquare(xc, k - 1 - 1);
        ysq[k - 1] = cauchysquare(yc, k - 1 - 1);
        xy[k - 1] = cauchyprod(xc, yc, k - 1);
        xyz[k - 1] = cauchyprod(xy, zc, k - 1);
        xsqz[k - 1] = cauchyprod(xsq, zc, k - 1);
        ysqz[k - 1] = cauchyprod(ysq, zc, k - 1);

        // added for C22, S22
        xysq[k - 1] = cauchyprod(xy, yc, k - 1);
        xsqy[k - 1] = cauchyprod(xy, xc, k - 1);
        ycub[k - 1] = cauchypower(yc, ycub, k - 2 - 1, 3);
        xcub[k - 1] = cauchypower(xc, xcub, k - 2 - 1, 3);

        // added for J3
        f[k - 1] = cauchypower(rsq, f, k - 2 - 1, -4.5);
        ze[k - 1] = cauchyprod(zc, e, k - 1);
        zcub[k - 1] = cauchypower(zc, zcub, k - 2 - 1, 3);
        zcubf[k - 1] = cauchyprod(zcub, f, k - 1);
        z4[k - 1] = cauchypower(zc, z4, k - 2 - 1, 4);
        z4f[k - 1] = cauchyprod(z4, f, k - 1);

        // x component
        gx = -(G * M * cauchyprod(xc, b, k - 1) +  
        J2 * ((3 / 2) * cauchyprod(xc, d, k - 1) -              // J2
        (3 / 2) * 5 * cauchyprod(xc, zsqe, k - 1)) +            // J2
        + 3 * C12 * (cauchyprod(zc, d, k - 1) - 5 * cauchyprod(xsqz, e, k - 1)) -       // 21
        15 * S12 * cauchyprod(xyz, e, k - 1)                                            // 21
        + 3 * C22 * (2 * cauchyprod(xc, d, k - 1) - 5 * cauchyprod(xcub, e, k - 1) +    // 22
        5 * cauchyprod(xysq, e, k - 1)) + 6 * S22 * (cauchyprod(yc, d, k - 1)           // 22
        - 5 * cauchyprod(xsqy, e, k - 1)) +                 // 22
        J3 * ((35 / 2) * cauchyprod(xc, zcubf, k - 1) +     // J3
        - (15 / 2) * cauchyprod(xc, ze, k - 1)));           // J3

        // y component
        gy = -(G * M * cauchyprod(yc, b, k - 1) +  
        J2 * ((3 / 2) * cauchyprod(yc, d, k - 1) -          //J2
        3 / 2) * 5 * cauchyprod(yc, zsqe, k - 1)) +         // J2
        3 * S12 * (cauchyprod(zc, d, k - 1) - 5 * cauchyprod(ysqz, e, k - 1)) -         // 21
        15 * C12 * cauchyprod(xyz, e, k - 1)                                            // 21
        + 3 * C22 * (-2 * cauchyprod(yc, d, k - 1) + 5 * cauchyprod(ycub, e, k - 1) -   // 22
        5 * cauchyprod(xsqy, e, k - 1)) + 6 * S22 * (cauchyprod(xc, d, k - 1)           // 22
        - 5 * cauchyprod(xysq, e, k - 1)) +                 // 22
        J3 * ((35 / 2) * cauchyprod(yc, zcubf, k - 1) +     // J3
        - (15 / 2) * cauchyprod(yc, ze, k - 1));           // J3


        // z component
        gz = -(G * M * cauchyprod(zc, b, k - 1) +  
        J2 * ((3 / 2) * 3 * cauchyprod(zc, d, k - 1) -      // J2
        (3 / 2) * 5 * cauchyprod(zc, zsqe, k - 1))          // J2
        + 3 * C12 * (cauchyprod(xc, d, k - 1) - 5 * cauchyprod(xc, zsqe, k - 1))     // 21
        + 3 * S12 * (cauchyprod(yc, d, k - 1) - 5 * cauchyprod(yc, zsqe, k - 1))     // 21
        + 15 * C22 * (cauchyprod(xsqz, e, k - 1) - cauchyprod(ysqz, e, k - 1))       // 22
        + 30 * S22 * cauchyprod(xyz, e, k - 1)          // 22
        + J3 * (-(35 / 2) * z4f[k - 1] +                // J3
        15 * zsqe[k - 1] - (3 / 2) * d[k - 1]));        // J3
    }


}