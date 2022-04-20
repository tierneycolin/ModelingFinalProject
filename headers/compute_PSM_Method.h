// Compute PSM Method
// Colin Tierney
// 
// Chooses which type of PSM Method to use for compute trajectory files


//libaries
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iomanip>

using namespace std;

void computeMethodPSM(int method, vector<double>& xc, vector<double>& yc, vector<double>& zc, vector<double>& vxc, vector<double>& vyc,
    vector<double>& vzc, int k, int mDeg, double tol, int maxStep, double& h, int ns, vector<int>& degM)
{
    //Fixed Time Step Fixed Order
    //if method == 1
    if (method == 1)
    {
        degM[ns] = mDeg;
        double maxValo = 0;
    }

    //Adaptive Time Step Fixed Order
    //if method == 2
    else if (method == 2)
    {
        //define variables
        double maxVal = 0;

        //for syntax; to find abs value below
        double absVar[6] = { xc[k], yc[k], zc[k], vxc[k], vyc[k], vzc[k] };

        //calculate abs value of each number in the array
        for (int i = 0; i <= 5; i++)
            absVar[i] = abs(absVar[i]);
        
        //find largest value in the array
        for (int i = 0; i <= 5; i++)
        {
            if (absVar[i] > maxVal)
                maxVal = absVar[i];
        }

        //find the abs value of the following
        double absV1 = pow((tol / (2.0 * maxVal)), (1/ static_cast <double> (k)));
        double absV2 = abs(maxStep);

        //find the smaller element and put the value in h
        h = min(absV1, absV2);

        //place k into degM at place ns
        degM[ns] = k;

        //for formatting the function; not used here but used in ASAO function
        double maxValo = 0;
    }

    //Adaptive Step Adaptive Order
    //if method == 3
    else if (method == 3)
    {
        //define variables
        double flag;
        double maxValo = 0, maxVal = 0;

        //for syntax; to find abs value below
        double absVar[6] = { xc[k], yc[k], zc[k], vxc[k], vyc[k], vzc[k] };

        //calculate abs value of each number in the array
        for (int i = 0; i <= 5; i++)
            absVar[i] = abs(absVar[i]);

        //find largest value in the array
        for (int i = 0; i <= 5; i++)
        {
            if (absVar[i] > maxVal)
                maxVal = absVar[i];
        }

        //find the abs value of the following
        double absV1 = pow((tol / (2.0 * maxVal)), (1 / static_cast <double> (k)));
        double absV2 = abs(maxStep);

        //find the smaller element and put the value in h
        h = min(absV1, absV2);

        //place k into degM at place ns
        degM[ns] = k;

        if (h == maxStep || maxVal * h < 0.1 * tol * maxValo || maxVal * pow(h, k) < 0.1 * tol)
        {
            degM[ns] = k;
            flag = 1;
        }

        //store variables
        maxValo = maxVal;
        degM[ns] = k;
        flag = 0;
    }
}
