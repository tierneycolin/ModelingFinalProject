// Cauchy Power Function
// Colin Tierney
//
//
// The function "<double> cauchyPower()"
// is the function for solving for the cauchy power product of a set of numbers


#include <iostream>
#include <vector>

using namespace std;

double cauchypower(vector <double>& a, vector <double>& b, int n, double p)
{
    //declare/initialize vector
    double cauchyPower = 0;

    //solve for cauchy power 
    for (int j = 0; j <= n - 1; j++)
    {
        cauchyPower = (j * (p + 1) + p - n) * a[j + 2] * b[n - j + 1] + cauchyPower;
    }
    cauchyPower = ((n + 1) * p * a[n + 2] * b[1] + cauchyPower) / ((n + 1) * a[1]);

    //return vector of new values to main function
    return cauchyPower;
}