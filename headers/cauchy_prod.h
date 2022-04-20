// Cauchy Product Algorithm
// Colin Tierney
//
//
// The function "<double> cauchyProd()"
// is the function for solving for the cauchy product of two sets of numbers


#include <iostream>
#include <vector>

using namespace std;

double cauchyprod(const vector<double>& a, const vector<double>& b, int n)
{
    //declare and initialize vector up to "n" with zeros
    double cp = 0;

    //initialize variable x
    double x = 0;

    //solve for cauchy product
    for (int i = 0; i <= n - 1; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            x += a[j + 1] * b[n - j];
        }
        cp = x;
        x = 0;
    }

    //return vector of coefficients to main
    return cp;
}