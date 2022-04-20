// Cauchy Square Algorithm
// Colin Tierney
//
//
// The function "<double> cauchySquare()"
// is the function for solving for cauchy square of an array of numbers


#include <iostream>
#include <vector>

using namespace std;

double cauchysquare(vector <double>& a, int n)
{
    //declare return variable
    double cPS = 0;

    //solve for cauchy square product
    if (n % 2 == 0)
    {
        for (int j = 0; j <= n / 2 - 1; j++)
        {
            cPS = a[j + 1] * a[n - j + 1] + cPS;
        }
        cPS = 2 * cPS + a[n / 2 + 1] * a[n / 2 + 1];
    }
    else
    {
        for (int j = 0; j <= (n - 1) / 2; j++)
        {
            cPS = a[j + 1] * a[n - j + 1] + cPS;
        }
        cPS = 2 * cPS;
    }

    //return value of cPS to where function was called
    return cPS;
}
