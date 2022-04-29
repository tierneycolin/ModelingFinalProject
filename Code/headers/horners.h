// Horner's Method Algorithm
// Colin Tierney
//
//
// The function "double horners()"
// is the function for solving for the horner's method for a specified array/vector of numbers
//
// Test numbers come from the "PSM_coeffs.xlsx" file


#include <iostream>
#include <vector>

using namespace std;

double horners(double t, vector <double>& a, int n)
{
    //declare variable
    double h;

    //intilialize h for horners method
    h = a[n - 1] + a[n] * t;

    //solves using horners method and puts values into a[]
    for (int j = 1; j <= n - 2; j++)
    {
        h = a[n - j - 1] + t * h;
    }

    return h;
}