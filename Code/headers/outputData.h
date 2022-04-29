// Chooses how to output the data for each different solving Method
// 
// Triton Systems inc. By Colin Tierney

//libaries
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iomanip>

using namespace std;

void fileClose(ofstream& f)
{
	f.close();
}

void outputData(ofstream& f, int method, vector<double>& xPSM, vector<double>& yPSM, vector<double>& zPSM, vector<double>& vxPSM, vector<double>& vyPSM,
	vector<double>& vzPSM, vector<double>& xc, vector<double>& yc, vector<double>& zc, vector<double>& vxc, vector<double>& vyc, vector<double>& vzc,
	vector <double>& psmOrder, vector <double>& timeStep, int ns, int endTime, vector<double>& ta, int mDeg, vector <int>& degM, double& h, int dt, int choice)
{
	//formatting variables
	int w = 10;
	int TG = 161;

	//FSFO Output
	if (method == 1)
	{
		//TGx file size comparison output
		if (choice == 1)
		{
			//output time index
			if (ns * dt <= dt)
				f << TG << ",";
			else
				f << TG + dt * (ns - 1) << ",";

			//output x component Maclaurin Coefficients
			for (int i = 1; i <= mDeg + 1; i++)
			{
				f << setprecision(w) << xc[i] << ",";
			}
			//output y component Maclaurin Coefficients
			for (int i = 1; i <= mDeg + 1; i++)
			{
				f << setprecision(w) << yc[i] << ",";
			}
			//output z component Maclaurin Coefficients
			for (int i = 1; i <= mDeg + 1; i++)
			{
				f << setprecision(w) << zc[i] << ",";
			}
			f << mDeg << endl;
		}

		//output choice
		if (choice == 2)
		{
			//COLUMN LABELS
			if (ns == 1)
			{
				f << "X Coordinate" << ",";
				//f << "Y Coordinate" << ",";
				//f << "Z Coordinate" << ",";
				//f << "X Velocity" << ",";
				//f << "Y Velocity" << ",";
				//f << "Z Velocity" << ",";
				f << "PSM Order" << ",";
				f << "Time Step";
				f << endl;
			}

			//OUTPUTS
			f << setprecision(w) << xPSM[ns] << ",";			//x coord
			//f << setprecision(w) << yPSM[ns] << ",";			//y coord
			//f << setprecision(w) << zPSM[ns] << ",";			//z coord
			//f << setprecision(w) << vxPSM[ns] << ",";		//x velocity
			//f << setprecision(w) << vyPSM[ns] << ",";		//y velocity
			//f << setprecision(w) << vzPSM[ns] << ",";		//z velocity
			f << setprecision(w) << mDeg << ",";				//psm order at every index
			f << setprecision(w) << h;							//timestep at every index
			f << endl;
		}
	}

	//ASFO Output
	else if (method == 2)
	{
		//TGx file size comparison output
		if (choice == 1)
		{
			//output time index
			//if (ns * h <= h)
				//f << TG << ",";
			//else
			f << TG + h * (ns - 1) << ",";

			//output x component Maclaurin Coefficients
			for (int i = 1; i <= mDeg + 1; i++)
			{
				f << setprecision(w) << xc[i] << ",";
			}
			//output y component Maclaurin Coefficients
			for (int i = 1; i <= mDeg + 1; i++)
			{
				f << setprecision(w) << yc[i] << ",";
			}
			//output z component Maclaurin Coefficients
			for (int i = 1; i <= mDeg + 1; i++)
			{
				f << setprecision(w) << zc[i] << ",";
			}
			f << mDeg << endl;
		}

		//output choice
		if (choice == 2)
		{
			//COLUMN LABELS
			if (ns == 1)
			{
				f << "X Coordinate" << ",";
				//f << "Y Coordinate" << ",";
				//f << "Z Coordinate" << ",";
				//f << "X Velocity" << ",";
				//f << "Y Velocity" << ",";
				//f << "Z Velocity" << ",";
				f << "PSM Order" << ",";
				f << "Time Step" << ",";
				f << endl;
			}

			//OUTPUTS
			f << setprecision(w) << xPSM[ns] << ",";			//x coord
			//f << setprecision(w) << yPSM[ns] << ",";			//y coord
			//f << setprecision(w) << zPSM[ns] << ",";			//z coord
			//f << setprecision(w) << vxPSM[ns] << ",";		//x velocity
			//f << setprecision(w) << vyPSM[ns] << ",";		//y velocity
			//f << setprecision(w) << vzPSM[ns] << ",";		//z velocity
			f << setprecision(w) << mDeg << ",";				//psm order at every point
			f << setprecision(w) << h << ",";					//timestep at every point
			f << endl;
		}
	}

	//ASAO Ouput
	else if (method == 3)
	{
		//COLUMN LABELS
		if (ns == 1)
		{
			f << "X Coordinate" << ",";
			//f << "Y Coordinate" << ",";
			//f << "Z Coordinate" << ",";
			//f << "X Velocity" << ",";
			//f << "Y Velocity" << ",";
			//f << "Z Velocity" << ",";
			f << "PSM Order" << ",";
			f << "Time Step" << ",";
			f << endl;
		}

		//OUTPUTS
		f << setprecision(w) << xPSM[ns] << ",";			//x coord
		//f << setprecision(w) << yPSM[ns] << ",";			//y coord
		//f << setprecision(w) << zPSM[ns] << ",";			//z coord
		//f << setprecision(w) << vxPSM[ns] << ",";		//x velocity
		//f << setprecision(w) << vyPSM[ns] << ",";		//y velocity
		//f << setprecision(w) << vzPSM[ns] << ",";		//z velocity
		f << setprecision(w) << psmOrder[ns] << ",";		//psm order at every point
		f << setprecision(w) << timeStep[ns] << ",";		//timestep at every point
		f << endl;
	}
	if (ta[ns] > endTime)
		fileClose(f);
}

