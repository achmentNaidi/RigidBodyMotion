#pragma once

#include <iostream>
#include <vector>
#include "qualityCheck.h"
#include "DataStructures.h"
#include <iomanip>
#include <cmath>

using namespace std;

class Numerics
{
public:
	Numerics();
	Numerics(double **coorp, double **coor, vector<int> logfr, vector<int> ndeg, vector<int> jaret, vector<int> nu,
			 int ns, int np);
	void Solver(int Iterations);
	~Numerics();

private:
	//	----------------------------//
	/**/ double **dxi, **coorp, **coor;		 /**/
	/**/ vector<int> logfr, ndeg, jaret, nu; /**/
	/**/ vector<double> cM;					 /**/
	//	----------------------------//

	int ndim, ns, np, mID, neiTot, iter, inei;
	double theta, cc, ss, sx, sy, dx, dy;
	double aTerm, bTerm, xnei, ynei, xneiP, yneiP,
		xrms, yrms, xErrMax, yErrMax, xi, yi, xi2r, yi2r;
	double fun, dfun, dxOld, dyOld, thOld;
	double tiny = 1e-25;

	int kountNeis(int nodeID);
	void globalToLocal(int ndim, vector<double> &cM, int nodeID);
	void localToGlobal(int ndim, vector<double> &cM, int nodeID);
	bool converged(double dxOld, double dyOld, double thOld,
				   double dx, double dy, double theta);
	bool isPeriodic(int mID);
};
