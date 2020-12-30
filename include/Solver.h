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
	Numerics(double **coorp, double **coor, int **iper, vector<int> logfr, vector<int> ndeg, vector<int> jaret, vector<int> &nu,
			 int ns, int np, double pitch, double nper, int isperiph);
	void Solver(int Iterations);
	~Numerics();

private:
	//	----------------------------//
	/**/ double **dxi, **coorp, **coor;		 /**/
	/**/ vector<int> logfr, ndeg, jaret, nu; /**/
	/**/ vector<double> cM;
	/**/ int **iper;
	//	----------------------------//

	int ndim, ns, np, mID, neiTot, iter, inei, isperiph;
	double theta, cc, ss, sx, sy, dx, dy, pitch, nper;
	double aTerm, bTerm, xnei, ynei, xneiP, yneiP,
		xrms, yrms, xErrMax, yErrMax, xi, yi, xi2r, yi2r;
	double fun, dfun, dxOld, dyOld, thOld;
	double tiny = 1e-25;

	int kountNeis(int nodeID);
	void globalToLocal(int ndim, vector<double> &cM, int nodeID);
	void localToGlobal(int ndim, vector<double> &cM, int nodeID);
	void addpitch(int mid, int lgfr);
	void removepitch(int mid, int lgfr);
	void update_periodic_nodes(int mid);
	bool converged(double dxOld, double dyOld, double thOld,
				   double dx, double dy, double theta);
	int getQuarter(double x, double y);
	double getAngle(double x, double y);
	void addangle(int mid, double x, double y);
};
