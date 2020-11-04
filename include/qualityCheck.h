#pragma once
#include "DataStructures.h"
#include <string>
#include <iomanip>
#include <algorithm>
#include <cmath>

using namespace std;

class qualityCheck
{
public:
	qualityCheck();
	qualityCheck(double **coorp_, vector<int> nu_, int np_);
	bool invertedElements2D(double **coorp);
	void meshQuality2D(double **coor, vector<int> &nu);
	void vecProd(double v1x, double v1y, double v1z, double v2x,
				 double v2y, double v2z, double &xre, double &yre, double &zre);
	void shapeMetric(double **cl, int ntet, string eleFile);
	void aspectRatio(double **cl, int ntet, string eleFile);

	void Jacobian(double **cl, int ntet, string eleFile);
	~qualityCheck();

private:
	int np;
	vector<int> nu;
	double **coorp;
};
