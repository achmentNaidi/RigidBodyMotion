#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;

struct centroid
{
	double x;
	double y;
};

class geomHandle
{
public:
	geomHandle(double **coor, vector<int> logfr, int ns);

	void wingBending(double alpha);

	void wingTorsionBending(double alpha);

	void translateGeometry(double dx, double dy, double dz);

	void rotateGeometry(double u, double v, double w);

	double **get_movement()
	{
		return coorp;
	}

	centroid centerOfGravity(vector<double> x, vector<double> y);

	void read2dGeometry(string fileNmame, vector<double> &x_g, vector<double> &y_g);

	void rotate2D(vector<double> &x, vector<double> &y, double cx, double cy, double angle);

	void getGeometry();

	void overwriteDefFile(vector<double> xg, vector<double> yg, string fileName);

	~geomHandle();

private:
	double **coor, **coorp;
	vector<double> x_g, y_g;
	vector<int> logfr;
	int ns;
};
