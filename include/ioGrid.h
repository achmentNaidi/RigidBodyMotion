#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include "DataStructures.h"
#include <iomanip>
#include <string>
#include <map>

using namespace std;

struct patch
{
	int Nodes;
	vector<double> patch_node_ids;
	string BCType;
};

class ioGrids
{
public:
	ioGrids();

	void readInitial2D_Structured(string fileName);
	void write2D_Structured(gridInfo2D_Structured gI, string gridName);

	void readInitialMesh(string eleFile, string nodFile);				  // 2D MESH
	void readInitialMesh(string eleFile, string hybFile, string nodFile); // 3D MESH
	void read(string name);

	void readDeformed2D_Unstructured(string defFile);
	void write_mesh(int ndim, double **coorp, double **coor, vector<int> &logfr, int ns, string gridName);
	void write_Un3D(double **coorp, string name);
	void vtk_graphics_3D_unstr(double **coorp);
	void Stop(string message, string file, int line);
	void Stop(string message);
	bool readPatches(string filename);
	void setListOfPatches(ifstream nod[]);
	void readPatchesFile(string filename);

	gridInfo2D_Structured getInfo2D_Structured() { return gInfo_st; }
	gridInfo2D_Unstructured getInfo2D_Unstructured() { return gInfo_un; }
	gridInfo3D_Unstructured getInfo3D_Unstructured() { return gInfo3D_un; }
	gridInfo getInfo() { return gInfo; }
	~ioGrids();

private:
	void stringIndexing2D();

	string fileName;
	map<string, patch> patces;
	gridInfo2D_Structured gInfo_st;
	gridInfo2D_Unstructured gInfo_un;
	gridInfo3D_Unstructured gInfo3D_un;
	gridInfo gInfo;
};