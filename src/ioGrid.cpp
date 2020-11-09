#include "ioGrid.h"

/*=====================================================
  =====================================================
  Created by Alexandros Eskantar and
			 Ahmet-Triantafilos Naidi
  Master students of Computational Mechanics,
  National Technichal University of Athens. All 
  rights reserved to Professor Kyriakos C. Giannakoglou
  and the Lab. of  THERMAL TURBOMACHINES PARALLEL CFD 
						& OPTIMIZATION UNIT
  ======================================================
  ======================================================
*/

ioGrids::ioGrids() {}

void ioGrids::readDeformed2D_Unstructured(string defFile)
{
	int dum, ns2;

	std::ifstream def;
	//Check if .ele file exists
	if (def.fail())
	{

		cout << "Non-existent def file" << endl;
		exit;
	}

	gInfo.coorp = matrix<double>(3, gInfo.ns);
	//Read.def file
	def.open("input_data/" + defFile + ".def");
	def >> ns2;
	if (gInfo.ns != ns2)
	{
		cout << "Incompatible derfomation file" << endl;
	}
	for (int i = 0; i < gInfo.ns; i++)
		def >> dum;
	for (int i = 0; i < gInfo.ns; i++)
	{
		def >> gInfo.coorp[0][i];
	}

	for (int i = 0; i < gInfo.ns; i++)
	{
		def >> gInfo.coorp[1][i];
	}

	for (int i = 0; i < gInfo.ns; i++)
	{
		gInfo.coorp[2][i] = 0;
	}

	def.close();
}

void ioGrids::read(string name)
{
	ostringstream fileName;

	fileName << "input_data/" + name;
	// >> .nod file (allocations)
	///////////////////////////////////////////////////////////////

	int ns;
	int _ASCII_ = 0, _BINARY_ = 1;

	// Connectivity
	///////////////////////////////////////////////////////////////
	std::ifstream hyb[2];
	std::ifstream ele[2];

	hyb[_ASCII_].open((fileName.str() + ".hyb").c_str(), std::ios::in);
	ele[_ASCII_].open((fileName.str() + ".ele").c_str(), std::ios::in);

	hyb[_BINARY_].open((fileName.str() + "_BINARY.hyb").c_str(), std::ios::in | std::ios::binary);
	ele[_BINARY_].open((fileName.str() + "_BINARY.ele").c_str(), std::ios::in | std::ios::binary);

	bool foundHyb = false;
	bool foundEle = false;
	if (hyb[_ASCII_].is_open() || hyb[_BINARY_].is_open())
		foundHyb = true;
	if (ele[_ASCII_].is_open() || ele[_BINARY_].is_open())
		foundEle = true;

	if (!foundEle && !foundHyb)
		Stop("ERROR: There is no mesh connectivity file", __FILE__, __LINE__);
	if (foundEle && foundHyb)
		Stop("ERROR: Both *.hyb and *.ele files are located in the mesh folder ", __FILE__, __LINE__);

	if (hyb[_ASCII_].is_open() || ele[_ASCII_].is_open())
		cout << "#"
			 << " Reading mesh connectivity in ASCII format ... " << endl;
	if (hyb[_BINARY_].is_open() || ele[_BINARY_].is_open())
	{
		std::cout << "#"
				  << " Reading mesh connectivity in BINARY format ... " << endl;
	}

	if (hyb[_ASCII_].is_open() ||
		hyb[_BINARY_].is_open())
	{
		if (hyb[_ASCII_].is_open())
		{
			hyb[_ASCII_] >> gInfo.np;
			hyb[_ASCII_] >> gInfo.nq;
			hyb[_ASCII_] >> gInfo.ntet;
			hyb[_ASCII_] >> gInfo.npyr;
			hyb[_ASCII_] >> gInfo.npri;
			hyb[_ASCII_] >> gInfo.nhex;
		}
		else
		{
			hyb[_BINARY_].read((char *)&(gInfo.np), sizeof(int));
			hyb[_BINARY_].read((char *)&(gInfo.nq), sizeof(int));
			hyb[_BINARY_].read((char *)&(gInfo.ntet), sizeof(int));
			hyb[_BINARY_].read((char *)&(gInfo.npyr), sizeof(int));
			hyb[_BINARY_].read((char *)&(gInfo.npri), sizeof(int));
			hyb[_BINARY_].read((char *)&(gInfo.nhex), sizeof(int));
		}
		if (gInfo.ntet < 0)
			Stop("Invalid number of tetrahedra in .hyb file\n", __FILE__, __LINE__);
		if (gInfo.npyr < 0)
			Stop("Invalid number of pyramids   in .hyb file\n", __FILE__, __LINE__);
		if (gInfo.npri < 0)
			Stop("Invalid number of prisms     in .hyb file\n", __FILE__, __LINE__);
		if (gInfo.nhex < 0)
			Stop("Invalid number of hexahedra  in .hyb file\n", __FILE__, __LINE__);

		const int nentries = gInfo.np * 3 +
							 gInfo.nq * 4 +
							 gInfo.ntet * 4 +
							 gInfo.npyr * 5 +
							 gInfo.npri * 6 +
							 gInfo.nhex * 8;

		gInfo.nu.resize(nentries + 1);
		// nu[0] = 0;

		if (hyb[_ASCII_].is_open())
		{
			int ind = 1;
			if (gInfo.np > 0)
				for (int ip = 1; ip <= 3 * gInfo.np; ip++, ind++)
					hyb[_ASCII_] >> gInfo.nu[ind];
			if (gInfo.nq > 0)
				for (int iq = 1; iq <= 3 * gInfo.np; iq++, ind++)
					hyb[_ASCII_] >> gInfo.nu[ind];
			if (gInfo.ntet > 0)
				for (int itet = 1; itet <= 4 * gInfo.ntet; itet++, ind++)
					hyb[_ASCII_] >> gInfo.nu[ind];
			if (gInfo.npyr > 0)
				for (int ipyr = 1; ipyr <= 5 * gInfo.npyr; ipyr++, ind++)
					hyb[_ASCII_] >> gInfo.nu[ind];
			if (gInfo.npri > 0)
				for (int ipri = 1; ipri <= 6 * gInfo.npri; ipri++, ind++)
					hyb[_ASCII_] >> gInfo.nu[ind];
			if (gInfo.nhex > 0)
				for (int ihex = 1; ihex <= 8 * gInfo.nhex; ihex++, ind++)
					hyb[_ASCII_] >> gInfo.nu[ind];
			// if (!istream_operational(hyb[_ASCII_]))
			// 	Stop("EOF/bad data in .hyb file\n");

			hyb[_ASCII_].close();
		}
		else
		{
			hyb[_BINARY_].read((char *)&(gInfo.nu[1]), nentries * sizeof(int));
			hyb[_BINARY_].close();
		}
	}

	if (ele[_ASCII_].is_open() ||
		ele[_BINARY_].is_open())
	{
		if (ele[_ASCII_].is_open())
		{
			ele[_ASCII_] >> gInfo.np;
			ele[_ASCII_] >> gInfo.nq;
			ele[_ASCII_] >> gInfo.ntet;
			ele[_ASCII_] >> gInfo.npyr;
			ele[_ASCII_] >> gInfo.npri;
			ele[_ASCII_] >> gInfo.nhex;
		}
		else
		{
			ele[_BINARY_].read((char *)&(gInfo.np), sizeof(int));
			ele[_BINARY_].read((char *)&(gInfo.nq), sizeof(int));
			ele[_BINARY_].read((char *)&(gInfo.ntet), sizeof(int));
			ele[_BINARY_].read((char *)&(gInfo.npyr), sizeof(int));
			ele[_BINARY_].read((char *)&(gInfo.npri), sizeof(int));
			ele[_BINARY_].read((char *)&(gInfo.nhex), sizeof(int));
		}
		if (gInfo.ntet < 0)
			Stop("Invalid number of tetrahedra in .ele file\n");

		const int nentries = gInfo.np * 3 +
							 gInfo.nq * 4 +
							 gInfo.ntet * 4 +
							 gInfo.npyr * 5 +
							 gInfo.npri * 6 +
							 gInfo.nhex * 8;

		// nu.host_alloc(nentries + 1, "SubMesh::Connectivity");
		// nu[0] = 0;
		gInfo.nu.resize(nentries + 1);
		if (ele[_ASCII_].is_open())
		{
			int ind = 1;
			if (gInfo.np > 0)
				for (int ip = 1; ip <= 3 * gInfo.np; ip++, ind++)
					ele[_ASCII_] >> gInfo.nu[ind];
			if (gInfo.nq > 0)
				for (int iq = 1; iq <= 3 * gInfo.np; iq++, ind++)
					ele[_ASCII_] >> gInfo.nu[ind];
			if (gInfo.ntet > 0)
				for (int itet = 1; itet <= 4 * gInfo.ntet; itet++, ind++)
					ele[_ASCII_] >> gInfo.nu[ind];
			if (gInfo.npyr > 0)
				for (int ipyr = 1; ipyr <= 5 * gInfo.npyr; ipyr++, ind++)
					ele[_ASCII_] >> gInfo.nu[ind];
			if (gInfo.npri > 0)
				for (int ipri = 1; ipri <= 6 * gInfo.npri; ipri++, ind++)
					ele[_ASCII_] >> gInfo.nu[ind];
			if (gInfo.nhex > 0)
				for (int ihex = 1; ihex <= 8 * gInfo.nhex; ihex++, ind++)
					ele[_ASCII_] >> gInfo.nu[ind];
			// if (!istream_operational(ele[_ASCII_]))
			// 	Stop("EOF/bad data in .ele \n");

			ele[_ASCII_].close();
		}
		else
		{
			ele[_BINARY_].read((char *)&(gInfo.nu[1]), nentries * sizeof(int));
			ele[_BINARY_].close();
		}
	}
	gInfo.nall = gInfo.np + gInfo.nq + gInfo.ntet + gInfo.npyr + gInfo.npri + gInfo.nhex;

	// read .nod and patch file
	ifstream nod[2];

	nod[_ASCII_].open((fileName.str() + ".nod").c_str(), std::ios::in);
	nod[_BINARY_].open((fileName.str() + "_BINARY.nod").c_str(), std::ios::in | std::ios::binary);

	if (!nod[_ASCII_].is_open() && !nod[_BINARY_].is_open())
	{
		Stop("Cannot open coordinates file", __FILE__, __LINE__);
	}

	if (nod[_ASCII_].is_open())
	{
		nod[_ASCII_] >> ns;
		gInfo.ns = ns;
	}
	else
		nod[_BINARY_].read((char *)&(ns), sizeof(int));
	gInfo.ns = ns;

	if (ns <= 0)
	{
		Stop("Invalid number of nodes in .nod file", __FILE__, __LINE__);
	}

	//.patch file (list of patches)
	///////////////////////////////////////////////////////////////
	if (!readPatches(fileName.str()))
	{ //  in case of old format ("logfr") set list of patches :
		setListOfPatches(nod);
	}
	else
	{
		readPatchesFile(fileName.str());
	}

	gInfo.coor = matrix<double>(3, gInfo.ns);

	// .nod file (coordinates)
	///////////////////////////////////////////////////////////////

	if (nod[_ASCII_].is_open())
	{
		for (int is = 1; is <= ns; is++)
			nod[_ASCII_] >> setprecision(9) >> gInfo.coor[0][is - 1];
		// coor[1][0] = 0.;

		for (int is = 1; is <= ns; is++)
			nod[_ASCII_] >> setprecision(9) >> gInfo.coor[1][is - 1];
		// coor[2][0] = 0.;
		for (int is = 1; is <= ns; is++)
		{
			if (gInfo.np > 0 || gInfo.nq > 0)
			{
				gInfo.coor[2][is - 1] = 0;
			}
			else
			{
				nod[_ASCII_] >> setprecision(9) >> gInfo.coor[2][is - 1];
			}
		}

		// if (!istream_operational(nod[_ASCII_]))
		// 	Stop("EOF/bad data in nod file\n");

		nod[_ASCII_].close();
	}
	else if (nod[_BINARY_].is_open())
	{
		nod[_BINARY_].read((char *)&(gInfo.coor[0][0]), ns * sizeof(double) - 1);
		nod[_BINARY_].read((char *)&(gInfo.coor[1][0]), ns * sizeof(double) - 1);
		if (gInfo.np != 0 || gInfo.nq != 0)
		{
			nod[_BINARY_].read((char *)&(gInfo.coor[2][0]), ns * sizeof(double) - 1);
		}
		nod[_BINARY_].close();
	}

	// ########### SCALE MESH #############
	// for (int is = 1; is <= ns; is++)
	//   { //  scale mesh :
	// 	  coor[0][is - 1] *= Global::input.mesh->scale;
	// 	  coor[1][is - 1] *= Global::input.mesh->scale;
	// 	  coor[2][is - 1] *= Global::input.mesh->scale;
	//   }
}

void ioGrids::write_mesh(int ndim, double **coorp, double **coor, vector<int> &logfr, int ns, string gridName)
{
	std::ofstream nod;
	nod.open("output_data/" + gridName);
	nod << ns << endl;
	for (int i = 1; i <= ns; i++)
	{
		nod << logfr[i] << " ";
	}
	nod << endl;
	for (size_t i = 0; i < ndim; i++)
	{
		for (int j = 0; j < ns; j++)
		{
			nod << coorp[i][j] << " ";
		}
		nod << endl;
	}

	/*for (int i = 0; i < ns; i++)
	{
		nod << coorp[0][i] << "	" << coorp[1][i] << " " << coorp[2][i] << endl;
	}
	nod.close();*/
}

void ioGrids::write_Un3D(double **coorp, string name)
{
	std::ofstream grid;
	grid.open("output_data/" + name);
	// grid << endl
	// 	 << endl;
	// grid << gInfo.ntet << endl;
	// grid << gInfo.ns << endl;

	for (int i = 0; i < gInfo.ns; i++)
	{
		grid << coorp[0][i] << " " << coorp[1][i] << " " << coorp[2][i] << endl;
	}
	grid.close();
}

void ioGrids::vtk_graphics_3D_unstr(double **coorp)
{
	const int ln = 1e6;
	int ntet, npop, ns, nt;
	// vector<int> nds, nelem;
	// vector<double> x, y, z;
	// string line;
	cout << "-------------------------------------" << endl;
	cout << "Writing VTK file..." << endl;
	// cin >> fileName;
	// ifstream inData;
	// inData.open("output_data/" + fileName);
	// getline(inData, line);
	// getline(inData, line);
	// inData >> ntet;
	// inData >> ns;
	if (gInfo.ns > ln)
	{
		cout << "INCREASE L" << endl;
		exit;
	}
	cout << "DIMENSIONS : [ Total Nodes " << gInfo.ns << " ]" << endl;
	cout << "-------------------------------------" << endl;
	// x.resize(ns);
	// y.resize(ns);
	// z.resize(ns);
	// nds.resize(ns);
	npop = 4 * ntet;
	// nelem.resize(npop);

	// for (int i = 0; i < ns; i++)
	// {
	// 	inData >> x[i] >> y[i] >> z[i];
	// }

	ofstream vtk;
	vtk.open("output_data/out3D_unstr.vtk");

	vtk << "# vtk DataFile Version 3.1" << endl;
	vtk << fileName << endl;
	vtk << "ASCII" << endl;
	vtk << "DATASET UNSTRUCTURED_GRID" << endl;
	vtk << "POINTS " << gInfo.ns << " DOUBLE" << endl;

	for (int i = 0; i < gInfo.ns; i++)
	{
		vtk << setprecision(9) << coorp[0][i] << "		" << coorp[1][i] << "		" << coorp[2][i] << endl;
	}

	int nsize = 5 * gInfo.ntet;
	int nnum1 = 4;
	int ntype1 = 10;
	vtk << "CELLS " << gInfo.ntet << " " << nsize << endl;
	// Open .ele file
	// cout << "-------------------------------------" << endl;
	// cout << " ENTER .ele FILE " << endl;
	// cin >> fileName;
	// ifstream elefile;
	// elefile.open("input_data/" + fileName);
	// getline(elefile, line);
	// elefile >> nt;

	// for (int i = 0; i < npop; i++)
	// {
	// 	elefile >> nelem[i];
	// }
	int kounter = 1;
	for (int i = 0; i < gInfo.ntet; i++)
	{
		vtk << nnum1 << " " << gInfo.nu[kounter] - 1 << " " << gInfo.nu[kounter + 1] - 1
			<< " " << gInfo.nu[kounter + 2] - 1 << " " << gInfo.nu[kounter + 3] - 1 << endl;
		kounter += 4;
	}

	vtk << "CELL_TYPES " << gInfo.ntet << endl;
	for (int i = 0; i < gInfo.ntet; i++)
	{
		vtk << ntype1 << endl;
	}
}

void ioGrids::vtk_graphics_2D_unstr(double **coorp, string name)
{
	const int ln = 1e6;
	int np, npop, ns, nt;
	// vector<int> nds, nelem;
	// vector<double> x, y, z;
	// string line;
	cout << "-------------------------------------" << endl;
	cout << "Writing VTK file..." << endl;
	// cin >> fileName;
	// ifstream inData;
	// inData.open("output_data/" + fileName);
	// getline(inData, line);
	// getline(inData, line);
	// inData >> ntet;
	// inData >> ns;
	if (gInfo.ns > ln)
	{
		cout << "INCREASE L" << endl;
		exit;
	}
	cout << "DIMENSIONS : [ Total Nodes " << gInfo.ns << " ]" << endl;
	cout << "-------------------------------------" << endl;
	// x.resize(ns);
	// y.resize(ns);
	// z.resize(ns);
	// nds.resize(ns);
	npop = 4 * np;
	// nelem.resize(npop);

	// for (int i = 0; i < ns; i++)
	// {
	// 	inData >> x[i] >> y[i] >> z[i];
	// }

	ofstream vtk;
	name = "output_data/" + name + ".vtk";
	vtk.open(name);

	vtk << "# vtk DataFile Version 3.1" << endl;
	vtk << "vtk output" << endl;
	vtk << "ASCII" << endl;
	vtk << "DATASET UNSTRUCTURED_GRID" << endl;
	vtk << "POINTS " << gInfo.ns << " DOUBLE" << endl;

	for (int i = 0; i < gInfo.ns; i++)
	{
		vtk << setprecision(9) << coorp[0][i] << "		" << coorp[1][i] << " " << coorp[2][i] << endl;
	}

	int nsize = 4 * gInfo.np;
	int nnum1 = 3;
	int ntype1 = 5;
	vtk << "CELLS " << gInfo.np << " " << nsize << endl;

	int kounter = 1;
	for (int i = 0; i < gInfo.np; i++)
	{
		vtk << nnum1 << " " << gInfo.nu[kounter] - 1 << " " << gInfo.nu[kounter + 1] - 1
			<< " " << gInfo.nu[kounter + 2] - 1 << endl;
		kounter += 3;
	}

	vtk << "CELL_TYPES " << gInfo.np << endl;
	for (int i = 0; i < gInfo.np; i++)
	{
		vtk << ntype1 << endl;
	}
}

void ioGrids::Stop(string message, string file, int line)
{
	cout << message << file << line << endl;
	exit;
}

void ioGrids::Stop(string message)
{
	cout << message << endl;
	exit;
}

bool ioGrids::readPatches(string filename)
{
	ifstream patch;
	filename = "input_data/" + filename;
	patch.open(filename);

	bool foundPatch = false;
	if (patch.is_open())
		foundPatch = true;
	if (!foundPatch)
	{
		cout << "Patch file does not exist!!" << endl
			 << endl;
		// Stop("ERROR: There is no mesh patch file", __FILE__, __LINE__);
	}
	return foundPatch;
}

void ioGrids::readPatchesFile(string filename)
{
	int nbPatches;
	string line, patchName;

	ifstream file;
	filename = "input_data/" + filename + ".patch";
	file.open(filename);

	file >> nbPatches;

	patch Patch;

	int counter = 1;
	while (counter <= nbPatches)
	{

		file >> patchName;
		line = patchName;
		while (line != "}")
		{
			file >> line;
			if (line == "BCType")
			{
				file >> Patch.BCType;
			}
			if (line == "Nodes")
			{
				file >> Patch.Nodes;
				Patch.patch_node_ids.resize(Patch.Nodes + 1);
				for (int i = 1; i <= Patch.Nodes; i++)
				{
					file >> Patch.patch_node_ids[i];
				}
			}
		}
		patces[patchName] = Patch;
		counter++;
	}
}

void ioGrids::setListOfPatches(ifstream nod[])
{
	int ASCII = 0, BINARY = 1;
	gInfo.logfr.resize(gInfo.ns + 1);
	for (unsigned int i = 1; i <= gInfo.ns; i++)
		nod[ASCII] >> gInfo.logfr[i];
}

void ioGrids::checkPitch(double **coorp_, int **iper, double pitch)
{

	ofstream difPer;
	difPer.open("output_data/pitchDiference.dat");
	for (int i = 0; i < 121; i++)
	{
		difPer << setprecision(9) << coorp_[0][iper[1][i] - 1] << " " << abs(coorp_[1][iper[1][i] - 1] - coorp_[1][iper[0][i] - 1] - pitch) << endl;
	}
}
ioGrids::~ioGrids() {}