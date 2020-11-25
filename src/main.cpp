#include <iostream>
#include <math.h>
#include "DataStructures.h"
#include "ioGrid.h"
#include "qualityCheck.h"
#include "Solver.h"
#include "geomHandle.h"

using namespace std;
int main()
{
	int dim;

	cout << "==== CHOOSE 2D (2) or 3D (3) to solve ====" << endl;

	cin >> dim;

	switch (dim)
	{
	case 2:
	{
		// 1. CREATE GRID INFO STRUCT
		// gridInfo2D_Unstructured ginfo_un;
		gridInfo gInfo;
		ioGrids ionGrid;

		string filename = "box_in_box";
		// 2. READ INITIAL MESH

		ionGrid.read(filename);

		gInfo = ionGrid.getInfo();

		// 3. INITIAL MESH QUALITY CHECK
		qualityCheck gQuality(gInfo.coor, gInfo.nu, gInfo.np);
		gQuality.meshQuality2D(gInfo.coor, gInfo.nu);

		// 4. MOVE GEOMETRY'S BOUNDARY
		// ionGrid.readDeformed2D_Unstructured(filename);

		// gInfo = ionGrid.getInfo();

		// 5. CREATE DATA STRUCTURE OF UNSTRUCTURED GRID
		DataStructures dataStructure2D(gInfo.ns, gInfo.np, gInfo.nq, gInfo.coor, gInfo.logfr, gInfo.nu);
		dataStructure2D.Create2D();

		gInfo = ionGrid.getInfo();
		gInfo.pitch = dataStructure2D.get_pitch();
		gInfo.nper = dataStructure2D.get_nper();
		ionGrid.checkPitch(gInfo.coor, dataStructure2D.get_iper(), gInfo.pitch, gInfo.nper);

		geomHandle GeometryHandle(gInfo.coor, gInfo.logfr, gInfo.ns);
		// GeometryHandle.wingBending(20); // alpha = 0.1
		// GeometryHandle.wingTorsionBending(1); // alpha = 0.05
		GeometryHandle.rotateGeometry(0., 0., 10.0);
		gInfo.coorp = GeometryHandle.get_movement();
		ionGrid.vtk_graphics_2D_unstr(gInfo.coorp, filename + "_final");

		// 6. ADAPTATION OF INITIAL MESH WITH RBM METHOD
		// Numerics newtonRaphson2D(ginfo_un.coorp, ginfo_un.coor, ginfo_un.logfr, dataStructure2D.get_ndeg(),
		// 	dataStructure2D.get_jaret(), ginfo_un.nu, ginfo_un.ns, ginfo_un.np);

		Numerics newtonRaphson2D(gInfo.coorp, gInfo.coor, dataStructure2D.get_iper(), dataStructure2D.get_logfr(), dataStructure2D.get_ndeg(),
								 dataStructure2D.get_jaret(), gInfo.nu, gInfo.ns, gInfo.np, dataStructure2D.get_pitch(), gInfo.nper);

		newtonRaphson2D.Solver(3000);
		// gInfo = ionGrid.getInfo2D_Unstructured();
		// gInfo = ionGrid.getInfo();
		gInfo.pitch = dataStructure2D.get_pitch();
		gInfo.nper = dataStructure2D.get_nper();

		// 7. MESH QUALITY OF DEFORMED GEOMETRY MESH
		//gQuality.meshQuality2D(gInfo.coorp, gInfo.nu);

		// 8. WRITE new.nod file
		//ionGrid.write_mesh(dim, gInfo.coorp, gInfo.coor, gInfo.logfr, gInfo.ns, "siev.nod");
		// ionGrid.write_Un3D(gInfo.coorp, "sievGnu.dat");
		ionGrid.vtk_graphics_2D_unstr(gInfo.coorp, filename + "_final");
		ionGrid.vtk_graphics_2D_unstr(gInfo.coor, filename + "_initial");
		ionGrid.checkPitch(gInfo.coorp, dataStructure2D.get_iper(), gInfo.pitch, gInfo.nper);
		break;
	}
	case 3:
	{
		// 1. CREATE GRID INFO STRUCT
		gridInfo Mesh3D;
		ioGrids ioGrid;

		string filename = "mesh1";

		// 2. READ INITIAL MESH
		// ioGrid.readInitialMesh("wing.ele", "NONE", "wing.nod");
		ioGrid.read(filename);
		Mesh3D = ioGrid.getInfo();

		//ioGrid.vtk_graphics_3D_unstr();

		// 3. INITIAL MESH QUALITY CHECK
		qualityCheck MeshQuality;
		MeshQuality.shapeMetric(Mesh3D.coor, Mesh3D.ntet, filename);
		//MeshQuality.aspectRatio(Mesh3D.coor, Mesh3D.ntet, "wing");
		// MeshQuality.Jacobian(Mesh3D.coor, Mesh3D.ntet, filename);

		// 4. CREATE DATA STRUCTURE OF UNSTRUCTURED GRID
		DataStructures DS3D(Mesh3D.logfr, Mesh3D.ns, Mesh3D.coor,
							Mesh3D.ntet, Mesh3D.npyr, Mesh3D.npri, Mesh3D.nhex, Mesh3D.nall, Mesh3D.nu);
		DS3D.Create3D();

		// 5. MOVE GEOMETRY'S BOUNDARY

		geomHandle GeometryHandle(Mesh3D.coor, Mesh3D.logfr, Mesh3D.ns);
		// GeometryHandle.wingBending(20); // alpha = 0.1
		// GeometryHandle.wingTorsionBending(1); // alpha = 0.05
		GeometryHandle.rotateGeometry(15., 15., 10.0);
		GeometryHandle.translateGeometry(0.0, 0.25, -0.0);

		Mesh3D.coorp = GeometryHandle.get_movement();

		// ioGrid.vtk_graphics_3D_unstr(Mesh3D.coorp);
		// 6. ADAPTATION OF INITIAL MESH WITH RBM METHOD
		Numerics NewtonRaphson3D(Mesh3D.coorp, Mesh3D.coor, DS3D.get_iper(), DS3D.get_logfr(), DS3D.get_ndeg(), DS3D.get_jaret(),
								 Mesh3D.nu, Mesh3D.ns, 0, DS3D.get_pitch(), Mesh3D.nper);

		NewtonRaphson3D.Solver(3000);

		// 7. MESH QUALITY OF DEFORMED GEOMETRY MESH
		//MeshQuality.meshQuality3D(Mesh3D.coorp, Mesh3D.ntet, "wing");
		MeshQuality.aspectRatio(Mesh3D.coorp, Mesh3D.ntet, filename);

		// 8. WRITE new.nod file
		ioGrid.write_mesh(dim, Mesh3D.coorp, Mesh3D.coor, Mesh3D.logfr, Mesh3D.ns, "new.nod");

		// 9. WRITE .VTK file for paraview
		// ioGrid.write_Un3D(Mesh3D.coorp, "gr1d");
		ioGrid.vtk_graphics_3D_unstr(Mesh3D.coorp);
		break;
	}
	default:
		break;
	}

	std::cin.get();
	return 0;
}
