#include "Solver.h"
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

Numerics::Numerics() {}

Numerics::Numerics(double **coorp_, double **coor_, vector<int> logfr_, vector<int> ndeg_, vector<int> jaret_,
				   vector<int> nu_, int ns_, int np_) : coorp(coorp_), coor(coor_), logfr(logfr_), ndeg(ndeg_), jaret(jaret_), nu(nu_), ns(ns_), np(np_) {}

void Numerics::Solver(int Iterations)
{
	bool testInt;

	qualityCheck triCheck(coorp, nu, np);

	int kIterOut = 0, iterMax;

	double u, v, w, ccu, ssu, ccv, ssv, ccw, ssw, **R,
		sx, sy, sz, dz, AuTerm, BuTerm, AvTerm, BvTerm,
		AwTerm, BwTerm, znei, zneiP, defx, defy, defz,
		funu, dfunu, funv, dfunv, funw, dfunw, dzOld, uOld,
		vOld, wOld, zrms, zErrMax, zi2r, zi, Zexist;

	int lime, lime1, kextra;

	ndim = 3;
	dxi = matrix<double>(ndim, ns);
	R = matrix<double>(ndim, ndim);

	// Read maxiter
	// ------------
	cout << "ENTER MAX NUMBER OF ITERATIONS" << endl;
	cout << " (DEFAULT = 1000 ITERATIONS)  " << endl;
	cout << "------------------------------ " << endl;
	cin >> lime;
	int maxiter = 1000;
	if (lime != ' ')
		cin >> maxiter;
	do
	{
		kIterOut++;
		for (int i = 0; i < ndim; i++)
		{
			for (int j = 0; j < ns; j++)
			{
				dxi[i][j] = coorp[i][j];
			}
		}
		/*
	     Sweep all internal grid nodes and compute deformations & rotations
	     ------------------------------------------------------------------
		*/
		for (unsigned int is = 1; is <= ns; is++)
		{
			if (logfr[is] != 0)
				continue; //  only internal nodes are of interest
			mID = is;	  // current mID node

			/*
		 Count neighbouring nodes
		 ------------------------
		*/
			neiTot = kountNeis(mID);
			globalToLocal(ndim, cM, mID);

			/*
			 ========================================
		     Compute Deformations and Rotations in 3D 
		     ========================================
			*/
			u = 0.;
			v = 0.;
			w = 0.;
			iter = 0;
			iterMax = 3;

			do
			{
				iter++;
				ccu = cos(u);
				ssu = sin(u);
				ccv = cos(v);
				ssv = sin(v);
				ccw = cos(w);
				ssw = sin(w);

				/*
				      Creating R table
				      ----------------
				*/
				R[0][0] = ccv * ccw;
				R[0][1] = -(ccu * ssw) + (ccw * ssu * ssv);
				R[0][2] = (ssu * ssw) + (ccu * ccw * ssv);
				R[1][0] = ccv * ssw;
				R[1][1] = (ccu * ccw) + (ssu * ssv * ssw);
				R[1][2] = -(ccw * ssu) + (ssv * ssw * ccu);
				R[2][0] = -ssv;
				R[2][1] = ccv * ssu;
				R[2][2] = ccu * ccv;

				//  Update X, Y, Z
				sx = 0.;
				sy = 0.;
				sz = 0.;

				for (unsigned int k = ndeg[mID - 1] + 1; k <= ndeg[mID]; k++)
				{
					inei = jaret[k];
					sx += coorp[0][inei - 1] - R[0][0] * coor[0][inei - 1] - R[0][1] * coor[1][inei - 1] - R[0][2] * coor[2][inei - 1];
					sy += coorp[1][inei - 1] - R[1][0] * coor[0][inei - 1] - R[1][1] * coor[1][inei - 1] - R[1][2] * coor[2][inei - 1];
					sz += coorp[2][inei - 1] - R[2][0] * coor[0][inei - 1] - R[2][1] * coor[1][inei - 1] - R[2][2] * coor[2][inei - 1];
				}
				dx = sx / double(neiTot);
				dy = sy / double(neiTot);
				dz = sz / double(neiTot);

				//  Update u, v, w
				AuTerm = 0.;
				BuTerm = 0.;
				AvTerm = 0.;
				BvTerm = 0.;
				AwTerm = 0.;
				BwTerm = 0.;

				for (unsigned int k = ndeg[mID - 1] + 1; k <= ndeg[mID]; k++)
				{
					inei = jaret[k];
					xnei = coor[0][inei - 1];
					ynei = coor[1][inei - 1];
					znei = coor[2][inei - 1];
					xneiP = coorp[0][inei - 1];
					yneiP = coorp[1][inei - 1];
					zneiP = coorp[2][inei - 1];

					defx = dx - xneiP;
					defy = dy - yneiP;
					defz = dz - zneiP;

					AuTerm = AuTerm - znei * (-defx * ssv * ccw + defy * ssv * ssw + defz * ccv) -
							 ynei * (defx * ssw + defy * ccw);
					BuTerm = BuTerm + ynei * (-defx * ssv * ccw + defy * ssv * ssw + defz * ccv) -
							 znei * (defx * ssw + defy * ccw);

					AvTerm = AvTerm - xnei * (ccw * defx - ssw * defy) - defz * (ccu * znei + ssu * ynei);
					BvTerm = BvTerm + (ccu * znei + ssu * ynei) * (ccw * defx - ssw * defy) + xnei * defz;

					AwTerm = AwTerm - (ccv * xnei + ssv * ssu * ynei + ssv * ccu * znei) * defx -
							 (ccu * ynei - ssu * znei) * defy;
					BwTerm = BwTerm + (ccv * xnei + ssv * ssu * ynei + ssv * ccu * znei) * defy -
							 (ccu * ynei - ssu * znei) * defx;
				}
				/*
				      Newton - Raphson
 			    	  ----------------
				*/
				funu = ssu * AuTerm + ccu * BuTerm;
				dfunu = ccu * AuTerm - ssu * BuTerm;
				//	if (abs(dfunu < -1e+61)) break;
				u += -funu / dfunu;

				funv = ssv * AvTerm + ccv * BvTerm;
				dfunv = ccv * AvTerm - ssv * BvTerm;
				//		if (abs(dfunv < -1e+40)) break;
				v += -funv / dfunv;

				funw = ssw * AwTerm + ccw * BwTerm;
				dfunw = ccw * AwTerm - ssw * BwTerm;
				//		if (abs(dfunw < -1e+40)) break;
				w += -funw / dfunw;

				dxOld = dx;
				dyOld = dy;
				dzOld = dz;
				uOld = u;
				vOld = v;
				wOld = w;

			} while (iter <= 3);

			/*
				 Update coordinates of point M
				 -----------------------------
			*/
			coorp[0][mID - 1] = cM[0] + dx;
			coorp[1][mID - 1] = cM[1] + dy;
			coorp[2][mID - 1] = cM[2] + dz;

			/*
			Local to Global w.r.t. point M
			------------------------------
			*/
			localToGlobal(ndim, cM, mID);
		}
		xrms = 0.;
		yrms = 0.;
		zrms = 0.;

		xErrMax = -1.e10;
		yErrMax = -1.e10;
		zErrMax = -1.e10;

		for (int is = 1; is <= ns; is++)
		{
			xi = dxi[0][is - 1] - coorp[0][is - 1];
			yi = dxi[1][is - 1] - coorp[1][is - 1];
			zi = dxi[2][is - 1] - coorp[2][is - 1];
			xi2r = sqrt(xi * xi);
			yi2r = sqrt(yi * yi);
			zi2r = sqrt(zi * zi);
			if (xi2r > xErrMax)
				xErrMax = xi2r;
			if (yi2r > yErrMax)
				yErrMax = yi2r;
			if (zi2r > zErrMax)
				zErrMax = zi2r;
			xrms += xi * xi;
			yrms += yi * yi;
			zrms += zi * zi;
		}
		xrms = sqrt(xrms / double(ns));
		yrms = sqrt(yrms / double(ns));
		zrms = sqrt(zrms / double(ns));
		Zexist = zrms;
		xrms = log10(xrms);
		yrms = log10(yrms);
		zrms = log10(zrms);
		xErrMax = log10(xErrMax);
		yErrMax = log10(yErrMax);
		zErrMax = log10(zErrMax);

		if (kIterOut % 100 == 0 || kIterOut <= 10)
		{
			if (Zexist == 0)
			{
				cout << setprecision(10) << kIterOut << " " << xrms << " " << yrms << endl;
				if (kIterOut > 10)
					testInt = triCheck.invertedElements2D(coorp);
			}
			else if (kIterOut > 10)
			{
				cout << setprecision(10) << kIterOut << " " << xrms << " " << yrms << " " << zrms << endl;
			}
		}

		if (kIterOut == maxiter)
		{
			cout << "ENTER ADDITIONAL ITERATIONS (0=STOP)" << endl;
			cout << "(DEFAULT = 1000 ITERATIONS" << endl;
			cout << "(IF NEGATIVE PRINTS TEMPORARY FILE FOR SAFETY)" << endl;
			cout << "----------------------------------------------" << endl;
			cin >> lime1;
			maxiter += 1000;
			if (lime1 != ' ')
				cin >> kextra;
			if (kextra <= 0)
				goto L30;
			else
				maxiter += -1000 + kextra;
		}

	} while (kIterOut < Iterations);
L30:;
}

int Numerics::kountNeis(int nodeId)
{

	int neiTot = 0;

	for (int i = ndeg[nodeId - 1] + 1; i <= ndeg[nodeId]; i++)
	{
		neiTot++;
	}

	return neiTot;
}

void Numerics::globalToLocal(int ndim, vector<double> &cM, int nodeID)
{

	int inei;
	cM.resize(ndim);

	for (int id = 0; id < ndim; id++)
	{
		cM[id] = coor[id][nodeID - 1]; // exw valei -1 sto coor[id][nodeID]
	}
	//     Global to Local w.r.t. point M
	for (int k = ndeg[nodeID - 1] + 1; k <= ndeg[nodeID]; k++)
	{
		inei = jaret[k];
		for (int id = 0; id < ndim; id++)
		{
			coor[id][inei - 1] = coor[id][inei - 1] - cM[id];	// exw valei -1 sto inei
			coorp[id][inei - 1] = coorp[id][inei - 1] - cM[id]; // exw valei -1 sto inei
		}
	}
}

void Numerics::localToGlobal(int ndim, vector<double> &cM, int nodeID)
{

	int inei;
	for (int k = ndeg[nodeID - 1] + 1; k <= ndeg[nodeID]; k++)
	{
		inei = jaret[k];
		for (int id = 0; id < ndim; id++)
		{
			coor[id][inei - 1] = coor[id][inei - 1] + cM[id];	// exw valei -1 sto inei
			coorp[id][inei - 1] = coorp[id][inei - 1] + cM[id]; // exw valei -1 sto inei
		}
	}
}

bool Numerics::converged(double dxOld, double dyOld, double thOld, double dx, double dy, double theta)
{

	double kconvT = 0;
	bool convergence = false;

	if (abs(dxOld - dx) < 1.e-14)
	{
		kconvT++;
	}
	if (abs(dyOld - dy) < 1.e-14)
	{
		kconvT++;
	}
	if (abs(thOld - theta) < 1.e-14)
	{
		kconvT++;
	}

	if (kconvT == 3)
	{
		convergence = true;
	}

	return convergence;
}

Numerics::~Numerics() {}