#include "qualityCheck.h"

qualityCheck::qualityCheck(){};
qualityCheck::qualityCheck(double **coorp_, vector<int> nu_, int np_)
	: coorp(coorp_), nu(nu_), np(np_) {}

bool qualityCheck::invertedElements2D(double **cl)
{
	bool invElnts = false;
	int is1, is2, is3, ind = 0;
	double dx12, dx13, dy12, dy13, criterion;
	double dxa, dya, dza = 0.;
	int kountInv = 0; // number of inverted triangles
	ind = 0;
	for (int ip = 0; ip < np; ip++)
	{
		is1 = nu[ind + 1] - 1; // exw valei -1
		is2 = nu[ind + 2] - 1; // exw valei -1
		is3 = nu[ind + 3] - 1; // exw valei -1
		dx12 = cl[0][is2] - cl[0][is1];
		dx13 = cl[0][is3] - cl[0][is1];
		dy12 = cl[1][is2] - cl[1][is1];
		dy13 = cl[1][is3] - cl[1][is1];
		vecProd(dx12, dy12, 0., dx13, dy13, 0., dxa, dya, dza);
		criterion = dza * .5;
		if (criterion < 0.)
		{
			kountInv++;
		}
		ind = ind + 3;
	}
	if (kountInv > 0)
	{
		cout << "Inverted triangles : " << kountInv << endl;
		invElnts = true;
	}
	return invElnts;
}

void qualityCheck::vecProd(double v1x, double v1y, double v1z, double v2x,
						   double v2y, double v2z, double &xre, double &yre, double &zre)
{

	xre = v1y * v2z - v1z * v2y;
	yre = v1z * v2x - v1x * v2z;
	zre = v1x * v2y - v1y * v2x;
}

void qualityCheck::meshQuality2D(double **coor, vector<int> &nu)
{
	double qMin, qMax, qMean, qLoc, qSTD;
	double riza3 = sqrt(3.0);
	double wr11 = 1.;
	double wr12 = -1. / riza3;
	double wr21 = 0.;
	double wr22 = 2. / riza3;
	double a11, a12, a21, a22, s11, s12, s21, s22;
	double aja, trace;
	int is1, is2, is3;
	vector<int> nod;
	nod.resize(3);
	vector<double> qT, qual;
	qT.resize(3);
	qual.resize(np);
	int ind = 0;
	int kountInv = 0; // number of inverted triangles

	for (unsigned int ip = 0; ip < np; ip++)
	{
		nod[0] = nu[ind + 1];
		nod[1] = nu[ind + 2];
		nod[2] = nu[ind + 3];
		is1 = nod[0];
		is2 = nod[1];
		is3 = nod[2];
		ind += 3;

		qT[0] = 100.;
		qT[1] = 100.;
		qT[2] = 100.; // set a high (bad) value
					  //
					  //		the three nodes of the triangle
		for (unsigned int j = 0; j < 3; j++)
		{
			if (j == 1)
			{
				is1 = nod[0];
				is2 = nod[1];
				is3 = nod[2];
			}
			if (j == 2)
			{
				is1 = nod[0];
				is2 = nod[1];
				is3 = nod[2];
			}

			a11 = coor[0][is2 - 1] - coor[0][is1 - 1];
			a12 = coor[0][is3 - 1] - coor[0][is1 - 1];
			a21 = coor[1][is2 - 1] - coor[1][is1 - 1];
			a22 = coor[1][is3 - 1] - coor[1][is1 - 1];

			s11 = a11 * wr11 + a12 * wr21;
			s12 = a11 * wr12 + a12 * wr22;
			s21 = a21 * wr11 + a22 * wr21;
			s22 = a21 * wr12 + a22 * wr22;

			aja = s11 * s22 - s12 * s21;
			trace = s11 * s11 + s12 * s12 + s21 * s21 + s22 * s22;

			if (aja < 0)
			{
				kountInv++;
				exit;
			}
			qT[j] = aja / trace;
		}

		qual[ip] = 2. * (qT[0] + qT[1] + qT[2]) / 3.;
	}

	qMin = 100.;
	qMax = -100.;
	qMean = 0.;
	for (unsigned int ip = 0; ip < np; ip++)
	{
		qLoc = qual[ip];
		qMean = qMean + qLoc;
		if (qMin > qLoc)
			qMin = qLoc;
		if (qMax < qLoc)
			qMax = qLoc;
	}
	double a, b;
	a = qual[0];
	b = qual[np - 1];
	qMean = qMean / np;

	qSTD = 0.;
	for (unsigned int ip = 0; ip < np; ip++)
	{
		qSTD = qSTD + (qual[ip] - qMean) * (qual[ip] - qMean);
	}
	qSTD = sqrt(qSTD / np);

	cout << " Quality Metrics  " << endl;
	cout << " ---------------- " << endl;
	cout << " Min: " << qMin << "	Max: " << qMax << endl;
	cout << " Mean: " << qMean << "	std: " << qSTD << endl;
	cout << "  ------  ------  ------  ------  ------ " << endl;
}

void qualityCheck::shapeMetric(double **cl, int ntet, string eleFile)
{

	string line;
	int is0, is1, is2, is3;
	int nod[4];
	double x0, x1, x2, x3, y0, y1, y2, y3, z0, z1, z2, z3,
		deta, s11, s22, s33, s12, s23, s13, sum1, sum2, anum,
		qMin, qMax, qMean, qLoc, den;
	vector<double> qual;
	qual.resize(ntet + 1);

	eleFile += ".ele";
	ifstream elefile;

	elefile.open("input_data/" + eleFile);

	if (elefile.fail())
	{
		cout << " Non-existent .ele file " << eleFile << endl;
		exit;
	}
	getline(elefile, line);
	int kountInv = 0; // number of inverted tetrahedrons
	int kount = 0;
	double riza2 = sqrt(2.);

	ofstream qualityFile;
	qualityFile.open("output_data/Mesh_Quality");
	qualityFile << endl
				<< endl;

	for (int ip = 1; ip <= ntet; ip++)
	{
		elefile >> is0 >> is1 >> is2 >> is3;
		x0 = cl[0][is0 - 1];
		y0 = cl[1][is0 - 1];
		z0 = cl[2][is0 - 1];
		x1 = cl[0][is1 - 1];
		y1 = cl[1][is1 - 1];
		z1 = cl[2][is1 - 1];
		x2 = cl[0][is2 - 1];
		y2 = cl[1][is2 - 1];
		z2 = cl[2][is2 - 1];
		x3 = cl[0][is3 - 1];
		y3 = cl[1][is3 - 1];
		z3 = cl[2][is3 - 1];

		deta = (x1 - x0) * ((y2 - y0) * (z3 - z0) - (y3 - y0) * (z2 - z0)) -
			   (x2 - x0) * ((y1 - y0) * (z3 - z0) - (y3 - y0) * (z1 - z0)) +
			   (x3 - x0) * ((y1 - y0) * (z2 - z0) - (y2 - y0) * (z1 - z0));

		s11 = pow(x1 - x0, 2) + pow(y1 - y0, 2.) + pow(z1 - z0, 2);
		s22 = pow(x2 - x0, 2) + pow(y2 - y0, 2.) + pow(z2 - z0, 2);
		s33 = pow(x3 - x0, 2) + pow(y3 - y0, 2.) + pow(z3 - z0, 2);
		;

		s12 = (x1 - x0) * (x2 - x0) + (y1 - y0) * (y2 - y0) + (z1 - z0) * (z2 - z0);
		s23 = (x2 - x0) * (x3 - x0) + (y2 - y0) * (y3 - y0) + (z2 - z0) * (z3 - z0);
		s13 = (x1 - x0) * (x3 - x0) + (y1 - y0) * (y3 - y0) + (z1 - z0) * (z3 - z0);

		sum1 = s11 + s22 + s33;
		sum2 = s12 + s23 + s13;

		anum = 3. * pow(deta * riza2, 2. / 3.);
		den = ((3. / 2.) * sum1) - sum2;

		qual[ip] = anum / den;
		qualityFile << qual[ip] << " " << ip << endl;

		if (qual[ip] < 1.e-10)
		{
			kountInv++;
		}
	}
	qMin = 100.;
	qMax = -100.;
	qMean = 0.;

	for (int ip = 1; ip <= ntet; ip++)
	{
		qLoc = qual[ip];
		qMean = qMean + qLoc;
		qMin = min(qMin, qLoc);
		qMax = max(qMax, qLoc);
	}
	qMean = qMean / double(ntet);

	double qSTD = 0.;

	for (int ip = 1; ip <= ntet; ip++)
	{
		qSTD += pow(qual[ip] - qMean, 2.);
	}
	qSTD = sqrt(qSTD / double(ntet));

	cout << "Quality Metrics  :: SHAPE " << endl;
	cout << "---------------- " << endl;
	cout << setprecision(10) << "Min,  Max   :" << qMin << "  " << qMax << endl;
	cout << setprecision(10) << "Mean,  STD  :" << qMean << "  " << qSTD << endl;
	cout << "  ------  ------  ------  ------  ------ " << endl;
	cout << " Inverted   :" << kountInv << endl;
	return;
}

void qualityCheck::Jacobian(double **cl, int ntet, string eleFile)
{
	string line;
	int is0, is1, is2, is3;
	int nod[4];
	double x0, x1, x2, x3, y0, y1, y2, y3, z0, z1, z2, z3,
		deta, qMin, qMax, qMean, qLoc, den;
	vector<double> qual;
	qual.resize(ntet + 1);

	eleFile += ".ele";
	ifstream elefile;
	elefile.open("input_data/" + eleFile);

	if (elefile.fail())
	{
		cout << " Non-existent .ele file " << eleFile << endl;
		exit;
	}
	getline(elefile, line);
	int kountInv = 0; // number of inverted tetrahedrons
	int kount = 0;
	double riza2 = sqrt(2.);

	ofstream qualityFile;
	qualityFile.open("output_data/Mesh_Quality");
	qualityFile << endl
				<< endl;

	for (int ip = 1; ip <= ntet; ip++)
	{
		elefile >> is0 >> is1 >> is2 >> is3;
		x0 = cl[0][is0 - 1];
		y0 = cl[1][is0 - 1];
		z0 = cl[2][is0 - 1];
		x1 = cl[0][is1 - 1];
		y1 = cl[1][is1 - 1];
		z1 = cl[2][is1 - 1];
		x2 = cl[0][is2 - 1];
		y2 = cl[1][is2 - 1];
		z2 = cl[2][is2 - 1];
		x3 = cl[0][is3 - 1];
		y3 = cl[1][is3 - 1];
		z3 = cl[2][is3 - 1];

		deta = (x1 - x0) * ((y2 - y0) * (z3 - z0) - (y3 - y0) * (z2 - z0)) -
			   (x2 - x0) * ((y1 - y0) * (z3 - z0) - (y3 - y0) * (z1 - z0)) +
			   (x3 - x0) * ((y1 - y0) * (z2 - z0) - (y2 - y0) * (z1 - z0));

		qual[ip] = deta;
	}

	qMin = 100.;
	qMax = -100.;
	qMean = 0.;

	for (int ip = 1; ip <= ntet; ip++)
	{
		qLoc = qual[ip];
		qMean = qMean + qLoc;
		qMin = min(qMin, qLoc);
		qMax = max(qMax, qLoc);
	}
	qMean = qMean / double(ntet);

	double qSTD = 0.;

	for (int ip = 1; ip <= ntet; ip++)
	{
		qSTD += pow(qual[ip] - qMean, 2.);
	}
	qSTD = sqrt(qSTD / double(ntet));

	cout << "Quality Metrics  :: JACOBIAN " << endl;
	cout << "---------------- " << endl;
	cout << setprecision(10) << "Min,  Max   :" << qMin << "  " << qMax << endl;
	cout << setprecision(10) << "Mean,  STD  :" << qMean << "  " << qSTD << endl;
	cout << "  ------  ------  ------  ------  ------ " << endl;
}

void qualityCheck::aspectRatio(double **cl, int ntet, string eleFile)
{
	string line;
	int is0, is1, is2, is3;
	int nod[4];
	double x0, x1, x2, x3, y0, y1, y2, y3, z0, z1, z2, z3,
		deta, s11, s22, s33, s12, s23, s13, sum1, sum2, anum,
		qMin, qMax, qMean, qLoc, den;
	vector<double> qual;
	double ab[3], ac[3], ad[3], bc[3], bd[3], cd[3],
		Nabc[3], Nabd[3], Nacd[3], Nbcd[3], a, hab, hac, had,
		hbc, hbd, hcd, hmax, r, NABC, NABD, NACD, NBCD, H[6],
		*maxVal;

	qual.resize(ntet + 1);

	eleFile += ".ele";
	ifstream elefile;

	elefile.open("input_data/" + eleFile);

	if (elefile.fail())
	{
		cout << " Non-existent .ele file " << eleFile << endl;
		exit;
	}
	getline(elefile, line);
	int kountInv = 0; // number of inverted tetrahedrons
	int kount = 0;
	double riza2 = sqrt(2.);

	ofstream qualityFile;
	qualityFile.open("output_data/aspect_ratio_tet");
	qualityFile << endl
				<< endl;

	for (int ip = 1; ip <= ntet; ip++)
	{
		elefile >> is0 >> is1 >> is2 >> is3;
		x0 = cl[0][is0 - 1];
		y0 = cl[1][is0 - 1];
		z0 = cl[2][is0 - 1];
		x1 = cl[0][is1 - 1];
		y1 = cl[1][is1 - 1];
		z1 = cl[2][is1 - 1];
		x2 = cl[0][is2 - 1];
		y2 = cl[1][is2 - 1];
		z2 = cl[2][is2 - 1];
		x3 = cl[0][is3 - 1];
		y3 = cl[1][is3 - 1];
		z3 = cl[2][is3 - 1];

		// VECTORS OF TETRAHEDRON
		ab[0] = x1 - x0;
		ab[1] = y1 - y0;
		ab[2] = z1 - z0;
		ac[0] = x2 - x0;
		ac[1] = y2 - y0;
		ac[2] = z2 - z0;
		ad[0] = x3 - x0;
		ad[1] = y3 - y0;
		ad[2] = z3 - z0;
		bc[0] = x2 - x1;
		bc[1] = y2 - y1;
		bc[2] = z2 - z1;
		cd[0] = x3 - x2;
		cd[1] = y3 - y2;
		cd[2] = z3 - z2;
		bd[0] = x3 - x1;
		bd[1] = y3 - y1;
		bd[2] = z3 - z1;

		a = ab[0] * (ac[1] * ad[2] - ac[2] * ad[1]) + ab[1] * (ac[2] * ad[0] - ac[0] * ad[2]) +
			ab[2] * (ac[0] * ad[1] - ac[1] * ad[0]);

		// NORMALS TO FACES
		Nabc[0] = (ab[1] * ac[2] - ab[2] * ac[1]);
		Nabc[1] = (ab[2] * ac[0] - ab[0] * ac[2]);
		Nabc[2] = (ab[0] * ac[1] - ab[1] * ac[0]);

		NABC = sqrt(Nabc[0] * Nabc[0] + Nabc[1] * Nabc[1] + Nabc[2] * Nabc[2]);

		Nabd[0] = (ab[1] * ad[2] - ab[2] * ad[1]);
		Nabd[1] = (ab[2] * ad[0] - ab[0] * ad[2]);
		Nabd[2] = (ab[0] * ad[1] - ab[1] * ad[0]);

		NABD = sqrt(Nabd[0] * Nabd[0] + Nabd[1] * Nabd[1] + Nabd[2] * Nabd[2]);

		Nacd[0] = (ac[1] * ad[2] - ac[2] * ad[1]);
		Nacd[1] = (ac[2] * ad[0] - ac[0] * ad[2]);
		Nacd[2] = (ac[0] * ad[1] - ac[1] * ad[0]);

		NACD = sqrt(Nacd[0] * Nacd[0] + Nacd[1] * Nacd[1] + Nacd[2] * Nacd[2]);

		Nbcd[0] = (bc[1] * bd[2] - bc[2] * bd[1]);
		Nbcd[1] = (bc[2] * bd[0] - bc[0] * bd[2]);
		Nbcd[2] = (bc[0] * bd[1] - bc[1] * bd[0]);

		NBCD = sqrt(Nbcd[0] * Nbcd[0] + Nbcd[1] * Nbcd[1] + Nbcd[2] * Nbcd[2]);

		// LENGTHS of EDGES
		hab = sqrt(ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2]);
		hac = sqrt(ac[0] * ac[0] + ac[1] * ac[1] + ac[2] * ac[2]);
		had = sqrt(ad[0] * ad[0] + ad[1] * ad[1] + ad[2] * ad[2]);
		hbc = sqrt(bc[0] * bc[0] + bc[1] * bc[1] + bc[2] * bc[2]);
		hbd = sqrt(bd[0] * bd[0] + bd[1] * bd[1] + bd[2] * bd[2]);
		hcd = sqrt(cd[0] * cd[0] + cd[1] * cd[1] + cd[2] * cd[2]);

		H[0] = hab;
		H[1] = hac;
		H[2] = had;
		H[3] = hbc;
		H[4] = hbd;
		H[5] = hcd;

		maxVal = max_element(H, H + 6);
		hmax = *maxVal;

		r = abs(a) / (NABC + NABD + NACD + NBCD);

		qual[ip] = hmax / (2 * sqrt(6) * r);

		qualityFile << qual[ip] << " " << ip << endl;
	}

	qMin = 100.;
	qMax = -100.;
	qMean = 0.;

	for (int ip = 1; ip <= ntet; ip++)
	{
		qLoc = qual[ip];
		qMean = qMean + qLoc;
		qMin = min(qMin, qLoc);
		qMax = max(qMax, qLoc);
	}
	qMean = qMean / double(ntet);

	double qSTD = 0.;

	for (int ip = 1; ip <= ntet; ip++)
	{
		qSTD += pow(qual[ip] - qMean, 2.);
	}
	qSTD = sqrt(qSTD / double(ntet));

	cout << "Quality Metrics  :: ASPECT RATIO " << endl;
	cout << "---------------- " << endl;
	cout << setprecision(10) << "Min,  Max   :" << qMin << "  " << qMax << endl;
	cout << setprecision(10) << "Mean,  STD  :" << qMean << "  " << qSTD << endl;
	cout << "  ------  ------  ------  ------  ------ " << endl;
}

qualityCheck::~qualityCheck(){};