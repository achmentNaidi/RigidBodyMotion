#include "../include/geomHandle.h"
#include "../include/DataStructures.h"

geomHandle::geomHandle(double **coor_, vector<int> logfr_, int ns_)
    : coor(coor_), logfr(logfr_), ns(ns_) {}

void geomHandle::rotate2D(vector<double> &x, vector<double> &y, double cx, double cy, double angle)
{
    double xlocal, ylocal;
    angle = angle * 3.141592653589 / 180.;

    for (int i = 0; i < x.size(); i++)
    {
        xlocal = x[i] - cx;
        ylocal = y[i] - cy;
        x[i] = xlocal * cos(angle) - ylocal * sin(angle) + cx;
        y[i] = xlocal * sin(angle) + ylocal * cos(angle) + cy;
    }

    ofstream rotatedGeo;
    rotatedGeo.open("RotatedGeo.dat");

    rotatedGeo << cx << "   " << cy << endl;

    for (int i = 0; i < x.size(); i++)
    {
        rotatedGeo << x[i] << "  " << y[i] << endl;
    }
    rotatedGeo.close();
}

void geomHandle::read2dGeometry(string fileName, vector<double> &x_g, vector<double> &y_g)
{
    int id;
    string line;
    vector<double> p;
    double xi, yi, cx = 0, cy = 0;
    ifstream geometry;
    geometry.open(fileName);

    if (geometry.is_open())
    {
        id = 0;
        while (getline(geometry, line))
        {
            cout << line << '\n';

            geometry >> xi >> yi;

            x_g.push_back(xi);
            y_g.push_back(yi);
            id++;
        }

        geometry.close();
    }
    else
        cout << "Unable to open file";
}

centroid geomHandle::centerOfGravity(vector<double> x_g, vector<double> y_g)
{
    centroid p;
    p.x = 0;
    p.y = 0;
    double area = 0;
    for (int i = 0; i < x_g.size() - 1; i++)
    {
        area += 0.5 * (x_g[i] * y_g[i + 1] - x_g[i + 1] * y_g[i]);
    }

    for (int i = 0; i < x_g.size() - 1; i++)
    {
        p.x += (x_g[i] + x_g[i + 1]) * (x_g[i] * y_g[i + 1] - x_g[i + 1] * y_g[i]) / 6 / area;
        p.y += (y_g[i] + y_g[i + 1]) * (x_g[i] * y_g[i + 1] - x_g[i + 1] * y_g[i]) / 6 / area;
    }
    return p;
}

double **geomHandle::wingBending(double alpha)
{
    // alpha = bending coefficient
    double **coorp = matrix<double>(3, ns);
    for (int is = 1; is <= ns; is++)
    {
        if (logfr[is] == 3 || logfr[is] == 2)
        {
            coorp[0][is - 1] = coor[0][is - 1];
            coorp[1][is - 1] = coor[1][is - 1] + alpha * coor[2][is - 1] * coor[2][is - 1];
            coorp[2][is - 1] = coor[2][is - 1];
        }
        else
        {
            coorp[0][is - 1] = coor[0][is - 1];
            coorp[1][is - 1] = coor[1][is - 1];
            coorp[2][is - 1] = coor[2][is - 1];
        }
    }
    return coorp;
}

double **geomHandle::wingTorsionBending(double alpha)
{

    /*
          Wing Bending & Torsion
          ---------------------
    */

    // alpha = bending coefficient
    double **coorp = matrix<double>(3, ns);
    double phi, xr, yr = 5.;
    for (int is = 1; is <= ns; is++)
    {
        if (logfr[is] == 3 || logfr[is] == 2)
        {
            coorp[0][is - 1] = coor[0][is - 1];
            coorp[1][is - 1] = coor[1][is - 1] + alpha * coor[2][is - 1] * coor[2][is - 1];
            coorp[2][is - 1] = coor[2][is - 1];

            xr = 5.25 + 0.45 * coor[2][is - 1];
            phi = alpha * pow(coor[2][is - 1], 2);

            coorp[0][is - 1] = (coor[0][is - 1] - xr) * cos(phi) -
                               (coor[1][is - 1] - yr) * sin(phi) + xr;
            coorp[1][is - 1] = (coor[1][is - 1] - yr) * cos(phi) -
                               (coor[0][is - 1] - xr) * sin(phi) + yr + phi;
            coorp[2][is - 1] = coor[2][is - 1];
        }
        else
        {
            coorp[0][is - 1] = coor[0][is - 1];
            coorp[1][is - 1] = coor[1][is - 1];
            coorp[2][is - 1] = coor[2][is - 1];
        }
    }
    return coorp;
}

void geomHandle::getGeometry()
{
    ofstream geometry;
    geometry.open("geom2D.dat");

    for (int i = 0; i < ns; i++)
    {
        if (logfr[i] == 3)
        {
            geometry << coor[0][i] << "	" << coor[1][i] << endl;
        }
    }
}

void geomHandle::overwriteDefFile(vector<double> xg, vector<double> yg, string fileName)
{
    fileName = fileName + ".def";
    ifstream indef;
    indef.open(fileName);

    int nodes, counter;
    vector<double> x, y;
    vector<int> logfr;
    indef >> nodes;

    logfr.resize(nodes);
    x.resize(nodes);
    y.resize(nodes);
    for (int i = 0; i < nodes; i++)
    {
        indef >> logfr[i];
    }

    for (int i = 0; i < nodes; i++)
    {
        indef >> x[i];
    }

    for (int i = 0; i < nodes; i++)
    {
        indef >> y[i];
    }

    indef.close();

    ofstream outdef;
    outdef.open(fileName);
    counter = 0;
    for (int i = 0; i < nodes; i++)
    {

        if (logfr[i] == 3)
        {
            x[i] = xg[counter];
            y[i] = yg[counter];
            counter++;
        }
    }

    outdef << nodes << endl;

    for (int i = 0; i < nodes; i++)
    {
        outdef << logfr[i] << endl;
        ;
    }
    for (int i = 0; i < nodes; i++)
    {
        outdef << x[i] << endl;
        ;
    }

    for (int i = 0; i < nodes; i++)
    {
        outdef << y[i] << endl;
    }

    outdef.close();
}

geomHandle::~geomHandle() {}