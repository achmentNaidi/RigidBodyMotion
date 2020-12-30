#include "geomHandle.h"

#include "DataStructures.h"

geomHandle::geomHandle(double **coor_, vector<int> logfr_, int ns_)
    : coor(coor_), logfr(logfr_), ns(ns_) {
    coorp = matrix<double>(3, ns);
    for (int i = 0; i < ns; i++) {
        coorp[0][i] = coor[0][i];
        coorp[1][i] = coor[1][i];
        coorp[2][i] = coor[2][i];
    }
}

void geomHandle::rotate2D(vector<double> &x, vector<double> &y, double cx, double cy, double angle) {
    double xlocal, ylocal;
    angle = angle * 3.141592653589 / 180.;

    for (int i = 0; i < x.size(); i++) {
        xlocal = x[i] - cx;
        ylocal = y[i] - cy;
        x[i] = xlocal * cos(angle) - ylocal * sin(angle) + cx;
        y[i] = xlocal * sin(angle) + ylocal * cos(angle) + cy;
    }

    ofstream rotatedGeo;
    rotatedGeo.open("RotatedGeo.dat");

    rotatedGeo << cx << "   " << cy << endl;

    for (int i = 0; i < x.size(); i++) {
        rotatedGeo << x[i] << "  " << y[i] << endl;
    }
    rotatedGeo.close();
}

void geomHandle::read2dGeometry(string fileName, vector<double> &x_g, vector<double> &y_g) {
    int id;
    string line;
    vector<double> p;
    double xi, yi, cx = 0, cy = 0;
    ifstream geometry;
    geometry.open(fileName);

    if (geometry.is_open()) {
        id = 0;
        while (getline(geometry, line)) {
            cout << line << '\n';

            geometry >> xi >> yi;

            x_g.push_back(xi);
            y_g.push_back(yi);
            id++;
        }

        geometry.close();
    } else
        cout << "Unable to open file";
}

centroid geomHandle::centerOfGravity(vector<double> x_g, vector<double> y_g) {
    centroid p;
    p.x = 0;
    p.y = 0;
    double area = 0;
    for (int i = 0; i < x_g.size() - 1; i++) {
        area += 0.5 * (x_g[i] * y_g[i + 1] - x_g[i + 1] * y_g[i]);
    }

    for (int i = 0; i < x_g.size() - 1; i++) {
        p.x += (x_g[i] + x_g[i + 1]) * (x_g[i] * y_g[i + 1] - x_g[i + 1] * y_g[i]) / 6 / area;
        p.y += (y_g[i] + y_g[i + 1]) * (x_g[i] * y_g[i + 1] - x_g[i + 1] * y_g[i]) / 6 / area;
    }
    return p;
}

void geomHandle::wingBending(double alpha) {
    // alpha = bending coefficient
    for (int is = 1; is <= ns; is++) {
        if (logfr[is] == 3 || logfr[is] == 2) {
            coorp[0][is - 1] = coorp[0][is - 1];
            coorp[1][is - 1] = coorp[1][is - 1] + alpha * coorp[2][is - 1] * coorp[2][is - 1];
            coorp[2][is - 1] = coorp[2][is - 1];
        } else {
            coorp[0][is - 1] = coorp[0][is - 1];
            coorp[1][is - 1] = coorp[1][is - 1];
            coorp[2][is - 1] = coorp[2][is - 1];
        }
    }
}

void geomHandle::wingTorsionBending(double alpha) {
    /*
          Wing Bending & Torsion
          ---------------------
    */

    // alpha = bending coefficient
    double phi, xr, yr = 5.;
    for (int is = 1; is <= ns; is++) {
        if (logfr[is] == 3 || logfr[is] == 2) {
            coorp[0][is - 1] = coorp[0][is - 1];
            coorp[1][is - 1] = coorp[1][is - 1] + alpha * coorp[2][is - 1] * coorp[2][is - 1];
            coorp[2][is - 1] = coorp[2][is - 1];

            xr = 5.25 + 0.45 * coorp[2][is - 1];
            phi = alpha * pow(coorp[2][is - 1], 2);

            coorp[0][is - 1] = (coorp[0][is - 1] - xr) * cos(phi) -
                               (coorp[1][is - 1] - yr) * sin(phi) + xr;
            coorp[1][is - 1] = (coorp[1][is - 1] - yr) * cos(phi) -
                               (coorp[0][is - 1] - xr) * sin(phi) + yr + phi;
            coorp[2][is - 1] = coorp[2][is - 1];
        } else {
            coorp[0][is - 1] = coorp[0][is - 1];
            coorp[1][is - 1] = coorp[1][is - 1];
            coorp[2][is - 1] = coorp[2][is - 1];
        }
    }
}

void geomHandle::translateGeometry(double dx, double dy, double dz) {
    double phi, xr, yr = 5.;
    for (int is = 1; is <= ns; is++) {
        if (logfr[is] == 3 || logfr[is] == 2) {
            coorp[0][is - 1] = coorp[0][is - 1] + dx;
            coorp[1][is - 1] = coorp[1][is - 1] + dy;
            coorp[2][is - 1] = coorp[2][is - 1] + dz;
        } else {
            coorp[0][is - 1] = coorp[0][is - 1];
            coorp[1][is - 1] = coorp[1][is - 1];
            coorp[2][is - 1] = coorp[2][is - 1];
        }
    }
}

void geomHandle::rotateGeometry(double u, double v, double w) {
    // Rotation point X,Y,Z

    double **R = matrix<double>(3, 3);
    u = u * 3.141592653589 / 180.;
    v = v * 3.141592653589 / 180.;
    w = w * 3.141592653589 / 180.;

    double ccu = cos(u);
    double ssu = sin(u);
    double ccv = cos(v);
    double ssv = sin(v);
    double ccw = cos(w);
    double ssw = sin(w);

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

    for (int inei = 1; inei <= ns; inei++) {
        if (logfr[inei] == 3) {
            coorp[0][inei - 1] = R[0][0] * coor[0][inei - 1] + R[0][1] * coor[1][inei - 1] + R[0][2] * coor[2][inei - 1];
            coorp[1][inei - 1] = R[1][0] * coor[0][inei - 1] + R[1][1] * coor[1][inei - 1] + R[1][2] * coor[2][inei - 1];
            coorp[2][inei - 1] = R[2][0] * coor[0][inei - 1] + R[2][1] * coor[1][inei - 1] + R[2][2] * coor[2][inei - 1];
        } else {
            coorp[0][inei - 1] = coor[0][inei - 1];
            coorp[1][inei - 1] = coor[1][inei - 1];
            coorp[2][inei - 1] = coor[2][inei - 1];
        }
    }
}

void geomHandle::getGeometry() {
    ofstream geometry;
    geometry.open("defShape.dat");

    for (int i = 0; i < ns; i++) {
        if (logfr[i] == 3) {
            geometry << coor[0][i] << "	" << coor[1][i] << "\t" << coor[2][i] << endl;
        }
    }
    geometry.close();
}

void geomHandle::overwriteDefFile(vector<double> xg, vector<double> yg, string fileName) {
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
    for (int i = 0; i < nodes; i++) {
        indef >> logfr[i];
    }

    for (int i = 0; i < nodes; i++) {
        indef >> x[i];
    }

    for (int i = 0; i < nodes; i++) {
        indef >> y[i];
    }

    indef.close();

    ofstream outdef;
    outdef.open(fileName);
    counter = 0;
    for (int i = 0; i < nodes; i++) {
        if (logfr[i] == 3) {
            x[i] = xg[counter];
            y[i] = yg[counter];
            counter++;
        }
    }

    outdef << nodes << endl;

    for (int i = 0; i < nodes; i++) {
        outdef << logfr[i] << endl;
        ;
    }
    for (int i = 0; i < nodes; i++) {
        outdef << x[i] << endl;
        ;
    }

    for (int i = 0; i < nodes; i++) {
        outdef << y[i] << endl;
    }

    outdef.close();
}

geomHandle::~geomHandle() {}