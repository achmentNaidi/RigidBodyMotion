#ifndef __DATASTRUCTURES_H__
#define __DATASTRUCTURES_H__

#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>

#include "qualityCheck.h"

using namespace std;

struct gridInfo2D_Structured {
    int Imax = 0, Jmax = 0, bountaryStart = 0, bountaryEnd = 0;
    int totalNodes;
    int **idMatrix = 0;
    vector<double> x0, y0;
    vector<double> xf, yf;
};

//  STAY HERE
struct gridInfo {
    vector<int> nu, logfr;
    double **coor, **coorp, pitch;
    int nall, np, nq, ntet, npyr, npri, nhex;
    int ns, nper;
};

class DataStructures {
   public:
    DataStructures();
    DataStructures(int ns, int np, int nq, double **coor_, vector<int> &logfr_, vector<int> &nu);
    DataStructures(vector<int> &logfr_, int ns_, double **coor_, int ntet,
                   int npyr, int npri, int nhex, int nall, vector<int> &nu);
    void Create2D();

    void Create3D();
    void setPeriodicityConstants(int mpar, int kaxial, int isperiph);
    vector<int> get_ndeg() { return ndeg; }
    vector<int> get_jaret() { return jaret; }
    int get_nbseg() { return nbseg; }
    vector<int> get_logfr() { return logfr; }
    double get_pitch() { return pitch; }
    int **get_iper() { return iper; }
    int get_nper() { return nper; }

    int mpar, kaxial, isperiph, nper;

    ~DataStructures();

   private:
    //================//
    // INTEGER PARAMETERS

    const int maxseg = 5000000, nvmaxall = 6000000, maxbfac = 100000,
              maxsegpe = 100000, maxlist = 10, maxnod = 50000, maxpyr = 4 * maxnod,
              maxnu = 300000, maxbseg = 5000;

    // REAL ALLOCATABLE ARRAYS

    double **coor, **vnocl, **vnofac;
    vector<double> cell, vol;

    // INTEGER ALLOCATABLE ARRAYS

    vector<int> logfr, nu, nusg, ndeg, listn, listbf;
    int **iper, **ipersg;

    // INTEGER ARRAYS

    int **nubo, **nubf;
    vector<int> listnp1, listnp2, ibfc2te,
        listbf1, listbf2, ibsg2te, jaret, ibsg2tr;

    // REAL VARIABLES

    double pint, pext, pitch;

    // INTEGER VARIABLES

    int np, nq, nbseg, nvseg;  // 2D ONLY
    int ns, nseg, ntet, npyr, npri, nhex, nall, nbfac;
    int npersg;
    int max_str, nvmaxall_out;

    // PRIVATE FUNCTIONS
    void resizeVectors(int maxseg, int nvmxall, int maxbfac, int maxsegpe, int maxlist, int nall);
    void vecProd(double v1x, double v1y, double v1z, double v2x,
                 double v2y, double v2z, double &xre, double &yre, double &zre);
    void filistn();
    void analnode(int kkk0, int &mposa, int &m100, int &m010, int &m001);
    void setperio();
    void setperiph();
    void numfaces();
    void volumes();
    void filistbf();
    void perfac();
    int analface(int lgf10, int lgf20, int lgf30, int lgf40, int ktipos);
    void calvnocl();
    void numsegs3D();
    void numsegs2D();
    void perseg();
    void fjaret3D();
    void fjaret2D();
    void fjaret_el();
    void matrix_mult(int nvar, double **xx, double **yy, double **zz, int mvar1, int mvar2, int mvar3);
    void rotation_matrix(int nvar, int iop, double cm, double sm, double **rra, double **rrb);
    void virtualPeriodicNeighbours();
};

template <class myType>
myType **matrix(int rows, int columns);
#endif  // __DATASTRUCTURES_H__