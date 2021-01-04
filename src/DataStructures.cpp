#include "DataStructures.h"

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

DataStructures::DataStructures() {
}

//===================2D===========================//

DataStructures::DataStructures(int ns_, int np_, int nq_, double **coor_, vector<int> &logfr_, vector<int> &nu_)
    : ns(ns_), np(np_), nq(nq_), nu(nu_), logfr(logfr_), coor(coor_) {
    listnp1.resize(maxlist + 1);
    listnp2.resize(maxlist + 1);
}

void DataStructures::Create2D() {
    if (np + nq > maxpyr) {
        cout << "geom2d: Increase maxpyr at DataStructures.h" << endl;
        exit(0);
    }
    if (ns > maxnod) {
        cout << "geom2d: Increase maxnod at DataStructures.h" << endl;
    }

    //	--------------------------------------------------
    //	Numerate the segments(nubo, nusg, ibsg2tr)
    //	 -------------------------------------------------
    if (isPeriodic) {
        if (isperiph == 1) {
            setperiph();  //  peripheral cascade
            cout << " setperiph completed...." << endl;
        } else {
            setperio();  //  linear cascade
            cout << " setperio completed...." << endl;
        }
        numsegs2D();  // provisionally uses jaret, ndeg
        cout << "numsegs completed...." << endl;
        fjaret2D();
        cout << "fjaret completed...." << endl;
        virtualPeriodicNeighbours();
    } else {
        numsegs2D();  // provisionally uses jaret, ndeg
        cout << "numsegs completed...." << endl;
        fjaret2D();
        cout << "fjaret completed...." << endl;
    }

    //	  ------
    //	  Ending
    //	  ------
    cout << " Number of nodes                : " << ns << endl;
    cout << " Number of segments			 : " << nseg << endl;
    cout << " Number of triangles			 : " << np << endl;
    cout << " Number of quadrangles			 : " << nq << endl;
    cout << " Number of boundary segments	 : " << nbseg << endl;

    cout << " Geom2d completed ..." << endl;
}

void DataStructures::numsegs2D() {
    int is, kpoi_1, is1, is2, iprmid, kount, k1, k2,
        is3, is4, ivseg1, ivseg2;
    int ind = 0, ind2 = 0;
    int **iquadseg;
    iquadseg = matrix<int>(4, 2);  //side enumeration for quadrangles
    iquadseg[0][0] = 2;
    iquadseg[1][0] = 3;
    iquadseg[2][0] = 4;
    iquadseg[3][0] = 1;
    iquadseg[0][1] = 3;
    iquadseg[1][1] = 4;
    iquadseg[2][1] = 1;
    iquadseg[3][1] = 2;
    vector<int> nod;
    nod.resize(3);
    ndeg.resize(ns + 1);
    nusg.resize(maxnu);
    jaret.resize(nvmaxall);
    ibsg2tr.resize(maxbseg);
    int **iauxn = matrix<int>(3, 2);  //opposite nodes of segments for triangles
    iauxn[0][0] = 2;
    iauxn[1][0] = 1;
    iauxn[2][0] = 1;
    iauxn[0][1] = 3;
    iauxn[1][1] = 3;
    iauxn[2][1] = 2;

    nubo = matrix<int>(2, maxseg);
    nseg = 0;
    nbseg = 0;
    //
    //	Only in this subroutine nu turns -ve and
    //	jaret is provisionally the triangles around a node
    for (int i = 0; i < nvmaxall; i++)
        jaret[i] = 0;
    //
    //  ----------------------------------------------------------------
    //	Find the elements around a node (ATTENTION: in bounds,#elements
    //									 around a node is different from
    //									 #segments around a node!!)
    //	----------------------------------------------------------------
    for (int i = 1; i <= np; i++) {
        for (int j = 0; j < 3; j++) {
            ind += 1;
            is = nu[ind];
            ndeg[is] = ndeg[is] + 1;
        }
    }
    for (int i = np + 1; i <= np + nq; i++) {
        for (int j = 0; j < 4; j++) {
            ind += 1;
            is = nu[ind];
            ndeg[is] = ndeg[is] + 1;
        }
    }
    //
    for (int k = 1; k <= ns; k++) {
        if (ndeg[k] == 0) {
            cout << endl;
            cout << "ERROR: hanging node, index: " << k << endl;
            exit(0);
        }
        ndeg[k] = ndeg[k] + ndeg[k - 1];  // build index
        jaret[ndeg[k]] = ndeg[k - 1];     // provisory
    }
    if (2 * ndeg[ns] > nvmaxall) {
        cout << "numsegs:Increase nvmaxall at DataStructures.h-->" << 2 * ndeg[ns] << "++" << endl;
        exit(0);
    }
    //
    ind = 0;
    for (int i = 1; i <= np; i++) {
        for (int j = 1; j <= 3; j++) {
            ind++;
            is = nu[ind];
            jaret[ndeg[is]] = jaret[ndeg[is]] + 1;
            kpoi_1 = jaret[ndeg[is]];
            jaret[kpoi_1] = i;  // Element
        }
    }
    for (int i = np + 1; i <= np + nq; i++) {
        for (int j = 1; j <= 4; j++) {
            ind++;
            is = nu[ind];
            jaret[ndeg[is]] = jaret[ndeg[is]] + 1;
            kpoi_1 = jaret[ndeg[is]];
            jaret[kpoi_1] = i;  // Element
        }
    }

    //
    //	-----------------------------------
    //	Numerate the segments and find NUBO
    //	-----------------------------------
    ind = 0;
    //std::vector<int>::iterator it;
    //it = nu.begin();
    //nu.insert(it,0);
    for (int ip = 1; ip <= np; ip++) {
        for (int m = 0; m < 3; m++) {
            //
            //			Each segment is identidief by its opposite node, in the triangle:
            ind++;
            if (nu[ind] < 0) {
                continue;
            }  // this segment has been treated before
            //
            //			Else, this is a NEW segment, then:
            ind2 = (ip - 1) * 3;
            is1 = abs(nu[ind2 + iauxn[m][0]]);
            is2 = abs(nu[ind2 + iauxn[m][1]]);

            if (nseg + 2 * nq > maxseg) {
                cout << "numsegs: Increase maxseg at DataStructures.h" << endl;
                exit(0);
            }
            if (nbseg > maxbseg) {
                cout << "numsegs: Increase maxbseg at DataStructures.h" << endl;
                exit(0);
            }
            ibsg2tr[nbseg] = nseg;
            nubo[0][nseg] = is1;
            nubo[1][nseg] = is2;

            nusg[ind] = nseg;
            nu[ind] = -nu[ind];  // attention, nu turn -ve
            //
            //          Find adjacent segments to avoid duplicating.
            //			Also mark the segments that lie over the boundary
            for (int ik = ndeg[is1 - 1] + 1; ik <= ndeg[is1]; ik++) {
                iprmid = jaret[ik];
                //
                if (iprmid == ip)
                    continue;
                if (iprmid <= np) {
                    kount = 0;
                    ind2 = (iprmid - 1) * 3;
                    for (int m1 = 0; m1 < 3; m1++) {
                        k1 = abs(nu[ind2 + m1 + 1]);
                        k2 = 1;
                        if (k1 == is1 || k1 == is2) {
                            k2 = 0;
                        }
                        if (k2 != 0) {
                            kount++;
                            nod[0] = m1 + 1;
                        }
                    }
                    if (kount == 2) {
                        continue;
                    }
                    if (kount == 0) {
                        cout << "numsegs Error, with data !" << endl;
                        exit(0);
                    }

                    nu[ind2 + nod[0]] = -nu[ind2 + nod[0]];  // for this also,nu turns-ve
                    int check;
                    check = nu[ind2 + nod[0]];
                    nusg[ind2 + nod[0]] = nseg;
                    nbseg = nbseg - 1;  // This was a boundary segment
                    break;
                } else {
                    kount = 0;
                    ind2 = np * 3 + (iprmid - np - 1) * 4;
                    for (int m1 = 1; m1 <= 4; m1++) {
                        k1 = abs(nu[ind2 + m1]);
                        k2 = 1;
                        if (k1 == is1 || k1 == is2) {
                            k2 = 0;
                        }
                        if (k2 != 0) {
                            nod[kount] = m1;
                            kount++;
                        }
                    }
                    if (kount == 3) {
                        continue;
                    }
                    if (kount <= 1) {
                        cout << "numsegs: Error, with data !" << endl;
                        exit(0);
                    }
                    if (nod[0] == 1 && nod[1] == 4) {
                        nod[1] = 1;
                    }
                    nu[ind2 + nod[1]] = -nu[ind2 + nod[1]];
                    nusg[ind2 + nod[1]] = nseg;
                    nbseg = nbseg - 1;  // This was not a boundary segment
                    break;
                }
            }
            nseg++;
            nbseg++;
        }
    }
    //
    for (int ip = np + 1; ip <= np + nq; ip++)  // Loop on quandrangle
    {
        for (int m = 0; m < 4; m++)  //Loop on the four nodes of each quadrangle
        {
            ind++;
            //			Each segment is identified according to the enumeration in iquadseg:
            if (nu[ind] < 0) {
                continue;
            }
            //			Else, this is a NEW segment, then:
            ind2 = np * 3 + (ip - np - 1) * 4;
            is1 = abs(nu[ind2 + iquadseg[m][0]]);
            is2 = abs(nu[ind2 + iquadseg[m][1]]);
            if (nseg + 2 * nq > maxseg) {
                cout << "numsegs: Increase maxseg at DataStructures.h " << endl;  // EULER
                exit(0);
            }
            if (nbseg > maxbseg) {
                cout << "numsegs: Increase maxbseg at DataStructures.h " << endl;
            }
            ibsg2tr[nbseg] = nseg;
            nubo[0][nseg] = is1;
            nubo[1][nseg] = is2;

            nusg[ind] = nseg;
            nu[ind] = -nu[ind];  // attention, nu turns -ve
            //
            //			Find adjacent segment to avoid duplicating.
            //			Also mark the segments that lie over the boundary.
            for (int ik = ndeg[is1 - 1] + 1; ik <= ndeg[is1]; ik++) {
                iprmid = jaret[ik];
                //
                if (iprmid == ip) {
                    continue;
                }  // skip the quadrangle at hand
                if (iprmid <= np) {
                    kount = 0;
                    ind2 = (iprmid - 1) * 3;
                    for (int m1 = 1; m1 <= 3; m1++)  // Loop on the 3 nodes of the triangle
                    {
                        k1 = abs(nu[ind2 + m1]);
                        k2 = 1;
                        if (k1 == is1 || k1 == is2) {
                            k2 = 0;
                        }
                        if (k2 != 0) {
                            kount++;
                            nod[0] = m1;
                        }
                    }
                    if (kount == 2) {
                        continue;
                    }
                    if (kount == 0) {
                        cout << "numsegs: Error, with data !" << endl;
                        exit(0);
                    }
                    nu[ind2 + nod[0]] = -nu[ind2 + nod[0]];  // for this also, nu turns-ve
                    nusg[ind2 + nod[0]] = nseg;
                    nbseg = nbseg - 1;  // This was not a boundary segment
                    break;
                } else {
                    kount = 0;
                    ind2 = np * 3 + (iprmid - np - 1) * 4;
                    for (int m1 = 1; m1 <= 4; m1++)  // Loop on the 4 nodes of the quadrangle
                    {
                        k1 = abs(nu[ind2 + m1]);
                        k2 = 1;
                        if (k1 == is1 || k1 == is2) {
                            k2 = 0;
                        }
                        if (k2 != 0) {
                            nod[kount] = m1;
                            kount++;
                        }
                    }
                    if (kount == 3) {
                        continue;
                    }
                    if (kount <= 1) {
                        cout << "numsegs: Error, with data !" << endl;
                        exit(0);
                    }
                    if (nod[0] == 1 && nod[1] == 4) {
                        nod[1] = 1;
                    }
                    nu[ind2 + nod[1]] = -nu[ind2 + nod[1]];  //for this also, nu turn-ve
                    int check;
                    check = nu[ind2 + nod[1]];
                    nusg[ind2 + nod[1]] = nseg;
                    nbseg = nbseg - 1;  // This was a boundary segment
                    break;
                }
            }
            nseg++;
            nbseg++;
        }
    }
    //
    //     -------------------------
    //		Reconstitute positive NU
    //     -------------------------
    kount = 0;
    ind = 0;
    for (int ip = 1; ip <= np; ip++) {
        for (int i = 1; i <= 3; i++) {
            nu[ind + i] = -nu[ind + i];
            if (nu[ind + i] < 0) {
                cout << ind + i << endl;
                cout << "tr:" << ip << " Nd:" << nu[ind + 1] << " " << nu[ind + 2] << " " << nu[ind + 3] << endl;
                kount++;
            }
        }
        ind = ind + 3;
    }
    for (int ip = 1; ip <= nq; ip++) {
        for (int i = 1; i <= 4; i++) {
            nu[ind + i] = -nu[ind + i];
            if (nu[ind + i] < 0) {
                cout << "quad: " << ip << " Nd:" << nu[ind + 1] << nu[ind + 2] << nu[ind + 3] << nu[ind + 4] << endl;
                kount++;
            }
        }
        ind = ind + 4;
    }

    if (kount > 0) {
        cout << "numsegs: Error with data !" << endl;
        exit(0);
    }
    //
    //     ------------------------------
    //	   Build NUBO of VIRTUAL segments
    //	   ------------------------------
    ind = np * 3;
    nvseg = 0;
    for (int ip = np + 1; ip <= np + nq; ip++) {
        is1 = nu[ind + 1];
        is2 = nu[ind + 2];
        is3 = nu[ind + 3];
        is4 = nu[ind + 4];
        ivseg1 = nseg + nvseg;
        ivseg2 = nseg + nvseg + 1;
        nubo[0][ivseg1] = is3;
        nubo[1][ivseg1] = is1;
        nubo[0][ivseg2] = is4;
        nubo[1][ivseg2] = is2;
        ind = ind + 4;
        nvseg = nvseg + 2;
    }
}

void DataStructures::fjaret2D() {
    int inod1, inod2, kpoi_1, kpoi_2, max_str;
    //
    //	Jaret: Serial Storage of Jaret

    for_each(ndeg.begin(), ndeg.end(), [](int &i) { i = 0; });
    for_each(jaret.begin(), jaret.end(), [](int &i) { i = 0; });

    //
    //	Find # of nodes - segments around each node
    //	ATT: Jaret also includes neighbours due to virtual segments

    for (int iseg = 0; iseg < nseg + nvseg; iseg++) {
        inod1 = nubo[0][iseg];
        inod2 = nubo[1][iseg];
        ndeg[inod1] = ndeg[inod1] + 1;
        ndeg[inod2] = ndeg[inod2] + 1;
    }
    //
    for (int k = 1; k <= ns; k++) {
        ndeg[k] = ndeg[k - 1] + ndeg[k];  // build index
        jaret[ndeg[k]] = ndeg[k - 1];     // provisory
    }
    //
    // Calculate Max Length of Jaret(String), Only Node Storing
    for (int iseg = 0; iseg < nseg + nvseg; iseg++) {
        inod1 = nubo[0][iseg];
        inod2 = nubo[1][iseg];
        jaret[ndeg[inod1]] = jaret[ndeg[inod1]] + 1;
        kpoi_1 = jaret[ndeg[inod1]];
        jaret[kpoi_1] = inod2;
        jaret[ndeg[inod2]] = jaret[ndeg[inod2]] + 1;
        kpoi_2 = jaret[ndeg[inod2]];
        jaret[kpoi_2] = inod1;
    }
    //
    //	Create Extra Segments
    max_str = ndeg[ns];  // Length of String when all
    // JARET(Nodes) are Stored
    //
    if (2 * ndeg[ns] > nvmaxall) {
        cout << "fjaret:Increase nvmaxall at DataStructures.h-->" << 2 * ndeg[ns] << "++" << endl;
    }
    for (int k = 1; k <= ns; k++) {
        jaret[ndeg[k] + max_str] = max_str + ndeg[k - 1];
    }
    //
    for (int iseg = 0; iseg < nseg + nvseg; iseg++) {
        inod1 = nubo[0][iseg];
        inod2 = nubo[1][iseg];
        jaret[ndeg[inod1] + max_str] = jaret[ndeg[inod1] + max_str] + 1;
        kpoi_1 = jaret[ndeg[inod1] + max_str];
        jaret[kpoi_1] = iseg;
        jaret[ndeg[inod2] + max_str] = jaret[ndeg[inod2] + max_str] + 1;
        kpoi_2 = jaret[ndeg[inod2] + max_str];
        jaret[kpoi_2] = iseg + 1;
    }
}

//===================3D===========================//

DataStructures::DataStructures(vector<int> &logfr_, int ns_, double **coor_, int ntet_, int npyr_,
                               int npri_, int nhex_, int nall_, vector<int> &nu_)
    : logfr(logfr_), ns(ns_), coor(coor_), ntet(ntet_), npyr(npyr_), npri(npri_), nhex(nhex_), nall(nall_), nu(nu_) {
    resizeVectors(maxseg, nvmaxall, maxbfac, maxseg, maxlist, nall);
}

void DataStructures::Create3D() {
    /*-------------------------- -
	Make Lists of Nodes(listn)
	-------------------------- */
    filistn();
    cout << "fillstn completed...." << endl;

    /*------------------------------------------------
	Find the Correspondence of Periodic nodes(iper)
	// ------------------------------------------------*/
    if (isPeriodic) {
        if (isperiph == 1) {
            setperiph();  //  peripheral cascade
            cout << " setperiph completed...." << endl;
        } else {
            setperio();  //  linear cascade
            cout << " setperio completed...." << endl;
        }
        if (nper > 0) {
            cout << " Pint = " << pint;
            cout << " Pext = " << pext;
        }
        numsegs3D();
        virtualPeriodicNeighbours();
        // nper = 0;
        // calvnocl();

        cout << " numsegs completed...." << endl;
    } else {
        numsegs3D();
        virtualPeriodicNeighbours();
        // nper = 0;
        // calvnocl();

        cout << " numsegs completed...." << endl;
    }

    /*----------------------------------------------------
	Calculate + ve Element volumes(vol), may reorder nu
	------------------------------------------------------*/

    // volumes();
    // cout << " volumes completed...." << endl;

    /*-------------------------------------------------- -
	Numerate the faces, provisionally uses jaret(nubf)
	-------------------------------------------------- -*/
    // numfaces();
    // cout << " numfaces completed...." << endl;

    /*
	-------------------------------------------------------
	!     Fill List of boundary faces, painted (listbf), (vnofac)
	!     -------------------------------------------------------
	*/
    // filistbf();
    // cout << " filistbf completed...." << endl;

    // /*
    // -----------------------------------------------
    // !     Numerate the segments (nusg, jaret, ndeg, nubo)
    // !     -----------------------------------------------
    // */
    // if (nper > 0)
    // {
    // 	perseg();
    // 	cout << " perseg completed...." << endl;
    // }

    /*
	-------------------------------------------------
	!     Calculate VNOCL for the segments and cell volume
	!     -------------------------------------------------
	*/
    /*
	--------------------------------------------------
	!     Find pairs of periodic segments (ipersg)
	!     --------------------------------------------------
	*/

    // cout << " calvnocl completed...." << endl;

    /*
	------
	!     Ending
	!     ------
	*/

    cout << " Number of nodes                  : " << ns << endl;
    //	cout << " Number of segments               : " << nseg << endl;
    cout << " Number of tetrahedra              : " << ntet << endl;
    cout << " Number of pyramids                : " << npyr << endl;
    cout << " Number of prisms                  : " << npri << endl;
    cout << " Number of hexahedra               : " << nhex << endl;
    //	cout << " Number of boundary faces          : " << nbfac << endl;

    return;
}

void DataStructures::resizeVectors(int maxseg, int nvmxall, int maxbfac,
                                   int maxsegpe, int maxlist, int nall) {
    nubo = matrix<int>(2, maxseg);
    jaret.resize(nvmxall + 1);
    ibfc2te.resize(maxbfac + 1);
    nubf = matrix<int>(4, maxbfac);
    listnp1.resize(maxlist + 1);
    listnp2.resize(maxlist + 1);
    listbf1.resize(maxlist + 1);
    listbf2.resize(maxlist + 1);
    ibsg2te.resize(maxlist + 1);
    vol.resize(nall + 1);
    cell.resize(ns + 1);
    ndeg.resize(2 * ns + 2);
    nusg.resize(ntet * 6 + npyr * 8 + npri * 9 + nhex * 12 + 1);
}

void DataStructures::filistn() {
    int mposa;
    vector<int> md;
    md.resize(3);

    int kount = 0, kounter;

    for (unsigned int i = 1; i <= ns; i++) {
        if (logfr[i] < 10) {
            kount++;
        }  //appears once
        else if (logfr[i] < 100) {
            kount += 2;
        }  //apears twice
        else {
            kount += 3;
        }  // apears thrice
    }

    listn.resize(kount + 1);

    cout << "Nodes are listed as follows: " << endl;

    for_each(listnp1.begin(), listnp1.end(), [](int i) -> int { return i = 1; });

    for (unsigned int ktype = 0; ktype <= 6; ktype++) {
        if (ktype == 0) {
            kounter = 0;
            listnp1[0] = 1;
        } else {
            listnp1[ktype] = listnp2[ktype - 1] + 1;
            kounter = listnp2[ktype - 1];
        }

        for (unsigned int i = 1; i <= ns; i++) {
            if (logfr[i] == 0) {
                if (ktype == 0) {
                    kounter++;
                    listn[kounter] = i;
                }
            } else {
                analnode(logfr[i], mposa, md[2], md[1], md[0]);

                for (unsigned int mva = 0; mva <= mposa; mva++) {
                    if (md[mva] == ktype) {
                        kounter++;
                        listn[kounter] = i;
                    }
                }
            }
        }
        listnp2[ktype] = kounter;
        if (listnp2[ktype] < listnp1[ktype]) {
            listnp2[ktype] = listnp2[ktype - 1];
        } else {
            cout << " Nodes of type " << ktype << " :between :" << listnp1[ktype] << " " << listnp2[ktype] << endl;
        }

        continue;
    }
    cout << " ------  ------  ------  ------  ------ " << endl;
    return;
}

void DataStructures::analnode(int kkk0, int &mposa, int &m100,
                              int &m010, int &m001) {
    mposa = 2;

    m100 = kkk0 / 100;
    if (m100 == 0) {
        mposa = 1;
    }
    kkk0 -= m100 * 100;
    m010 = kkk0 / 10;
    if (m010 == 0) {
        mposa = 0;
    }
    m001 = kkk0 - m010 * 10;
}

void DataStructures::setPeriodicityConstants(bool isPeriodic_, int mpar_, int kaxial_, int isperiph_) {
    isPeriodic = isPeriodic_;
    mpar = mpar_;
    kaxial = kaxial_;
    isperiph = isperiph_;
}

void DataStructures::setperio() {
    double dyy, dzz, ds;
    int is, js;
    nper = 0;    //default: no periodic pairs
    pitch = 0.;  // default: zero pitch

    if (listnp2[1] <= listnp1[1] + 1) {
        return;
    }  // exit if no periodic nodes

    int nper_exp = (listnp2[1] - listnp1[1] + 1) / 2;

    iper = matrix<int>(2, nper_exp);

    /*if (kaxial != 3)
	{
	cout << "KAXIAL not eq to 3 in peripheral case" << endl;
	cout << "It will be adjusted accordingly" << endl;
	kaxial = 3;
	}*/

    // "Angles" used in peripheral cascade are taken on default values

    pint = 1.;
    pext = 0.;

    // Periodic nodes to be matched based on D(COORDS)<EPSILON

    double eps = 1.e-4;
    cout << "Periodic eps: " << eps;

    // Matching Periodic Nodes, one by one

    for (unsigned int idum = listnp1[1]; idum <= listnp2[1]; idum++)  // Sweep all nodes with LOGFR=1
    {
        is = listn[idum];
        if (is < 0) {
            listn[idum] = -is;
            continue;
        }  // pair already treated

        for (unsigned int idum1 = idum + 1; idum1 <= listnp2[1]; idum1++) {
            js = listn[idum1];
            if (js < 0) {
                continue;
            }
            dyy = abs(coor[0][is - 1] - coor[0][js - 1]);
            dzz = abs(coor[2][is - 1] - coor[2][js - 1]);
            if (dyy < eps && dzz < eps) {
                listn[idum1] = -listn[idum1];
                nper++;
                iper[0][nper - 1] = is;
                iper[1][nper - 1] = js;
                goto Line10;
            }
        }
        cout << "Warning: No pair found for node " << is << ", " << coor[0][is - 1] << endl;
    Line10:;
    }
    // If sequential run (not a slave) all LOGFR=1 nodes should be matched

    cout << " Pairs of Periodic nodes = " << nper << endl;
    if (2 * nper != listnp2[1] - listnp1[1] + 1) {
        cout << "ERROR in SETPERIO" << endl;
        exit;
    }

    // Reorder IPER so as IPER(1) be the node with min_X

    double dsmax = -9.e28;
    double dsmin = 9.e28;

    for (unsigned int kper = 0; kper < nper; kper++) {
        is = iper[0][kper];
        js = iper[1][kper];
        if (coor[1][is - 1] > coor[1][js - 1]) {
            iper[0][kper] = js;
            iper[1][kper] = is;
            is = iper[0][kper];  // renew is,js for ds
            js = iper[1][kper];
        }
        ds = coor[1][js - 1] - coor[1][is - 1];
        if (ds > dsmax)
            dsmax = ds;
        if (ds < dsmin)
            dsmin = ds;
    }
    if (dsmax * dsmin < 0. && nper > 0) {
        cout << " Wrong Orientation in SETPERIO" << endl;
        exit;
    }
    /* Prints and saves the pitch, only if Sequential Run
	otherwise it will be sent by master
	--------------------------------------------------*/
    mpar = 10;
    if (mpar == 0) {
        cout << "Found PITCH between: " << dsmin << " " << dsmax << endl;
    } else {
        pitch = 0.5 * (dsmin + dsmax);
    }
    return;
}

void DataStructures::setperiph() {
    int is, js;
    double x1, y1, z1, r1, x2, y2, z2, r2, drr, dzz, r12,
        pint, pint_1, pext_1, pext, phi, phideg, check;
    nper = 0;
    pitch = 0.;

    if (listnp2[1] <= listnp1[1]) {
        return;
    }  // exit if no periodic nodes

    int nper_exp = (listnp2[1] - listnp1[1] + 1) / 2;

    iper = matrix<int>(2, nper_exp);

    if (kaxial != 3) {
        cout << "KAXIAL not eq to 3 in peripheral case" << endl;
        cout << "It will be adjusted accordingly" << endl;
        kaxial = 3;
    }

    double phimax = -9.e28;
    double phimin = 9.e28;

    // Periodic nodes to be matched based on D(COORDS) + RADIUS < EPSILON

    double eps = 1.e-04;
    cout << " Periodic eps: " << eps << endl;

    // Matching Periodic Nodes, one by one

    for (int idum = listnp1[1]; idum <= listnp2[1]; idum++)  //10
    {
        is = listn[idum];
        if (is < 0) {
            listn[idum] = -is;  // pair already treated
            continue;
        }
        x1 = coor[0][is - 1];
        y1 = coor[1][is - 1];
        z1 = coor[2][is - 1];
        r1 = sqrt(x1 * x1 + y1 * y1);
        for (int idum1 = idum + 1; idum1 <= listnp2[1]; idum1++)  //20
        {
            js = listn[idum1];
            if (js < 0)
                continue;
            x2 = coor[0][js - 1];
            y2 = coor[1][js - 1];
            z2 = coor[2][js - 1];
            r2 = sqrt(x2 * x2 + y2 * y2);
            drr = abs(r1 - r2);  // Same R?
            dzz = abs(z1 - z2);  // Same Z?
            if (drr < eps && dzz < eps) {
                listn[idum1] = -listn[idum1];
                nper++;
                iper[0][nper - 1] = is;
                iper[1][nper - 1] = js;
                r12 = r1 * r2;
                if (nper == 1 && mpar == 0)  // first periodic, serial
                {
                    pint = (x1 * x2 + y1 * y2) / r12;  // cos(pitch) LAW
                    pext = (x1 * y2 - x2 * y1) / r12;  // sin(pitch) LAW
                    phi = abs(acos(pint));
                    if (phi > phimax)
                        phimax = phi;
                    if (phi < phimin)
                        phimin = phi;
                    goto LINE10;
                } else {
                    pint_1 = (x1 * x2 + y1 * y2) / r12;  // cos(pitch)
                    pext_1 = (x1 * y2 - x2 * y1) / r12;  // sin(pitch)
                    check = pext * pext_1;
                    if (check < 0) {
                        iper[0][nper - 1] = js;
                        iper[1][nper - 1] = is;
                    }
                    phi = abs(acos(pint_1));
                    if (phi > phimax)
                        phimax = phi;
                    if (phi < phimin)
                        phimin = phi;
                    goto LINE10;
                }
            }
        }
        cout << "Warning: No pair found for node " << is << ", z " << coor[2][is - 1] << endl;
    LINE10:;
    }
    // If sequential run(not a slave) all LOGFR = 1 nodes should be matched
    //------------------------------------------------------------------ -
    if (mpar == 0) {
        cout << "Found PITCH between (rad) :" << phimin << " " << phimax << endl;
        if (abs(phimax - phimin) > eps) {
            cout << "ERROR:Large variation of pitch at peripheral" << endl;
            exit;
        } else {
            pitch = (phimin + phimax) * 0.5;  // pitch angle in radians
            phideg = pitch * 180. / (4. * atan(1.));
            cout << " PERIPHERAL CASCADE angle(deg) = " << phideg << endl;
        }
    }

    return;
}

void DataStructures::volumes() {
    qualityCheck gcheck;
    int pyra[4][4] = {{1, 2, 3, 4}, {2, 3, 4, 1}, {3, 4, 1, 2}, {5, 5, 5, 5}};
    int prism[4][10] = {{1, 3, 6, 5, 2, 3, 1, 5, 4, 2}, {2, 1, 5, 4, 1, 2, 2, 4, 1, 5}, {3, 2, 4, 6, 4, 5, 3, 6, 3, 6}, {6, 5, 1, 3, 6, 4, 4, 2, 5, 1}};
    int hexa[4][32] = {{1, 2, 3, 4, 5, 8, 7, 6, 2, 4, 1, 3, 1, 2, 3, 4, 5, 8, 7, 6, 5, 3, 4, 2, 1, 2, 3, 4, 5, 6, 7, 8},
                       {2, 3, 4, 1, 8, 7, 6, 5, 1, 3, 4, 2, 2, 3, 4, 1, 8, 7, 6, 5, 1, 7, 8, 6, 2, 3, 4, 1, 8, 5, 6, 7},
                       {3, 4, 1, 2, 7, 6, 5, 8, 5, 7, 8, 6, 3, 4, 1, 2, 7, 6, 5, 8, 4, 8, 5, 7, 4, 1, 2, 3, 6, 7, 8, 5},
                       {7, 8, 5, 6, 3, 2, 1, 4, 8, 6, 7, 5, 5, 6, 7, 8, 1, 4, 3, 2, 6, 2, 3, 1, 5, 6, 7, 8, 1, 2, 3, 4}};

    double dx12, dx13, dx14, dy12, dy13, dy14, dz12, dz13, dz14, dxa, dya, dza;
    double dx15, dy15, dz15, dx24, dy24, dz24, ajac;
    int is1, is2, is3, is4, is5, is6, is7, is8, kountj, kountv;
    int js1, js2, js3, js4;
    double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, x12, x13, x14,
        y12, y13, y14, z12, z13, z14, xre, yre, zre, ajaco;
    double dx521, dy521, dz521, dx31, dy31, dz31, dx54, dy54, dz54, vo1,
        dx51, dy51, dz51, dx52, dy52, dz52, vo2;
    double dx7281, dy7281, dz7281, dx61, dy61, dz61, dx74, dy74, dz74,
        dx81, dy81, dz81, dx72, dy72, dz72, dx75, dy75, dz75, dx7461, dy7461, dz7461,
        dx7531, dy7531, dz7531, vo3;
    double dx45, dy45, dz45, dx46, dy46, dz46, dx62, dy62, dz62,
        dx6351, dy6351, dz6351, dx64, dy64, dz64, dx53, dy53, dz53;

    double us6 = 1. / 6.;
    double us12 = 1. / 12.;
    int ind = 0;

    // this is for tetrahedra
    double volmin = 1e28;
    int kount = 0;

    for (unsigned int ip = 1; ip <= ntet; ip++) {
        // chexk orientation ana calculate volume
        is1 = nu[ind + 1];
        is2 = nu[ind + 2];
        is3 = nu[ind + 3];
        is4 = nu[ind + 4];
        dx12 = coor[0][is2 - 1] - coor[0][is1 - 1];
        dx13 = coor[0][is3 - 1] - coor[0][is1 - 1];
        dx14 = coor[0][is4 - 1] - coor[0][is1 - 1];
        dy12 = coor[1][is2 - 1] - coor[1][is1 - 1];
        dy13 = coor[1][is3 - 1] - coor[1][is1 - 1];
        dy14 = coor[1][is4 - 1] - coor[1][is1 - 1];
        dz12 = coor[2][is2 - 1] - coor[2][is1 - 1];
        dz13 = coor[2][is3 - 1] - coor[2][is1 - 1];
        dz14 = coor[2][is4 - 1] - coor[2][is1 - 1];

        gcheck.vecProd(dx12, dy12, dz12, dx13, dy13, dz13, dxa, dya, dza);

        vol[ip] = (dx14 * dxa + dy14 * dya + dz14 * dza) * us6;
        if (vol[ip] <= 0.)
        //  uncomment if I want to allow reordering
        {
            vol[ip] = -vol[ip];
            nu[ind + 2] = is3;  // reorder
            nu[ind + 3] = is2;  // reorder
            kount++;

            //  uncomment to stop if invalid element is found

            /*cout << " Invalid tetrahedral found in volumes" << endl;
			cout << " Tetrahedral no. " << ip << endl;
			cout << "Nodes no. " << is1 << " " << is2 << " " << is3 << " " << is4 << endl;
			exit;*/
        }
        if (vol[ip] < volmin)
            volmin = vol[ip];
        ind += 4;
    }
    if (ntet > 0) {
        cout << " min vol at tetrahedra    :  " << volmin << endl;
        cout << " reorderings at tetrahedra:  " << kount << endl;
    }
    // this is for pyramids
    volmin = 1.e28;
    kount = 0;                                                   //reorderings
    kountj = 0;                                                  //negative jacobians
    for (unsigned int ip = ntet + 1; ip <= ntet + npyr; ip++) {  // check orientation
        is1 = nu[ind + 1];
        is2 = nu[ind + 2];
        is3 = nu[ind + 3];
        is4 = nu[ind + 4];
        is5 = nu[ind + 5];
        dx13 = coor[0][is3 - 1] - coor[0][is1 - 1];
        dx24 = coor[0][is4 - 1] - coor[0][is2 - 1];
        dx15 = coor[0][is5 - 1] - coor[0][is1 - 1];
        dy13 = coor[1][is3 - 1] - coor[1][is1 - 1];
        dy24 = coor[1][is4 - 1] - coor[1][is2 - 1];
        dy15 = coor[1][is5 - 1] - coor[1][is1 - 1];
        dz13 = coor[2][is3 - 1] - coor[2][is1 - 1];
        dz24 = coor[2][is4 - 1] - coor[2][is2 - 1];
        dz15 = coor[2][is5 - 1] - coor[2][is1 - 1];

        gcheck.vecProd(dx13, dy13, dz13, dx24, dy24, dz24, dxa, dya, dza);
        ajac = (dx15 * dxa + dy15 * dya + dz15 * dza);

        if (ajac <= 0.)
        //  uncomment if I want to allow reordering
        {
            /*
			nu[ind + 4] = is2; // reorder
			nu[ind + 2] = is4; // reorder
			is2 = nu[ind + 2];
			is4 = nu[ind + 4];
			kount++;
			*/

            //  uncomment to stop if invalid element is found

            cout << " Invalid pyramid found in volumes" << endl;
            cout << " Pyramid no. " << ip << endl;
            cout << "Nodes no. " << is1 << ", " << is2 << ", " << is3 << ", "
                 << is4 << ", " << is5 << endl;
            exit;
        }

        // check jacobian
        for (unsigned int i = 1; i <= 4; i++) {
            js1 = nu[ind + pyra[0][i - 1]];
            js2 = nu[ind + pyra[1][i - 1]];
            js3 = nu[ind + pyra[2][i - 1]];
            js4 = nu[ind + pyra[3][i - 1]];

            x1 = coor[0][js1 - 1];
            y1 = coor[1][js1 - 1];
            z1 = coor[2][js1 - 1];
            x2 = coor[0][js2 - 1];
            y2 = coor[1][js2 - 1];
            z2 = coor[2][js2 - 1];
            x3 = coor[0][js3 - 1];
            y3 = coor[1][js3 - 1];
            z3 = coor[2][js3 - 1];
            x4 = coor[0][js4 - 1];
            y4 = coor[1][js4 - 1];
            z4 = coor[2][js4 - 1];

            x12 = x2 - x1;
            y12 = y2 - y1;
            z12 = z2 - z1;
            x13 = x3 - x1;
            y13 = y3 - y1;
            z13 = z3 - z1;
            x14 = x4 - x1;
            y14 = y4 - y1;
            z14 = z4 - z1;

            gcheck.vecProd(x12, y12, z12, x13, y13, z13, xre, yre, zre);
            ajaco = xre * x14 + yre * y14 + zre * z14;
            if (ajaco < 0.) {
                cout << "Negative jacobian at pyramid element: " << ip << endl;
                cout << "Nodes: " << js1 << ", " << js2 << ", " << js3 << ", " << js4 << endl;
                cout << "Jacobian: " << ajaco << endl;
                kountj++;
            }
        }

        // calculate volume
        dx521 = 2. * coor[0][is5 - 1] - coor[0][is2 - 1] - coor[0][is1 - 1];
        dx54 = coor[0][is5 - 1] - coor[0][is4 - 1];
        dx31 = coor[0][is3 - 1] - coor[0][is1 - 1];
        dy521 = 2. * coor[1][is5 - 1] - coor[1][is2 - 1] - coor[1][is1 - 1];
        dy54 = coor[1][is5 - 1] - coor[1][is4 - 1];
        dy31 = coor[1][is3 - 1] - coor[1][is1 - 1];
        dz521 = 2. * coor[2][is5 - 1] - coor[2][is2 - 1] - coor[2][is1 - 1];
        dz54 = coor[2][is5 - 1] - coor[2][is4 - 1];
        dz31 = coor[2][is3 - 1] - coor[2][is1 - 1];

        gcheck.vecProd(dx521, dy521, dz521, dx54, dy54, dz54, xre, yre, zre);
        vo1 = dx31 * xre + dy31 * yre + dz31 * zre;

        dx52 = coor[0][is5 - 1] - coor[0][is2 - 1];
        dx51 = coor[0][is5 - 1] - coor[0][is1 - 1];
        dx31 = coor[0][is3 - 1] - coor[0][is1 - 1];
        dy52 = coor[1][is5 - 1] - coor[1][is2 - 1];
        dy51 = coor[1][is5 - 1] - coor[1][is1 - 1];
        dy31 = coor[1][is3 - 1] - coor[1][is1 - 1];
        dz52 = coor[2][is5 - 1] - coor[2][is2 - 1];
        dz51 = coor[2][is5 - 1] - coor[2][is1 - 1];
        dz31 = coor[2][is3 - 1] - coor[2][is1 - 1];

        gcheck.vecProd(dx52, dy52, dz52, dx51, dy51, dz51, xre, yre, zre);
        vo2 = dx31 * xre + dy31 * yre + dz31 * zre;

        vol[ip] = (vo1 + vo2) * us12;
        if (vol[ip] <= 0.) {
            cout << "Negative vol Pyramid" << endl;
            cout << "El = " << ip << " Vol = " << vol[ip] << endl;
            exit;
        }
        if (vol[ip] < volmin)
            volmin = vol[ip];
        ind += 5;
    }
    if (npyr > 0) {
        cout << " min vol at pyramids    :" << volmin << endl;
        cout << " reorderings at pyramids:" << kount << endl;
    }
    // this is for prisms
    volmin = 1.e28;
    kount = 0;   // reorderings
    kountj = 0;  // negative jacobians
    kountv = 0;  // negative volumes
    for (unsigned int ip = ntet + npyr + 1; ip <= ntet + npyr + npri; ip++)
    // check orientation
    {
        is1 = nu[ind + 1];
        is2 = nu[ind + 2];
        is3 = nu[ind + 3];
        is4 = nu[ind + 4];
        is5 = nu[ind + 5];
        is6 = nu[ind + 6];
        dx12 = coor[0][is2 - 1] - coor[0][is1 - 1];
        dx13 = coor[0][is3 - 1] - coor[0][is1 - 1];
        dx14 = coor[0][is4 - 1] - coor[0][is1 - 1];
        dx45 = coor[0][is5 - 1] - coor[0][is4 - 1];
        dx46 = coor[0][is6 - 1] - coor[0][is4 - 1];
        dy12 = coor[1][is2 - 1] - coor[1][is1 - 1];
        dy13 = coor[1][is3 - 1] - coor[1][is1 - 1];
        dy14 = coor[1][is4 - 1] - coor[1][is1 - 1];
        dy45 = coor[1][is5 - 1] - coor[1][is4 - 1];
        dy46 = coor[1][is6 - 1] - coor[1][is4 - 1];
        dz12 = coor[2][is2 - 1] - coor[2][is1 - 1];
        dz13 = coor[2][is3 - 1] - coor[2][is1 - 1];
        dz14 = coor[2][is4 - 1] - coor[2][is1 - 1];
        dz45 = coor[2][is5 - 1] - coor[2][is4 - 1];
        dz46 = coor[2][is6 - 1] - coor[2][is4 - 1];

        gcheck.vecProd(dx12, dy12, dz12, dx13, dy13, dz13, dxa, dya, dza);
        vo1 = (dx14 * dxa + dy14 * dya + dz14 * dza);
        gcheck.vecProd(dx45, dy45, dz45, dx46, dy46, dz46, dxa, dya, dza);
        vo2 = (dx14 * dxa + dy14 * dya + dz14 * dza);

        if (vo1 < 0.)
        //uncomment if I want to allow reordering
        {
            nu[ind + 2] = is3;  //1nd reorder
            nu[ind + 3] = is2;  //1nd reorder
            is2 = nu[ind + 2];
            is3 = nu[ind + 3];
            kount++;

            //uncomment to stop if invalid element is found
            /*
			cout << "Invalid prism found in volumes" << endl;
			cout << " Prism no. "<< ip << endl;
			cout << "Nodes no. " << is1 << ", " << is2 << ", " << is3 << ", "
			<< is4 << ", " << is5 <<", " << is6 << endl;
			exit;
			*/
        }
        if (vo2 < 0.)
        //uncomment if I want to allow reordering
        {
            nu[ind + 5] = is6;  //2st reorder
            nu[ind + 6] = is5;  //2st reorder
            is5 = nu[ind + 5];
            is6 = nu[ind + 6];
            kount++;

            //uncomment to stop if invalid element is found
            /*
			cout << "Invalid prism found in volumes" << endl;
			cout << " Prism no. "<< ip << endl;
			cout << "Nodes no. " << is1 << ", " << is2 << ", " << is3 << ", "
			<< is4 << ", " << is5 <<", " << is6 << endl;
			exit;
			*/
        }

        // check jacobian

        for (unsigned int i = 1; i <= 10; i++) {
            js1 = nu[ind + prism[0][i - 1]];
            js2 = nu[ind + prism[1][i - 1]];
            js3 = nu[ind + prism[2][i - 1]];
            js4 = nu[ind + prism[3][i - 1]];

            x1 = coor[0][js1 - 1];
            y1 = coor[1][js1 - 1];
            z1 = coor[2][js1 - 1];
            x2 = coor[0][js2 - 1];
            y2 = coor[1][js2 - 1];
            z2 = coor[2][js2 - 1];
            x3 = coor[0][js3 - 1];
            y3 = coor[1][js3 - 1];
            z3 = coor[2][js3 - 1];
            x4 = coor[0][js4 - 1];
            y4 = coor[1][js4 - 1];
            z4 = coor[2][js4 - 1];

            x12 = x2 - x1;
            y12 = y2 - y1;
            z12 = z2 - z1;
            x13 = x3 - x1;
            y13 = y3 - y1;
            z13 = z3 - z1;
            x14 = x4 - x1;
            y14 = y4 - y1;
            z14 = z4 - z1;

            gcheck.vecProd(x12, y12, z12, x13, y13, z13, xre, yre, zre);
            ajaco = xre * x14 + yre * y14 + zre * z14;
            if (ajaco < 0.) {
                cout << "Negative jacobian at pyramid element: " << ip << endl;
                cout << "Nodes: " << js1 << ", " << js2 << ", " << js3 << ", " << js4 << endl;
                cout << "Jacobian: " << ajaco << endl;
                kountj++;
            }
        }

        //calculate volume
        dx31 = coor[0][is3 - 1] - coor[0][is1 - 1];
        dx62 = coor[0][is6 - 1] - coor[0][is2 - 1];
        dx6351 = coor[0][is6 - 1] - coor[0][is3 - 1] + coor[0][is5 - 1] - coor[0][is1 - 1];
        dy31 = coor[1][is3 - 1] - coor[1][is1 - 1];
        dy62 = coor[1][is6 - 1] - coor[1][is2 - 1];
        dy6351 = coor[1][is6 - 1] - coor[1][is3 - 1] + coor[1][is5 - 1] - coor[1][is1 - 1];
        dz31 = coor[2][is3 - 1] - coor[2][is1 - 1];
        dz62 = coor[2][is6 - 1] - coor[2][is2 - 1];
        dz6351 = coor[2][is6 - 1] - coor[2][is3 - 1] + coor[2][is5 - 1] - coor[2][is1 - 1];

        gcheck.vecProd(dx31, dy31, dz31, dx62, dy62, dz62, xre, yre, zre);
        vo1 = dx6351 * xre + dy6351 * yre + dz6351 * zre;

        dx64 = coor[0][is6 - 1] - coor[0][is4 - 1];
        dx61 = coor[0][is6 - 1] - coor[0][is1 - 1];
        dx53 = coor[0][is5 - 1] - coor[0][is3 - 1];
        dy64 = coor[1][is6 - 1] - coor[1][is4 - 1];
        dy61 = coor[1][is6 - 1] - coor[1][is1 - 1];
        dy53 = coor[1][is5 - 1] - coor[1][is3 - 1];
        dz64 = coor[2][is6 - 1] - coor[2][is4 - 1];
        dz61 = coor[2][is6 - 1] - coor[2][is1 - 1];
        dz53 = coor[2][is5 - 1] - coor[2][is3 - 1];

        gcheck.vecProd(dx64, dy64, dz64, dx61, dy61, dz61, xre, yre, zre);
        vo2 = dx53 * xre + dy53 * yre + dz53 * zre;

        dx64 = coor[0][is6 - 1] - coor[0][is4 - 1];
        dx62 = coor[0][is6 - 1] - coor[0][is2 - 1];
        dx51 = coor[0][is5 - 1] - coor[0][is1 - 1];
        dy64 = coor[1][is6 - 1] - coor[1][is4 - 1];
        dy62 = coor[1][is6 - 1] - coor[1][is2 - 1];
        dy51 = coor[1][is5 - 1] - coor[1][is1 - 1];
        dz64 = coor[2][is6 - 1] - coor[2][is4 - 1];
        dz62 = coor[2][is6 - 1] - coor[2][is2 - 1];
        dz51 = coor[2][is5 - 1] - coor[2][is1 - 1];

        gcheck.vecProd(dx64, dy64, dz64, dx62, dy62, dz62, xre, yre, zre);
        vo3 = dx51 * xre + dy51 * yre + dz51 * zre;

        vol[ip] = (vo1 + vo2 + vo3) * us12;
        if (vol[ip] < volmin)
            volmin = vol[ip];
        if (vol[ip] <= 0.) {
            cout << "Negative vol Prisms" << endl;
            cout << "El = " << ip << " Vol = " << vol[ip] << endl;
            exit;
        }
        ind += 6;
    }
    if (npri > 0) {
        cout << " min vol at prisms    :" << volmin << endl;
        cout << " reordering prisms    :" << kount << endl;
        cout << " negative jacobians at prisms: " << kountj << endl;
    }

    // this is for hexahedra

    volmin = 1.e28;
    kount = 0;   // reorderings
    kountj = 0;  // negative jacobians
    kountv = 0;  // negative volumes

    for (unsigned int ip = ntet + npyr + npri + 1; ip <= ntet + npyr + npri + nhex; ip++) {
        is1 = nu[ind + 1];
        is2 = nu[ind + 2];
        is3 = nu[ind + 3];
        is4 = nu[ind + 4];
        is5 = nu[ind + 5];
        is6 = nu[ind + 6];
        is7 = nu[ind + 7];
        is8 = nu[ind + 8];

        // check jacobian
        for (unsigned int i = 1; i <= 32; i++) {
            js1 = nu[ind + hexa[0][i - 1]];
            js2 = nu[ind + hexa[1][i - 1]];
            js3 = nu[ind + hexa[2][i - 1]];
            js4 = nu[ind + hexa[3][i - 1]];

            x1 = coor[0][js1 - 1];
            y1 = coor[1][js1 - 1];
            z1 = coor[2][js1 - 1];
            x2 = coor[0][js2 - 1];
            y2 = coor[1][js2 - 1];
            z2 = coor[2][js2 - 1];
            x3 = coor[0][js3 - 1];
            y3 = coor[1][js3 - 1];
            z3 = coor[2][js3 - 1];
            x4 = coor[0][js4 - 1];
            y4 = coor[1][js4 - 1];
            z4 = coor[2][js4 - 1];

            x12 = x2 - x1;
            y12 = y2 - y1;
            z12 = z2 - z1;
            x13 = x3 - x1;
            y13 = y3 - y1;
            z13 = z3 - z1;
            x14 = x4 - x1;
            y14 = y4 - y1;
            z14 = z4 - z1;

            gcheck.vecProd(x12, y12, z12, x13, y13, z13, xre, yre, zre);
            ajaco = xre * x14 + yre * y14 + zre * z14;
            if (ajaco < 0.)
                if (ajaco < 0.) {
                    cout << "Negative jacobian at pyramid element: " << ip << endl;
                    cout << "Nodes: " << js1 << ", " << js2 << ", " << js3 << ", " << js4 << endl;
                    cout << "Jacobian: " << ajaco << endl;
                    kountj++;
                }
        }
        // calculate volume
        dx7281 = coor[0][is7 - 1] - coor[0][is2 - 1] + coor[0][is8 - 1] - coor[0][is1 - 1];
        dx74 = coor[0][is7 - 1] - coor[0][is4 - 1];
        dx31 = coor[0][is3 - 1] - coor[0][is1 - 1];
        dy7281 = coor[1][is7 - 1] - coor[1][is2 - 1] + coor[1][is8 - 1] - coor[1][is1 - 1];
        dy74 = coor[1][is7 - 1] - coor[1][is4 - 1];
        dy31 = coor[1][is3 - 1] - coor[1][is1 - 1];
        dz7281 = coor[2][is7 - 1] - coor[2][is2 - 1] + coor[2][is8 - 1] - coor[2][is1 - 1];
        dz74 = coor[2][is7 - 1] - coor[2][is4 - 1];
        dz31 = coor[2][is3 - 1] - coor[2][is1 - 1];

        gcheck.vecProd(dx7281, dy7281, dz7281, dx74, dy74, dz74, xre, yre, zre);
        vo1 = dx31 * xre + dy31 * yre + dz31 * zre;

        dx81 = coor[0][is8 - 1] - coor[0][is1 - 1];
        dx7461 = coor[0][is7 - 1] - coor[0][is4 - 1] + coor[0][is1 - 1];
        dx75 = coor[0][is7 - 1] - coor[0][is5 - 1];
        dy81 = coor[1][is8 - 1] - coor[1][is1 - 1];
        dy7461 = coor[1][is7 - 1] - coor[1][is4 - 1] + coor[1][is1 - 1];
        dy75 = coor[1][is7 - 1] - coor[1][is5 - 1];
        dz81 = coor[2][is8 - 1] - coor[2][is1 - 1];
        dz7461 = coor[2][is7 - 1] - coor[2][is4 - 1] + coor[2][is1 - 1];
        dz75 = coor[2][is7 - 1] - coor[2][is5 - 1];

        gcheck.vecProd(dx81, dy81, dz81, dx7461, dy7461, dz7461, xre, yre, zre);
        vo2 = dx75 * xre + dy75 * yre + dz75 * zre;

        dx72 = coor[0][is7 - 1] - coor[0][is2 - 1];
        dx61 = coor[0][is6 - 1] - coor[0][is1 - 1];
        dx7531 = coor[0][is7 - 1] - coor[0][is5 - 1] + coor[0][is3 - 1] - coor[0][is1 - 1];
        dy72 = coor[1][is7 - 1] - coor[1][is2 - 1];
        dy61 = coor[1][is6 - 1] - coor[1][is1 - 1];
        dy7531 = coor[1][is7 - 1] - coor[1][is5 - 1] + coor[1][is3 - 1] - coor[1][is1 - 1];
        dz72 = coor[2][is7 - 1] - coor[2][is2 - 1];
        dz61 = coor[2][is6 - 1] - coor[2][is1 - 1];
        dz7531 = coor[2][is7 - 1] - coor[2][is5 - 1] + coor[2][is3 - 1] - coor[2][is1 - 1];

        gcheck.vecProd(dx72, dy72, dz72, dx61, dy61, dz61, xre, yre, zre);
        vo3 = dx7531 * xre + dy7531 * yre + dz7531 * zre;

        vol[ip] = (vo1 + vo2 + vo3) * us12;
        if (vol[ip] < volmin)
            volmin = vol[ip];
        if (vol[ip] <= 0.) {
            cout << "Negative vol Hexahedral" << endl;
            cout << "El = " << ip << " Vol = " << vol[ip] << endl;
            exit;
        }
        ind += 8;
    }
    if (nhex > 0) {
        cout << " min vol at hexahedra    :" << volmin << endl;
        cout << " ATT: No checks for reorderings in hexahedrals" << endl;
        cout << " negative jacobians at hexahedra: " << kountj << endl;
    }
    return;
}

void DataStructures::numfaces() {
    //  nodes of the faces, tetrahedra
    int iaux1[4][3] = {{2, 3, 4}, {1, 4, 3}, {1, 2, 4}, {1, 3, 2}};
    //  nodes of the faces, pyramids
    int iaux2[5][4] = {{1, 2, 5, 0}, {3, 4, 5, 0}, {4, 1, 5, 0}, {2, 3, 5, 0}, {1, 4, 3, 2}};
    //  nodes of the faces, prisms
    int iaux3[5][4] = {{1, 2, 5, 4}, {1, 3, 2, 0}, {4, 5, 6, 0}, {2, 3, 6, 5}, {1, 4, 6, 3}};
    //  nodes of the faces, hexahedra
    int iaux4[6][4] = {{2, 3, 7, 6}, {1, 5, 8, 4}, {1, 2, 6, 5}, {5, 6, 7, 8}, {1, 4, 3, 2}, {3, 4, 8, 7}};

    int is1, is2, is3, is4;
    int ielem, ind2, kk, ndif, ind2f;
    int k1, k2, kount, indp;
    vector<int> nod;
    nod.resize(5);

    /*
	Only in this routine NU turns - ve and
	JARET stores, provisionally, the pyramids around the node
	*/
    int mnodes, ind;
    nbfac = 0;
    fjaret_el();

    /*

	=======================
	FIND the boundary faces
	=======================

	*/
    ind = 0;  // node index

    /*
	!     Loop on elements
	!     -----------------------------
	*/

    for (unsigned int ip = 1; ip <= ntet + npyr + npri + nhex; ip++)
    //---------------------------- -
    {
        if (ip <= ntet)
            mnodes = 4;  // find the kind of the elem.
        else if (ip <= ntet + npyr)
            mnodes = 5;
        else if (ip <= ntet + npyr + npri)
            mnodes = 6;
        else
            mnodes = 8;

        for (unsigned int m = 1; m <= mnodes; m++)  // Loop on the nodes of each element
        {
            if (nu[ind + m] < 0)
                continue;
            if (m == 6 && mnodes == 6) {
                nu[ind + m] = -nu[ind + m];  // turn negative, no face attached to it
                continue;
            }
            if (m > 6 && mnodes == 8) {
                nu[ind + m] = -nu[ind + m];  // turn negative, no face attached to it
                continue;
            }

            //Else, this is a NEW face, then:
            if (ip <= ntet) {
                is1 = abs(nu[ind + iaux1[m - 1][0]]);
                is2 = abs(nu[ind + iaux1[m - 1][1]]);
                is3 = abs(nu[ind + iaux1[m - 1][2]]);
                is4 = 0;
            } else if (ip <= ntet + npyr) {
                is1 = abs(nu[ind + iaux2[m - 1][0]]);
                is2 = abs(nu[ind + iaux2[m - 1][1]]);
                is3 = abs(nu[ind + iaux2[m - 1][2]]);
                if (iaux2[m - 1][3] == 0)
                    is4 = 0;
                else
                    is4 = abs(nu[ind + iaux2[m - 1][3]]);
            } else if (ip <= ntet + npyr + npri) {
                is1 = abs(nu[ind + iaux3[m - 1][0]]);
                is2 = abs(nu[ind + iaux3[m - 1][1]]);
                is3 = abs(nu[ind + iaux3[m - 1][2]]);
                if (iaux3[m - 1][3] == 0)
                    is4 = 0;
                else
                    is4 = abs(nu[ind + iaux3[m - 1][3]]);
            } else {
                is1 = abs(nu[ind + iaux4[m - 1][0]]);
                is2 = abs(nu[ind + iaux4[m - 1][1]]);
                is3 = abs(nu[ind + iaux4[m - 1][2]]);
                is4 = abs(nu[ind + iaux4[m - 1][3]]);
            }

            nbfac++;

            if (nbfac > maxbfac) {
                cout << " numfaces tet: increase MAXBFAC" << endl;
                exit;
            }

            nubf[0][nbfac - 1] = is1;
            nubf[1][nbfac - 1] = is2;
            nubf[2][nbfac - 1] = is3;
            nubf[3][nbfac - 1] = is4;

            nu[ind + m] = -nu[ind + m];  // attention, nu turns -ve
            ibfc2te[nbfac] = ip;

            /*
			Find adjacent face to avoid duplicating.
			Also mark the faces that lie over the boundary.
			*/
            int ds;
            for (unsigned int ik = ndeg[is1 - 1] + 1; ik <= ndeg[is1]; ik++) {
            LOOP42:
                if (ik > ndeg[is1])
                    break;
                ds = ndeg[is1];
                ielem = jaret[ik];
                if (ielem == ip)
                    continue;       // skip the element at hand
                if (ielem <= ntet)  // find the kind of the adjacent elem.
                {
                    if (mnodes == 8)
                        continue;  // hex can't be adjacent to tet
                    if (is4 > 0)
                        continue;  // quad face can't belong to tet
                    kk = 4;
                    ind2 = (ielem - 1) * 4;
                    ind2f = ind2;
                    ndif = 2;
                } else if (ielem <= ntet + npyr) {
                    kk = 5;
                    ind2 = ntet * 4 + (ielem - ntet - 1) * 5;
                    ind2f = ind2;
                    ndif = 3;
                    if (is4 > 0)
                        ndif = 2;
                } else if (ielem <= ntet + npyr + npri) {
                    kk = 6;
                    ind2 = ntet * 4 + npyr * 5 + (ielem - ntet - npyr - 1) * 6;
                    ind2f = ntet * 4 + npyr * 5 + (ielem - ntet - npyr - 1) * 5;
                    ndif = 4;
                    if (is4 > 0)
                        ndif = 3;
                } else {
                    if (mnodes == 4)
                        continue;  // hex can't be adjacent to tet
                    if (is4 == 0)
                        continue;  // tri face can't belong to hex
                    kk = 8;
                    ind2 = ntet * 4 + npyr * 5 + npri * 6 +
                           (ielem - ntet - npyr - npri - 1) * 8;
                    ind2f = ntet * 4 + npyr * 5 + npri * 5 +  // LATHOS??? to 5
                            (ielem - ntet - npyr - npri - 1) * 6;
                    ndif = 5;
                }
                kount = 0;
                for (int m1 = 1; m1 <= kk; m1++)  // Loop on the kk nodes of the element
                {
                    k1 = abs(nu[ind2 + m1]);
                    k2 = 1;
                    if (k1 == is1 || k1 == is2 || k1 == is3 || k1 == is4)
                        k2 = 0;
                    if (k2 != 0) {
                        kount++;
                        if (kount == ndif) {
                            ik++;
                            goto LOOP42;
                        }
                        nod[kount] = m1;
                    }
                }
                if (kount < 1) {
                    cout << " numfaces: Error, with data !" << endl;
                    exit;
                }
                if (kk == 5 && nod[1] == 3)
                    nod[kount] = 1;
                if (kk == 6 && nod[1] == 3)
                    nod[kount] = 1;
                if (kk == 6 && nod[1] == 4)
                    nod[kount] = 2;
                if (nod[kount] == 8 && nod[1] == 1)
                    nod[kount] = 1;
                if (nod[kount] == 7 && nod[1] == 2)
                    nod[kount] = 2;
                if (nod[kount] == 8 && nod[1] == 3)
                    nod[kount] = 3;
                if (nod[kount] == 8 && nod[1] == 5)
                    nod[kount] = 5;
                nu[ind2 + nod[kount]] = -nu[ind2 + nod[kount]];  // for this also,nu turns - ve

                nbfac--;  // This was not a boundary face
                break;
            }
        }
        ind += mnodes;
    }

    /*Allocate mem for the other boundary face arrays
	============================================== =*/
    vnofac = matrix<double>(3, nbfac);
    listbf.resize(nbfac + 1);
    /*vnofac = new double*[3];
	for (unsigned int i = 0; i < 3; i++)
	{
	vnofac[i] = new double[nbfac];
	}*/

    /*Reconstitute positive NU
	========================*/

    indp = 0;
    // for tetrahedra
    for (unsigned int i = 1; i <= ntet; i++) {
        for (unsigned int j = 1; j <= 4; j++) {
            indp++;
            nu[indp] = -nu[indp];
            if (nu[indp] <= 0) {
                cout << " Numfaces: Error with data at tetrahedra" << endl;
                // exit;
            }
        }
    }
    // for pyramids
    for (unsigned int i = 1; i <= npyr; i++) {
        for (unsigned int j = 1; j <= 5; j++) {
            indp++;
            nu[indp] = -nu[indp];
            if (nu[indp] <= 0) {
                cout << " Numfaces: Error with data at pyramids" << endl;
                // exit;
            }
        }
    }
    // for prisms
    for (unsigned int i = 1; i <= npri; i++) {
        for (unsigned int j = 1; j <= 6; j++) {
            indp++;
            nu[indp] = -nu[indp];
            if (nu[indp] <= 0) {
                cout << " Numfaces: Error with data at prisms" << endl;
                // exit;
            }
        }
    }

    // for hexaedra
    for (unsigned int i = 1; i <= nhex; i++) {
        for (unsigned int j = 1; j <= 8; j++) {
            indp++;
            nu[indp] = -nu[indp];
            if (nu[indp] <= 0) {
                cout << " Numfaces: Error with data at hexaedra" << endl;
                // exit;
            }
        }
    }
    return;
}

void DataStructures::filistbf() {
    qualityCheck gcheck;
    vector<int> itype;
    itype.resize(nbfac + 1);
    int is1, is2, is3, is4, logfr1, logfr2, logfr3, logfr4,
        ityp, ifac, i1, i2, i3, i4;
    double x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4,
        resx, resy, resz;

    cout << " Bound. faces are listed as follows: " << endl;

    for (unsigned int i = 0; i <= maxlist; i++) {
        listbf1[i] = 1;
        listbf2[i] = 0;
    }
    int kount, kounter, kount1 = 0, kount2 = 0, kount3 = 0, kount4 = 0, kount5 = 0, kount6 = 0;

    //     Sweep boundary faces & and analyse them

    for (unsigned int i = 1; i <= nbfac; i++) {
        is1 = nubf[0][i - 1];
        is2 = nubf[1][i - 1];
        is3 = nubf[2][i - 1];
        is4 = nubf[3][i - 1];
        logfr1 = logfr[is1];
        x1 = coor[0][is1 - 1];
        y1 = coor[1][is1 - 1];
        z1 = coor[2][is1 - 1];
        logfr2 = logfr[is2];
        x2 = coor[0][is2 - 1];
        y2 = coor[1][is2 - 1];
        z2 = coor[2][is2 - 1];
        logfr3 = logfr[is3];
        x3 = coor[0][is3 - 1];
        y3 = coor[1][is3 - 1];
        z3 = coor[2][is3 - 1];
        if (is4 == 0) {
            logfr4 = -9;  // 3 node face
            x4 = y4 = z4 = 0.;
        } else {
            logfr4 = logfr[is4];
            x4 = coor[0][is4 - 1];
            y4 = coor[1][is4 - 1];
            z4 = coor[2][is4 - 1];
        }
        ityp = 0;
        ityp = analface(logfr1, logfr2, logfr3, logfr4, ityp);

        if (ityp == 1)
            kount1++;
        else if (ityp == 2)
            kount2++;
        else if (ityp == 3)
            kount3++;
        else if (ityp == 4)
            kount4++;
        else if (ityp == 5)
            kount5++;
        else if (ityp == 6)
            kount6++;
        else {
            cout << " Unknown boundary face type returned " << ityp << endl;
            cout << " Nodes :" << is1 << ", " << is2 << ", " << is3 << ", " << is4 << endl;
            cout << " x-coor :" << x1 << ", " << x2 << ", " << x3 << ", " << x4 << endl;
            cout << " y-coor :" << y1 << ", " << y2 << ", " << y3 << ", " << y4 << endl;
            cout << " z-coor :" << z1 << ", " << z2 << ", " << z3 << ", " << z4 << endl;
            exit;
        }
        itype[i] = ityp;
    }
    cout << kount1 << " faces are PERIODIC" << endl;
    cout << kount2 << " faces are SYMMETRY" << endl;
    cout << kount3 << " faces are SOLID WALL" << endl;
    cout << kount4 << " faces are INLET" << endl;
    cout << kount5 << " faces are OUTLET" << endl;
    cout << kount6 << " faces are ROTATING WALL" << endl;
    cout << "-------------------------------------------" << endl;
    cout << " TOTAL NUMBER OF BOUNDARY FACES: " << nbfac << endl;
    cout << "-------------------------------------------" << endl;

    // Check if all faces are analysed

    kount = kount1 + kount2 + kount3 + kount4 + kount5 + kount6;
    if (kount != nbfac) {
        cout << " error in filistbf " << endl;
        exit;
    }

    //      FILL THE LIST
    //     ================
    for (unsigned int ktype = 1; ktype <= 6; ktype++) {
        if (ktype == 1) {
            kounter = 0;
            listbf1[1] = 1;
        } else {
            listbf1[ktype] = listbf2[ktype - 1] + 1;
            kounter = listbf2[ktype - 1];
        }

        for (unsigned int ifac = 1; ifac <= nbfac; ifac++) {
            ityp = itype[ifac];
            if (ityp == ktype) {
                kounter++;
                listbf[kounter] = ifac;
            }
        }

        listbf2[ktype] = kounter;
        if (listbf2[ktype] < listbf1[ktype]) {
            listbf2[ktype] = listbf2[ktype - 1];
        } else {
            cout << " Faces of type " << ktype << " :between :" << listbf1[ktype] << " " << listbf2[ktype] << endl;
        }
    }
    cout << " ------  ------  ------  ------  ------ " << endl;
    if (kounter > nbfac) {
        cout << " Filistbf: ERROR" << endl;
        exit;
    }
    // STOP AT LINE 1255
    if (listnp2[1] > listnp1[1]) {
        perfac();
        cout << " Perfac completed...." << endl;  // correct periodic lists - find special faces
    }
    /* ----------------------------------------------------
	Find Normal Vectors for Boundary faces (outwards)
	----------------------------------------------------
	*/
    for (unsigned int ibfac = listbf1[1]; ibfac <= listbf2[6]; ibfac++) {
        ifac = listbf[ibfac];
        // node1
        i1 = nubf[0][ifac - 1];
        x1 = coor[0][i1 - 1];
        y1 = coor[1][i1 - 1];
        z1 = coor[2][i1 - 1];
        // node2
        i2 = nubf[1][ifac - 1];
        x2 = coor[0][i2 - 1];
        y2 = coor[1][i2 - 1];
        z2 = coor[2][i2 - 1];
        // node3
        i3 = nubf[2][ifac - 1];
        x3 = coor[0][i3 - 1];
        y3 = coor[1][i3 - 1];
        z3 = coor[2][i3 - 1];
        // node4
        i4 = nubf[3][ifac - 1];

        if (i4 == 0)
        //  The normal/area of this 3-node face is:
        {
            gcheck.vecProd(x2 - x1, y2 - y1, z2 - z1, x3 - x1, y3 - y1, z3 - z1,
                           resx, resy, resz);
            vnofac[0][ibfac - 1] = resx * 0.5;
            vnofac[1][ibfac - 1] = resy * 0.5;
            vnofac[2][ibfac - 1] = resz * 0.5;
        } else {
            x4 = coor[0][i4 - 1];
            y4 = coor[1][i4 - 1];
            z4 = coor[2][i4 - 1];
            // The normal/area (approx.) of this 4-node face is:
            gcheck.vecProd(x3 - x1, y3 - y1, z3 - z1, x4 - x2, y4 - y2, z4 - z2,
                           resx, resy, resz);

            vnofac[0][ibfac - 1] = resx * 0.5;
            vnofac[1][ibfac - 1] = resy * 0.5;
            vnofac[2][ibfac - 1] = resz * 0.5;
        }
    }
    return;
}

void DataStructures::perfac() {
    /*
	locates boundary faces falsely marked as periodic by
	inding correspondance of periodic faces
	(so it is also a useful check for the periodicity)
	provisionally uses jaret
	*/

    vector<int> ik, jk;
    ik.resize(5);
    jk.resize(5);
    int inod1, inod2, ifac, jfac, kou, kou2, kount;
    int lgf10, lgf20, lgf30, lgf40, ktipos, ik1, ik2, ik3, ik4,
        iempty, irep, jdum2;

    for (unsigned int i = 1; i <= ns; i++)
        jaret[i] = -1;

    for (unsigned int i = 1; i <= nper; i++) {
        inod1 = iper[0][i - 1];
        inod2 = iper[1][i - 1];
        jaret[inod1] = inod2;
        jaret[inod2] = inod1;
    }
    int nperfac = 0;
    int iexit = listbf2[1];  // listbf2 will be altered
    // so keep it here for the do-enddo

    for (unsigned int idum = listbf1[1]; idum <= iexit; idum++) {
    LINE99:
        ifac = listbf[idum];
        if (ifac < 0)  // pair already treated
        {
            listbf[idum] = ifac;  // restore it
            continue;
        }
        if (idum > listbf2[1])
            exit;  // past all periodic faces
        ik[1] = nubf[0][ifac - 1];
        ik[2] = nubf[1][ifac - 1];
        ik[3] = nubf[2][ifac - 1];
        ik[4] = nubf[3][ifac - 1];
        kou = 4;
        if (ik[4] == 0)
            kou = 3;
        // pick its periodic nodes
        for (unsigned int i = 1; i <= kou; i++) {
            ik[i] = jaret[ik[i]];
            if (ik[i] == -1) {
                cout << "Problem with periodic nodes in perfac" << endl;
                exit;
            }
        }

        for (unsigned int jdum = idum + 1; jdum <= listbf2[1]; jdum++) {
            jfac = listbf[jdum];  // bf numbering
            if (listbf[jfac] < 0)
                continue;  // already matched
            jk[1] = nubf[0][jfac - 1];
            jk[2] = nubf[1][jfac - 1];
            jk[3] = nubf[2][jfac - 1];
            jk[4] = nubf[3][jfac - 1];
            kou2 = 4;
            if (jk[4] == 0)
                kou2 = 3;
            if (kou != kou2)
                continue;  // different face type
            kount = 0;
            for (unsigned int m1 = 1; m1 <= kou; m1++) {
                for (unsigned int m2 = 1; m2 <= kou; m2++) {
                    if (ik[m1] - jk[m2] == 0)
                        kount++;
                }
            }
            if (kount == kou) {
                listbf[jdum] = -listbf[jdum];
                nperfac += 1;
                exit;
            }
            jdum2 = jdum;
        }
        if (listbf[jdum2] < 0)
            continue;
        // check if it is a 'special' face and compress lists accordingly
        ik1 = nubf[0][ifac - 1];
        ik2 = nubf[1][ifac - 1];
        ik3 = nubf[2][ifac - 1];
        ik4 = nubf[3][ifac - 1];

        cout << endl;
        cout << "No periodic found for face " << ifac << endl;
        cout << " nodes(f3d numbering): " << ik1 << ", " << ik2
             << ", " << ik3 << ", " << ik4 << endl;
        cout << " Checking if it is a ''special'' face" << endl;
        cout << endl;
        lgf10 = logfr[ik1];
        lgf20 = logfr[ik2];
        lgf30 = logfr[ik3];
        if (ik4 == 0) {
            if (lgf10 == lgf20 && lgf10 == lgf30 && lgf10 == 12) {
                ktipos = 2;
            } else if (lgf10 == lgf20 && lgf10 == lgf30 && lgf10 == 13) {
                ktipos = 3;
            } else if (lgf10 == lgf20 && lgf10 == lgf30 && lgf10 == 16) {
                ktipos = 6;
            } else {
                cout << " Unknown boundary face type returned " << endl;
                cout << " Nodes : " << ik1 << ", " << ik2
                     << ", " << ik3 << ", " << endl;
                cout << " Logfrs: " << lgf10 << ", " << lgf20 << ", " << lgf30 << endl;
                exit;
            }
        } else  // quad face
        {
            lgf40 = logfr[ik4];
            if (lgf10 == lgf20 && lgf10 == lgf30 &&
                lgf10 == lgf40 && lgf10 == 12) {
                ktipos = 2;
            } else if (lgf10 == lgf20 && lgf10 == lgf30 &&
                       lgf10 == lgf40 && lgf10 == 13) {
                ktipos = 3;
            } else if (lgf10 == lgf20 && lgf10 == lgf30 &&
                       lgf10 == lgf40 && lgf10 == 16) {
                ktipos = 6;
            } else {
                cout << " Unknown boundary face type returned " << endl;
                cout << " Nodes : " << ik1 << ", " << ik2
                     << ", " << ik3 << ", " << ik4 << endl;
                cout << " Logfrs: " << lgf10 << ", " << lgf20 << ", "
                     << lgf30 << ", " << lgf40 << endl;
                exit;
            }
        }
        iempty = idum;  //  first empty pos - erased face
        for (unsigned int ipos = 1; ipos <= ktipos; ipos++) {
            irep = listbf2[ipos];
            listbf[iempty] = listbf[irep];
            iempty = listbf2[ipos];
            if (ipos != 1)
                listbf1[ipos] = listbf1[ipos] - 1;
            if (ipos != ktipos)
                listbf2[ipos] = listbf2[ipos] - 1;
        }
        listbf[iempty] = ifac;  // put face to its correct position
        cout << "OK" << endl
             << endl;
        goto LINE99;  // repeat for current position, it is now another face
    }
    if (2 * nperfac != listbf2[1] - listbf1[1] + 1) {
        cout << "error in perfac" << endl;
        exit;
    }
    return;
}

int DataStructures::analface(int lgf10, int lgf20, int lgf30, int lgf40, int ktipos) {
    vector<int> kdig, lgf;
    kdig.resize(7);
    lgf.resize(5);
    int kkk, m010, m100, m001;
    //
    int nnode = 4;  // default 4 node (quadrangle) face
    if (lgf40 < 0)  // else 3 node (triangle) face
    {
        nnode = 3;
        if (lgf10 == lgf20 && lgf10 == lgf30 && lgf10 > 10) {
            //fast 'obvious' classification
            if (lgf10 == 34)
                ktipos = 4;
            if (lgf10 == 35)
                ktipos = 5;
            if (lgf10 == 13 || lgf10 == 16)
                ktipos = 1;
            if (lgf10 == 24)
                ktipos = 4;
            if (lgf10 == 25)
                ktipos = 5;
            return ktipos;
        }
    }
    if (lgf10 == lgf20 && lgf10 == lgf30 && lgf10 == lgf40 && lgf10 > 10) {
        //fast 'obvious' classification
        if (lgf10 == 34)
            ktipos = 4;
        if (lgf10 == 35)
            ktipos = 5;
        if (lgf10 == 13 || lgf10 == 16)
            ktipos = 1;
        if (lgf10 == 24)
            ktipos = 4;
        if (lgf10 == 25)
            ktipos = 5;
        return ktipos;
    }
    //
    //	zeroing kdig
    for (int i = 1; i <= 6; i++) {
        kdig[i] = 0;
    }
    //
    lgf[1] = lgf10;
    lgf[2] = lgf20;
    lgf[3] = lgf30;
    lgf[4] = lgf40;
    //
    // Sweep the nodes of the face
    for (int m = 1; m <= nnode; m++) {
        kkk = lgf[m];
        m100 = kkk / 100;
        if (m100 != 0)
            kdig[m100] = kdig[m100] + 1;
        kkk = kkk - m100 * 100;
        m010 = kkk / 10;
        if (m010 != 0)
            kdig[m010] = kdig[m010] + 1;
        m001 = kkk - m010 * 10;
        if (m001 != 0)
            kdig[m001] = kdig[m001] + 1;
    }
    ktipos = -999;

    ktipos = 0;
    //	find face type
    for (int i = 1; i < 6; i++) {
        if (kdig[i] == nnode) {
            if (ktipos != 0) {
                cout << "SPECIAL FACE OF TYPE	" << ktipos << "AND " << i << endl;
            }
            ktipos = i;
        }
    }
    //
    return ktipos;
}

void DataStructures::calvnocl() {
    qualityCheck gcheck;
    double po1, po2, signp, vnox1, vnoy1, vnoz1, vnox2, vnoy2, vnoz2,
        conx1, cony1, conx2, cony2;
    int kk1, kk2, iseg1, iseg2;
    int ind, ind2, inod, iseg, inod1, inod2, isegloc, isegloc1, isegloc2, kk,
        inodloc1, inodloc2;
    double barx, bary, barz, barfx, barfy, barfz, xmseg, ymseg, zmseg, resx,
        resy, resz, desirx, desiry, desirz, pinner, aind, coe, is, js, cm, sm;

    int iaux1[6][2] = {{1, 2}, {1, 3}, {1, 4},  // elements type 1 --> tetrahedra
                       {2, 3},
                       {2, 4},
                       {3, 4}};

    int iaux2[8][2] = {{1, 2}, {1, 4}, {1, 5},  // elements type 2 --> pyramids
                       {2, 3},
                       {2, 5},
                       {3, 4},
                       {3, 5},
                       {4, 5}};

    int iaux3[9][2] = {{1, 2}, {1, 3}, {1, 4},  // elements type 3 --> prisms
                       {2, 3},
                       {2, 5},
                       {3, 6},
                       {4, 5},
                       {4, 6},
                       {5, 6}};

    int iaux4[12][2] = {{1, 2}, {1, 4}, {1, 5},  // elements type 4 --> hexahedra
                        {2, 3},
                        {2, 6},
                        {3, 4},
                        {3, 7},
                        {4, 8},
                        {5, 6},
                        {5, 8},
                        {6, 7},
                        {7, 8}};

    int iauxs1[4][3] = {{4, 5, 6}, {2, 3, 6},  // elements type 1 --> tetrahedra
                        {1, 3, 5},
                        {1, 2, 4}};

    int iauxs2[5][4] = {{1, 3, 5, 0}, {6, 7, 8, 0}, {2, 3, 8, 0},  // elements type 2 --> pyramis
                        {4, 5, 7, 0},
                        {1, 2, 4, 6}};

    int iauxs3[5][4] = {{1, 3, 5, 7}, {1, 2, 4, 0}, {7, 8, 9, 0},  // elements type 3 --> prisms
                        {4, 5, 6, 9},
                        {2, 3, 6, 8}};

    int iauxs4[6][4] = {{4, 5, 7, 11}, {2, 3, 8, 10}, {1, 3, 5, 9},  // elements type 4 --> hexahedra
                        {9, 10, 11, 12},
                        {1, 2, 4, 6},
                        {6, 7, 8, 12}};
    //
    double us6 = 1. / 6.;
    double us8 = 1. / 8.;

    // Initialization
    for (int i = 1; i <= nseg; i++) {
        vnocl[0][i - 1] = 0.;
        vnocl[1][i - 1] = 0.;
        vnocl[2][i - 1] = 0.;
    }
    //
    //
    // VNOCL always points NUBO1
    ind = 0;   // nodal index
    ind2 = 0;  // edge index
    //  Sweep tetrahedra
    for (int ip = 1; ip <= ntet; ip++)  // 10
    {
        // Find the Tetrahedron Barycentre
        barx = 0.;
        bary = 0.;
        barz = 0.;
        for (int k = 1; k <= 4; k++) {
            inod = nu[ind + k];
            barx += coor[0][inod - 1];  // exw valei -1
            bary += coor[1][inod - 1];  // exw valei -1
            barz += coor[2][inod - 1];  // exw valei -1
        }
        barx *= 0.25;
        bary *= 0.25;
        barz *= 0.25;
        //
        for (int ifac = 1; ifac <= 4; ifac++)  // 11 // the four faces of the tetrahedron
        {
            //			Find the Face Barycentre
            barfx = 0.;
            barfy = 0.;
            barfz = 0.;
            for (int k = 1; k <= 3; k++)  // three segments
            {
                iseg = nusg[ind2 + iauxs1[ifac - 1][k - 1]];
                inod1 = nubo[0][iseg - 1];
                inod2 = nubo[1][iseg - 1];
                barfx = barfx + coor[0][inod1 - 1] + coor[0][inod2 - 1];
                barfy = barfy + coor[1][inod1 - 1] + coor[1][inod2 - 1];
                barfz = barfz + coor[2][inod1 - 1] + coor[2][inod2 - 1];
            }
            barfx *= us6;
            barfy *= us6;
            barfz *= us6;
            //
            //			Sweep the 3 segments and scatter-added VNOCL,
            //			cell contributions
            for (int isegdum = 1; isegdum <= 3; isegdum++) {
                isegloc = iauxs1[ifac - 1][isegdum - 1];
                iseg = nusg[ind2 + isegloc];
                inod1 = nubo[0][iseg - 1];
                inod2 = nubo[1][iseg - 1];
                inodloc1 = iaux1[isegloc - 1][0];
                inodloc2 = iaux1[isegloc - 1][1];
                xmseg = 0.5 * (coor[0][inod1 - 1] + coor[0][inod2 - 1]);
                ymseg = 0.5 * (coor[1][inod1 - 1] + coor[1][inod2 - 1]);
                zmseg = 0.5 * (coor[2][inod1 - 1] + coor[2][inod2 - 1]);
                gcheck.vecProd(barx - xmseg, bary - ymseg, barz - zmseg,
                               barfx - xmseg, barfy - ymseg, barfz - zmseg,
                               resx, resy, resz);
                desirx = coor[0][inod1 - 1] - xmseg;
                desiry = coor[1][inod1 - 1] - ymseg;
                desirz = coor[2][inod1 - 1] - zmseg;
                pinner = resx * desirx + resy * desiry + resz * desirz;
                aind = 1.;
                if (nu[ind + inodloc1] == inod2)
                    aind = -1.;
                if (pinner > 0.) {
                    // segment - based scatter-adding for convective fluxes
                    vnocl[0][iseg - 1] = vnocl[0][iseg - 1] + 0.5 * resx;
                    vnocl[1][iseg - 1] = vnocl[1][iseg - 1] + 0.5 * resy;
                    vnocl[2][iseg - 1] = vnocl[2][iseg - 1] + 0.5 * resz;
                    // element-based scatter-adding for diffusive fluxes
                } else {
                    //change sign
                    // segment-based scatter-adding for convective fluxes
                    vnocl[0][iseg - 1] = vnocl[0][iseg - 1] - 0.5 * resx;
                    vnocl[1][iseg - 1] = vnocl[1][iseg - 1] - 0.5 * resy;
                    vnocl[2][iseg - 1] = vnocl[2][iseg - 1] - 0.5 * resz;
                    // element-based scatter-adding for diffusive fluxes
                    pinner = -pinner;  // it is now a +ve volume
                }
                cell[inod1] = cell[inod1] + pinner * us6;
                cell[inod2] = cell[inod2] + pinner * us6;
            }  // isegdum
        }
        ind += 4;
        ind2 += 6;
    }
    //	-----------
    //
    //	Sweep pyramids
    for (int ip = ntet + 1; ip <= ntet + npyr; ip++)  //20
    {
        // Find the Pyramid Barycentre
        barx = 0.;
        bary = 0.;
        barz = 0.;
        for (int k = 1; k < 5; k++) {
            inod = nu[ind + k];
            barx = barx + coor[0][inod - 1];
            bary = bary + coor[1][inod - 1];
            barz = barz + coor[2][inod - 1];
        }
        barx *= 0.2;
        bary *= 0.2;
        barz *= 0.2;

        for (int ifac = 1; ifac <= 5; ifac++)  //21 // the five faces
        {
            // Find the Face Barycentre
            barfx = 0.;
            barfy = 0.;
            barfz = 0.;
            kk = 4;
            if (iauxs2[ifac - 1][3] == 0)
                kk = 3;  // it is a triangle face
            for (int k = 1; k <= kk; k++) {
                iseg = nusg[ind2 + iauxs2[ifac - 1][k - 1]];
                inod1 = nubo[0][iseg - 1];
                inod2 = nubo[1][iseg - 1];
                barfx = barfx + coor[0][inod1 - 1] + coor[0][inod2 - 1];
                barfy = barfy + coor[1][inod1 - 1] + coor[1][inod2 - 1];
                barfz = barfz + coor[2][inod1 - 1] + coor[2][inod2 - 1];
            }
            coe = 2. * kk;
            barfx /= coe;
            barfy /= coe;
            barfz /= coe;
            //
            //		Sweep the kk segments and scatter-added VNOCL
            //		cell contributions
            for (int isegdum = 1; isegdum <= kk; isegdum++) {
                isegloc = iauxs2[ifac - 1][isegdum - 1];
                iseg = nusg[ind2 + isegloc];
                inod1 = nubo[0][iseg - 1];
                inod2 = nubo[1][iseg - 1];
                inodloc1 = iaux2[isegloc - 1][0];
                inodloc2 = iaux2[isegloc - 1][1];
                xmseg = 0.5 * (coor[0][inod1 - 1] + coor[0][inod2 - 1]);
                ymseg = 0.5 * (coor[1][inod1 - 1] + coor[1][inod2 - 1]);
                zmseg = 0.5 * (coor[2][inod1 - 1] + coor[2][inod2 - 1]);
                gcheck.vecProd(barx - xmseg, bary - ymseg, barz - zmseg,
                               barfx - xmseg, barfy - ymseg, barfz - zmseg,
                               resx, resy, resz);
                desirx = coor[0][inod1 - 1] - xmseg;
                desiry = coor[1][inod1 - 1] - ymseg;
                desirz = coor[2][inod1 - 1] - zmseg;
                pinner = resx * desirx + resy * desiry + resz * desirz;
                aind = 1.;
                if (nu[ind + inodloc1] == inod2)
                    aind = -1.;
                if (pinner > 0.) {
                    // segment-based scatter-adding for convective fluxes
                    vnocl[0][iseg - 1] = vnocl[0][iseg - 1] + 0.5 * resx;
                    vnocl[1][iseg - 1] = vnocl[1][iseg - 1] + 0.5 * resy;
                    vnocl[2][iseg - 1] = vnocl[2][iseg - 1] + 0.5 * resz;
                } else  // change sign
                {
                    // segment-based scatter-adding for convective fluxes
                    vnocl[0][iseg - 1] = vnocl[0][iseg - 1] - 0.5 * resx;
                    vnocl[1][iseg - 1] = vnocl[1][iseg - 1] - 0.5 * resy;
                    vnocl[2][iseg - 1] = vnocl[2][iseg - 1] - 0.5 * resz;
                    pinner = -pinner;
                }
                cell[inod1] = cell[inod1] + pinner * us6;
                cell[inod2] = cell[inod2] + pinner * us6;
            }
        }
        ind += 5;
        ind2 += 8;
    }
    //-----------------
    //
    // Sweep prisms
    for (int ip = ntet + npyr + 1; ip <= ntet + npyr + npri; ip++)  //30
    {
        // Find the Prism Barycentre
        barx = 0.;
        bary = 0.;
        barz = 0.;
        for (int k = 1; k <= 6; k++) {
            inod = nu[ind + k];
            barx = barx + coor[0][inod - 1];
            bary = bary + coor[1][inod - 1];
        }
        barx *= us6;
        bary *= us6;
        barz *= us6;

        for (int ifac = 1; ifac <= 5; ifac++)  // 31
        {
            // Find the Face Barycentre
            barfx = 0.;
            barfy = 0.;
            barfz = 0.;
            kk = 4;
            if (iauxs3[ifac - 1][3] == 0)
                kk = 3;                    // it is a triangle face
            for (int k = 1; k <= kk; k++)  // for the kk segments of the face
            {
                iseg = nusg[ind2 + iauxs3[ifac - 1][k - 1]];
                inod1 = nubo[0][iseg - 1];
                inod2 = nubo[1][iseg - 1];
                barfx = barfx + coor[0][inod1 - 1] + coor[0][inod2 - 1];
                barfy = barfy + coor[1][inod1 - 1] + coor[1][inod2 - 1];
                barfz = barfz + coor[2][inod1 - 1] + coor[2][inod2 - 1];
            }
            coe = 2. * kk;
            barfx /= coe;
            barfy /= coe;
            barfz /= coe;
            //
            //		Sweep the kk segmens and scatter-added VNOCL,
            //		cell contributions
            for (int isegdum = 1; isegdum <= kk; isegdum++) {
                isegloc = iauxs3[ifac - 1][isegdum - 1];
                iseg = nusg[ind2 + isegloc];
                inod1 = nubo[0][iseg - 1];
                inod2 = nubo[1][iseg - 1];
                inodloc1 = iaux3[isegloc - 1][0];
                inodloc2 = iaux3[isegloc - 1][1];
                xmseg = 0.5 * (coor[0][inod1 - 1] + coor[0][inod2 - 1]);
                ymseg = 0.5 * (coor[1][inod1 - 1] + coor[1][inod2 - 1]);
                zmseg = 0.5 * (coor[2][inod1 - 1] + coor[2][inod2 - 1]);
                gcheck.vecProd(barx - xmseg, bary - ymseg, barz - zmseg,
                               barfx - xmseg, barfy - ymseg, barfz - zmseg,
                               resx, resy, resz);
                desirx = coor[0][inod1 - 1] - xmseg;
                desiry = coor[1][inod1 - 1] - ymseg;
                desirz = coor[2][inod1 - 1] - zmseg;
                pinner = resx * desirx + resy * desiry + resz * desirz;
                aind = 1.;
                if (nu[ind + inodloc1] == inod2)
                    aind = -1.;
                if (pinner > 0.) {
                    // segment-based scatter-adding for convective fluxes
                    vnocl[0][iseg - 1] = vnocl[0][iseg - 1] + 0.5 * resx;
                    vnocl[1][iseg - 1] = vnocl[1][iseg - 1] + 0.5 * resy;
                    vnocl[2][iseg - 1] = vnocl[2][iseg - 1] + 0.5 * resz;
                } else {
                    // segment-based scatter-adding for convective fluxes
                    vnocl[0][iseg - 1] = vnocl[0][iseg - 1] - 0.5 * resx;
                    vnocl[1][iseg - 1] = vnocl[1][iseg - 1] - 0.5 * resy;
                    vnocl[2][iseg - 1] = vnocl[2][iseg - 1] - 0.5 * resz;
                    pinner = -pinner;  // it is now a +ve volume
                }
                cell[inod1] = cell[inod1] + pinner * us6;
                cell[inod2] = cell[inod2] + pinner * us6;
            }  // isegdum
        }
        ind += 6;
        ind2 += 9;
    }
    //-------------------
    //
    // Sweep hexahedra
    for (int ip = ntet + npyr + npri + 1; ip <= ntet + npyr + npri + nhex; ip++)  //40
    {
        // Find the Hexahedron Barycentre
        barx = 0.;
        bary = 0.;
        barz = 0.;
        for (int k = 1; k <= 8; k++) {
            inod = nu[ind + k];
            barx = barx + coor[0][inod - 1];
            bary = bary + coor[1][inod - 1];
            barz = barz + coor[2][inod - 1];
        }
        barx *= us8;
        bary *= us8;
        barz *= us8;

        for (int ifac = 1; ifac <= 6; ifac++)  //41 // the six faces
        {
            // Find the Face Barycentre
            barfx = 0.;
            barfy = 0.;
            barfz = 0.;
            for (int k = 1; k <= 4; k++)  // for the 4 segments of the face
            {
                iseg = nusg[ind2 + iauxs4[ifac - 1][k - 1]];
                inod1 = nubo[0][iseg - 1];
                inod2 = nubo[1][iseg - 1];
                barfx = barfx + coor[0][inod1 - 1] + coor[0][inod2 - 1];
                barfy = barfy + coor[1][inod1 - 1] + coor[1][inod2 - 1];
                barfz = barfz + coor[2][inod1 - 1] + coor[2][inod2 - 1];
            }
            barfx *= us8;
            barfy *= us8;
            barfz *= us8;
            //
            //	Sweep the kk segments and scatter-added VNOCL
            //	cell contributions
            for (int isegdum = 1; isegdum <= 4; isegdum++) {
                isegloc = iauxs4[ifac - 1][isegdum - 1];
                iseg = nusg[ind2 + isegloc];
                inod1 = nubo[0][iseg - 1];
                inod2 = nubo[1][iseg - 1];
                inodloc1 = iaux4[isegloc - 1][0];
                inodloc2 = iaux4[isegloc - 1][1];
                xmseg = 0.5 * (coor[0][inod1 - 1] + coor[0][inod2 - 1]);
                ymseg = 0.5 * (coor[1][inod1 - 1] + coor[1][inod2 - 1]);
                zmseg = 0.5 * (coor[2][inod1 - 1] + coor[2][inod2 - 1]);
                gcheck.vecProd(barx - xmseg, bary - ymseg, barz - zmseg,
                               barfx - xmseg, barfy - ymseg, barfz - zmseg,
                               resx, resy, resz);
                desirx = coor[0][inod1 - 1] - xmseg;
                desiry = coor[1][inod1 - 1] - ymseg;
                desirz = coor[2][inod1 - 1] - zmseg;
                pinner = resx * desirx + resy * desiry + resz * desirz;
                aind = 1.;
                if (nu[ind + inodloc1] == inod2)
                    aind = -1.;
                if (pinner > 0.) {
                    // segment-based scatter-adding for convective fluxes
                    vnocl[0][iseg - 1] = vnocl[0][iseg - 1] + 0.5 * resx;
                    vnocl[1][iseg - 1] = vnocl[1][iseg - 1] + 0.5 * resy;
                    vnocl[2][iseg - 1] = vnocl[2][iseg - 1] + 0.5 * resz;
                } else  // change sign
                {
                    //segment-based scatter-adding for convective fluxes
                    vnocl[0][iseg - 1] = vnocl[0][iseg - 1] - 0.5 * resx;
                    vnocl[1][iseg - 1] = vnocl[1][iseg - 1] - 0.5 * resy;
                    vnocl[2][iseg - 1] = vnocl[2][iseg - 1] - 0.5 * resz;
                    pinner = -pinner;  // it is a +ve volume
                }
                cell[inod1] = cell[inod1] + pinner * us6;
                cell[inod2] = cell[inod2] + pinner * us6;
            }  // isegdum
        }
        ind += 8;
        ind2 += 12;
    }
    // ----------------
    //
    // Periodicity for array CELL
    if (nper > 0) {
        for (int i = 1; i <= nper; i++) {
            is = iper[0][i - 1];
            js = iper[1][i - 1];
            cell[is] = cell[is] + cell[js];
            cell[js] = cell[is];
        }
    }
    //
    //	Fix vnocl for periodic segments
    //	---------------------
    if (nper > 0)
    //	---------------------
    {
        if (isperiph == 1) {
            cm = pint;
            sm = pext;
        } else {
            cm = 1.;
            sm = 0.;
        }
        for (int ipsg = 1; ipsg <= npersg; ipsg++) {
            iseg1 = ipersg[0][ipsg - 1];
            iseg2 = ipersg[1][ipsg - 1];
            iseg1 = abs(iseg1);
            iseg2 = abs(iseg2);
            kk1 = 1;
            kk2 = 1;
            po1 = 0.5 * kk1;
            po2 = 0.5 * kk2;
            signp = ipersg[0][ipsg - 1];
            signp = copysign(1., signp);
            vnox1 = vnocl[0][iseg1 - 1] * signp;  // vnocl1 oriented as iseg2
            vnoy1 = vnocl[1][iseg1 - 1] * signp;
            vnoz1 = vnocl[2][iseg1 - 1] * signp;
            vnox2 = vnocl[0][iseg2 - 1] * signp;  // vnocl2 oriented as iseg1
            vnoy2 = vnocl[1][iseg2 - 1] * signp;
            vnoz2 = vnocl[2][iseg2 - 1] * signp;
            //
            conx1 = vnox1 * cm - vnoy1 * sm;  // 1 --> 2 Turning(oriented)
            cony1 = vnox1 * sm + vnoy1 * cm;
            conx2 = vnox2 * cm + vnoy2 * sm;  // 1 <-- 2 Turning(oriented)
            cony2 = -vnox2 * sm + vnoy2 * cm;
            //
            vnocl[0][iseg1 - 1] = po1 * (vnocl[0][iseg1 - 1] + conx2);
            vnocl[1][iseg1 - 1] = po1 * (vnocl[1][iseg1 - 1] + cony2);
            vnocl[2][iseg1 - 1] = po1 * (vnocl[2][iseg1 - 1] + vnoz2);
            vnocl[0][iseg2 - 1] = po2 * (vnocl[0][iseg2 - 1] + conx1);
            vnocl[1][iseg2 - 1] = po2 * (vnocl[1][iseg2 - 1] + cony1);
            vnocl[2][iseg2 - 1] = po2 * (vnocl[2][iseg2 - 1] + vnoz1);
            //
            //		turn them positive for further use
            ipersg[0][ipsg - 1] = iseg1;
            ipersg[1][ipsg - 1] = iseg2;
        }  // ipsg
           //
           // -----
    }
    // -----
    //
    return;
}

void DataStructures::numsegs3D() {
    int ind, inds, is1, is2, iseg, i1, i2, is_min, is_max, isaid,
        kpoin_1, kpoin_2;

    int iaux1[6][2] = {{1, 2}, {1, 3}, {1, 4},  // elements type 1 --> tetrahedra
                       {2, 3},
                       {2, 4},
                       {3, 4}};
    int iaux2[8][2] = {{1, 2}, {1, 4}, {1, 5}, {2, 3},  // elements type 2 --> pyramids
                       {2, 5},
                       {3, 4},
                       {3, 5},
                       {4, 5}};
    int iaux3[9][2] = {{1, 2}, {1, 3}, {1, 4}, {2, 3},  // elements type 3 --> prisms
                       {2, 5},
                       {3, 6},
                       {4, 5},
                       {4, 6},
                       {5, 6}};
    int iaux4[12][2] = {{1, 2}, {1, 4}, {1, 5}, {2, 3},  // elements type 4 --> hexahedra
                        {2, 6},
                        {3, 4},
                        {3, 7},
                        {4, 8},
                        {5, 6},
                        {5, 8},
                        {6, 7},
                        {7, 8}};
    //
    fjaret3D();
    //
    cout << "JARET length (Only Nodes)=	" << max_str;
    //
    //	Auxiliary storage in the second part of JARET
    for (int i = 1; i <= ns; i++) {
        jaret[ndeg[i] + max_str] = ndeg[i - 1];  // Provisory Pointer is the
    }                                            // Last Seat of Previous node
    // ---------------------------------
    // Compute JARET (edges), NUSG, NUBO
    // ---------------------------------
    for (int j = 1; j <= maxseg; j++) {
        for (int i = 1; i <= 2; i++) {
            nubo[i - 1][j - 1] = 0;  // init set to 0
        }
    }

    nseg = 0;
    ind = 0;   // node index
    inds = 0;  // edge index
    //	Sweep all tetrahedra
    for (int ip = 1; ip <= ntet; ip++)  //10
    {
        for (int isg = 1; isg <= 6; isg++)  // 11 // sweep tetrahedra segments
        {
            is1 = nu[ind + iaux1[isg - 1][0]];
            is2 = nu[ind + iaux1[isg - 1][1]];
            for (int kv = ndeg[is1 - 1] + 1; kv <= ndeg[is1]; kv++)  // neis of 1
            {
                iseg = jaret[max_str + kv];
                if (iseg == 0 || iseg > nseg)
                    continue;
                if (iseg > maxseg) {
                    cout << "numsegs 1: increase MAXSEG" << endl;
                    exit(0);
                }
                i1 = nubo[0][iseg - 1];
                i2 = nubo[1][iseg - 1];
                is_min = min(is1, is2);
                is_max = max(is1, is2);
                if (is_min == i1 && is_max == i2)  // This Segments Exists
                {
                    isaid = iseg;
                    goto LINE12;
                }
            }
            //		this is a new segment (no-name):
            nseg = nseg + 1;
            if (nseg > maxseg) {
                cout << "numsegs 2: increase MAXSEG" << endl;
                exit(0);
            }
            isaid = nseg;
            kpoin_1 = jaret[max_str + ndeg[is1]] + 1;
            jaret[max_str + ndeg[is1]] = kpoin_1;
            jaret[kpoin_1 + max_str] = nseg;
            kpoin_2 = jaret[max_str + ndeg[is2]] + 1;
            jaret[max_str + ndeg[is2]] = kpoin_2;
            jaret[kpoin_2 + max_str] = nseg;
            nubo[0][nseg - 1] = min(is1, is2);
            nubo[1][nseg - 1] = max(is1, is2);
        LINE12:
            nusg[inds + isg] = isaid;
        }
        ind += 4;
        inds += 6;
    }
    //
    //	Sweep all pyramids
    for (int ip = ntet + 1; ip <= ntet + npyr; ip++)  //20
    {
        for (int isg = 1; isg <= 8; isg++)  // 21 //Sweep segments of the pyramids
        {
            is1 = nu[ind + iaux2[isg - 1][0]];
            is2 = nu[ind + iaux2[isg - 1][1]];
            for (int kv = ndeg[is1 - 1]; kv <= ndeg[is1]; kv++)  // neis of 1
            {
                iseg = jaret[max_str + kv];
                if (iseg == 0. || iseg > nseg)
                    continue;
                if (iseg == maxseg) {
                    cout << "numsegs 3: increase MAXSEG" << endl;
                    exit(0);
                }
                i1 = nubo[0][iseg - 1];
                i2 = nubo[1][iseg - 1];
                is_min = min(is1, is2);
                is_max = max(is1, is2);
                if (is_min == i1 && is_max == i2)  // This Segments Exists
                {
                    isaid = iseg;  //jaret[max_str+kv]
                    goto LINE22;
                }
            }
            // this is a new segment (no-name):
            nseg = nseg + 1;
            if (nseg > maxseg) {
                cout << "numsegs 4: increase MAXSEG" << endl;
                exit(0);
            }
            isaid = nseg;
            kpoin_1 = jaret[max_str + ndeg[is1]] + 1;
            jaret[max_str + ndeg[is1]] = kpoin_1;
            jaret[kpoin_1 + max_str] = nseg;
            kpoin_2 = jaret[max_str + ndeg[is2]] + 1;
            jaret[max_str + ndeg[is2]] = kpoin_2;
            jaret[kpoin_2 + max_str] = nseg;
            nubo[0][nseg - 1] = min(is1, is2);
            nubo[1][nseg - 1] = max(is1, is2);
        LINE22:
            nusg[inds + isg] = isaid;
        }
        ind += 5;
        inds += 8;
    }
    //
    //	Sweep all prisms
    for (int ip = ntet + npyr + 1; ip <= ntet + npyr + npri; ip++)  //30
    {
        for (int isg = 1; isg <= 9; isg++)  //31
        {
            is1 = nu[ind + iaux3[isg - 1][0]];
            is2 = nu[ind + iaux3[isg - 1][1]];
            for (int kv = ndeg[is1 - 1] + 1; kv <= ndeg[is1]; kv++)  // neis of 1
            {
                iseg = jaret[max_str + kv];
                if (iseg == 0 || iseg > nseg)
                    continue;
                if (iseg > maxseg) {
                    cout << "numsegs 5: increase MAXSEG" << endl;
                    exit(0);
                }
                i1 = nubo[0][iseg - 1];
                i2 = nubo[1][iseg - 1];
                is_min = min(is1, is2);
                is_max = max(is1, is2);
                if (is_min == i1 && is_max == i2)  // This Segment Exists
                {
                    isaid = iseg;
                    goto LINE32;
                }
            }
            // this is a new segment (no-name):
            nseg = nseg + 1;
            if (nseg > maxseg) {
                cout << "numsegs 6: increase MAXSEG" << endl;
                exit(0);
            }
            isaid = nseg;
            kpoin_1 = jaret[max_str + ndeg[is1]] + 1;
            jaret[max_str + ndeg[is1]] = kpoin_1;
            jaret[kpoin_1 + max_str] = nseg;
            kpoin_2 = jaret[max_str + ndeg[is2]] + 1;
            jaret[max_str + ndeg[is2]] = kpoin_2;
            jaret[kpoin_2 + max_str] = nseg;
            nubo[0][nseg - 1] = min(is1, is2);
            nubo[1][nseg - 1] = max(is1, is2);
        LINE32:
            nusg[inds + isg] = isaid;
        }
        ind += 6;
        inds += 9;
    }
    //
    //	Sweep all hexahedra
    for (int ip = ntet + npyr + npri + 1; ip <= ntet + npyr + npri + nhex; ip++)  //40
    {
        for (int isg = 1; isg <= 12; isg++)  // 41 // Sweep segments of the hexahedron
        {
            is1 = nu[ind + iaux4[isg - 1][0]];
            is2 = nu[ind + iaux4[isg - 1][1]];
            for (int kv = ndeg[is1 - 1] + 1; kv <= ndeg[is1]; kv++) {
                iseg = jaret[max_str + kv];
                if (iseg == 0 || iseg > nseg)
                    continue;
                if (iseg > maxseg) {
                    cout << "numsegs 7: increase MAXSEG";
                    exit(0);
                }
                i1 = nubo[0][iseg - 1];
                i2 = nubo[1][iseg - 1];
                is_min = min(is1, is2);
                is_max = max(is1, is2);
                if (is_min == i1 && is_max == i2)  // This Segment Exists
                {
                    isaid = iseg;  // jaret[max_str+kv]
                    goto LINE42;
                }
            }
            // This is a new segment (no-name):
            nseg = nseg + 1;
            if (nseg > maxseg) {
                cout << "numsegs 8: increase MAXSEG" << endl;
                exit(0);
            }
            isaid = nseg;
            kpoin_1 = jaret[max_str + ndeg[is1]] + 1;
            jaret[max_str + ndeg[is1]] = kpoin_1;
            jaret[kpoin_1 + max_str] = nseg;
            kpoin_2 = jaret[max_str + ndeg[is2]] + 1;
            jaret[max_str + ndeg[is2]] = kpoin_2;
            jaret[kpoin_2 + max_str] = nseg;
            nubo[0][nseg - 1] = min(is1, is2);
            nubo[1][nseg - 1] = max(is1, is2);
        LINE42:
            nusg[inds + isg] = isaid;
        }
        ind += 8;
        inds += 12;
    }

    cout << endl
         << " ---------------------------------" << endl;
    cout << " TOTAL NUMBER OF SEGMENTS:	" << nseg << endl;
    cout << " ---------------------------------" << endl;
    //
    //	Allocate the other segment-based arrays
    vnocl = matrix<double>(3, nseg);

    return;
}

void DataStructures::perseg() {
    double inod1, inod2, signs;
    int kount, ifac, is1, is2, is3, is4, ip, kk, ind2, iseg2,
        ndif, kou, k1, k2, nse, inds, iseglo, iseg, npersg_exp,
        inod3, inod4, kkk, iseg1, mflag1, n1, n2;
    int iauxs1[4][3] = {{4, 5, 6}, {2, 3, 6},  // elements type 1 --> tetrahedra
                        {1, 3, 5},
                        {1, 2, 4}};
    int iauxs2[5][4] = {{1, 3, 5, 0}, {6, 7, 8, 0}, {2, 3, 8, 0},  // elements type 2 --> pyramids
                        {4, 5, 7, 0},
                        {1, 2, 4, 6}};
    int iauxs3[5][4] = {{1, 3, 5, 7}, {1, 2, 4, 0}, {7, 8, 9, 0},  // elements type 3 --> prisms
                        {4, 5, 6, 9},
                        {2, 3, 6, 8}};
    int iauxs4[6][4] = {{4, 5, 7, 11}, {2, 3, 8, 10}, {1, 3, 5, 9},  // elements type 4 --> hexahedra
                        {9, 10, 11, 12},
                        {1, 2, 4, 6},
                        {6, 7, 8, 12}};
    vector<int> marksgper, nod, icom;
    marksgper.resize(maxseg + 1);
    nod.resize(5);
    icom.resize(ns + 1);
    //
    //	--------------------------------------------
    //	FIND THE CORRESPONDANCE OF PERIODIC SEGMENTS
    //	--------------------------------------------
    //
    //	init icom and marksgper
    for (int i = 1; i <= ns; i++) {
        icom[i] = -1;
    }
    for (int i = 1; i <= nseg; i++) {
        marksgper[i] = 0;
    }
    //
    //	fill icom
    for (int i = 1; i <= nper; i++) {
        inod1 = iper[0][i - 1];
        inod2 = iper[1][i - 1];
        icom[inod1] = inod2;
        icom[inod2] = inod1;
    }
    //
    npersg = 0;
    kount = 0;

    // FIND THE INDEX OF THE PERIODIC SEGMENTS FROM
    // THE CORRESPONDING PERIODIC FACE
    //	==============================
    for (int ii = listbf1[1]; ii <= listbf2[1]; ii++) {
        ifac = listbf[ii];
        is1 = nubf[0][ifac - 1];
        is2 = nubf[1][ifac - 1];
        is3 = nubf[2][ifac - 1];
        is4 = nubf[3][ifac - 1];
        ip = ibfc2te[ifac];

        if (ip <= ntet)  // find the kind of the adjacent elem.
        {
            kk = 4;
            ind2 = (ip - 1) * 4;  // node index
            ndif = 2;
        } else if (ip <= ntet + npyr) {
            kk = 5;
            ind2 = ntet * 4 + (ip - ntet - 1) * 5;  // node index
            ndif = 3;
            if (is4 != 0)
                ndif = 2;
        } else if (ip <= ntet + npyr + npri) {
            kk = 6;
            ind2 = ntet * 4 + npyr * 5 + (ip - ntet - npyr - 1) * 6;  // node index
            ndif = 4;
        } else {
            kk = 8;
            ind2 = ntet * 4 + npyr * 5 + npri * 6 + (ip - ntet - npyr - npri - 1) * 8;  // node index
            ndif = 5;
        }
        kou = 0;
        for (int m1 = 1; m1 <= kk; m1++)  // Loop on the kk nodes of the element
        {
            k1 = nu[ind2 + m1];
            k2 = 1;
            if (k1 == is1 || k1 == is2 || k1 == is3 || k1 == is4)
                k2 = 0;
            if (k2 != 0) {
                kou++;
                if (kou == ndif) {
                    cout << "Problem ndif in perseg" << endl;
                    exit(0);
                }
                nod[kou] = m1;
            }
        }
        if (kou < 1) {
            cout << "Perseg: Error, with data!" << endl;
            exit(0);
        }
        if (kk == 5 && nod[1] == 3)
            nod[kou] = 1;
        if (kk == 6 && nod[1] == 3)
            nod[kou] = 1;
        if (kk == 6 && nod[1] == 4)
            nod[kou] = 2;
        if (nod[kou] == 8 && nod[1] == 1)
            nod[kou] = 1;
        if (nod[kou] == 7 && nod[1] == 2)
            nod[kou] = 2;
        if (nod[kou] == 8 && nod[1] == 3)
            nod[kou] = 3;
        if (nod[kou] == 8 && nod[1] == 5)
            nod[kou] = 5;

        nse = kk + 1 - ndif;  //# periodic segments
        // REGISTER ALL PERIODIC SEGMENTS TO IBSG2TE -ONCE-
        // ------------------------------------------------
        if (ip <= ntet) {
            inds = (ip - 1) * 6;  // segment index
            for (int ise = 1; ise <= nse; ise++) {
                iseglo = iauxs1[nod[kou] - 1][ise - 1];
                iseg = nusg[inds + iseglo];
                if (marksgper[iseg] != 1)  // non-dupe
                {
                    kount++;
                    if (kount > maxsegpe) {
                        cout << "perseg tet: increase MAXSEGPE" << endl;
                        exit(0);
                    }
                    ibsg2te[kount] = iseg;
                    marksgper[iseg] = 1;
                }
            }
        } else if (ip <= ntet + npyr) {
            inds = ntet * 6 + (ip - ntet - 1) * 8;  // segment index
            for (int ise = 1; ise <= nse; ise++) {
                iseglo = iauxs2[nod[kou] - 1][ise - 1];
                iseg = nusg[inds + iseglo];
                if (marksgper[iseg] != 1)  // non-dupe
                {
                    kount++;
                    if (kount > maxsegpe) {
                        cout << "perseg pyr: increase MAXSEGPE" << endl;
                        exit(0);
                    }
                    ibsg2te[kount] = iseg;
                    marksgper[iseg] = 1;
                }
            }
        } else if (ip <= ntet + npyr + npri) {
            inds = ntet * 6 + npyr * 8 + (ip - ntet - npyr - 1) * 9;  // segment index
            for (int ise = 1; ise <= nse; ise++) {
                iseglo = iauxs3[nod[kou] - 1][ise - 1];
                iseg = nusg[inds + iseglo];
                if (marksgper[iseg] != 1)  // non-dupe
                {
                    kount++;
                    if (kount > maxsegpe) {
                        cout << "perseg pri: increase MASEGPE" << endl;
                        exit(0);
                    }
                    ibsg2te[kount] = iseg;
                    marksgper[iseg] = 1;
                }
            }
        } else  // hexahedra
        {
            inds = ntet * 6 + npyr * 8 + npri * 9 + (ip - ntet - npyr - npri - 1) * 12;  // segment index
            for (int ise = 1; ise <= nse; ise++) {
                iseglo = iauxs4[nod[kou] - 1][ise - 1];
                iseg = nusg[inds + iseglo];
                if (marksgper[iseg] != 1)  // non-dupe
                {
                    kount++;
                    if (kount > maxsegpe) {
                        cout << "perseg hex: increase MAXSEGPE" << endl;
                        exit(0);
                    }
                    ibsg2te[kount] = iseg;
                    marksgper[iseg] = 1;
                }
            }
        }
        // =====
    }  // end of loop on periodic faces
    // =====

    cout << " -----------------------------------------" << endl;
    cout << " TOTAL NUMBER OF PERIODIC SEGMENTS:	" << kount << endl;
    cout << " -----------------------------------------" << endl;

    npersg_exp = kount / 2;
    ipersg = matrix<int>(2, npersg_exp);

    //	MATCH THE PERIODIC SEGMENTS
    //	---------------------------
    for (int i = 1; i <= kount; i++)  // 112
    {
        iseg = ibsg2te[i];
        if (iseg < 0)
            continue;  // already treated
        inod1 = nubo[0][iseg - 1];
        inod2 = nubo[1][iseg - 1];
        inod1 = icom[inod1];
        inod2 = icom[inod2];
        for (int j = i + 1; j <= kount; j++)  // 111
        {
            iseg2 = ibsg2te[j];
            if (iseg2 < 0)
                continue;
            inod3 = nubo[0][iseg2 - 1];
            inod4 = nubo[1][iseg2 - 1];
            if (inod3 == inod1 || inod3 == inod2) {
                if (inod4 == inod1 || inod4 == inod2)  // corresponding seg.
                {
                    ibsg2te[j] = -iseg2;
                    npersg++;

                    if (inod3 == inod1) {
                        ipersg[0][npersg - 1] = iseg;
                        ipersg[1][npersg - 1] = iseg2;
                    } else {
                        ipersg[0][npersg - 1] = -iseg;
                        ipersg[1][npersg - 1] = -iseg2;
                    }
                    break;
                }
            }
        }
    }
    //
    // CHECK IF NODES OF IPERSG(1,...) CORRESPOND TO IPER(1,...)
    // IF NOT SWAP IPERSG(1,...) -><- IPERSG(2,...)
    // -----------------------------
    for (int ipsg = 1; ipsg <= npersg; ipsg++)  //200
    {
        kkk = ipersg[0][ipsg - 1];
        signs = copysign(1., kkk);
        iseg1 = abs(ipersg[0][ipsg - 1]);
        iseg2 = abs(ipersg[1][ipsg - 1]);
        mflag1 = 0;
        n1 = nubo[0][iseg1 - 1];
        for (int ii = 1; ii <= nper; ii++) {
            k1 = iper[0][ii - 1];
            if (n1 == k1) {
                mflag1++;
                break;
            }
        }
        n2 = nubo[1][iseg1 - 1];
        for (int ii = 1; ii <= nper; ii++) {
            k2 = iper[0][ii - 1];
            if (n2 == k2) {
                mflag1++;
                break;
            }
        }
        if (mflag1 == 0) {
            ipersg[0][ipsg - 1] = iseg2 * signs;
            ipersg[1][ipsg - 1] = iseg1 * signs;
        } else if (mflag1 == 1) {
            cout << "STOP in perseg: mixed iper i,j" << endl;
            exit(0);
        }
    }
    //
    icom.clear();
    return;
}

void DataStructures::fjaret3D() {
    int maxindex = 10000, nvmaxall_out, move, index, ielem, iaux, ind,
        knod, is1, isn, nei;
    vector<int> newnod;
    newnod.resize(maxindex + 1);
    // iaux[i][j] : j=nodes connected to i with edge in element
    int iaux1[4][3] = {{2, 3, 4}, {1, 3, 4},  // tetrahedra
                       {1, 2, 4},
                       {1, 2, 3}};
    int iaux2[5][4] = {{2, 4, 5, 0}, {1, 3, 5, 0}, {2, 4, 5, 0},  // pyramids
                       {1, 3, 5, 0},
                       {1, 2, 3, 4}};
    int iaux3[6][3] = {{2, 3, 4}, {1, 3, 5}, {1, 2, 6},  // prisms
                       {1, 5, 6},
                       {2, 4, 6},
                       {3, 4, 5}};
    int iaux4[8][3] = {{2, 4, 5}, {1, 3, 6}, {2, 4, 7},  // hexahedra
                       {1, 3, 8},
                       {1, 6, 8},
                       {2, 5, 7},
                       {3, 6, 8},
                       {4, 5, 7}};
    //
    fjaret_el();
    //
    nvmaxall_out = ndeg[ns];
    //	Move ndeg & Jaret
    move = ns + 1;  // nodes plus 0-position
    //
    if (2 * ndeg[ns] > nvmaxall) {
        cout << " fjaret 3: Increase nvmaxall ---- > " << endl;
        cout << " Current nvmaxall: " << nvmaxall << endl;
        exit(0);
    }
    //
    nvmaxall_out = max(nvmaxall_out, 2 * ndeg[ns]);
    cout << "nvmaxall_out =	" << nvmaxall_out << endl;
    //
    // Move jaret with elements to new positions
    for (int is = 0; is <= ns; is++) {
        ndeg[is + move] = ndeg[is] + nvmaxall / 2;  // split jaret equally
    }
    for (int is = 1; is <= ns; is++) {
        for (int ii = ndeg[is - 1] + 1; ii <= ndeg[is]; ii++) {
            jaret[ii + nvmaxall / 2] = jaret[ii];
        }
    }
    //
    for (int is = 1; is <= ns; is++)  // 40 // SWEEP NODES
    {
        index = 0;
        // In newnod Store the NEW nei Nodes
        for (int ii = 1; ii <= maxindex; ii++) {
            newnod[ii] = 0;
        }
        for (int i = ndeg[is - 1 + move] + 1; i <= ndeg[is + move]; i++)  // 41
        {
            ielem = jaret[i];  // SWEEP ELEMENTS AROUND NODE
            if (ielem <= ntet) {
                knod = 4;  // it's a tetrahedron
                ind = (ielem - 1) * 4;
                iaux = 3;
            } else if (ielem <= ntet + npyr) {
                knod = 5;  // it's a pyramid
                ind = ntet * 4 + (ielem - ntet - 1) * 5;
                iaux = 4;
            } else if (ielem <= ntet + npyr + npri) {
                knod = 6;  // it's a prism
                ind = ntet * 4 + npyr * 5 + (ielem - ntet - npyr - 1) * 6;
                iaux = 3;
            } else {
                knod = 8;  // it's a hexahedron
                ind = ntet * 4 + npyr * 5 + npri * 6 + (ielem - ntet - npyr - npri - 1) * 8;
                iaux = 3;
            }

            for (int j = 1; j <= knod; j++)  // 42 // SWEEP NODES OF ELEMENT AROUND NODE
            {
                is1 = nu[ind + j];
                if (is1 == is) {
                    for (int kk = 1; kk <= iaux; kk++)  // 43 //
                    {
                    LOOP43:
                        if (kk > iaux) {
                            break;
                        }
                        if (knod == 4)  // tetrahedron
                        {
                            isn = iaux1[j - 1][kk - 1];
                        } else if (knod == 5)  // pyramid
                        {
                            isn = iaux2[j - 1][kk - 1];
                            if (isn == 0) {
                                break;
                            }
                        } else if (knod == 6)  // prism
                        {
                            isn = iaux3[j - 1][kk - 1];
                        } else  // hexahedron
                        {
                            isn = iaux4[j - 1][kk - 1];
                        }
                        for (int ij = 1; ij <= index; ij++) {
                            nei = newnod[ij];
                            if (nei == nu[ind + isn]) {
                                kk++;
                                goto LOOP43;
                            }
                        }
                        index++;
                        if (index > maxindex) {
                            cout << "fjaret:increase dimension of newnod" << endl;
                            exit(0);
                        }
                        newnod[index] = nu[ind + isn];
                        ndeg[is] = ndeg[is - 1] + index;  // counter of jaret
                        jaret[ndeg[is]] = nu[ind + isn];
                    }
                }
            }
        }
    }
    //
    // Zero neww positions
    for (int i = ns + 1; i <= 2 * ns + 1; i++) {
        ndeg[i] = 0;
    }

    max_str = ndeg[ns];

    // Zero remaining jaret
    int i;
    if (2 * ndeg[ns] > nvmaxall) {
        i = ndeg[ns] + max_str;
        cout << "fjaret: Increase nvmaxall --- >	" << i << endl;
        cout << "Current nvmaxall:	" << nvmaxall;
    }
    for (int i = ndeg[ns] + 1; i <= 2 * ndeg[ns]; i++) {
        jaret[i] = 0;
    }

    return;
}

void DataStructures::fjaret_el() {
    int is, kpoi_1;

    for (unsigned int i = 0; i <= ns; i++)
        ndeg[i] = 0;
    for (unsigned int i = 1; i <= nvmaxall; i++)
        jaret[i] = 0;

    int ind = 0;
    // Sweep tetrahedra
    for (unsigned int i = 1; i <= ntet; i++) {
        for (unsigned int j = 1; j <= 4; j++) {
            ind++;
            is = nu[ind];
            ndeg[is] = ndeg[is] + 1;
        }
    }
    // Sweep Pyramids
    for (unsigned int i = ntet + 1; i <= ntet + npyr; i++) {
        for (unsigned int j = 1; j <= 5; j++) {
            ind++;
            is = nu[ind];
            ndeg[is] = ndeg[is] + 1;
        }
    }
    // Sweep Prisms
    for (unsigned int i = ntet + npyr + 1; i <= ntet + npyr + npri; i++) {
        for (unsigned int j = 1; j <= 6; j++) {
            ind++;
            is = nu[ind];
            ndeg[is] = ndeg[is] + 1;
        }
    }
    // Sweep hexahedra
    for (unsigned int i = ntet + npyr + npri + 1; i <= ntet + npyr + npri + nhex; i++) {
        for (unsigned int j = 1; j <= 8; j++) {
            ind++;
            is = nu[ind];
            ndeg[is] = ndeg[is] + 1;
        }
    }
    // Precaution
    int istr = 0;
    for (unsigned int k = 1; k <= ns; k++) {
        istr += ndeg[k];
    }
    if (istr > nvmaxall) {
        cout << " fjaret_el:Increase nvmaxall--> " << istr << endl;
        cout << " Current nvmaxall: " << nvmaxall << endl;
        exit;
    }

    for (unsigned int k = 1; k <= ns; k++) {
        if (ndeg[k] == 0) {
            cout << endl;
            cout << " ERROR: hanging node, index: " << k << endl;
            exit;
        }
        ndeg[k] += ndeg[k - 1];
        jaret[ndeg[k]] = ndeg[k - 1];  // provisory
    }

    ind = 0;
    // Sweep tetrahedra
    for (unsigned int i = 1; i <= ntet; i++) {
        for (unsigned int j = 1; j <= 4; j++) {
            ind++;
            is = nu[ind];
            jaret[ndeg[is]] = jaret[ndeg[is]] + 1;
            kpoi_1 = jaret[ndeg[is]];
            jaret[kpoi_1] = i;  // Tetrahedra
        }
    }

    // Sweep Pyramids
    for (unsigned int i = ntet + 1; i <= ntet + npyr; i++) {
        for (unsigned int j = 1; j <= 5; j++) {
            ind++;
            is = nu[ind];
            jaret[ndeg[is]] = jaret[ndeg[is]] + 1;
            kpoi_1 = jaret[ndeg[is]];
            jaret[kpoi_1] = i;  // Pyramids
        }
    }

    // Sweep prisms
    for (unsigned int i = ntet + npyr + 1; i <= ntet + npyr + npri; i++) {
        for (unsigned int j = 1; j <= 6; j++) {
            ind++;
            is = nu[ind];
            jaret[ndeg[is]] = jaret[ndeg[is]] + 1;
            kpoi_1 = jaret[ndeg[is]];
            jaret[kpoi_1] = i;  // Prisms
        }
    }
    // Sweep hexahedra
    for (unsigned int i = ntet + npyr + npri + 1; i <= ntet + npyr + npri + nhex; i++) {
        for (unsigned int j = 1; j <= 8; j++) {
            ind++;
            is = nu[ind];
            jaret[ndeg[is]] = jaret[ndeg[is]] + 1;
            kpoi_1 = jaret[ndeg[is]];
            jaret[kpoi_1] = i;  // Hexahedra
        }
    }
    return;
}

void DataStructures::matrix_mult(int nvar, double **xx, double **yy, double **zz,
                                 int mvar1, int mvar2, int mvar3) {
    for (unsigned int i = 0; i < mvar1; i++) {
        for (unsigned int j = 0; j < mvar3; j++) {
            zz[i][j] = 0.;
            for (unsigned int k = 0; k < mvar2; k++) {
                zz[i][j] = zz[i][j] + zz[i][k] * yy[k][j];
            }
        }
    }
    return;
}

void DataStructures::rotation_matrix(int nvar, int iop, double cm, double sm, double **rra, double **rrb) {
    //rra = matrix_double(nvar, nvar);
    rra = matrix<double>(nvar, nvar);
    rrb = matrix<double>(nvar, nvar);
    //rrb = matrix_double(nvar, nvar);
    for (unsigned int i = 0; i < nvar; i++) {
        for (unsigned int j = 0; j < nvar; j++) {
            rra[i][j] = 0.;
            rrb[i][j] = 0.;
        }
        rra[i][i] = 1.;
        rrb[i][i] = 1.;
    }

    if (iop > nvar - 1) {
        cout << "rotation_matrix called with invalid iop" << endl;
        exit;
    }

    rra[iop][iop] = cm;
    rra[iop][iop + 1] = -sm;
    rra[iop + 1][iop] = sm;
    rra[iop + 1][iop + 1] = cm;
    rrb[iop][iop] = cm;
    rrb[iop][iop + 1] = -sm;
    rrb[iop + 1][iop] = sm;
    rrb[iop + 1][iop + 1] = cm;
    return;
}

void DataStructures::virtualPeriodicNeighbours() {
    vector<int> neisUP, neisDown;
    vector<int> newInsert;
    int ifriend, inode, ivpseg1, ivpseg2, ifriendnei, inodenei;
    for (int ip = 0; ip < nper; ip++) {
        inode = iper[1][ip];
        ifriend = iper[0][ip];
        int counter = 0;
        for (unsigned int k = ndeg[ifriend - 1] + 1; k <= ndeg[ifriend]; k++) {
            ifriendnei = jaret[k];

            if (logfr[ifriendnei] == 0 && logfr[ifriend] == 1) {
                counter++;
                newInsert.push_back(ifriendnei);
                neisDown.push_back(ifriendnei);
            }
        }

        if (counter == 0) {
            continue;
        }

        // keep old jaret from inode + 1 and after

        vector<int>::iterator position;
        position = jaret.begin();
        advance(position, ndeg[inode] + 1);
        insert_iterator<vector<int>> jaret_inserter(jaret, position);

        // insert new virtual neighbours of inode in jaret after inode
        int kounter = 0;
        for (int neighbour : newInsert) {
            *jaret_inserter = neighbour;
        }

        // ndeg_new  = ndeg_old + number of new neighbours
        for (int id = inode; id < ndeg.size(); id++) {
            ndeg[id] += counter;
        }

        newInsert.clear();

        if (logfr[ifriend] == 1)
            logfr[ifriend] = 11;
    }

    for (int x : neisDown) {
        logfr[x] = 110;
    }

    neisDown.clear();
}

template <typename myType>
myType **matrix(int rows, int columns) {
    myType **matr = new myType *[rows];
    for (unsigned int i = 0; i < rows; i++) {
        matr[i] = new myType[columns];
    }
    return matr;
}

DataStructures::~DataStructures() {
}