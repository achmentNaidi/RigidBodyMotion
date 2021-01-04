#include <cmath>
#include <numeric>
#include <vector>
#

using namespace std;

struct point {
   public:
    double x;
    double y;
    double z = 0;
    point operator-(const point&) const;
    point operator+(const point&) const;
    point operator/(int) const;
};

point point::operator-(const point& p1) const {
    point tem;
    tem = *this;
    tem.x -= p1.x;
    tem.y -= p1.y;
    tem.z -= p1.z;
    return tem;
}

point point::operator/(int divisor) const {
    point tem;
    tem = *this;
    tem.x /= divisor;
    tem.y /= divisor;
    tem.z /= divisor;
    return tem;
}

point point::operator+(const point& p1) const {
    point tem;
    tem = *this;
    tem.x += p1.x;
    tem.y += p1.y;
    tem.z += p1.z;
    return tem;
}

class deformAlgebraic {
   private:
    double** coorp;
    double** coor;
    int ns;
    double Ldef, a, b, alpha;
    vector<int> logfr;
    double wi(const point& ri, const point& r);
    double coumputeLdef(vector<point>&);

   public:
    deformAlgebraic(double** initial, double** deformed, vector<int>& flags, int nodes)
        : coor(initial),
          coorp(deformed),
          logfr(flags),
          ns(nodes) {}
    double** IDW(double a, double b);
    double twonorm(point);
    vector<point> getBoundary(double** coordinates);
    double computeAlpha(vector<point>& si);

    ~deformAlgebraic() {}
};

inline double deformAlgebraic::wi(const point& ri, const point& r) {
    double p = 2.;
    return (pow(Ldef / twonorm(r - ri), a) +
            pow(alpha * Ldef / twonorm(r - ri), b));
    // return pow(1.0 / twonorm(r - ri), p);
}

double deformAlgebraic::coumputeLdef(vector<point>& bPoints) {
    // bounding box of deformed geometry
    double xmin = 100000, ymin = 100000, xmax = -10000, ymax = -10000;
    for (point p : bPoints) {
        xmin = min(xmin, p.x);
        ymin = min(ymin, p.y);
        xmax = max(xmax, p.x);
        ymax = max(ymax, p.y);
    }

    point p2{xmin - xmax, ymax - ymin};
    return twonorm(p2);
}

double** deformAlgebraic::IDW(double a, double b) {
    this->a = a;
    this->b = b;
    // get all boundary nodes coordinates (deformed,initial)
    vector<point> rb_def = getBoundary(coorp);
    vector<point> rb_ini = getBoundary(coor);
    vector<point> si(rb_def.size());  // deformation of boundary node

    for (int i = 0; i < si.size(); i++) {
        si[i] = rb_def[i] - rb_ini[i];
    }

    // 1. COMPUTE LDEF
    Ldef = coumputeLdef(rb_def);
    // alpha = computeAlpha(si);
    alpha = 0.1;

    point r;  // Internal node coordinates
    double w;

    for (int i = 0; i < ns; i++) {
        // for every internal node
        if (logfr[i + 1] == 0) {
            double sumx = 0, sumy = 0, sumz = 0;
            double sumw = 0;
            r = {coor[0][i], coor[1][i], coor[2][i]};
            //
            for (int j = 0; j < rb_def.size(); j++) {
                // compute deformations
                // compute weights for internal node i
                w = wi(rb_def[j], r);
                // compute sums (x,y,z)
                sumx += w * si[j].x;
                sumy += w * si[j].y;
                sumz += w * si[j].z;
                sumw += w;
            }
            // IDW formula for every internal node
            coorp[0][i] += sumx / sumw;
            coorp[1][i] += sumy / sumw;
            coorp[2][i] += sumz / sumw;
        }
    }
    return coorp;
}

inline double deformAlgebraic::twonorm(point p) {
    return sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
}

vector<point> deformAlgebraic::getBoundary(double** coordinates) {
    vector<point> ri;
    for (int i = 0; i < ns; i++) {
        if (logfr[i + 1] != 0) {
            ri.push_back({coordinates[0][i], coordinates[1][i], coordinates[2][i]});
        }
    }
    return ri;
}

double deformAlgebraic::computeAlpha(vector<point>& si) {
    point s_mean;
    double eta = 5;
    for (point p : si) {
        s_mean = s_mean + p;
    }
    s_mean = s_mean / si.size();

    double a_tem = -1000000;
    for (int i = 0; i < si.size(); i++) {
        a_tem = max(a_tem, twonorm(si[i] - s_mean));
    }
    return (eta / Ldef) * a_tem;
}
