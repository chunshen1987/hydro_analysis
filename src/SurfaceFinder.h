#ifndef SRC_SurfaceFinder_H_
#define SRC_SurfaceFinder_H_

#include "./Hydroinfo_h5.h"
#include "./ParameterReader.h"

using namespace std;

class SurfaceFinder {
 private:
    HydroinfoH5 *hydroinfo_ptr;
    ParameterReader *paraRdr;
    double T_cut;

 public:
    SurfaceFinder(HydroinfoH5* hydroinfo_ptr_in, ParameterReader* paraRdr_in);
    SurfaceFinder(HydroinfoH5* hydroinfo_ptr_in, ParameterReader* paraRdr_in,
                  double T_cut_in);
    ~SurfaceFinder();

    bool check_intersect(double T_cut, double tau, double x, double y,
                         double dt, double dx, double dy, double ***cube);
    int Find_full_hypersurface();
};

#endif  // SRC_SurfaceFinder_H_
