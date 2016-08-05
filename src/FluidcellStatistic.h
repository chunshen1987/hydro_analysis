#ifndef SRC_FLUIDCELLSTATISTIC_H_
#define SRC_FLUIDCELLSTATISTIC_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include <fstream>

#include "./Hydroinfo_h5.h"
#include "./ParameterReader.h"

using namespace std;

class FluidcellStatistic {
 private:
    int nbin;
    HydroinfoH5 *hydroinfo_ptr;
    ParameterReader *paraRdr;
    double t_i, t_f, df;
    double *tempature_ptr;
    double *fluidcellvolume_ptr;
    double *flowVeclocity_ptr;
    double *energyDensity_ptr;

public:
    FluidcellStatistic(HydroinfoH5* hydroinfo_ptr_in,
                       ParameterReader* paraRdr_in);
    ~FluidcellStatistic();
    void checkFreezeoutSurface(double Tdec);
    double getTemperature(double tau_local, double x_local, double y_local);
    void checkMomentumAnisotropy(double tau_local);
    void calculateAvgandStdtemperature(double tau_local);
    void calculateAvgandStdflowvelocity(double tau_local);
    void Countcellvolume(double temp_local, double volume);
    void OutputCellvolume(string filename);
    void outputTempasTauvsX();
    void calAvgVvsT(double temp_local, double volume, double velocity,
                    double ed);
    void outputAvgV(string filename);
    void outputKnudersonNumberasTauvsX();
    void outputinverseReynoldsNumberasTauvsX();
    void outputBulkinverseReynoldsNumberasTauvsX();
    double calculate_spacetime_4volume(double T_cut);
    double calculate_average_tau(double T_cut);
    double calculate_hypersurface_3volume(double T_cut);
};

#endif  // SRC_FLUIDCELLSTATISTIC_H_
