#ifndef FLUIDCELLSTATISTIC_H
#define FLUIDCELLSTATISTIC_H

#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>
#include<fstream>

#include "Hydroinfo_h5.h"

using namespace std;

class FluidcellStatistic
{
    private:
        int nbin;
        HydroinfoH5 *hydroinfo_ptr;
        double t_i, t_f, df;
        double *tempature_ptr;
        double *fluidcellvolume_ptr;
        double *flowVeclocity_ptr;
        double *energyDensity_ptr;


    public:
        FluidcellStatistic(HydroinfoH5* hydroinfo_ptr_in);
        ~FluidcellStatistic();
        
        void checkFreezeoutSurface(double Tdec);
        double getTemperature(double tau_local, double x_local, double y_local);
        void checkMomentumAnisotropy(double tau_local);
        void calculateAvgandStdtemperature(double tau_local);
        void calculateAvgandStdflowvelocity(double tau_local);
        void Countcellvolume(double temp_local, double volume);
        void OutputCellvolume(string filename);
        void outputTempasTauvsX();
        void calAvgVvsT(double temp_local, double volume, double velocity, double ed);
        void outputAvgV(string filename);

};

#endif
