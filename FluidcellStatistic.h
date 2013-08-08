#ifndef FLUIDCELLSTATISTIC_H
#define FLUIDCELLSTATISTIC_H

#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>
#include<fstream>

using namespace std;

class FluidcellStatistic
{
    private:
        int nbin;
        double t_i, t_f, df;
        double *tempature_ptr;
        double *fluidcellvolume_ptr;
        double *flowVeclocity_ptr;
        double *energyDensity_ptr;


    public:
        FluidcellStatistic();
        ~FluidcellStatistic();

        void Countcellvolume(double temp_local, double volume);
        void OutputCellvolume(string filename);
        void calAvgV(double temp_local, double volume, double velocity, double ed);
        void outputAvgV(string filename);

};

#endif
