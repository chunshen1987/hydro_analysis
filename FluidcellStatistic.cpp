#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>
#include<fstream>

#include "FluidcellStatistic.h"

using namespace std;

FluidcellStatistic::FluidcellStatistic()
{
    nbin = 400;
    t_i = 0.1e0;
    t_f = 0.5e0;
    df = (t_f - t_i)/(nbin);

    tempature_ptr = new double [nbin];
    fluidcellvolume_ptr = new double [nbin];
    flowVeclocity_ptr = new double [nbin];
    energyDensity_ptr = new double [nbin];

    for(int i=0; i<nbin; i++)
    {
        tempature_ptr[i] = t_i + i * df + df / 2;
        fluidcellvolume_ptr[i] = 0.0e0;
        flowVeclocity_ptr[i] = 0.0e0;
        energyDensity_ptr[i] = 0.0e0;
    }
}


FluidcellStatistic::~FluidcellStatistic()
{
    delete [] tempature_ptr;
    delete [] fluidcellvolume_ptr;
    delete [] energyDensity_ptr;
    delete [] flowVeclocity_ptr;
}

void FluidcellStatistic::Countcellvolume(double temp_local, double volume)
{
   int temp_idx;
   temp_idx = (int) ((temp_local - t_i)/df);
   if(temp_idx > nbin) 
   {
       cout << "FluidcellStatistic::Countcellvolume Warning: temperature is exceed the upper boundary!" << endl;
       cout << "Temperature = " << temp_local << " GeV. Temp_max = " << t_f 
             << " GeV. " << endl;
   }
   if(temp_idx < 0) 
   {
       cout << "FluidcellStatistic::Countcellvolume Warning: temperature is exceed the lower boundary!" << endl;
       cout << "Temperature = " << temp_local << " GeV. Temp_min = " << t_i 
             << " GeV. " << endl;
   }
   fluidcellvolume_ptr[temp_idx] += volume; 

   return;
}

void FluidcellStatistic::calAvgV(double temp_local, double volume, double velocity, double ed)
{
   double gamma = 1./sqrt(1. - velocity*velocity);
   
   int temp_idx = (int) ((temp_local - t_i)/df);
   if(temp_idx > nbin) 
   {
       cout << "FluidcellStatistic::Countcellvolume Warning: temperature is exceed the upper boundary!" << endl;
       cout << "Temperature = " << temp_local << " GeV. Temp_max = " << t_f 
             << " GeV. " << endl;
   }
   if(temp_idx < 0) 
   {
       cout << "FluidcellStatistic::Countcellvolume Warning: temperature is exceed the lower boundary!" << endl;
       cout << "Temperature = " << temp_local << " GeV. Temp_min = " << t_i 
             << " GeV. " << endl;
   }

   energyDensity_ptr[temp_idx] += volume*gamma*ed; 
   flowVeclocity_ptr[temp_idx] += volume*gamma*ed*velocity; 

   return;

}

void FluidcellStatistic::OutputCellvolume(string filename)
{
    ostringstream filename_stream;
    filename_stream << filename << ".dat";
    ofstream OpCellvolume(filename_stream.str().c_str());

    for(int i=0; i<nbin ; i++)
       OpCellvolume << scientific << setw(15) << setprecision(8)
                        << tempature_ptr[i] << "   " 
                        << fluidcellvolume_ptr[i] << endl;

    OpCellvolume.close();
    return;
}

void FluidcellStatistic::outputAvgV(string filename)
{
    double eps = 1e-10;
    ostringstream filename_stream;
    filename_stream << filename << ".dat";
    ofstream OpAvgV(filename_stream.str().c_str());

    for(int i=0; i<nbin ; i++)
    {
       double avgV;
       if(energyDensity_ptr[i] < eps)
          avgV = 0.0;
       else
          avgV = flowVeclocity_ptr[i]/energyDensity_ptr[i];
       OpAvgV << scientific << setw(15) << setprecision(8)
                        << tempature_ptr[i] << "   " 
                        << avgV << endl;
    }

    OpAvgV.close();
    return;
}
