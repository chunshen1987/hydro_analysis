/////////////////////////////////////////////////////////////////////////
//                      hydrodynamics analysis
//                          photon emission 
//
//              author: Chun Shen <shen@mps.ohio-state.edu>
//
//  This program calculates the photon emission from the relativistic
//  heavy ion collision. It reads in viscous hydrodynamics results in 
//  OSCAR format, which is achieved by the functions in readframe.cpp.
//  
//
//  To do in the future:
/////////////////////////////////////////////////////////////////////////

#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>

#include "Hydroinfo_h5.h"
#include "Stopwatch.h"
#include "FluidcellStatistic.h"

using namespace std;


int main()
{
  Stopwatch sw;

  sw.tic();
  HydroinfoH5* hydroinfo_ptr = new HydroinfoH5("results/JetData.h5", 500, 0);   //hydro data file pointer

  double T_dec = 0.12;
  //get grid information from OSCAR header
  double grid_t0 = 0.6;
  double grid_x0, grid_y0;
  grid_x0 = -13.0;
  grid_y0 = -13.0;
  double grid_dt, grid_dx, grid_dy;
  grid_dt = 0.04;
  grid_dx = 0.2;
  grid_dy = 0.2;
  int ntime = (int)((hydroinfo_ptr->getHydrogridTaumax() - grid_t0)/grid_dt) + 1;
  int nx = (int)(abs(2*grid_x0)/grid_dx) + 1;
  int ny = (int)(abs(2*grid_y0)/grid_dy) + 1;
  double volume;  //volume of the cell
  volume = grid_dt*grid_dx*grid_dy;
  
  FluidcellStatistic Fluidcellanalysis;

  fluidCell* fluidCellptr = new fluidCell();
  double e_local, temp_local, vx_local, vy_local, vz_local, tau_local, x_local, y_local;

  for(int itime=0;itime<ntime;itime++) //loop over time evolution
  {
    tau_local = grid_t0 + itime*grid_dt;
    for(int i=0;i<nx;i++) //loops over the transverse plane
    {
      x_local = grid_x0 + i*grid_dx;
      for(int j=0;j<ny;j++)
      {
        y_local = grid_y0 + j*grid_dy;
        hydroinfo_ptr->getHydroinfo(tau_local, x_local, y_local, fluidCellptr);
        temp_local = fluidCellptr->temperature;
        if(temp_local > T_dec)
        {
           e_local = fluidCellptr->ed;
           vx_local = fluidCellptr->vx;
           vy_local = fluidCellptr->vy;
           vz_local = 0.0;
           double v_perp = sqrt(vx_local*vx_local + vy_local*vy_local);
        
           Fluidcellanalysis.Countcellvolume(temp_local, tau_local*volume);
           Fluidcellanalysis.calAvgV(temp_local, tau_local*volume, v_perp, e_local);
        }
      }
    }
    cout<<"frame "<< itime << " : ";
    cout<<" tau = " << setw(4) << setprecision(3) << tau_local;
    cout<<"done!" <<endl ;
  }

  Fluidcellanalysis.OutputCellvolume("TvsCellvolume");
  Fluidcellanalysis.outputAvgV("TvsAvgV");

  sw.toc();
  cout << "totally takes : " << sw.takeTime() << " seconds." << endl;

  return(0);
}


