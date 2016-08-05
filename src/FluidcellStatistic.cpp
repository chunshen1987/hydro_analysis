// Copyright @ Chun Shen 2014
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include <cstdlib>

#include "./FluidcellStatistic.h"
#include "./SurfaceFinder.h"

using namespace std;

FluidcellStatistic::FluidcellStatistic(HydroinfoH5* hydroinfo_ptr_in,
                                       ParameterReader* paraRdr_in) {
    hydroinfo_ptr = hydroinfo_ptr_in;
    paraRdr = paraRdr_in;

    nbin = 400;
    t_i = 0.1e0;
    t_f = 0.5e0;
    df = (t_f - t_i)/(nbin);

    tempature_ptr = new double [nbin];
    fluidcellvolume_ptr = new double [nbin];
    flowVeclocity_ptr = new double [nbin];
    energyDensity_ptr = new double [nbin];

    for (int i=0; i<nbin; i++) {
        tempature_ptr[i] = t_i + i * df + df / 2;
        fluidcellvolume_ptr[i] = 0.0e0;
        flowVeclocity_ptr[i] = 0.0e0;
        energyDensity_ptr[i] = 0.0e0;
    }
}


FluidcellStatistic::~FluidcellStatistic() {
    delete [] tempature_ptr;
    delete [] fluidcellvolume_ptr;
    delete [] energyDensity_ptr;
    delete [] flowVeclocity_ptr;
}

double FluidcellStatistic::getTemperature(double tau_local, double x_local,
                                          double y_local) {
    fluidCell* fluidCellptr = new fluidCell();
    hydroinfo_ptr->getHydroinfo(tau_local, x_local, y_local, fluidCellptr);
    double temp_local = fluidCellptr->temperature;
    delete fluidCellptr;
    return(temp_local);
}

void FluidcellStatistic::checkFreezeoutSurface(double Tdec) {
    double grid_t0, grid_x0, grid_y0;
    grid_t0 = 0.6;
    grid_x0 = -13.0;
    grid_y0 = -13.0;
    double grid_dt, grid_dx, grid_dy;
    grid_dt = 0.04;
    grid_dx = 0.2;
    grid_dy = 0.2;
    int ntime = (int)((12.0 - grid_t0)/grid_dt) + 1;
    int nx = (int)(abs(2*grid_x0)/grid_dx) + 1;
    int ny = (int)(abs(2*grid_y0)/grid_dy) + 1;

    double tau_local;
    fluidCell* fluidCellptr = new fluidCell();
    ofstream output;
    output.open("results/checkFreezeoutSurface.dat", std::ofstream::app);

    for (int itime = 0; itime < ntime; itime++) {
        // loop over time evolution
        tau_local = grid_t0 + itime*grid_dt;
        for (int i = 0; i < nx; i++) {
            // loops over the transverse plane
            double x_local = grid_x0 + i*grid_dx;
            for (int j = 0; j < ny; j++) {
                double y_local = grid_y0 + j*grid_dy;
                hydroinfo_ptr->getHydroinfo(tau_local, x_local, y_local,
                                            fluidCellptr);
                double temp_local = fluidCellptr->temperature;
                if (fabs(temp_local - Tdec) < 0.001) {
                    output << tau_local << "   " << x_local << "   "
                           << y_local << "   "
                           << sqrt(x_local*x_local + y_local*y_local) << endl;
                }
            }
        }
    }
    output.close();
    delete fluidCellptr;
    return;
}

void FluidcellStatistic::calculateAvgandStdtemperature(double tau_local) {
    double avgTemp, stdTemp;
    double T_dec = 0.120;
    double grid_x0, grid_y0;
    grid_x0 = -13.0;
    grid_y0 = -13.0;
    double grid_dx, grid_dy;
    grid_dx = 0.2;
    grid_dy = 0.2;
    int nx = (int)(abs(2*grid_x0)/grid_dx) + 1;
    int ny = (int)(abs(2*grid_y0)/grid_dy) + 1;

    int numCell = 0;
    double tempSum = 0.0;
    double tempSumsq = 0.0;
    double tempMax = 0.0;
    fluidCell* fluidCellptr = new fluidCell();
    for (int i = 0; i < nx; i++) {
        // loops over the transverse plane
        double x_local = grid_x0 + i*grid_dx;
        for (int j = 0; j < ny; j++) {
            double y_local = grid_y0 + j*grid_dy;
            hydroinfo_ptr->getHydroinfo(tau_local, x_local, y_local,
                                        fluidCellptr);
            double temp_local = fluidCellptr->temperature;
            if (temp_local > T_dec) {
                numCell += 1;
                tempSum += temp_local;
                tempSumsq += temp_local*temp_local;
                if (tempMax < temp_local)
                    tempMax = temp_local;
            }
        }
    }
    if (numCell == 0) {
        avgTemp = 0.0;
        stdTemp = 0.0;
    } else {
        avgTemp = tempSum/numCell;
        stdTemp = sqrt((tempSumsq - tempSum*tempSum/numCell)/(numCell - 1));
    }

    ofstream output;
    output.open("results/tau_vs_T.dat", std::ofstream::app);
    output << tau_local << "   " << avgTemp << "   " << stdTemp 
           << "   " << tempMax << endl;
    delete fluidCellptr;
    return;
}

void FluidcellStatistic::calculateAvgandStdflowvelocity(double tau_local) {
    double avgV, stdV, avgVsurf;
    double T_dec = 0.120;
    double grid_x0, grid_y0;
    grid_x0 = -13.0;
    grid_y0 = -13.0;
    double grid_dx, grid_dy;
    grid_dx = 0.2;
    grid_dy = 0.2;
    int nx = (int)(abs(2*grid_x0)/grid_dx) + 1;
    int ny = (int)(abs(2*grid_y0)/grid_dy) + 1;

    double edSum = 0.0;
    double edSurfsum = 0.0;
    double vSurfsum = 0.0;
    double VSum = 0.0;
    double VSumsq = 0.0;
    fluidCell* fluidCellptr = new fluidCell();
    for (int i = 0; i < nx; i++) {
        // loops over the transverse plane
        double x_local = grid_x0 + i*grid_dx;
        for (int j = 0; j < ny; j++) {
            double y_local = grid_y0 + j*grid_dy;
            hydroinfo_ptr->getHydroinfo(tau_local, x_local, y_local,
                                        fluidCellptr);
            double temp_local = fluidCellptr->temperature;
            if (temp_local > T_dec) {
                double e_local = fluidCellptr->ed;
                double vx_local = fluidCellptr->vx;
                double vy_local = fluidCellptr->vy;
                double v_perp = sqrt(vx_local*vx_local + vy_local*vy_local);
                double gamma = 1./sqrt(1. - v_perp*v_perp);
                edSum += gamma*e_local;
                VSum += gamma*e_local*v_perp;
                VSumsq += gamma*e_local*v_perp*v_perp;
                if (fabs(temp_local - T_dec) < 0.001) {
                    vSurfsum += gamma*e_local*v_perp;
                    edSurfsum += gamma*e_local;
                }
            }
        }
    }
    if (fabs(edSum) < 1e-10) {
        avgV = 0.0;
        stdV = 0.0;
        vSurfsum += 0.0;
    } else {
        avgV = VSum/edSum;
        stdV = sqrt(VSumsq/edSum - avgV*avgV);
        avgVsurf = vSurfsum/edSurfsum;
    }

    ofstream output;
    output.open("results/tau_vs_v.dat", std::ofstream::app);
    output << tau_local << "   " << avgV << "   " << stdV 
           << "   " << avgVsurf << endl;
    return;
}

void FluidcellStatistic::checkMomentumAnisotropy(double tau_local) {
    double T_dec = 0.120;
    double grid_x0, grid_y0;
    grid_x0 = -13.0;
    grid_y0 = -13.0;
    double grid_dx, grid_dy;
    grid_dx = 0.2;
    grid_dy = 0.2;
    int nx = (int)(abs(2*grid_x0)/grid_dx) + 1;
    int ny = (int)(abs(2*grid_y0)/grid_dy) + 1;

    double TxxSum = 0.0;
    double TyySum = 0.0;
    double epsP = 0.0;
    fluidCell* fluidCellptr = new fluidCell();
    for (int i = 0; i < nx; i++) {
        // loops over the transverse plane
        double x_local = grid_x0 + i*grid_dx;
        for(int j = 0; j < ny; j++) {
            double y_local = grid_y0 + j*grid_dy;
            hydroinfo_ptr->getHydroinfo(tau_local, x_local, y_local,
                                        fluidCellptr);
            double temp_local = fluidCellptr->temperature;
            if (temp_local > T_dec) {
                double e_local = fluidCellptr->ed;
                double p_local = fluidCellptr->pressure;
                double vx_local = fluidCellptr->vx;
                double vy_local = fluidCellptr->vy;
                double v_perp = sqrt(vx_local*vx_local + vy_local*vy_local);
                double gamma = 1./sqrt(1. - v_perp*v_perp);
                double ux_local = gamma*vx_local;
                double uy_local = gamma*vy_local;
                double pi_xx = fluidCellptr->pi[1][1];
                double pi_yy = fluidCellptr->pi[2][2];
                TxxSum += (e_local+p_local)*ux_local*ux_local + p_local + pi_xx;
                TyySum += (e_local+p_local)*uy_local*uy_local + p_local + pi_yy;
            }
        }
    }
    if (fabs(TxxSum)+ fabs(TyySum) < 1e-10)
        epsP = 0.0;
    else
        epsP = (TxxSum - TyySum)/(TxxSum + TyySum);

    ofstream output;
    output.open("results/tau_vs_epsP.dat", std::ofstream::app);
    output << tau_local-0.6 << "   " << epsP << endl;
    return;
}

void FluidcellStatistic::Countcellvolume(double temp_local, double volume) {
   int temp_idx;
   temp_idx = (int) ((temp_local - t_i)/df);
   if (temp_idx > nbin)  {
       cout << "FluidcellStatistic::Countcellvolume Warning: "
            << "temperature is exceed the upper boundary!" << endl;
       cout << "Temperature = " << temp_local << " GeV. Temp_max = " << t_f 
             << " GeV. " << endl;
   }
   if (temp_idx < 0)  {
       cout << "FluidcellStatistic::Countcellvolume Warning: "
            << "temperature is exceed the lower boundary!" << endl;
       cout << "Temperature = " << temp_local << " GeV. Temp_min = " << t_i 
             << " GeV. " << endl;
   }
   fluidcellvolume_ptr[temp_idx] += volume; 

   return;
}

void FluidcellStatistic::calAvgVvsT(double temp_local, double volume,
                                    double velocity, double ed) {
   double gamma = 1./sqrt(1. - velocity*velocity);
   
   int temp_idx = (int) ((temp_local - t_i)/df);
   if (temp_idx > nbin) {
       cout << "FluidcellStatistic::Countcellvolume Warning: "
            << "temperature is exceed the upper boundary!" << endl;
       cout << "Temperature = " << temp_local << " GeV. Temp_max = " << t_f 
            << " GeV. " << endl;
   }
   if (temp_idx < 0) {
       cout << "FluidcellStatistic::Countcellvolume Warning: "
            << "temperature is exceed the lower boundary!" << endl;
       cout << "Temperature = " << temp_local << " GeV. Temp_min = " << t_i 
             << " GeV. " << endl;
   }

   energyDensity_ptr[temp_idx] += volume*gamma*ed; 
   flowVeclocity_ptr[temp_idx] += volume*gamma*ed*velocity; 

   return;
}

void FluidcellStatistic::OutputCellvolume(string filename) {
    ostringstream filename_stream;
    filename_stream << filename << ".dat";
    ofstream OpCellvolume(filename_stream.str().c_str());

    for (int i=0; i<nbin ; i++) {
        OpCellvolume << scientific << setw(15) << setprecision(8)
                     << tempature_ptr[i] << "   " 
                     << fluidcellvolume_ptr[i] << endl;
    }
    OpCellvolume.close();
    return;
}

void FluidcellStatistic::outputAvgV(string filename) {
    double eps = 1e-10;
    ostringstream filename_stream;
    filename_stream << filename << ".dat";
    ofstream OpAvgV(filename_stream.str().c_str());

    for (int i=0; i<nbin ; i++) {
        double avgV;
        if (energyDensity_ptr[i] < eps)
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

void FluidcellStatistic::outputTempasTauvsX() {
    double grid_t0 = hydroinfo_ptr->getHydrogridTau0();
    double grid_x0 = hydroinfo_ptr->getHydrogridX0();
    double grid_y0 = hydroinfo_ptr->getHydrogridY0();
    double grid_dt = 0.04;
    double grid_dx = 0.2;
    double grid_dy = 0.2;

    int ntime = static_cast<int>(
            (hydroinfo_ptr->getHydrogridTaumax() - grid_t0)/grid_dt) + 1;
    int nx = static_cast<int>(fabs(2.*grid_x0)/grid_dx) + 1;
    int ny = static_cast<int>(fabs(2.*grid_y0)/grid_dy) + 1;

    double tau_local;
    fluidCell* fluidCellptr = new fluidCell();
    ofstream output;
    output.open("results/TempasTauvsX.dat");

    for (int itime = 0; itime < ntime; itime++) {
        // loop over time evolution
        tau_local = grid_t0 + itime*grid_dt;
        for (int i = 0; i < nx; i++) {
            // loops over the transverse plane
            double x_local = grid_x0 + i*grid_dx;
            double y_local = 0.0;
            hydroinfo_ptr->getHydroinfo(tau_local, x_local, y_local,
                                        fluidCellptr);
            double temp_local = fluidCellptr->temperature;
            if (temp_local > 0.05)
                output << temp_local << "   " ;
            else
                output << 0.0 << "   " ;
        }
        output << endl;
    }
    output.close();
    delete fluidCellptr;
    return;
}


void FluidcellStatistic::outputKnudersonNumberasTauvsX() {
    double hbarC = 0.19733;
    double grid_t0 = hydroinfo_ptr->getHydrogridTau0();
    double grid_x0 = hydroinfo_ptr->getHydrogridX0();
    double grid_y0 = hydroinfo_ptr->getHydrogridY0();

    double grid_dt = 0.04;
    double grid_dx = 0.2;
    double grid_dy = 0.2;
    int ntime = static_cast<int>(
            (hydroinfo_ptr->getHydrogridTaumax() - grid_t0)/grid_dt) + 1;
    int nx = static_cast<int>(fabs(2.*grid_x0)/grid_dx) + 1;
    int ny = static_cast<int>(fabs(2.*grid_y0)/grid_dy) + 1;
   
    double eps = 1e-15;
    double MAX = 1000.;

    fluidCell* fluidCellptr = new fluidCell();
    fluidCell* fluidCellptrt1 = new fluidCell();
    fluidCell* fluidCellptrt2 = new fluidCell();
    fluidCell* fluidCellptrx1 = new fluidCell();
    fluidCell* fluidCellptrx2 = new fluidCell();
    fluidCell* fluidCellptry1 = new fluidCell();
    fluidCell* fluidCellptry2 = new fluidCell();

    ofstream output;
    output.open("results/KnudsenNumberasTauvsX.dat");

    for (int itime = 0; itime < ntime; itime++) {
        // loop over time evolution
        double tau_local = grid_t0 + itime*grid_dt;
        for (int i = 0; i < nx; i++) {
            // loops over the transverse plane
            double x_local = grid_x0 + i*grid_dx;
            double y_local = 0.0;
            double grid_dy = grid_dx;
            hydroinfo_ptr->getHydroinfo(tau_local, x_local, y_local,
                                        fluidCellptr);
            hydroinfo_ptr->getHydroinfo(tau_local-grid_dt, x_local, y_local,
                                        fluidCellptrt1);
            hydroinfo_ptr->getHydroinfo(tau_local+grid_dt, x_local, y_local,
                                        fluidCellptrt2);
            hydroinfo_ptr->getHydroinfo(tau_local, x_local-grid_dx, y_local,
                                        fluidCellptrx1);
            hydroinfo_ptr->getHydroinfo(tau_local, x_local+grid_dx, y_local,
                                        fluidCellptrx2);
            hydroinfo_ptr->getHydroinfo(tau_local, x_local, y_local-grid_dy,
                                        fluidCellptry1);
            hydroinfo_ptr->getHydroinfo(tau_local, x_local, y_local+grid_dy,
                                        fluidCellptry2);
       
            double u0 = 1./sqrt(1. - fluidCellptr->vx*fluidCellptr->vx
                                + fluidCellptr->vy*fluidCellptr->vy);
            double u0t1 = 1./sqrt(1. - fluidCellptrt1->vx*fluidCellptrt1->vx
                                  + fluidCellptrt1->vy*fluidCellptrt1->vy);
            double u0t2 = 1./sqrt(1. - fluidCellptrt2->vx*fluidCellptrt2->vx
                                  + fluidCellptrt2->vy*fluidCellptrt2->vy);
            double u1x1 = (fluidCellptrx1->vx
                           /sqrt(1. - fluidCellptrx1->vx*fluidCellptrx1->vx
                                 + fluidCellptrx1->vy*fluidCellptrx1->vy));
            double u1x2 = (fluidCellptrx2->vx
                           /sqrt(1. - fluidCellptrx2->vx*fluidCellptrx2->vx
                                 + fluidCellptrx2->vy*fluidCellptrx2->vy));
            double u2y1 = (fluidCellptry1->vy
                           /sqrt(1. - fluidCellptry1->vx*fluidCellptry1->vx
                                 + fluidCellptry1->vy*fluidCellptry1->vy));
            double u2y2 = (fluidCellptry2->vy
                           /sqrt(1. - fluidCellptry2->vx*fluidCellptry2->vx
                                 + fluidCellptry2->vy*fluidCellptry2->vy));

            double d0u0 = (u0t2 - u0t1)/2./grid_dt;
            double d1u1 = (u1x2 - u1x1)/2./grid_dx;
            double d2u2 = (u2y2 - u2y1)/2./grid_dy;
            double theta = (d0u0 + d1u1 + d2u2 + u0/tau_local);

            double eta_s = 0.08;
            double L_micro = (5.*eta_s
                              /(fabs(fluidCellptr->temperature) + eps))*hbarC;
            double L_macro = 1/(fabs(theta) + eps);
            double Knudsen = L_micro/L_macro;

            output << Knudsen << "    ";
        }
        output << endl;
    }
    output.close();
    delete fluidCellptr;
    delete fluidCellptrt1;
    delete fluidCellptrt2;
    delete fluidCellptrx1;
    delete fluidCellptrx2;
    delete fluidCellptry1;
    delete fluidCellptry2;

    return;
}

void FluidcellStatistic::outputinverseReynoldsNumberasTauvsX() {
    double grid_t0 = hydroinfo_ptr->getHydrogridTau0();
    double grid_x0 = hydroinfo_ptr->getHydrogridX0();
    double grid_y0 = hydroinfo_ptr->getHydrogridY0();
    double grid_dt = 0.04;
    double grid_dx = 0.2;
    double grid_dy = 0.2;
    int ntime = static_cast<int>(
            (hydroinfo_ptr->getHydrogridTaumax() - grid_t0)/grid_dt) + 1;
    int nx = static_cast<int>(fabs(2.*grid_x0)/grid_dx) + 1;
    int ny = static_cast<int>(fabs(2.*grid_y0)/grid_dy) + 1;

    double MAX = 1000.;

    fluidCell* fluidCellptr = new fluidCell();
    ofstream output;
    output.open("results/inverseReynoldsNumberasTauvsX.dat");

    for (int itime = 0; itime < ntime; itime++) {
        // loop over time evolution
        double tau_local = grid_t0 + itime*grid_dt;
        for(int i = 0; i < nx; i++) {
            // loops over the transverse plane
            double x_local = grid_x0 + i*grid_dx;
            double y_local = 0.0;
            hydroinfo_ptr->getHydroinfo(tau_local, x_local, y_local,
                                        fluidCellptr);

            double pi2 =   fluidCellptr->pi[0][0]*fluidCellptr->pi[0][0] 
                         + fluidCellptr->pi[1][1]*fluidCellptr->pi[1][1]
                         + fluidCellptr->pi[2][2]*fluidCellptr->pi[2][2]
                         + fluidCellptr->pi[3][3]*fluidCellptr->pi[3][3]
                         - 2.*(  fluidCellptr->pi[0][1]*fluidCellptr->pi[0][1]
                               + fluidCellptr->pi[0][2]*fluidCellptr->pi[0][2]
                               + fluidCellptr->pi[0][3]*fluidCellptr->pi[0][3])
                         + 2.*(  fluidCellptr->pi[1][2]*fluidCellptr->pi[1][2]
                               + fluidCellptr->pi[1][3]*fluidCellptr->pi[1][3]
                               + fluidCellptr->pi[2][3]*fluidCellptr->pi[2][3]);
       
            double inverseReynold;

            if (pi2 >= 0)
                inverseReynold = sqrt(pi2)/fluidCellptr->pressure;
            else
                inverseReynold = MAX;

            output << inverseReynold << "    " ;
        }
        output << endl;
    }
    output.close();
    delete fluidCellptr;
    return;
}

void FluidcellStatistic::outputBulkinverseReynoldsNumberasTauvsX() {
    double grid_t0 = hydroinfo_ptr->getHydrogridTau0();
    double grid_x0 = hydroinfo_ptr->getHydrogridX0();
    double grid_y0 = hydroinfo_ptr->getHydrogridY0();

    double grid_dt = 0.04;
    double grid_dx = 0.2;
    double grid_dy = 0.2;

    int ntime = static_cast<int>(
            (hydroinfo_ptr->getHydrogridTaumax() - grid_t0)/grid_dt) + 1;
    int nx = static_cast<int>(fabs(2.*grid_x0)/grid_dx) + 1;
    int ny = static_cast<int>(fabs(2.*grid_y0)/grid_dy) + 1;

    fluidCell* fluidCellptr = new fluidCell();
    ofstream output;
    output.open("results/inverseReynoldsNumberasTauvsX.dat");

    for (int itime = 0 ; itime < ntime; itime++) {
        // loop over time evolution
        double tau_local = grid_t0 + itime*grid_dt;
        for (int i = 0; i < nx; i++) {
            // loops over the transverse plane
            double x_local = grid_x0 + i*grid_dx;
            double y_local = 0.0;
            hydroinfo_ptr->getHydroinfo(tau_local, x_local, y_local,
                                        fluidCellptr);
            double inverseReynold;
            inverseReynold = fabs(fluidCellptr->bulkPi)/fluidCellptr->pressure;
            output << inverseReynold << "    " ;
        }
        output << endl;
    }
    output.close();
    delete fluidCellptr;
    return;
}

void FluidcellStatistic::analysis_hydro_volume_for_photon(double T_cut) {
    double V_3 = calculate_hypersurface_3volume(T_cut);
    double V_4 = calculate_spacetime_4volume(T_cut);
    double average_tau = calculate_average_tau(T_cut);
    stringstream output;
    output << "volume_info_for_photon_Tcut_" << T_cut << ".dat";
    ofstream of(output.str().c_str());
    of << "# V_4  <tau>  V_3" << endl;
    of << scientific << setw(18) << setprecision(8)
       << V_4 << "  " << average_tau << "  " << V_3 << endl;
    of.close();
}

double FluidcellStatistic::calculate_spacetime_4volume(double T_cut) {
    // this function calculates the space-time 4 volume of the medium
    // inside a give temperature T_cut [GeV]
    // the output volume V_4 is in [fm^4]
    // deta = 1
   
    // first get hydro grid information
    double grid_t0 = hydroinfo_ptr->getHydrogridTau0();
    double grid_tmax = hydroinfo_ptr->getHydrogridTaumax();
    double grid_x0 = - hydroinfo_ptr->getHydrogridXmax();
    double grid_y0 = - hydroinfo_ptr->getHydrogridYmax();
    double grid_dt = 0.1;
    double grid_dx = 0.2;
    double grid_dy = 0.2;
    int ntime = static_cast<int>((grid_tmax - grid_t0)/grid_dt) + 1;
    int nx = static_cast<int>(fabs(2.*grid_x0)/grid_dx) + 1;
    int ny = static_cast<int>(fabs(2.*grid_y0)/grid_dy) + 1;

    fluidCell* fluidCellptr = new fluidCell();

    double volume = 0.0;
    for (int itime = 0; itime < ntime; itime++) {
        // loop over time evolution
        double tau_local = grid_t0 + itime*grid_dt;
        double volume_element = tau_local*grid_dt*grid_dx*grid_dy;
        for (int i = 0; i < nx; i++) {
            // loops over the transverse plane
            double x_local = grid_x0 + i*grid_dx;
            for (int j = 0; j < ny; j++) {
                double y_local = grid_y0 + j*grid_dy;
                hydroinfo_ptr->getHydroinfo(tau_local, x_local, y_local,
                                            fluidCellptr);
                double T_local = fluidCellptr->temperature;  // GeV
                if (T_local > T_cut) {
                    volume += volume_element;
                }
            }
        }
    }
    delete fluidCellptr;
    return(volume);
}

double FluidcellStatistic::calculate_average_tau(double T_cut) {
    // this function calculates the average <\tau> of the medium
    // inside a give temperature T_cut [GeV]
    // the output <tau> is in [fm]
   
    // first get hydro grid information
    double grid_t0 = hydroinfo_ptr->getHydrogridTau0();
    double grid_tmax = hydroinfo_ptr->getHydrogridTaumax();
    double grid_x0 = - hydroinfo_ptr->getHydrogridXmax();
    double grid_y0 = - hydroinfo_ptr->getHydrogridYmax();
    double grid_dt = 0.1;
    double grid_dx = 0.2;
    double grid_dy = 0.2;
    int ntime = static_cast<int>((grid_tmax - grid_t0)/grid_dt) + 1;
    int nx = static_cast<int>(fabs(2.*grid_x0)/grid_dx) + 1;
    int ny = static_cast<int>(fabs(2.*grid_y0)/grid_dy) + 1;

    fluidCell* fluidCellptr = new fluidCell();

    double average_tau = 0.0;
    double volume = 0.0;
    for (int itime = 0; itime < ntime; itime++) {
        // loop over time evolution
        double tau_local = grid_t0 + itime*grid_dt;
        double volume_element = tau_local*grid_dt*grid_dx*grid_dy;
        for (int i = 0; i < nx; i++) {
            // loops over the transverse plane
            double x_local = grid_x0 + i*grid_dx;
            for (int j = 0; j < ny; j++) {
                double y_local = grid_y0 + j*grid_dy;
                hydroinfo_ptr->getHydroinfo(tau_local, x_local, y_local,
                                            fluidCellptr);
                double T_local = fluidCellptr->temperature;  // GeV
                if (T_local > T_cut) {
                    volume += volume_element;
                    average_tau += tau_local*volume_element;
                }
            }
        }
    }
    average_tau /= volume;
    delete fluidCellptr;
    return(average_tau);
}

double FluidcellStatistic::calculate_hypersurface_3volume(double T_cut) {
    // this function calculates the surface area of an isothermal hyper-surface
    // at a give temperature T_cut [GeV]
    // the output volume V_3 = int u^\mu d^3\sigma_\mu is in [fm^3]
    // deta = 1
   
    // first find the hyper-surface
    SurfaceFinder surfFinder(hydroinfo_ptr, paraRdr, T_cut);
    surfFinder.Find_full_hypersurface();

    // load the hyper-surface file
    ifstream surf("hyper_surface_2+1d.dat");
    if (!surf.is_open()) {
        cout << "FluidcellStatistic::calculate_hypersurface_3volume:"
             << "can not open the file hyper_surface_2+1d.dat!" << endl;
        exit(1);
    }
    
    double volume = 0.0;
    double tau, x, y, da0, da1, da2, vx, vy, T;
    surf >> tau >> x >> y >> da0 >> da1 >> da2 >> T >> vx >> vy;
    while (!surf.eof()) {
        double u0 = 1./sqrt(1. - vx*vx - vy*vy);
        double ux = u0*vx;
        double uy = u0*vy;
        volume += tau*(u0*da0 + ux*da1 + uy*da2);
        surf >> tau >> x >> y >> da0 >> da1 >> da2 >> T >> vx >> vy;
    }
    return(volume);
}

