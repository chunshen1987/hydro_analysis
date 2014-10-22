#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>
#include<fstream>

#include "SurfaceFinder.h"
#include "cornelius.h"

using namespace std;

SurfaceFinder::SurfaceFinder(HydroinfoH5* hydroinfo_ptr_in, ParameterReader* paraRdr_in)
{
    hydroinfo_ptr = hydroinfo_ptr_in;
    paraRdr = paraRdr_in;

    T_cut = paraRdr->getVal("T_cut");
}

SurfaceFinder::~SurfaceFinder()
{

}

bool SurfaceFinder::check_intersect(double T_cut, double tau, double x, double y, double dt, double dx, double dy, double ***cube)
{
    fluidCell *fluidCellptr = new fluidCell();
    bool intersect = true;

    double tau_low = tau - dt/2.;
    double tau_high = tau + dt/2.;
    double x_left = x - dx/2.;
    double x_right = x + dx/2.;
    double y_left = y - dy/2.;
    double y_right = y + dy/2.;

    hydroinfo_ptr->getHydroinfo(tau_low, x_left, y_left, fluidCellptr);
    cube[0][0][0] = fluidCellptr->temperature;
    hydroinfo_ptr->getHydroinfo(tau_low, x_left, y_right, fluidCellptr);
    cube[0][0][1] = fluidCellptr->temperature;
    hydroinfo_ptr->getHydroinfo(tau_low, x_right, y_left, fluidCellptr);
    cube[0][1][0] = fluidCellptr->temperature;
    hydroinfo_ptr->getHydroinfo(tau_low, x_right, y_right, fluidCellptr);
    cube[0][1][1] = fluidCellptr->temperature;
    hydroinfo_ptr->getHydroinfo(tau_high, x_left, y_left, fluidCellptr);
    cube[1][0][0] = fluidCellptr->temperature;
    hydroinfo_ptr->getHydroinfo(tau_high, x_left, y_right, fluidCellptr);
    cube[1][0][1] = fluidCellptr->temperature;
    hydroinfo_ptr->getHydroinfo(tau_high, x_right, y_left, fluidCellptr);
    cube[1][1][0] = fluidCellptr->temperature;
    hydroinfo_ptr->getHydroinfo(tau_high, x_right, y_right, fluidCellptr);
    cube[1][1][1] = fluidCellptr->temperature;

    if((T_cut - cube[0][0][0])*(cube[1][1][1] - T_cut) < 0.0)
        if((T_cut - cube[0][1][0])*(cube[1][0][1] - T_cut) < 0.0)
            if((T_cut - cube[0][1][1])*(cube[1][0][0] - T_cut) < 0.0)
                if((T_cut - cube[0][0][1])*(cube[1][1][0] - T_cut) < 0.0)
                    intersect = false;

    delete fluidCellptr;
    return(intersect);
}

int SurfaceFinder::Find_full_hypersurface()
{
    ofstream output;
    output.open("hyper_surface_2+1d.dat");

    double grid_t0 = paraRdr->getVal("grid_t0");
    double grid_x0 = paraRdr->getVal("grid_x0");
    double grid_y0 = paraRdr->getVal("grid_y0");

    double grid_dt = paraRdr->getVal("grid_dt");
    double grid_dx = paraRdr->getVal("grid_dx");
    double grid_dy = paraRdr->getVal("grid_dy");

    int dim = 3;
    double *lattice_spacing = new double [dim];
    lattice_spacing[0] = grid_dt;
    lattice_spacing[1] = grid_dx;
    lattice_spacing[2] = grid_dy;

    Cornelius* cornelius_ptr = new Cornelius();
    cornelius_ptr->init(dim, T_cut, lattice_spacing);
  
    int ntime = (int)((hydroinfo_ptr->getHydrogridTaumax() - grid_t0)/grid_dt);
    int nx = (int)(abs(2*grid_x0)/grid_dx);
    int ny = (int)(abs(2*grid_y0)/grid_dy);

    double ***cube = new double** [2];
    for(int i = 0; i < 2; i++)
    {
        cube[i] = new double* [2];
        for(int j = 0; j < 2; j++)
        {
            cube[i][j] = new double [2];
            for(int k = 0; k < 2; k++)
                cube[i][j][k] = 0.0;
        }
    }
    
    fluidCell *fluidCellptr = new fluidCell();
  

    ofstream output_test1("decdat2.dat");
    ofstream output_test2("surface.dat");
    for(int itime = 0; itime < ntime; itime++) //loop over time evolution
    {
        double tau_local = grid_t0 + (itime + 0.5)*grid_dt;
        for(int i = 0; i < nx; i++) //loops over the transverse plane
        {
            double x_local = grid_x0 + (i + 0.5)*grid_dx;
            for(int j = 0; j < ny; j++)
            {
                double y_local = grid_y0 + (j + 0.5)*grid_dy;
                bool intersect = check_intersect(T_cut, tau_local, x_local, y_local, grid_dt, grid_dx, grid_dy, cube);
                if(intersect)
                {
                    cornelius_ptr->find_surface_3d(cube);
                    for(int isurf = 0; isurf < cornelius_ptr->get_Nelements(); isurf++)
                    {
                        double tau_center = cornelius_ptr->get_centroid_elem(isurf, 0) + tau_local - grid_dt/2.;
                        double x_center = cornelius_ptr->get_centroid_elem(isurf, 1) + x_local - grid_dx/2.;
                        double y_center = cornelius_ptr->get_centroid_elem(isurf, 2) + y_local - grid_dy/2.;

                        double da_tau = cornelius_ptr->get_normal_elem(isurf, 0);
                        double da_x = cornelius_ptr->get_normal_elem(isurf, 1);
                        double da_y = cornelius_ptr->get_normal_elem(isurf, 2);
                       
                        hydroinfo_ptr->getHydroinfo(tau_center, x_center, y_center, fluidCellptr);

                        output << scientific << setw(18) << setprecision(8) 
                               << tau_center << "   " << x_center << "   " << y_center << "   " 
                               << da_tau << "   " << da_x << "   " << da_y << "   " 
                               << fluidCellptr->temperature << "   " << fluidCellptr->vx << "   " << fluidCellptr->vy 
                               << endl;
                        
                        output_test1 << tau_center << "   " << da_tau << "    " << da_x << "   " << da_y << "   "
                                     << fluidCellptr->vx << "   " << fluidCellptr->vy << "   " << 0.18 << "   " 
                                     << 0.0 << "   " << 0.11956 << "   " << 0.0 << "   " << 0.0 << "   " 
                                     << 0.026813212 << "   " << fluidCellptr->pi[3][3] << "   " << fluidCellptr->pi[0][0] << "   "
                                     << fluidCellptr->pi[0][1] << "   " << fluidCellptr->pi[0][2] << "   " 
                                     << fluidCellptr->pi[1][1] << "   " << fluidCellptr->pi[1][2] << "   " 
                                     << fluidCellptr->pi[2][2] << "   " << 0.0 << endl;

                        output_test2 << tau_center << "   " << tau_center << "   " 
                                     << x_center << "   " << y_center << "   " 
                                     << x_center << "   " << x_center << "   " << x_center << endl;

                    }
                }
            }
        }
    }
    output.close();
    output_test1.close();
    output_test2.close();
    
    delete fluidCellptr;
    delete cornelius_ptr;
    delete [] lattice_spacing;
    for(int i = 0; i < 2; i++)
    {
        for(int j = 0; j < 2; j++)
            delete [] cube[i][j];
        delete [] cube[i];
    }
    delete [] cube;
    return 0;
}
