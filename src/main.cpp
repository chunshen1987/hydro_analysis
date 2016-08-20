/////////////////////////////////////////////////////////////////////////
//                      hydrodynamics analysis
//
//              author: Chun Shen <chunshen@physics.mcgill.ca>
//              Copyright: Chun Shen 2014
//
//  This program load hydrodynamic evolution files and perform various
//  kinds of analysis
//  
//
//  To do in the future:
/////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "./Hydroinfo_h5.h"
#include "./Stopwatch.h"
#include "./FluidcellStatistic.h"
#include "./ParameterReader.h"
#include "./SurfaceFinder.h"

using namespace std;

int main(int argc, char *argv[]) {
    ParameterReader *paraRdr = new ParameterReader();
    paraRdr->readFromFile("parameters.dat");
    paraRdr->readFromArguments(argc, argv);
    paraRdr->echo();

    int load_viscous = paraRdr->getVal("load_viscous_info");

    Stopwatch sw;

    sw.tic();
    // hydro data file pointer
    HydroinfoH5* hydroinfo_ptr = new HydroinfoH5("JetData.h5", 500,
                                                 load_viscous);

    FluidcellStatistic fluidcellanalysis(hydroinfo_ptr, paraRdr);
    double T_cut = paraRdr->getVal("T_cut");
    fluidcellanalysis.analysis_hydro_volume_for_photon(T_cut);
    fluidcellanalysis.output_temperature_vs_avg_utau();

    // construct freeze-out hyper-surface
    // SurfaceFinder* surface_ptr = new SurfaceFinder(hydroinfo_ptr, paraRdr);
    // surface_ptr->Find_full_hypersurface();

    sw.toc();
    cout << "totally takes : " << sw.takeTime() << " seconds." << endl;

    return(0);
}


