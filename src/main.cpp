/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "Pipeline.hpp"
#include "FreeFree.hpp"

int main(int argc, char ** argv)
{
    /*string filename= "/mnt/f/work/free/gaunt_cube_with_axes_and_Z.fits";
    CFreeFree free;
    free.loadGauntFITS(filename);
    free.printParameters();*/
    
    /*float Z = 6.3f;
    float log_gam2 = -4.5f;
    float log_u = -10.7f;/**/

    /*float Z = 36.0f;
    float log_gam2 = -6.0f;
    float log_u = -16.0f;/**/
    
    /*float Z = 25.0f;
    float log_gam2 = 0.0f;
    float log_u = 0.0f;/**/

    /*float value = free.interpolate(Z, log_gam2, log_u);
    std::cout << "Interpolated Gaunt factor: " << value << std::endl;/**/
    
    
    CPipeline pipeline;
    if(!pipeline.Init(argc, argv))
        return 0;

    pipeline.Run();

    return 0;
}
