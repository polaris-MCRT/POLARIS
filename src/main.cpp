//#include "Synchrotron.h"
#include "Pipeline.h"

int main(int argc, char** argv)
{      
    /*CSynchrotron syn;
    syn_param p_syn;
    
    double n_e=1.94677151093791051e-6;
    double B=2.828427044752808e-05*1e4;
    double T_e=0.0;
    double l=con_c/408e9;
    double theta=0.5480284076203128;
    
    
    p_syn=syn.get_Thermal_Parameter(n_e,T_e,l,B,theta);
    
    cout << "kappa_Q: " << p_syn.kappa_Q << endl; 
    cout << "kappa_V: " << p_syn.kappa_V << endl; 
    return 0; /**/
    
    CPipeline pipeline;
    if(!pipeline.Init(argc, argv))
        return 0;

    pipeline.Run();

    return 0;
}