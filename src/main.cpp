#include "Pipeline.h"
//#include "Grid.h"
//#include "OcTree.h"

int main(int argc, char ** argv)
{
    /*CGridOcTree tree;
    
    parameters param;
    omp_set_num_threads(7);
    string str_input = "/home/ilion/1/reisslst/polaris projects/Michael/grid_2species_lvl21_root7_403_upd.dat";
    string str_output = "/home/ilion/1/reisslst/polaris projects/Michael/grid_2species_lvl21_root7_403_upd_sub_7.dat";
    
    param.setCommand(CMD_DUST_EMISSION);
    param.setPathGrid(str_output);
    tree.loadGridFromBinrayFile(param);
    tree.printParameters();
    tree.printPhysicalParameters();
    
    tree.reduceBinrayFile(str_input, str_output,9);
    return 0;/**/
        
    CPipeline pipeline;
    if(!pipeline.Init(argc, argv))
        return 0;

    pipeline.Run();

    return 0;
}