#include "Pipeline.h"
#include "Grid.h"
#include "OcTree.h"


int main(int argc, char ** argv)
{
    /*CGridOcTree tree;
    
    parameters param;
    
    string str_input = "/home/ilion/1/reisslst/polaris projects/Michael/grid_2species_lvl21_root7_403.dat";
    string str_output = "/home/ilion/1/reisslst/polaris projects/Michael/grid_level8.dat";
    
    tree.reduceBinrayFile(str_input, str_output,9);
    return 0;/**/
        
    CPipeline pipeline;
    if(!pipeline.Init(argc, argv))
        return 0;

    pipeline.Run();

    return 0;
}