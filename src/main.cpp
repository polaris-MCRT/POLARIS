#include "Pipeline.h"
#include "OPIATE.h"

int main(int argc, char ** argv)
{
    COpiateDataBase op;
    
    op.readEmissivityData("/mnt/c/Users/Stefan/Documents/work/opiate/test.fits");
    op.readExtincitonData("/mnt/c/Users/Stefan/Documents/work/opiate/test_ext.fits");
    
    cout << "\n";
    cout << op.biListIndexDataSearch(8) << "\n" << flush;
    
    op.findIndexByName("MOL50");
    
    cout << flush;
    return 0;
    
    CPipeline pipeline;
    if(!pipeline.Init(argc, argv))
        return 0;

    pipeline.Run();

    return 0;
}