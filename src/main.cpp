#include "Pipeline.h"
#include "OPIATE.h"

int main(int argc, char ** argv)
{
    OPIATE op;
    
    op.readFitsData("/home/ilion/1/reisslst/polaris projects/OPIATE/data/test.fits");
    
    cout << "\n";
    cout << op.biListIndexDataSearch(200) << "\n" << flush;
    
    return 0;
    
    CPipeline pipeline;
    if(!pipeline.Init(argc, argv))
        return 0;

    pipeline.Run();

    return 0;
}