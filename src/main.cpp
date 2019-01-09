#include "Pipeline.h"

int main(int argc, char ** argv)
{
    CPipeline pipeline;
    if(!pipeline.Init(argc, argv))
        return 0;

    pipeline.Run();

    return 0;
}