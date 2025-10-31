/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                      Copyright (C) 2018 - 2026 Stefan Reissl                      *
************************************************************************************/
#include "Pipeline.hpp"

int main(int argc, char ** argv)
{
    CPipeline pipeline;
    if(!pipeline.Init(argc, argv))
        return 0;

    pipeline.Run();

    return 0;
}
