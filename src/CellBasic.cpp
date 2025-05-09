/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "CellBasic.hpp"
#include "GridBasic.hpp"

bool cell_basic::isValid()
{
    return data != 0;
}

double cell_basic::getData(uint i) const
{
    if(data == 0)
        return 0;

    if(i == MAX_UINT)
        return 0;

    //if(i > 4138637280)
    //    return 0;

    return data[i];
}

void cell_basic::setData(uint i, double d)
{
    #pragma omp atomic write
    data[i] = d;
}

void cell_basic::updateData(uint i, double d)
{
    #pragma omp atomic update
    data[i] += d;
}

void cell_basic::convertData(uint i, double c)
{
    #pragma omp atomic update
    data[i] *= c;
}

void cell_basic::resize(uint size)
{
    if(data != 0)
        delete[] data;

    data = new double[size];
    for(uint i = 0; i < size; i++)
        data[i] = 0;
}

void cell_basic::setID(uint _id)
{
    id = _id;
}

ulong cell_basic::getUniqueID() const
{
    return ulong(id);
}

void cell_basic::updateID(uint _id)
{
    id += _id;
}
