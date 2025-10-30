/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CELL_VORONOI_H
#define CELL_VORONOI_H

#include "CellBasic.hpp"

class cell_vo : public cell_basic
{
public:
    cell_vo()
    {
        id = 0;
        data = 0;
        nr_neighbors = 0;
        neighbors = 0;
        volume = -1;
        ref_length=0;
    }

    ~cell_vo()
    {
        if(data != 0)
            delete[] data;

        data = 0;

        if(neighbors != 0)
            delete[] neighbors;

        neighbors = 0;
    }

    void initNeighbors(short nr);

    void setNeighbor(uint pos, int id);

    void setCenter(double cx, double cy, double cz);

    void setVolume(double v);

    Vector3D getCenter() const;

    double getX();

    double getY();

    double getZ();

    double getVolume() const;
    
    double getRefLength() const;
    
    int getNeighborID(uint pos);

    ushort getNrOfNeighbors();

private:
    Vector3D center;
    ushort nr_neighbors;
    int * neighbors;
    double volume;
    double ref_length;
};

#endif /* CELL_VORONOI_H */
