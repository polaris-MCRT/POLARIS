/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CGRID_CYLINCDRICAL_H
#define CGRID_CYLINCDRICAL_H

#include "GridBasic.hpp"
#include "Typedefs.hpp"
#include "Vector3D.hpp"
#include "CellCylindrical.hpp"
#include "Photon.hpp"

class CGridCylindrical : public CGridBasic
{
public:
    CGridCylindrical(void)
    {
        // geometry & class-specific members
        Rmin = 0;
        Rmax = 1;
        Zmax = 1;
        N_r = 0;
        N_z = 0;
        N_ph = 0;
        log_factorR  = 0;
        log_factorPh = 0;
        log_factorZ  = 0;

        grid_cells   = 0;
        center_cells = 0;

        listR = 0;
        listPh = 0;
        listZ = 0;
    }

    ~CGridCylindrical()
    {
        if(grid_cells != 0)
        {
            for(uint i_r = 0; i_r < N_r; i_r++)
            {
                // cout << "Cleaning memory for cylindrical grid file: "
                //      << float(100.0 * double(i_r) / double(N_r)) << "      \r" << flush;

                for(uint i_ph = 0; i_ph < N_ph[i_r]; i_ph++)
                {
                    for(uint i_z = 0; i_z < N_z; i_z++)
                    {
                        delete grid_cells[i_r][i_ph][i_z];
                        grid_cells[i_r][i_ph][i_z] = 0;
                    }

                    delete[] grid_cells[i_r][i_ph];
                    grid_cells[i_r][i_ph] = 0;
                }

                delete[] grid_cells[i_r];
                grid_cells[i_r] = 0;
            }
            delete[] grid_cells;
        }

        if(center_cells != 0)
        {
            for(uint i_z = 0; i_z < N_z; i_z++)
            {
                delete center_cells[i_z];
                center_cells[i_z] = 0;
            }
            delete[] center_cells;
        }

        if(cell_list != 0)
        {
            delete[] cell_list;
            cell_list = 0;
        }

        if(listR != 0)
        {
            delete[] listR;
            listR = 0;
        }

        if(listPh != 0)
        {
            for(uint i_r = 0; i_r < N_r; i_r++)
                delete[] listPh[i_r];
            delete[] listPh;
            listPh = 0;
        }

        if(listZ != 0)
        {
            for(uint i_r = 0; i_r < N_r; i_r++)
                delete[] listZ[i_r];
            delete[] listZ;
            listZ = 0;
        }

        if(N_ph != 0)
        {
            delete[] N_ph;
            N_ph = 0;
        }

        // cout << CLR_LINE << flush;
    }


    bool goToNextCellBorder(photon_package * pp);
    bool updateShortestDistance(photon_package * pp);

    Vector3D getCenter(const cell_basic & cell) const;

    bool next(photon_package * pp);

    bool findStartingPoint(photon_package * pp);

    void getLengths(uint bins, double & step_xy, double & off_xy);

    bool createCellList();

    double getVolume(const cell_basic & cell) const;

    /*
    This routine rotates a vector with direction "dir" from the current position of the photon package
    into the center of the current cell or vice versa (inv = true).
    Useful if the complete model space is symmetrical with respect to a coordinate,
    e.g. one phi cell, star at (0,0,0), magnetic field along z-axis.
    */
    Vector3D rotateToCenter(const photon_package & pp, Vector3D dir, bool inv, bool phi_only) const;

    bool positionPhotonInGrid(photon_package * pp);

    bool createArtificialGrid(string path);

    bool saveBinaryGridFile(string filename);

    bool loadGridFromBinaryFile(parameters & param, uint _data_len);
    bool saveBinaryGridFile(string filename, ushort id, ushort data_size);

    bool loadGridFromBinaryFile(parameters & param);

    void clear();

    void printParameters();

    bool getPolarRTGridParameter(double max_len,
                                 double pixel_width,
                                 uint max_subpixel_lvl,
                                 dlist & _listR,
                                 uint & N_polar_r,
                                 uint *& N_polar_ph);

private:
    double Rmin, Rmax, Zmax;
    uint N_r, N_z;
    uint * N_ph;
    double log_factorR, log_factorPh, log_factorZ;

    double * listR;
    double ** listPh;
    double ** listZ;

    cell_cyl **** grid_cells;
    cell_cyl ** center_cells;

    bool isInside(const Vector3D & pos) const;

    void setRndPositionInCell(photon_package * pp, CRandomGenerator * rand_gen);
};

#endif /* CGRID_CYLINCDRICAL_H */
