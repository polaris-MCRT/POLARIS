/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CGRID_SPHERICAL_H
#define CGRID_SPHERICAL_H

#include "GridBasic.hpp"
#include "Typedefs.hpp"
#include "Vector3D.hpp"
#include "CellSpherical.hpp"
#include "Photon.hpp"

class CGridSpherical : public CGridBasic
{
public:
    CGridSpherical(void)
    {
        // spherical grid geometry & class-specific members
        Rmin = 0;
        Rmax = 1;
        N_r  = 4;
        N_ph = 4;
        N_th = 3;
        log_factorR  = 0;
        log_factorPh = 0;
        log_factorTh = 0;

        grid_cells  = 0;
        center_cell = 0;

        listR = 0;
        listPh = 0;
        listTh = 0;
    }

    ~CGridSpherical()
    {
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
            delete[] listPh;
            listPh = 0;
        }

        if(listTh != 0)
        {
            delete[] listTh;
            listTh = 0;
        }

        if(center_cell != 0)
        {
            delete center_cell;
            center_cell = 0;
        }

        if(grid_cells != 0)
        {
            for(uint i_r = 0; i_r < N_r; i_r++)
            {
                // cout << "Cleaning memory for spherical grid file: "
                //      << float(100.0 * double(i_r) / double(N_r)) << "      \r" << flush;

                for(uint i_ph = 0; i_ph < N_ph; i_ph++)
                {
                    for(uint i_th = 0; i_th < N_th; i_th++)
                    {
                        delete grid_cells[i_r][i_ph][i_th];
                        grid_cells[i_r][i_ph][i_th] = 0;
                    }

                    delete[] grid_cells[i_r][i_ph];
                    grid_cells[i_r][i_ph] = 0;
                }

                delete[] grid_cells[i_r];
                grid_cells[i_r] = 0;
            }
            delete[] grid_cells;
            grid_cells = 0;
        }

        // cout << CLR_LINE << flush;
    }

    

    bool goToNextCellBorder(photon_package * pp);
    bool updateShortestDistance(photon_package * pp);

    Vector3D getCenter(const cell_basic & cell) const;

    void setRndPositionInCell(photon_package * pp, CRandomGenerator * rand_gen);

    bool next(photon_package * pp);

    /*
    void getBoundingPoints(Vector3D & p_min, Vector3D & p_max);

    void getBoundingPoints(cell_basic * cell, Vector3D & p_min, Vector3D & p_max);
    */

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
    double Rmin, Rmax;
    uint N_r, N_ph, N_th;
    double log_factorR, log_factorPh, log_factorTh;

    double * listR;
    double * listPh;
    double * listTh;

    cell_sp **** grid_cells;
    cell_sp * center_cell;

    bool isInside(const Vector3D & pos) const;
};

#endif /* CGRID_SPHERICAL_H */
