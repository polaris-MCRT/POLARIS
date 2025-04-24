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
        basic_path = 0;
        buffer_size = 0;

        max_cells = 0;
        max_value = 0;
        max_data = 0;

        /*min_delta = 0;
        max_delta = 0;*/

        min_mach = 1e300;
        max_mach = 0;

        min_mag = 0;
        max_mag = 0;

        min_vel = 0;
        max_vel = 0;

        min_len = 0;
        max_len = 0;

        min_gas_temp = 0;
        max_gas_temp = 0;

        min_dust_temp = 0;
        max_dust_temp = 0;

        min_gas_dens = 0;
        max_gas_dens = 0;

        min_dust_dens = 0;
        max_dust_dens = 0;

        aalg_min = 1e300;
        aalg_max = 0;

        min_larm_limit = 0;
        max_larm_limit = 0;

        min_pres = 0;
        max_pres = 0;

        line_counter = 0;
        char_counter = 0;
        ru[0] = '|';
        ru[1] = '/';
        ru[2] = '-';
        ru[3] = '\\';

        conv_length_in_SI = 1;
        conv_dens_in_SI = 1;
        conv_Bfield_in_SI = 1;
        conv_Vfield_in_SI = 1;

        nrOfPlotPoints = 1000;
        nrOfPlotVectors = 1000;
        maxPlotLines = 3;

        cell_list = 0;

        data_offset = 6;
        dataID = 0;
        data_len = 0;

        total_gas_mass = 0;
        mu = 2.0;

        nrOfDensRatios = 0;
        nrOfOpiateIDs = 0;

        // data_pos_gd = MAX_UINT;
        // data_pos_dd = MAX_UINT;
        // data_pos_td = MAX_UINT;
        data_pos_tg = MAX_UINT;
        data_pos_mx = MAX_UINT;
        data_pos_my = MAX_UINT;
        data_pos_mz = MAX_UINT;
        data_pos_vx = MAX_UINT;
        data_pos_vy = MAX_UINT;
        data_pos_vz = MAX_UINT;
        data_pos_px = MAX_UINT;
        data_pos_py = MAX_UINT;
        data_pos_pz = MAX_UINT;
        // data_pos_aalg = MAX_UINT;
        data_pos_amin = MAX_UINT;
        data_pos_amax = MAX_UINT;
        data_pos_size_param = MAX_UINT;
        data_pos_ra = MAX_UINT;

        data_pos_vt = MAX_UINT;
        data_pos_pda = MAX_UINT;

        data_pos_op = MAX_UINT;

        data_pos_n_th = MAX_UINT;
        data_pos_T_e = MAX_UINT;
        data_pos_n_cr = MAX_UINT;
        data_pos_g_min = MAX_UINT;
        data_pos_g_max = MAX_UINT;
        data_pos_p = MAX_UINT;

        pos_GasSpecRatios = 0;
        pos_OpiateIDS = 0;

        rot_angle1 = 0;
        rot_angle2 = 0;

        plt_gas_dens = false;
        plt_mol_dens = false;
        plt_dust_dens = false;
        plt_gas_temp = false;
        plt_dust_temp = false;
        plt_mag = false;
        plt_vel = false;
        plt_rat = false;
        //plt_delta = false;
        plt_larm = false;
        //plt_mach = false;
        plt_dust_id = false;

        total_volume = 0;
        cell_volume = 0;

        Rmin = 0;
        Rmax = 1;
        N_r = 4;
        N_ph = 4;
        N_th = 3;
        log_factorR = 0;
        log_factorPh = 0;
        log_factorTh = 0;

        grid_cells = 0;
        center_cell = 0;

        cell_list = 0;

        listR = 0;
        listPh = 0;
        listTh = 0;

        line_counter = 0;
        char_counter = 0;

        rot_angle1 = 0;
        rot_angle2 = 0;
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

    bool writePlotFiles(string path, parameters & param);

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
