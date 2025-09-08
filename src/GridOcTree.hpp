/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CGRID_OCTREE_H
#define CGRID_OCTREE_H

#include "GridBasic.hpp"
#include "RandomGenerator.hpp"
#include "Matrix2D.hpp"
#include "Typedefs.hpp"
#include "Vector3D.hpp"
#include "CellOcTree.hpp"
#include "Photon.hpp"

class CGridOcTree : public CGridBasic
{
public:
    CGridOcTree(void)
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
        //data_pos_ra = MAX_UINT;

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

        cell_oc_root = 0;
        cell_oc_pos = 0;
        rec_counter = 0;
        max_data = data_offset;
        max_level = 0;

        line_counter = 0;
        char_counter = 0;

        nrOfDensRatios = 0;

        rot_angle1 = 0;
        rot_angle2 = 0;
    }

    ~CGridOcTree()
    {
        if(cell_oc_root == 0)
            return;

        clear(cell_oc_root);

        if(cell_list != 0)
        {
            delete[] cell_list;
            cell_list = 0;
        }

        cout << CLR_LINE;
    }

    // begin IO functions
    bool writePlotFiles(string path, parameters & param);

    void goToRoot();

    bool nextLowLevelCell();
    bool nextLowLevelCell(cell_basic * cell);
    // end   IO functions
    bool goToNextCellBorder(photon_package * pp);
    bool updateShortestDistance(photon_package * pp);

    // void deleteSubCells(slist & source);

    bool reduceBinaryFile(string in_filename, string out_filename, uint tr_level);
    bool reduceLevelOfBinaryFile(cell_oc * cell, uint tr_level);

    Vector3D getCenter(const cell_basic & cell) const;

    Vector3D getMidplaneCenter(cell_basic * cell);

    bool createCellList();

    bool findMatchingCell(photon_package * pp);

    bool next(photon_package * pp);

    bool findStartingPoint(photon_package * pp);

    void getLengths(uint bins, double & step_xy, double & off_xy);

    double getVolume(const cell_basic & cell) const;

    bool positionPhotonInGrid(photon_package * pp);

    const cell_oc * getTopLevelCell() const;

    const cell_oc * getCurrentCell() const;

    void printParameters();

    bool createArtificialGrid(string path);

    void createNextLevel(cell_oc * cell);

    bool saveBinaryGridFile(string filename);

    bool loadGridFromBinaryFile(parameters & param, uint data_len);

    bool saveBinaryGridFile(string filename, ushort id, ushort data_size);

    bool loadGridFromBinaryFile(parameters & param);

    void nextBinaryDataCell(ofstream & file_stream, cell_oc * cell, uint data_size);

    void clear();

    void goNextLevelDown(photon_package * pp);

    bool createTree(cell_oc * parent,
                    double _x_min,
                    double _y_min,
                    double _z_min,
                    double _length,
                    uint _level);

    bool initiateTreeFromFile(uint _nx,
                              uint _max_level,
                              double _fa,
                              double _length,
                              string str_dens,
                              string str_temp,
                              string str_magx,
                              string str_magy,
                              string str_magz);

private:
    void clear(cell_oc * cell);
    cell_oc * cell_oc_root;
    cell_oc * cell_oc_pos;

    Matrix3D datdens;
    Matrix3D dattemp;
    Matrix3D datmx;
    Matrix3D datmy;
    Matrix3D datmz;

    uint nx;
    // uint ny;
    // uint nz;
    double f_min;
    double f_max;
    double factor;
    uint treelevel_counter;
    uint tagged_cells;

    uint rec_counter;
    double max_level;

    void plotNextDataPoint(ofstream * file_streams, cell_oc * cell, uint level);
    void plotNextDataVector(ofstream * file_streams, cell_oc * cell, uint level);
    void plotNextGridCell(ofstream * grid_streams, cell_oc * cell, uint level);

    void createBoundingCell();
    void createBoundingCell(cell_oc * cell);

    void goNextLevelUp(photon_package * pp);

    bool isInside(const Vector3D & pos, const cell_basic & _cell) const;

    bool isInside(const Vector3D & pos) const;

    // bool isInside(photon_package * pp, Vector3D & pos);

    void setRndPositionInCell(photon_package * pp, CRandomGenerator * rand_gen);
};

#endif /* CGRID_OCTREE_H */
