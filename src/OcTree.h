#pragma once
#include <iostream>
#include <string>

#include "Grid.h"
#include "Matrix2D.h"
#include "Typedefs.h"
#include "Vector.h"
#include "Cell.h"
#include "Photon.h"

class parameters;

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

        min_delta = 0;
        max_delta = 0;

        min_mach = 0;
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
        plt_delta = false;
        plt_larm = false;
        plt_mach = false;
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

    void goToRoot()
    {
        cell_oc_pos = cell_oc_root;
    }

    bool nextLowLevelCell();
    bool nextLowLevelCell(cell_basic * cell);
    // end   IO functions
    bool goToNextCellBorder(photon_package * pp);
    bool updateShortestDistance(photon_package * pp);

    // void deleteSubCells(slist & source);

    bool reduceBinaryFile(string in_filename, string out_filename, uint tr_level);
    bool reduceLevelOfBinaryFile(cell_oc * cell, uint tr_level);

    Vector3D getCenter(const cell_basic & cell) const
    {
        Vector3D center;
        const cell_oc * tmp_cell = (const cell_oc *)&cell;

        center.setX(tmp_cell->getXmin() + 0.5 * tmp_cell->getLength());
        center.setY(tmp_cell->getYmin() + 0.5 * tmp_cell->getLength());
        center.setZ(tmp_cell->getZmin() + 0.5 * tmp_cell->getLength());

        return center;
    }

    Vector3D getMidplaneCenter(cell_basic * cell)
    {
        Vector3D center;
        cell_oc * tmp_cell = (cell_oc *)cell;

        center.setX(tmp_cell->getXmin() + 0.5 * tmp_cell->getLength());
        center.setY(tmp_cell->getYmin() + 0.5 * tmp_cell->getLength());
        center.setZ(0);

        return center;
    }

    double maxLength()
    {
        return cell_oc_root->getLength();
    }

    bool createCellList()
    {
        if(max_cells == 0)
        {
            cout << "\nERROR: OcTree grid contains no cells!" << endl;
            cout << "       Cell list cannot be created!" << endl;
            return false;
        }

        cell_list = new cell_basic *[max_cells];
        ulong pos_counter = 0;
        goToRoot();
        cout << CLR_LINE;
        cout << "-> Creating cell list    : 0 [%]           \r";

        while(nextLowLevelCell())
        {
#pragma warning(suppress : 6386)
            /*if(pos_counter>=2097152)
                return true;*/

            cell_list[pos_counter] = (cell_basic *)cell_oc_pos;
            cell_oc_pos->setUniqueID(pos_counter);

            pos_counter++;
            // if(pos_counter % 15000 == 0)
            //     cout << "-> Creating cell list     : " << 100.0 * float(pos_counter) / float(max_cells)
            //          << " [%]        \r" << flush;
        }

        // cout << CLR_LINE;
        // cout << "- Creating cell list                   : done          \n" << flush;
        return true;
    }

    bool findMatchingCell(photon_package * pp)
    {
        Vector3D pos = pp->getPosition();

        if(pos.X() < cell_oc_root->getXmin() || pos.Y() < cell_oc_root->getYmin() ||
           pos.Z() < cell_oc_root->getZmin())
            return false;

        if(pos.X() > cell_oc_root->getXmax() || pos.Y() > cell_oc_root->getYmax() ||
           pos.Z() > cell_oc_root->getZmax())
            return false;

        goNextLevelUp(pp);
        goNextLevelDown(pp);

        return true;
    }

    bool next(photon_package * pp)
    {
        if(!findMatchingCell(pp))
            return false;

        if(!goToNextCellBorder(pp))
            return false;

        return true;
    };

    bool findStartingPoint(photon_package * pp);

    void getLengths(uint bins, double & step_xy, double & off_xy)
    {
        step_xy = (cell_oc_root->getLength()) / double(bins);

        off_xy = step_xy / 2.0;
    }

    double getVolume(const cell_basic & cell) const
    {
        const cell_oc * cell_pos = (const cell_oc *)&cell;

        double volume = cell_pos->getLength();
        volume = volume * volume * volume;

        return volume;
    }

    bool positionPhotonInGrid(photon_package * pp)
    {
        pp->setPositionCell(cell_oc_root);
        return findMatchingCell(pp);
    }

    const cell_oc * getTopLevelCell() const
    {
        return cell_oc_root;
    }

    const cell_oc * getCurrentCell() const
    {
        return cell_oc_pos;
    }

    double getMaxLength()
    {
        return cell_oc_root->getLength();
    }

    void printParameters(parameters & param);

    bool createArtificialGrid(string path);

    void createNextLevel(cell_oc * cell);

    bool saveBinaryGridFile(string filename)
    {
        return saveBinaryGridFile(filename, GRID_ID_OCT, data_offset);
    };

    bool loadGridFromBinaryFile(parameters & param, uint data_len);

    bool saveBinaryGridFile(string filename, ushort id, ushort data_size);

    bool loadGridFromBinaryFile(parameters & param)
    {
        return loadGridFromBinaryFile(param, 0);
    };

    void nextBinaryDataCell(ofstream & file_stream, cell_oc * cell, uint data_size);

    void clear()
    {
        line_counter = 0;
        char_counter = 0;
        clear(cell_oc_root);
        cell_oc_root = 0;
        cell_oc_pos = 0;
        cout << "Final cleanup                                : done" << endl;
    }

    void goNextLevelDown(photon_package * pp)
    {
        cell_oc * tmp_cell = (cell_oc *)pp->getPositionCell();

        if(tmp_cell->getChildren() != 0)
        {
            Vector3D pos = pp->getPosition();
            Vector3D center = getCenter(*pp->getPositionCell());
            double x_mid = center.X();
            double y_mid = center.Y();
            double z_mid = center.Z();

            if(pos.Z() < z_mid) // z 0 1 2 3
            {
                if(pos.Y() < y_mid) // y 0 1
                {
                    if(pos.X() < x_mid) // x 0
                        tmp_cell = tmp_cell->getChild(0);
                    else
                        // x 1
                        tmp_cell = tmp_cell->getChild(1);
                }
                else // y 2 3
                {
                    if(pos.X() < x_mid) // x 2
                        tmp_cell = tmp_cell->getChild(2);
                    else // x 3
                        tmp_cell = tmp_cell->getChild(3);
                }
            }
            else // z 4 5 6 7
            {
                if(pos.Y() < y_mid) // y 4 5
                {
                    if(pos.X() < x_mid) // x 4
                        tmp_cell = tmp_cell->getChild(4);
                    else // x 5
                        tmp_cell = tmp_cell->getChild(5);
                }
                else // y 6 7
                {
                    if(pos.X() < x_mid) // x 6
                        tmp_cell = tmp_cell->getChild(6);
                    else // x 7
                        tmp_cell = tmp_cell->getChild(7);
                }
            }

            pp->setPositionCell(tmp_cell);
            goNextLevelDown(pp);
        }
    }

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

    void goNextLevelUp(photon_package * pp)
    {
        cell_oc * tmp_cell = (cell_oc *)pp->getPositionCell();

        if(tmp_cell == 0)
            return;

        if(tmp_cell->getLevel() == 0)
            return;

        tmp_cell = tmp_cell->getParent();

        pp->setPositionCell(tmp_cell);

        if(!isInside(pp->getPosition(), *tmp_cell))
            goNextLevelUp(pp);
    }

    bool isInside(const Vector3D & pos, const cell_basic & _cell) const
    {
        const cell_oc * tmp_cell = (const cell_oc *)&_cell;

        if(tmp_cell->getXmin() > pos.X())
            return false;

        if(tmp_cell->getYmin() > pos.Y())
            return false;

        if(tmp_cell->getZmin() > pos.Z())
            return false;

        if(tmp_cell->getXmax() < pos.X())
            return false;

        if(tmp_cell->getYmax() < pos.Y())
            return false;

        if(tmp_cell->getZmax() < pos.Z())
            return false;

        return true;
    }

    bool isInside(const Vector3D & pos) const
    {
        if(pos.X() < cell_oc_root->getXmin() || pos.Y() < cell_oc_root->getYmin() ||
           pos.Z() < cell_oc_root->getZmin())
            return false;

        if(pos.X() > cell_oc_root->getXmax() || pos.Y() > cell_oc_root->getYmax() ||
           pos.Z() > cell_oc_root->getZmax())
            return false;

        return true;
    }

    // bool isInside(photon_package * pp, Vector3D & pos)
    // {
    //     cell_oc * cell = (cell_oc *)pp->getPositionCell();
    //     if(pos.X() < cell->getXmin() || pos.Y() < cell->getYmin() || pos.Z() < cell->getZmin())
    //         return false;

    //     if(pos.X() > cell->getXmax() || pos.Y() > cell->getYmax() || pos.Z() > cell->getZmax())
    //         return false;

    //     return true;
    // }

    void setRndPositionInCell(photon_package * pp, CRandomGenerator * rand_gen)
    {
        Vector3D pos;
        cell_oc * tmp_cell = (cell_oc *)pp->getPositionCell();
        double x, dx, y, dy, z, dz;

        x = tmp_cell->getXmin();
        dx = tmp_cell->getXmax() - x;
        y = tmp_cell->getYmin();
        dy = tmp_cell->getYmax() - y;
        z = tmp_cell->getZmin();
        dz = tmp_cell->getZmax() - z;

        double rnd_x = rand_gen->getRND();
        double rnd_y = rand_gen->getRND();
        double rnd_z = rand_gen->getRND();

        pos = Vector3D(x + rnd_x * dx, y + rnd_y * dy, z + rnd_z * dz);
        pp->setPosition(pos);
    }
};
