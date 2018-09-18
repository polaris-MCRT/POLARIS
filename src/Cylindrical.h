#include "typedefs.h"
#include "Vector.h"
#include "chelper.h"
#include "MathFunctions.h"
#include "Grid.h"
#include "Matrix2D.h"
#include "Source.h"

class CGridCylindrical: public CGridBasic
{
public:

    CGridCylindrical(void)
    {
        basic_path = 0;
        buffer_size = 0;

        max_cells = 0;
        max_value = 0;
        max_data = 0;

        min_delta = 1e300;
        max_delta = 0;

        min_mach = 1e300;
        max_mach = 0;

        min_mag = 1e300;
        max_mag = 0;

        min_vel = 1e300;
        max_vel = 0;

        min_len = 1e300;
        max_len = 0;

        min_gas_temp = 1e300;
        max_gas_temp = 0;

        min_dust_temp = 1e300;
        max_dust_temp = 0;

        min_gas_dens = 1e300;
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

        nrOfGnuPoints = 1000;
        nrOfGnuVectors = 1000;
        maxGridLines = 3;

        cell_list = 0;

        data_offset = 6;
        dataID = 0;
        data_len = 0;

        total_gas_mass = 0;
        mu = 2.0;

        nrOfDensRatios = 0;
        nrOfOpiateIDs = 0;

        //data_pos_gd = MAX_UINT;
        //data_pos_dd = MAX_UINT;
        //data_pos_td = MAX_UINT;
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
        data_pos_aalg = MAX_UINT;
        data_pos_amin = MAX_UINT;
        data_pos_amax = MAX_UINT;
        data_pos_eq = MAX_UINT;
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

        plt_gas_dens=false;
        plt_dust_dens=false;
        plt_gas_temp=false;
        plt_dust_temp=false;
        plt_mag=false;
        plt_vel=false;
        plt_rat=false;
        plt_delta=false;
        plt_larm=false;
        plt_mach=false;
        plt_dust_id=false;

        total_volume = 0;
        cell_volume = 0;

        Rmin = 0;
        Rmax = 1;
        Zmax = 1;
        N_r = 0;
        N_ph = 0;
        N_z = 0;
        log_factorR = 0;
        log_factorPh = 0;
        log_factorZ = 0;

        grid_cells = 0;
        center_cells = 0;

        cell_list = 0;

        listR = 0;
        listPh = 0;
        listZ = 0;

        line_counter = 0;
        char_counter = 0;

        rot_angle1 = 0;
        rot_angle2 = 0;

        turbulent_velocity = 0;
    }

    ~CGridCylindrical()
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

        if(grid_cells != 0)
        {
            for(uint i_r = 0; i_r < N_r; i_r++)
            {
                cout << "Cleaning memory for cylindrical grid file : " <<
                        float(100.0 * double(i_r) / double(N_r)) << "      \r" << flush;

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
        }

        if(center_cells != 0)
        {
            for(uint i_z = 0; i_z < N_z; i_z++)
            {
                delete center_cells[i_z];
                center_cells[i_z] = 0;
            }
        }

        cout << CLR_LINE << flush;
    }

    bool writeGNUPlotFiles(string path, parameter & param);

    bool goToNextCellBorder(photon_package * pp);
    bool updateShortestDistance(photon_package * pp);

    Vector3D getCenter(cell_basic * cell)
    {
        Vector3D center;
        cell_cyl * tmp_cell = (cell_cyl *) cell;

        if(tmp_cell->getRID() == MAX_UINT)
        {
            double z = listZ[0][tmp_cell->getZID()];
            double dz = listZ[0][tmp_cell->getZID() + 1] - z;
            return Vector3D(0, 0, z + 0.5 * dz);
        }

        double r = listR[tmp_cell->getRID()];
        double dr = listR[tmp_cell->getRID() + 1] - r;
        double ph = listPh[tmp_cell->getRID()][tmp_cell->getPhID()];
        double dph = listPh[tmp_cell->getRID()][tmp_cell->getPhID() + 1] - ph;
        double z = listZ[tmp_cell->getRID()][tmp_cell->getZID()];
        double dz = listZ[tmp_cell->getRID()][tmp_cell->getZID() + 1] - z;

        double sin_ph = sin(ph + 0.5 * dph);
        double cos_ph = cos(ph + 0.5 * dph);

        center = (r + 0.5 * dr) * Vector3D(cos_ph, sin_ph, 0);
        center.setZ(z + 0.5 * dz);

        return center;
    }

    double maxLength()
    {
        return max(Rmax, Zmax);
    }

    bool next(photon_package * pp)
    {
        if(!positionPhotonInGrid(pp))
            return false;

        if(!goToNextCellBorder(pp))
            return false;

        return true;
    };

    /*void getBoundingPoints(Vector3D & p_min, Vector3D & p_max)
    {
        p_min.set(cell_oc_root->x_min, cell_oc_root->y_min,
                cell_oc_root->z_min);
        p_max.set(cell_oc_root->x_max, cell_oc_root->y_max,
                cell_oc_root->z_max);
    }

    void getBoundingPoints(cell_basic * cell, Vector3D & p_min,
            Vector3D & p_max)
    {
        cell_oc * curr_cell = (cell_oc*) cell;
        p_min.set(curr_cell->x_min, curr_cell->y_min, curr_cell->z_min);
        p_max.set(curr_cell->x_max, curr_cell->y_max, curr_cell->z_max);
    }*/

    bool findStartingPoint(photon_package * pp);

    void getLengths(uint bins, double & step_xy, double & off_xy)
    {
        step_xy = 2 * Rmax / double(bins);
        off_xy = step_xy / 2.0;
    }

    bool createCellList()
    {
        if(max_cells == 0)
        {
            cout << "ERROR: Cylindrical grid contains no cells!" << endl;
            cout << "       Cell list cannot be created!" << endl;
            return false;
        }

        cell_list = new cell_basic *[max_cells];
        ulong pos_counter = 0;

        cout << CLR_LINE;
        cout << "-> Creating cell list    : 0 %           \r" << flush;

        for(uint i_r = 0; i_r < N_r; i_r++)
        {
            cout << "-> Creating cell list     : "
                    << 100.0 * float(i_r) / float(N_r)
                    << " %        \r" << flush;

            for(uint i_ph = 0; i_ph < N_ph[i_r]; i_ph++)
            {
                for(uint i_z = 0; i_z < N_z; i_z++)
                {
                    cell_list[pos_counter] = (cell_basic*) grid_cells[i_r][i_ph][i_z];
                    pos_counter++;
                }
            }
        }

        for(uint i_z = 0; i_z < N_z; i_z++)
        {
            cell_list[pos_counter] = (cell_basic*) center_cells[i_z];
            pos_counter++;
        }

        cout << CLR_LINE;
        //cout << "- Creating of cell list         : done          \n" << flush;
        return true;
    }

    double getVolume(cell_basic * cell)
    {
        cell_cyl * cell_pos = (cell_cyl*) cell;

        double volume = 0;

        if(cell_pos->getRID() == MAX_UINT)
        {
            double z1 = listZ[0][cell_pos->getZID()];
            double z2 = listZ[0][cell_pos->getZID() + 1];
            volume = PI * Rmin * Rmin * (z2 - z1);
        }
        else
        {
            double r1 = listR[cell_pos->getRID()];
            double r2 = listR[cell_pos->getRID() + 1];
            double ph1 = listPh[cell_pos->getRID()][cell_pos->getPhID()];
            double ph2 = listPh[cell_pos->getRID()][cell_pos->getPhID() + 1];
            double z1 = listZ[cell_pos->getRID()][cell_pos->getZID()];
            double z2 = listZ[cell_pos->getRID()][cell_pos->getZID() + 1];

            volume = 0.5 * (ph2 - ph1)* (r2 * r2 - r1 * r1)* (z2 - z1);
        }

        return volume;
    }

    double getVolume(photon_package * pp)
    {
        cell_basic * cell_pos = pp->getPositionCell();
        return getVolume(cell_pos);
    }

    Vector3D rotateToCenter(photon_package * pp, Vector3D dir, bool inv=false)
    {
        cell_cyl * cell_pos = (cell_cyl*) pp->getPositionCell();
        double phi = pp->getPosition().getPhiCoord();

        double phi_center = 0;
        if(cell_pos->getRID() != MAX_UINT)
            phi_center = 0.5 * (listPh[cell_pos->getRID()][cell_pos->getPhID()] +
                listPh[cell_pos->getRID()][cell_pos->getPhID() + 1]);

        double dph = phi_center - phi;
        if(inv)
            dph *= -1;

        dir.cart2cyl();
        dir.setPhi(dir.Phi() + dph);
        dir.cyl2cart();
        return dir;
    }

    /*double getMinArea(photon_package * pp)
    {
        //tbd
        cell_oc * cell_pos = (cell_oc *) pp->getPositionCell();
        return 0;
    }*/

    bool positionPhotonInGrid(photon_package * pp);

    double getMaxLength()
    {
        return max(Rmax, Zmax);
    }

    bool createArtificialGrid(string path);

    bool saveBinaryGridFile(string filename)
    {
        return saveBinaryGridFile(filename, GRID_ID_CYL, data_offset);
    }

    bool loadGridFromBinrayFile(parameter & param, uint _data_len);
    bool saveBinaryGridFile(string filename, ushort id, ushort data_size);

    bool loadGridFromBinrayFile(parameter & param)
    {
        return loadGridFromBinrayFile(param, 0);
    };

    void clear()
    {
        line_counter = 0;
        char_counter = 0;
        cout << "Final cleanup                                : done" << endl;
    }

    void printParameter();

    bool getPolarRTGridParameter(double max_len, double pixel_width, uint max_subpixel_lvl, 
            dlist &_listR, uint &N_polar_r, uint * &N_polar_ph)
    {
        uint subpixel_multiplier = pow(2, max_subpixel_lvl);

        // Calculate additional rings for the center
        uint N_r_center = uint(ceil(listR[0] / (listR[1] - listR[0])));

        for(uint i_r = 0; i_r <= N_r_center; i_r++)
            _listR.push_back(listR[0] * (i_r / double(N_r_center)));

        for(uint i_r = 1; i_r <= N_r; i_r++)
        {
            double r0 = _listR[_listR.size() - 2];
            double r1 = _listR[_listR.size() - 1];
            double r2 = listR[i_r];
            if((r2 - r1) < 5.0 * (r1 - r0))
                for(int i_subpixel = 1; i_subpixel <= subpixel_multiplier; i_subpixel++)
                    _listR.push_back(r1 + (r2 - r1) * i_subpixel / double(subpixel_multiplier));
            else
            {
                uint N_r_sub = uint(ceil((r2 - r1) / (5.0 * (r1 - r0))));
                for(int i_r_sub = 1; i_r_sub <= N_r_sub; i_r_sub++)
                    for(int i_subpixel = 1; i_subpixel <= subpixel_multiplier; i_subpixel++)
                        _listR.push_back(r1 + (r2 - r1) * i_r_sub / double(N_r_sub) * 
                            i_subpixel / double(subpixel_multiplier));
            }

            // break if sidelength is smaller than full grid
            if(_listR.back() > max_len)
            {
                _listR.pop_back();
                _listR.push_back(max_len);
                break;
            }
        }

        if(_listR.back() <= max_len)
        {
            // Calculate additional rings for the outer rings
            uint N_r_outer = uint(ceil((max_len - listR[N_r]) / (listR[N_r] - listR[N_r - 1])));

            for(uint i_r = 1; i_r <= N_r_outer; i_r++)
                for(int i_subpixel = 1; i_subpixel <= subpixel_multiplier; i_subpixel++)
                    _listR.push_back(listR[N_r] + (max_len - listR[N_r]) * i_r / double(N_r_outer) *
                         i_subpixel / double(subpixel_multiplier));
        }

        // Set total size of the radial cells
        N_polar_r = _listR.size() - 1;

       // Calc the number of phi background grid pixel
        N_polar_ph = new uint[N_polar_r];
        for(uint i_r = 0; i_r < N_polar_r; i_r++)
            N_polar_ph[i_r] = uint(ceil(PIx2 * _listR[i_r + 1] / 
                min(pixel_width, (_listR[i_r + 1] - _listR[i_r]))));

        return true;
    }

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

    bool isInside(Vector3D & pos, cell_basic * _cell)
    {
        cell_cyl * cell = (cell_cyl *) _cell;

        Vector3D tmp_pos = pos.getSphericalCoord();

        if(cell->getRID() == MAX_UINT)
        {
            double sq_r = pos.X() * pos.X() + pos.Y() * pos.Y();
            if(sq_r > Rmin * Rmin)
                return false;

            if(tmp_pos.Z() < -Zmax)
                return false;

            if(tmp_pos.Z() > Zmax)
                return false;
        }

        double r1 = listR[cell->getRID()];
        double r2 = listR[cell->getRID() + 1];
        double ph1 = listPh[cell->getRID()][cell->getPhID()];
        double ph2 = listPh[cell->getRID()][cell->getPhID() + 1];
        double z1 = listZ[cell->getRID()][cell->getZID()];
        double z2 = listZ[cell->getRID()][cell->getZID() + 1];

        if(tmp_pos.R() < r1)
            return false;

        if(tmp_pos.R() > r2)
            return false;

        if(tmp_pos.Phi() < ph1)
            return false;

        if(tmp_pos.Phi() > ph2)
            return false;

        if(tmp_pos.Z() < z1)
            return false;

        if(tmp_pos.Z() > z2)
            return false;

        return true;
    }

    bool isInside(Vector3D & pos)
    {
        double sq_r = pos.X() * pos.X() + pos.Y() * pos.Y();
        if(sq_r > Rmax * Rmax)
            return false;

        if(pos.Z() < -Zmax)
            return false;

        if(pos.Z() > Zmax)
            return false;

        return true;
    }

    bool isInside(photon_package * pp, Vector3D & pos)
    {
        return isInside(pos, pp->getPositionCell());
    }

    void setRndPositionInCell(photon_package * pp)
    {
        Vector3D pos;
        cell_cyl * tmp_cell = (cell_cyl *) pp->getPositionCell();
        double r1, r2, ph1, ph2, z1, z2;

        double rnd_r = pp->getRND();
        double rnd_ph = pp->getRND();
        double rnd_z = pp->getRND();

        if(tmp_cell->getZID() == MAX_UINT)
            return;
        else if(tmp_cell->getRID() == MAX_UINT)
        {
            r1 = 0;
            r2 = listR[0];
            ph1 = listPh[0][0];
            ph2 = listPh[0][N_ph[0]];
            z1 = listZ[0][0];
            z2 = listZ[0][N_z];
        }
        else
        {
            r1 = listR[tmp_cell->getRID()];
            r2 = listR[tmp_cell->getRID() + 1];
            ph1 = listPh[tmp_cell->getRID()][tmp_cell->getPhID()];
            ph2 = listPh[tmp_cell->getRID()][tmp_cell->getPhID() + 1];
            z1 = listZ[tmp_cell->getRID()][tmp_cell->getZID()];
            z2 = listZ[tmp_cell->getRID()][tmp_cell->getZID() + 1];
        }

        double sin_ph = sin(ph1 + rnd_ph * (ph2 - ph1));
        double cos_ph = cos(ph1 + rnd_ph * (ph2 - ph1));

        pos = (pow(pow(r1, 2) + rnd_r * (pow(r2, 2) - pow(r1, 2)), 1.0 / 2.0)
                * Vector3D(cos_ph, sin_ph, 0))
                + Vector3D(0, 0, z1 + rnd_z * (z2 - z1));

        pp->setPosition(pos);
    }
};
