#pragma once
#include "Grid.h"
#include "MathFunctions.h"
#include "Matrix2D.h"
#include "Source.h"
#include "Vector.h"
#include "chelper.h"
#include "typedefs.h"

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

        turbulent_velocity = 0;
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
                cout << "Cleaning memory for spherical grid file: "
                     << float(100.0 * double(i_r) / double(N_r)) << "      \r" << flush;

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

        cout << CLR_LINE << flush;
    }

    bool writeGNUPlotFiles(string path, parameters & param);

    bool goToNextCellBorder(photon_basic * pp);
    bool updateShortestDistance(photon_package * pp);

    Vector3D getCenter(cell_basic * cell)
    {
        Vector3D center;
        cell_sp * tmp_cell = (cell_sp *)cell;

        if(tmp_cell->getRID() == MAX_UINT)
            return center;

        double r = listR[tmp_cell->getRID()];
        double dr = listR[tmp_cell->getRID() + 1] - r;
        double ph = listPh[tmp_cell->getPhID()];
        double dph = listPh[tmp_cell->getPhID() + 1] - ph;
        double th = listTh[tmp_cell->getThID()];
        double dth = listTh[tmp_cell->getThID() + 1] - th;

        double sin_th = sin(th + 0.5 * dth);
        double cos_th = cos(th + 0.5 * dth);
        double sin_ph = sin(ph + 0.5 * dph);
        double cos_ph = cos(ph + 0.5 * dph);

        center = (r + 0.5 * dr) * Vector3D(sin_th * cos_ph, sin_th * sin_ph, cos_th);

        return center;
    }

    void setRndPositionInCell(photon_basic * pp)
    {
        Vector3D pos;
        cell_sp * tmp_cell = (cell_sp *)pp->getPositionCell();
        double r1, r2, ph1, ph2, th1, th2;

        double rnd_r = pp->getRND();
        double rnd_ph = pp->getRND();
        double rnd_th = pp->getRND();

        if(tmp_cell->getRID() == MAX_UINT)
        {
            r1 = 0;
            r2 = listR[0];
            ph1 = listPh[0];
            ph2 = listPh[N_ph];
            th1 = listTh[0];
            th2 = listTh[N_th];
        }
        else
        {
            r1 = listR[tmp_cell->getRID()];
            r2 = listR[tmp_cell->getRID() + 1];
            ph1 = listPh[tmp_cell->getPhID()];
            ph2 = listPh[tmp_cell->getPhID() + 1];
            th1 = listTh[tmp_cell->getThID()];
            th2 = listTh[tmp_cell->getThID() + 1];
        }

        double cos_th = cos(th2) + rnd_th * (cos(th1) - cos(th2));
        double sin_th = sin(acos(cos_th));
        double sin_ph = sin(ph1 + rnd_ph * (ph2 - ph1));
        double cos_ph = cos(ph1 + rnd_ph * (ph2 - ph1));

        pos = pow(pow(r1, 3) + rnd_r * (pow(r2, 3) - pow(r1, 3)), 1.0 / 3.0) *
              Vector3D(sin_th * cos_ph, sin_th * sin_ph, cos_th);

        pp->setPosition(pos);
    }

    double maxLength()
    {
        return Rmax;
    }
    
    bool next(photon_basic * pp)
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

    bool findStartingPoint(photon_basic * pp);

    void getLengths(uint bins, double & step_xy, double & off_xy)
    {
        step_xy = 2 * Rmax / double(bins);
        off_xy = step_xy / 2.0;
    }

    bool createCellList()
    {
        if(max_cells == 0)
        {
            cout << "\nERROR: Spherical grid contains no cells!" << endl;
            cout << "       Cell list cannot be created!" << endl;
            return false;
        }

        cell_list = new cell_basic *[max_cells];
        ulong pos_counter = 0;

        cout << CLR_LINE;
        cout << "-> Creating cell list    : 0.0 [%]           \r" << flush;

        for(uint i_r = 0; i_r < N_r; i_r++)
        {
            cout << "-> Creating cell list     : " << 100.0 * float(i_r) / float(N_r) << " %        \r"
                 << flush;

            for(uint i_ph = 0; i_ph < N_ph; i_ph++)
            {
                for(uint i_th = 0; i_th < N_th; i_th++)
                {
                    cell_list[pos_counter] = (cell_basic *)grid_cells[i_r][i_ph][i_th];
                    pos_counter++;
                }
            }
        }

        cell_list[pos_counter] = (cell_basic *)center_cell;

        // cout << CLR_LINE;
        // cout << "- Creating of cell list                : done          \n" << flush;
        return true;
    }

    double getVolume(cell_basic * cell)
    {
        cell_sp * cell_pos = (cell_sp *)cell;

        if(cell_pos->getRID() == MAX_UINT)
        {
            return 4.0 / 3.0 * PI * Rmin * Rmin * Rmin;
        }

        double r1 = listR[cell_pos->getRID()];
        double r2 = listR[cell_pos->getRID() + 1];
        double ph1 = listPh[cell_pos->getPhID()];
        double ph2 = listPh[cell_pos->getPhID() + 1];
        double th1 = listTh[cell_pos->getThID()];
        double th2 = listTh[cell_pos->getThID() + 1];

        double volume = (ph1 - ph2) * (r1 * r1 * r1 - r2 * r2 * r2) * (cos(th1) - cos(th2)) / 3.0;

        return volume;
    }

    double getVolume(photon_basic * pp)
    {
        cell_basic * cell_pos = pp->getPositionCell();

        return getVolume(cell_pos);
    }

    double getSolidAngle(cell_basic * cell)
    {
        // Returns the solid angle for flux scaling needed for time-dependent RT
        
        cell_sp * cell_pos = (cell_sp *)cell;
        
        double dphi = 0;
        double dtheta = 0;
        
        double solang = 0;

        if(cell_pos->getRID() == MAX_UINT)
        {
            return 4.0 * PI;
        }
        
        double ph1 = listPh[cell_pos->getPhID()];
        double ph2 = listPh[cell_pos->getPhID() + 1];
        double th1 = listTh[cell_pos->getThID()];
        double th2 = listTh[cell_pos->getThID() + 1];
        
        dphi = ph2 - ph1;
        dtheta = th2 - th1;

        solang = dtheta * dphi / (4 * PI);

        return solang;
    }
    
    Vector3D rotateToCenter(photon_basic * pp, Vector3D dir, bool inv, bool phi_only)
    {
        cell_sp * cell_pos = (cell_sp *)pp->getPositionCell();
        Vector3D pos = pp->getPosition().getSphericalCoord();

        double phi_center = cell_pos->getRID() == MAX_UINT
                                ? 0
                                : 0.5 * (listPh[cell_pos->getPhID()] + listPh[cell_pos->getPhID() + 1]);
        dir.rot(Vector3D(0, 0, 1), inv ? pos.Phi() - phi_center : phi_center - pos.Phi());

        if(!phi_only)
        {
            double theta_center = cell_pos->getRID() == MAX_UINT
                                      ? PI2
                                      : 0.5 * (listTh[cell_pos->getThID()] + listTh[cell_pos->getThID() + 1]);

            Vector3D n = Vector3D(dir.Y(), -dir.X(), 0);
            n.normalize();
            dir.rot(n, inv ? pos.Theta() - theta_center : theta_center - pos.Theta());
        }

        return dir;
    }

    bool positionPhotonInGrid(photon_basic * pp);

    double getMaxLength()
    {
        return Rmax;
    }

    bool createArtificialGrid(string path);

    bool saveBinaryGridFile(string filename)
    {
        return saveBinaryGridFile(filename, GRID_ID_SPH, data_offset);
    }

    bool loadGridFromBinrayFile(parameters & param, uint _data_len);
    bool saveBinaryGridFile(string filename, ushort id, ushort data_size);

    bool loadGridFromBinrayFile(parameters & param)
    {
        return loadGridFromBinrayFile(param, 0);
    };

    void clear()
    {
        line_counter = 0;
        char_counter = 0;
        cout << "Final cleanup                                : done" << endl;
    }

    void printParameters();

    bool getPolarRTGridParameter(double max_len,
                                 double pixel_width,
                                 uint max_subpixel_lvl,
                                 dlist & _listR,
                                 uint & N_polar_r,
                                 uint *& N_polar_ph)
    {
        uint subpixel_multiplier = pow(2, max_subpixel_lvl);

        // Add same number of cell in the inner region
        uint N_r_inner = min(uint(ceil(subpixel_multiplier * listR[0] / listR[N_r] * N_r)), N_r);
        for(uint i_r = 0; i_r <= N_r_inner; i_r++)
            _listR.push_back(listR[0] * (i_r / double(N_r_inner)));

        for(uint i_r = 1; i_r <= N_r; i_r++)
        {
            double r1 = _listR[_listR.size() - 1];
            double r2 = listR[i_r];
            uint N_r_sub = min(subpixel_multiplier, uint(ceil((r2 - r1) * 5.0 / pixel_width)));

            for(int i_r_sub = 1; i_r_sub <= N_r_sub; i_r_sub++)
                _listR.push_back(r1 + (r2 - r1) * i_r_sub / double(N_r_sub));

            // break if sidelength is smaller than full grid
            if(_listR.back() > max_len)
            {
                _listR.pop_back();
                _listR.push_back(max_len);
                break;
            }
        }

        if(_listR.back() < max_len)
        {
            // Calculate additional rings for the outer rings
            uint N_r_outer = uint(ceil((max_len - listR[N_r]) / (listR[N_r] - listR[N_r - 1])));

            for(uint i_r = 1; i_r <= N_r_outer; i_r++)
                _listR.push_back(listR[N_r] + (max_len - listR[N_r]) * i_r / double(N_r_outer));
        }

        // Set total size of the radial cells
        N_polar_r = _listR.size() - 1;

        // Calc the number of phi background grid pixel
        N_polar_ph = new uint[N_polar_r];
        for(uint i_r = 0; i_r < N_polar_r; i_r++)
        {
            N_polar_ph[i_r] = min(uint(subpixel_multiplier * 360.0),
                                  uint(ceil(PIx2 * _listR[i_r + 1] / (_listR[i_r + 1] - _listR[i_r]))));
        }

        return true;
    }

  private:
    double Rmin, Rmax;
    uint N_r, N_ph, N_th;
    double log_factorR, log_factorPh, log_factorTh;

    double * listR;
    double * listPh;
    double * listTh;

    cell_sp **** grid_cells;
    cell_sp * center_cell;

    bool isInside(Vector3D & pos, cell_basic * _cell)
    {
        cell_sp * cell = (cell_sp *)_cell;

        if(cell->getRID() == MAX_UINT)
        {
            if(pos.sq_length() > Rmin * Rmin)
                return false;
        }

        double r1 = listR[cell->getRID()];
        double r2 = listR[cell->getRID() + 1];
        double ph1 = listPh[cell->getPhID()];
        double ph2 = listPh[cell->getPhID() + 1];
        double th1 = listTh[cell->getThID()];
        double th2 = listTh[cell->getThID() + 1];

        Vector3D tmp_pos = pos.getSphericalCoord();

        if(tmp_pos.R() < r1)
            return false;

        if(tmp_pos.R() > r2)
            return false;

        if(tmp_pos.Phi() < ph1)
            return false;

        if(tmp_pos.Phi() > ph2)
            return false;

        if(tmp_pos.Theta() < th1)
            return false;

        if(tmp_pos.Theta() > th2)
            return false;

        return true;
    }

    bool isInside(Vector3D & pos)
    {
        if(pos.sq_length() > Rmax * Rmax)
            return false;

        return true;
    }

    bool isInside(photon_basic * pp, Vector3D & pos)
    {
        return isInside(pos, pp->getPositionCell());
    }
};
