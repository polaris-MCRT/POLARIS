#pragma once
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include <string>

#include "Parameters.h"
#include "Photon.h"
#include "Vector.h"
#include "Cell.h"
#include "MathFunctions.h"
#include "Stokes.h"
#include "Typedefs.h"

#ifndef CGRIDBASIC
#define CGRIDBASIC

// Additional Structure
struct VelFieldInterp
{
    spline vel_field;
    bool zero_vel_field;
    Vector3D start_pos;
};

// Additional Structures
class MagFieldInfo
{
    public:
        MagFieldInfo()
        {
            cos_theta=0;
            sin_theta=0;
            cos_2_phi=0;
            sin_2_phi=0;
        }


        double cos_theta;
        double sin_theta;
        double cos_2_phi;
        double sin_2_phi;
        Vector3D mag_field;
};

class LineBroadening
{
    public:
    LineBroadening()
    {
        gauss_a=0;
        voigt_a=0;
    }

    double gauss_a;
    double voigt_a;
};

class CGridBasic
{
  public:
    CGridBasic()
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

        aalg_min = 0;
        aalg_max = 0;

        a_min_min = 0;
        a_min_max = 0;

        a_max_min = 0;
        a_max_max = 0;

        size_param_min = 0;
        size_param_max = 0;

        dust_id_min = 0;
        dust_id_max = 0;

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

        dust_is_mass_density = false;
        gas_is_mass_density = false;
        velocity_field_needed = false;
        spec_length_as_vector = false;

        nrOfPlotPoints = 1000;
        nrOfPlotVectors = 1000;
        maxPlotLines = 3;

        cell_list = 0;

        data_offset = 6;
        dataID = 0;
        data_len = 0;

        total_gas_mass = 0;
        mu = 0;

        nrOfDensRatios = 0;
        nrOfOpiateIDs = 0;

        nr_densities = 1;
        size_gd_list = 0;
        size_dd_list = 0;
        multi_temperature_entries = 0;
        stochastic_temperature_entries = 0;
        nr_mixtures = 0;

        nr_dust_temp_sizes = 0;
        nr_stochastic_sizes = 0;
        nr_stochastic_temps = 0;
        size_skip = 0;

        level_to_pos = 0;
        line_to_pos = 0;

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
        data_pos_amin = MAX_UINT;
        data_pos_amax = MAX_UINT;
        data_pos_size_param = MAX_UINT;
        data_pos_ra = MAX_UINT;
        data_pos_id = MAX_UINT;

        data_pos_vt = MAX_UINT;
        data_pos_pda = MAX_UINT;

        pos_GasSpecRatios = 0;
        pos_OpiateIDS = 0;

        rot_angle1 = 0;
        rot_angle2 = 0;

        nr_rad_field_comp = 1;

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
        plt_amin = false;
        plt_amax = false;
        plt_size_param = false;
        plt_rad_field1 = false;
        plt_u_rad = false;
        plt_g_zero1 = false;
        plt_n_th = false;
        plt_T_e = false;
        plt_n_cr = false;
        plt_g_min = false;
        plt_g_max = false;
        plt_p = false;

        plt_avg_dir = false;
        plt_avg_th = false;

        total_volume = 0;
        cell_volume = 0;

        buffer_gas_dens = 0;
        buffer_mol_dens = 0;
        buffer_dust_dens = 0;
        buffer_gas_temp = 0;
        buffer_dust_temp = 0;
        buffer_rat = 0;
        buffer_delta = 0;
        buffer_mag = 0;
        buffer_mag_x = 0;
        buffer_mag_y = 0;
        buffer_mag_z = 0;
        buffer_vel = 0;
        buffer_vel_x = 0;
        buffer_vel_y = 0;
        buffer_vel_z = 0;
        buffer_larm = 0;
        buffer_mach = 0;
        buffer_dust_mixture = 0;
        buffer_dust_amin = 0;
        buffer_dust_amax = 0;
        buffer_dust_size_param = 0;
        buffer_rad_field = 0;
        buffer_g_zero1 = 0;
        buffer_u_rad = 0;
        buffer_n_th = 0;
        buffer_T_e = 0;
        buffer_n_cr = 0;
        buffer_g_min = 0;
        buffer_g_max = 0;
        buffer_p = 0;

        buffer_avg_dir = 0;
        buffer_avg_th = 0;

        turbulent_velocity = 0;

        wl_list.resize(WL_STEPS);
        CMathFunctions::LogList(WL_MIN, WL_MAX, wl_list, 10);

        CextMeanTab = 0;
        CabsMeanTab = 0;
        CscaMeanTab = 0;
        numberDensityTab = 0;
        totalCellEmissionTab = 0;
        chosen_wID = 0;
        max_wavelengths = 0;
    }

    virtual ~CGridBasic(void)
    {
        if(cell_list != 0)
        {
            delete[] cell_list;
            cell_list = 0;
        }

        if(pos_GasSpecRatios != 0)
        {
            delete[] pos_GasSpecRatios;
            pos_GasSpecRatios = 0;
        }

        if(pos_OpiateIDS != 0)
        {
            delete[] pos_OpiateIDS;
            pos_OpiateIDS = 0;
        }

        if(nr_dust_temp_sizes != 0)
            delete[] nr_dust_temp_sizes;

        if(nr_stochastic_sizes != 0)
            delete[] nr_stochastic_sizes;

        if(nr_stochastic_temps != 0)
            delete[] nr_stochastic_temps;

        if(size_skip != 0)
            delete[] size_skip;

        if(CextMeanTab != 0)
        {
            for(uint wID = 0; wID < max_wavelengths; wID++)
                delete[] CextMeanTab[wID];
            delete[] CextMeanTab;
            CextMeanTab = 0;
        }

        if(CabsMeanTab != 0)
        {
            for(uint wID = 0; wID < max_wavelengths; wID++)
                delete[] CabsMeanTab[wID];
            delete[] CabsMeanTab;
            CabsMeanTab = 0;
        }

        if(CscaMeanTab != 0)
        {
            for(uint wID = 0; wID < max_wavelengths; wID++)
                delete[] CscaMeanTab[wID];
            delete[] CscaMeanTab;
            CscaMeanTab = 0;
        }

        if(numberDensityTab != 0)
        {
            delete[] numberDensityTab;
            numberDensityTab = 0;
        }

        if(totalCellEmissionTab != 0)
        {
            delete[] totalCellEmissionTab;
            totalCellEmissionTab = 0;
        }
    }

    void printPhysicalParameters();

    void resetGridValues()
    {
        max_cells = 0;
        max_value = 0;
        max_data = 0;

        delta0 = 8.28e23 * 2.5e-12 * 1e8 * 1e-6 * 1e6;
        larm_f = 4.1e-21;

        max_gas_dens = -1e300;
        min_gas_dens = 1e300;

        max_dust_dens = -1e300;
        min_dust_dens = 1e300;

        max_gas_temp = -1e300;
        min_gas_temp = 1e300;

        max_dust_temp = -1e300;
        min_dust_temp = 1e300;

        max_larm_limit = -1e300;
        min_larm_limit = 1e300;

        max_delta = -1e300;
        min_delta = 1e300;

        max_mach = -1e300;
        min_mach = 1e300;

        aalg_max = -1e300;
        aalg_min = 1e300;

        a_min_min = 1e300;
        a_min_max = -1e300;

        a_max_min = 1e300;
        a_max_max = -1e300;

        size_param_min = 1e300;
        size_param_max = -1e300;

        dust_id_min = MAX_UINT;
        dust_id_max = 0;

        max_pres = -1e300;
        min_pres = 1e300;

        max_vel = -1e300;
        min_vel = 1e300;

        max_mag = -1e300;
        min_mag = 1e300;

        max_value = 0;

        max_len = -1e300;
        min_len = 1e300;

        min_n_th = 1e300;
        min_T_e = 1e300;
        min_n_cr = 1e300;
        min_g_min = 1e300;
        min_g_max = 1e300;
        min_p = 1e300;

        max_n_th = -1e300;
        max_T_e = -1e300;
        max_n_cr = -1e300;
        max_g_min = -1e300;
        max_g_max = -1e300;
        max_p = -1e300;

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
        data_pos_amin = MAX_UINT;
        data_pos_amax = MAX_UINT;
        data_pos_size_param = MAX_UINT;
        data_pos_ra = MAX_UINT;
        data_pos_id = MAX_UINT;

        data_pos_vt = MAX_UINT;
        data_pos_pda = MAX_UINT;

        data_pos_op = MAX_UINT;

        data_pos_n_th = MAX_UINT;
        data_pos_T_e = MAX_UINT;
        data_pos_n_cr = MAX_UINT;
        data_pos_g_min = MAX_UINT;
        data_pos_g_max = MAX_UINT;
        data_pos_p = MAX_UINT;

        data_pos_avg_th = MAX_UINT;
        data_pos_avg_dir = MAX_UINT;

        nr_rad_field_comp = 1;

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
        plt_amin = false;
        plt_amax = false;
        plt_size_param = false;
        plt_rad_field1 = false;
        plt_u_rad = false;
        plt_g_zero1 = false;
        plt_n_th = false;
        plt_T_e = false;
        plt_n_cr = false;
        plt_g_min = false;
        plt_g_max = false;
        plt_p = false;

        plt_avg_dir = false;
        plt_avg_th = false;

        total_volume = 0;
        cell_volume = 0;

        buffer_gas_dens = 0;
        buffer_mol_dens = 0;
        buffer_dust_dens = 0;
        buffer_gas_temp = 0;
        buffer_dust_temp = 0;
        buffer_rat = 0;
        buffer_delta = 0;
        buffer_mag = 0;
        buffer_mag_x = 0;
        buffer_mag_y = 0;
        buffer_mag_z = 0;
        buffer_vel = 0;
        buffer_vel_x = 0;
        buffer_vel_y = 0;
        buffer_vel_z = 0;
        buffer_larm = 0;
        buffer_mach = 0;
        buffer_dust_mixture = 0;
        buffer_dust_amin = 0;
        buffer_dust_amax = 0;
        buffer_dust_size_param = 0;
        buffer_rad_field = 0;
        buffer_g_zero1 = 0;
        buffer_u_rad = 0;
        buffer_n_th = 0;
        buffer_T_e = 0;
        buffer_n_cr = 0;
        buffer_g_min = 0;
        buffer_g_max = 0;
        buffer_p = 0;

        buffer_avg_dir = 0;
        buffer_avg_th = 0;

        turbulent_velocity = 0;

        CextMeanTab = 0;
        CabsMeanTab = 0;
        CscaMeanTab = 0;
        numberDensityTab = 0;
        totalCellEmissionTab = 0;
        chosen_wID = 0;
        max_wavelengths = 0;
    }

    double getCextMeanTab(uint cellID) const;
    double getCabsMeanTab(uint cellID) const;
    double getCscaMeanTab(uint cellID) const;
    double getNumberDensityTab(uint cellID) const;
    double getTotalCellEmissionTab(uint cellID) const;

    void setCextMeanTab(double Cext, uint cellID);
    void setCabsMeanTab(double Cext, uint cellID);
    void setCscaMeanTab(double Cext, uint cellID);
    void setNumberDensityTab(double density, uint cellID);
    void setTotalCellEmissionTab(double cell_emission, uint cellID);

    void setWaveID(uint wID);
    void initPreCalcTables(uint nr_used_wavelengths);

    double getTurbulentVelocity(cell_basic * cell)
    {
        if(turbulent_velocity > 0)
            return turbulent_velocity;
        else if(hasTurbulentVelocity())
            return cell->getData(data_pos_vt);
        return 0;
    }

    double getTurbulentVelocity(photon_package * pp)
    {
        return getTurbulentVelocity(pp->getPositionCell());
    }

    void updateDataRange(cell_basic * cell);

    void updateVelocity(cell_basic * cell, parameters & param)
    {
        if(param.getIsSpeedOfSound() && data_pos_tg != MAX_UINT)
        {
            double speed_of_sound = sqrt((con_kB * cell->getData(data_pos_tg)) / (mu * m_H));
            double vx_tmp = cell->getData(data_pos_vx);
            cell->setData(data_pos_vx, vx_tmp * speed_of_sound);
            double vy_tmp = cell->getData(data_pos_vy);
            cell->setData(data_pos_vy, vy_tmp * speed_of_sound);
            double vz_tmp = cell->getData(data_pos_vz);
            cell->setData(data_pos_vz, vz_tmp * speed_of_sound);
        }
    }

    uint getDataOffset()
    {
        return data_offset;
    }

    uint getDataID()
    {
        return dataID;
    }

    bool hasVelocityField()
    {
        return (data_pos_vx != MAX_UINT && data_pos_vy != MAX_UINT && data_pos_vz != MAX_UINT);
    }

    bool hasTurbulentVelocity()
    {
        return (data_pos_vt != MAX_UINT);
    }

    virtual Vector3D getCenter(const cell_basic & cell) const = 0;

    Vector3D getCenter(const photon_package & pp) const
    {
        return getCenter(*pp.getPositionCell());
    }

    uint getDataLength()
    {
        return data_len;
    }

    virtual double maxLength() = 0;

    ulong getMaxDataCells()
    {
        return max_cells;
    }

    uint getDataSize()
    {
        return max_data;
    }

    virtual void printParameters() = 0;

    virtual bool goToNextCellBorder(photon_package * pp) = 0;

    virtual bool updateShortestDistance(photon_package * pp)
    {
        return false;
    }

    void setOrientation(Vector3D n1, Vector3D n2, double _rot_angle1, double _rot_angle2)
    {
        rot_angle1 = _rot_angle1;
        rot_angle2 = _rot_angle2;

        ex.set(1, 0, 0);
        ey.set(0, 1, 0);
        ez.set(0, 0, 1);

        double cos_a = cos(rot_angle1);
        double sin_a = sin(_rot_angle1);

        ex.rot(n1, cos_a, sin_a);
        ey.rot(n1, cos_a, sin_a);
        ez.rot(n1, cos_a, sin_a);

        cos_a = cos(rot_angle2);
        sin_a = sin(rot_angle2);

        ex.rot(n2, cos_a, sin_a);
        ey.rot(n2, cos_a, sin_a);
        ez.rot(n2, cos_a, sin_a);

        ex.normalize();
        ey.normalize();
        ez.normalize();

        cout << "grid: " << ex << ey << ez << endl;
    }

    void getMagFieldInfo(const photon_package & pp, MagFieldInfo * mfo) const
    {
        // Get the magnetic field from grid
        mfo->mag_field = getMagField(pp);

        // Get the theta and phi angle from the magnetic field direction
        double theta = getThetaMag(pp);
        double phi = getPhiMag(pp);

        // Calculate the sine and cosine including double angles
        mfo->cos_theta = cos(theta);
        mfo->sin_theta = sin(theta);
        mfo->cos_2_phi = cos(2.0 * phi);
        mfo->sin_2_phi = sin(2.0 * phi);
    }

    virtual bool findStartingPoint(photon_package * pp) = 0;

    virtual void getLengths(uint bins, double & step_xy, double & off_xy) = 0;

    virtual bool positionPhotonInGrid(photon_package * pp) = 0;

    double getMinLength()
    {
        return min_len;
    }

    double getMaxLength()
    {
        return max_len;
    }

    virtual bool getPolarRTGridParameter(double max_len,
                                         double pixel_width,
                                         uint max_subpixel_lvl,
                                         dlist & _listR,
                                         uint & N_polar_r,
                                         uint *& N_polar_ph)
    {
        return false;
    }

    virtual bool createCellList() = 0;

    cell_basic * getCellFromIndex(ulong i)
    {
        return cell_list[i];
    }

    void setSIConversionFactors(parameters & param)
    {
        mu = param.getMu();
        conv_length_in_SI = param.getSIConvLength();

        delta0 = param.getDelta0();
        larm_f = param.getLarmF();

        conv_dens_in_SI = abs(param.getSIConvDH());
        conv_Bfield_in_SI = param.getSIConvBField();
        conv_Vfield_in_SI = param.getSIConvVField();
    }

    virtual bool writePlotFiles(string path, parameters & param) = 0;

    void setDataSize(uint sz)
    {
        max_data = data_offset + sz;
    }

    void setDustInformation(uint _nr_mixtures,
                            uint * _nr_dust_temp_sizes,
                            uint * _nr_stochastic_sizes,
                            uint * _nr_stochastic_temps)
    {
        nr_mixtures = _nr_mixtures;
        nr_dust_temp_sizes = _nr_dust_temp_sizes;
        nr_stochastic_sizes = _nr_stochastic_sizes;
        nr_stochastic_temps = _nr_stochastic_temps;
    }

    void setGasInformation(uint ** _level_to_pos, uint *** _line_to_pos)
    {
        level_to_pos = _level_to_pos;
        line_to_pos = _line_to_pos;
    }

    void setVelocityFieldNeeded(bool val)
    {
        velocity_field_needed = val;
    }

    void setDataOffset(uint off)
    {
        data_offset = off;
    }

    void setDataID(uint id)
    {
        dataID = id;
    }

    void setGasDensity(photon_package * pp, double dens)
    {
        setGasDensity(pp->getPositionCell(), dens);
    }

    void setGasDensity(photon_package * pp, uint i_density, double dens)
    {
        setGasDensity(pp->getPositionCell(), i_density, dens);
    }

    virtual bool isInside(const Vector3D & pos) const = 0;

    void setSpecLengthAsVector(bool val)
    {
        spec_length_as_vector = val;
    }

    bool specLengthIsVector()
    {
        return spec_length_as_vector;
    }

    void updateSpecLength(photon_package * pp, double len)
    {
        cell_basic * cell = pp->getPositionCell();
        if(spec_length_as_vector)
        {
            uint data_pos = data_offset + 4 * pp->getDustWavelengthID();
            Vector3D e_dir = len * rotateToCenter(*pp, pp->getDirection());
            cell->updateData(data_pos + 0, len);
            cell->updateData(data_pos + 1, e_dir.X());
            cell->updateData(data_pos + 2, e_dir.Y());
            cell->updateData(data_pos + 3, e_dir.Z());
        }
        else
        {
            uint data_pos = data_offset + pp->getDustWavelengthID();
            cell->updateData(data_pos, len);
        }
    }

    void updateSpecLength(cell_basic * cell, uint i_offset, StokesVector stokes) const
    {
        stokes /= getVolume(*cell);
        uint data_pos = data_offset + 4 * i_offset;
        cell->updateData(data_pos + 0, stokes.I());
        cell->updateData(data_pos + 1, stokes.Q());
        cell->updateData(data_pos + 2, stokes.U());
        cell->updateData(data_pos + 3, stokes.V());
    }

    inline double getSpecLength(const cell_basic & cell, uint wID) const
    {
#ifdef CAMPS_BENCHMARK
        // To perform Camps et. al (2015) benchmark.
        double res = 0, wavelength = wl_list[wID], mult = 1e6;
        res = mult * CMathFunctions::mathis_isrf(wavelength);
        // res = 2.99e-14 * CMathFunctions::planck(wl_list[wID], 9000.0);
        return PIx4 * res * getVolume(cell);
#else
        if(spec_length_as_vector)
            return cell.getData(data_offset + 4 * wID + 0);
        else
            return cell.getData(data_offset + wID);
#endif
    }

    inline double getSpecLength(const photon_package & pp, uint wID) const
    {
        return getSpecLength(*pp.getPositionCell(), wID);
    }

    void getSpecLength(const cell_basic & cell, uint wID, double * us, Vector3D * e_dir) const
    {
        uint data_pos = data_offset + 4 * wID;

        *us = cell.getData(data_pos + 0);
        e_dir->setX(cell.getData(data_pos + 1));
        e_dir->setY(cell.getData(data_pos + 2));
        e_dir->setZ(cell.getData(data_pos + 3));
    }

    void saveRadiationField()
    {
#pragma omp parallel for schedule(dynamic)
        for(long c_i = 0; c_i < long(max_cells); c_i++)
        {
            cell_basic * cell = cell_list[c_i];
            double inv_vol = 1 / getVolume(*cell);
            for(uint wID = 0; wID < WL_STEPS; wID++)
            {
                cell->convertData(data_offset + 4 * wID + 0, inv_vol);
                cell->convertData(data_offset + 4 * wID + 1, inv_vol);
                cell->convertData(data_offset + 4 * wID + 2, inv_vol);
                cell->convertData(data_offset + 4 * wID + 3, inv_vol);
            }
        }

        for(uint wID = 0; wID < WL_STEPS; wID++)
        {
            data_ids.push_back(GRIDrad);
            data_ids.push_back(GRIDradx);
            data_ids.push_back(GRIDrady);
            data_ids.push_back(GRIDradz);
        }
        data_offset += 4 * WL_STEPS;
    }

    double getRadiationField(const cell_basic & cell, uint wID) const
    {
#ifdef CAMPS_BENCHMARK
        // To perform Camps et. al (2015) benchmark.
        double res = 0, wavelength = wl_list[wID], mult = 1e6;
        res = mult * CMathFunctions::mathis_isrf(wavelength);
        // res = 2.99e-14 * CMathFunctions::planck(wl_list[wID], 9000.0);
        return PIx4 * res;
#else
        // If the radiation field is needed after temp calculation, use the SpecLength
        // instead
        if(data_pos_rf_list.empty())
            return getSpecLength(cell, wID) / getVolume(cell);
        return cell.getData(data_pos_rf_list[wID]);
#endif
    }

    double getRadiationField(const photon_package & pp, uint wID) const
    {
        return getRadiationField(*pp.getPositionCell(), wID);
    }

    double getRadiationFieldX(const cell_basic & cell, uint wID) const
    {
        // If the radiation field is needed after temp calculation, use the SpecLength
        // instead
        if(data_pos_rx_list.empty())
        {
            if(spec_length_as_vector)
                return cell.getData(data_offset + 4 * wID + 1) / getVolume(cell);
            else
                return 0;
        }
        return cell.getData(data_pos_rx_list[wID]);
    }

    double getRadiationFieldX(const photon_package & pp, uint wID) const
    {
        return getRadiationFieldX(*pp.getPositionCell(), wID);
    }

    double getRadiationFieldY(const cell_basic & cell, uint wID) const
    {
        // If the radiation field is needed after temp calculation, use the SpecLength
        // instead
        if(data_pos_ry_list.empty())
        {
            if(spec_length_as_vector)
                return cell.getData(data_offset + 4 * wID + 2) / getVolume(cell);
            else
                return 0;
        }
        return cell.getData(data_pos_ry_list[wID]);
    }

    double getRadiationFieldY(const photon_package & pp, uint wID) const
    {
        return getRadiationFieldY(*pp.getPositionCell(), wID);
    }

    double getRadiationFieldZ(const cell_basic & cell, uint wID) const
    {
        // If the radiation field is needed after temp calculation, use the SpecLength
        // instead
        if(data_pos_rz_list.empty())
        {
            if(spec_length_as_vector)
                return cell.getData(data_offset + 4 * wID + 3) / getVolume(cell);
            else
                return 0;
        }
        return cell.getData(data_pos_rz_list[wID]);
    }

    double getRadiationFieldZ(const photon_package & pp, uint wID) const
    {
        return getRadiationFieldZ(*pp.getPositionCell(), wID);
    }

    void getRadiationField(const photon_package & pp, uint w, double * us, Vector3D * e_dir) const
    {
        // Init variables and get current cell
        Vector3D tmp_dir;
        const cell_basic & cell = *pp.getPositionCell();

        // Get radiation field strength and direction from cell
        if(data_pos_rf_list.empty())
        {
            if(spec_length_as_vector)
            {
                // Get SpecLength instead if no radiation field in grid
                getSpecLength(cell, w, us, &tmp_dir);
                *us /= getVolume(cell);
            }
            else
                cout << "ERROR: This should not happen" << endl;
        }
        else
        {
            *us = cell.getData(data_pos_rf_list[w]);
            tmp_dir.setX(cell.getData(data_pos_rx_list[w]));
            tmp_dir.setY(cell.getData(data_pos_ry_list[w]));
            tmp_dir.setZ(cell.getData(data_pos_rz_list[w]));
        }

        // Rotate vector from cell center to position
        *e_dir = rotateToCenter(pp, tmp_dir, true);

        // Normalize the radiation field vector
        e_dir->normalize();
    }

    StokesVector getStokesFromRadiationField(const photon_package & pp, uint i_offset) const
    {
        // Init variables
        StokesVector scattering_stokes;
        const cell_basic & cell = *pp.getPositionCell();
        uint data_pos = data_offset + 4 * i_offset;

        scattering_stokes.setI(cell.getData(data_pos + 0));
        scattering_stokes.setQ(cell.getData(data_pos + 1));
        scattering_stokes.setU(cell.getData(data_pos + 2));
        scattering_stokes.setV(cell.getData(data_pos + 3));

        // Rotate vector from cell center to position
        Vector3D rot_dir = rotateToCenter(pp, getCenter(pp), true);

        // Get rotation angle to rotate back into the map/detector frame
        double phi_map =
            getAnglePhi(pp.getEX(), pp.getEY(), rot_dir) - getAnglePhi(pp.getEX(), pp.getEY(), getCenter(pp));

        // Rotate Stokes Vector to be in agreement with the detector plane
        scattering_stokes.rot(phi_map);

        return scattering_stokes;
    }

    void getRadiationFieldInterp(const photon_package & pp,
                                 double wavelength,
                                 double * us,
                                 Vector3D * e_dir) const
    {
        // Do not interpolate if outside of wavelength list
        if(wl_list.back() < wavelength || wl_list.front() > wavelength)
        {
            uint w = 0;
            if(wl_list.back() < wavelength)
                w = wl_list.size() - 1;

            getRadiationField(pp, w, us, e_dir);

            return;
        }

        // Init variables and get current cell
        Vector3D tmp_dir;
        const cell_basic & cell = *pp.getPositionCell();

        // Get wavelength indices from radiation field calculation
        uint wID1 = CMathFunctions::biListIndexSearch(wavelength, wl_list);
        uint wID2 = wID1 + 1;

        // Interpolate radiation field strength and direction
        *us = CMathFunctions::interpolate(wl_list[wID1],
                                          wl_list[wID2],
                                          cell.getData(data_pos_rf_list[wID1]),
                                          cell.getData(data_pos_rf_list[wID2]),
                                          wavelength);
        tmp_dir.setX(CMathFunctions::interpolate(wl_list[wID1],
                                                 wl_list[wID2],
                                                 cell.getData(data_pos_rx_list[wID1]),
                                                 cell.getData(data_pos_rx_list[wID2]),
                                                 wavelength));
        tmp_dir.setY(CMathFunctions::interpolate(wl_list[wID1],
                                                 wl_list[wID2],
                                                 cell.getData(data_pos_ry_list[wID1]),
                                                 cell.getData(data_pos_ry_list[wID2]),
                                                 wavelength));
        tmp_dir.setZ(CMathFunctions::interpolate(wl_list[wID1],
                                                 wl_list[wID2],
                                                 cell.getData(data_pos_rz_list[wID1]),
                                                 cell.getData(data_pos_rz_list[wID2]),
                                                 wavelength));

        // Rotate vector to cell center
        *e_dir = rotateToCenter(pp, tmp_dir, true);

        // Normalize the radiation field vector
        e_dir->normalize();
    }

    double getGZero(const cell_basic & cell) const
    {
        // Init variables
        const double wl1 = 9.1165e-08, wl2 = 2.06640e-07;
        double g_zero = 0;

        // If the radiation field is needed after temp calculation, use the SpecLength
        // instead
        if(spec_length_as_vector)
        {
            for(uint w = 1; w < WL_STEPS; w++)
            {
                double rad_field_1 = getRadiationField(cell, w - 1);
                double rad_field_2 = getRadiationField(cell, w);
                double mult = 0;
                if(wl_list[w] > wl1 && wl_list[w + 1] < wl2)
                    mult = 1;
                else if(wl_list[w] < wl1 && wl_list[w + 1] > wl2)
                    mult = (wl2 - wl1) / (wl_list[w + 1] - wl_list[w]);
                else if(wl_list[w] < wl2 && wl_list[w + 1] > wl2)
                    mult = (wl2 - wl_list[w]) / (wl_list[w + 1] - wl_list[w]);
                else if(wl_list[w] < wl1 && wl_list[w + 1] > wl1)
                    mult = (wl_list[w + 1] - wl1) / (wl_list[w + 1] - wl_list[w]);
                g_zero += mult * ((wl_list[w] - wl_list[w - 1]) * rad_field_1 +
                                  0.5 * (wl_list[w] - wl_list[w - 1]) * (rad_field_2 - rad_field_1));
            }
        }
        else
        {
            for(uint w = 0; w < WL_STEPS; w++)
            {
                if(wl_list[w] > wl1 && wl_list[w + 1] < wl2)
                    g_zero += getSpecLength(cell, w) / getVolume(cell);
                else if(wl_list[w] < wl1 && wl_list[w + 1] > wl2)
                    g_zero += getSpecLength(cell, w) * (wl2 - wl1) / (wl_list[w + 1] - wl_list[w]) /
                              getVolume(cell);
                else if(wl_list[w] < wl2 && wl_list[w + 1] > wl2)
                    g_zero += getSpecLength(cell, w) * (wl2 - wl_list[w]) / (wl_list[w + 1] - wl_list[w]) /
                              getVolume(cell);
                else if(wl_list[w] < wl1 && wl_list[w + 1] > wl1)
                    g_zero += getSpecLength(cell, w) * (wl_list[w + 1] - wl1) /
                              (wl_list[w + 1] - wl_list[w]) / getVolume(cell);
            }
        }

        g_zero /= 1.7836e-06;
        return g_zero;
    }

    double getGZero(const photon_package & pp) const
    {
        return getGZero(*pp.getPositionCell());
    }

    double getUrad(const cell_basic & cell) const
    {
        double u_rad = 0;

        if(spec_length_as_vector)
        {
            for(uint w = 1; w < WL_STEPS; w++)
            {
                double rad_field_1 = getRadiationField(cell, w - 1) / con_c;
                double rad_field_2 = getRadiationField(cell, w) / con_c;

                u_rad += ((wl_list[w] - wl_list[w - 1]) * rad_field_1 +
                          0.5 * (wl_list[w] - wl_list[w - 1]) * (rad_field_2 - rad_field_1));
            }
        }

        return u_rad / (8.64e-14);
    }

    double getUrad(const photon_package & pp) const
    {
        return getUrad(*pp.getPositionCell());
    }

    double getMu() const
    {
        return mu;
    }

    void setRelOutsidePosition(photon_package * pp, double tx, double ty, double tz)
    {
        pp->setPosition(tx * ex + ty * ey + tz * ez);
    }

    void setRelOutsidePosition(photon_package * pp, double tx, double ty)
    {
        pp->setPosition(tx * ex + ty * ey - max_len * ez);
    }

    void setRelDirection(photon_package * pp)
    {
        pp->setEX(ex);
        pp->setEY(ey);
        pp->setEZ(ez);
    }

    virtual void setRndPositionInCell(photon_package * pp)
    {
        pp->setPosition(Vector3D(0, 0, 0));
    }

    Vector3D rotateToCenter(const photon_package & pp, bool inv = false, bool phi_only = false) const
    {
        return rotateToCenter(pp, pp.getDirection(), inv, phi_only);
    }

    virtual Vector3D rotateToCenter(const photon_package & pp,
                                    Vector3D dir,
                                    bool inv = false,
                                    bool phi_only = false) const
    {
        return dir;
    }

    void setDustTemperature(cell_basic * cell, uint i_density, uint a, double temp)
    {
        if(!data_pos_dt_list.empty())
        {
            uint id = a + nr_densities;
            for(uint i = 0; i < i_density; i++)
                id += size_skip[i];
            cell->setData(data_pos_dt_list[id], temp);
        }
    }

    void setDustTemperature(cell_basic * cell, uint i_density, double temp)
    {
        cell->setData(data_pos_dt_list[i_density], temp);
    }

    void setDustTemperature(cell_basic * cell, double temp)
    {
        for(uint i_density = 0; i_density < data_pos_dt_list.size(); i_density++)
            cell->setData(data_pos_dt_list[i_density], temp);
    }

    void setDustTemperature(photon_package * pp, uint i_density, uint a, double temp)
    {
        setDustTemperature(pp->getPositionCell(), i_density, a, temp);
    }

    void setDustTemperature(photon_package * pp, uint i_density, double temp)
    {
        setDustTemperature(pp->getPositionCell(), i_density, temp);
    }

    void setDustTemperature(photon_package * pp, double temp)
    {
        setDustTemperature(pp->getPositionCell(), temp);
    }

    void setDustTempProbability(cell_basic * cell, uint i_density, uint a, uint t, double temp)
    {
        uint id = a * nr_stochastic_temps[i_density] + t;
        for(uint i = 0; i < i_density; i++)
            id += nr_stochastic_sizes[i] * nr_stochastic_temps[i];
        cell->setData(data_offset + id, temp);
    }

    void setPDAValue(cell_basic * cell, double val)
    {
        cell->setData(data_pos_pda, val);
    }

    double getMagMax()
    {
        return max_mag;
    }

    void setGasTemperature(cell_basic * cell, double temp)
    {
        cell->setData(data_pos_tg, temp);
    }

    void setElectronTemperature(cell_basic * cell, double temp)
    {
        if(data_pos_T_e != MAX_UINT)
            cell->setData(data_pos_T_e, temp);
    }

    void setThermalElectronDensity(cell_basic * cell, double dens)
    {
        if(data_pos_n_th != MAX_UINT)
            cell->setData(data_pos_n_th, dens);
    }

    void setCRElectronDensity(cell_basic * cell, double dens)
    {
        if(data_pos_n_cr != MAX_UINT)
            cell->setData(data_pos_n_cr, dens);
    }

    void setGammaMin(cell_basic * cell, double g_min)
    {
        if(data_pos_g_min != MAX_UINT)
            cell->setData(data_pos_g_min, g_min);
    }

    void setGammaMax(cell_basic * cell, double g_max)
    {
        if(data_pos_g_max != MAX_UINT)
            cell->setData(data_pos_g_max, g_max);
    }

    void setPowerLawIndex(cell_basic * cell, double p)
    {
        if(data_pos_p != MAX_UINT)
            cell->setData(data_pos_p, p);
    }

    void setAvgTheta(cell_basic * cell, double phi)
    {
        if(data_pos_avg_th != MAX_UINT)
            cell->setData(data_pos_avg_th, phi);
    }

    void setAvgDir(cell_basic * cell, double dir)
    {
        if(data_pos_avg_dir != MAX_UINT)
            cell->setData(data_pos_avg_dir, dir);
    }

    void setDustChoiceID(cell_basic * cell, uint dust_id)
    {
        if(data_pos_id != MAX_UINT)
            cell->setData(data_pos_id, dust_id);
    }

    void setGasDensity(cell_basic * cell, double dens)
    {
        cell->setData(data_pos_gd_list[0], dens);
    }

    void setGasDensity(cell_basic * cell, uint i_density, double dens)
    {
        cell->setData(data_pos_gd_list[i_density], dens);
    }

    double getQBOffset(const cell_basic & cell, uint i_density) const
    {
        if(data_pos_dt_list.size() == 1)
            return cell.getData(data_pos_dt_list[0]);
        else if(data_pos_dt_list.size() > i_density)
            return cell.getData(data_pos_dt_list[i_density]);
        else
            return 0;
    }

    double getQBOffset(const photon_package & pp, uint i_density) const
    {
        return getQBOffset(*pp.getPositionCell(), i_density);
    }

    double getQBOffset(const cell_basic & cell, uint i_density, uint a) const
    {
        if(!data_pos_dt_list.empty())
        {
            uint id = a + nr_densities;
            for(uint i = 0; i < i_density; i++)
                id += size_skip[i];
            return cell.getData(data_pos_dt_list[id]);
        }
        else
            return 0;
    }

    double getQBOffset(const photon_package & pp, uint i_density, uint a) const
    {
        return getQBOffset(*pp.getPositionCell(), i_density, a);
    }

    void setQBOffset(cell_basic * cell, uint i_density, uint a, double temp)
    {
        if(!data_pos_dt_list.empty())
        {
            uint id = a + nr_densities;
            for(uint i = 0; i < i_density; i++)
                id += size_skip[i];
            cell->setData(data_pos_dt_list[id], temp);
        }
    }

    void setQBOffset(cell_basic * cell, uint i_density, double temp)
    {
        if(data_pos_dt_list.size() == 1)
            cell->setData(data_pos_dt_list[0], temp);
        else if(data_pos_dt_list.size() > 1)
            cell->setData(data_pos_dt_list[i_density], temp);
    }

    uint getNrAlignedRadii()
    {
        return data_pos_aalg_list.size();
    }

    double getAlignedRadius(const cell_basic & cell, uint i_density) const
    {
        if(data_pos_aalg_list.size() == 1)
            return cell.getData(data_pos_aalg_list[0]);
        else if(data_pos_aalg_list.size() > i_density)
            return cell.getData(data_pos_aalg_list[i_density]);
        else
            return 0;
    }

    double getAlignedRadius(const photon_package & pp, uint i_density) const
    {
        return getAlignedRadius(*pp.getPositionCell(), i_density);
    }

    void setAlignedRadius(cell_basic * cell, uint i_density, double _a_alg)
    {
        cell->setData(data_pos_aalg_list[i_density], _a_alg);
    }

    double getMinGrainRadius(const cell_basic & cell) const
    {
        if(data_pos_amin != MAX_UINT)
            return cell.getData(data_pos_amin);
        return 0;
    }

    double getMinGrainRadius(const photon_package & pp) const
    {
        return getMinGrainRadius(*pp.getPositionCell());
    }

    double getMaxGrainRadius(const cell_basic & cell) const
    {
        if(data_pos_amax != MAX_UINT)
            return cell.getData(data_pos_amax);
        return 0;
    }

    double getMaxGrainRadius(const photon_package & pp) const
    {
        return getMaxGrainRadius(*pp.getPositionCell());
    }

    double getGrainSizeParam(const cell_basic & cell) const
    {
        if(data_pos_amax != MAX_UINT)
            return cell.getData(data_pos_size_param);
        return 0;
    }

    double getGrainSizeParam(const photon_package & pp) const
    {
        return getGrainSizeParam(*pp.getPositionCell());
    }

    uint getDustChoiceID(const photon_package & pp) const
    {
        return getDustChoiceID(*pp.getPositionCell());
    }

    uint getDustChoiceID(const cell_basic & cell) const
    {
        if(data_pos_id != MAX_UINT)
            return uint(cell.getData(data_pos_id));
        else
            return 0;
    }

    void getLineBroadening(const photon_package & pp, uint i_trans, LineBroadening * line_broadening) const
    {
        line_broadening->gauss_a = getGaussA(pp);
        line_broadening->voigt_a = getVoigtA(pp, i_trans);
    }

    double getGaussA(const cell_basic & cell) const
    {
        return cell.getData(data_offset);
    }

    double getGaussA(const photon_package & pp) const
    {
        return getGaussA(*pp.getPositionCell());
    }

    double getVoigtA(const cell_basic & cell, uint i_line) const
    {
        return cell.getData(data_offset + 1 + i_line);
    }

    double getVoigtA(const photon_package & pp, uint i_line) const
    {
        return getVoigtA(*pp.getPositionCell(), i_line);
    }

    double getLvlPop(const cell_basic & cell, uint i_lvl, uint i_sublvl = 0) const
    {
        if(level_to_pos[i_lvl][i_sublvl] != MAX_UINT)
            return cell.getData(data_offset + level_to_pos[i_lvl][i_sublvl]);

        return 0;
    }

    double getLvlPop(const photon_package & pp, uint i_lvl, uint i_sublvl = 0) const
    {
        return getLvlPop(*pp.getPositionCell(), i_lvl, i_sublvl);
    }

    void setVelocityField(cell_basic * cell, const Vector3D & vel)
    {
        cell->setData(data_pos_vx, vel.X());
        cell->setData(data_pos_vy, vel.Y());
        cell->setData(data_pos_vz, vel.Z());
    }

    void setVelocityField(photon_package * pp, const Vector3D & vel)
    {
        pp->getPositionCell()->setData(data_pos_vx, vel.X());
        pp->getPositionCell()->setData(data_pos_vy, vel.Y());
        pp->getPositionCell()->setData(data_pos_vz, vel.Z());
    }

    uint getCellID(cell_basic * cell)
    {
        return cell->getUniqueID();
    }

    uint validateDataPositions(parameters & param);
    uint getDataIdsOffset(parameters & param);

    void setLvlPopLower(cell_basic * cell, uint i_line, uint i_sublvl, double lvl_lower)
    {
        cell->setData(data_offset + line_to_pos[i_line][0][i_sublvl], lvl_lower);
    }

    void setLvlPopLower(photon_package * pp, uint i_line, uint i_sublvl, double lvl_lower)
    {
        return setLvlPopLower(pp->getPositionCell(), i_line, i_sublvl, lvl_lower);
    }

    void setLvlPopUpper(cell_basic * cell, uint i_line, uint i_sublvl, double lvl_upper)
    {
        cell->setData(data_offset + line_to_pos[i_line][1][i_sublvl], lvl_upper);
    }

    void setLvlPopUpper(photon_package * pp, uint i_line, uint i_sublvl, double lvl_upper)
    {
        return setLvlPopUpper(pp->getPositionCell(), i_line, i_sublvl, lvl_upper);
    }

    void setLvlPop(cell_basic * cell, uint i_lvl, uint i_sublvl, double lvl_pop)
    {
        if(level_to_pos[i_lvl][i_sublvl] != MAX_UINT)
            cell->setData(data_offset + level_to_pos[i_lvl][i_sublvl], lvl_pop);
    }

    void setLvlPop(photon_package * pp, uint i_lvl, uint i_sublvl, double lvl_pop)
    {
        return setLvlPop(pp->getPositionCell(), i_lvl, i_sublvl, lvl_pop);
    }

    void addToZeemanSublevel(cell_basic * cell, uint i_lvl, uint i_sublvl, double lvl_pop)
    {
        if(level_to_pos[i_lvl][1 + i_sublvl] != MAX_UINT)
            cell->setData(data_offset + level_to_pos[i_lvl][1 + i_sublvl], lvl_pop);
    }

    void addToZeemanSublevel(photon_package * pp, uint i_lvl, uint i_sublvl, double lvl_pop)
    {
        return addToZeemanSublevel(pp->getPositionCell(), i_lvl, i_sublvl, lvl_pop);
    }

    void setLineBroadening(cell_basic * cell, uint i_line_broad, const LineBroadening & line_broadening)
    {
        cell->setData(data_offset + 1 + i_line_broad, line_broadening.voigt_a);
    }

    void setGaussA(cell_basic * cell, double gauss_a)
    {
        // Gauss_a only has to be set once
        cell->setData(data_offset, gauss_a);
    }

    /*void assignOpiateID(cell_basic * cell)
    {
        if(opiate == 0)
            return;

        if(nrOfOpiateIDs != 4)
            return;

        double dx = getPathLength(cell) / conv_length_in_SI;
        double dens = getGasDensity(cell); ///conv_dens_in_SI;
        double temp = getGasTemperature(cell);
        double flih = getOpiateIDParameter(cell, 0);
        double fli2 = getOpiateIDParameter(cell, 1);
        double fluv = getOpiateIDParameter(cell, 2);
        double flge = getOpiateIDParameter(cell, 3);

        dx = log10(dx);
        dx = round(dx * 1000) / 1000;

        dens = log10(dens / 1.67e-24);
        dens = round(dens * 10) / 10;

        temp = log10(temp);

        temp = round(temp * 10) / 10;

        if(temp < 3)
            temp = 3;

        if(temp > 5)
            temp = 5;


        flih = (flih) - 2.0 * dx;
        if(flih < 0)
            flih = -99;

        flih = round(flih * 10) / 10;


        fli2 = (fli2) - 2.0 * dx;
        if(fli2 < 0)
            fli2 = -99;

        fli2 = round(fli2 * 10) / 10;


        fluv = (fluv) - 2.0 * dx;
        if(fluv < 0)
            fluv = -99;

        fluv = round(fluv * 10) / 10;


        flge = log10(flge);

        if(flge < -4)
            flge = -99;

        flge = round(flge * 10) / 10;

        opiate_ids res = opiate->findID(dx, dens, temp, flge, fluv, flih, fli2);
        uint id;

        if(res.uiniqD == MAX_UINT)
            id = MAX_UINT;
        else
            id = res.uiniqD;

        setOpiateID(cell, id);
    };*/

    uint getOpiateID(cell_basic * cell)
    {
        uint pos = pos_OpiateIDS[0];
        return uint(cell->getData(pos));
    }

    uint getOpiateID(const photon_package * pp)
    {
        uint pos = pos_OpiateIDS[0];
        return uint(pp->getPositionCell()->getData(pos));
    }

    uint getOpiateID(photon_package pp)
    {
        uint pos = pos_OpiateIDS[0];
        return uint(pp.getPositionCell()->getData(pos));
    }

    void setOpiateID(cell_basic * cell, uint id)
    {
        uint pos = pos_OpiateIDS[0];
        cell->setData(pos, double(id));
    }

    double getOpiateTestData(photon_package * pp)
    {
        uint pos = pos_OpiateIDS[1];

        return pp->getPositionCell()->getData(pos);
    }

    double getOpiateTestData(cell_basic * cell)
    {
        uint pos = pos_OpiateIDS[1];
        return uint(cell->getData(pos));
    }

    void setOpiateTestData(cell_basic * cell, double val)
    {
        uint pos = pos_OpiateIDS[1];
        cell->setData(pos, val);
    }

    bool fillGridWithOpiateData(uint col_id);

    double getElectronTemperature(const photon_package & pp) const
    {
        return getElectronTemperature(*pp.getPositionCell());
    }

    double getElectronTemperature(const cell_basic & cell) const
    {
        if(data_pos_T_e != MAX_UINT)
            return cell.getData(data_pos_T_e);

        return 0;
    }

    double getThermalElectronDensity(const photon_package & pp) const
    {
        return getThermalElectronDensity(*pp.getPositionCell());
    }

    double getThermalElectronDensity(const cell_basic & cell) const
    {
        if(data_pos_n_th != MAX_UINT)
            return cell.getData(data_pos_n_th);

        return 0;
    }

    double getCRElectronDensity(const photon_package & pp) const
    {
        return getCRElectronDensity(*pp.getPositionCell());
    }

    double getCRElectronDensity(const cell_basic & cell) const
    {
        if(data_pos_n_cr != MAX_UINT)
            return cell.getData(data_pos_n_cr);

        return 0;
    }

    double getGammaMin(const photon_package & pp) const
    {
        return getGammaMin(*pp.getPositionCell());
    }

    double getGammaMin(const cell_basic & cell) const
    {
        if(data_pos_g_min != MAX_UINT)
            return cell.getData(data_pos_g_min);

        return 0;
    }

    double getGammaMax(const photon_package & pp) const
    {
        return getGammaMax(*pp.getPositionCell());
    }

    double getGammaMax(const cell_basic & cell) const
    {
        if(data_pos_g_max != MAX_UINT)
            return cell.getData(data_pos_g_max);

        return 0;
    }

    double getPowerLawIndex(const photon_package & pp) const
    {
        return getPowerLawIndex(*pp.getPositionCell());
    }

    double getPowerLawIndex(const cell_basic & cell) const
    {
        if(data_pos_p != MAX_UINT)
            return cell.getData(data_pos_p);

        return 0;
    }

    double getAvgTheta(const photon_package & pp) const
    {
        return getAvgTheta(*pp.getPositionCell());
    }

    double getAvgTheta(const cell_basic & cell) const
    {
        if(data_pos_avg_th != MAX_UINT)
            return cell.getData(data_pos_avg_th);

        return 0;
    }

    double getAvgDir(const photon_package & pp) const
    {
        return getAvgDir(*pp.getPositionCell());
    }

    double getAvgDir(const cell_basic & cell) const
    {
        if(data_pos_avg_dir != MAX_UINT)
            return cell.getData(data_pos_avg_dir);

        return 0;
    }

    double getDustTemperature(const cell_basic & cell, uint i_density, uint a) const
    {
        if(!data_pos_dt_list.empty())
        {
            uint id = a + nr_densities;
            for(uint i = 0; i < i_density; i++)
                id += size_skip[i];
            return cell.getData(data_pos_dt_list[id]);
        }
        else
            return 0;
    }

    double getDustTemperature(const photon_package & pp, uint i_density, uint a) const
    {
        return getDustTemperature(*pp.getPositionCell(), i_density, a);
    }

    double getDustTemperature(const cell_basic & cell, uint i_density) const
    {
        if(data_pos_dt_list.size() == 1)
            return cell.getData(data_pos_dt_list[0]);
        else if(data_pos_dt_list.size() > i_density)
            return cell.getData(data_pos_dt_list[i_density]);
        else
            return 0;
    }

    double getDustTemperature(const photon_package & pp, uint i_density) const
    {
        return getDustTemperature(*pp.getPositionCell(), i_density);
    }

    double getDustTemperature(const cell_basic & cell) const
    {
        double sum = 0;
        for(uint i_density = 0; i_density < nr_densities; i_density++)
            sum += getDustTemperature(cell, i_density) * getRelativeDustDensity(cell, i_density);
        return sum;
    }

    double getDustTemperature(const photon_package & pp) const
    {
        return getDustTemperature(*pp.getPositionCell());
    }

    double getDustTempProbability(const cell_basic & cell, uint i_density, uint a, uint t) const
    {
        uint id = a * nr_stochastic_temps[i_density] + t;
        for(uint i = 0; i < i_density; i++)
            id += nr_stochastic_sizes[i] * nr_stochastic_temps[i];
        return cell.getData(data_offset + id);
    }

    double getDustTempProbability(const photon_package & pp, uint i_density, uint a, uint t) const
    {
        return getDustTempProbability(*pp.getPositionCell(), i_density, a, t);
    }

    double getPDAValue(const cell_basic & cell) const
    {
        return cell.getData(data_pos_pda);
    }

    double getGasTemperature(const photon_package & pp) const
    {
        return getGasTemperature(*pp.getPositionCell());
    }

    double getGasTemperature(const cell_basic & cell) const
    {
        return max(TEMP_MIN, cell.getData(data_pos_tg));
    }

    void setPlaneParameter(uint plane_index,
                           double xy_step,
                           double off_xy,
                           double z_step,
                           double off_z,
                           double shift_z,
                           int j,
                           int k,
                           int l,
                           double & tx,
                           double & ty,
                           double & tz)
    {
        switch(plane_index)
        {
            case PROJ_XY:
                if(j != 0)
                {
                    double sg = CMathFunctions::sgn(j);
                    tx = double(j) * xy_step - sg * off_xy;
                }
                else
                    tx = numeric_limits<double>::min();
                if(k != 0)
                {
                    double sg = CMathFunctions::sgn(k);
                    ty = double(k) * xy_step - sg * off_xy;
                }
                else
                    ty = numeric_limits<double>::min();
                if(l != 0)
                {
                    double sg = CMathFunctions::sgn(l);
                    tz = double(l) * z_step - sg * off_z + shift_z;
                }
                else
                    tz = numeric_limits<double>::min();
                tz += shift_z;
                break;

            case PROJ_XZ:
                if(j != 0)
                {
                    double sg = CMathFunctions::sgn(j);
                    tx = double(j) * xy_step - sg * off_xy;
                }
                else
                    tx = numeric_limits<double>::min();
                if(l != 0)
                {
                    double sg = CMathFunctions::sgn(l);
                    ty = double(l) * z_step - sg * off_z + shift_z;
                }
                else
                    ty = numeric_limits<double>::min();
                ty += shift_z;
                if(k != 0)
                {
                    double sg = CMathFunctions::sgn(k);
                    tz = double(k) * xy_step - sg * off_xy;
                }
                else
                    tz = numeric_limits<double>::min();
                break;

            case PROJ_YZ:
                if(l != 0)
                {
                    double sg = CMathFunctions::sgn(l);
                    tx = double(l) * z_step - sg * off_z + shift_z;
                }
                else
                    tx = numeric_limits<double>::min();
                tx += shift_z;
                if(j != 0)
                {
                    double sg = CMathFunctions::sgn(j);
                    ty = double(j) * xy_step - sg * off_xy;
                }
                else
                    ty = numeric_limits<double>::min();
                if(k != 0)
                {
                    double sg = CMathFunctions::sgn(k);
                    tz = double(k) * xy_step - sg * off_xy;
                }
                else
                    tz = numeric_limits<double>::min();
                break;

            default:
                break;
        }
    }

    void fillMidplaneBuffer(double tx, double ty, double tz, uint i_cell)
    {
        photon_package pp = photon_package();
        pp.setPosition(Vector3D(tx, ty, tz));
        if(positionPhotonInGrid(&pp))
        {
            if(plt_gas_dens)
            {
                buffer_gas_dens[i_cell][0] = getGasDensity(pp);
                // Do it only once if only one gas distribution is defined
                if(nr_densities > 1 && size_gd_list == nr_densities)
                    for(uint i_density = 0; i_density < nr_densities; i_density++)
                        buffer_gas_dens[i_cell][i_density + 1] = getGasDensity(pp, i_density);
            }
            if(plt_mol_dens)
            {
                for(uint i_density = 0; i_density < nrOfDensRatios; i_density++)
                        buffer_mol_dens[i_cell][i_density] = getGasDensity(pp, i_density)*getCellAbundance(pp, i_density);
            }
            if(plt_dust_dens)
            {
                buffer_dust_dens[i_cell][0] = getDustDensity(pp);
                // Do it only once if only one dust distribution is defined
                if(nr_densities > 1 && size_dd_list == nr_densities)
                    for(uint i_density = 0; i_density < nr_densities; i_density++)
                        buffer_dust_dens[i_cell][i_density + 1] = getDustDensity(pp, i_density);
            }
            if(plt_gas_temp)
                buffer_gas_temp[i_cell] = getGasTemperature(pp);
            if(plt_dust_temp)
            {
                buffer_dust_temp[i_cell][0] = getDustTemperature(pp);
                // Do it only once if only one dust temperatures is defined
                if(nr_densities > 1 && data_pos_dt_list.size() >= nr_densities)
                    for(uint i_density = 0; i_density < nr_densities; i_density++)
                        buffer_dust_temp[i_cell][i_density + 1] = getDustTemperature(pp, i_density);
            }
            if(plt_rat)
                for(uint i_density = 0; i_density < data_pos_aalg_list.size(); i_density++)
                    buffer_rat[i_cell][i_density] = getAlignedRadius(pp, i_density);
            if(plt_delta)
            {
                double field = getMagField(pp).length();
                double Td = getDustTemperature(pp);
                double Tg = getGasTemperature(pp);
                double dens = getGasDensity(pp);
                double delta = CMathFunctions::calc_delta(field, Td, Tg, dens);
                buffer_delta[i_cell] = delta;
            }
            if(plt_mag)
            {
                Vector3D mag_field = getMagField(pp);
                buffer_mag[i_cell] = mag_field.length();
                buffer_mag_x[i_cell] = mag_field.X();
                buffer_mag_y[i_cell] = mag_field.Y();
                buffer_mag_z[i_cell] = mag_field.Z();
            }
            if(plt_vel)
            {
                Vector3D vel_field = getVelocityField(pp);
                buffer_vel[i_cell] = vel_field.length();
                buffer_vel_x[i_cell] = vel_field.X();
                buffer_vel_y[i_cell] = vel_field.Y();
                buffer_vel_z[i_cell] = vel_field.Z();
            }
            if(plt_larm)
            {
                double field = getMagField(pp).length();
                double Td = getDustTemperature(pp);
                double Tg = getGasTemperature(pp);
                double dens = getGasDensity(pp);
                double a_limit = CMathFunctions::calc_larm_limit(field, Td, Tg, dens, 0.5, 4.1e-19);
                buffer_larm[i_cell] = a_limit;
            }
            if(plt_mach)
            {
                Vector3D vel_field = getVelocityField(pp);
                double gas_temp = getGasTemperature(pp);
                double mach = vel_field.length() / sqrt(con_kB * gas_temp / (mu * m_H));
                buffer_mach[i_cell] = mach;
            }
            if(plt_dust_id)
                buffer_dust_mixture[i_cell] = getDustChoiceID(pp);
            if(plt_amin)
                buffer_dust_amin[i_cell] = getMinGrainRadius(pp);
            if(plt_amax)
                buffer_dust_amax[i_cell] = getMaxGrainRadius(pp);
            if(plt_size_param)
                buffer_dust_size_param[i_cell] = getGrainSizeParam(pp);
            if(plt_rad_field1)
                for(uint i_comp = 0; i_comp < nr_rad_field_comp; i_comp++)
                {
                    for(uint wID = 0; wID < WL_STEPS; wID++)
                    {
                        double val = 0;
                        switch(i_comp)
                        {
                            default:
                                val = getRadiationField(pp, wID);
                                break;

                            case 1:
                                val = getRadiationFieldX(pp, wID);
                                break;

                            case 2:
                                val = getRadiationFieldY(pp, wID);
                                break;

                            case 3:
                                val = getRadiationFieldZ(pp, wID);
                                break;
                        }
                        buffer_rad_field[i_cell][wID][i_comp] = val;
                    }
                }
            if(plt_g_zero1)
                buffer_g_zero1[i_cell] = getGZero(pp);

            if(plt_u_rad)
                buffer_u_rad[i_cell] = getUrad(pp);

            if(plt_n_th)
                buffer_n_th[i_cell] = getThermalElectronDensity(pp);
            if(plt_T_e)
                buffer_T_e[i_cell] = getElectronTemperature(pp);
            if(plt_n_cr)
                buffer_n_cr[i_cell] = getCRElectronDensity(pp);
            if(plt_g_min)
                buffer_g_min[i_cell] = getGammaMin(pp);
            if(plt_g_max)
                buffer_g_max[i_cell] = getGammaMax(pp);
            if(plt_p)
                buffer_p[i_cell] = getPowerLawIndex(pp);
            if(plt_avg_dir)
                buffer_avg_dir[i_cell] = getAvgDir(pp);
            if(plt_avg_th)
                buffer_avg_th[i_cell] = getAvgTheta(pp);
        }
        else
        {
            if(plt_gas_dens)
            {
                buffer_gas_dens[i_cell][0] = 0;
                if(nr_densities > 1 && size_gd_list == nr_densities)
                    for(uint i_density = 1; i_density <= nr_densities; i_density++)
                        buffer_gas_dens[i_cell][i_density] = 0;
            }
            if(plt_mol_dens)
            {
                for(uint i_density = 0; i_density < nrOfDensRatios; i_density++)
                        buffer_mol_dens[i_cell][i_density] = getCellAbundance(pp, i_density);
            }
            if(plt_dust_dens)
            {
                buffer_dust_dens[i_cell][0] = 0;
                if(nr_densities > 1 && size_dd_list == nr_densities)
                    for(uint i_density = 1; i_density <= nr_densities; i_density++)
                        buffer_dust_dens[i_cell][i_density] = 0;
            }
            if(plt_gas_temp)
                buffer_gas_temp[i_cell] = 0;
            if(plt_dust_temp)
            {
                buffer_dust_temp[i_cell][0] = 0;
                if(nr_densities > 1 && data_pos_dt_list.size() >= nr_densities)
                    for(uint i_density = 1; i_density <= nr_densities; i_density++)
                        buffer_dust_temp[i_cell][i_density] = 0;
            }
            if(plt_rat)
            {
                for(uint i_density = 0; i_density < data_pos_aalg_list.size(); i_density++)
                    buffer_rat[i_cell][i_density] = 0;
            }
            if(plt_delta)
                buffer_delta[i_cell] = 0;
            if(plt_mag)
            {
                buffer_mag[i_cell] = 0;
                buffer_mag_x[i_cell] = 0;
                buffer_mag_y[i_cell] = 0;
                buffer_mag_z[i_cell] = 0;
            }
            if(plt_vel)
            {
                buffer_vel[i_cell] = 0;
                buffer_vel_x[i_cell] = 0;
                buffer_vel_y[i_cell] = 0;
                buffer_vel_z[i_cell] = 0;
            }
            if(plt_larm)
                buffer_larm[i_cell] = 0;
            if(plt_mach)
                buffer_mach[i_cell] = 0;
            if(plt_dust_id)
                buffer_dust_mixture[i_cell] = 0;
            if(plt_amin)
                buffer_dust_amin[i_cell] = 0;
            if(plt_amax)
                buffer_dust_amax[i_cell] = 0;
            if(plt_size_param)
                buffer_dust_size_param[i_cell] = 0;
            if(plt_rad_field1)
                for(uint i_comp = 0; i_comp < nr_rad_field_comp; i_comp++)
                    for(uint wID = 0; wID < WL_STEPS; wID++)
                        buffer_rad_field[i_cell][wID][i_comp] = 0;
            if(plt_g_zero1)
                buffer_g_zero1[i_cell] = 0;
            if(plt_u_rad)
                buffer_u_rad[i_cell] = 0;
            if(plt_n_th)
                buffer_n_th[i_cell] = 0;
            if(plt_T_e)
                buffer_T_e[i_cell] = 0;
            if(plt_n_cr)
                buffer_n_cr[i_cell] = 0;
            if(plt_g_min)
                buffer_g_min[i_cell] = 0;
            if(plt_g_max)
                buffer_g_max[i_cell] = 0;
            if(plt_p)
                buffer_p[i_cell] = 0;
            if(plt_avg_dir)
                buffer_avg_dir[i_cell] = 0;
            if(plt_avg_th)
                buffer_avg_th[i_cell] = 0;
        }
    }

    bool writeSpecialLines(string path);
    bool writeAMIRAFiles(string path, parameters & param, uint bins);
    bool writeMidplaneFits(string data_path, parameters & param, uint bins, bool all = false);

    void updateMidplaneString(char * str_1, char * str_2, uint counter)
    {
#ifdef WINDOWS
        sprintf_s(str_1, "MIDPLANE%i", counter);
        sprintf_s(str_2, "quantity of %i. image", counter);
#else
        sprintf(str_1, "MIDPLANE%i", counter);
        sprintf(str_2, "quantity of %i. image", counter);
#endif
    }

    string getDensityString(string quantity, uint counter)
    {
        char str_char[256];
#ifdef WINDOWS
        sprintf_s(str_char, quantity.c_str(), counter);
#else
        sprintf(str_char, quantity.c_str(), counter);
#endif
        string tmp_str(str_char);
        return tmp_str;
    }

    virtual bool saveBinaryGridFile(string filename) = 0;

    virtual bool saveBinaryGridFile(string filename, ushort id, ushort data_size) = 0;

    virtual bool loadGridFromBinrayFile(parameters & param) = 0;

    virtual bool loadGridFromBinrayFile(parameters & param, uint data_len) = 0;

    virtual double getVolume(const cell_basic & cell) const = 0;

    double getVolume(const photon_package & pp) const
    {
        return getVolume(*pp.getPositionCell());
    }

    inline double getGasDensity(const cell_basic & cell) const
    {
        double sum = 0;
        for(uint i_density = 0; i_density < size_gd_list; i_density++)
            sum += getGasDensity(cell, i_density);
        return sum;
    }

    inline double getGasDensity(const cell_basic & cell, uint i_density) const
    {
        if(size_gd_list > i_density)
            return cell.getData(data_pos_gd_list[i_density]);
        else
            return 0;
    }

    inline double getGasDensity(const photon_package & pp) const
    {
        return getGasDensity(*pp.getPositionCell());
    }

    inline double getGasDensity(const photon_package & pp, uint i_density) const
    {
        return getGasDensity(*pp.getPositionCell(), i_density);
    }

    inline double getGasNumberDensity(const cell_basic & cell) const
    {
        double sum = 0;
        for(uint i_density = 0; i_density < size_gd_list; i_density++)
            sum += cell.getData(data_pos_gd_list[i_density]);
        if(gas_is_mass_density)
            sum /= (mu * m_H);
        return sum;
    }

    inline double getGasNumberDensity(const cell_basic & cell, uint i_density) const
    {
        double dens = 0;
        if(size_gd_list > i_density)
            dens = cell.getData(data_pos_gd_list[i_density]);
        if(gas_is_mass_density)
            dens /= (mu * m_H);
        return dens;
    }

    inline double getGasNumberDensity(const photon_package & pp) const
    {
        return getGasNumberDensity(*pp.getPositionCell());
    }

    inline double getGasNumberDensity(const photon_package & pp, uint i_density) const
    {
        return getGasNumberDensity(*pp.getPositionCell(), i_density);
    }

    inline double getGasMassDensity(const cell_basic & cell) const
    {
        double sum = 0;
        for(uint i_density = 0; i_density < size_gd_list; i_density++)
            sum += cell.getData(data_pos_gd_list[i_density]);
        if(!gas_is_mass_density)
            sum *= (mu * m_H);
        return sum;
    }

    inline double getGasMassDensity(const cell_basic & cell, uint i_density) const
    {
        double dens = 0;
        if(size_gd_list > i_density)
            dens = cell.getData(data_pos_gd_list[i_density]);
        if(!gas_is_mass_density)
            dens *= (mu * m_H);
        return dens;
    }

    inline double getGasMassDensity(const photon_package & pp) const
    {
        return getGasMassDensity(*pp.getPositionCell());
    }

    inline double getGasMassDensity(const photon_package & pp, uint i_density) const
    {
        return getGasMassDensity(*pp.getPositionCell(), i_density);
    }

    bool useDustChoice()
    {
        if(size_gd_list > 1 || size_dd_list > 1)
            return false;
        return true;
    }

    bool useConstantGrainSizes()
    {
        if(data_pos_amin != MAX_UINT || data_pos_amax != MAX_UINT || data_pos_size_param != MAX_UINT)
            return false;
        return true;
    }

    bool useDustDensities()
    {
        if(size_dd_list > 1)
            return true;
        return false;
    }

    void setDustDensity(cell_basic * cell, double val)
    {
        if(size_dd_list > 0)
            for(uint i_density = 0; i_density < size_dd_list; i_density++)
                cell->setData(data_pos_dd_list[i_density], val);
        else
            for(uint i_density = 0; i_density < size_gd_list; i_density++)
                cell->setData(data_pos_gd_list[i_density], val);
    }

    void setDustDensity(cell_basic * cell, uint i_density, double val)
    {
        if(size_dd_list > 0)
            cell->setData(data_pos_dd_list[i_density], val);
        else
            cell->setData(data_pos_gd_list[i_density], val);
    }

    void setDustDensity(photon_package * pp, double val)
    {
        setDustDensity(pp->getPositionCell(), val);
    }

    void setDustDensity(photon_package * pp, uint i_density, double val)
    {
        setDustDensity(pp->getPositionCell(), i_density, val);
    }

    double getDustDensity(const cell_basic & cell) const
    {
        double sum = 0;
        if(size_dd_list > 0)
        {
            for(uint i_density = 0; i_density < size_dd_list; i_density++)
                sum += cell.getData(data_pos_dd_list[i_density]);
            return sum;
        }
        else
            return 0;
    }

    double getDustDensity(const cell_basic & cell, uint i_density) const
    {
        if(size_dd_list > i_density)
            return cell.getData(data_pos_dd_list[i_density]);
        else
            return 0;
    }

    double getDustDensity(const photon_package & pp) const
    {
        return getDustDensity(*pp.getPositionCell());
    }

    double getDustDensity(const photon_package & pp, uint i_density) const
    {
        return getDustDensity(*pp.getPositionCell(), i_density);
    }

    double getRelativeDustDensity(const cell_basic & cell, uint i_density) const
    {
        if(getDustDensity(cell) != 0)
            return getDustDensity(cell, i_density) / getDustDensity(cell);
        else if(getGasDensity(cell) != 0)
            return getGasDensity(cell, i_density) / getGasDensity(cell);
        else
            return 0;
    }

    void adjustDustDensity(cell_basic * cell, uint i_density, double factor)
    {
        if(size_dd_list > 0)
        {
            double dust_dens = getDustDensity(*cell, i_density);
            cell->setData(data_pos_dd_list[i_density], dust_dens * factor);
        }
        else
        {
            double gas_dens = getGasDensity(*cell, i_density);
            cell->setData(data_pos_gd_list[i_density], gas_dens * factor);
        }
    }

    virtual bool positionPhotonInGridTest(photon_package * pp)
    {
        return false;
    };

    bool isVelocityFieldAvailable()
    {
        if(data_pos_vx == MAX_UINT || data_pos_vy == MAX_UINT || data_pos_vz == MAX_UINT)
            return false;
        return true;
    }

    bool isTurbulentVelocityAvailable()
    {
        if(data_pos_vt == MAX_UINT)
            return false;
        return true;
    }

    Vector3D getVelocityField(const photon_package & pp) const
    {
        if(data_pos_vx == MAX_UINT || data_pos_vy == MAX_UINT || data_pos_vz == MAX_UINT)
            return Vector3D();

        const cell_basic & cell = *pp.getPositionCell();
        Vector3D tmp_dir(cell.getData(data_pos_vx), cell.getData(data_pos_vy), cell.getData(data_pos_vz));
        // Rotate vector from cell center to position
        return rotateToCenter(pp, tmp_dir, true, true);
    }

    Vector3D getVelocityField(const cell_basic & cell) const
    {
        if(data_pos_vx == MAX_UINT || data_pos_vy == MAX_UINT || data_pos_vz == MAX_UINT)
            return Vector3D();

        return Vector3D(cell.getData(data_pos_vx), cell.getData(data_pos_vy), cell.getData(data_pos_vz));
    }

    double getCellAbundance(const photon_package & pp, uint id) const
    {
        if(id > nrOfDensRatios - 1)
            return 0;

        return getCellAbundance(*pp.getPositionCell(), id);
    }

    double getCellAbundance(const cell_basic & cell, uint id) const
    {
        if(id > nrOfDensRatios - 1)
            return 0;

        uint pos = pos_GasSpecRatios[id];
        return cell.getData(pos);
    }

    double getOpiateIDParameter(cell_basic * cell, uint id)
    {
        if(id > nrOfOpiateIDs - 1)
            return 0;

        uint pos = pos_OpiateIDS[id];

        return cell->getData(pos);
    }

    Vector3D getMagField(const cell_basic & cell) const
    {
        return Vector3D(cell.getData(data_pos_mx), cell.getData(data_pos_my), cell.getData(data_pos_mz));
    }

    Vector3D getMagField(const photon_package & pp) const
    {
        const cell_basic & cell = *pp.getPositionCell();
        Vector3D tmp_dir(cell.getData(data_pos_mx), cell.getData(data_pos_my), cell.getData(data_pos_mz));
        // Rotate vector from cell center to position
        return rotateToCenter(pp, tmp_dir, true, true);
    }

    void setMagField(cell_basic * cell, const Vector3D & mag)
    {
        cell->setData(data_pos_mx, mag.X());
        cell->setData(data_pos_my, mag.Y());
        cell->setData(data_pos_mz, mag.Z());
    }

    virtual bool next(photon_package * pp) = 0;

    double getThetaMag(const photon_package & pp) const
    {
        return getAngleTheta(pp.getDirection(), getMagField(pp));
    }

    double getPhiMag(const photon_package & pp) const
    {
        // 0 deg are in the e_y direction
        return getAnglePhi(pp.getEX(), pp.getEY(), getMagField(pp)) - PI2;
    }

    double getTheta(const cell_basic & cell, Vector3D & dir) const
    {
        return getAngleTheta(dir, getMagField(cell));
    }

    double getThetaPhoton(const photon_package & pp, Vector3D & dir) const
    {
        return getAngleTheta(pp.getDirection(), dir);
    }

    bool getDustIsMassDensity() const
    {
        return dust_is_mass_density;
    }

    bool getGasIsMassDensity() const
    {
        return gas_is_mass_density;
    }

    bool isRadiationFieldAvailable() const
    {
        if(data_pos_rx_list.empty() || data_pos_ry_list.empty() || data_pos_rz_list.empty() ||
           data_pos_rf_list.empty())
            return false;
        return true;
    }

    double getTotalGasMass() const
    {
        return total_gas_mass;
    };

    void setDustTemperatureRange(double _min_dust_temp, double _max_dust_temp)
    {
        max_dust_temp = _max_dust_temp;
        min_dust_temp = _min_dust_temp;
    }

    void setalignedRadiusRange(double a_min, double a_max)
    {
        aalg_min = a_min;
        aalg_max = a_max;
    }

    bool doPDA(parameters & param, uint pda_id);

    uint getTemperatureFieldInformation() const
    {
        // Check which kind of temperature calculation the grid supports
        if(multi_temperature_entries > nr_densities && data_pos_dt_list.size() == multi_temperature_entries)
            return TEMP_FULL;
        else if(data_pos_dt_list.size() == nr_densities)
            return TEMP_EFF;
        else if(data_pos_dt_list.size() == 1)
            return TEMP_SINGLE;
        else if(data_pos_dt_list.size() == 0)
            return TEMP_EMPTY;
        else if(data_pos_dt_list.size() == stochastic_temperature_entries)
            return TEMP_STOCH;
        else
            return MAX_UINT;
    }

    bool setDataPositionsVariable()
    {
        nrOfDensRatios = 0;

        for(uint i = 0; i < data_offset; i++)
        {
            switch(data_ids[i])
            {
                case GRIDgas_dens:
                    if(!data_pos_gd_list.empty())
                    {
                        if(data_pos_id != MAX_UINT)
                        {
                            cout << "\nERROR: Multiple densities and dust choices cannot "
                                    "be combined!"
                                 << endl;
                            return false;
                        }
                        if(gas_is_mass_density == true)
                        {
                            cout << "\nERROR: Gas number densities cannot be combined "
                                    "with gas mass densities!"
                                 << endl;
                            return false;
                        }
                    }
                    data_pos_gd_list.push_back(i);
                    gas_is_mass_density = false;
                    break;

                case GRIDgas_mdens:
                    if(!data_pos_gd_list.empty())
                    {
                        if(data_pos_id != MAX_UINT)
                        {
                            cout << "\nERROR: Multiple densities and dust choices cannot "
                                    "be combined!"
                                 << endl;
                            return false;
                        }
                        if(gas_is_mass_density == false)
                        {
                            cout << "\nERROR: Gas mass densities cannot be combined with "
                                    "gas number densities!"
                                 << endl;
                            return false;
                        }
                    }
                    data_pos_gd_list.push_back(i);
                    gas_is_mass_density = true;
                    break;

                case GRIDdust_dens:
                    if(!data_pos_dd_list.empty())
                    {
                        if(data_pos_id != MAX_UINT)
                        {
                            cout << "\nERROR: Multiple densities and dust choices cannot "
                                    "be combined!"
                                 << endl;
                            return false;
                        }
                        if(dust_is_mass_density == true)
                        {
                            cout << "\nERROR: Dust number densities cannot be combined "
                                    "with dust mass densities!"
                                 << endl;
                            return false;
                        }
                    }
                    data_pos_dd_list.push_back(i);
                    dust_is_mass_density = false;
                    break;

                case GRIDdust_mdens:
                    if(!data_pos_dd_list.empty())
                    {
                        if(data_pos_id != MAX_UINT)
                        {
                            cout << "\nERROR: Multiple densities and dust choices cannot "
                                    "be combined!"
                                 << endl;
                            return false;
                        }
                        if(dust_is_mass_density == false)
                        {
                            cout << "\nERROR: Dust mass densities cannot be combined "
                                    "with dust number densities!"
                                 << endl;
                            return false;
                        }
                    }
                    data_pos_dd_list.push_back(i);
                    dust_is_mass_density = true;
                    break;

                case GRIDdust_temp:
                    data_pos_dt_list.push_back(i);
                    break;

                case GRIDgas_temp:
                    if(data_pos_tg != MAX_UINT)
                    {
                        cout << "\nERROR: Grid ID " << GRIDgas_temp << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_tg = i;
                    break;

                case GRIDmx:
                    if(data_pos_mx != MAX_UINT)
                    {
                        cout << "\nERROR: Grid ID " << GRIDmx << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_mx = i;
                    break;

                case GRIDmy:
                    if(data_pos_my != MAX_UINT)
                    {
                        cout << "\nERROR: Grid ID " << GRIDmy << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_my = i;
                    break;

                case GRIDmz:
                    if(data_pos_mz != MAX_UINT)
                    {
                        cout << "\nERROR: Grid ID " << GRIDmz << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_mz = i;
                    break;

                case GRIDvx:
                    if(data_pos_vx != MAX_UINT)
                    {
                        cout << "\nERROR: Grid ID " << GRIDvx << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_vx = i;
                    break;

                case GRIDvy:
                    if(data_pos_vy != MAX_UINT)
                    {
                        cout << "\nERROR: Grid ID " << GRIDvy << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_vy = i;
                    break;

                case GRIDvz:
                    if(data_pos_vz != MAX_UINT)
                    {
                        cout << "\nERROR: Grid ID " << GRIDvz << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_vz = i;
                    break;

                case GRIDpx:
                    if(data_pos_px != MAX_UINT)
                    {
                        cout << "\nERROR: Grid ID " << GRIDpx << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_px = i;
                    break;

                case GRIDpy:
                    if(data_pos_py != MAX_UINT)
                    {
                        cout << "\nERROR: Grid ID " << GRIDpy << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_py = i;
                    break;

                case GRIDpz:
                    if(data_pos_pz != MAX_UINT)
                    {
                        cout << "\nERROR: Grid ID " << GRIDpz << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_pz = i;
                    break;

                case GRIDa_alg:
                    data_pos_aalg_list.push_back(i);
                    break;

                case GRIDa_min:
                    if(data_pos_amin != MAX_UINT)
                    {
                        cout << "\nERROR: Grid ID " << GRIDa_min << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_amin = i;
                    break;

                case GRIDa_max:
                    if(data_pos_amax != MAX_UINT)
                    {
                        cout << "\nERROR: Grid ID " << GRIDa_max << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_amax = i;
                    break;

                case GRIDq:
                    if(data_pos_size_param != MAX_UINT)
                    {
                        cout << "\nERROR: Grid ID " << GRIDq << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_size_param = i;
                    break;

                case GRIDv_turb:
                    if(data_pos_vt != MAX_UINT)
                    {
                        cout << "\nERROR: Grid ID " << GRIDvx << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_vt = i;
                    break;

                case GRIDn_th:
                    if(data_pos_n_th != MAX_UINT)
                    {
                        cout << "\nERROR: Grid ID " << GRIDn_th << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_n_th = i;
                    break;

                case GRIDT_e:
                    if(data_pos_T_e != MAX_UINT)
                    {
                        cout << "\nERROR: Grid ID " << GRIDT_e << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_T_e = i;
                    break;

                case GRIDn_cr:
                    if(data_pos_n_cr != MAX_UINT)
                    {
                        cout << "\nERROR: Grid ID " << GRIDn_cr << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_n_cr = i;
                    break;

                case GRIDg_min:
                    if(data_pos_g_min != MAX_UINT)
                    {
                        cout << "\nERROR: Grid ID " << GRIDg_min << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_g_min = i;
                    break;

                case GRIDg_max:
                    if(data_pos_g_max != MAX_UINT)
                    {
                        cout << "\nERROR: Grid ID " << GRIDg_max << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_g_max = i;
                    break;

                case GRIDp:
                    if(data_pos_p != MAX_UINT)
                    {
                        cout << "\nERROR: Grid ID " << GRIDp << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_p = i;
                    break;

                case GRIDavg_dir:
                    if(data_pos_avg_dir != MAX_UINT)
                    {
                        cout << "\nERROR: Grid ID " << GRIDavg_dir << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_avg_dir = i;
                    break;

                case GRIDavg_th:
                    if(data_pos_avg_th != MAX_UINT)
                    {
                        cout << "\nERROR: Grid ID " << GRIDavg_th << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_avg_th = i;
                    break;

                case GRIDratio:
                    nrOfDensRatios++;
                    break;

                case GRIDopiate:
                    nrOfOpiateIDs++;
                    break;

                case GRIDdust_id:
                    if(data_pos_gd_list.size() > 1 || data_pos_dd_list.size() > 1)
                    {
                        cout << "\nERROR: Multiple densities and dust choices cannot be "
                                "combined!"
                             << endl;
                        return false;
                    }
                    if(data_pos_id != MAX_UINT)
                    {
                        cout << "\nERROR: Grid ID " << GRIDdust_id << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_id = i;
                    break;

                case GRIDrad:
                    data_pos_rf_list.push_back(i);
                    break;

                case GRIDradx:
                    data_pos_rx_list.push_back(i);
                    break;

                case GRIDrady:
                    data_pos_ry_list.push_back(i);
                    break;

                case GRIDradz:
                    data_pos_rz_list.push_back(i);
                    break;

                default:
                    cout << "\nERROR: Unknown data IDs!" << endl;
                    cout << "         IDs have to be between " << minGRID << " and " << maxGRID << "!"
                         << endl;
                    return false;
            }
        }
        size_gd_list = data_pos_gd_list.size();
        size_dd_list = data_pos_dd_list.size();

        if(size_gd_list == 0)
        {
            cout << "\nERROR: Grid requires a gas density! " << endl;
            return false;
        }

        pos_GasSpecRatios = new uint[nrOfDensRatios];
        uint pos_counter = 0;

        for(uint i = 0; i < data_offset; i++)
        {
            if(data_ids[i] == GRIDratio)
            {
                pos_GasSpecRatios[pos_counter] = i;
                // cout << pos_counter << "\t" << pos_ration[pos_counter] << endl;
                pos_counter++;
            }
        }

        pos_OpiateIDS = new uint[nrOfOpiateIDs];
        pos_counter = 0;

        for(uint i = 0; i < data_offset; i++)
        {
            if(data_ids[i] == GRIDopiate)
            {
                pos_OpiateIDS[pos_counter] = i;
                // cout << pos_counter << "\t" << pos_ration[pos_counter] << endl;
                pos_counter++;
            }
        }

        return true;
    }

    bool createCompatibleTree()
    {
        switch(dataID)
        {
            case 0:
                if(data_offset != 6)
                {
                    cout << "\nERROR: A grid ID of 0 requires a data length of 6!" << endl;
                    cout << "       Cannot create a compatible octree!" << endl;
                    return false;
                }

                data_pos_gd_list.push_back(0);
                data_pos_dt_list.push_back(1);
                data_pos_tg = 2;

                data_pos_mx = 3;
                data_pos_my = 4;
                data_pos_mz = 5;

                data_ids[0] = GRIDgas_dens;
                data_ids[1] = GRIDdust_temp;
                data_ids[2] = GRIDgas_temp;

                data_ids[3] = GRIDmx;
                data_ids[4] = GRIDmy;
                data_ids[5] = GRIDmz;

                break;

            case 1:
                if(data_offset != 7)
                {
                    cout << "\nERROR: A grid ID of 1 requires a data length of 7!" << endl;
                    cout << "       Cannot create a compatible octree!" << endl;
                    return false;
                }

                data_pos_gd_list.push_back(0);
                data_pos_dt_list.push_back(1);
                data_pos_tg = 2;

                data_pos_mx = 3;
                data_pos_my = 4;
                data_pos_mz = 5;

                data_pos_aalg_list.push_back(6);

                data_ids[0] = GRIDgas_dens;
                data_ids[1] = GRIDdust_temp;
                data_ids[2] = GRIDgas_temp;

                data_ids[3] = GRIDmx;
                data_ids[4] = GRIDmy;
                data_ids[5] = GRIDmz;

                data_ids[6] = GRIDa_alg;

                break;

            case 6:
                if(data_offset != 9)
                {
                    cout << "\nERROR: A grid ID of 6 requires a data length of 9!" << endl;
                    cout << "       Cannot create a compatible octree!" << endl;
                    return false;
                }

                data_pos_gd_list.push_back(0);
                data_pos_dt_list.push_back(1);
                data_pos_tg = 2;

                data_pos_mx = 3;
                data_pos_my = 4;
                data_pos_mz = 5;

                data_pos_vx = 6;
                data_pos_vy = 7;
                data_pos_vz = 8;

                data_ids[0] = GRIDgas_dens;
                data_ids[1] = GRIDdust_temp;
                data_ids[2] = GRIDgas_temp;

                data_ids[3] = GRIDmx;
                data_ids[4] = GRIDmy;
                data_ids[5] = GRIDmz;

                data_ids[6] = GRIDvx;
                data_ids[7] = GRIDvy;
                data_ids[8] = GRIDvz;
                break;

            case 7:
                if(data_offset != 10)
                {
                    cout << "\nERROR: A grid ID of 7 requires a data length of 10!" << endl;
                    cout << "       Cannot create a compatible octree!" << endl;
                    return false;
                }

                data_pos_gd_list.push_back(0);
                data_pos_dt_list.push_back(1);
                data_pos_tg = 2;

                data_pos_mx = 3;
                data_pos_my = 4;
                data_pos_mz = 5;

                data_pos_vx = 6;
                data_pos_vy = 7;
                data_pos_vz = 8;

                data_pos_aalg_list.push_back(9);

                data_ids[0] = GRIDgas_dens;
                data_ids[1] = GRIDdust_temp;
                data_ids[2] = GRIDgas_temp;

                data_ids[3] = GRIDmx;
                data_ids[4] = GRIDmy;
                data_ids[5] = GRIDmz;

                data_ids[6] = GRIDvx;
                data_ids[7] = GRIDvy;
                data_ids[8] = GRIDvz;

                data_ids[9] = GRIDa_alg;

                break;

            default:
                cout << "\nERROR: Grid ID = " << dataID << " doesn't match any of the old octree IDs!"
                     << endl;
                cout << "       Cannot create a compatible octree!" << endl;
                return false;
        }

        dataID = GRID_ID_OCT;

        return true;
    }

    uint CheckSynchrotron(parameters & param)
    {
        if(data_pos_n_th == MAX_UINT && data_pos_n_th == MAX_UINT)
        {
            cout << "\nERROR: Neither thermal electrons nor CR electrons are defined in "
                    "grid file!                                            \n"
                 << endl;
            cout << "       No SYNCHROTRON calculation possible." << endl;
            return MAX_UINT;
        }

        if(data_pos_mx == MAX_UINT)
        {
            cout << "\nERROR: Grid contains no magnetic Bx component!" << endl;
            cout << "       No SYNCHROTRON calculation possible." << endl;
            return MAX_UINT;
        }
        if(data_pos_my == MAX_UINT)
        {
            cout << "\nERROR: Grid contains no magnetic By component!" << endl;
            cout << "       No SYNCHROTRON calculation possible." << endl;
            return MAX_UINT;
        }
        if(data_pos_mz == MAX_UINT)
        {
            cout << "\nERROR: Grid contains no magnetic Bz component!" << endl;
            cout << "       No SYNCHROTRON calculation possible." << endl;
            return MAX_UINT;
        }

        if(data_pos_n_th == MAX_UINT)
        {
            cout << "\nWARNING: Grid contains no thermal electron component!" << endl;
            cout << "         Only CR SYNCHROTRON calculation possible." << endl;
        }
        else
        {
            if(data_pos_T_e != MAX_UINT)
            {
                cout << "\nHINT: Grid contains a electron temperature component!" << endl;
                cout << "      This component is currently ignored!          " << endl;
            }
        }

        if(data_pos_n_cr == MAX_UINT)
        {
            cout << "\nWARNING: Grid contains no thermal electron component!         " << endl;
            cout << "         Only Fraraday RM calculations possible." << endl;
        }
        else
        {
            if(data_pos_g_min == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no gamma_min component!" << endl;
                cout << "       No SYNCHROTRON calculation possible." << endl;
                return MAX_UINT;
            }

            if(data_pos_g_max == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no gamma_max component!" << endl;
                cout << "       No SYNCHROTRON calculation possible." << endl;
                return MAX_UINT;
            }

            if(data_pos_p == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no electron power-law index p component!" << endl;
                cout << "       No SYNCHROTRON calculation possible." << endl;
                return MAX_UINT;
            }
        }

        return 0;
    }

    uint CheckOpiate(parameters & param)
    {
        if(data_pos_tg == MAX_UINT)
        {
            cout << "\nERROR: Grid contains no gas temperature!" << endl;
            cout << "       No OPIATE calculation possible." << endl;
            return MAX_UINT;
        }
        return 0;
    }

    uint CheckTemp(parameters & param, uint & tmp_data_offset)
    {
        uint extra_temp_entries = 0;
        if(getTemperatureFieldInformation() == MAX_UINT)
        {
            cout << "\nERROR: The grid does not include the correct information for "
                    "temperature calculations"
                 << endl;
            cout << "       No dust temperature calculation possible (full_dust_temp or "
                    "stochastic heating?)."
                 << endl;
            return MAX_UINT;
        }
        else
        {
            // Calculate the entries for the temperature that have to be added
            if(param.getDustTempMulti())
                extra_temp_entries = multi_temperature_entries;
            else if(param.getStochasticHeatingMaxSize() > 0 && !param.getSaveRadiationField())
                extra_temp_entries = stochastic_temperature_entries;
            else
                extra_temp_entries = nr_densities;

            // Entries that are already in the grid do not need to be added
            if(getTemperatureFieldInformation() == TEMP_SINGLE)
                extra_temp_entries--;
            else if(getTemperatureFieldInformation() == TEMP_EFF)
                extra_temp_entries -= nr_densities;
        }

        // Add entries to grid
        for(uint i_entries = 0; i_entries < extra_temp_entries; i_entries++)
        {
            data_pos_dt_list.push_back(data_offset + tmp_data_offset);
            data_ids.push_back(GRIDdust_temp);
            tmp_data_offset++;
        }

        if(param.getSaveRadiationField())
            if(data_pos_rx_list.size() != 0 || data_pos_ry_list.size() != 0 || data_pos_rz_list.size() != 0 ||
               data_pos_rf_list.size() != 0)
            {
                cout << "\nERROR: The grid includes partial/broken information about a "
                        "radiation field!"
                     << endl;
                cout << "       No dust temperature calculation possible." << endl;
                return MAX_UINT;
            }

        if(data_pos_tg == MAX_UINT)
        {
            if(param.getAdjTgas() != 0)
            {
                data_pos_tg = data_offset + tmp_data_offset;
                data_ids.push_back(GRIDgas_temp);
                tmp_data_offset++;
                // cout << "Create entries for gas temperature   : done" << endl;
            }
            else
            {
                param.setAdjTgas(1.0);
                data_pos_tg = data_offset + tmp_data_offset;
                data_ids.push_back(GRIDgas_temp);
                tmp_data_offset++;
                cout << SEP_LINE;
                cout << "\nHINT: No gas temperature found in grid." << endl;
                cout << "    Add entry and set gas temperature to dust temperature after "
                        "calculation!"
                     << endl;
                cout << SEP_LINE;
            }
        }
        return 0;
    }

    uint CheckRat(parameters & param, uint & tmp_data_offset)
    {
        if(data_pos_aalg_list.empty())
        {
            for(uint i_density = 0; i_density < nr_densities; i_density++)
            {
                data_pos_aalg_list.push_back(data_offset + tmp_data_offset);
                data_ids.push_back(GRIDa_alg);
                tmp_data_offset++;
            }
        }

        if(data_pos_avg_dir == MAX_UINT)
        {
            data_pos_avg_dir = data_offset + tmp_data_offset;
            data_ids.push_back(GRIDavg_dir);
            tmp_data_offset++;
        }

        if(data_pos_avg_th == MAX_UINT)
        {
            data_pos_avg_th = data_offset + tmp_data_offset;
            data_ids.push_back(GRIDavg_th);
            tmp_data_offset++;
        }

        if(getTemperatureFieldInformation() == TEMP_EMPTY)
        {
            cout << "\nERROR: Grid contains no dust temperature!" << endl;
            cout << "       No RAT calculation possible." << endl;
            return MAX_UINT;
        }
        else if(getTemperatureFieldInformation() == MAX_UINT)
        {
            cout << "\nERROR: The grid does not include the information for temperature "
                    "calculations"
                 << endl;
            cout << "       No RAT calculation possible." << endl;
            return MAX_UINT;
        }

        if(data_pos_tg == MAX_UINT)
        {
            cout << "\nERROR: Grid contains no gas temperature!" << endl;
            cout << "       No RAT calculation possible." << endl;
            return MAX_UINT;
        }
        if(data_pos_mx == MAX_UINT)
        {
            cout << "\nWARNING: Grid contains no magnetic Bx component!" << endl;
            cout << "         No follow up calculations possible." << endl;
        }
        if(data_pos_my == MAX_UINT)
        {
            cout << "\nWARNING: Grid contains no magnetic By component!" << endl;
            cout << "         No follow up calculations possible." << endl;
        }
        if(data_pos_mz == MAX_UINT)
        {
            cout << "\nWARNING: Grid contains no magnetic Bz component!" << endl;
            cout << "         No follow up calculations possible." << endl;
        }
        return 0;
    }

    uint CheckDustEmission(parameters & param)
    {
        // Check if stochastic heating temperatures are saved in grid
        if(data_pos_dt_list.size() > nr_densities && data_pos_dt_list.size() < multi_temperature_entries)
            stochastic_temperature_entries = data_pos_dt_list.size();

        if(getTemperatureFieldInformation() == TEMP_EMPTY)
        {
            cout << "\nERROR: Grid contains no dust temperature!" << endl;
            cout << "       No dust emission possible." << endl;
            return MAX_UINT;
        }
        else if(getTemperatureFieldInformation() == MAX_UINT)
        {
            cout << "\nERROR: The grid does not include the information for temperature "
                    "calculations"
                 << endl;
            cout << "       No dust emission possible." << endl;
            return MAX_UINT;
        }

        if(getTemperatureFieldInformation() == TEMP_STOCH)
            param.setStochasticHeatingMaxSize(0.0);

        if(param.getStochasticHeatingMaxSize())
        {
            if(data_pos_rf_list.size() != WL_STEPS)
            {
                cout << "\nERROR: The grid includes partial/no information about a "
                        "radiation field!"
                     << endl;
                cout << "       No dust emission with stochastic heating possible." << endl;
                return MAX_UINT;
            }
        }

        if(!data_pos_rf_list.empty() && data_pos_rf_list.size() != WL_STEPS)
        {
            cout << "\nERROR: The grid includes partial/no information about a radiation "
                    "field!"
                 << endl;
            cout << "       No dust emission possible." << endl;
            return MAX_UINT;
        }

        if(param.getAlign() != 0 && !param.getAligPA())
        {
            if(data_pos_tg == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no gas temperature!" << endl;
                cout << "       No dust emission with aligned dust grains possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_mx == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no magnetic Bx component!" << endl;
                cout << "       No dust emission with aligned dust grains possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_my == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no magnetic By component!" << endl;
                cout << "       No dust emission with aligned dust grains possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_mz == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no magnetic Bz component!" << endl;
                cout << "       No dust emission with aligned dust grains possible." << endl;
                return MAX_UINT;
            }
        }

        if(param.getAligGOLD())
        {
            if(data_pos_vx == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no velocity vx component!" << endl;
                cout << "        No dust emission with GOLD alignment possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_vy == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no velocity vy component!" << endl;
                cout << "        No dust emission with GOLD alignment possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_vz == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no velocity vz component!" << endl;
                cout << "        No dust emission with GOLD alignment possible." << endl;
                return MAX_UINT;
            }
        }

        if(param.getAligRAT())
        {
            if(data_pos_aalg_list.empty())
            {
                cout << "\nERROR: Grid contains no minimum alignment radius for RATs!" << endl;
                cout << "        No dust emission with RAT alignment possible." << endl;
                return MAX_UINT;
            }
            else if(data_pos_aalg_list.size() != 1 && data_pos_aalg_list.size() != nr_densities)
            {
                cout << "\nERROR: Grid contains not the correct amount of minimum alignment radii for "
                        "RATs!"
                     << endl;
                cout << "        No dust emission with RAT alignment possible." << endl;
                return MAX_UINT;
            }
        }
        return 0;
    }

    uint CheckDustScattering(parameters & param)
    {
        if(param.getAlign() != 0 && !param.getAligPA())
        {
            if(data_pos_tg == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no gas temperature!           " << endl;
                cout << "       No dust scattering calculations with aligned dust grains "
                        "possible."
                     << endl;
                return MAX_UINT;
            }
            if(data_pos_mx == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no magnetic Bx component!     " << endl;
                cout << "       No dust scattering calculations with aligned dust grains "
                        "possible."
                     << endl;
                return MAX_UINT;
            }
            if(data_pos_my == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no magnetic By component!     " << endl;
                cout << "       No dust scattering calculations with aligned dust grains "
                        "possible."
                     << endl;
                return MAX_UINT;
            }
            if(data_pos_mz == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no magnetic Bz component!     " << endl;
                cout << "       No dust scattering calculations with aligned dust grains "
                        "possible."
                     << endl;
                return MAX_UINT;
            }
        }

        if(param.getAligGOLD())
        {
            if(data_pos_vx == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no velocity vx component!  " << endl;
                cout << "        No dust scattering calculations with GOLD alignment "
                        "possible."
                     << endl;
                return MAX_UINT;
            }
            if(data_pos_vy == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no velocity vy component!  " << endl;
                cout << "        No dust scattering calculations with GOLD alignment "
                        "possible."
                     << endl;
                return MAX_UINT;
            }
            if(data_pos_vz == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no velocity vz component!  " << endl;
                cout << "        No dust scattering calculations with GOLD alignment "
                        "possible."
                     << endl;
                return MAX_UINT;
            }
        }

        if(param.getAligRAT())
        {
            if(data_pos_aalg_list.empty())
            {
                cout << "\nERROR: Grid contains no minimum alignment radius for RATs!" << endl;
                cout << "        No dust scattering calculations with RAT alignment "
                        "possible."
                     << endl;
                return MAX_UINT;
            }
            else if(data_pos_aalg_list.size() != 1 && data_pos_aalg_list.size() != nr_densities)
            {
                cout << "\nERROR: Grid contains not the correct amount of minimum alignment radii for "
                        "RATs!"
                     << endl;
                cout << "        No dust scattering calculations with RAT alignment "
                        "possible."
                     << endl;
                return MAX_UINT;
            }
        }

        return 0;
    }

    uint CheckRadiationForce(parameters & param)
    {
        if(getTemperatureFieldInformation() == TEMP_EMPTY)
        {
            cout << "\nERROR: Grid contains no dust temperature!" << endl;
            cout << "       No FORCE calculation possible." << endl;
            return MAX_UINT;
        }
        else if(getTemperatureFieldInformation() == MAX_UINT)
        {
            cout << "\nERROR: The grid does not include the information for temperature "
                    "calculations"
                 << endl;
            cout << "       No FORCE calculation possible." << endl;
            return MAX_UINT;
        }

        if(!data_pos_rf_list.empty() && data_pos_rf_list.size() != WL_STEPS)
        {
            cout << "\nERROR: The grid includes partial/no information about a radiation "
                    "field!"
                 << endl;
            cout << "       No FORCE calculation possible." << endl;
            return MAX_UINT;
        }

        return 0;
    }

    uint CheckLineEmission(parameters & param)
    {
        if(param.getTotalNrOfDustComponents() != 0)
        {
            if(getTemperatureFieldInformation() == TEMP_EMPTY)
            {
                cout << "\nERROR: Grid contains no dust temperature!" << endl;
                cout << "       No line transfer including dust emission possible." << endl;
                return MAX_UINT;
            }
            else if(getTemperatureFieldInformation() == MAX_UINT)
            {
                cout << "\nERROR: The grid does not include the information for "
                        "temperature calculations"
                     << endl;
                cout << "       No line transfer including dust emission possible." << endl;
                return MAX_UINT;
            }
        }

        if(data_pos_tg == MAX_UINT && !param.isGasSpeciesLevelPopMC())
        {
            cout << "\nERROR: Grid contains no gas temperature!" << endl;
            cout << "       No line transfer with possible.  " << endl;
            return MAX_UINT;
        }

        if(velocity_field_needed)
        {
            if(data_pos_mx == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no magnetic Bx component!" << endl;
                cout << "       No line transfer possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_my == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no magnetic By component!" << endl;
                cout << "       No line transfer possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_mz == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no magnetic Bz component!" << endl;
                cout << "       No line transfer possible." << endl;
                return MAX_UINT;
            }
        }

        // Velocity field should be simply zero then.
        // if(param.getKeplerStarMass() == 0)
        // {
        //     if(data_pos_vx == MAX_UINT)
        //     {
        //         cout << "\nERROR: Grid contains no velocity vx component!" << endl;
        //         cout << "        No line transfer possible." << endl;
        //         return MAX_UINT;
        //     }
        //     if(data_pos_vy == MAX_UINT)
        //     {
        //         cout << "\nERROR: Grid contains no velocity vy component!" << endl;
        //         cout << "        No line transfer possible." << endl;
        //         return MAX_UINT;
        //     }
        //     if(data_pos_vz == MAX_UINT)
        //     {
        //         cout << "\nERROR: Grid contains no velocity vz component!" << endl;
        //         cout << "        No line transfer possible." << endl;
        //         return MAX_UINT;
        //     }
        // }
        return 0;
    }

    uint CheckProbing(parameters & param)
    {
        if(getTemperatureFieldInformation() == TEMP_EMPTY)
        {
            cout << "\nERROR: Grid contains no dust temperature!" << endl;
            cout << "       No LOS analysis with aligned dust grains possible." << endl;
            return MAX_UINT;
        }
        else if(getTemperatureFieldInformation() == MAX_UINT)
        {
            cout << "\nERROR: The grid does not include the information for temperature "
                    "calculations"
                 << endl;
            cout << "       No LOS analysis with aligned dust grains possible." << endl;
            return MAX_UINT;
        }

        if(!data_pos_rf_list.empty() && data_pos_rf_list.size() != WL_STEPS)
        {
            cout << "\nERROR: The grid includes partial/no information about a radiation "
                    "field!"
                 << endl;
            cout << "       No LOS analysis with aligned dust grains possible." << endl;
            return MAX_UINT;
        }

        if(param.getAlign() != 0 && !param.getAligPA())
        {
            if(data_pos_tg == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no gas temperature!" << endl;
                cout << "       No LOS analysis with aligned dust grains possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_mx == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no magnetic Bx component!" << endl;
                cout << "       No LOS analysis with aligned dust grains possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_my == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no magnetic By component!" << endl;
                cout << "       No LOS analysis with aligned dust grains possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_mz == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no magnetic Bz component!" << endl;
                cout << "       No LOS analysis with aligned dust grains possible." << endl;
                return MAX_UINT;
            }
        }

        if(param.getAligGOLD())
        {
            if(data_pos_vx == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no velocity vx component!" << endl;
                cout << "        No LOS analysis with GOLD alignment possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_vy == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no velocity vy component!" << endl;
                cout << "        No LOS analysis with GOLD alignment possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_vz == MAX_UINT)
            {
                cout << "\nERROR: Grid contains no velocity vz component!" << endl;
                cout << "        No LOS analysis with GOLD alignment possible." << endl;
                return MAX_UINT;
            }
        }
        return 0;
    }

  protected:
    // uint grid_type;
    ulong max_cells;
    uint max_data;
    Vector3D meanBdir;
    Vector3D meanVdir;

    char * basic_path;

    double max_gas_dens;
    double min_gas_dens;

    double max_dust_dens;
    double min_dust_dens;

    double max_gas_temp;
    double min_gas_temp;

    double max_dust_temp;
    double min_dust_temp;

    double max_larm_limit;
    double min_larm_limit;

    double min_delta;
    double max_delta;

    double max_mach;
    double min_mach;

    double aalg_min;
    double aalg_max;

    double a_min_min;
    double a_min_max;
    double a_max_min;
    double a_max_max;

    double size_param_min;
    double size_param_max;

    uint dust_id_min;
    uint dust_id_max;

    double min_pres;
    double max_pres;

    double max_vel;
    double min_vel;

    double max_mag;
    double min_mag;

    double max_value;

    double max_len;
    double min_len;

    double min_n_th;
    double min_T_e;
    double min_n_cr;
    double min_g_min;
    double min_g_max;
    double min_p;

    double max_n_th;
    double max_T_e;
    double max_n_cr;
    double max_g_min;
    double max_g_max;
    double max_p;

    double conv_length_in_SI, conv_dens_in_SI;
    double conv_Bfield_in_SI, conv_Vfield_in_SI;

    double total_gas_mass;
    double mu;

    dlist wl_list;

    int line_counter;
    uint char_counter;
    unsigned char ru[4];

    uint nrOfPlotPoints, nrOfPlotVectors;
    uint maxPlotLines;

    Vector3D ex, ey, ez;

    cell_basic ** cell_list;
    uint buffer_size;

    uint data_offset;
    uint dataID;
    uint data_len;

    uint nrOfDensRatios;
    uint nrOfOpiateIDs;

    uint nr_mixtures;
    uint nr_densities;
    uint size_gd_list;
    uint size_dd_list;
    uint multi_temperature_entries;
    uint stochastic_temperature_entries;
    uint * nr_dust_temp_sizes;
    uint * nr_stochastic_sizes;
    uint * nr_stochastic_temps;
    uint * size_skip;

    uint ** level_to_pos;
    uint *** line_to_pos;

    uilist data_pos_gd_list;
    uilist data_pos_dd_list;
    uilist data_pos_dt_list;
    uint data_pos_tg;
    uint data_pos_mx;
    uint data_pos_my;
    uint data_pos_mz;
    uint data_pos_vx;
    uint data_pos_vy;
    uint data_pos_vz;
    uint data_pos_px;
    uint data_pos_py;
    uint data_pos_pz;
    uilist data_pos_aalg_list;
    uint data_pos_amin;
    uint data_pos_amax;
    uint data_pos_size_param;
    uint data_pos_ra;
    uint data_pos_id;

    uint data_pos_vt;
    uint data_pos_pda;
    uint data_pos_op;

    uint data_pos_n_th;
    uint data_pos_T_e;
    uint data_pos_n_cr;
    uint data_pos_g_min;
    uint data_pos_g_max;
    uint data_pos_p;

    uint nr_rad_field_comp;

    uilist data_pos_rx_list;
    uilist data_pos_ry_list;
    uilist data_pos_rz_list;
    uilist data_pos_rf_list;

    uint data_pos_avg_th;
    uint data_pos_avg_dir;

    double turbulent_velocity;

    uslist data_ids;
    uint * pos_GasSpecRatios;
    uint * pos_OpiateIDS;

    double rot_angle1, rot_angle2;

    bool plt_gas_dens;
    bool plt_mol_dens;
    bool plt_dust_dens;
    bool plt_gas_temp;
    bool plt_dust_temp;
    bool plt_mag;
    bool plt_vel;
    bool plt_rat;
    bool plt_delta;
    bool plt_larm;
    bool plt_mach;
    bool plt_dust_id;
    bool plt_amin;
    bool plt_amax;
    bool plt_size_param;
    bool plt_rad_field1;
    bool plt_u_rad;
    bool plt_g_zero1;
    bool plt_n_th;
    bool plt_T_e;
    bool plt_n_cr;
    bool plt_g_min;
    bool plt_g_max;
    bool plt_p;

    bool plt_avg_dir;
    bool plt_avg_th;

    bool dust_is_mass_density, gas_is_mass_density;
    bool velocity_field_needed;
    bool spec_length_as_vector;

    double delta0;
    double larm_f;

    double total_volume;
    double cell_volume;

    double ** buffer_gas_dens;
    double ** buffer_mol_dens;
    double ** buffer_dust_dens;
    double * buffer_gas_temp;
    double ** buffer_dust_temp;
    double ** buffer_rat;
    double * buffer_delta;
    double * buffer_mag;
    double * buffer_mag_x;
    double * buffer_mag_y;
    double * buffer_mag_z;
    double * buffer_vel;
    double * buffer_vel_x;
    double * buffer_vel_y;
    double * buffer_vel_z;
    double * buffer_larm;
    double * buffer_mach;
    double * buffer_dust_mixture;
    double * buffer_dust_amin;
    double * buffer_dust_amax;
    double * buffer_dust_size_param;
    double *** buffer_rad_field;
    double * buffer_g_zero1;
    double * buffer_u_rad;
    double * buffer_n_th;
    double * buffer_T_e;
    double * buffer_n_cr;
    double * buffer_g_min;
    double * buffer_g_max;
    double * buffer_p;

    double * buffer_avg_th;
    double * buffer_avg_dir;

    double ** CextMeanTab;
    double ** CabsMeanTab;
    double ** CscaMeanTab;
    double * numberDensityTab;
    double * totalCellEmissionTab;
    uint chosen_wID;
    uint max_wavelengths;
};

#endif
