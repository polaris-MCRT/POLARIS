#pragma once
#include "Vector.h"
#include "chelper.h"

#ifndef CGRIDBASIC
#define CGRIDBASIC

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

        nrOfGnuPoints = 1000;
        nrOfGnuVectors = 1000;
        maxGridLines = 3;

        cell_list = 0;

        data_offset = 6;
        dataID = 0;
        data_len = 0;

        total_gas_mass = 0;
        mu = 0;

        nrOfDensRatios = 0;
        nrOfOpiateIDs = 0;

        nr_densities = 1;
        nr_mixtures = 0;

        multi_temperature_entries = 0;
        stochastic_entries = 0;
        nr_dust_sizes = 0;
        nr_stochastic_sizes = 0;
        nr_stochastic_temps = 0;

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
        data_pos_id = MAX_UINT;

        data_pos_vt = MAX_UINT;
        data_pos_pda = MAX_UINT;

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
        plt_amin=false;
        plt_amax=false;
        plt_rad_field = false;
        plt_g_zero = false;

        total_volume = 0;
        cell_volume = 0;

        buffer_gas_dens = 0;
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
        buffer_rad_field = 0;
        buffer_g_zero = 0;

        wl_list.resize(WL_STEPS);
        CMathFunctions::LogList(WL_MIN, WL_MAX, wl_list, 10); 
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
    }

    void printPhysicalParameter();

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

        max_pres = -1e300;
        min_pres = 1e300;

        max_vel = -1e300;
        min_vel = 1e300;

        max_mag = -1e300;
        min_mag = 1e300;

        max_value = 0;

        max_len = -1e300;
        min_len = 1e300;

        min_n_th=1e300;
        min_T_e=1e300;
        min_n_cr=1e300;
        min_g_min=1e300;
        min_g_max=1e300;
        min_p=1e300;

        max_n_th=-1e300;
        max_T_e=-1e300;
        max_n_cr=-1e300;
        max_g_min=-1e300;
        max_g_max=-1e300;
        max_p=-1e300;

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
        plt_amin=false;
        plt_amax=false;
        plt_rad_field = false;
        plt_g_zero = false;

        total_volume = 0;
        cell_volume = 0;

        buffer_gas_dens = 0;
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
        buffer_rad_field = 0;
        buffer_g_zero = 0;
    }

    double getTurbulentVelocity(cell_basic * cell)
    {
        if(turbulent_velocity > 0)
            return turbulent_velocity;

        if(data_pos_vt != MAX_UINT && turbulent_velocity < 0)
            return cell->getData(data_pos_vt);

        return 0;
    }

    double getTurbulentVelocity(photon_package * pp)
    {
        if(turbulent_velocity > 0)
            return turbulent_velocity;

        if(data_pos_vt != MAX_UINT && turbulent_velocity < 0)
            return pp->getPositionCell()->getData(data_pos_vt);

        return 0;
    }

    void updateDataRange(cell_basic * cell);

    void updateVelocity(cell_basic * cell, parameter & param)
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
        return(data_pos_vx != MAX_UINT && data_pos_vy != MAX_UINT &&
                data_pos_vz != MAX_UINT);

    }

    bool hasTurbulentVelocity()
    {
        return(data_pos_vt != MAX_UINT);
    }

    virtual Vector3D getCenter(cell_basic * cell)
    {
        return Vector3D(0, 0, 0);
    }

    uint getDataLength()
    {
        return data_len;
    }

    virtual double maxLength()
    {
        return 0;
    }

    ulong getMaxDataCells()
    {
        return max_cells;
    }

    uint getDataSize()
    {
        return max_data;
    }

    virtual void printParameter()
    {
    };

    virtual bool goToNextCellBorder(photon_package * pp)
    {
        return false;
    }

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

    void getMagFieldOrientation(photon_package * pp, Vector3D & mag_field, double & cos_theta, double & sin_theta,
            double & cos_2_phi, double & sin_2_phi)
    {
        // Get the magnetic field from grid
        mag_field = getMagField(pp);

        // Get the theta and phi angle from the magnetic field direction
        double theta = getThetaMag(pp);
        double phi = abs(getPhiMag(pp));

        // Calculate the sine and cosine including double angles
        cos_theta = cos(theta);
        sin_theta = sin(theta);
        cos_2_phi = cos(2.0 * phi);
        sin_2_phi = sin(2.0 * phi);
    }

    virtual bool findStartingPoint(photon_package * pp)
    {
        return false;
    }

    virtual void getLengths(uint bins, double & step_xy, double & off_xy)
    {
    }

    virtual bool positionPhotonInGrid(photon_package * pp)
    {
        return false;
    }

    double getMinLength()
    {
        return min_len;
    }

    double getMaxLength()
    {
        return max_len;
    }

    virtual bool getPolarRTGridParameter(double max_len, double pixel_width, uint max_subpixel_lvl,
        double * &_listR, uint &N_polar_r, uint * &N_polar_ph)
    {
        return false;
    }

    virtual bool createCellList()
    {
        return false;
    }

    cell_basic * getCellFromIndex(ulong i)
    {
        return cell_list[i];
    }

    void setSIConversionFactors(parameter & param)
    {
        mu = param.getMu();
        conv_length_in_SI = param.getSIConvLength();

        delta0 = param.getDelta0();
        larm_f = param.getLarmF();

        conv_dens_in_SI = abs(param.getSIConvDH());
        conv_Bfield_in_SI = param.getSIConvBField();
        conv_Vfield_in_SI = param.getSIConvVField();
        mass_fraction = param.getMassFraction();
    }

    //overloaded functions
    virtual bool findMatchingCell(photon_package * pp)
    {
        return false;
    }

    virtual bool writeGNUPlotFiles(string path, parameter & param)
    {
        return false;
    }

    void setDataSize(uint sz)
    {
        max_data = data_offset + sz;
    }

    void setDustInformation(uint _nr_mixtures, uint * _nr_dust_sizes,
            uint * _nr_stochastic_sizes, uint * _nr_stochastic_temps)
    {
        nr_mixtures = _nr_mixtures;
        nr_dust_sizes = _nr_dust_sizes;
        nr_stochastic_sizes = _nr_stochastic_sizes;
        nr_stochastic_temps = _nr_stochastic_temps;
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
        cell_basic * cell = pp->getPositionCell();
        setGasDensity(cell, dens);
    }

    void setGasDensity(photon_package * pp, uint i_density, double dens)
    {
        cell_basic * cell = pp->getPositionCell();
        setGasDensity(cell, i_density, dens);
    }

    virtual bool isInside(photon_package * pp, Vector3D & pos)
    {
        return false;
    }

    virtual bool isInside(Vector3D & pos)
    {
        return false;
    }

    void setSpecLengthAsVector(bool val)
    {
        spec_length_as_vector = val;
    }

    void updateSpecLength(photon_package * pp, double len)
    {
        cell_basic * cell = pp->getPositionCell();
        uint data_pos;
        if(spec_length_as_vector)
        {
            data_pos = data_offset + 4 * pp->getWavelengthID();
            Vector3D e_dir = len * rotateToCenter(pp, pp->getDirection(), false);
            cell->updateData(data_pos + 0, len);
            cell->updateData(data_pos + 1, e_dir.X());
            cell->updateData(data_pos + 2, e_dir.Y());
            cell->updateData(data_pos + 3, e_dir.Z());
        }
        else
        {
            data_pos = data_offset + pp->getWavelengthID();
            cell->updateData(data_pos, len);
        }
    }

    double getSpecLength(cell_basic * cell, uint wID)
    {
#ifdef CAMPS_BENCHMARK
        // To perform Camps et. al (2015) benchmark.
        double res = 0, wavelength = wl_list[wID], mult = 1.0e6;
        res = mult * CMathFunctions::mathis_isrf(wavelength);
        //res = 2.99e-14 * CMathFunctions::planck(wl_list[wID], 9000.0);
        return PIx4 * res * getVolume(cell);
#else
        if(spec_length_as_vector)
            return cell->getData(data_offset + 4 * wID + 0);
        else
            return cell->getData(data_offset + wID);
#endif
    }

    double getSpecLength(photon_package * pp, uint wID)
    {
        cell_basic * cell = pp->getPositionCell();
        return getSpecLength(cell, wID);
    }

    void getSpecLength(cell_basic * cell, uint wID, double & us, Vector3D & e_dir)
    {
        uint data_pos = data_offset + 4 * wID;

        us = cell->getData(data_pos + 0);
        e_dir.setX(cell->getData(data_pos + 1));
        e_dir.setY(cell->getData(data_pos + 2));
        e_dir.setZ(cell->getData(data_pos + 3));
    }

    void saveRadiationField()
    {
#pragma omp parallel for schedule(dynamic)
        for(long c_i = 0; c_i < long(max_cells); c_i++)
        {
            cell_basic * cell =  cell_list[c_i];
            double inv_vol = 1 / getVolume(cell);
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

    double getRadiationField(cell_basic * cell, uint wID)
    {
#ifdef CAMPS_BENCHMARK
        // To perform Camps et. al (2015) benchmark.
        double res = 0, wavelength = wl_list[wID], mult = 1e6;
        res = mult * CMathFunctions::mathis_isrf(wavelength);
        //res = 2.99e-14 * CMathFunctions::planck(wl_list[wID], 9000.0);
        return PIx4 * res;
#else
        // If the radiation field is needed after temp calculation, use the SpecLength instead
        if(data_pos_rf_list.empty())
            return getSpecLength(cell, wID) / getVolume(cell);
        return cell->getData(data_pos_rf_list[wID]);
#endif
    }

    double getRadiationField(photon_package * pp, uint wID)
    {
        cell_basic * cell = pp->getPositionCell();
        return getRadiationField(cell, wID);
    }

    void getRadiationField(photon_package * pp, uint w, double & us, Vector3D & e_dir)
    {
        // Init variables and get current cell
        Vector3D tmp_dir;
        cell_basic * cell = pp->getPositionCell();

        // Get radiation field strength and direction from cell
        us = cell->getData(data_pos_rf_list[w]);
        tmp_dir.setX(cell->getData(data_pos_rx_list[w]));
        tmp_dir.setY(cell->getData(data_pos_ry_list[w]));
        tmp_dir.setZ(cell->getData(data_pos_rz_list[w]));

        // Rotate vector to cell center
        e_dir = rotateToCenter(pp, tmp_dir, true);

        // Normalize the radiation field vector
        e_dir.normalize();
    }

    void getRadiationFieldInterp(photon_package * pp, double wavelength, double & us, Vector3D & e_dir)
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
        cell_basic * cell = pp->getPositionCell();

        // Get wavelength indices from radiation field calculation
        uint wID1 = CMathFunctions::biListIndexSearch(wavelength, wl_list);
        uint wID2 = wID1 + 1;

        // Interpolate radiation field strength and direction
        us = CMathFunctions::interpolate(wl_list[wID1], wl_list[wID2],
                cell->getData(data_pos_rf_list[wID1]), cell->getData(data_pos_rf_list[wID2]), wavelength);
        tmp_dir.setX(CMathFunctions::interpolate(wl_list[wID1], wl_list[wID2],
                cell->getData(data_pos_rx_list[wID1]), cell->getData(data_pos_rx_list[wID2]), wavelength));
        tmp_dir.setY(CMathFunctions::interpolate(wl_list[wID1], wl_list[wID2],
                cell->getData(data_pos_ry_list[wID1]), cell->getData(data_pos_ry_list[wID2]), wavelength));
        tmp_dir.setZ(CMathFunctions::interpolate(wl_list[wID1], wl_list[wID2],
                cell->getData(data_pos_rz_list[wID1]), cell->getData(data_pos_rz_list[wID2]), wavelength));

        // Rotate vector to cell center
        e_dir = rotateToCenter(pp, tmp_dir, true);

        // Normalize the radiation field vector
        e_dir.normalize();
    }

    double getGZero(cell_basic * cell)
    {
        // Init variables
        const double wl1 = 9.1165e-08, wl2 = 2.06640e-07;
        double g_zero;

        // If the radiation field is needed after temp calculation, use the SpecLength instead
        dlist rad_field(WL_STEPS);
        for(uint wID = 0; wID < WL_STEPS; wID++)
            rad_field[wID] = getRadiationField(cell, wID);

        g_zero = max(0.0, CMathFunctions::full_integ(wl_list, 
            rad_field, wl1, wl2) / 1.7836e-06);
        return g_zero;
    }

    double getGZero(photon_package * pp)
    {
        cell_basic * cell = pp->getPositionCell();
        return getGZero(cell);
    }

    double getMassFraction()
    {
        return mass_fraction;
    }

    double getMu()
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

    virtual Vector3D rotateToCenter(photon_package * pp, Vector3D dir, bool inv)
    {
        return dir;
    }

    void setDustTemperature(cell_basic * cell, uint i_density, uint a, double temp)
    {
        uint id = a + nr_densities;
        for(uint i = 0; i < i_density; i++)
            id += nr_dust_sizes[i];
        cell->setData(data_pos_dt_list[id], temp);
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
        cell_basic * cell = pp->getPositionCell();
        setDustTemperature(cell, i_density, a, temp);
    }

    void setDustTemperature(photon_package * pp, uint i_density, double temp)
    {
        cell_basic * cell = pp->getPositionCell();
        setDustTemperature(cell, i_density, temp);
    }

    void setDustTemperature(photon_package * pp, double temp)
    {
        cell_basic * cell = pp->getPositionCell();
        setDustTemperature(cell, temp);
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

    double getQBOffset(cell_basic * cell, uint i_density)
    {
        return cell->getData(data_pos_dt_list[i_density]);
    }

    double getQBOffset(photon_package * pp, uint i_density)
    {
        cell_basic * cell = pp->getPositionCell();
        return getQBOffset(cell, i_density);
    }

    double getQBOffset(cell_basic * cell, uint i_density, uint a)
    {
        uint id = a + nr_densities;
        for(uint i = 0; i < i_density; i++)
            id += nr_dust_sizes[i];
        return cell->getData(data_pos_dt_list[id]);
    }

    double getQBOffset(photon_package * pp, uint i_density, uint a)
    {
        cell_basic * cell = pp->getPositionCell();
        return getQBOffset(cell, i_density, a);
    }

    void setQBOffset(cell_basic * cell, uint i_density, uint a, double temp)
    {
        uint id = a + nr_densities;
        for(uint i = 0; i < i_density; i++)
            id += nr_dust_sizes[i];
        cell->setData(data_pos_dt_list[id], temp);
    }

    void setQBOffset(cell_basic * cell, uint i_density, double temp)
    {
        cell->setData(data_pos_dt_list[i_density], temp);
    }

    double getAlignedRadius(cell_basic *cell)
    {
        return cell->getData(data_pos_aalg);
    }

    double getAlignedRadius(photon_package * pp)
    {
        cell_basic * cell = pp->getPositionCell();
        return getAlignedRadius(cell);
    }

    void setAlignedRadius(cell_basic * cell, double _a_alg)
    {
        cell->setData(data_pos_aalg, _a_alg);
    }

    double getMinGrainRadius(cell_basic *cell)
    {
        if(data_pos_amin != MAX_UINT)
            return cell->getData(data_pos_amin);
        return 0;
    }

    double getMinGrainRadius(photon_package * pp)
    {
        cell_basic * cell = pp->getPositionCell();
        return getMinGrainRadius(cell);
    }

    double getMaxGrainRadius(cell_basic *cell)
    {
        if(data_pos_amax != MAX_UINT)
            return cell->getData(data_pos_amax);
        return 0;
    }

    double getMaxGrainRadius(photon_package * pp)
    {
        cell_basic * cell = pp->getPositionCell();
        return getMaxGrainRadius(cell);
    }

    uint getDustChoiceID(photon_package * pp)
    {
        if(data_pos_id != MAX_UINT)
            return uint(pp->getPositionCell()->getData(data_pos_id));
        else
            return 0;
    }

    uint getDustChoiceID(cell_basic * cell)
    {
        if(data_pos_id != MAX_UINT)
            return uint(cell->getData(data_pos_id));
        else
            return 0;
    }

    double getLvlPopLower(photon_package * pp, uint i_line)
    {
        return pp->getPositionCell()->getData(data_offset + 6 * i_line);
    }

    double getLvlPopUpper(photon_package * pp, uint i_line)
    {
        return pp->getPositionCell()->getData(data_offset + 6 * i_line + 1);
    }

    double getGaussA(photon_package * pp, uint i_line)
    {
        return pp->getPositionCell()->getData(data_offset + 6 * i_line + 2);
    }

    double getGamma(photon_package * pp, uint i_line)
    {
        return pp->getPositionCell()->getData(data_offset + 6 * i_line + 3);
    }

    double getDopplerWidth(photon_package * pp, uint i_line)
    {
        return pp->getPositionCell()->getData(data_offset + 6 * i_line + 4);
    }

    double getVoigtA(photon_package * pp, uint i_line)
    {
        return pp->getPositionCell()->getData(data_offset + 6 * i_line + 5);
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
        return cell->getID();
    }

    uint validateDataPositions(parameter & param);
    uint getDataIdsOffset(parameter & param);

    void setAlignedPos(uint pos)
    {
        data_pos_aalg = pos;
    }

    void setLvlPopLower(cell_basic *cell, uint i_line, double lvl_lower)
    {
        cell->setData(data_offset + 6 * i_line, lvl_lower);
    }

    void setLvlPopUpper(cell_basic *cell, uint i_line, double lvl_upper)
    {
        cell->setData(data_offset + 6 * i_line + 1, lvl_upper);
    }

    void setLineBroadening(cell_basic *cell, uint i_line,
            double gauss_a, double Gamma, double doppler_width, double voigt_a)
    {
        cell->setData(data_offset + 6 * i_line + 2, gauss_a);
        cell->setData(data_offset + 6 * i_line + 3, Gamma);
        cell->setData(data_offset + 6 * i_line + 4, doppler_width);
        cell->setData(data_offset + 6 * i_line + 5, voigt_a);
    }

    virtual double getPathLength(cell_basic * cell)
    {
        return 0;
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

    double getElectronTemperature(photon_package * pp)
    {
        if(data_pos_T_e!=MAX_UINT)
            return pp->getPositionCell()->getData(data_pos_T_e);

        return 0;
    }

    double getThermalElectronDensity(photon_package * pp)
    {
        if(data_pos_n_th!=MAX_UINT)
            return pp->getPositionCell()->getData(data_pos_n_th);

        return 0;
    }

    double getCRElectronDensity(photon_package * pp)
    {
        if(data_pos_n_cr!=MAX_UINT)
            return pp->getPositionCell()->getData(data_pos_n_cr);

        return 0;
    }

    double getGammaMin(photon_package * pp)
    {
        if(data_pos_g_min!=MAX_UINT)
            return pp->getPositionCell()->getData(data_pos_g_min);

        return 0;
    }

    double getGammaMax(photon_package * pp)
    {
        if(data_pos_g_max!=MAX_UINT)
            return pp->getPositionCell()->getData(data_pos_g_max);

        return 0;
    }

    double getPowerLawIndex(photon_package * pp)
    {
        if(data_pos_p!=MAX_UINT)
            return pp->getPositionCell()->getData(data_pos_p);

        return 0;
    }

    double getDustTemperature(cell_basic * cell, uint i_density, uint a)
    {
        if(data_pos_dt_list.size() > 0)
        {
            uint id = a + nr_densities;
            for(uint i = 0; i < i_density; i++)
                id += nr_dust_sizes[i];
            return cell->getData(data_pos_dt_list[id]);
        }
        else
            return 0;
    }

    double getDustTemperature(photon_package * pp, uint i_density, uint a)
    {
        cell_basic * cell = pp->getPositionCell();
        return getDustTemperature(cell, i_density, a);
    }

    double getDustTemperature(cell_basic * cell, uint i_density)
    {
        if(data_pos_dt_list.size() > 0)
            return cell->getData(data_pos_dt_list[i_density]);
        else
            return 0;
    }

    double getDustTemperature(photon_package * pp, uint i_density)
    {
        cell_basic * cell = pp->getPositionCell();
        return getDustTemperature(cell, i_density);
    }

    double getDustTemperature(cell_basic * cell)
    {
        double sum = 0;
        for(uint i_density = 0; i_density < nr_densities; i_density++)
            sum += getDustTemperature(cell, i_density) * getRelativeDustDensity(cell, i_density);
        return sum;
    }

    double getDustTemperature(photon_package * pp)
    {
        cell_basic * cell = pp->getPositionCell();
        return getDustTemperature(cell);
    }

    double getDustTempProbability(cell_basic * cell, uint i_density, uint a, uint t)
    {
        uint id = a * nr_stochastic_temps[i_density] + t;
        for(uint i = 0; i < i_density; i++)
            id += nr_stochastic_sizes[i] * nr_stochastic_temps[i];
        return cell->getData(data_offset + id);
    }

    double getDustTempProbability(photon_package * pp, uint i_density, uint a, uint t)
    {
        cell_basic * cell = pp->getPositionCell();
        return getDustTempProbability(cell, i_density, a, t);
    }

    double getPDAValue(cell_basic * cell)
    {
        return cell->getData(data_pos_pda);
    }

    double getGasTemperature(photon_package * pp)
    {
        cell_basic * cell = pp->getPositionCell();
        return getGasTemperature(cell);
    }

    double getGasTemperature(cell_basic * cell)
    {
        return cell->getData(data_pos_tg);
    }

    void setPlaneParameter(uint plane_index, double xy_step, double off_xy, double z_step, double off_z,
        int j, int k, int l, double & tx, double & ty, double & tz)
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
                    tx = 0.1;
                if(k != 0)
                {
                    double sg = CMathFunctions::sgn(k);
                    ty = double(k) * xy_step - sg * off_xy;
                }
                else
                    ty = 0.1;
                if(l != 0)
                {
                    double sg = CMathFunctions::sgn(l);
                    tz = double(l) * z_step - sg * off_z;
                }
                else
                    tz = 0.1;
                break;

            case PROJ_XZ:
                if(j != 0)
                {
                    double sg = CMathFunctions::sgn(j);
                    tx = double(j) * xy_step - sg * off_xy;
                }
                else
                    tx = 0.1;
                if(l != 0)
                {
                    double sg = CMathFunctions::sgn(l);
                    ty = double(l) * z_step - sg * off_z;
                }
                else
                    ty = 0.1;
                if(k != 0)
                {
                    double sg = CMathFunctions::sgn(k);
                    tz = double(k) * xy_step - sg * off_xy;
                }
                else
                    tz = 0.1;
                break;

            case PROJ_YZ:
                if(l != 0)
                {
                    double sg = CMathFunctions::sgn(l);
                    tx = double(l) * z_step - sg * off_z;
                }
                else
                    tx = 0.1;
                if(j != 0)
                {
                    double sg = CMathFunctions::sgn(j);
                    ty = double(j) * xy_step - sg * off_xy;
                }
                else
                    ty = 0.1;
                if(k != 0)
                {
                    double sg = CMathFunctions::sgn(k);
                    tz = double(k) * xy_step - sg * off_xy;
                }
                else
                    tz = 0.1;
                break;

            default:
                break;
        }
    }

    void fillMidplaneBuffer(double tx, double ty, double tz, uint i_cell)
    {

        photon_package * pp = new photon_package;
        pp->setPosition(Vector3D(tx, ty, tz));
        if(positionPhotonInGrid(pp))
        {
            uint id = 0;
            if(plt_gas_dens)
            {
                buffer_gas_dens[i_cell][0] = getGasDensity(pp);
                // Do it only once if only one gas distribution is defined
                if(nr_densities > 1 && data_pos_gd_list.size() >= nr_densities)
                    for(uint i_density = 0; i_density < nr_densities; i_density++)
                        buffer_gas_dens[i_cell][i_density + 1] = getGasDensity(pp, i_density);
            }
            if(plt_dust_dens)
            {
                buffer_dust_dens[i_cell][0] = getDustDensity(pp);
                // Do it only once if only one dust distribution is defined
                if(nr_densities > 1 && data_pos_dd_list.size() >= nr_densities)
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
                buffer_rat[i_cell] = getAlignedRadius(pp);
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
                double gas_temp = getGasTemperature(pp);
                double field = getMagField(pp).length();
                double Td = getDustTemperature(pp);
                double Tg = getGasTemperature(pp);
                double dens = getGasDensity(pp);
                double a_limit = CMathFunctions::calc_larm_limit(field, Td, Tg, dens, 0.5,4.1e-19);
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
            if(plt_rad_field)
                for(uint wID = 0; wID < WL_STEPS; wID++)
                    buffer_rad_field[i_cell][wID] = getRadiationField(pp, wID);
            if(plt_g_zero)
                buffer_g_zero[i_cell] = getGZero(pp);
        }
        else
        {
            if(plt_gas_dens)
            {
                buffer_gas_dens[i_cell][0] = 0;
                if(nr_densities > 1 && data_pos_gd_list.size() >= nr_densities)
                    for(uint i_density = 1; i_density <= nr_densities; i_density++)
                        buffer_gas_dens[i_cell][i_density] = 0;
            }
            if(plt_dust_dens)
            {
                buffer_dust_dens[i_cell][0] = 0;
                if(nr_densities > 1 && data_pos_dd_list.size() >= nr_densities)
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
                buffer_rat[i_cell] = 0;
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
            if(plt_rad_field)
                for(uint wID = 0; wID < WL_STEPS; wID++)
                    buffer_rad_field[i_cell][wID] = 0;
            if(plt_g_zero)
                buffer_g_zero[i_cell] = 0;
        }
        delete pp;
    }

    bool writeSpecialLines(string path);
    bool writeAMIRAFiles(string path, parameter & param, uint bins);
    bool writeMidplaneFits(string data_path, parameter & param, uint bins, bool all=false);
    bool writeMidplaneFits1(string data_path, parameter & param, uint bins, bool all=false);

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

    virtual bool saveBinaryGridFile(string filename)
    {
        return false;
    }

    virtual bool saveBinaryGridFile(string filename, ushort id, ushort data_size)
    {
        return false;
    }

    virtual bool loadGridFromBinrayFile(parameter & param)
    {
        return false;
    }

    virtual bool loadGridFromBinrayFile(parameter & param, uint data_len)
    {
        return false;
    }

    virtual double getVolume(cell_basic * cell)
    {
        return 0;
    }

    virtual double getVolume(photon_package * pp)
    {
        return 0;
    }

    virtual Vector3D getLowerBoundary(cell_basic * cell)
    {
        return Vector3D(0, 0, 0);
    }

    virtual Vector3D getLowerBoundary(photon_package * pp)
    {
        return Vector3D(0, 0, 0);
    }

    virtual Vector3D getUpperBoundary(cell_basic * cell)
    {
        return Vector3D(0, 0, 0);
    }

    virtual Vector3D getUpperBoundary(photon_package * pp)
    {
        return Vector3D(0, 0, 0);
    }

    double getGasDensity(cell_basic * cell)
    {
        double sum = 0;
        for(uint i_density = 0; i_density < data_pos_gd_list.size(); i_density++)
            sum += getGasDensity(cell, i_density);
        return sum;
    }

    double getGasDensity(cell_basic * cell, uint i_density)
    {
        if(data_pos_gd_list.size() > 0)
            return cell->getData(data_pos_gd_list[i_density]);
        else
            return 0;
    }

    double getGasDensity(photon_package * pp)
    {
        cell_basic * cell = pp->getPositionCell();
        return getGasDensity(cell);
    }

    double getGasDensity(photon_package * pp, uint i_density)
    {
        cell_basic * cell = pp->getPositionCell();
        return getGasDensity(cell, i_density);
    }

    double getGasNumberDensity(cell_basic * cell)
    {
        double sum = 0;
        for(uint i_density = 0; i_density < data_pos_gd_list.size(); i_density++)
            sum += cell->getData(data_pos_gd_list[i_density]);
        if(gas_is_mass_density)
            sum /= (mu * m_H);
        return sum;
    }

    double getGasNumberDensity(cell_basic * cell, uint i_density)
    {
        double dens = cell->getData(data_pos_gd_list[i_density]);
        if(gas_is_mass_density)
            dens /= (mu * m_H);
        return dens;
    }

    double getGasNumberDensity(photon_package * pp)
    {
        cell_basic * cell = pp->getPositionCell();
        return getGasNumberDensity(cell);
    }

    double getGasNumberDensity(photon_package * pp, uint i_density)
    {
        cell_basic * cell = pp->getPositionCell();
        return getGasNumberDensity(cell, i_density);
    }

    double getGasMassDensity(cell_basic * cell)
    {
        double sum = 0;
        for(uint i_density = 0; i_density < data_pos_gd_list.size(); i_density++)
            sum += cell->getData(data_pos_gd_list[i_density]);
        if(!gas_is_mass_density)
            sum *= (mu * m_H);
        return sum;
    }

    double getGasMassDensity(cell_basic * cell, uint i_density)
    {
        double dens = cell->getData(data_pos_gd_list[i_density]);
        if(!gas_is_mass_density)
            dens *= (mu * m_H);
        return dens;
    }

    double getGasMassDensity(photon_package * pp)
    {
        cell_basic * cell = pp->getPositionCell();
        return getGasMassDensity(cell);
    }

    double getGasMassDensity(photon_package * pp, uint i_density)
    {
        cell_basic * cell = pp->getPositionCell();
        return getGasMassDensity(cell, i_density);
    }

    bool useDustChoice()
    {
        if(data_pos_gd_list.size() > 1 || data_pos_dd_list.size() > 1)
            return false;
        return true;
    }

    bool useConstantGrainSizes()
    {
        if(data_pos_amin != MAX_UINT || data_pos_amax != MAX_UINT)
            return false;
        return true;
    }

    bool useDustDensities()
    {
        if(data_pos_dd_list.size() > 1)
            return true;
        return false;
    }

    double getDustDensity(cell_basic * cell)
    {
        double sum = 0;
        if(data_pos_dd_list.size() > 0)
        {
            for(uint i_density = 0; i_density < data_pos_dd_list.size(); i_density++)
                sum += cell->getData(data_pos_dd_list[i_density]);
            return sum;
        }
        else
            return 0;
    }

    double getDustDensity(cell_basic * cell, uint i_density)
    {
        if(data_pos_dd_list.size() > 0)
            return cell->getData(data_pos_dd_list[i_density]);
        else
            return 0;
    }

    double getDustDensity(photon_package * pp)
    {
        cell_basic * cell = pp->getPositionCell();
        return getDustDensity(cell);
    }

    double getDustDensity(photon_package * pp, uint i_density)
    {
        cell_basic * cell = pp->getPositionCell();
        return getDustDensity(cell, i_density);
    }

    double getRelativeDustDensity(cell_basic * cell, uint i_density)
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
        if(data_pos_dd_list.size() > 0)
        {
            double dust_dens = getDustDensity(cell, i_density);
            cell->setData(data_pos_dd_list[i_density], dust_dens * factor);
        }
        else
        {
            double gas_dens = getGasDensity(cell, i_density);
            cell->setData(data_pos_gd_list[i_density], gas_dens * factor);
        }
    }

    virtual bool positionPhotonInGridTest(photon_package * pp)
    {
        return false;
    };

    Vector3D getVelocityField(photon_package * pp)
    {
        return Vector3D(pp->getPositionCell()->getData(data_pos_vx),
            pp->getPositionCell()->getData(data_pos_vy),
            pp->getPositionCell()->getData(data_pos_vz));
    }

    Vector3D getVelocityField(cell_basic * cell)
    {
        return Vector3D(cell->getData(data_pos_vx),
                cell->getData(data_pos_vy),
                cell->getData(data_pos_vz));
    }

    virtual void setCrossSections(cross_sections & cs)
    {

    }

    virtual void setID(uint id)
    {

    }

    double getCellAbundance(photon_package * pp, uint id)
    {
        cell_basic * cell = pp->getPositionCell();
        return getCellAbundance(cell, id);
    }

    double getCellAbundance(cell_basic * cell, uint id)
    {
        if(id > nrOfDensRatios - 1)
            return 0;

        uint pos = pos_GasSpecRatios[id];
        return cell->getData(pos);
    }

    double getOpiateIDParameter(cell_basic * cell, uint id)
    {
        if(id > nrOfOpiateIDs - 1)
            return 0;

        uint pos = pos_OpiateIDS[id];

        return cell->getData(pos);
    }

    Vector3D getMagField(cell_basic * cell)
    {
        return Vector3D(cell->getData(data_pos_mx), cell->getData(data_pos_my), cell->getData(data_pos_mz));
    }

    Vector3D getMagField(photon_package * pp)
    {
        return Vector3D(pp->getPositionCell()->getData(data_pos_mx),
            pp->getPositionCell()->getData(data_pos_my),
            pp->getPositionCell()->getData(data_pos_mz));
    }

    void setMagField(cell_basic * cell, const Vector3D & mag)
    {
        cell->setData(data_pos_mx, mag.X());
        cell->setData(data_pos_my, mag.Y());
        cell->setData(data_pos_mz, mag.Z());
    }

    virtual bool next(photon_package * pp)
    {
        return false;
    }

    double getThetaMag(photon_package * pp)
    {
        return getAngleTheta(pp->getDirection(), getMagField(pp));
    }

    double getPhiMag(photon_package * pp)
    {
        return getAnglePhi(pp->getEX(), pp->getEY(), getMagField(pp)) - PI2;
    }

    double getTheta(cell_basic * cell, Vector3D & dir)
    {
        return getAngleTheta(dir, getMagField(cell));
    }

    double getThetaPhoton(photon_package * pp, Vector3D & dir)
    {
        return getAngleTheta(pp->getDirection(), dir);
    }

    bool getDustIsMassDensity()
    {
        return dust_is_mass_density;
    }

    bool getGasIsMassDensity()
    {
        return gas_is_mass_density;
    }

    bool getRadiationFieldAvailable()
    {
        if(data_pos_rx_list.empty() || data_pos_ry_list.empty() ||
                data_pos_rz_list.empty() || data_pos_rf_list.empty())
            return false;
        return true;
    }

    double getTotalGasMass()
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

    bool doPDA(parameter & param, uint pda_id);

    uint getTemperatureFieldInformation()
    {
        // Check which kind of temperature calculation the grid supports
        if(data_pos_dt_list.size() == multi_temperature_entries)
            return TEMP_FULL;
        else if(data_pos_dt_list.size() == nr_densities)
            return TEMP_EFF;
        else if(data_pos_dt_list.size() == 1)
            return TEMP_SINGLE;
        else if(data_pos_dt_list.size() == 0)
            return TEMP_EMPTY;
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
                    if(data_pos_gd_list.size() > 0)
                    {
                        if(data_pos_id != MAX_UINT)
                        {
                            cout << "ERROR: Multiple densities and dust choices cannot be combined!" << endl;
                            return false;
                        }
                        if(gas_is_mass_density == true)
                        {
                            cout << "ERROR: Gas number densities cannot be combined with gas mass densities!" << endl;
                            return false;
                        }
                    }
                    data_pos_gd_list.push_back(i);
                    gas_is_mass_density = false;
                    break;

                case GRIDgas_mdens:
                    if(data_pos_gd_list.size() > 0)
                    {
                        if(data_pos_id != MAX_UINT)
                        {
                            cout << "ERROR: Multiple densities and dust choices cannot be combined!" << endl;
                            return false;
                        }
                        if(gas_is_mass_density == false)
                        {
                            cout << "ERROR: Gas mass densities cannot be combined with gas number densities!" << endl;
                            return false;
                        }
                    }
                    data_pos_gd_list.push_back(i);
                    gas_is_mass_density = true;
                    break;

                case GRIDdust_dens:
                    if(data_pos_dd_list.size() > 0)
                    {
                        if(data_pos_id != MAX_UINT)
                        {
                            cout << "ERROR: Multiple densities and dust choices cannot be combined!" << endl;
                            return false;
                        }
                        if(dust_is_mass_density == true)
                        {
                            cout << "ERROR: Dust number densities cannot be combined with dust mass densities!" << endl;
                            return false;
                        }
                    }
                    data_pos_dd_list.push_back(i);
                    dust_is_mass_density = false;
                    break;

                case GRIDdust_mdens:
                    if(data_pos_dd_list.size() > 0)
                    {
                        if(data_pos_id != MAX_UINT)
                        {
                            cout << "ERROR: Multiple densities and dust choices cannot be combined!" << endl;
                            return false;
                        }
                        if(dust_is_mass_density == false)
                        {
                            cout << "ERROR: Dust mass densities cannot be combined with dust number densities!" << endl;
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
                        cout << "ERROR: Grid ID " << GRIDgas_temp << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_tg = i;
                    break;

                case GRIDmx:
                    if(data_pos_mx != MAX_UINT)
                    {
                        cout << "ERROR: Grid ID " << GRIDmx << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_mx = i;
                    break;

                case GRIDmy:
                    if(data_pos_my != MAX_UINT)
                    {
                        cout << "ERROR: Grid ID " << GRIDmy << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_my = i;
                    break;

                case GRIDmz:
                    if(data_pos_mz != MAX_UINT)
                    {
                        cout << "ERROR: Grid ID " << GRIDmz << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_mz = i;
                    break;

                case GRIDvx:
                    if(data_pos_vx != MAX_UINT)
                    {
                        cout << "ERROR: Grid ID " << GRIDvx << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_vx = i;
                    break;

                case GRIDvy:
                    if(data_pos_vy != MAX_UINT)
                    {
                        cout << "ERROR: Grid ID " << GRIDvy << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_vy = i;
                    break;

                case GRIDvz:
                    if(data_pos_vz != MAX_UINT)
                    {
                        cout << "ERROR: Grid ID " << GRIDvz << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_vz = i;
                    break;

                case GRIDpx:
                    if(data_pos_px != MAX_UINT)
                    {
                        cout << "ERROR: Grid ID " << GRIDpx << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_px = i;
                    break;

                case GRIDpy:
                    if(data_pos_py != MAX_UINT)
                    {
                        cout << "ERROR: Grid ID " << GRIDpy << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_py = i;
                    break;

                case GRIDpz:
                    if(data_pos_pz != MAX_UINT)
                    {
                        cout << "ERROR: Grid ID " << GRIDpz << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_pz = i;
                    break;

                case GRIDa_alg:
                    if(data_pos_aalg != MAX_UINT)
                    {
                        cout << "ERROR: Grid ID " << GRIDa_alg << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_aalg = i;
                    break;

                case GRIDa_min:
                    if(data_pos_amin != MAX_UINT)
                    {
                        cout << "ERROR: Grid ID " << GRIDa_min << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_amin = i;
                    break;

                case GRIDa_max:
                    if(data_pos_amax != MAX_UINT)
                    {
                        cout << "ERROR: Grid ID " << GRIDa_max << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_amax = i;
                    break;

                case GRIDq:
                    if(data_pos_eq != MAX_UINT)
                    {
                        cout << "ERROR: Grid ID " << GRIDq << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_eq = i;
                    break;

                case GRIDv_turb:
                    if(data_pos_vt != MAX_UINT)
                    {
                        cout << "ERROR: Grid ID " << GRIDvx << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_vt = i;
                    break;

                case GRIDn_th:
                    if(data_pos_n_th != MAX_UINT)
                    {
                        cout << "ERROR: Grid ID " << GRIDn_th << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_n_th = i;
                    break;

                case GRIDT_e:
                    if(data_pos_T_e != MAX_UINT)
                    {
                        cout << "ERROR: Grid ID " << GRIDT_e << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_T_e = i;
                    break;

                case GRIDn_cr:
                    if(data_pos_n_cr != MAX_UINT)
                    {
                        cout << "ERROR: Grid ID " << GRIDn_cr << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_n_cr = i;
                    break;

                case GRIDg_min:
                    if(data_pos_g_min != MAX_UINT)
                    {
                        cout << "ERROR: Grid ID " << GRIDg_min << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_g_min = i;
                    break;

                 case GRIDg_max:
                    if(data_pos_g_max != MAX_UINT)
                    {
                        cout << "ERROR: Grid ID " << GRIDg_max << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_g_max = i;
                    break;

                case GRIDp:
                    if(data_pos_p != MAX_UINT)
                    {
                        cout << "ERROR: Grid ID " << GRIDp << " can be set only once!" << endl;
                        return false;
                    }

                    data_pos_p = i;
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
                        cout << "ERROR: Multiple densities and dust choices cannot be combined!" << endl;
                        return false;
                    }
                    if(data_pos_id != MAX_UINT)
                    {
                        cout << "ERROR: Grid ID " << GRIDdust_id << " can be set only once!" << endl;
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
                    cout << "ERROR: Unknown data IDs!" << endl;
                    cout << "  IDs have to be between " << minGRID << " and " << maxGRID << "!" << endl;
                    return false;
            }
        }

        if(data_pos_gd_list.size() == 0)
        {
            cout << "ERROR: Grid requires a gas density! " << endl;
            return false;
        }

        pos_GasSpecRatios = new uint[nrOfDensRatios];
        uint pos_counter = 0;

        for(uint i = 0; i < data_offset; i++)
        {
            if(data_ids[i] == GRIDratio)
            {
                pos_GasSpecRatios[pos_counter] = i;
                //cout << pos_counter << "\t" << pos_ration[pos_counter] << endl;
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
                //cout << pos_counter << "\t" << pos_ration[pos_counter] << endl;
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
                    cout << "ERROR: A grid ID of 0 requires a data length of 6!" << endl;
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
                    cout << "ERROR: A grid ID of 1 requires a data length of 7!" << endl;
                    cout << "       Cannot create a compatible octree!" << endl;
                    return false;
                }

                data_pos_gd_list.push_back(0);
                data_pos_dt_list.push_back(1);
                data_pos_tg = 2;

                data_pos_mx = 3;
                data_pos_my = 4;
                data_pos_mz = 5;

                data_pos_aalg = 6;

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
                    cout << "ERROR: A grid ID of 6 requires a data length of 9!" << endl;
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
                    cout << "ERROR: A grid ID of 7 requires a data length of 10!" << endl;
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

                data_pos_aalg = 9;

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
                cout << "ERROR: Grid ID = " << dataID << " doesn't match any of the old octree IDs!" << endl;
                cout << "       Cannot create a compatible octree!" << endl;
                return false;
        }

        dataID = GRID_ID_OCT;

        return true;
    }

    uint SynchrotronCheck(parameter & param)
    {
        if(data_pos_n_th == MAX_UINT && data_pos_n_th == MAX_UINT)
        {
            cout << "ERROR: Neither thermal electrons nor CR electrons are defined in grid file!" << endl;
            cout << "       No SYNCHROTRON calculation possible." << endl;
            return MAX_UINT;
        }

        if(data_pos_mx == MAX_UINT)
        {
            cout << "ERROR: Grid contains no magnetic Bx component!" << endl;
            cout << "       No SYNCHROTRON calculation possible." << endl;
            return MAX_UINT;
        }
        if(data_pos_my == MAX_UINT)
        {
            cout << "ERROR: Grid contains no magnetic By component!" << endl;
            cout << "       No SYNCHROTRON calculation possible." << endl;
            return MAX_UINT;
        }
        if(data_pos_mz == MAX_UINT)
        {
            cout << "ERROR: Grid contains no magnetic Bz component!" << endl;
            cout << "       No SYNCHROTRON calculation possible." << endl;
            return MAX_UINT;
        }

        if(data_pos_n_th == MAX_UINT)
        {
            cout << "WARNING: Grid contains no thermal electron component!" << endl;
            cout << "         Only CR SYNCHROTRON calculation possible." << endl;
        }
        else
        {
            if(data_pos_tg == MAX_UINT)
            {
                cout << "ERROR: Grid contains no gas temperature!" << endl;
                cout << "       No SYNCHROTRON calculation possible." << endl;
                return MAX_UINT;
            }

            if(data_pos_T_e == MAX_UINT)
            {
                cout << "WARNING: Grid contains no electron temperature component!" << endl;
                cout << "         Electron temperature will be set to equal gas temperature!" << endl;
                data_pos_T_e = data_pos_tg;
            }
        }

        if(data_pos_n_cr == MAX_UINT)
        {
            cout << "ERROR: Grid contains no CR electron component!" << endl;
            cout << "       No SYNCHROTRON calculation possible." << endl;
            return MAX_UINT;
        }

        if(data_pos_g_min == MAX_UINT)
        {
            cout << "ERROR: Grid contains no gamma_min component!" << endl;
            cout << "       No SYNCHROTRON calculation possible." << endl;
            return MAX_UINT;
        }

        if(data_pos_g_min == MAX_UINT)
        {
            cout << "ERROR: Grid contains no gamma_min component!" << endl;
            cout << "       No SYNCHROTRON calculation possible." << endl;
            return MAX_UINT;
        }

        if(data_pos_g_max == MAX_UINT)
        {
            cout << "ERROR: Grid contains no gamma_max component!" << endl;
            cout << "       No SYNCHROTRON calculation possible." << endl;
            return MAX_UINT;
        }

        if(data_pos_p == MAX_UINT)
        {
            cout << "ERROR: Grid contains no electron power-law index p component!" << endl;
            cout << "       No SYNCHROTRON calculation possible." << endl;
            return MAX_UINT;
        }
        return 0;
    }

    uint OpiateCheck(parameter & param)
    {
        if(data_pos_tg == MAX_UINT)
        {
            cout << "ERROR: Grid contains no gas temperature!" << endl;
            cout << "       No OPIATE calculation possible." << endl;
            return MAX_UINT;
        }
        return 0;
    }

    uint TempCheck(parameter & param, uint & tmp_data_offset)
    {
        uint extra_temp_entries = 0;
        if(getTemperatureFieldInformation() == TEMP_EMPTY)
        {
            if(param.getDustTempMulti())
                extra_temp_entries = multi_temperature_entries;
            else
                extra_temp_entries = nr_densities;
        }
        else if(getTemperatureFieldInformation() == TEMP_SINGLE)
        {
            if(param.getDustTempMulti())
                extra_temp_entries = multi_temperature_entries - 1;
            else
                extra_temp_entries = nr_densities - 1;
        }
        else if(getTemperatureFieldInformation() == TEMP_EFF)
        {
            if(param.getDustTempMulti())
                extra_temp_entries = multi_temperature_entries - nr_densities;
            else
                extra_temp_entries = 0;
        }
        else if(getTemperatureFieldInformation() == MAX_UINT)
        {
            cout << "ERROR: The grid does not include the correct information for temperature calculations" << endl;
            cout << "       No dust temperature calculation possible (full_dust_temp or stochastic heating?)." << endl;
            return MAX_UINT;
        }
        for(uint i_entries = 0; i_entries < extra_temp_entries; i_entries++)
        {
            data_pos_dt_list.push_back(data_offset + tmp_data_offset);
            data_ids.push_back(GRIDdust_temp);
            tmp_data_offset++;
        }

        if(param.getSaveRadiationField())
            if(data_pos_rx_list.size() != 0 || data_pos_ry_list.size() != 0 ||
                    data_pos_rz_list.size() != 0 || data_pos_rf_list.size() != 0)
            {
                cout << "ERROR: The grid includes partial/broken information about a radiation field!" << endl;
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
                //cout << "Create entries for gas temperature   : done" << endl;
            }
            else
            {
                param.setAdjTgas(1.0);
                data_pos_tg = data_offset + tmp_data_offset;
                data_ids.push_back(GRIDgas_temp);
                tmp_data_offset++;
                cout << SEP_LINE;
                cout << "HINT: No gas temperature found in grid." << endl;
                cout << "    Add entry and set gas temperature to dust temperature after calculation!" << endl;
                cout << SEP_LINE;
            }
        }
        return 0;
    }

    uint RatCheck(parameter & param, uint & tmp_data_offset)
    {
        if(data_pos_aalg == MAX_UINT)
        {
            data_pos_aalg = data_offset + tmp_data_offset;
            data_ids.push_back(GRIDa_alg);
            tmp_data_offset++;
            // cout << "Create entries for RAT alignment     : done" << endl;
        }

        if(getTemperatureFieldInformation() == TEMP_EMPTY)
        {
            cout << "ERROR: Grid contains no dust temperature!" << endl;
            cout << "       No RAT calculation possible." << endl;
            return MAX_UINT;
        }
        else if(getTemperatureFieldInformation() == MAX_UINT || getTemperatureFieldInformation() == TEMP_SINGLE)
        {
            cout << "ERROR: The grid does not include the information for temperature calculations" << endl;
            cout << "       No RAT calculation possible." << endl;
            return MAX_UINT;
        }

        if(data_pos_tg == MAX_UINT)
        {
            cout << "ERROR: Grid contains no gas temperature!" << endl;
            cout << "       No RAT calculation possible." << endl;
            return MAX_UINT;
        }
        if(data_pos_mx == MAX_UINT)
        {
            cout << "WARNING: Grid contains no magnetic Bx component!" << endl;
            cout << "         No follow up calculations possible." << endl;
        }
        if(data_pos_my == MAX_UINT)
        {
            cout << "WARNING: Grid contains no magnetic By component!" << endl;
            cout << "         No follow up calculations possible." << endl;
        }
        if(data_pos_mz == MAX_UINT)
        {
            cout << "WARNING: Grid contains no magnetic Bz component!" << endl;
            cout << "         No follow up calculations possible." << endl;
        }
        return 0;
    }

    uint DustEmissionCheck(parameter & param)
    {
        if(getTemperatureFieldInformation() == TEMP_EMPTY)
        {
            cout << "ERROR: Grid contains no dust temperature!" << endl;
            cout << "       No dust emission possible." << endl;
            return MAX_UINT;
        }
        else if(getTemperatureFieldInformation() == MAX_UINT || getTemperatureFieldInformation() == TEMP_SINGLE)
        {
            cout << "ERROR: The grid does not include the information for temperature calculations" << endl;
            cout << "       No dust emission possible." << endl;
            return MAX_UINT;
        }

        if(param.getStochasticHeatingMaxSize())
        {
            if(data_pos_rf_list.size() != WL_STEPS)
            {
               cout << "ERROR: The grid includes partial/no information about a radiation field!" << endl;
               cout << "       No dust emission with stochastic heating possible." << endl;
               return MAX_UINT;
           }
        }

        if(data_pos_rf_list.size() > 0 && data_pos_rf_list.size() != WL_STEPS)
        {
           cout << "ERROR: The grid includes partial/no information about a radiation field!" << endl;
           cout << "       No dust emission possible." << endl;
           return MAX_UINT;
       }

        if(param.getAlign() != 0 && !param.getAligPA())
        {
            if(data_pos_tg == MAX_UINT)
            {
                cout << "ERROR: Grid contains no gas temperature!" << endl;
                cout << "       No dust emission with aligned dust grains possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_mx == MAX_UINT)
            {
                cout << "ERROR: Grid contains no magnetic Bx component!" << endl;
                cout << "       No dust emission with aligned dust grains possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_my == MAX_UINT)
            {
                cout << "ERROR: Grid contains no magnetic By component!" << endl;
                cout << "       No dust emission with aligned dust grains possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_mz == MAX_UINT)
            {
                cout << "ERROR: Grid contains no magnetic Bz component!" << endl;
                cout << "       No dust emission with aligned dust grains possible." << endl;
                return MAX_UINT;
            }
        }

        if(param.getAligGOLD())
        {
            if(data_pos_vx == MAX_UINT)
            {
                cout << "ERROR: Grid contains no velocity vx component!" << endl;
                cout << "        No dust emission with GOLD alignment possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_vy == MAX_UINT)
            {
                cout << "ERROR: Grid contains no velocity vy component!" << endl;
                cout << "        No dust emission with GOLD alignment possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_vz == MAX_UINT)
            {
                cout << "ERROR: Grid contains no velocity vz component!" << endl;
                cout << "        No dust emission with GOLD alignment possible." << endl;
                return MAX_UINT;
            }
        }

        if(param.getAligRAT())
        {
            if(data_pos_aalg == MAX_UINT)
            {
                cout << "ERROR: Grid contains no minimum alignment radius for RATs!" << endl;
                cout << "        No dust emission with RAT alignment possible." << endl;
                return MAX_UINT;
            }
        }
        return 0;
    }

    uint DustScatteringCheck(parameter & param)
    {
        if(param.getAlign() != 0 && !param.getAligPA())
        {
            if(data_pos_tg == MAX_UINT)
            {
                cout << "ERROR: Grid contains no gas temperature!           " << endl;
                cout << "       No dust scattering calculations with aligned dust grains possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_mx == MAX_UINT)
            {
                cout << "ERROR: Grid contains no magnetic Bx component!     " << endl;
                cout << "       No dust scattering calculations with aligned dust grains possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_my == MAX_UINT)
            {
                cout << "ERROR: Grid contains no magnetic By component!     " << endl;
                cout << "       No dust scattering calculations with aligned dust grains possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_mz == MAX_UINT)
            {
                cout << "ERROR: Grid contains no magnetic Bz component!     " << endl;
                cout << "       No dust scattering calculations with aligned dust grains possible." << endl;
                return MAX_UINT;
            }
        }

        if(param.getAligGOLD())
        {
            if(data_pos_vx == MAX_UINT)
            {
                cout << "ERROR: Grid contains no velocity vx component!  " << endl;
                cout << "        No dust scattering calculations with GOLD alignment possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_vy == MAX_UINT)
            {
                cout << "ERROR: Grid contains no velocity vy component!  " << endl;
                cout << "        No dust scattering calculations with GOLD alignment possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_vz == MAX_UINT)
            {
                cout << "ERROR: Grid contains no velocity vz component!  " << endl;
                cout << "        No dust scattering calculations with GOLD alignment possible." << endl;
                return MAX_UINT;
            }
        }

        if(param.getAligRAT())
        {
            if(data_pos_aalg == MAX_UINT)
            {
                cout << "ERROR: Grid contains no minimum alignment radius for RATs!" << endl;
                cout << "        No dust scattering calculations with RAT alignment possible." << endl;
                return MAX_UINT;
            }
        }

        return 0;
    }

    uint RadiationForceCheck(parameter & param)
    {
        if(getTemperatureFieldInformation() == TEMP_EMPTY)
        {
            cout << "ERROR: Grid contains no dust temperature!" << endl;
            cout << "       No FORCE calculation possible." << endl;
            return MAX_UINT;
        }
        else if(getTemperatureFieldInformation() == MAX_UINT || getTemperatureFieldInformation() == TEMP_SINGLE)
        {
            cout << "ERROR: The grid does not include the information for temperature calculations" << endl;
            cout << "       No FORCE calculation possible." << endl;
            return MAX_UINT;
        }

        if(data_pos_rf_list.size() > 0 && data_pos_rf_list.size() != WL_STEPS)
        {
           cout << "ERROR: The grid includes partial/no information about a radiation field!" << endl;
           cout << "       No FORCE calculation possible." << endl;
           return MAX_UINT;
       }

        return 0;
    }

    uint LineEmissionCheck(parameter & param)
    {
        if(param.getTotalNrOfDustComponents() != 0)
        {
            if(getTemperatureFieldInformation() == TEMP_EMPTY)
            {
                cout << "ERROR: Grid contains no dust temperature!" << endl;
                cout << "       No line transfer including dust emission possible." << endl;
                return MAX_UINT;
            }
            else if(getTemperatureFieldInformation() == MAX_UINT || getTemperatureFieldInformation() == TEMP_SINGLE)
            {
                cout << "ERROR: The grid does not include the information for temperature calculations" << endl;
                cout << "       No line transfer including dust emission possible." << endl;
                return MAX_UINT;
            }
        }

        if(data_pos_tg == MAX_UINT)
        {
            cout << "ERROR: Grid contains no gas temperature!" << endl;
            cout << "       No line transfer with possible.  " << endl;
            return MAX_UINT;
        }

        if(velocity_field_needed)
        {
            if(data_pos_mx == MAX_UINT)
            {
                cout << "ERROR: Grid contains no magnetic Bx component!" << endl;
                cout << "       No line transfer possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_my == MAX_UINT)
            {
                cout << "ERROR: Grid contains no magnetic By component!" << endl;
                cout << "       No line transfer possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_mz == MAX_UINT)
            {
                cout << "ERROR: Grid contains no magnetic Bz component!" << endl;
                cout << "       No line transfer possible." << endl;
                return MAX_UINT;
            }
        }

        if(param.getKeplerStarMass() == 0)
        {
            if(data_pos_vx == MAX_UINT)
            {
                cout << "ERROR: Grid contains no velocity vx component!" << endl;
                cout << "        No line transfer possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_vy == MAX_UINT)
            {
                cout << "ERROR: Grid contains no velocity vy component!" << endl;
                cout << "        No line transfer possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_vz == MAX_UINT)
            {
                cout << "ERROR: Grid contains no velocity vz component!" << endl;
                cout << "        No line transfer possible." << endl;
                return MAX_UINT;
            }
        }
        return 0;
    }

    uint ProbingCheck(parameter & param)
    {
        if(getTemperatureFieldInformation() == TEMP_EMPTY)
        {
            cout << "ERROR: Grid contains no dust temperature!" << endl;
            cout << "       No LOS analysis with aligned dust grains possible." << endl;
            return MAX_UINT;
        }
        else if(getTemperatureFieldInformation() == MAX_UINT || getTemperatureFieldInformation() == TEMP_SINGLE)
        {
            cout << "ERROR: The grid does not include the information for temperature calculations" << endl;
            cout << "       No LOS analysis with aligned dust grains possible." << endl;
            return MAX_UINT;
        }

        if(data_pos_rf_list.size() > 0 && data_pos_rf_list.size() != WL_STEPS)
        {
           cout << "ERROR: The grid includes partial/no information about a radiation field!" << endl;
           cout << "       No LOS analysis with aligned dust grains possible." << endl;
           return MAX_UINT;
       }

        if(param.getAlign() != 0 && !param.getAligPA())
        {
            if(data_pos_tg == MAX_UINT)
            {
                cout << "ERROR: Grid contains no gas temperature!" << endl;
                cout << "       No LOS analysis with aligned dust grains possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_mx == MAX_UINT)
            {
                cout << "ERROR: Grid contains no magnetic Bx component!" << endl;
                cout << "       No LOS analysis with aligned dust grains possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_my == MAX_UINT)
            {
                cout << "ERROR: Grid contains no magnetic By component!" << endl;
                cout << "       No LOS analysis with aligned dust grains possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_mz == MAX_UINT)
            {
                cout << "ERROR: Grid contains no magnetic Bz component!" << endl;
                cout << "       No LOS analysis with aligned dust grains possible." << endl;
                return MAX_UINT;
            }
        }

        if(param.getAligGOLD())
        {
            if(data_pos_vx == MAX_UINT)
            {
                cout << "ERROR: Grid contains no velocity vx component!" << endl;
                cout << "        No LOS analysis with GOLD alignment possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_vy == MAX_UINT)
            {
                cout << "ERROR: Grid contains no velocity vy component!" << endl;
                cout << "        No LOS analysis with GOLD alignment possible." << endl;
                return MAX_UINT;
            }
            if(data_pos_vz == MAX_UINT)
            {
                cout << "ERROR: Grid contains no velocity vz component!" << endl;
                cout << "        No LOS analysis with GOLD alignment possible." << endl;
                return MAX_UINT;
            }
        }
        return 0;
    }

protected:
    //uint grid_type;
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
    double mass_fraction;

    double total_gas_mass;
    double mu;

    dlist wl_list;

    int line_counter;
    uint char_counter;
    unsigned char ru[4];

    uint nrOfGnuPoints, nrOfGnuVectors;
    uint maxGridLines;

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
    uint multi_temperature_entries;
    uint stochastic_entries;
    uint * nr_dust_sizes;
    uint * nr_stochastic_sizes;
    uint * nr_stochastic_temps;

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
    uint data_pos_aalg;
    uint data_pos_amin;
    uint data_pos_amax;
    uint data_pos_eq;
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

    uilist data_pos_rx_list;
    uilist data_pos_ry_list;
    uilist data_pos_rz_list;
    uilist data_pos_rf_list;

    double turbulent_velocity;

    uslist data_ids;
    uint * pos_GasSpecRatios;
    uint * pos_OpiateIDS;

    double rot_angle1, rot_angle2;

    bool plt_gas_dens;
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
    bool plt_rad_field;
    bool plt_g_zero;

    bool dust_is_mass_density, gas_is_mass_density;
    bool velocity_field_needed;
    bool spec_length_as_vector;

    double delta0;
    double larm_f;

    double total_volume;
    double cell_volume;

    double ** buffer_gas_dens;
    double ** buffer_dust_dens;
    double * buffer_gas_temp;
    double ** buffer_dust_temp;
    double * buffer_rat;
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
    double ** buffer_rad_field;
    double * buffer_g_zero;
};

#endif
