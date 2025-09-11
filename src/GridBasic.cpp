/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include <valarray>
#include "CCfits/FITS.h"
#include "CCfits/FITSUtilT.h"
#include "CCfits/FitsError.h"
#include "CCfits/KeyData.h"
#include "CCfits/PHDU.h"
#include "CCfits/PHDUT.h"
#include "fitsio.h"
#include "GridBasic.hpp"

void CGridBasic::resetGridValues()
{
    max_cells = 0;
    max_value = 0;
    max_data = 0;

    delta0 = 8.28e23 * 2.5e-12 * 1e8 * 1e-6 * 1e6;
    // This would always reset the value to the value given below. I.e. the command "<larm_f> value" in the command file would effectively be ignored, hence it needs to be commented out
    // larm_f = 4.1e-21;

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

    min_dust_aalg = -1e300;
    max_dust_aalg =  1e300;
    
    min_dust_akrat = -1e300;
    max_dust_akrat =  1e300;

    min_dust_amin = -1e300;
    max_dust_amin =  1e300;
    
    min_dust_amax = -1e300;
    max_dust_amax =  1e300;

    min_dust_size_param = -1e300;
    max_dust_size_param =  1e300;
    
    min_ion_dens = -1e300;
    max_ion_dens =  1e300;
    
    min_ion_Z = -1e300;
    max_ion_Z =  1e300;

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
    //data_pos_amin = MAX_UINT;
    //data_pos_amax = MAX_UINT;
    //data_pos_size_param = MAX_UINT;
    
    
    //data_pos_ra = MAX_UINT;
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
    //plt_delta = false;
    plt_larm = false;
    //plt_mach = false;
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
    //buffer_delta = 0;
    buffer_mag = 0;
    buffer_mag_x = 0;
    buffer_mag_y = 0;
    buffer_mag_z = 0;
    buffer_vel = 0;
    buffer_vel_x = 0;
    buffer_vel_y = 0;
    buffer_vel_z = 0;
    buffer_larm = 0;
    //buffer_mach = 0;
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
    max_wavelengths = 0;
}

double CGridBasic::getCextMeanTab(uint cellID, uint wID) const
{
    if(CextMeanTab != 0)
        return CextMeanTab[wID][cellID];
    return MAX_DOUBLE;
}

double CGridBasic::getCabsMeanTab(uint cellID, uint wID) const
{
    if(CabsMeanTab != 0)
        return CabsMeanTab[wID][cellID];
    return MAX_DOUBLE;
}

double CGridBasic::getCscaMeanTab(uint cellID, uint wID) const
{
    if(CscaMeanTab != 0)
        return CscaMeanTab[wID][cellID];
    return MAX_DOUBLE;
}

double CGridBasic::getNumberDensityTab(uint cellID) const
{
    if(numberDensityTab != 0)
        return numberDensityTab[cellID];
    return MAX_DOUBLE;
}

double CGridBasic::getTotalCellEmissionTab(uint cellID) const
{
    if(totalCellEmissionTab != 0)
        return totalCellEmissionTab[cellID];
    return MAX_DOUBLE;
}

void CGridBasic::setCextMeanTab(double Cext, uint cellID, uint wID)
{
    CextMeanTab[wID][cellID] = Cext;
}

void CGridBasic::setCabsMeanTab(double Cabs, uint cellID, uint wID)
{
    CabsMeanTab[wID][cellID] = Cabs;
}

void CGridBasic::setCscaMeanTab(double Csca, uint cellID, uint wID)
{
    CscaMeanTab[wID][cellID] = Csca;
}

void CGridBasic::setNumberDensityTab(double number_density, uint cellID)
{
    numberDensityTab[cellID] = number_density;
}

void CGridBasic::setTotalCellEmissionTab(double cell_emission, uint cellID)
{
    totalCellEmissionTab[cellID] = cell_emission;
}

void CGridBasic::initPreCalcTables(uint nr_used_wavelengths)
{
    max_wavelengths = nr_used_wavelengths;
    CextMeanTab = new double *[max_wavelengths];
    CabsMeanTab = new double *[max_wavelengths];
    CscaMeanTab = new double *[max_wavelengths];

    for(uint wID = 0; wID < max_wavelengths; wID++)
    {
        CextMeanTab[wID] = new double[max_cells];
        fill(CextMeanTab[wID], CextMeanTab[wID] + max_cells, MAX_DOUBLE);
        CabsMeanTab[wID] = new double[max_cells];
        fill(CabsMeanTab[wID], CabsMeanTab[wID] + max_cells, MAX_DOUBLE);
        CscaMeanTab[wID] = new double[max_cells];
        fill(CscaMeanTab[wID], CscaMeanTab[wID] + max_cells, MAX_DOUBLE);
    }

    numberDensityTab = new double[max_cells];
    fill(numberDensityTab, numberDensityTab + max_cells, MAX_DOUBLE);

    totalCellEmissionTab = new double[max_cells];
    fill(totalCellEmissionTab, totalCellEmissionTab + max_cells, MAX_DOUBLE);
}

void CGridBasic::updateDataRange(cell_basic * cell)
{
    double gas_dens = 0;
    double dust_dens = 0;
    double mx = 0;
    double my = 0;
    double mz = 0;
    double vx = 0;
    double vy = 0;
    double vz = 0;
    double a_alg = 0;
    double dust_temp = 0;
    double gas_temp = 0;
    double mach = 0;
    double delta = 0;
    double a_limit = 0;

    if(size_gd_list > 0)
    {
        for(uint i_dens = 0; i_dens < size_gd_list; i_dens++)
            cell->convertData(data_pos_gd_list[i_dens], conv_dens_in_SI);
        gas_dens = getGasDensity(*cell);

        if(gas_dens > max_gas_dens)
            max_gas_dens = gas_dens;
        if(gas_dens < min_gas_dens)
            min_gas_dens = gas_dens;
    }

    if(size_dd_list > 0)
    {
        for(uint i_dens = 0; i_dens < size_dd_list; i_dens++)
            cell->convertData(data_pos_dd_list[i_dens], conv_dens_in_SI);

        dust_dens = getDustDensity(*cell);

        if(dust_dens > max_dust_dens)
            max_dust_dens = dust_dens;
        if(dust_dens < min_dust_dens)
            min_dust_dens = dust_dens;
    }

    if(!data_pos_dt_list.empty())
    {
        dust_temp = cell->getData(data_pos_dt_list[0]);
        // to do if conversion is implemented

        if(dust_temp > max_dust_temp)
            max_dust_temp = dust_temp;
        if(dust_temp < min_dust_temp)
            min_dust_temp = dust_temp;
    }

    if(data_pos_tg != MAX_UINT)
    {
        gas_temp = cell->getData(data_pos_tg);
        // to do if conversion is implemented
        if(gas_temp > max_gas_temp)
            max_gas_temp = gas_temp;
        if(gas_temp < min_gas_temp)
            min_gas_temp = gas_temp;
    }

    if(data_pos_mx != MAX_UINT)
    {
        cell->convertData(data_pos_mx, conv_Bfield_in_SI);
        mx = cell->getData(data_pos_mx);
    }

    if(data_pos_my != MAX_UINT)
    {
        cell->convertData(data_pos_my, conv_Bfield_in_SI);
        my = cell->getData(data_pos_my);
    }

    if(data_pos_mz != MAX_UINT)
    {
        cell->convertData(data_pos_mz, conv_Bfield_in_SI);
        mz = cell->getData(data_pos_mz);
    }

    if(data_pos_vx != MAX_UINT)
    {
        cell->convertData(data_pos_vx, conv_Vfield_in_SI);
        vx = cell->getData(data_pos_vx);
    }

    if(data_pos_vy != MAX_UINT)
    {
        cell->convertData(data_pos_vy, conv_Vfield_in_SI);
        vy = cell->getData(data_pos_vy);
    }

    if(data_pos_vz != MAX_UINT)
    {
        cell->convertData(data_pos_vz, conv_Vfield_in_SI);
        vz = cell->getData(data_pos_vz);
    }

    if(!data_pos_aalg_list.empty())
    {
        for(uint i = 0; i < data_pos_aalg_list.size(); i++)
        {
            a_alg = cell->getData(data_pos_aalg_list[i]);

            if(a_alg > float(aalg_max))
                aalg_max = (double)a_alg;

            if(a_alg < float(aalg_min))
                aalg_min = (double)a_alg;
        }
    }

    if(!data_pos_amin_list.empty())
    {
        for(uint i = 0; i < data_pos_amin_list.size(); i++)
        {
            double a_min = cell->getData(data_pos_amin_list[i]);

            if(a_min > float(a_min_max))
                a_min_max = a_min;

            if(a_min < float(a_min_min))
                a_min_min = a_min;
        }
    }

    if(!data_pos_amax_list.empty())
    {
        for(uint i = 0; i < data_pos_amax_list.size(); i++)
        {
            double a_max = cell->getData(data_pos_amax_list[i]);

            if(a_max > float(a_max_max))
                a_max_max = a_max;

            if(a_max < float(a_max_min))
                a_max_min = a_max;
        }
    }

    if(!data_pos_size_param_list.empty())
    {
        for(uint i = 0; i < data_pos_size_param_list.size(); i++)
        {
            double size_param = cell->getData(data_pos_size_param_list[i]);

            if(size_param > float(size_param_max))
                size_param_max = size_param;

            if(size_param < float(size_param_min))
                size_param_min = size_param;
        }
    }

    if(data_pos_id != MAX_UINT)
    {
        uint dust_id = cell->getData(data_pos_id);

        if(dust_id > float(dust_id_max))
            dust_id_max = dust_id;

        if(dust_id < float(dust_id_min))
            dust_id_min = dust_id;
    }

    // data positions for synchrotron
    if(data_pos_n_th != MAX_UINT)
    {
        cell->convertData(data_pos_n_th, conv_dens_in_SI);
        double data = cell->getData(data_pos_n_th);

        if(data < min_n_th)
            min_n_th = data;

        if(data > max_n_th)
            max_n_th = data;
    }

    if(data_pos_T_e != MAX_UINT)
    {
        double data = cell->getData(data_pos_T_e);

        if(data < min_T_e)
            min_T_e = data;

        if(data > max_T_e)
            max_T_e = data;
    }

    if(data_pos_n_cr != MAX_UINT)
    {
        cell->convertData(data_pos_n_cr, conv_dens_in_SI);
        double data = cell->getData(data_pos_n_cr);

        if(data < min_n_cr)
            min_n_cr = data;

        if(data > max_n_cr)
            max_n_cr = data;
    }

    if(data_pos_g_min != MAX_UINT)
    {
        double data = cell->getData(data_pos_g_min);

        if(data < min_g_min)
            min_g_min = data;

        if(data > max_g_min)
            max_g_min = data;
    }

    if(data_pos_g_max != MAX_UINT)
    {
        double data = cell->getData(data_pos_g_max);

        if(data < min_g_max)
            min_g_max = data;

        if(data > max_g_max)
            max_g_max = data;
    }

    if(data_pos_p != MAX_UINT)
    {
        double data = cell->getData(data_pos_p);

        if(data < min_p)
            min_p = data;

        if(data > max_p)
            max_p = data;
    }

    double Bfield = sqrt(mx * mx + my * my + mz * mz);
    double Vfield = sqrt(vx * vx + vy * vy + vz * vz);

    if(Bfield > 0)
    {
        if(dust_temp * gas_temp * gas_dens >= 0)
        {
            //delta = CMathFunctions::calc_delta(Bfield, dust_temp, gas_temp, gas_dens) * delta0;
            a_limit = CMathFunctions::calc_larm_limit(Bfield, dust_temp, gas_temp, gas_dens, 0.5, larm_f);
        }
    }
    else
    {
        Bfield = 0;
        //delta = 0;
        a_limit = 0;
    }

    /*if(delta > max_delta)
        max_delta = delta;
    if(delta < min_delta)
        min_delta = delta;*/

    if(Bfield > max_mag)
        max_mag = Bfield;
    if(Bfield < min_mag)
        min_mag = Bfield;

    meanBdir += Vector3D(mx, my, mz);

    if(a_limit > max_larm_limit)
        max_larm_limit = a_limit;
    if(a_limit < min_larm_limit)
        min_larm_limit = a_limit;

    if(Vfield >= 0)
    {
        if(gas_temp > 0)
            mach = Vfield / sqrt(con_kB * gas_temp / (mu * m_H));
        else
            mach = 0;
    }

    if(Vfield > max_vel)
        max_vel = Vfield;
    if(Vfield < min_vel)
        min_vel = Vfield;
    if(mach > max_mach)
        max_mach = mach;
    if(mach < min_mach)
        min_mach = mach;

    meanVdir += Vector3D(vx, vy, vz);
}

bool CGridBasic::fillGridWithOpiateData(uint col_id)
{
    /* uint cell_count = 0;
     uint found_count = 0;
     //#pragma omp parallel for schedule(dynamic)
     for(long i_cell = 0; i_cell < long(max_cells); i_cell++)
     {
         cell_count++;

         cell_basic * cell = cell_list[i_cell];

         uint id = getOpiateID(cell);
         double val = 0;

         if(id != MAX_UINT)
         {
             val = opiate->getData(id, col_id);
             setOpiateTestData(cell, val);
             found_count++;
         }
         else
             setOpiateTestData(cell, 0);
     }

     cout << CLR_LINE;
     cout << " - " << found_count << " of " << max_cells << " cells match with the OPIATE
 paramter file." << endl;
     */
    return true;
}

uint CGridBasic::validateDataPositions(parameters & param)
{
    uint tmp_data_offset = 0;

    cout << CLR_LINE;

    if(size_gd_list == 0)
    {
        cout << ERROR_LINE << "Grid contains no gas (number) density!" << endl;
        cout << "       No RT calculations possible!" << endl;
        return MAX_UINT;
    }

    if(nr_mixtures > 0 && (param.isTemperatureSimulation() || param.isRatSimulation() ||
                           param.getCommand() == CMD_DUST_EMISSION || param.getCommand() == CMD_LINE_EMISSION ||
                           param.getCommand() == CMD_FORCE || param.getCommand() == CMD_PROBING))
    {
        // Get Number of temperature fields for temperature calculation
        if(size_dd_list > 0)
            nr_densities = size_dd_list;
        else
            nr_densities = size_gd_list;

        // Precalculate the number of temperature entries, if the grid has a
        // temperature for each grain size or stochastically heated grains
        for(uint i_density = 0; i_density < nr_densities; i_density++)
        {
            multi_temperature_entries += nr_dust_temp_sizes[i_density] + 1;
            stochastic_temperature_entries += nr_stochastic_sizes[i_density] + 1;
        }

        // Check for a valid combination between densities and dust mixtures
        if(nr_densities > 1 && nr_mixtures < nr_densities)
        {
            cout << ERROR_LINE << "Amount of densities in the grid (" << nr_densities
                 << ") does not fit with the defined dust mixtures (" << nr_mixtures << ")!\n"
                 << "(Use a grid with only one density distribution or define more/less "
                    "dust mixtures!)"
                 << endl;
            return MAX_UINT;
        }

        // Init list to know how many dust sizes are used per dust component
        size_skip = new uint[nr_densities];

        // Calculate the entries for the temperature that have to be added
        if(param.getDustTempMulti())
            for(uint i_density = 0; i_density < nr_densities; i_density++)
                size_skip[i_density] = nr_dust_temp_sizes[i_density];
        else if(param.getStochasticHeatingMaxSize() > 0 && !param.getSaveRadiationField())
            for(uint i_density = 0; i_density < nr_densities; i_density++)
                size_skip[i_density] = nr_stochastic_sizes[i_density];
        else
            for(uint i_density = 0; i_density < nr_densities; i_density++)
                size_skip[i_density] = 1;
    }

    switch(param.getCommand())
    {
        case CMD_SYNCHROTRON:
            if(CheckSynchrotron(param) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_OPIATE:
            if(CheckOpiate(param) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_TEMP:
            if(CheckTemp(param, tmp_data_offset) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_TEMP_RAT:
            if(CheckTemp(param, tmp_data_offset) == MAX_UINT)
                return MAX_UINT;

            if(CheckRat(param, tmp_data_offset) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_RAT:
            if(CheckRat(param, tmp_data_offset) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_DUST_EMISSION:
            if(CheckDustEmission(param) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_DUST_SCATTERING:
            if(CheckDustScattering(param) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_FORCE:
            if(CheckRadiationForce(param) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_LINE_EMISSION:
            if(CheckLineEmission(param) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_PROBING:
            if(CheckProbing(param) == MAX_UINT)
                return MAX_UINT;
            break;

        default:
            cout << ERROR_LINE << "Command is unknown!" << endl;
            return MAX_UINT;
    }

    return tmp_data_offset;
}

void CGridBasic::printPhysicalParameters()
{
    cout << "- Volume (total, cells, diff)   : " << total_volume << " [m^3], " << cell_volume << " [m^3], "
         << float(100.0 * abs(total_volume - cell_volume) / max(total_volume, cell_volume)) << " [%]" << endl;
    cout << "- Total gas mass                : " << total_gas_mass / M_sun << " [M_sun], " << total_gas_mass
         << " [kg]" << endl;
    cout << "- Grid length         (min,max) : [" << min_len << ", " << max_len << "] [m]" << endl;
    if(gas_is_mass_density)
        cout << "- Gas mass density    (min,max) : [" << min_gas_dens << ", " << max_gas_dens << "] [kg m^-3]"
             << endl;
    else
        cout << "- Gas number density  (min,max) : [" << min_gas_dens << ", " << max_gas_dens << "] [m^-3]"
             << endl;
    if(size_dd_list > 0)
    {
        if(dust_is_mass_density)
            cout << "- Dust mass density   (min,max) : [" << min_dust_dens << ", " << max_dust_dens
                 << "] [kg m^-3]" << endl;
        else
            cout << "- Dust number density (min,max) : [" << min_dust_dens << ", " << max_dust_dens
                 << "] [m^-3]" << endl;
    }
    if(data_pos_tg != MAX_UINT)
        cout << "- Gas temperature     (min,max) : [" << min_gas_temp << ", " << max_gas_temp << "] [K]"
             << endl;
    else
        cout << "- Gas temperature     (min,max) : none" << endl;

    if(!data_pos_dt_list.empty())
        cout << "- Dust temperature    (min,max) : [" << min_dust_temp << ", " << max_dust_temp << "] [K]"
             << endl;
    else
        cout << "- Dust temperature    (min,max) : none" << endl;

    if(data_pos_mx != MAX_UINT)
    {
        meanBdir.normalize();
        cout << "- Magnetic field      (min,max) : [" << min_mag << ", " << max_mag << "] [T]" << endl;
        cout << "- Mean direction      (norm.)   : X: " << meanBdir.X() << " Y: " << meanBdir.Y()
             << " Z: " << meanBdir.Z() << endl;
        //cout << "- Delta0              (min,max) : [" << min_delta << ", " << max_delta << "] [m]" << endl;
        cout << "- Larm. limit         (min,max) : [" << min_larm_limit << ", " << max_larm_limit << "] [m]"
             << endl;
    }
    else
        cout << "- Magnetic field      (min,max) : none" << endl;

    if(!data_pos_aalg_list.empty())
        cout << "- a_alig              (min,max) : [" << aalg_min << ", " << aalg_max << "] [m]" << endl;

    if(data_pos_vx != MAX_UINT)
    {
        meanVdir.normalize();
        cout << "- Velocity field      (min,max) : [" << min_vel << ", " << max_vel << "] [m/s]" << endl;
        cout << "- Mean direction      (norm.)   : X: " << meanVdir.X() << " Y: " << meanVdir.Y()
             << " Z: " << meanVdir.Z() << endl;
        cout << "- Mach number         (min,max) : [" << min_mach << ", " << max_mach << "]" << endl;
    }

    if(data_pos_amin_list.size()>0)
        cout << "- Minimum grain size  (min,max) : [" << a_min_min << ", " << a_min_max << "] [m]" << endl;

    if(data_pos_amax_list.size()>0)
        cout << "- Maximum grain size  (min,max) : [" << a_max_min << ", " << a_max_max << "] [m]" << endl;

    if(data_pos_size_param_list.size()>0)
        cout << "- Dust size parameter (min,max) : [" << size_param_min << ", " << size_param_max << "]"
             << endl;

    if(data_pos_id != MAX_UINT)
        cout << "- Dust mixture ID     (min,max) : [" << dust_id_min << ", " << dust_id_max << "]" << endl;

    if(data_pos_n_cr != MAX_UINT)
    {
        cout << "- CR el. density      (min,max) : [" << min_n_cr << "; " << max_n_cr << "] [m^-3]" << endl;

        if(data_pos_g_min != MAX_UINT)
            cout << "- Gamma_min           (min,max) : [" << min_g_min << "; " << max_g_min << "]" << endl;

        if(data_pos_g_max != MAX_UINT)
            cout << "- Gamma_max           (min,max) : [" << min_g_max << "; " << max_g_max << "]" << endl;

        if(data_pos_p != MAX_UINT)
            cout << "- El. energy index p  (min,max) : [" << min_p << "; " << max_p << "]" << endl;
    }
    else
        cout << "- CR el. density      (min,max) : none   " << endl;

    if(data_pos_n_th != MAX_UINT)
    {
        cout << "- Therm. el. density  (min,max) : [" << min_n_th << "; " << max_n_th << "] [m]" << endl;

        if(data_pos_T_e != MAX_UINT)
        {
            if(min_T_e == 1e300)
                cout << "- Electron temperature          : same as dust temperature" << endl;
            else
                cout << "- Electron temp.      (min,max) : [" << min_T_e << "; " << max_T_e << "] [K]"
                     << endl;
        }
    }
    else
        cout << "- Therm. el. density  (min,max) : none" << endl;

    if(nrOfOpiateIDs > 0 || nrOfDensRatios > 0)
    {
        cout << SEP_LINE;
        cout << "Additional grid data:" << endl;
    }

    if(nrOfDensRatios > 0)
    {
        cout << "- Density. ratio IDs: ";
        cout << 1 << ":" << pos_GasSpecRatios[0];

        for(uint i = 1; i < nrOfDensRatios; i++)
            cout << ", " << i + 1 << ":" << pos_GasSpecRatios[i];

        cout << endl;
    }

    if(nrOfOpiateIDs > 0)
    {
        cout << "- Unique OPIATE IDs : ";
        cout << 1 << ":" << pos_OpiateIDS[0];

        for(uint i = 1; i < nrOfOpiateIDs; i++)
            cout << ", " << i + 1 << ":" << pos_OpiateIDS[i];

        cout << endl;
    }

    if(data_pos_op != UINT_MAX)
        cout << " - Unique OPIATE IDs" << endl;
}

bool CGridBasic::writeMidplaneFits(string data_path, parameters & param, uint bins, bool all)
{
    bool res = true;

    if(bins == 0)
        return res;

    int cmd = param.getCommand();


    cout << CLR_LINE;
    cout << " -> Allocating memory for plotting midplane files ...             \r" << flush;

    if(all)
    {
        plt_gas_dens = (size_gd_list > 0) && param.isInPlotList(GRIDgas_dens);
        plt_mol_dens = (nrOfDensRatios>0 && param.isInPlotList(GRIDratio) );
        plt_dust_dens = (size_dd_list > 0) && param.isInPlotList(GRIDdust_dens);
        plt_gas_temp = (data_pos_tg != MAX_UINT) && param.isInPlotList(GRIDgas_temp);

        plt_mag = (data_pos_mx != MAX_UINT) && (data_pos_my != MAX_UINT) && (data_pos_my != MAX_UINT) &&
                  param.isInPlotList(GRIDmx) && param.isInPlotList(GRIDmy) && param.isInPlotList(GRIDmz);

        plt_vel = (data_pos_vx != MAX_UINT) && (data_pos_vy != MAX_UINT) && (data_pos_vz != MAX_UINT) &&
                  param.isInPlotList(GRIDvx) && param.isInPlotList(GRIDvy) && param.isInPlotList(GRIDvz);

        //plt_delta = plt_gas_temp && plt_mag && (!data_pos_dt_list.empty());
        plt_larm = plt_gas_temp && plt_mag && (!data_pos_dt_list.empty());
        //plt_mach = plt_vel && plt_gas_temp;

        plt_dust_id = (data_pos_id != MAX_UINT);

        plt_amin = (data_pos_amin_list.size()>0) && param.isInPlotList(GRIDa_min);
        plt_amax = (data_pos_amax_list.size()>0) && param.isInPlotList(GRIDa_max);
        plt_size_param = (data_pos_size_param_list.size()>0) && param.isInPlotList(GRIDq);

        plt_n_th = (data_pos_n_th != MAX_UINT) && param.isInPlotList(GRIDn_th);
        plt_T_e = (data_pos_T_e != MAX_UINT) && param.isInPlotList(GRIDT_e);
        plt_n_cr = (data_pos_n_cr != MAX_UINT) && param.isInPlotList(GRIDn_cr);
        plt_g_min = (data_pos_g_min != MAX_UINT) && param.isInPlotList(GRIDg_min);
        plt_g_max = (data_pos_g_max != MAX_UINT) && param.isInPlotList(GRIDg_max);
        plt_p = (data_pos_p != MAX_UINT) && param.isInPlotList(GRIDp);

        if(cmd != CMD_RAT && cmd != CMD_TEMP_RAT)
        {
            plt_rat = (!data_pos_aalg_list.empty()) && param.isInPlotList(GRIDa_alg);
            plt_avg_th = (data_pos_avg_th != MAX_UINT) && param.isInPlotList(GRIDavg_th);
            plt_avg_dir = (data_pos_avg_dir != MAX_UINT) && param.isInPlotList(GRIDavg_dir);
        }

        if(cmd != CMD_TEMP && cmd != CMD_TEMP_RAT)
            plt_dust_temp = (!data_pos_dt_list.empty()) && param.isInPlotList(GRIDdust_temp);

        // if(getRadiationFieldAvailable())
        {
            switch(param.getWriteRadiationField())
            {
                default:
                    plt_u_rad = false;
                    plt_rad_field1 = false;
                    break;

                case 1:
                    plt_u_rad = (cmd == CMD_RAT || cmd == CMD_TEMP_RAT);
                    plt_rad_field1 = false;
                    break;

                case 2:
                    plt_u_rad = false;
                    plt_rad_field1 = true;
                    break;

                case 3:
                    plt_u_rad = false;
                    plt_rad_field1 = true;
                    nr_rad_field_comp = 4;
                    break;
            }

            if(param.getWriteGZero())
                plt_g_zero1 = true;
        }
    }
    else
    {
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
        plt_amin = false;
        plt_amax = false;
        plt_size_param = false;
        plt_rad_field1 = false;
        plt_g_zero1 = false;
        plt_u_rad = false;
        plt_n_th = false;
        plt_T_e = false;
        plt_n_cr = false;
        plt_g_min = false;
        plt_g_max = false;
        plt_p = false;

        plt_avg_th = false;
        plt_avg_dir = false;

        if(cmd == CMD_TEMP || cmd == CMD_TEMP_RAT)
        {
            if(param.getAdjTgas() > 0)
                plt_gas_temp = param.isInPlotList(GRIDgas_temp);

            plt_dust_temp = param.isInPlotList(GRIDdust_temp);
        }

        if(cmd == CMD_RAT || cmd == CMD_TEMP_RAT)
        {
            plt_rat = param.isInPlotList(GRIDa_alg);
            plt_avg_th = param.isInPlotList(GRIDavg_th);
            plt_avg_dir = param.isInPlotList(GRIDavg_dir);
        }

        switch(param.getWriteRadiationField())
        {
            default:
                plt_u_rad = false;
                plt_rad_field1 = false;
                break;

            case 1:
                plt_u_rad = (cmd == CMD_RAT || cmd == CMD_TEMP_RAT);
                plt_rad_field1 = false;
                break;

            case 2:
                plt_u_rad = false;
                plt_rad_field1 = true;

                if(!spec_length_as_vector)
                {
                    cout << INFO_LINE << "The full radiation field can only be saved if it was used by the simulation" << endl;
                    cout << "  (when saving the radiation field in the grid or calculating RATs)!" << endl;
                }
                break;

            case 3:
                plt_u_rad = false;
                plt_rad_field1 = true;
                nr_rad_field_comp = 4;
                break;
        }

        if(param.getWriteGZero())
            plt_g_zero1 = true;
    }

    //+ uint(plt_delta) + uint(plt_mach)
    uint nr_parameters = uint(plt_gas_dens) + uint(plt_dust_dens) + uint(plt_gas_temp) + uint(plt_dust_temp) +
                         4 * uint(plt_mag) + 4 * uint(plt_vel) + uint(plt_rat) + 
                         uint(plt_larm) + uint(plt_dust_id) +
                         uint(plt_rad_field1) * nr_rad_field_comp * WL_STEPS + uint(plt_g_zero1) +
                         uint(plt_u_rad) + uint(plt_n_th) + uint(plt_T_e) + uint(plt_n_cr) + uint(plt_g_min) +
                         uint(plt_g_max) + uint(plt_p) + uint(plt_avg_th) + uint(plt_avg_dir);

    if(nr_parameters == 0)
        return res;

    if(plt_gas_dens)
        if(nr_densities > 1 && size_gd_list >= nr_densities)
            nr_parameters += nr_densities;
    if(plt_dust_dens)
        if(nr_densities > 1 && size_dd_list >= nr_densities)
            nr_parameters += nr_densities;
    if(plt_dust_temp)
        if(nr_densities > 1 && data_pos_dt_list.size() >= nr_densities)
            nr_parameters += nr_densities;
    if(plt_mol_dens)
        nr_parameters += nrOfDensRatios;

    long naxis = 4;
    long naxes[4] = { uint(bins), uint(bins), 3, nr_parameters };
    uint per_max = 3 * bins * bins;

    double max_midplane_len = (max_len / param.getMidplaneZoom());

    dlist midplane_3d_param = param.getMidplane3dParams();
    double z_step, off_z, shift_z = 0;
    uint plane_3d = 0;
    if(midplane_3d_param.size() == 4)
    {
        plane_3d = midplane_3d_param[0];

        if(midplane_3d_param[1] != 0)
        {
            naxes[2] = uint(midplane_3d_param[1]);
            per_max = bins * bins * midplane_3d_param[1];
        }
        else
        {
            naxes[2] = uint(bins);
            per_max = bins * bins * bins;
        }

        if(midplane_3d_param[2] != 0 || midplane_3d_param[3] != 0)
        {
            z_step = (midplane_3d_param[3] - midplane_3d_param[2]) / double(naxes[2]);
            off_z = 0.5 * z_step;
            shift_z = (midplane_3d_param[3] + midplane_3d_param[2]) / 2.0;
        }
        else
        {
            z_step = max_midplane_len / double(naxes[2]);
            off_z = 0.5 * z_step;
        }
    }
    else
    {
        z_step = max_midplane_len / double(bins);
        off_z = 0.5 * z_step;
    }

    double xy_step = max_midplane_len / double(bins);
    double off_xy = 0.5 * xy_step;
    int b_limit_z, b_limit_xy;

    if(naxes[2] % 2)
    {
        b_limit_z = (naxes[2] - 1) / 2;
        off_z = 0;
    }
    else
        b_limit_z = naxes[2] / 2;

    if(naxes[0] % 2)
    {
        b_limit_xy = (naxes[0] - 1) / 2;
        off_xy = 0;
    }
    else
        b_limit_xy = naxes[0] / 2;

    ullong per_counter = 0;

    // auto_ptr<CCfits::FITS> pFits(0);
    unique_ptr<CCfits::FITS> pFits;

    try
    {
        string path_out = data_path + "midplane" + FITS_COMPRESS_EXT;
        if(midplane_3d_param.size() == 4)
            path_out = data_path + "midplane_3d" + FITS_COMPRESS_EXT;
        remove(path_out.c_str());
        pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
    }
    catch(CCfits::FITS::CantCreate)
    {
        return false;
    }

    long nelements = bins * bins;

    valarray<double> array_gas_dens(nelements);
    valarray<double> array_mol_dens(nelements);
    valarray<double> array_dust_dens(nelements);
    valarray<double> array_gas_temp(nelements);
    valarray<double> array_dust_temp(nelements);
    valarray<double> array_rat(nelements);
    //valarray<double> array_delta(nelements);
    valarray<double> array_mag(nelements);
    valarray<double> array_mag_x(nelements);
    valarray<double> array_mag_y(nelements);
    valarray<double> array_mag_z(nelements);
    valarray<double> array_vel(nelements);
    valarray<double> array_vel_x(nelements);
    valarray<double> array_vel_y(nelements);
    valarray<double> array_vel_z(nelements);
    valarray<double> array_larm(nelements);
    valarray<double> array_mach(nelements);
    valarray<double> array_dust_mixture(nelements);
    valarray<double> array_amin(nelements);
    valarray<double> array_amax(nelements);
    valarray<double> array_size_param(nelements);
    valarray<double> array_rad_field(nelements);
    valarray<double> array_g_zero1(nelements);
    valarray<double> array_u_rad(nelements);
    valarray<double> array_n_th(nelements);
    valarray<double> array_T_e(nelements);
    valarray<double> array_n_cr(nelements);
    valarray<double> array_g_min(nelements);
    valarray<double> array_g_max(nelements);
    valarray<double> array_p(nelements);

    valarray<double> array_avg_th(nelements);
    valarray<double> array_avg_dir(nelements);

    if(plt_gas_dens)
    {
        buffer_gas_dens = new double *[nelements];
        for(long i_cell = 0; i_cell < nelements; i_cell++)
        {
            // +1 for the average/sum of the quantity, but only if multiple quantities are
            // in the grid
            if(nr_densities > 1 && size_gd_list == nr_densities)
                buffer_gas_dens[i_cell] = new double[nr_densities + 1];
            else
                buffer_gas_dens[i_cell] = new double[nr_densities];
        }
    }

    if(plt_mol_dens)
    {
        buffer_mol_dens = new double *[nelements];
        for(long i_cell = 0; i_cell < nelements; i_cell++)
        {
            buffer_mol_dens[i_cell] = new double[nrOfDensRatios];
        }
    }

    if(plt_dust_dens)
    {
        buffer_dust_dens = new double *[nelements];
        for(long i_cell = 0; i_cell < nelements; i_cell++)
        {
            // +1 for the average/sum of the quantity, but only if multiple quantities are
            // in the grid
            if(nr_densities > 1 && size_dd_list == nr_densities)
                buffer_dust_dens[i_cell] = new double[nr_densities + 1];
            else
                buffer_dust_dens[i_cell] = new double[nr_densities];
        }
    }
    if(plt_gas_temp)
        buffer_gas_temp = new double[nelements];
    if(plt_dust_temp)
    {
        buffer_dust_temp = new double *[nelements];
        for(long i_cell = 0; i_cell < nelements; i_cell++)
        {
            // +1 for the average/sum of the quantity, but only if multiple quantities are
            // in the grid
            if(nr_densities > 1 && data_pos_dt_list.size() >= nr_densities)
                buffer_dust_temp[i_cell] = new double[nr_densities + 1];
            else
                buffer_dust_temp[i_cell] = new double[nr_densities];
        }
    }
    if(plt_rat)
    {
        buffer_rat = new double *[nelements];
        for(long i_cell = 0; i_cell < nelements; i_cell++)
            buffer_rat[i_cell] = new double[data_pos_aalg_list.size()];
    }
        
    //if(plt_delta)
    //    buffer_delta = new double[nelements];
        
    if(plt_mag)
    {
        buffer_mag = new double[nelements];
        buffer_mag_x = new double[nelements];
        buffer_mag_y = new double[nelements];
        buffer_mag_z = new double[nelements];
    }
    if(plt_vel)
    {
        buffer_vel = new double[nelements];
        buffer_vel_x = new double[nelements];
        buffer_vel_y = new double[nelements];
        buffer_vel_z = new double[nelements];
    }
    if(plt_larm)
        buffer_larm = new double[nelements];
    //if(plt_mach)
    //    buffer_mach = new double[nelements];
    if(plt_dust_id)
        buffer_dust_mixture = new double[nelements];
    if(plt_amin)
        buffer_dust_amin = new double[nelements];
    if(plt_amax)
        buffer_dust_amax = new double[nelements];
    if(plt_size_param)
        buffer_dust_size_param = new double[nelements];
    if(plt_rad_field1)
    {
        buffer_rad_field = new double **[nelements];
        for(long i_cell = 0; i_cell < nelements; i_cell++)
        {
            buffer_rad_field[i_cell] = new double *[WL_STEPS];
            for(uint wID = 0; wID < WL_STEPS; wID++)
                buffer_rad_field[i_cell][wID] = new double[nr_rad_field_comp];
        }
    }
    if(plt_g_zero1)
        buffer_g_zero1 = new double[nelements];
    if(plt_u_rad)
        buffer_u_rad = new double[nelements];
    if(plt_n_th)
        buffer_n_th = new double[nelements];
    if(plt_T_e)
        buffer_T_e = new double[nelements];
    if(plt_n_cr)
        buffer_n_cr = new double[nelements];
    if(plt_g_min)
        buffer_g_min = new double[nelements];
    if(plt_g_max)
        buffer_g_max = new double[nelements];
    if(plt_p)
        buffer_p = new double[nelements];
    if(plt_avg_th)
        buffer_avg_th = new double[nelements];
    if(plt_avg_dir)
        buffer_avg_dir = new double[nelements];

    vector<long> fpixel(4);

    fpixel[0] = 1;
    fpixel[1] = 1;
    fpixel[2] = 0;

    if(midplane_3d_param.size() == 4)
    {
        for(int l = -b_limit_z; l <= b_limit_z; l++)
        {
            if(l == 0 && naxes[2] % 2 == 0)
                continue;

            fpixel[2]++;

            #pragma omp parallel for schedule(dynamic)
            for(long i_cell = 0; i_cell < nelements; i_cell++)
            {
                int j = (i_cell % bins);

                int k = i_cell / bins - b_limit_xy;

                j -= b_limit_xy;

                if(bins % 2 == 0)
                    if(k > -1)
                        k++;

                if(bins % 2 == 0)
                    if(j > -1)
                        j++;

                double tx, ty, tz;
                setPlaneParameter(plane_3d, xy_step, off_xy, z_step, off_z, shift_z, j, k, l, tx, ty, tz);

                fillMidplaneBuffer(tx, ty, tz, i_cell);

                per_counter++;
            }

            fpixel[3] = 0;
            if(plt_gas_dens)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_gas_dens[i_cell] = buffer_gas_dens[i_cell][0];
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_gas_dens);

                if(nr_densities > 1 && size_gd_list >= nr_densities)
                    for(uint i_density = 0; i_density < nr_densities; i_density++)
                    {
                        for(long i_cell = 0; i_cell < nelements; i_cell++)
                            array_gas_dens[i_cell] = buffer_gas_dens[i_cell][i_density + 1];
                        fpixel[3]++;
                        pFits->pHDU().write(fpixel, nelements, array_gas_dens);
                    }
            }
            if(plt_mol_dens)
            {
                for(uint i_density = 0; i_density < nrOfDensRatios; i_density++)
                {
                    for(long i_cell = 0; i_cell < nelements; i_cell++)
                        array_mol_dens[i_cell] = buffer_mol_dens[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_mol_dens);
                }
            }
            if(plt_dust_dens)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_dust_dens[i_cell] = buffer_dust_dens[i_cell][0];
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_dust_dens);

                if(nr_densities > 1 && size_dd_list >= nr_densities)
                    for(uint i_density = 0; i_density < nr_densities; i_density++)
                    {
                        for(long i_cell = 0; i_cell < nelements; i_cell++)
                            array_dust_dens[i_cell] = buffer_dust_dens[i_cell][i_density + 1];
                        fpixel[3]++;
                        pFits->pHDU().write(fpixel, nelements, array_dust_dens);
                    }
            }
            if(plt_gas_temp)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_gas_temp[i_cell] = buffer_gas_temp[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_gas_temp);
            }
            if(plt_dust_temp)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_dust_temp[i_cell] = buffer_dust_temp[i_cell][0];
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_dust_temp);

                if(nr_densities > 1 && data_pos_dt_list.size() >= nr_densities)
                    for(uint i_density = 0; i_density < nr_densities; i_density++)
                    {
                        for(long i_cell = 0; i_cell < nelements; i_cell++)
                            array_dust_temp[i_cell] = buffer_dust_temp[i_cell][i_density + 1];
                        fpixel[3]++;
                        pFits->pHDU().write(fpixel, nelements, array_dust_temp);
                    }
            }
            if(plt_rat)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(long i_cell = 0; i_cell < nelements; i_cell++)
                        array_rat[i_cell] = buffer_rat[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_rat);
                }
            }
            
            /*if(plt_delta)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_delta[i_cell] = buffer_delta[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_delta);
            }*/
            
            if(plt_mag)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                {
                    array_mag[i_cell] = buffer_mag[i_cell];
                    array_mag_x[i_cell] = buffer_mag_x[i_cell];
                    array_mag_y[i_cell] = buffer_mag_y[i_cell];
                    array_mag_z[i_cell] = buffer_mag_z[i_cell];
                }

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_mag);
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_mag_x);
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_mag_y);
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_mag_z);
            }
            if(plt_vel)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                {
                    array_vel[i_cell] = buffer_vel[i_cell];
                    array_vel_x[i_cell] = buffer_vel_x[i_cell];
                    array_vel_y[i_cell] = buffer_vel_y[i_cell];
                    array_vel_z[i_cell] = buffer_vel_z[i_cell];
                }

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_vel);
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_vel_x);
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_vel_y);
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_vel_z);
            }
            if(plt_larm)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_larm[i_cell] = buffer_larm[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_larm);
            }
            
            /*if(plt_mach)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_mach[i_cell] = buffer_mach[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_mach);
            }*/
            
            if(plt_dust_id)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_dust_mixture[i_cell] = buffer_dust_mixture[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_dust_mixture);
            }
            if(plt_amin)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_amin[i_cell] = buffer_dust_amin[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_amin);
            }
            if(plt_amax)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_amax[i_cell] = buffer_dust_amax[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_amax);
            }
            if(plt_size_param)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_size_param[i_cell] = buffer_dust_size_param[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_size_param);
            }
            if(plt_rad_field1)
            {
                for(uint i_comp = 0; i_comp < nr_rad_field_comp; i_comp++)
                    for(uint wID = 0; wID < WL_STEPS; wID++)
                    {
                        for(long i_cell = 0; i_cell < nelements; i_cell++)
                            array_rad_field[i_cell] = buffer_rad_field[i_cell][wID][i_comp];

                        fpixel[3]++;
                        pFits->pHDU().write(fpixel, nelements, array_rad_field);
                    }
            }
            if(plt_g_zero1)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_g_zero1[i_cell] = buffer_g_zero1[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_g_zero1);
            }
            if(plt_u_rad)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_u_rad[i_cell] = buffer_u_rad[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_u_rad);
            }
            if(plt_n_th)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_n_th[i_cell] = buffer_n_th[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_n_th);
            }
            if(plt_T_e)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_T_e[i_cell] = buffer_T_e[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_T_e);
            }
            if(plt_n_cr)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_n_cr[i_cell] = buffer_n_cr[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_n_cr);
            }
            if(plt_g_min)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_g_min[i_cell] = buffer_g_min[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_g_min);
            }
            if(buffer_g_max)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_g_max[i_cell] = buffer_g_max[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_g_max);
            }
            if(plt_p)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_p[i_cell] = buffer_p[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_p);
            }
            if(plt_avg_th)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_avg_th[i_cell] = buffer_avg_th[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_avg_th);
            }
            if(plt_avg_dir)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_avg_dir[i_cell] = buffer_avg_dir[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_avg_dir);
            }
        }
    }
    else
    {
        for(int i = 1; i <= 3; i++)
        {
            fpixel[2] = i;

            #pragma omp parallel for schedule(dynamic)
            for(long i_cell = 0; i_cell < nelements; i_cell++)
            {
                int j = (i_cell % bins);

                int k = i_cell / bins - b_limit_xy;

                j -= b_limit_xy;

                if(bins % 2 == 0)
                    if(k > -1)
                        k++;

                if(bins % 2 == 0)
                    if(j > -1)
                        j++;

                double tx, ty, tz;

                setPlaneParameter(i, xy_step, off_xy, 0, 0, 0, j, k, 0, tx, ty, tz);

                fillMidplaneBuffer(tx, ty, tz, i_cell);

                per_counter++;
            }

            fpixel[3] = 0;
            if(plt_gas_dens)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_gas_dens[i_cell] = buffer_gas_dens[i_cell][0];
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_gas_dens);

                if(nr_densities > 1 && size_gd_list >= nr_densities)
                    for(uint i_density = 0; i_density < nr_densities; i_density++)
                    {
                        for(long i_cell = 0; i_cell < nelements; i_cell++)
                            array_gas_dens[i_cell] = buffer_gas_dens[i_cell][i_density + 1];
                        fpixel[3]++;
                        pFits->pHDU().write(fpixel, nelements, array_gas_dens);
                    }
            }
            if(plt_mol_dens)
            {
                for(uint i_density = 0; i_density < nrOfDensRatios; i_density++)
                {
                    for(long i_cell = 0; i_cell < nelements; i_cell++)
                        array_mol_dens[i_cell] = buffer_mol_dens[i_cell][i_density];
                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_mol_dens);
                }
            }
            if(plt_dust_dens)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_dust_dens[i_cell] = buffer_dust_dens[i_cell][0];
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_dust_dens);

                if(nr_densities > 1 && size_dd_list >= nr_densities)
                    for(uint i_density = 0; i_density < nr_densities; i_density++)
                    {
                        for(long i_cell = 0; i_cell < nelements; i_cell++)
                            array_dust_dens[i_cell] = buffer_dust_dens[i_cell][i_density + 1];
                        fpixel[3]++;
                        pFits->pHDU().write(fpixel, nelements, array_dust_dens);
                    }
            }
            if(plt_gas_temp)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_gas_temp[i_cell] = buffer_gas_temp[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_gas_temp);
            }
            if(plt_dust_temp)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_dust_temp[i_cell] = buffer_dust_temp[i_cell][0];
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_dust_temp);

                if(nr_densities > 1 && data_pos_dt_list.size() >= nr_densities)
                    for(uint i_density = 0; i_density < nr_densities; i_density++)
                    {
                        for(long i_cell = 0; i_cell < nelements; i_cell++)
                            array_dust_temp[i_cell] = buffer_dust_temp[i_cell][i_density + 1];
                        fpixel[3]++;
                        pFits->pHDU().write(fpixel, nelements, array_dust_temp);
                    }
            }
            if(plt_rat)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(long i_cell = 0; i_cell < nelements; i_cell++)
                        array_rat[i_cell] = buffer_rat[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_rat);
                }
            }
            
            /*if(plt_delta)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_delta[i_cell] = buffer_delta[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_delta);
            }*/
            
            if(plt_mag)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                {
                    array_mag[i_cell] = buffer_mag[i_cell];
                    array_mag_x[i_cell] = buffer_mag_x[i_cell];
                    array_mag_y[i_cell] = buffer_mag_y[i_cell];
                    array_mag_z[i_cell] = buffer_mag_z[i_cell];
                }

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_mag);
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_mag_x);
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_mag_y);
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_mag_z);
            }
            if(plt_vel)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                {
                    array_vel[i_cell] = buffer_vel[i_cell];
                    array_vel_x[i_cell] = buffer_vel_x[i_cell];
                    array_vel_y[i_cell] = buffer_vel_y[i_cell];
                    array_vel_z[i_cell] = buffer_vel_z[i_cell];
                }

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_vel);
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_vel_x);
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_vel_y);
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_vel_z);
            }
            if(plt_larm)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_larm[i_cell] = buffer_larm[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_larm);
            }
            
            /*if(plt_mach)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_mach[i_cell] = buffer_mach[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_mach);
            }*/
            
            if(plt_dust_id)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_dust_mixture[i_cell] = buffer_dust_mixture[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_dust_mixture);
            }
            if(plt_amin)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_amin[i_cell] = buffer_dust_amin[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_amin);
            }
            if(plt_amax)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_amax[i_cell] = buffer_dust_amax[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_amax);
            }
            if(plt_size_param)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_size_param[i_cell] = buffer_dust_size_param[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_size_param);
            }
            if(plt_rad_field1)
            {
                for(uint i_comp = 0; i_comp < nr_rad_field_comp; i_comp++)
                    for(uint wID = 0; wID < WL_STEPS; wID++)
                    {
                        for(long i_cell = 0; i_cell < nelements; i_cell++)
                            array_rad_field[i_cell] = buffer_rad_field[i_cell][wID][i_comp];

                        fpixel[3]++;
                        pFits->pHDU().write(fpixel, nelements, array_rad_field);
                    }
            }
            if(plt_g_zero1)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_g_zero1[i_cell] = buffer_g_zero1[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_g_zero1);
            }
            if(plt_u_rad)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_u_rad[i_cell] = buffer_u_rad[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_u_rad);
            }
            if(plt_n_th)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_n_th[i_cell] = buffer_n_th[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_n_th);
            }
            if(plt_T_e)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_T_e[i_cell] = buffer_T_e[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_T_e);
            }
            if(plt_n_cr)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_n_cr[i_cell] = buffer_n_cr[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_n_cr);
            }
            if(plt_g_min)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_g_min[i_cell] = buffer_g_min[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_g_min);
            }
            if(buffer_g_max)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_g_max[i_cell] = buffer_g_max[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_g_max);
            }
            if(plt_p)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_p[i_cell] = buffer_p[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_p);
            }
            if(plt_avg_th)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_avg_th[i_cell] = buffer_avg_th[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_avg_th);
            }
            if(plt_avg_dir)
            {
                for(long i_cell = 0; i_cell < nelements; i_cell++)
                    array_avg_dir[i_cell] = buffer_avg_dir[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_avg_dir);
            }
        }
    }

    double bin_width = max_midplane_len / bins;
    double first_pix_val = -max_midplane_len / 2.0 + (bin_width / 2.0);

    // Grid
    pFits->pHDU().addKey("CTYPE1", "PARAM", "type of unit 1");
    pFits->pHDU().addKey("CRVAL1", first_pix_val, "value of axis 1");
    pFits->pHDU().addKey("CRPIX1", 1, "pixel where CRVAL1 is defined ");
    pFits->pHDU().addKey("CDELT1", bin_width, "delta of axis 1");
    pFits->pHDU().addKey("CUNIT1", "m", "unit of axis 1");

    // Alternatively as AU grid
    pFits->pHDU().addKey("CTYPE1B", "PARAM", "type of unit 1");
    pFits->pHDU().addKey("CRVAL1B", first_pix_val / con_AU, "value of axis 1");
    pFits->pHDU().addKey("CRPIX1B", 1, "pixel where CRVAL1 is defined ");
    pFits->pHDU().addKey("CDELT1B", bin_width / con_AU, "delta of axis 1");
    pFits->pHDU().addKey("CUNIT1B", "AU", "unit of axis 1");

    // Alternatively as pc grid
    pFits->pHDU().addKey("CTYPE1C", "PARAM", "type of unit 1");
    pFits->pHDU().addKey("CRVAL1C", first_pix_val / con_pc, "value of axis 1");
    pFits->pHDU().addKey("CRPIX1C", 1, "pixel where CRVAL1 is defined ");
    pFits->pHDU().addKey("CDELT1C", bin_width / con_pc, "delta of axis 1");
    pFits->pHDU().addKey("CUNIT1C", "pc", "unit of axis 1");

    // Grid
    pFits->pHDU().addKey("CTYPE2", "PARAM", "type of unit 2");
    pFits->pHDU().addKey("CRVAL2", first_pix_val, "value of axis 2");
    pFits->pHDU().addKey("CRPIX2", 1, "pixel where CRVAL2 is defined ");
    pFits->pHDU().addKey("CDELT2", bin_width, "delta of axis 2");
    pFits->pHDU().addKey("CUNIT2", "m", "unit of axis 2");

    // Alternatively as AU grid
    pFits->pHDU().addKey("CTYPE2B", "PARAM", "type of unit 2");
    pFits->pHDU().addKey("CRVAL2B", first_pix_val / con_AU, "value of axis 2");
    pFits->pHDU().addKey("CRPIX2B", 1, "pixel where CRVAL2 is defined ");
    pFits->pHDU().addKey("CDELT2B", bin_width / con_AU, "delta of axis 2");
    pFits->pHDU().addKey("CUNIT2B", "AU", "unit of axis 2");

    // Alternatively as pc grid
    pFits->pHDU().addKey("CTYPE2C", "PARAM", "type of unit 2");
    pFits->pHDU().addKey("CRVAL2C", first_pix_val / con_pc, "value of axis 2");
    pFits->pHDU().addKey("CRPIX2C", 1, "pixel where CRVAL2 is defined ");
    pFits->pHDU().addKey("CDELT2C", bin_width / con_pc, "delta of axis 2");
    pFits->pHDU().addKey("CUNIT2C", "pc", "unit of axis 2");
    if(midplane_3d_param.size() == 4)
    {
        double bin_width_z = z_step;
        double first_pix_val_z = shift_z - (z_step * double(naxes[2])) / 2.0 + (bin_width_z / 2.0);

        // Grid
        pFits->pHDU().addKey("CTYPE3", "PARAM", "type of unit 3");
        pFits->pHDU().addKey("CRVAL3", first_pix_val_z, "value of axis 3");
        pFits->pHDU().addKey("CRPIX3", 1, "pixel where CRVAL3 is defined ");
        pFits->pHDU().addKey("CDELT3", bin_width_z, "delta of axis 3");
        pFits->pHDU().addKey("CUNIT3", "m", "unit of axis 3");

        // Alternatively as AU grid
        pFits->pHDU().addKey("CTYPE3B", "PARAM", "type of unit 3");
        pFits->pHDU().addKey("CRVAL3B", first_pix_val_z / con_AU, "value of axis 3");
        pFits->pHDU().addKey("CRPIX3B", 1, "pixel where CRVAL3 is defined ");
        pFits->pHDU().addKey("CDELT3B", bin_width_z / con_AU, "delta of axis 3");
        pFits->pHDU().addKey("CUNIT3B", "AU", "unit of axis 3");

        // Alternatively as pc grid
        pFits->pHDU().addKey("CTYPE3C", "PARAM", "type of unit 3");
        pFits->pHDU().addKey("CRVAL3C", first_pix_val_z / con_pc, "value of axis 3");
        pFits->pHDU().addKey("CRPIX3C", 1, "pixel where CRVAL3 is defined ");
        pFits->pHDU().addKey("CDELT3C", bin_width_z / con_pc, "delta of axis 3");
        pFits->pHDU().addKey("CUNIT3C", "pc", "unit of axis 3");

        // Quantities
        pFits->pHDU().addKey("CTYPE4", "PARAM", "type of unit 4");
        pFits->pHDU().addKey("CRVAL4", 1, "value of axis 4");
        pFits->pHDU().addKey("CRPIX4", 1, "pixel where CRVAL4 is defined ");
        pFits->pHDU().addKey("CDELT4", 1, "delta of axis 4");
        pFits->pHDU().addKey("CUNIT4", "see MIDPLANEX", "unit of axis 4");
    }

    uint counter = 0;
    char str_1[1024];
    char str_2[1024];
    if(plt_gas_dens)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        if(nr_densities > 1 && size_gd_list >= nr_densities)
        {
            if(gas_is_mass_density)
                pFits->pHDU().addKey(str_1, "total_gas_mass_density [kg/m^3]", str_2);
            else
                pFits->pHDU().addKey(str_1, "total_gas_number_density [m^-3]", str_2);
            for(uint i_density = 1; i_density <= nr_densities; i_density++)
            {
                counter++;
                updateMidplaneString(str_1, str_2, counter);
                string str_3;
                if(gas_is_mass_density)
                    str_3 = getDensityString("gas_mass_density_%i [kg/m^3]", i_density);
                else
                    str_3 = getDensityString("gas_number_density_%i [m^-3]", i_density);
                pFits->pHDU().addKey(str_1, str_3, str_2);
            }
        }
        else
        {
            if(gas_is_mass_density)
                pFits->pHDU().addKey(str_1, "gas_mass_density [kg/m^3]", str_2);
            else
                pFits->pHDU().addKey(str_1, "gas_number_density [m^-3]", str_2);
        }
    }
    if(plt_mol_dens)
    {
        updateMidplaneString(str_1, str_2, counter);

        for(uint i_density = 1; i_density <= nrOfDensRatios; i_density++)
        {
            counter++;
            updateMidplaneString(str_1, str_2, counter);
            string str_3;
            if(gas_is_mass_density)
                str_3 = getDensityString("mol_mass_density_%i [kg/m^3]", i_density);
            else
                str_3 = getDensityString("mol_number_density_%i [m^-3]", i_density);
            pFits->pHDU().addKey(str_1, str_3, str_2);
        }
    }
    if(plt_dust_dens)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        if(nr_densities > 1 && size_dd_list >= nr_densities)
        {
            if(dust_is_mass_density)
                pFits->pHDU().addKey(str_1, "total_dust_mass_density [kg/m^3]", str_2);
            else
                pFits->pHDU().addKey(str_1, "total_dust_number_density [m^-3]", str_2);
            for(uint i_density = 1; i_density <= nr_densities; i_density++)
            {
                counter++;
                updateMidplaneString(str_1, str_2, counter);
                string str_3;
                if(dust_is_mass_density)
                    str_3 = getDensityString("dust_mass_density_%i [kg/m^3]", i_density);
                else
                    str_3 = getDensityString("dust_number_density_%i [m^-3]", i_density);
                pFits->pHDU().addKey(str_1, str_3, str_2);
            }
        }
        else
        {
            if(dust_is_mass_density)
                pFits->pHDU().addKey(str_1, "dust_mass_density [kg/m^3]", str_2);
            else
                pFits->pHDU().addKey(str_1, "dust_number_density [m^-3]", str_2);
        }
    }
    if(plt_gas_temp)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "gas_temperature [K]", str_2);
    }
    if(plt_dust_temp)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        if(nr_densities > 1)
        {
            pFits->pHDU().addKey(str_1, "average_dust_temperature [K]", str_2);
            for(uint i_density = 1; i_density <= nr_densities; i_density++)
            {
                counter++;
                updateMidplaneString(str_1, str_2, counter);
                string str_3 = getDensityString("dust_temperature_%i [K]", i_density);
                pFits->pHDU().addKey(str_1, str_3, str_2);
            }
        }
        else
            pFits->pHDU().addKey(str_1, "dust_temperature [K]", str_2);
    }
    if(plt_rat)
    {
        if(nr_densities > 1)
            for(uint i_density = 1; i_density <= nr_densities; i_density++)
            {
                counter++;
                updateMidplaneString(str_1, str_2, counter);
                string str_3 = getDensityString("rat_aalig_%i [m]", i_density);
                pFits->pHDU().addKey(str_1, str_3, str_2);
            }
        else
        {
            counter++;
            updateMidplaneString(str_1, str_2, counter);
            pFits->pHDU().addKey(str_1, "rat_aalig [m]", str_2);
        }
    }
    /*if(plt_delta)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "delta [m]", str_2);
    }*/
    if(plt_mag)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "mag_total [T]", str_2);
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "mag_x [T]", str_2);
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "mag_y [T]", str_2);
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "mag_z [T]", str_2);
    }
    if(plt_vel)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "vel_total [m/s]", str_2);
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "vel_x [m/s]", str_2);
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "vel_y [m/s]", str_2);
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "vel_z [m/s]", str_2);
    }
    if(plt_larm)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "larm [m]", str_2);
    }
    
    /*if(plt_mach)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "mach number", str_2);
    }*/
    
    if(plt_dust_id)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "dust_mixture [index]", str_2);
    }
    if(plt_amin)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "a_min [m]", str_2);
    }
    if(plt_amax)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "a_max [m]", str_2);
    }
    if(plt_size_param)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "size param [value]", str_2);
    }
    if(plt_rad_field1)
    {
        for(uint i_comp = 0; i_comp < nr_rad_field_comp; i_comp++)
        {
            for(uint wID = 0; wID < WL_STEPS; wID++)
            {
                counter++;
                updateMidplaneString(str_1, str_2, counter);
                char str_3[1024];
                switch(i_comp)
                {
                    default:
#ifdef WINDOWS
                        sprintf_s(str_3, "rad_field [W/m/m^2] (%.3e [m])", wl_list[wID]);
#else
                        sprintf(str_3, "rad_field [W/m/m^2] (%.3e [m])", wl_list[wID]);
#endif
                        break;

                    case 1:
#ifdef WINDOWS
                        sprintf_s(str_3, "rad_field_x [W/m/m^2] (%.3e [m])", wl_list[wID]);
#else
                        sprintf(str_3, "rad_field_x [W/m/m^2] (%.3e [m])", wl_list[wID]);
#endif
                        break;

                    case 2:
#ifdef WINDOWS
                        sprintf_s(str_3, "rad_field_y [W/m/m^2] (%.3e [m])", wl_list[wID]);
#else
                        sprintf(str_3, "rad_field_y [W/m/m^2] (%.3e [m])", wl_list[wID]);
#endif
                        break;

                    case 3:
#ifdef WINDOWS
                        sprintf_s(str_3, "rad_field_z [W/m/m^2] (%.3e [m])", wl_list[wID]);
#else
                        sprintf(str_3, "rad_field_z [W/m/m^2] (%.3e [m])", wl_list[wID]);
#endif
                        break;
                }
                pFits->pHDU().addKey(str_1, string(str_3), str_2);
            }
        }
    }
    if(plt_g_zero1)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "G_0 (dustem)", str_2);
    }
    if(plt_u_rad)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "u_rad/u_isrf", str_2);
    }
    if(plt_n_th)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "therm_el_density [m^-3]", str_2);
    }
    if(plt_T_e)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "electron temperature [K]", str_2);
    }
    if(plt_n_cr)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "cr_el_density [m^-3]", str_2);
    }
    if(plt_g_min)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "gamma_min", str_2);
    }
    if(plt_g_max)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "gamma_max", str_2);
    }
    if(plt_p)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "syn_p", str_2);
    }
    if(plt_avg_th)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "avg. RAT cos(theta)", str_2);
    }
    if(plt_avg_dir)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "avg. RAT aniso. (gamma)", str_2);
    }


    cout << CLR_LINE;
    cout << "Memory cleanup of the plotting arrays ...     \r" << flush;
    // Free memory of pointer arrays
    if(plt_gas_dens)
    {
        for(long i_cell = 0; i_cell < nelements; i_cell++)
            delete[] buffer_gas_dens[i_cell];
        delete[] buffer_gas_dens;
    }
    if(plt_mol_dens)
    {
        for(long i_cell = 0; i_cell < nelements; i_cell++)
            delete[] buffer_mol_dens[i_cell];
        delete[] buffer_mol_dens;
    }
    if(plt_dust_dens)
    {
        for(long i_cell = 0; i_cell < nelements; i_cell++)
            delete[] buffer_dust_dens[i_cell];
        delete[] buffer_dust_dens;
    }
    if(plt_gas_temp)
        delete[] buffer_gas_temp;
    if(plt_dust_temp)
    {
        for(long i_cell = 0; i_cell < nelements; i_cell++)
            delete[] buffer_dust_temp[i_cell];
        delete[] buffer_dust_temp;
    }
    if(plt_rat)
        delete[] buffer_rat;
        
    /*if(plt_delta)
        delete[] buffer_delta;*/
        
    if(plt_mag)
    {
        delete[] buffer_mag;
        delete[] buffer_mag_x;
        delete[] buffer_mag_y;
        delete[] buffer_mag_z;
    }
    if(plt_vel)
    {
        delete[] buffer_vel;
        delete[] buffer_vel_x;
        delete[] buffer_vel_y;
        delete[] buffer_vel_z;
    }
    if(plt_larm)
        delete[] buffer_larm;
    //if(plt_mach)
    //    delete[] buffer_mach;
    if(plt_dust_id)
        delete[] buffer_dust_mixture;
    if(plt_amin)
        delete[] buffer_dust_amin;
    if(plt_amax)
        delete[] buffer_dust_amax;
    if(plt_size_param)
        delete[] buffer_dust_size_param;
    if(plt_rad_field1)
    {
        for(long i_cell = 0; i_cell < nelements; i_cell++)
        {
            for(uint wID = 0; wID < WL_STEPS; wID++)
                delete[] buffer_rad_field[i_cell][wID];
            delete[] buffer_rad_field[i_cell];
        }
        delete[] buffer_rad_field;
    }
    if(plt_g_zero1)
        delete[] buffer_g_zero1;
    if(plt_u_rad)
        delete[] buffer_u_rad;
    if(plt_n_th)
        delete[] buffer_n_th;
    if(plt_T_e)
        delete[] buffer_T_e;
    if(plt_n_cr)
        delete[] buffer_n_cr;
    if(plt_g_min)
        delete[] buffer_g_min;
    if(plt_g_max)
        delete[] buffer_g_max;
    if(plt_p)
        delete[] buffer_p;
    if(plt_avg_th)
        delete[] buffer_avg_th;
    if(plt_avg_dir)
        delete[] buffer_avg_dir;

    cout << CLR_LINE;
    cout << "- Writing of midplane files     : done" << endl;

    return res;
}

bool CGridBasic::getPolarRTGridParameterWorker(double max_len,
                           double pixel_width,
                           uint max_subpixel_lvl,
                           dlist & _listR,
                           uint & N_polar_r,
                           uint *& N_polar_ph,
                           const uint &N_r,
                           const double *listR
    )
{
    uint subpixel_multiplier = pow(2, max_subpixel_lvl);

    // Add polar detector pixels in the inner region to obtain resolution specified by max_subpixel_lvl
    // inner grid cell diameter is 2.*listR[0]
    uint N_r_inner = uint(ceil(subpixel_multiplier * 2.0 * listR[0] / pixel_width)); 

    for(uint i_r = 0; i_r <= N_r_inner; i_r++)
        _listR.push_back(listR[0] * (i_r / double(N_r_inner)));

    for(uint i_r = 1; i_r <= N_r; i_r++)
    {
        double r1 = _listR[_listR.size() - 1];
        double r2 = listR[i_r];

        // if sidelength is smaller than full grid, only consider visible grid
        if(r2 > max_len)
            r2 = max_len;

        // r2 - r1 is width of current grid cell's ring
        uint N_r_sub = uint(ceil(subpixel_multiplier * (r2 - r1) / pixel_width));

        for(uint i_r_sub = 1; i_r_sub <= N_r_sub; i_r_sub++)
            _listR.push_back(r1 + (r2 - r1) * i_r_sub / double(N_r_sub));

        // break if sidelength is smaller than full grid
        if(_listR.back() >= max_len)
        {
            _listR.pop_back();
            _listR.push_back(max_len);
            break;
        }
    }

    if(_listR.back() < max_len)
    {
        // Create additional outer rings with outermost grid cell radial distance
        // and store them in buffer to do subpixeling afterwards
        uint N_r_outer = uint(ceil((max_len - listR[N_r]) / (listR[N_r] - listR[N_r - 1])));
        std::vector<double> outerR_buffer;
            
        for(uint i_r = 1; i_r <= N_r_outer; i_r++)
            outerR_buffer.push_back(listR[N_r] + (max_len - listR[N_r]) * i_r / double(N_r_outer));

        // loop over equally spaced rings outside the grid and do subpixeling
        for (uint i_r = 0; i_r < N_r_outer; i_r++)
        {
            double r1 = _listR[_listR.size() - 1];
            double r2 = outerR_buffer[i_r];
            uint N_r_sub = uint(ceil(subpixel_multiplier * (r2 - r1) / pixel_width));  

            for(uint i_r_sub = 1; i_r_sub <= N_r_sub; i_r_sub++)   
                _listR.push_back(r1 + (r2 - r1) * i_r_sub / double(N_r_sub));
        }
    }

    // Set total size of the radial cells
    N_polar_r = _listR.size() - 1;

    // Compute the number of phi background grid pixel
    N_polar_ph = new uint[N_polar_r];
    for(uint i_r = 0; i_r < N_polar_r; i_r++)
    {
        N_polar_ph[i_r] = uint(ceil(PIx2 * _listR[i_r + 1] / (_listR[i_r + 1] - _listR[i_r])));
    }

    return true;
}

double CGridBasic::getTurbulentVelocity(cell_basic * cell)
{
    if(turbulent_velocity > 0)
        return turbulent_velocity;
    else if(hasTurbulentVelocity())
        return cell->getData(data_pos_vt);
    return 0;
}

double CGridBasic::getTurbulentVelocity(photon_package * pp)
{
    return getTurbulentVelocity(pp->getPositionCell());
}

void CGridBasic::updateVelocity(cell_basic * cell, parameters & param)
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

uint CGridBasic::getDataOffset()
{
    return data_offset;
}

uint CGridBasic::getDataID()
{
    return dataID;
}

bool CGridBasic::hasVelocityField()
{
    return (data_pos_vx != MAX_UINT && data_pos_vy != MAX_UINT && data_pos_vz != MAX_UINT);
}

bool CGridBasic::hasTurbulentVelocity()
{
    return (data_pos_vt != MAX_UINT);
}

Vector3D CGridBasic::getCenter(const photon_package & pp) const
{
    return getCenter(*pp.getPositionCell());
}

uint CGridBasic::getDataLength()
{
    return data_len;
}

ulong CGridBasic::getMaxDataCells()
{
    return max_cells;
}

uint CGridBasic::getDataSize()
{
    return max_data;
}

bool CGridBasic::updateShortestDistance(photon_package * pp)
{
    return false;
}

void CGridBasic::setOrientation(Vector3D n1, Vector3D n2, double _rot_angle1, double _rot_angle2)
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

void CGridBasic::getMagFieldInfo(const photon_package & pp, MagFieldInfo * mfo) const
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

double CGridBasic::getMinLength()
{
    return min_len;
}

double CGridBasic::getMaxLength()
{
    return max_len;
}

bool CGridBasic::getPolarRTGridParameter(double max_len,
                                        double pixel_width,
                                        uint max_subpixel_lvl,
                                        dlist & _listR,
                                        uint & N_polar_r,
                                        uint *& N_polar_ph)
{
    return false;
}

cell_basic * CGridBasic::getCellFromIndex(ulong i)
{
    return cell_list[i];
}

void CGridBasic::setSIConversionFactors(parameters & param)
{
    mu = param.getMu();
    conv_length_in_SI = param.getSIConvLength();

    delta0 = param.getDelta0();
    larm_f = param.getLarmF();

    conv_dens_in_SI = abs(param.getSIConvDH());
    conv_Bfield_in_SI = param.getSIConvBField();
    conv_Vfield_in_SI = param.getSIConvVField();
}

void CGridBasic::setDataSize(uint sz)
{
    max_data = data_offset + sz;
}

void CGridBasic::setDustInformation(uint _nr_mixtures,
                        uint * _nr_dust_temp_sizes,
                        uint * _nr_stochastic_sizes,
                        uint * _nr_stochastic_temps)
{
    nr_mixtures = _nr_mixtures;
    nr_dust_temp_sizes = _nr_dust_temp_sizes;
    nr_stochastic_sizes = _nr_stochastic_sizes;
    nr_stochastic_temps = _nr_stochastic_temps;
}

void CGridBasic::setGasInformation(uint ** _level_to_pos, uint *** _line_to_pos)
{
    level_to_pos = _level_to_pos;
    line_to_pos = _line_to_pos;
}

void CGridBasic::setVelocityFieldNeeded(bool val)
{
    velocity_field_needed = val;
}

void CGridBasic::setDataOffset(uint off)
{
    data_offset = off;
}

void CGridBasic::setDataID(uint id)
{
    dataID = id;
}

void CGridBasic::setGasDensity(photon_package * pp, double dens)
{
    setGasDensity(pp->getPositionCell(), dens);
}

void CGridBasic::setGasDensity(photon_package * pp, uint i_density, double dens)
{
    setGasDensity(pp->getPositionCell(), i_density, dens);
}

void CGridBasic::setSpecLengthAsVector(bool val)
{
    spec_length_as_vector = val;
}

bool CGridBasic::specLengthIsVector()
{
    return spec_length_as_vector;
}

void CGridBasic::updateSpecLength(photon_package * pp, double len)
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

void CGridBasic::updateSpecLength(cell_basic * cell, uint i_offset, StokesVector stokes) const
{
    stokes /= getVolume(*cell);
    uint data_pos = data_offset + 4 * i_offset;
    cell->updateData(data_pos + 0, stokes.I());
    cell->updateData(data_pos + 1, stokes.Q());
    cell->updateData(data_pos + 2, stokes.U());
    cell->updateData(data_pos + 3, stokes.V());
}

double CGridBasic::getSpecLength(const cell_basic & cell, uint wID) const
{
#if BENCHMARK == CAMPS
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

double CGridBasic::getSpecLength(const photon_package & pp, uint wID) const
{
    return getSpecLength(*pp.getPositionCell(), wID);
}

void CGridBasic::getSpecLength(const cell_basic & cell, uint wID, double * us, Vector3D * e_dir) const
{
    uint data_pos = data_offset + 4 * wID;

    *us = cell.getData(data_pos + 0);
    e_dir->setX(cell.getData(data_pos + 1));
    e_dir->setY(cell.getData(data_pos + 2));
    e_dir->setZ(cell.getData(data_pos + 3));
}

void CGridBasic::saveRadiationField()
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

double CGridBasic::getRadiationField(const cell_basic & cell, uint wID) const
{
#if BENCHMARK == CAMPS
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

double CGridBasic::getRadiationField(const photon_package & pp, uint wID) const
{
    return getRadiationField(*pp.getPositionCell(), wID);
}

double CGridBasic::getRadiationFieldX(const cell_basic & cell, uint wID) const
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

double CGridBasic::getRadiationFieldX(const photon_package & pp, uint wID) const
{
    return getRadiationFieldX(*pp.getPositionCell(), wID);
}

double CGridBasic::getRadiationFieldY(const cell_basic & cell, uint wID) const
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

double CGridBasic::getRadiationFieldY(const photon_package & pp, uint wID) const
{
    return getRadiationFieldY(*pp.getPositionCell(), wID);
}

double CGridBasic::getRadiationFieldZ(const cell_basic & cell, uint wID) const
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

double CGridBasic::getRadiationFieldZ(const photon_package & pp, uint wID) const
{
    return getRadiationFieldZ(*pp.getPositionCell(), wID);
}

void CGridBasic::getRadiationField(const photon_package & pp, uint w, double * us, Vector3D * e_dir) const
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
            cout << ERROR_LINE << "This should not happen" << endl;
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

StokesVector CGridBasic::getStokesFromRadiationField(const photon_package & pp, uint i_offset) const
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
        Vector3D::getAnglePhi(pp.getEX(), pp.getEY(), rot_dir) - Vector3D::getAnglePhi(pp.getEX(), pp.getEY(), getCenter(pp));

    // Rotate Stokes Vector to be in agreement with the detector plane
    scattering_stokes.rot(phi_map);

    return scattering_stokes;
}

void CGridBasic::getRadiationFieldInterp(const photon_package & pp,
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

double CGridBasic::getGZero(const cell_basic & cell) const
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

double CGridBasic::getGZero(const photon_package & pp) const
{
    return getGZero(*pp.getPositionCell());
}

double CGridBasic::getUrad(const cell_basic & cell) const
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

double CGridBasic::getUrad(const photon_package & pp) const
{
    return getUrad(*pp.getPositionCell());
}

double CGridBasic::getMu() const
{
    return mu;
}

void CGridBasic::setRelOutsidePosition(photon_package * pp, double tx, double ty, double tz)
{
    pp->setPosition(tx * ex + ty * ey + tz * ez);
}

void CGridBasic::setRelOutsidePosition(photon_package * pp, double tx, double ty)
{
    pp->setPosition(tx * ex + ty * ey - max_len * ez);
}

void CGridBasic::setRelDirection(photon_package * pp)
{
    pp->setEX(ex);
    pp->setEY(ey);
    pp->setEZ(ez);
}

void CGridBasic::setRndPositionInCell(photon_package * pp, CRandomGenerator * rand_gen)
{
    pp->setPosition(Vector3D(0, 0, 0));
}

Vector3D CGridBasic::rotateToCenter(const photon_package & pp, bool inv, bool phi_only) const
{
    return rotateToCenter(pp, pp.getDirection(), inv, phi_only);
}

Vector3D CGridBasic::rotateToCenter(const photon_package & pp,
                                Vector3D dir,
                                bool inv,
                                bool phi_only) const
{
    return dir;
}

void CGridBasic::setDustTemperature(cell_basic * cell, uint i_density, uint a, double temp)
{
    if(!data_pos_dt_list.empty())
    {
        uint id = a + nr_densities;
        for(uint i = 0; i < i_density; i++)
            id += size_skip[i];
        cell->setData(data_pos_dt_list[id], temp);
    }
}

void CGridBasic::setDustTemperature(cell_basic * cell, uint i_density, double temp)
{
    cell->setData(data_pos_dt_list[i_density], temp);
}

void CGridBasic::setDustTemperature(cell_basic * cell, double temp)
{
    for(uint i_density = 0; i_density < data_pos_dt_list.size(); i_density++)
        cell->setData(data_pos_dt_list[i_density], temp);
}

void CGridBasic::setDustTemperature(photon_package * pp, uint i_density, uint a, double temp)
{
    setDustTemperature(pp->getPositionCell(), i_density, a, temp);
}

void CGridBasic::setDustTemperature(photon_package * pp, uint i_density, double temp)
{
    setDustTemperature(pp->getPositionCell(), i_density, temp);
}

void CGridBasic::setDustTemperature(photon_package * pp, double temp)
{
    setDustTemperature(pp->getPositionCell(), temp);
}

void CGridBasic::setDustTempProbability(cell_basic * cell, uint i_density, uint a, uint t, double temp)
{
    uint id = a * nr_stochastic_temps[i_density] + t;
    for(uint i = 0; i < i_density; i++)
        id += nr_stochastic_sizes[i] * nr_stochastic_temps[i];
    cell->setData(data_offset + id, temp);
}

void CGridBasic::setPDAValue(cell_basic * cell, double val)
{
    cell->setData(data_pos_pda, val);
}

double CGridBasic::getMagMax()
{
    return max_mag;
}

void CGridBasic::setGasTemperature(cell_basic * cell, double temp)
{
    cell->setData(data_pos_tg, temp);
}

void CGridBasic::setElectronTemperature(cell_basic * cell, double temp)
{
    if(data_pos_T_e != MAX_UINT)
        cell->setData(data_pos_T_e, temp);
}

void CGridBasic::setThermalElectronDensity(cell_basic * cell, double dens)
{
    if(data_pos_n_th != MAX_UINT)
        cell->setData(data_pos_n_th, dens);
}

void CGridBasic::setCRElectronDensity(cell_basic * cell, double dens)
{
    if(data_pos_n_cr != MAX_UINT)
        cell->setData(data_pos_n_cr, dens);
}

void CGridBasic::setGammaMin(cell_basic * cell, double g_min)
{
    if(data_pos_g_min != MAX_UINT)
        cell->setData(data_pos_g_min, g_min);
}

void CGridBasic::setGammaMax(cell_basic * cell, double g_max)
{
    if(data_pos_g_max != MAX_UINT)
        cell->setData(data_pos_g_max, g_max);
}

void CGridBasic::setPowerLawIndex(cell_basic * cell, double p)
{
    if(data_pos_p != MAX_UINT)
        cell->setData(data_pos_p, p);
}

void CGridBasic::setAvgTheta(cell_basic * cell, double phi)
{
    if(data_pos_avg_th != MAX_UINT)
        cell->setData(data_pos_avg_th, phi);
}

void CGridBasic::setAvgDir(cell_basic * cell, double dir)
{
    if(data_pos_avg_dir != MAX_UINT)
        cell->setData(data_pos_avg_dir, dir);
}

void CGridBasic::setAvg_ux(cell_basic * cell, double ux)
{
    if(data_pos_avg_ux != MAX_UINT)
        cell->setData(data_pos_avg_ux, ux);
}

void CGridBasic::setAvg_uy(cell_basic * cell, double uy)
{
    if(data_pos_avg_uy != MAX_UINT)
        cell->setData(data_pos_avg_uy, uy);
}

void CGridBasic::setAvg_uz(cell_basic * cell, double uz)
{
    if(data_pos_avg_uz != MAX_UINT)
        cell->setData(data_pos_avg_uz, uz);
}

void CGridBasic::setIonDensity(cell_basic * cell, double n_i)
{
    if(data_pos_n_i != MAX_UINT)
        cell->setData(data_pos_n_i, n_i);
}

void CGridBasic::setIonCharge(cell_basic * cell, double Z)
{
    if(data_pos_Z != MAX_UINT)
        cell->setData(data_pos_Z, Z);
}

void CGridBasic::setDustChoiceID(cell_basic * cell, uint dust_id)
{
    if(data_pos_id != MAX_UINT)
        cell->setData(data_pos_id, dust_id);
}

void CGridBasic::setGasDensity(cell_basic * cell, double dens)
{
    cell->setData(data_pos_gd_list[0], dens);
}

void CGridBasic::setGasDensity(cell_basic * cell, uint i_density, double dens)
{
    cell->setData(data_pos_gd_list[i_density], dens);
}

double CGridBasic::getQBOffset(const cell_basic & cell, uint i_density) const
{
    if(data_pos_dt_list.size() == 1)
        return cell.getData(data_pos_dt_list[0]);
    else if(data_pos_dt_list.size() > i_density)
        return cell.getData(data_pos_dt_list[i_density]);
    else
        return 0;
}

double CGridBasic::getQBOffset(const photon_package & pp, uint i_density) const
{
    return getQBOffset(*pp.getPositionCell(), i_density);
}

double CGridBasic::getQBOffset(const cell_basic & cell, uint i_density, uint a) const
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

double CGridBasic::getQBOffset(const photon_package & pp, uint i_density, uint a) const
{
    return getQBOffset(*pp.getPositionCell(), i_density, a);
}

void CGridBasic::setQBOffset(cell_basic * cell, uint i_density, uint a, double temp)
{
    if(!data_pos_dt_list.empty())
    {
        uint id = a + nr_densities;
        for(uint i = 0; i < i_density; i++)
            id += size_skip[i];
        cell->setData(data_pos_dt_list[id], temp);
    }
}

void CGridBasic::setQBOffset(cell_basic * cell, uint i_density, double temp)
{
    if(data_pos_dt_list.size() == 1)
        cell->setData(data_pos_dt_list[0], temp);
    else if(data_pos_dt_list.size() > 1)
        cell->setData(data_pos_dt_list[i_density], temp);
}

uint CGridBasic::getNrAlignedRadii()
{
    return data_pos_aalg_list.size();
}

double CGridBasic::getAlignedRadius(const cell_basic & cell, uint i_density) const
{
    if(data_pos_aalg_list.size() == 1)
        return cell.getData(data_pos_aalg_list[0]);
    else if(data_pos_aalg_list.size() > i_density)
        return cell.getData(data_pos_aalg_list[i_density]);
    else
        return 0;
}

double CGridBasic::getAlignedRadius(const photon_package & pp, uint i_density) const
{
    return getAlignedRadius(*pp.getPositionCell(), i_density);
}

void CGridBasic::setAlignedRadius(cell_basic * cell, uint i_density, double _a_alg)
{
    cell->setData(data_pos_aalg_list[i_density], _a_alg);
}

double CGridBasic::getMinGrainRadius(const cell_basic & cell, uint i_density) const
{
    if(data_pos_amin_list.size() == 1)
        return cell.getData(data_pos_amin_list[0]);
    else if(data_pos_amin_list.size() > i_density)
        return cell.getData(data_pos_amin_list[i_density]);
    else
        return 0;
}

double CGridBasic::getMinGrainRadius(const photon_package & pp, uint i_density) const
{
    return getMinGrainRadius(*pp.getPositionCell(), i_density);
}

double CGridBasic::getMaxGrainRadius(const cell_basic & cell, uint i_density) const
{
    if(data_pos_amax_list.size() == 1)
        return cell.getData(data_pos_amax_list[0]);
    else if(data_pos_amax_list.size() > i_density)
        return cell.getData(data_pos_amax_list[i_density]);
    else
        return 0;
}

double CGridBasic::getMaxGrainRadius(const photon_package & pp, uint i_density) const
{
    return getMaxGrainRadius(*pp.getPositionCell(), i_density);
}

double CGridBasic::getGrainSizeParam(const cell_basic & cell, uint i_density) const
{
    if(data_pos_size_param_list.size() == 1)
        return cell.getData(data_pos_size_param_list[0]);
    else if(data_pos_size_param_list.size() > i_density)
        return cell.getData(data_pos_size_param_list[i_density]);
    else
        return 0;
}

double CGridBasic::getGrainSizeParam(const photon_package & pp, uint i_density) const
{
    return getGrainSizeParam(*pp.getPositionCell(), i_density);
}

uint CGridBasic::getDustChoiceID(const photon_package & pp) const
{
    return getDustChoiceID(*pp.getPositionCell());
}

uint CGridBasic::getDustChoiceID(const cell_basic & cell) const
{
    if(data_pos_id != MAX_UINT)
        return uint(cell.getData(data_pos_id));
    else
        return 0;
}

void CGridBasic::getLineBroadening(const photon_package & pp, uint i_trans, LineBroadening * line_broadening) const
{
    line_broadening->gauss_a = getGaussA(pp);
    line_broadening->voigt_a = getVoigtA(pp, i_trans);
}

double CGridBasic::getGaussA(const cell_basic & cell) const
{
    return cell.getData(data_offset);
}

double CGridBasic::getGaussA(const photon_package & pp) const
{
    return getGaussA(*pp.getPositionCell());
}

double CGridBasic::getVoigtA(const cell_basic & cell, uint i_line) const
{
    return cell.getData(data_offset + 1 + i_line);
}

double CGridBasic::getVoigtA(const photon_package & pp, uint i_line) const
{
    return getVoigtA(*pp.getPositionCell(), i_line);
}

double CGridBasic::getLvlPop(const cell_basic & cell, uint i_lvl, uint i_sublvl) const
{
    if(level_to_pos[i_lvl][i_sublvl] != MAX_UINT)
        return cell.getData(data_offset + level_to_pos[i_lvl][i_sublvl]);

    return 0;
}

double CGridBasic::getLvlPop(const photon_package & pp, uint i_lvl, uint i_sublvl) const
{
    return getLvlPop(*pp.getPositionCell(), i_lvl, i_sublvl);
}

void CGridBasic::setVelocityField(cell_basic * cell, const Vector3D & vel)
{
    cell->setData(data_pos_vx, vel.X());
    cell->setData(data_pos_vy, vel.Y());
    cell->setData(data_pos_vz, vel.Z());
}

void CGridBasic::setVelocityField(photon_package * pp, const Vector3D & vel)
{
    pp->getPositionCell()->setData(data_pos_vx, vel.X());
    pp->getPositionCell()->setData(data_pos_vy, vel.Y());
    pp->getPositionCell()->setData(data_pos_vz, vel.Z());
}

uint CGridBasic::getCellID(cell_basic * cell)
{
    return cell->getUniqueID();
}

void CGridBasic::setLvlPopLower(cell_basic * cell, uint i_line, uint i_sublvl, double lvl_lower)
{
    cell->setData(data_offset + line_to_pos[i_line][0][i_sublvl], lvl_lower);
}

void CGridBasic::setLvlPopLower(photon_package * pp, uint i_line, uint i_sublvl, double lvl_lower)
{
    return setLvlPopLower(pp->getPositionCell(), i_line, i_sublvl, lvl_lower);
}

void CGridBasic::setLvlPopUpper(cell_basic * cell, uint i_line, uint i_sublvl, double lvl_upper)
{
    cell->setData(data_offset + line_to_pos[i_line][1][i_sublvl], lvl_upper);
}

void CGridBasic::setLvlPopUpper(photon_package * pp, uint i_line, uint i_sublvl, double lvl_upper)
{
    return setLvlPopUpper(pp->getPositionCell(), i_line, i_sublvl, lvl_upper);
}

void CGridBasic::setLvlPop(cell_basic * cell, uint i_lvl, uint i_sublvl, double lvl_pop)
{
    if(level_to_pos[i_lvl][i_sublvl] != MAX_UINT)
        cell->setData(data_offset + level_to_pos[i_lvl][i_sublvl], lvl_pop);
}

void CGridBasic::setLvlPop(photon_package * pp, uint i_lvl, uint i_sublvl, double lvl_pop)
{
    return setLvlPop(pp->getPositionCell(), i_lvl, i_sublvl, lvl_pop);
}

void CGridBasic::addToZeemanSublevel(cell_basic * cell, uint i_lvl, uint i_sublvl, double lvl_pop)
{
    if(level_to_pos[i_lvl][1 + i_sublvl] != MAX_UINT)
        cell->setData(data_offset + level_to_pos[i_lvl][1 + i_sublvl], lvl_pop);
}

void CGridBasic::addToZeemanSublevel(photon_package * pp, uint i_lvl, uint i_sublvl, double lvl_pop)
{
    return addToZeemanSublevel(pp->getPositionCell(), i_lvl, i_sublvl, lvl_pop);
}

void CGridBasic::setLineBroadening(cell_basic * cell, uint i_line_broad, const LineBroadening & line_broadening)
{
    cell->setData(data_offset + 1 + i_line_broad, line_broadening.voigt_a);
}

void CGridBasic::setGaussA(cell_basic * cell, double gauss_a)
{
    // Gauss_a only has to be set once
    cell->setData(data_offset, gauss_a);
}

uint CGridBasic::getOpiateID(cell_basic * cell)
{
    uint pos = pos_OpiateIDS[0];
    return uint(cell->getData(pos));
}

uint CGridBasic::getOpiateID(const photon_package * pp)
{
    uint pos = pos_OpiateIDS[0];
    return uint(pp->getPositionCell()->getData(pos));
}

uint CGridBasic::getOpiateID(photon_package pp)
{
    uint pos = pos_OpiateIDS[0];
    return uint(pp.getPositionCell()->getData(pos));
}

void CGridBasic::setOpiateID(cell_basic * cell, uint id)
{
    uint pos = pos_OpiateIDS[0];
    cell->setData(pos, double(id));
}

double CGridBasic::getOpiateTestData(photon_package * pp)
{
    uint pos = pos_OpiateIDS[1];

    return pp->getPositionCell()->getData(pos);
}

double CGridBasic::getOpiateTestData(cell_basic * cell)
{
    uint pos = pos_OpiateIDS[1];
    return uint(cell->getData(pos));
}

void CGridBasic::setOpiateTestData(cell_basic * cell, double val)
{
    uint pos = pos_OpiateIDS[1];
    cell->setData(pos, val);
}

double CGridBasic::getElectronTemperature(const photon_package & pp) const
{
    return getElectronTemperature(*pp.getPositionCell());
}

double CGridBasic::getElectronTemperature(const cell_basic & cell) const
{
    if(data_pos_T_e != MAX_UINT)
        return cell.getData(data_pos_T_e);

    return 0;
}

double CGridBasic::getThermalElectronDensity(const photon_package & pp) const
{
    return getThermalElectronDensity(*pp.getPositionCell());
}

double CGridBasic::getThermalElectronDensity(const cell_basic & cell) const
{
    if(data_pos_n_th != MAX_UINT)
        return cell.getData(data_pos_n_th);

    return 0;
}

double CGridBasic::getCRElectronDensity(const photon_package & pp) const
{
    return getCRElectronDensity(*pp.getPositionCell());
}

double CGridBasic::getCRElectronDensity(const cell_basic & cell) const
{
    if(data_pos_n_cr != MAX_UINT)
        return cell.getData(data_pos_n_cr);

    return 0;
}

double CGridBasic::getGammaMin(const photon_package & pp) const
{
    return getGammaMin(*pp.getPositionCell());
}

double CGridBasic::getGammaMin(const cell_basic & cell) const
{
    if(data_pos_g_min != MAX_UINT)
        return cell.getData(data_pos_g_min);

    return 0;
}

double CGridBasic::getGammaMax(const photon_package & pp) const
{
    return getGammaMax(*pp.getPositionCell());
}

double CGridBasic::getGammaMax(const cell_basic & cell) const
{
    if(data_pos_g_max != MAX_UINT)
        return cell.getData(data_pos_g_max);

    return 0;
}

double CGridBasic::getPowerLawIndex(const photon_package & pp) const
{
    return getPowerLawIndex(*pp.getPositionCell());
}

double CGridBasic::getPowerLawIndex(const cell_basic & cell) const
{
    if(data_pos_p != MAX_UINT)
        return cell.getData(data_pos_p);

    return 0;
}

Vector3D CGridBasic::getAvg_u(const photon_package & pp) const
{
    return getAvg_u(*pp.getPositionCell());
}

Vector3D CGridBasic::getAvg_u(const cell_basic & cell) const
{
    double ux=0, uy=0, uz=0;

    if(data_pos_avg_ux != MAX_UINT)
        ux=cell.getData(data_pos_avg_ux);
    else
        return Vector3D(0,0,0);

    if(data_pos_avg_uy != MAX_UINT)
        uy=cell.getData(data_pos_avg_uy);
    else
        return Vector3D(0,0,0);    

    if(data_pos_avg_uz != MAX_UINT)
        uz=cell.getData(data_pos_avg_uz);
    else
        return Vector3D(0,0,0);    

    return Vector3D(ux,uy,uz);
}

double CGridBasic::getAvg_ux(const cell_basic & cell) const
{
    if(data_pos_avg_ux != MAX_UINT)
        return cell.getData(data_pos_avg_ux);

    return 0;
}

double CGridBasic::getAvg_ux(const photon_package & pp) const
{
    return getAvg_ux(*pp.getPositionCell());
}

double CGridBasic::getAvg_uy(const cell_basic & cell) const
{
    if(data_pos_avg_uy != MAX_UINT)
        return cell.getData(data_pos_avg_uy);

    return 0;
}

double CGridBasic::getAvg_uy(const photon_package & pp) const
{
    return getAvg_uy(*pp.getPositionCell());
}

double CGridBasic::getAvg_uz(const cell_basic & cell) const
{
    if(data_pos_avg_uz != MAX_UINT)
        return cell.getData(data_pos_avg_uz);

    return 0;
}

double CGridBasic::getAvg_uz(const photon_package & pp) const
{
    return getAvg_uz(*pp.getPositionCell());
}

double CGridBasic::getAvgTheta(const photon_package & pp) const
{
    return getAvgTheta(*pp.getPositionCell());
}

double CGridBasic::getAvgTheta(const cell_basic & cell) const
{
    if(data_pos_avg_th != MAX_UINT)
        return cell.getData(data_pos_avg_th);

    return 0;
}

double CGridBasic::getIonDensity(const photon_package & pp) const
{
    return getIonDensity(*pp.getPositionCell());
}

double CGridBasic::getIonDensity(const cell_basic & cell) const
{
    if(data_pos_n_i != MAX_UINT)
        return cell.getData(data_pos_n_i);

    return 0;
}

double CGridBasic::getIonCharge(const photon_package & pp) const
{
    return getIonCharge(*pp.getPositionCell());
}

double CGridBasic::getIonCharge(const cell_basic & cell) const
{
    if(data_pos_Z != MAX_UINT)
        return cell.getData(data_pos_Z);

    return 0;
}


double CGridBasic::getAvgDir(const photon_package & pp) const
{
    return getAvgDir(*pp.getPositionCell());
}

double CGridBasic::getAvgDir(const cell_basic & cell) const
{
    if(data_pos_avg_dir != MAX_UINT)
        return cell.getData(data_pos_avg_dir);

    return 0;
}

double CGridBasic::getDustTemperature(const cell_basic & cell, uint i_density, uint a) const
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

double CGridBasic::getDustTemperature(const photon_package & pp, uint i_density, uint a) const
{
    return getDustTemperature(*pp.getPositionCell(), i_density, a);
}

double CGridBasic::getDustTemperature(const cell_basic & cell, uint i_density) const
{
    if(data_pos_dt_list.size() == 1)
        return cell.getData(data_pos_dt_list[0]);
    else if(data_pos_dt_list.size() > i_density)
        return cell.getData(data_pos_dt_list[i_density]);
    else
        return 0;
}

double CGridBasic::getDustTemperature(const photon_package & pp, uint i_density) const
{
    return getDustTemperature(*pp.getPositionCell(), i_density);
}

double CGridBasic::getDustTemperature(const cell_basic & cell) const
{
    double sum = 0;
    for(uint i_density = 0; i_density < nr_densities; i_density++)
        sum += getDustTemperature(cell, i_density) * getRelativeDustDensity(cell, i_density);
    return sum;
}

double CGridBasic::getDustTemperature(const photon_package & pp) const
{
    return getDustTemperature(*pp.getPositionCell());
}

double CGridBasic::getDustTempProbability(const cell_basic & cell, uint i_density, uint a, uint t) const
{
    uint id = a * nr_stochastic_temps[i_density] + t;
    for(uint i = 0; i < i_density; i++)
        id += nr_stochastic_sizes[i] * nr_stochastic_temps[i];
    return cell.getData(data_offset + id);
}

double CGridBasic::getDustTempProbability(const photon_package & pp, uint i_density, uint a, uint t) const
{
    return getDustTempProbability(*pp.getPositionCell(), i_density, a, t);
}

double CGridBasic::getPDAValue(const cell_basic & cell) const
{
    return cell.getData(data_pos_pda);
}

double CGridBasic::getGasTemperature(const photon_package & pp) const
{
    return getGasTemperature(*pp.getPositionCell());
}

double CGridBasic::getGasTemperature(const cell_basic & cell) const
{
    return max(TEMP_MIN, cell.getData(data_pos_tg));
}

void CGridBasic::setPlaneParameter(uint plane_index,
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

void CGridBasic::fillMidplaneBuffer(double tx, double ty, double tz, uint i_cell)
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
        
        /*if(plt_delta)
        {
            double field = getMagField(pp).length();
            double Td = getDustTemperature(pp);
            double Tg = getGasTemperature(pp);
            double dens = getGasDensity(pp);
            double delta = CMathFunctions::calc_delta(field, Td, Tg, dens);
            buffer_delta[i_cell] = delta;
        }*/
        
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
        /*if(plt_mach)
        {
            Vector3D vel_field = getVelocityField(pp);
            double gas_temp = getGasTemperature(pp);
            double mach = vel_field.length() / sqrt(con_kB * gas_temp / (mu * m_H));
            buffer_mach[i_cell] = mach;
        }*/
        if(plt_dust_id)
            buffer_dust_mixture[i_cell] = getDustChoiceID(pp);
        if(plt_amin)
            buffer_dust_amin[i_cell] = getMinGrainRadius(pp,0);
        if(plt_amax)
            buffer_dust_amax[i_cell] = getMaxGrainRadius(pp,0);
        if(plt_size_param)
            buffer_dust_size_param[i_cell] = getGrainSizeParam(pp,0);
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
        
        //if(plt_delta)
        //    buffer_delta[i_cell] = 0;
        
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
        //if(plt_mach)
        //    buffer_mach[i_cell] = 0;
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

void CGridBasic::updateMidplaneString(char * str_1, char * str_2, uint counter)
{
#ifdef WINDOWS
    sprintf_s(str_1, "MIDPLANE%i", counter);
    sprintf_s(str_2, "quantity of %i. image", counter);
#else
    sprintf(str_1, "MIDPLANE%i", counter);
    sprintf(str_2, "quantity of %i. image", counter);
#endif
}

string CGridBasic::getDensityString(string quantity, uint counter)
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

double CGridBasic::getVolume(const photon_package & pp) const
{
    return getVolume(*pp.getPositionCell());
}

double CGridBasic::getGasDensity(const cell_basic & cell) const
{
    double sum = 0;
    for(uint i_density = 0; i_density < size_gd_list; i_density++)
        sum += getGasDensity(cell, i_density);
    return sum;
}

double CGridBasic::getGasDensity(const cell_basic & cell, uint i_density) const
{
    if(size_gd_list > i_density)
        return cell.getData(data_pos_gd_list[i_density]);
    else
        return 0;
}

double CGridBasic::getGasDensity(const photon_package & pp) const
{
    return getGasDensity(*pp.getPositionCell());
}

double CGridBasic::getGasDensity(const photon_package & pp, uint i_density) const
{
    return getGasDensity(*pp.getPositionCell(), i_density);
}

double CGridBasic::getGasNumberDensity(const cell_basic & cell) const
{
    double sum = 0;
    for(uint i_density = 0; i_density < size_gd_list; i_density++)
        sum += cell.getData(data_pos_gd_list[i_density]);
    if(gas_is_mass_density)
        sum /= (mu * m_H);
    return sum;
}

double CGridBasic::getGasNumberDensity(const cell_basic & cell, uint i_density) const
{
    double dens = 0;
    if(size_gd_list > i_density)
        dens = cell.getData(data_pos_gd_list[i_density]);
    if(gas_is_mass_density)
        dens /= (mu * m_H);
    return dens;
}

double CGridBasic::getGasNumberDensity(const photon_package & pp) const
{
    return getGasNumberDensity(*pp.getPositionCell());
}

double CGridBasic::getGasNumberDensity(const photon_package & pp, uint i_density) const
{
    return getGasNumberDensity(*pp.getPositionCell(), i_density);
}

double CGridBasic::getGasMassDensity(const cell_basic & cell) const
{
    double sum = 0;
    for(uint i_density = 0; i_density < size_gd_list; i_density++)
        sum += cell.getData(data_pos_gd_list[i_density]);
    if(!gas_is_mass_density)
        sum *= (mu * m_H);
    return sum;
}

double CGridBasic::getGasMassDensity(const cell_basic & cell, uint i_density) const
{
    double dens = 0;
    if(size_gd_list > i_density)
        dens = cell.getData(data_pos_gd_list[i_density]);
    if(!gas_is_mass_density)
        dens *= (mu * m_H);
    return dens;
}

double CGridBasic::getGasMassDensity(const photon_package & pp) const
{
    return getGasMassDensity(*pp.getPositionCell());
}

double CGridBasic::getGasMassDensity(const photon_package & pp, uint i_density) const
{
    return getGasMassDensity(*pp.getPositionCell(), i_density);
}

bool CGridBasic::useDustChoice()
{
    if(size_gd_list > 1 || size_dd_list > 1)
        return false;
    return true;
}

bool CGridBasic::useConstantGrainSizes()
{
    if(data_pos_amin_list.size() != 0 || data_pos_amax_list.size() != 0 || data_pos_size_param_list.size() != 0)
        return false;
    return true;
}

bool CGridBasic::useDustDensities()
{
    if(size_dd_list > 0)
        return true;
    return false;
}

void CGridBasic::setDustDensity(cell_basic * cell, double val)
{
    if(size_dd_list > 0)
        for(uint i_density = 0; i_density < size_dd_list; i_density++)
            cell->setData(data_pos_dd_list[i_density], val);
    else
        for(uint i_density = 0; i_density < size_gd_list; i_density++)
            cell->setData(data_pos_gd_list[i_density], val);
}

void CGridBasic::setDustDensity(cell_basic * cell, uint i_density, double val)
{
    if(size_dd_list > 0)
        cell->setData(data_pos_dd_list[i_density], val);
    else
        cell->setData(data_pos_gd_list[i_density], val);
}

void CGridBasic::setDustDensity(photon_package * pp, double val)
{
    setDustDensity(pp->getPositionCell(), val);
}

void CGridBasic::setDustDensity(photon_package * pp, uint i_density, double val)
{
    setDustDensity(pp->getPositionCell(), i_density, val);
}

double CGridBasic::getDustDensity(const cell_basic & cell) const
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

double CGridBasic::getDustDensity(const cell_basic & cell, uint i_density) const
{
    if(size_dd_list > i_density)
        return cell.getData(data_pos_dd_list[i_density]);
    else
        return 0;
}

double CGridBasic::getDustDensity(const photon_package & pp) const
{
    return getDustDensity(*pp.getPositionCell());
}

double CGridBasic::getDustDensity(const photon_package & pp, uint i_density) const
{
    return getDustDensity(*pp.getPositionCell(), i_density);
}

double CGridBasic::getRelativeDustDensity(const cell_basic & cell, uint i_density) const
{
    if(getDustDensity(cell) != 0)
        return getDustDensity(cell, i_density) / getDustDensity(cell);
    else if(getGasDensity(cell) != 0)
        return getGasDensity(cell, i_density) / getGasDensity(cell);
    else
        return 0;
}

void CGridBasic::adjustDustDensity(cell_basic * cell, uint i_density, double factor)
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

bool CGridBasic::positionPhotonInGridTest(photon_package * pp)
{
    return false;
}

bool CGridBasic::isVelocityFieldAvailable()
{
    if(data_pos_vx == MAX_UINT || data_pos_vy == MAX_UINT || data_pos_vz == MAX_UINT)
        return false;
    return true;
}

bool CGridBasic::isTurbulentVelocityAvailable()
{
    if(data_pos_vt == MAX_UINT)
        return false;
    return true;
}

Vector3D CGridBasic::getVelocityField(const photon_package & pp) const
{
    if(data_pos_vx == MAX_UINT || data_pos_vy == MAX_UINT || data_pos_vz == MAX_UINT)
        return Vector3D();

    const cell_basic & cell = *pp.getPositionCell();
    Vector3D tmp_dir(cell.getData(data_pos_vx), cell.getData(data_pos_vy), cell.getData(data_pos_vz));
    // Rotate vector from cell center to position
    return rotateToCenter(pp, tmp_dir, true, true);
}

Vector3D CGridBasic::getVelocityField(const cell_basic & cell) const
{
    if(data_pos_vx == MAX_UINT || data_pos_vy == MAX_UINT || data_pos_vz == MAX_UINT)
        return Vector3D();

    return Vector3D(cell.getData(data_pos_vx), cell.getData(data_pos_vy), cell.getData(data_pos_vz));
}

double CGridBasic::getCellAbundance(const photon_package & pp, uint id) const
{
    if(id > nrOfDensRatios - 1)
        return 0;

    return getCellAbundance(*pp.getPositionCell(), id);
}

double CGridBasic::getCellAbundance(const cell_basic & cell, uint id) const
{
    if(id > nrOfDensRatios - 1)
        return 0;

    uint pos = pos_GasSpecRatios[id];
    return cell.getData(pos);
}

double CGridBasic::getOpiateIDParameter(cell_basic * cell, uint id)
{
    if(id > nrOfOpiateIDs - 1)
        return 0;

    uint pos = pos_OpiateIDS[id];

    return cell->getData(pos);
}

Vector3D CGridBasic::getMagField(const cell_basic & cell) const
{
    return Vector3D(cell.getData(data_pos_mx), cell.getData(data_pos_my), cell.getData(data_pos_mz));
}

Vector3D CGridBasic::getMagField(const photon_package & pp) const
{
    const cell_basic & cell = *pp.getPositionCell();
    Vector3D tmp_dir(cell.getData(data_pos_mx), cell.getData(data_pos_my), cell.getData(data_pos_mz));
    // Rotate vector from cell center to position
    return rotateToCenter(pp, tmp_dir, true, true);
}

void CGridBasic::setMagField(cell_basic * cell, const Vector3D & mag)
{
    cell->setData(data_pos_mx, mag.X());
    cell->setData(data_pos_my, mag.Y());
    cell->setData(data_pos_mz, mag.Z());
}

//add some small value to avoid singularities
double CGridBasic::getThetaSync(const photon_package & pp) const
{
    return Vector3D::getAngleTheta(pp.getDirection(), getMagField(pp))+EPS_DOUBLE;
}

double CGridBasic::getThetaMag(const photon_package & pp) const
{
    return Vector3D::getAngleThetaOff(pp.getDirection(), getMagField(pp));
}

double CGridBasic::getPhiMag(const photon_package & pp) const
{
    // 0 deg are in the e_y direction
    return Vector3D::getAnglePhi(pp.getEX(), pp.getEY(), getMagField(pp)) - PI2;
}

double CGridBasic::getTheta(const cell_basic & cell, Vector3D & dir) const
{
    return Vector3D::getAngleTheta(dir, getMagField(cell));
}

double CGridBasic::getThetaPhoton(const photon_package & pp, Vector3D & dir) const
{
    return Vector3D::getAngleTheta(pp.getDirection(), dir);
}

bool CGridBasic::getDustIsMassDensity() const
{
    return dust_is_mass_density;
}

bool CGridBasic::getGasIsMassDensity() const
{
    return gas_is_mass_density;
}

bool CGridBasic::isRadiationFieldAvailable() const
{
    if(data_pos_rx_list.empty() || data_pos_ry_list.empty() || data_pos_rz_list.empty() ||
        data_pos_rf_list.empty())
        return false;
    return true;
}

double CGridBasic::getTotalGasMass() const
{
    return total_gas_mass;
}

void CGridBasic::setDustTemperatureRange(double _min_dust_temp, double _max_dust_temp)
{
    max_dust_temp = _max_dust_temp;
    min_dust_temp = _min_dust_temp;
}

void CGridBasic::setalignedRadiusRange(double a_min, double a_max)
{
    aalg_min = a_min;
    aalg_max = a_max;
}

uint CGridBasic::getTemperatureFieldInformation() const
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

bool CGridBasic::setDataPositionsVariable()
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
                        cout << ERROR_LINE << "Multiple densities and dust choices cannot "
                                "be combined!"
                                << endl;
                        return false;
                    }
                    if(gas_is_mass_density == true)
                    {
                        cout << ERROR_LINE << "Gas number densities cannot be combined "
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
                        cout << ERROR_LINE << "Multiple densities and dust choices cannot "
                                "be combined!"
                                << endl;
                        return false;
                    }
                    if(gas_is_mass_density == false)
                    {
                        cout << ERROR_LINE << "Gas mass densities cannot be combined with "
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
                        cout << ERROR_LINE << "Multiple densities and dust choices cannot "
                                "be combined!"
                                << endl;
                        return false;
                    }
                    if(dust_is_mass_density == true)
                    {
                        cout << ERROR_LINE << "Dust number densities cannot be combined "
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
                        cout << ERROR_LINE << "Multiple densities and dust choices cannot "
                                "be combined!"
                                << endl;
                        return false;
                    }
                    if(dust_is_mass_density == false)
                    {
                        cout << ERROR_LINE << "Dust mass densities cannot be combined "
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
                    cout << ERROR_LINE << "Grid ID " << GRIDgas_temp << " can be set only once!" << endl;
                    return false;
                }

                data_pos_tg = i;
                break;

            case GRIDmx:
                if(data_pos_mx != MAX_UINT)
                {
                    cout << ERROR_LINE << "Grid ID " << GRIDmx << " can be set only once!" << endl;
                    return false;
                }

                data_pos_mx = i;
                break;

            case GRIDmy:
                if(data_pos_my != MAX_UINT)
                {
                    cout << ERROR_LINE << "Grid ID " << GRIDmy << " can be set only once!" << endl;
                    return false;
                }

                data_pos_my = i;
                break;

            case GRIDmz:
                if(data_pos_mz != MAX_UINT)
                {
                    cout << ERROR_LINE << "Grid ID " << GRIDmz << " can be set only once!" << endl;
                    return false;
                }

                data_pos_mz = i;
                break;

            case GRIDvx:
                if(data_pos_vx != MAX_UINT)
                {
                    cout << ERROR_LINE << "Grid ID " << GRIDvx << " can be set only once!" << endl;
                    return false;
                }

                data_pos_vx = i;
                break;

            case GRIDvy:
                if(data_pos_vy != MAX_UINT)
                {
                    cout << ERROR_LINE << "Grid ID " << GRIDvy << " can be set only once!" << endl;
                    return false;
                }

                data_pos_vy = i;
                break;

            case GRIDvz:
                if(data_pos_vz != MAX_UINT)
                {
                    cout << ERROR_LINE << "Grid ID " << GRIDvz << " can be set only once!" << endl;
                    return false;
                }

                data_pos_vz = i;
                break;

            case GRIDpx:
                if(data_pos_px != MAX_UINT)
                {
                    cout << ERROR_LINE << "Grid ID " << GRIDpx << " can be set only once!" << endl;
                    return false;
                }

                data_pos_px = i;
                break;

            case GRIDpy:
                if(data_pos_py != MAX_UINT)
                {
                    cout << ERROR_LINE << "Grid ID " << GRIDpy << " can be set only once!" << endl;
                    return false;
                }

                data_pos_py = i;
                break;

            case GRIDpz:
                if(data_pos_pz != MAX_UINT)
                {
                    cout << ERROR_LINE << "Grid ID " << GRIDpz << " can be set only once!" << endl;
                    return false;
                }

                data_pos_pz = i;
                break;

            case GRIDa_alg:
                data_pos_aalg_list.push_back(i);
                break;

            case GRIDa_min:
                /*if(data_pos_amin != MAX_UINT)
                {
                    cout << ERROR_LINE << "Grid ID " << GRIDa_min << " can be set only once!" << endl;
                    return false;
                }

                data_pos_amin = i;*/
                data_pos_amin_list.push_back(i);
                break;

            case GRIDa_max:
                /*if(data_pos_amax != MAX_UINT)
                {
                    cout << ERROR_LINE << "Grid ID " << GRIDa_max << " can be set only once!" << endl;
                    return false;
                }

                data_pos_amax = i;*/
                data_pos_amax_list.push_back(i);
                break;

            case GRIDq:
                /*if(data_pos_size_param != MAX_UINT)
                {
                    cout << ERROR_LINE << "Grid ID " << GRIDq << " can be set only once!" << endl;
                    return false;
                }

                data_pos_size_param = i;*/
                data_pos_size_param_list.push_back(i);
                break;

            case GRIDv_turb:
                if(data_pos_vt != MAX_UINT)
                {
                    cout << ERROR_LINE << "Grid ID " << GRIDvx << " can be set only once!" << endl;
                    return false;
                }

                data_pos_vt = i;
                break;

            case GRIDn_th:
                if(data_pos_n_th != MAX_UINT)
                {
                    cout << ERROR_LINE << "Grid ID " << GRIDn_th << " can be set only once!" << endl;
                    return false;
                }

                data_pos_n_th = i;
                break;

            case GRIDT_e:
                if(data_pos_T_e != MAX_UINT)
                {
                    cout << ERROR_LINE << "Grid ID " << GRIDT_e << " can be set only once!" << endl;
                    return false;
                }

                data_pos_T_e = i;
                break;

            case GRIDn_cr:
                if(data_pos_n_cr != MAX_UINT)
                {
                    cout << ERROR_LINE << "Grid ID " << GRIDn_cr << " can be set only once!" << endl;
                    return false;
                }

                data_pos_n_cr = i;
                break;

            case GRIDg_min:
                if(data_pos_g_min != MAX_UINT)
                {
                    cout << ERROR_LINE << "Grid ID " << GRIDg_min << " can be set only once!" << endl;
                    return false;
                }

                data_pos_g_min = i;
                break;

            case GRIDg_max:
                if(data_pos_g_max != MAX_UINT)
                {
                    cout << ERROR_LINE << "Grid ID " << GRIDg_max << " can be set only once!" << endl;
                    return false;
                }

                data_pos_g_max = i;
                break;

            case GRIDp:
                if(data_pos_p != MAX_UINT)
                {
                    cout << ERROR_LINE << "Grid ID " << GRIDp << " can be set only once!" << endl;
                    return false;
                }

                data_pos_p = i;
                break;

            case GRIDavg_dir:
                if(data_pos_avg_dir != MAX_UINT)
                {
                    cout << ERROR_LINE << "Grid ID " << GRIDavg_dir << " can be set only once!" << endl;
                    return false;
                }

                data_pos_avg_dir = i;
                break;

            case GRIDavg_th:
                if(data_pos_avg_th != MAX_UINT)
                {
                    cout << ERROR_LINE << "Grid ID " << GRIDavg_th << " can be set only once!" << endl;
                    return false;
                }

                data_pos_avg_th = i;
                break;
                
            case GRIDavg_ux:
                if(data_pos_avg_ux != MAX_UINT)
                {
                    cout << "\nERROR: Grid ID " << GRIDavg_ux << " can be set only once!" << endl;
                    return false;
                }

                data_pos_avg_ux = i;
                break;  

            case GRIDavg_uy:
                if(data_pos_avg_uy != MAX_UINT)
                {
                    cout << "\nERROR: Grid ID " << GRIDavg_uy << " can be set only once!" << endl;
                    return false;
                }

                data_pos_avg_uy = i;
                break; 

            case GRIDavg_uz:
                if(data_pos_avg_uz != MAX_UINT)
                {
                    cout << "\nERROR: Grid ID " << GRIDavg_uz << " can be set only once!" << endl;
                    return false;
                }

                data_pos_avg_uz = i;
                break; 

            case GRID_akRAT:
                data_pos_dust_a_krat_list.push_back(i);
                break;

            case GRID_alarm:
                data_pos_dust_a_larm_list.push_back(i);
                break;                

                
            case GRID_ni:
                if(data_pos_n_i != MAX_UINT)
                {
                    cout << "\nERROR: Grid ID " << GRID_ni << " can be set only once!" << endl;
                    return false;
                }

                data_pos_n_i = i;
                break; 

            case GRID_Z:
                if(data_pos_Z != MAX_UINT)
                {
                    cout << "\nERROR: Grid ID " << GRID_Z << " can be set only once!" << endl;
                    return false;
                }

                data_pos_Z = i;
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
                    cout << ERROR_LINE << "Multiple densities and dust choices cannot be "
                            "combined!"
                            << endl;
                    return false;
                }
                if(data_pos_id != MAX_UINT)
                {
                    cout << ERROR_LINE << "Grid ID " << GRIDdust_id << " can be set only once!" << endl;
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
                cout << ERROR_LINE << "Unknown data IDs in grid file!" << endl;
                cout << "         IDs have to be between " << minGRID << " and " << maxGRID << "!"
                        << endl;
                return false;
        }
    }
    
    size_gd_list = data_pos_gd_list.size();
    size_dd_list = data_pos_dd_list.size();

    if(size_gd_list == 0)
    {
        cout << ERROR_LINE << "Grid requires a gas density! " << endl;
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

bool CGridBasic::createCompatibleTree()
{
    switch(dataID)
    {
        case 0:
            if(data_offset != 6)
            {
                cout << ERROR_LINE << "A grid ID of 0 requires a data length of 6!" << endl;
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
                cout << ERROR_LINE << "A grid ID of 1 requires a data length of 7!" << endl;
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
                cout << ERROR_LINE << "A grid ID of 6 requires a data length of 9!" << endl;
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
                cout << ERROR_LINE << "A grid ID of 7 requires a data length of 10!" << endl;
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
            cout << ERROR_LINE << "Grid ID = " << dataID << " doesn't match any of the old octree IDs!"
                    << endl;
            cout << "       Cannot create a compatible octree!" << endl;
            return false;
    }

    dataID = GRID_ID_OCT;

    return true;
}

uint CGridBasic::CheckSynchrotron(parameters & param)
{
    if(data_pos_n_th == MAX_UINT && data_pos_n_th == MAX_UINT)
    {
        cout << ERROR_LINE << "Neither thermal electrons nor CR electrons are defined in "
                "grid file!                                            \n"
                << endl;
        cout << "       No SYNCHROTRON calculation possible." << endl;
        return MAX_UINT;
    }

    if(data_pos_mx == MAX_UINT)
    {
        cout << ERROR_LINE << "Grid contains no magnetic Bx component!" << endl;
        cout << "       No SYNCHROTRON calculation possible." << endl;
        return MAX_UINT;
    }
    if(data_pos_my == MAX_UINT)
    {
        cout << ERROR_LINE << "Grid contains no magnetic By component!" << endl;
        cout << "       No SYNCHROTRON calculation possible." << endl;
        return MAX_UINT;
    }
    if(data_pos_mz == MAX_UINT)
    {
        cout << ERROR_LINE << "Grid contains no magnetic Bz component!" << endl;
        cout << "       No SYNCHROTRON calculation possible." << endl;
        return MAX_UINT;
    }

    if(data_pos_n_th == MAX_UINT)
    {
        cout << WARNING_LINE << "Grid contains no thermal electron component!" << endl;
        cout << "         Only CR SYNCHROTRON calculation possible." << endl;
    }
    else
    {
        if(data_pos_T_e != MAX_UINT)
        {
            cout << INFO_LINE << "Grid contains a electron temperature component!" << endl;
            cout << "      This component is currently ignored!          " << endl;
        }
    }

    if(data_pos_n_cr == MAX_UINT)
    {
        cout << WARNING_LINE << "Grid contains no thermal electron component!         " << endl;
        cout << "         Only Fraraday RM calculations possible." << endl;
    }
    else
    {
        if(data_pos_g_min == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no gamma_min component!" << endl;
            cout << "       No SYNCHROTRON calculation possible." << endl;
            return MAX_UINT;
        }

        if(data_pos_g_max == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no gamma_max component!" << endl;
            cout << "       No SYNCHROTRON calculation possible." << endl;
            return MAX_UINT;
        }

        if(data_pos_p == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no electron power-law index p component!" << endl;
            cout << "       No SYNCHROTRON calculation possible." << endl;
            return MAX_UINT;
        }
    }

    return 0;
}

uint CGridBasic::CheckOpiate(parameters & param)
{
    if(data_pos_tg == MAX_UINT)
    {
        cout << ERROR_LINE << "Grid contains no gas temperature!" << endl;
        cout << "       No OPIATE calculation possible." << endl;
        return MAX_UINT;
    }
    return 0;
}

uint CGridBasic::CheckTemp(parameters & param, uint & tmp_data_offset)
{
    uint extra_temp_entries = 0;
    if(getTemperatureFieldInformation() == MAX_UINT)
    {
        cout << ERROR_LINE << "The grid does not include the correct information for "
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
            cout << ERROR_LINE << "The grid includes partial/broken information about a "
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
            cout << INFO_LINE << "No gas temperature found in grid." << endl;
            cout << "    Add entry and set gas temperature to dust temperature after "
                    "calculation!"
                    << endl;
            cout << SEP_LINE;
        }
    }
    return 0;
}

uint CGridBasic::CheckRat(parameters & param, uint & tmp_data_offset)
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
    
    if(data_pos_dust_a_krat_list.empty())
    {
        for(uint i_density = 0; i_density < nr_densities; i_density++)
        {
            data_pos_dust_a_krat_list.push_back(data_offset + tmp_data_offset);
            data_ids.push_back(GRID_akRAT);
            tmp_data_offset++;
        }
    }

    if(data_pos_dust_a_larm_list.empty())
    {
        for(uint i_density = 0; i_density < nr_densities; i_density++)
        {
            data_pos_dust_a_larm_list.push_back(data_offset + tmp_data_offset);
            data_ids.push_back(GRID_alarm);
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
    
    if(data_pos_avg_ux == MAX_UINT)
    {
        data_pos_avg_ux = data_offset + tmp_data_offset;
        data_ids.push_back(GRIDavg_ux);
        tmp_data_offset++;
    }

    if(data_pos_avg_uy == MAX_UINT)
    {
        data_pos_avg_uy = data_offset + tmp_data_offset;
        data_ids.push_back(GRIDavg_uy);
        tmp_data_offset++;
    }

    if(data_pos_avg_uz == MAX_UINT)
    {
        data_pos_avg_uz = data_offset + tmp_data_offset;
        data_ids.push_back(GRIDavg_uz);
        tmp_data_offset++;
    }

    if(getTemperatureFieldInformation() == TEMP_EMPTY)
    {
        cout << ERROR_LINE << "Grid contains no dust temperature!" << endl;
        cout << "       No RAT calculation possible." << endl;
        return MAX_UINT;
    }
    else if(getTemperatureFieldInformation() == MAX_UINT)
    {
        cout << ERROR_LINE << "The grid does not include the information for temperature "
                "calculations"
                << endl;
        cout << "       No RAT calculation possible." << endl;
        return MAX_UINT;
    }

    if(data_pos_tg == MAX_UINT)
    {
        cout << ERROR_LINE << "Grid contains no gas temperature!" << endl;
        cout << "       No RAT calculation possible." << endl;
        return MAX_UINT;
    }
    if(data_pos_mx == MAX_UINT)
    {
        cout << WARNING_LINE << "Grid contains no magnetic Bx component!" << endl;
        cout << "         No follow up calculations possible." << endl;
    }
    if(data_pos_my == MAX_UINT)
    {
        cout << WARNING_LINE << "Grid contains no magnetic By component!" << endl;
        cout << "         No follow up calculations possible." << endl;
    }
    if(data_pos_mz == MAX_UINT)
    {
        cout << WARNING_LINE << "Grid contains no magnetic Bz component!" << endl;
        cout << "         No follow up calculations possible." << endl;
    }
    return 0;
}

uint CGridBasic::CheckDustEmission(parameters & param)
{
    // Check if stochastic heating temperatures are saved in grid
    if(data_pos_dt_list.size() > nr_densities && data_pos_dt_list.size() < multi_temperature_entries)
        stochastic_temperature_entries = data_pos_dt_list.size();

    if(getTemperatureFieldInformation() == TEMP_EMPTY)
    {
        cout << ERROR_LINE << "Grid contains no dust temperature!" << endl;
        cout << "       No dust emission possible." << endl;
        return MAX_UINT;
    }
    else if(getTemperatureFieldInformation() == MAX_UINT)
    {
        cout << ERROR_LINE << "The grid does not include the information for temperature "
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
            cout << ERROR_LINE << "The grid includes partial/no information about a "
                    "radiation field!"
                    << endl;
            cout << "       No dust emission with stochastic heating possible." << endl;
            return MAX_UINT;
        }
    }

    if(!data_pos_rf_list.empty() && data_pos_rf_list.size() != WL_STEPS)
    {
        cout << ERROR_LINE << "The grid includes partial/no information about a radiation "
                "field!"
                << endl;
        cout << "       No dust emission possible." << endl;
        return MAX_UINT;
    }

    if(param.getAlign() != 0 && !param.getAligPA())
    {
        if(data_pos_tg == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no gas temperature!" << endl;
            cout << "       No dust emission with aligned dust grains possible." << endl;
            return MAX_UINT;
        }
        if(data_pos_mx == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no magnetic Bx component!" << endl;
            cout << "       No dust emission with aligned dust grains possible." << endl;
            return MAX_UINT;
        }
        if(data_pos_my == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no magnetic By component!" << endl;
            cout << "       No dust emission with aligned dust grains possible." << endl;
            return MAX_UINT;
        }
        if(data_pos_mz == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no magnetic Bz component!" << endl;
            cout << "       No dust emission with aligned dust grains possible." << endl;
            return MAX_UINT;
        }
    }

    if(param.getAligGOLD())
    {
        if(data_pos_vx == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no velocity vx component!" << endl;
            cout << "        No dust emission with GOLD alignment possible." << endl;
            return MAX_UINT;
        }
        if(data_pos_vy == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no velocity vy component!" << endl;
            cout << "        No dust emission with GOLD alignment possible." << endl;
            return MAX_UINT;
        }
        if(data_pos_vz == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no velocity vz component!" << endl;
            cout << "        No dust emission with GOLD alignment possible." << endl;
            return MAX_UINT;
        }
    }

    if(param.getAligRAT())
    {
        if(data_pos_aalg_list.empty())
        {
            cout << ERROR_LINE << "Grid contains no minimum alignment radius for RATs!" << endl;
            cout << "        No dust emission with RAT alignment possible." << endl;
            return MAX_UINT;
        }
        else if(data_pos_aalg_list.size() != 1 && data_pos_aalg_list.size() != nr_densities)
        {
            cout << ERROR_LINE << "Grid contains not the correct amount of minimum alignment radii for "
                    "RATs!"
                    << endl;
            cout << "        No dust emission with RAT alignment possible." << endl;
            return MAX_UINT;
        }
    }
    return 0;
}

uint CGridBasic::CheckDustScattering(parameters & param)
{
    if(param.getAlign() != 0 && !param.getAligPA())
    {
        if(data_pos_tg == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no gas temperature!           " << endl;
            cout << "       No dust scattering calculations with aligned dust grains "
                    "possible."
                    << endl;
            return MAX_UINT;
        }
        if(data_pos_mx == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no magnetic Bx component!     " << endl;
            cout << "       No dust scattering calculations with aligned dust grains "
                    "possible."
                    << endl;
            return MAX_UINT;
        }
        if(data_pos_my == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no magnetic By component!     " << endl;
            cout << "       No dust scattering calculations with aligned dust grains "
                    "possible."
                    << endl;
            return MAX_UINT;
        }
        if(data_pos_mz == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no magnetic Bz component!     " << endl;
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
            cout << ERROR_LINE << "Grid contains no velocity vx component!  " << endl;
            cout << "        No dust scattering calculations with GOLD alignment "
                    "possible."
                    << endl;
            return MAX_UINT;
        }
        if(data_pos_vy == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no velocity vy component!  " << endl;
            cout << "        No dust scattering calculations with GOLD alignment "
                    "possible."
                    << endl;
            return MAX_UINT;
        }
        if(data_pos_vz == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no velocity vz component!  " << endl;
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
            cout << ERROR_LINE << "Grid contains no minimum alignment radius for RATs!" << endl;
            cout << "        No dust scattering calculations with RAT alignment "
                    "possible."
                    << endl;
            return MAX_UINT;
        }
        else if(data_pos_aalg_list.size() != 1 && data_pos_aalg_list.size() != nr_densities)
        {
            cout << ERROR_LINE << "Grid contains not the correct amount of minimum alignment radii for "
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

uint CGridBasic::CheckRadiationForce(parameters & param)
{
    if(getTemperatureFieldInformation() == TEMP_EMPTY)
    {
        cout << ERROR_LINE << "Grid contains no dust temperature!" << endl;
        cout << "       No FORCE calculation possible." << endl;
        return MAX_UINT;
    }
    else if(getTemperatureFieldInformation() == MAX_UINT)
    {
        cout << ERROR_LINE << "The grid does not include the information for temperature "
                "calculations"
                << endl;
        cout << "       No FORCE calculation possible." << endl;
        return MAX_UINT;
    }

    if(!data_pos_rf_list.empty() && data_pos_rf_list.size() != WL_STEPS)
    {
        cout << ERROR_LINE << "The grid includes partial/no information about a radiation "
                "field!"
                << endl;
        cout << "       No FORCE calculation possible." << endl;
        return MAX_UINT;
    }

    return 0;
}

uint CGridBasic::CheckLineEmission(parameters & param)
{
    if(param.getTotalNrOfDustComponents() != 0)
    {
        if(getTemperatureFieldInformation() == TEMP_EMPTY)
        {
            cout << ERROR_LINE << "Grid contains no dust temperature!" << endl;
            cout << "       No line transfer including dust emission possible." << endl;
            return MAX_UINT;
        }
        else if(getTemperatureFieldInformation() == MAX_UINT)
        {
            cout << ERROR_LINE << "The grid does not include the information for "
                    "temperature calculations"
                    << endl;
            cout << "       No line transfer including dust emission possible." << endl;
            return MAX_UINT;
        }
    }

    if(data_pos_tg == MAX_UINT && !param.isGasSpeciesLevelPopMC())
    {
        cout << ERROR_LINE << "Grid contains no gas temperature!" << endl;
        cout << "       No line transfer with possible.  " << endl;
        return MAX_UINT;
    }

    if(velocity_field_needed)
    {
        if(data_pos_mx == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no magnetic Bx component!" << endl;
            cout << "       No line transfer possible." << endl;
            return MAX_UINT;
        }
        if(data_pos_my == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no magnetic By component!" << endl;
            cout << "       No line transfer possible." << endl;
            return MAX_UINT;
        }
        if(data_pos_mz == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no magnetic Bz component!" << endl;
            cout << "       No line transfer possible." << endl;
            return MAX_UINT;
        }
    }

    // Velocity field should be simply zero then.
    // if(param.getKeplerStarMass() == 0)
    // {
    //     if(data_pos_vx == MAX_UINT)
    //     {
    //         cout << ERROR_LINE << "Grid contains no velocity vx component!" << endl;
    //         cout << "        No line transfer possible." << endl;
    //         return MAX_UINT;
    //     }
    //     if(data_pos_vy == MAX_UINT)
    //     {
    //         cout << ERROR_LINE << "Grid contains no velocity vy component!" << endl;
    //         cout << "        No line transfer possible." << endl;
    //         return MAX_UINT;
    //     }
    //     if(data_pos_vz == MAX_UINT)
    //     {
    //         cout << ERROR_LINE << "Grid contains no velocity vz component!" << endl;
    //         cout << "        No line transfer possible." << endl;
    //         return MAX_UINT;
    //     }
    // }
    return 0;
}

uint CGridBasic::CheckProbing(parameters & param)
{
    if(getTemperatureFieldInformation() == TEMP_EMPTY)
    {
        cout << ERROR_LINE << "Grid contains no dust temperature!" << endl;
        cout << "       No LOS analysis with aligned dust grains possible." << endl;
        return MAX_UINT;
    }
    else if(getTemperatureFieldInformation() == MAX_UINT)
    {
        cout << ERROR_LINE << "The grid does not include the information for temperature "
                "calculations"
                << endl;
        cout << "       No LOS analysis with aligned dust grains possible." << endl;
        return MAX_UINT;
    }

    if(!data_pos_rf_list.empty() && data_pos_rf_list.size() != WL_STEPS)
    {
        cout << ERROR_LINE << "The grid includes partial/no information about a radiation "
                "field!"
                << endl;
        cout << "       No LOS analysis with aligned dust grains possible." << endl;
        return MAX_UINT;
    }

    if(param.getAlign() != 0 && !param.getAligPA())
    {
        if(data_pos_tg == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no gas temperature!" << endl;
            cout << "       No LOS analysis with aligned dust grains possible." << endl;
            return MAX_UINT;
        }
        if(data_pos_mx == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no magnetic Bx component!" << endl;
            cout << "       No LOS analysis with aligned dust grains possible." << endl;
            return MAX_UINT;
        }
        if(data_pos_my == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no magnetic By component!" << endl;
            cout << "       No LOS analysis with aligned dust grains possible." << endl;
            return MAX_UINT;
        }
        if(data_pos_mz == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no magnetic Bz component!" << endl;
            cout << "       No LOS analysis with aligned dust grains possible." << endl;
            return MAX_UINT;
        }
    }

    if(param.getAligGOLD())
    {
        if(data_pos_vx == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no velocity vx component!" << endl;
            cout << "        No LOS analysis with GOLD alignment possible." << endl;
            return MAX_UINT;
        }
        if(data_pos_vy == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no velocity vy component!" << endl;
            cout << "        No LOS analysis with GOLD alignment possible." << endl;
            return MAX_UINT;
        }
        if(data_pos_vz == MAX_UINT)
        {
            cout << ERROR_LINE << "Grid contains no velocity vz component!" << endl;
            cout << "        No LOS analysis with GOLD alignment possible." << endl;
            return MAX_UINT;
        }
    }
    return 0;
}
