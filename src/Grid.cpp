#include "Grid.h"
#include <CCfits/CCfits>
#include <cmath>

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
    double px = 0;
    double py = 0;
    double pz = 0;
    double a_alg = 0;
    double dust_temp = 0;
    double gas_temp = 0;
    double mach = 0;
    double delta = 0;
    double a_limit = 0;

    if(data_pos_gd_list.size() > 0)
    {
        for(uint i_dens = 0; i_dens < data_pos_gd_list.size(); i_dens++)
            cell->convertData(data_pos_gd_list[i_dens], conv_dens_in_SI);
        gas_dens = getGasDensity(cell);

        if(gas_dens > max_gas_dens)
            max_gas_dens = gas_dens;
        if(gas_dens < min_gas_dens)
            min_gas_dens = gas_dens;

        if(data_pos_dd_list.size() == 0)
        {
            dust_dens = getGasMassDensity(cell) * getMassFraction();

            if(dust_dens > max_dust_dens)
                max_dust_dens = dust_dens;
            if(dust_dens < min_dust_dens)
                min_dust_dens = dust_dens;
        }
    }

    if(data_pos_dd_list.size() > 0)
    {
        for(uint i_dens = 0; i_dens < data_pos_dd_list.size(); i_dens++)
            cell->convertData(data_pos_dd_list[i_dens], conv_dens_in_SI);

        dust_dens = getDustDensity(cell);

        if(dust_dens > max_dust_dens)
            max_dust_dens = dust_dens;
        if(dust_dens < min_dust_dens)
            min_dust_dens = dust_dens;
    }

    if(!data_pos_dt_list.empty())
    {
        dust_temp = cell->getData(data_pos_dt_list[0]);
        //to do if conversion is implemented

        if(dust_temp > max_dust_temp)
            max_dust_temp = dust_temp;
        if(dust_temp < min_dust_temp)
            min_dust_temp = dust_temp;

    }

    if(data_pos_tg != MAX_UINT)
    {
        gas_temp = cell->getData(data_pos_tg);
        //to do if conversion is implemented
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

    if(data_pos_aalg != MAX_UINT)
    {
        a_alg = cell->getData(data_pos_aalg);

        if(a_alg > float(aalg_max))
            aalg_max = (double) a_alg;

        if(a_alg < float(aalg_min))
            aalg_min = (double) a_alg;
    }

    if(data_pos_amin != MAX_UINT)
    {
        double a_min = cell->getData(data_pos_amin);

        if(a_min > float(a_min_max))
            a_min_max = a_min;

        if(a_min < float(a_min_min))
            a_min_min = a_min;
    }

    if(data_pos_amax != MAX_UINT)
    {
        double a_max = cell->getData(data_pos_amax);

        if(a_max > float(a_max_max))
            a_max_max = a_max;

        if(a_max < float(a_max_min))
            a_max_min = a_max;
    }

    if(data_pos_id != MAX_UINT)
    {
        uint dust_id = cell->getData(data_pos_id);

        if(dust_id > float(dust_id_max))
            dust_id_max = dust_id;

        if(dust_id < float(dust_id_min))
            dust_id_min = dust_id;
    }

    //data positions for synchrotron
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
            delta = CMathFunctions::calc_delta(Bfield, dust_temp, gas_temp, gas_dens) * delta0;
            a_limit = CMathFunctions::calc_larm_limit(Bfield, dust_temp, gas_temp, gas_dens, 0.5, larm_f);
        }
    }
    else
    {
        Bfield = 0;
        delta = 0;
        a_limit = 0;
    }

    if(delta > max_delta)
        max_delta = delta;
    if(delta < min_delta)
        min_delta = delta;

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
 #pragma omp critical
         {
             if(cell_count % 500 == 0)
             {
                 cout << "- Filling grid with OPIATE data  : "
                         << 100 * float(cell_count) / float(max_cells) << " [%]               \r";
             }
         }

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
     cout << " - " << found_count << " of " << max_cells << " cells match with the OPIATE paramter file." << endl;
     */return true;
}

uint CGridBasic::validateDataPositions(parameter & param)
{
    uint tmp_data_offset = 0;

    cout << CLR_LINE;

    if(data_pos_gd_list.size() == 0)
    {
        cout << "ERROR: Grid contains no gas (number) density!" << endl;
        cout << "       No RT calculations possible!" << endl;
        return MAX_UINT;
    }

    if(nr_mixtures > 0 && (param.isMonteCarloSimulation() || param.getCommand() == CMD_DUST_EMISSION || 
            param.getCommand() == CMD_LINE_EMISSION || param.getCommand() == CMD_FORCE || 
            param.getCommand() == CMD_PROBING))
    {
        // Get Number of temperature fields for temperature calculation
        if(data_pos_dd_list.size() > 0)
            nr_densities = data_pos_dd_list.size();
        else
            nr_densities = data_pos_gd_list.size();

        // Precalculate the number of temperature entries, if the grid has a 
        // temperature for each grain size or stochastically heated grains
        for(uint i_density = 0; i_density < nr_densities; i_density++)
        {
            multi_temperature_entries += nr_dust_temp_sizes[i_density] + 1;
            stochastic_temperature_entries += nr_stochastic_sizes[i_density] + 1;
        }

        // Check for a valid combination between densities and dust mixtures
        if(nr_densities > 1 && nr_mixtures != nr_densities)
        {
            cout << "ERROR: Amount of densities in the grid (" << nr_densities
                << ") does not fit with the defined dust mixtures (" << nr_mixtures << ")!\n"
                << "(Use a grid with only one density distribution or define more/less dust mixtures!)" << endl;
            return MAX_UINT;
        }
    }

    switch(param.getCommand())
    {
        case CMD_SYNCHROTRON:
            if(SynchrotronCheck(param) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_OPIATE:
            if(OpiateCheck(param) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_TEMP:
            if(TempCheck(param, tmp_data_offset) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_TEMP_RAT:
            if(TempCheck(param, tmp_data_offset) == MAX_UINT)
                return MAX_UINT;

            if(RatCheck(param, tmp_data_offset) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_RAT:
            if(RatCheck(param, tmp_data_offset) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_DUST_EMISSION:
            if(DustEmissionCheck(param) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_DUST_SCATTERING:
            if(DustScatteringCheck(param) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_FORCE:
            if(RadiationForceCheck(param) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_LINE_EMISSION:
            if(LineEmissionCheck(param) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_PROBING:
            if(ProbingCheck(param) == MAX_UINT)
                return MAX_UINT;
            break;

        default:
            cout << "ERROR : Command is unknown!" << endl;
            return MAX_UINT;
    }

    return tmp_data_offset;
}

void CGridBasic::printPhysicalParameter()
{
    /*if(nrOfDensRatios > 0)
    {
        cout << "- Additional dens. ratio IDs    :";
        for(uint i = 0; i < nrOfDensRatios; i++)
            cout << " " << pos_GasSpecRatios[i];

        cout << endl;
    }*/

    cout << "- Volume (total, cells, diff)   : " << total_volume << " [m^3], " << cell_volume
            << " [m^3], " << float(100.0 * abs(total_volume - cell_volume) / max(total_volume, cell_volume))
            << " %" << endl;
    cout << "- Total gas mass                : " << total_gas_mass / M_sun << " [M_sun]" << endl;
    cout << "- Grid length         (min,max) : [" << min_len << ", " << max_len << "] [m]" << endl;
    if(gas_is_mass_density)
        cout << "- Gas mass density    (min,max) : [" << min_gas_dens << ", " << max_gas_dens << "] [kg m^-3]" << endl;
    else
        cout << "- Gas number density  (min,max) : [" << min_gas_dens << ", " << max_gas_dens << "] [m^-3]" << endl;
    if(!dust_is_mass_density && data_pos_dd_list.size() > 0)
        cout << "- Dust number density (min,max) : [" << min_dust_dens << ", " << max_dust_dens << "] [m^-3]" << endl;
    else
        cout << "- Dust mass density   (min,max) : [" << min_dust_dens << ", " << max_dust_dens
            << "] [kg m^-3]" << endl;
    if(data_pos_tg != MAX_UINT)
        cout << "- Gas temperature     (min,max) : [" << min_gas_temp << ", " << max_gas_temp << "] [K]" << endl;
    else
        cout << "- Gas temperature     (min,max) : none" << endl;

    if(!data_pos_dt_list.empty())
        cout << "- Dust temperature    (min,max) : [" << min_dust_temp << ", " << max_dust_temp << "] [K]" << endl;
    else
        cout << "- Dust temperature    (min,max) : none" << endl;

    if(data_pos_mx != MAX_UINT)
    {
        meanBdir.normalize();
        cout << "- Magnetic field      (min,max) : [" << min_mag << ", " << max_mag << "] [T]" << endl;
        cout << "- Mean direction      (norm.)   : X: " << meanBdir.X() << " Y: " << meanBdir.Y() << " Z: " << meanBdir.Z() << endl;
        cout << "- Delta0              (min,max) : [" << min_delta << ", " << max_delta << "] [m]" << endl;
        cout << "- Larm. limit         (min,max) : [" << min_larm_limit << ", " << max_larm_limit << "] [m]" << endl;
    }
    else
        cout << "- Magnetic field     (min,max) : none" << endl;

    if(data_pos_aalg != MAX_UINT)
        cout << "- a_alig              (min,max) : [" << aalg_min << ", " << aalg_max << "] [m]" << endl;

    if(data_pos_vx != MAX_UINT)
    {
        meanVdir.normalize();
        cout << "- Velocity field      (min,max) : [" << min_vel << ", " << max_vel << "] [m/s]" << endl;
        cout << "- Mean direction      (norm.)   : X: " << meanVdir.X() << " Y: " << meanVdir.Y() << " Z: " << meanVdir.Z() << endl;
        cout << "- Mach number         (min,max) : [" << min_mach << ", " << max_mach << "]" << endl;
    }

    if(data_pos_amin != MAX_UINT)
        cout << "- Minimum grain size  (min,max) : [" << a_min_min << ", " << a_min_max << "] [m]" << endl;

    if(data_pos_amax != MAX_UINT)
        cout << "- Maximum grain size  (min,max) : [" << a_max_min << ", " << a_max_max << "] [m]" << endl;

    if(data_pos_id != MAX_UINT)
        cout << "- Dust mixture ID     (min,max) : [" << dust_id_min << ", " << dust_id_max << "]" << endl;

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

    if(data_pos_n_cr != MAX_UINT)
        cout << "- CR el. density      (min,max) : [" << min_n_cr << "; " << max_n_cr << "] [m^-3]" << endl;

    if(data_pos_g_min != MAX_UINT)
        cout << "- Gamma_min           (min,max) : [" << min_g_min << "; " << max_g_min << "]" << endl;

    if(data_pos_g_max != MAX_UINT)
        cout << "- Gamma_max           (min,max) : [" << min_g_max << "; " << max_g_max << "]" << endl;

    if(data_pos_p != MAX_UINT)
        cout << "- El. energy index p  (min,max) : [" << min_p << "; " << max_p << "]" << endl;

    if(data_pos_n_th != MAX_UINT)
        cout << "- Therm. el. density  (min,max) : [" << min_n_th << "; " << max_n_th << "] [m]" << endl;

    if(data_pos_T_e != MAX_UINT)
    {
        if(min_T_e == 1e300)
            cout << "- Electron temperature          : same as dust temperature" << endl;
        else
            cout << "- Electron temp.      (min,max) : [" << min_T_e << "; " << max_T_e << "] [K]" << endl;
    }

    if(data_pos_op != UINT_MAX)
        cout << " - Unique OPIATE IDs" << endl;
};

bool CGridBasic::writeAMIRAFiles(string path, parameter & param, uint bins)
{
    if(bins == 0)
        return true;

    bool plt_gas_dens = (data_pos_gd_list.size() > 0);
    //bool plt_dust_dens = param.getPlot(plIDnd) && (data_pos_dd_list.size() > 0); //to do if dust denity is possible
    bool plt_gas_temp = (data_pos_tg != MAX_UINT);
    bool plt_dust_temp = (!data_pos_dt_list.empty());
    bool plt_mag = (data_pos_mx != MAX_UINT);
    bool plt_vel = (data_pos_vx != MAX_UINT);
    bool plt_rat = (data_pos_aalg != MAX_UINT);
    bool plt_delta = (data_pos_tg != MAX_UINT) && (data_pos_mx != MAX_UINT) && (!data_pos_dt_list.empty());
    bool plt_larm = (data_pos_tg != MAX_UINT) && (data_pos_mx != MAX_UINT) && (!data_pos_dt_list.empty());
    bool plt_mach = (data_pos_vx != MAX_UINT) && (data_pos_tg != MAX_UINT);

    uint per_counter = 0, per_max = bins * bins;
    string dens_filename = path + "gas_density.am";
    string dtemp_filename = path + "dust_temp.am";
    string gtemp_filename = path + "gas_temp.am";
    string magvec_filename = path + "mag_vec_field.am";
    string velvec_filename = path + "vel_vec_field.am";
    string magf_filename = path + "mag_field.am";
    string velf_filename = path + "vel_field.am";
    string a_filename = path + "aalig.am";
    string d_filename = path + "delta.am";

    ofstream dens_writer, rat_writer;
    ofstream gas_writer, dust_writer;
    ofstream magvec_writer, magf_writer;
    ofstream velvec_writer, velf_writer;
    ofstream d_writer;

    if(plt_gas_dens)
    {
        dens_writer.open(dens_filename.c_str(), ios::out);

        if(dens_writer.fail())
        {
            cout << "ERROR: Can't write to " << dens_filename << endl;
            return false;
        }
    }

    if(plt_gas_temp)
    {
        gas_writer.open(gtemp_filename.c_str(), ios::out);

        if(gas_writer.fail())
        {
            cout << "ERROR: Can't write to " << gtemp_filename
                    << endl;
            return false;
        }
    }

    if(plt_dust_temp)
    {
        dust_writer.open(dtemp_filename.c_str(), ios::out);
        if(dust_writer.fail())
        {
            cout << "ERROR: Can't write to " << dtemp_filename
                    << endl;
            return false;
        }
    }

    if(plt_rat)
    {
        rat_writer.open(a_filename.c_str(), ios::out);

        if(rat_writer.fail())
        {
            cout << "ERROR: Can't write to " << a_filename << endl;
            return false;
        }
    }

    if(plt_delta)
    {
        d_writer.open(d_filename.c_str(), ios::out);

        if(d_writer.fail())
        {
            cout << "ERROR: Can't write to " << d_filename << endl;
            return false;
        }
    }

    if(plt_mag)
    {
        magvec_writer.open(magvec_filename.c_str(), ios::out);

        if(magvec_writer.fail())
        {
            cout << "ERROR: Can't write to " << magvec_filename
                    << endl;
            return false;
        }

        magf_writer.open(magf_filename.c_str(), ios::out);

        if(magf_writer.fail())
        {
            cout << "ERROR: Can't write to " << magf_filename << endl;
            return false;
        }
    }

    if(plt_vel)
    {
        velvec_writer.open(velvec_filename.c_str(), ios::out);

        if(velvec_writer.fail())
        {
            cout << "ERROR: Can't write to " << velvec_filename
                    << endl;
            return false;
        }

        velf_writer.open(velf_filename.c_str(), ios::out);

        if(velf_writer.fail())
        {
            cout << "ERROR: Can't write to " << velf_filename << endl;
            return false;
        }
    }

    stringstream point_header, vec_header;
    point_header.str("");
    vec_header.str("");

    int b_limit = int(bins) / 2;
    double xyz_step = max_len / double(bins);

    double off_xyz = 0.5 * xyz_step;
    photon_package * pp = new photon_package;

    point_header << "# AmiraMesh 3D ASCII 2.0\n" << endl;
    point_header << "define Lattice " << bins << " " << bins << " " << bins;
    point_header << "\tParameters {" << endl;
    point_header << "Content \"" << bins << "x" << bins << "x" << bins
            << " float, uniform coordinates\"," << endl;
    point_header << "BoundingBox ";
    point_header << 0 << " "; //float(cell_oc_root->getXmin())
    point_header << 1 << " "; //float(cell_oc_root->getXmax())
    point_header << 0 << " "; //float(cell_oc_root->getXmin())
    point_header << 1 << " "; //float(cell_oc_root->getXmax())
    point_header << 0 << " "; //float(cell_oc_root->getXmin())
    point_header << 1 << " "; //float(cell_oc_root->getXmax())
    point_header << "," << endl;
    point_header << " CoordType \"uniform\"" << endl;
    point_header << "}" << endl;
    point_header << "Lattice { float Data } @1" << endl;
    point_header << "# Data section follows" << endl;
    point_header << "@1" << endl;

    vec_header << "# AmiraMesh 3D ASCII 2.0\n" << endl;
    vec_header << "define Lattice " << bins << " " << bins << " " << bins;
    vec_header << "\tParameters {" << endl;
    vec_header << "Content \"" << bins << "x" << bins << "x" << bins
            << " float[3], uniform coordinates\"," << endl;
    vec_header << "BoundingBox ";
    vec_header << 0 << " "; //float(cell_oc_root->getXmin())
    vec_header << 1 << " "; //float(cell_oc_root->getXmax())
    vec_header << 0 << " "; //float(cell_oc_root->getXmin())
    vec_header << 1 << " "; //float(cell_oc_root->getXmax())
    vec_header << 0 << " "; //float(cell_oc_root->getXmin())
    vec_header << 1 << " "; //float(cell_oc_root->getXmax())
    vec_header << "," << endl;
    vec_header << " CoordType \"uniform\"" << endl;
    vec_header << "}" << endl;
    vec_header << "Lattice { float[3] Data } @1" << endl;
    vec_header << "# Data section follows" << endl;
    vec_header << "@1" << endl;


    dens_writer << "# " << (min_len) << " " << (max_len);
    dens_writer << point_header.str();
    gas_writer << point_header.str();
    dust_writer << point_header.str();
    magf_writer << point_header.str();
    velf_writer << point_header.str();

    magvec_writer << vec_header.str();
    velvec_writer << vec_header.str();

    rat_writer << point_header.str();
    d_writer << point_header.str();

    Vector3D mag_field, vel_field;

    for(int z = -b_limit; z <= b_limit; z++)
    {
        if(z == 0)
            continue;

        for(int y = -b_limit; y <= b_limit; y++)
        {
            if(y == 0)
                continue;

            for(int x = -b_limit; x <= b_limit; x++)
            {
                if(x == 0)
                    continue;

                double sgx = CMathFunctions::sgn(x);
                double sgy = CMathFunctions::sgn(y);
                double sgz = CMathFunctions::sgn(z);
                double tx = double(x) * xyz_step - sgx * off_xyz;
                double ty = double(y) * xyz_step - sgy * off_xyz;
                double tz = double(z) * xyz_step - sgz * off_xyz;

                pp->setPosition(Vector3D(tx, ty, tz));

                if(!positionPhotonInGrid(pp))
                {
                    dens_writer << float(log10(min_gas_dens)) << endl;
                    gas_writer << uint(0) << endl;
                    dust_writer << uint(0) << endl;
                    magf_writer << float(log10(min_mag)) << endl;
                    velf_writer << float(log10(min_vel)) << endl;

                    magvec_writer << uint(0) << " " << uint(0) << " "
                            << uint(0) << endl;
                    velvec_writer << uint(0) << " " << uint(0) << " "
                            << uint(0) << endl;

                    d_writer << uint(0) << endl;

                    rat_writer << uint(0) << endl;
                }
                else
                {
                    dens_writer << float(log10(getGasDensity(pp))) << endl;
                    gas_writer << float((getGasTemperature(pp))) << endl;
                    dust_writer << float((getDustTemperature(pp))) << endl;

                    if(plt_mag)
                    {
                        mag_field = getMagField(pp);
                        magf_writer << float(log10(mag_field.length())) << endl;
                        mag_field.normalize();
                        //mag_field *= max_len;

                        magvec_writer << float(mag_field.X()) << " "
                                << float(mag_field.Y()) << " "
                                << float(mag_field.Z()) << endl;
                    }

                    if(plt_vel)
                    {
                        vel_field = getVelocityField(pp);
                        velf_writer << float(log10(vel_field.length())) << endl;

                        vel_field.normalize();
                        //vel_field *= max_len;

                        velvec_writer << float(vel_field.X()) << " "
                                << float(vel_field.Y()) << " "
                                << float(vel_field.Z()) << endl;
                    }

                    rat_writer << float(log10(getAlignedRadius(pp))) << endl;


                    if(plt_delta)
                    {
                        double field = getMagField(pp).length();
                        double Td = getDustTemperature(pp);
                        double Tg = getGasTemperature(pp);
                        double dens = getGasDensity(pp);
                        double delta = CMathFunctions::calc_delta(field, Td, Tg,
                                dens);
                        d_writer << float(log10(delta)) << endl;
                    }
                }
            }

            per_counter++;
            if(per_counter % 49 == 0)
                cout << " -> Writing AMIRA files:    [ "
                    << 100.0 * float(per_counter) / float(per_max)
                << " % ]                  \r";
        }
    }

    dens_writer.close();
    gas_writer.close();
    dust_writer.close();
    magvec_writer.close();
    magf_writer.close();

    velvec_writer.close();
    velf_writer.close();
    rat_writer.close();
    d_writer.close();

    //cout << "- Writing AMIRA files           : done" << endl;

    delete pp;
    return true;
}

bool CGridBasic::writeSpecialLines(string data_path)
{
    string x_filename = data_path + "center_lines_x.dat";
    string y_filename = data_path + "center_lines_y.dat";
    string z_filename = data_path + "center_lines_z.dat";

    ofstream x_writer, y_writer, z_writer;

    x_writer.open(x_filename.c_str(), ios::out);
    y_writer.open(y_filename.c_str(), ios::out);
    z_writer.open(z_filename.c_str(), ios::out);

    if(x_writer.fail())
    {
        cout << "ERROR: Can't write to " << x_filename << endl;
        return false;
    }

    if(y_writer.fail())
    {
        cout << "ERROR: Can't write to " << y_filename << endl;
        return false;
    }

    if(z_writer.fail())
    {
        cout << "ERROR: Can't write to " << z_filename << endl;
        return false;
    }

    x_writer.precision(8);
    x_writer << scientific;

    y_writer.precision(8);
    y_writer << scientific;

    z_writer.precision(8);
    z_writer << scientific;

    cout << " -> Writing lines : 0.0 %                   \r" << flush;

    photon_package * pp = new photon_package;

    //along z
    pp->setPosition(Vector3D(0, 0, -2.0 * max_len));
    pp->setDirection(Vector3D(0.0001, 0.0001, -1.00001).normalized());
    findStartingPoint(pp);

    z_writer << "ng\tTg\tTd\tmx\tmy\tmz\tvx\tvy\tvz\ta_alg" << endl;

    while(next(pp))
    {
        double pos = pp->getPosition().Z();
        double dens = getGasDensity(pp);
        double Tg = 0;
        double Td = 0;
        double a_alg = 0;
        double mx = 0, my = 0, mz = 0;
        double vx = 0, vy = 0, vz = 0;

        if(data_pos_tg != MAX_UINT)
            Tg = getGasTemperature(pp);

        if(!data_pos_dt_list.empty())
            Td = getDustTemperature(pp);

        if(data_pos_mx != MAX_UINT)
        {
            mx = getMagField(pp).X();
            my = getMagField(pp).Y();
            mz = getMagField(pp).Z();
        }

        if(data_pos_vx != MAX_UINT)
        {
            vx = getVelocityField(pp).X();
            vy = getVelocityField(pp).Y();
            vz = getVelocityField(pp).Z();
        }

        if(data_pos_aalg != MAX_UINT)
            a_alg = getAlignedRadius(pp);

        z_writer << pos << "\t" << dens << "\t" << Tg << "\t" << Td << "\t"
                << mx << "\t" << my << "\t" << mz << "\t"
                << vx << "\t" << vy << "\t" << vz << "\t" << a_alg << endl;
    }

    cout << " -> Writing lines : 33.3 %                   \r" << flush;

    pp->setPosition(Vector3D(0, -2.0 * max_len, 0));
    pp->setDirection(Vector3D(0.0001, -1.00001, 0.0001).normalized());
    findStartingPoint(pp);

    y_writer << "ng\tTg\tTd\tmx\tmy\tmz\tvx\tvy\tvz\ta_alg" << endl;

    while(next(pp))
    {
        double pos = pp->getPosition().Y();
        double dens = getGasDensity(pp);
        double Tg = 0;
        double Td = 0;
        double a_alg = 0;
        double mx = 0, my = 0, mz = 0;
        double vx = 0, vy = 0, vz = 0;

        if(data_pos_tg != MAX_UINT)
            Tg = getGasTemperature(pp);

        if(!data_pos_dt_list.empty())
            Td = getDustTemperature(pp);

        if(data_pos_mx != MAX_UINT)
        {
            mx = getMagField(pp).X();
            my = getMagField(pp).Y();
            mz = getMagField(pp).Z();
        }

        if(data_pos_vx != MAX_UINT)
        {
            vx = getVelocityField(pp).X();
            vy = getVelocityField(pp).Y();
            vz = getVelocityField(pp).Z();
        }

        if(data_pos_aalg != MAX_UINT)
            a_alg = getAlignedRadius(pp);

        y_writer << pos << "\t" << dens << "\t" << Tg << "\t" << Td << "\t"
                << mx << "\t" << my << "\t" << mz << "\t"
                << vx << "\t" << vy << "\t" << vz << "\t" << a_alg << endl;
    }

    cout << " -> Writing lines : 66.6 %                   \r" << flush;

    pp->setPosition(Vector3D(-2.0 * max_len, 0, 0));
    pp->setDirection(Vector3D(-1.00001, 0.0001, 0.0001).normalized());
    findStartingPoint(pp);

    x_writer << "ng\tTg\tTd\tmx\tmy\tmz\tvx\tvy\tvz\ta_alg" << endl;

    while(next(pp))
    {
        double pos = pp->getPosition().X();
        double dens = getGasDensity(pp);
        double Tg = 0;
        double Td = 0;
        double a_alg = 0;
        double mx = 0, my = 0, mz = 0;
        double vx = 0, vy = 0, vz = 0;

        if(data_pos_tg != MAX_UINT)
            Tg = getGasTemperature(pp);

        if(!data_pos_dt_list.empty())
            Td = getDustTemperature(pp);

        if(data_pos_mx != MAX_UINT)
        {
            mx = getMagField(pp).X();
            my = getMagField(pp).Y();
            mz = getMagField(pp).Z();
        }

        if(data_pos_vx != MAX_UINT)
        {
            vx = getVelocityField(pp).X();
            vy = getVelocityField(pp).Y();
            vz = getVelocityField(pp).Z();
        }

        if(data_pos_aalg != MAX_UINT)
            a_alg = getAlignedRadius(pp);

        x_writer << pos << "\t" << dens << "\t" << Tg << "\t" << Td << "\t"
                << mx << "\t" << my << "\t" << mz << "\t"
                << vx << "\t" << vy << "\t" << vz << "\t" << a_alg << endl;
    }

    x_writer.close();
    y_writer.close();
    z_writer.close();

    delete pp;

    cout << "- Writing lines                 : done" << endl;
    return true;
}

bool CGridBasic::writeMidplaneFits(string data_path, parameter & param, uint bins, bool all)
{
    bool res = true;

    if(bins == 0)
        return res;

    int cmd = param.getCommand();

    if(all)
    {
        plt_gas_dens = (!data_pos_gd_list.empty());
        plt_dust_dens = (!data_pos_dd_list.empty());
        plt_gas_temp = (data_pos_tg != MAX_UINT);
        plt_mag = (data_pos_mx != MAX_UINT);
        plt_vel = (data_pos_vx != MAX_UINT);
        plt_delta = (data_pos_tg != MAX_UINT) && (data_pos_mx != MAX_UINT) && (!data_pos_dt_list.empty());
        plt_larm = (data_pos_tg != MAX_UINT) && (data_pos_mx != MAX_UINT) && (!data_pos_dt_list.empty());
        plt_mach = (data_pos_vx != MAX_UINT) && (data_pos_tg != MAX_UINT);
        plt_dust_id = (data_pos_id != MAX_UINT);
        plt_amin = (data_pos_amin != MAX_UINT);
        plt_amax = (data_pos_amax != MAX_UINT);
        plt_n_th = (data_pos_n_th != MAX_UINT);
        plt_T_e = (data_pos_T_e != MAX_UINT);
        plt_n_cr = (data_pos_n_cr != MAX_UINT);
        plt_g_min = (data_pos_g_min != MAX_UINT);
        plt_g_max = (data_pos_g_max != MAX_UINT);
        plt_p = (data_pos_p != MAX_UINT);

        if(cmd != CMD_RAT && cmd != CMD_TEMP_RAT)
            plt_rat = (data_pos_aalg != MAX_UINT);

        if(cmd != CMD_TEMP && cmd != CMD_TEMP_RAT)
            plt_dust_temp = (!data_pos_dt_list.empty());
        
        if(getRadiationFieldAvailable())
        {
            if(param.getWriteRadiationField())
               plt_rad_field = true;
            else if(param.getWriteFullRadiationField())
            {
                plt_rad_field = true;
                nr_rad_field_comp = 4;
            }
            else if(param.getWriteGZero())
                plt_g_zero = true;
        }
    }
    else
    {
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
        plt_n_th = false;
        plt_T_e = false;
        plt_n_cr = false;
        plt_g_min = false;
        plt_g_max = false;
        plt_p = false;

        if(cmd == CMD_TEMP || cmd == CMD_TEMP_RAT)
        {
            if(param.getAdjTgas() > 0)
                plt_gas_temp = true;

            plt_dust_temp = true;
        }

        if(cmd == CMD_RAT || cmd == CMD_TEMP_RAT)
            plt_rat = true;

        if(param.getSaveRadiationField())
        {
            if(param.getWriteRadiationField())
                plt_rad_field = true;
            else if(param.getWriteFullRadiationField())
            {
                plt_rad_field = true;
                nr_rad_field_comp = 4;
            }
            else if(param.getWriteGZero())
                plt_g_zero = true;
        }
        else
            if(param.getWriteRadiationField() || param.getWriteGZero())
                cout << "HINT: The radiation field or G_0 can only be shown in the midplane fits files,\n"
                    << "      if the radiation field will be saved to the grid or is included in the grid!" << endl;
    }

    uint nr_parameter = uint(plt_gas_dens) + uint(plt_dust_dens) + uint(plt_gas_temp) + uint(plt_dust_temp)
            + 4 * uint(plt_mag) + 4 * uint(plt_vel) + uint(plt_rat) + uint(plt_delta) + uint(plt_larm) +
            uint(plt_mach) + uint(plt_dust_id) + uint(plt_rad_field) * nr_rad_field_comp * WL_STEPS + uint(plt_g_zero) + 
            uint(plt_n_th) + uint(plt_T_e) + uint(plt_n_cr) + uint(plt_g_min) + uint(plt_g_max) + uint(plt_p);

    if(nr_parameter == 0)
        return res;

    if(plt_gas_dens)
        if(nr_densities > 1 && data_pos_gd_list.size() >= nr_densities)
            nr_parameter += nr_densities;
    if(plt_dust_dens)
        if(nr_densities > 1 && data_pos_dd_list.size() >= nr_densities)
            nr_parameter += nr_densities;
    if(plt_dust_temp)
        if(nr_densities > 1 && data_pos_dt_list.size() >= nr_densities)
            nr_parameter += nr_densities;

    long naxis = 4;
    long naxes[4] = {uint(bins), uint(bins), 3, nr_parameter};
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

    uint per_counter = 0;

    auto_ptr<CCfits::FITS> pFits(0);
    //unique_ptr<CCfits::FITS> pFits;

    try
    {
        string path_out = data_path + "midplane.fits";
        if(midplane_3d_param.size() == 4)
            path_out = data_path + "midplane_3d.fits";
        remove(path_out.c_str());
        pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
    }
    catch(CCfits::FITS::CantCreate)
    {
        return false;
    }

    long nelements = bins * bins;

    valarray<double> array_gas_dens(nelements);
    valarray<double> array_dust_dens(nelements);
    valarray<double> array_gas_temp(nelements);
    valarray<double> array_dust_temp(nelements);
    valarray<double> array_rat(nelements);
    valarray<double> array_delta(nelements);
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
    valarray<double> array_rad_field(nelements);
    valarray<double> array_g_zero(nelements);
    valarray<double> array_n_th(nelements);
    valarray<double> array_T_e(nelements);
    valarray<double> array_n_cr(nelements);
    valarray<double> array_g_min(nelements);
    valarray<double> array_g_max(nelements);
    valarray<double> array_p(nelements);

    if(plt_gas_dens)
    {
        buffer_gas_dens = new double*[nelements];
        for(int i_cell = 0; i_cell < nelements; i_cell++)
        {
            // +1 for the average/sum of the quantity, but only if multiple quantities are in the grid
            if(nr_densities > 1 && data_pos_gd_list.size() >= nr_densities)
                buffer_gas_dens[i_cell] = new double[nr_densities + 1];
            else
                buffer_gas_dens[i_cell] = new double[nr_densities];
        }
    }
    if(plt_dust_dens)
    {
        buffer_dust_dens = new double*[nelements];
        for(int i_cell = 0; i_cell < nelements; i_cell++)
        {
            // +1 for the average/sum of the quantity, but only if multiple quantities are in the grid
            if(nr_densities > 1 && data_pos_dd_list.size() >= nr_densities)
                buffer_dust_dens[i_cell] = new double[nr_densities + 1];
            else
                buffer_dust_dens[i_cell] = new double[nr_densities];
        }
    }
    if(plt_gas_temp)
        buffer_gas_temp = new double[nelements];
    if(plt_dust_temp)
    {
        buffer_dust_temp = new double*[nelements];
        for(int i_cell = 0; i_cell < nelements; i_cell++)
        {
            // +1 for the average/sum of the quantity, but only if multiple quantities are in the grid
            if(nr_densities > 1 && data_pos_dt_list.size() >= nr_densities)
                buffer_dust_temp[i_cell] = new double[nr_densities + 1];
            else
                buffer_dust_temp[i_cell] = new double[nr_densities];
        }
    }
    if(plt_rat)
        buffer_rat = new double[nelements];
    if(plt_delta)
        buffer_delta = new double[nelements];
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
    if(plt_mach)
        buffer_mach = new double[nelements];
    if(plt_dust_id)
        buffer_dust_mixture = new double[nelements];
    if(plt_amin)
        buffer_dust_amin = new double[nelements];
    if(plt_amax)
        buffer_dust_amax = new double[nelements];
    if(plt_rad_field)
    {
        buffer_rad_field = new double**[nelements];
        for(int i_cell = 0; i_cell < nelements; i_cell++)
        {
            buffer_rad_field[i_cell] = new double*[WL_STEPS];
            for(int wID = 0; wID < WL_STEPS; wID++)
                buffer_rad_field[i_cell][wID] = new double[nr_rad_field_comp];
        }
    }
    if(plt_g_zero)
        buffer_g_zero = new double[nelements];    
    if(plt_n_th)
        buffer_n_th= new double[nelements];  
    if(plt_T_e)
        buffer_T_e= new double[nelements];  
    if(plt_n_cr)
        buffer_n_cr= new double[nelements];  
    if(plt_g_min)
        buffer_g_min= new double[nelements];  
    if(plt_g_max)
        buffer_g_max= new double[nelements];  
    if(plt_p)
        buffer_p= new double[nelements];  

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
            for(int i_cell = 0; i_cell < nelements; i_cell++)
            {
                int j = (i_cell % bins);

                int k = i_cell / bins - b_limit_xy;

                j -= b_limit_xy;


                if(bins % 2 == 0)
                    if(k>-1)
                        k++;

                if(bins % 2 == 0)
                    if(j>-1)
                        j++;

                double tx, ty, tz;
                setPlaneParameter(plane_3d, xy_step, off_xy, z_step, off_z, shift_z,
                    j, k, l, tx, ty, tz);

                fillMidplaneBuffer(tx, ty, tz, i_cell);

                per_counter++;
                if(per_counter % 220 == 0)
                {
#pragma omp critical
                    {
                        cout << " -> Writing 3D midplane file: [ " << 100.0 * float(per_counter) / float(per_max)
                                << " % ]                 \r";
                    }
                }
            }

            fpixel[3] = 0;
            if(plt_gas_dens)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_gas_dens[i_cell] = buffer_gas_dens[i_cell][0];
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_gas_dens);

                if(nr_densities > 1 && data_pos_gd_list.size() >= nr_densities)
                    for(uint i_density = 0; i_density < nr_densities; i_density++)
                    {
                        for(int i_cell = 0; i_cell < nelements; i_cell++)
                            array_gas_dens[i_cell] = buffer_gas_dens[i_cell][i_density + 1];
                        fpixel[3]++;
                        pFits->pHDU().write(fpixel, nelements, array_gas_dens);
                    }
            }
            if(plt_dust_dens)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_dust_dens[i_cell] = buffer_dust_dens[i_cell][0];
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_dust_dens);

                if(nr_densities > 1 && data_pos_dd_list.size() >= nr_densities)
                    for(uint i_density = 0; i_density < nr_densities; i_density++)
                    {
                        for(int i_cell = 0; i_cell < nelements; i_cell++)
                            array_dust_dens[i_cell] = buffer_dust_dens[i_cell][i_density + 1];
                        fpixel[3]++;
                        pFits->pHDU().write(fpixel, nelements, array_dust_dens);
                    }
            }
            if(plt_gas_temp)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_gas_temp[i_cell] = buffer_gas_temp[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_gas_temp);
            }
            if(plt_dust_temp)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_dust_temp[i_cell] = buffer_dust_temp[i_cell][0];
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_dust_temp);

                if(nr_densities > 1 && data_pos_dt_list.size() >= nr_densities)
                    for(uint i_density = 0; i_density < nr_densities; i_density++)
                    {
                        for(int i_cell = 0; i_cell < nelements; i_cell++)
                            array_dust_temp[i_cell] = buffer_dust_temp[i_cell][i_density + 1];
                        fpixel[3]++;
                        pFits->pHDU().write(fpixel, nelements, array_dust_temp);
                    }
            }
            if(plt_rat)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_rat[i_cell] = buffer_rat[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_rat);
            }
            if(plt_delta)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_delta[i_cell] = buffer_delta[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_delta);
            }
            if(plt_mag)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
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
                for(int i_cell = 0; i_cell < nelements; i_cell++)
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
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_larm[i_cell] = buffer_larm[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_larm);
            }
            if(plt_mach)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_mach[i_cell] = buffer_mach[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_mach);
            }
            if(plt_dust_id)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_dust_mixture[i_cell] = buffer_dust_mixture[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_dust_mixture);
            }
            if(plt_amin)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_amin[i_cell] = buffer_dust_amin[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_amin);
            }
                buffer_dust_amin = new double[nelements];
            if(plt_amax)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_amax[i_cell] = buffer_dust_amax[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_amax);
            }
            if(plt_rad_field)
            {
                for (int i_comp = 0; i_comp < nr_rad_field_comp; i_comp++)
                    for(int wID = 0; wID < WL_STEPS; wID++)
                    {
                        for(int i_cell = 0; i_cell < nelements; i_cell++)
                            array_rad_field[i_cell] = buffer_rad_field[i_cell][wID][i_comp];

                        fpixel[3]++;
                        pFits->pHDU().write(fpixel, nelements, array_rad_field);
                    }
            }
            if(plt_g_zero)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_g_zero[i_cell] = buffer_g_zero[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_g_zero);
            }
            if(plt_n_th)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_n_th[i_cell] = buffer_n_th[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_n_th);
            }  
            if(plt_T_e)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_T_e[i_cell] = buffer_T_e[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_T_e);
            }   
            if(plt_n_cr)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_n_cr[i_cell] = buffer_n_cr[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_n_cr);
            }    
            if(plt_g_min)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_g_min[i_cell] = buffer_g_min[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_g_min);
            }     
            if(buffer_g_max)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_g_max[i_cell] = buffer_g_max[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_g_max);
            }    
            if(plt_p)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_p[i_cell] = buffer_p[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_p);
            } 
        }
    }
    else
    {
        for(int i = 1; i <= 3; i++)
        {
            fpixel[2] = i;

#pragma omp parallel for schedule(dynamic)
            for(int i_cell = 0; i_cell < nelements; i_cell++)
            {
                int j = (i_cell % bins);

                int k = i_cell / bins - b_limit_xy;

                j -= b_limit_xy;

                if(bins % 2 == 0)
                    if(k>-1)
                        k++;

                if(bins % 2 == 0)
                    if(j>-1)
                        j++;

                double tx, ty, tz;

                setPlaneParameter(i, xy_step, off_xy, 0, 0, 0, j, k, 0, tx, ty, tz);

                fillMidplaneBuffer(tx, ty, tz, i_cell);

                per_counter++;
                if(per_counter % 220 == 0)
                {
#pragma omp critical
                    {
                        cout << " -> Writing midplane files: [ " << 100.0 * float(per_counter) / float(per_max)
                                << " % ]             \r";
                    }
                }
            }

            fpixel[3] = 0;
            if(plt_gas_dens)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_gas_dens[i_cell] = buffer_gas_dens[i_cell][0];
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_gas_dens);

                if(nr_densities > 1 && data_pos_gd_list.size() >= nr_densities)
                    for(uint i_density = 0; i_density < nr_densities; i_density++)
                    {
                        for(int i_cell = 0; i_cell < nelements; i_cell++)
                            array_gas_dens[i_cell] = buffer_gas_dens[i_cell][i_density + 1];
                        fpixel[3]++;
                        pFits->pHDU().write(fpixel, nelements, array_gas_dens);
                    }
            }
            if(plt_dust_dens)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_dust_dens[i_cell] = buffer_dust_dens[i_cell][0];
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_dust_dens);

                if(nr_densities > 1 && data_pos_dd_list.size() >= nr_densities)
                    for(uint i_density = 0; i_density < nr_densities; i_density++)
                    {
                        for(int i_cell = 0; i_cell < nelements; i_cell++)
                            array_dust_dens[i_cell] = buffer_dust_dens[i_cell][i_density + 1];
                        fpixel[3]++;
                        pFits->pHDU().write(fpixel, nelements, array_dust_dens);
                    }
            }
            if(plt_gas_temp)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_gas_temp[i_cell] = buffer_gas_temp[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_gas_temp);
            }
            if(plt_dust_temp)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_dust_temp[i_cell] = buffer_dust_temp[i_cell][0];
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_dust_temp);

                if(nr_densities > 1 && data_pos_dt_list.size() >= nr_densities)
                    for(uint i_density = 0; i_density < nr_densities; i_density++)
                    {
                        for(int i_cell = 0; i_cell < nelements; i_cell++)
                            array_dust_temp[i_cell] = buffer_dust_temp[i_cell][i_density + 1];
                        fpixel[3]++;
                        pFits->pHDU().write(fpixel, nelements, array_dust_temp);
                    }
            }
            if(plt_rat)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_rat[i_cell] = buffer_rat[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_rat);
            }
            if(plt_delta)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_delta[i_cell] = buffer_delta[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_delta);
            }
            if(plt_mag)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
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
                for(int i_cell = 0; i_cell < nelements; i_cell++)
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
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_larm[i_cell] = buffer_larm[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_larm);
            }
            if(plt_mach)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_mach[i_cell] = buffer_mach[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_mach);
            }
            if(plt_dust_id)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_dust_mixture[i_cell] = buffer_dust_mixture[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_dust_mixture);
            }
            if(plt_amin)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_amin[i_cell] = buffer_dust_amin[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_amin);
            }
                buffer_dust_amin = new double[nelements];
            if(plt_amax)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_amax[i_cell] = buffer_dust_amax[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_amax);
            }
            if(plt_rad_field)
            {
                for (int i_comp = 0; i_comp < nr_rad_field_comp; i_comp++)
                    for(int wID = 0; wID < WL_STEPS; wID++)
                    {
                        for(int i_cell = 0; i_cell < nelements; i_cell++)
                            array_rad_field[i_cell] = buffer_rad_field[i_cell][wID][i_comp];

                        fpixel[3]++;
                        pFits->pHDU().write(fpixel, nelements, array_rad_field);
                    }
            }
            if(plt_g_zero)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_g_zero[i_cell] = buffer_g_zero[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_g_zero);
            }
            if(plt_n_th)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_n_th[i_cell] = buffer_n_th[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_n_th);
            }      
            if(plt_T_e)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_T_e[i_cell] = buffer_T_e[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_T_e);
            }   
            if(plt_n_cr)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_n_cr[i_cell] = buffer_n_cr[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_n_cr);
            }     
            if(plt_g_min)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_g_min[i_cell] = buffer_g_min[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_g_min);
            }     
            if(buffer_g_max)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_g_max[i_cell] = buffer_g_max[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_g_max);
            }      
            if(plt_p)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_p[i_cell] = buffer_p[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_p);
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
        double first_pix_val_z = shift_z - (z_step *  double(naxes[2])) / 2.0 + (bin_width_z / 2.0);

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
        if(nr_densities > 1)
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
    if(plt_dust_dens)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        if(nr_densities > 1)
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
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "rat_aalig [m]", str_2);
    }
    if(plt_delta)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "delta [m]", str_2);
    }
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
    if(plt_mach)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "mach number", str_2);
    }
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
        buffer_dust_amin = new double[nelements];
    if(plt_amax)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "a_max [m]", str_2);
    }
    if(plt_rad_field)
    {
        for (int i_comp = 0; i_comp < nr_rad_field_comp; i_comp++)
        {
            for (int wID = 0; wID < WL_STEPS; wID++)
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
    if(plt_g_zero)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "G_0 (dustem)", str_2);
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

    // Free memory of pointer arrays
    if(plt_gas_dens)
        for(int i_cell = 0; i_cell < nelements; i_cell++)
            delete[] buffer_gas_dens[i_cell];
    if(plt_dust_dens)
        for(int i_cell = 0; i_cell < nelements; i_cell++)
            delete[] buffer_dust_dens[i_cell];
    if(plt_gas_temp)
        delete[] buffer_gas_temp;
    if(plt_dust_temp)
        for(int i_cell = 0; i_cell < nelements; i_cell++)
            delete[] buffer_dust_temp[i_cell];
    if(plt_rat)
        delete[] buffer_rat;
    if(plt_delta)
        delete[] buffer_delta;
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
    if(plt_mach)
        delete[] buffer_mach;
    if(plt_dust_id)
        delete[] buffer_dust_mixture;
    if(plt_amin)
        delete[] buffer_dust_amin;
    if(plt_amax)
        delete[] buffer_dust_amax;
    if(plt_rad_field)
    {
        for(int i_cell = 0; i_cell < nelements; i_cell++)
        {
            for(uint wID = 0; wID < WL_STEPS; wID++)
                delete[] buffer_rad_field[i_cell][wID];
            delete[] buffer_rad_field[i_cell];
        }
        delete[] buffer_rad_field;
    }
    if(plt_g_zero)
        delete[] buffer_g_zero;
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

    cout << CLR_LINE;
    //cout << "- Writing of midplane files     : done" << endl;

    return res;
}
