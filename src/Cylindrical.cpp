#include "typedefs.h"
#include "Cylindrical.h"
#include "MathFunctions.h"
#include "CommandParser.h"
#include <limits>

bool CGridCylindrical::loadGridFromBinrayFile(parameter & param, uint _data_len)
{
    ushort tmpID, tmpOffset;
    string filename = param.getPathGrid();
    bool b_center = false;
    bool additional_cells = false;

    uint r_counter = 0;
    uint ph_counter = 0;
    uint z_counter = 0;

    line_counter = 0;
    char_counter = 0;

    turbulent_velocity = param.getTurbulentVelocity();

    ifstream bin_reader(filename.c_str(), ios::in | ios::binary);

    if(bin_reader.fail())
    {
        cout << "ERROR: Can not load binary cylindrical grid file:" << endl;
        cout << filename << "\n" << endl;
        return false;
    }

    resetGridValues();

    min_len = 1e300;

    max_cells = 0;

    line_counter = 1;
    char_counter = 0;
    uint last_percentage = 0;

    bin_reader.read((char*) &tmpID, 2);
    bin_reader.read((char*) &tmpOffset, 2);

    dataID = tmpID;
    data_offset = (uint) tmpOffset;
    data_len = _data_len + data_offset;

    if(dataID == GRID_ID_CYL)
    {
        data_ids.resize(data_offset);

        for(uint i = 0; i < data_offset; i++)
        {
            ushort tmp_ids = 0;
            bin_reader.read((char*) &tmp_ids, 2);
            data_ids[i] = tmp_ids;
        }

        if(!setDataPositionsVariable())
            return false;
    }
    else
    {
        cout << "ERROR: A cylindrical grid requires an ID of \"" << GRID_ID_CYL << "\"!" << endl;
        return false;
    }

    uint tmp_data_offset = validateDataPositions(param);
    if(tmp_data_offset == MAX_UINT)
        return false;

    bin_reader.read((char*) &Rmin, 8);
    bin_reader.read((char*) &Rmax, 8);
    bin_reader.read((char*) &Zmax, 8);
    bin_reader.read((char*) &N_r, 2);
    bin_reader.read((char*) &N_ph, 2);
    bin_reader.read((char*) &N_z, 2);
    bin_reader.read((char*) &log_factorR, 8);
    bin_reader.read((char*) &log_factorPh, 8);
    bin_reader.read((char*) &log_factorZ, 8);

    grid_cells = new cell_cyl***[N_r];
    center_cells = new cell_cyl*[N_z];

    for(uint i_r = 0; i_r < N_r; i_r++)
    {
        grid_cells[i_r] = new cell_cyl**[N_ph];

        cout << "Allocating memory for cylindrical grid cells : " <<
                float(100.0 * double(i_r) / double(N_r)) << "      \r" << flush;

        for(uint i_ph = 0; i_ph < N_ph; i_ph++)
        {
            grid_cells[i_r][i_ph] = new cell_cyl*[N_z];

            for(uint i_z = 0; i_z < N_z; i_z++)
            {
                grid_cells[i_r][i_ph][i_z] = 0;
                center_cells[i_z] = 0;
            }
        }
    }

    cout << CLR_LINE;
    Rmin *= conv_length_in_SI;
    Rmax *= conv_length_in_SI;
    Zmax *= conv_length_in_SI;

    listR = new double[N_r + 1];
    listPhi = new double[N_ph + 1];
    listZ = new double*[N_r];

    if(log_factorR == 0)
    {
        listR[0] = Rmin;
        listR[N_r] = Rmax;
        for(uint i_r = 1; i_r < N_r; i_r++)
        {
            bin_reader.read((char*) &listR[i_r], 8);
            listR[i_r] *= conv_length_in_SI;
        }
    }
    else if(log_factorR == 1.0)
        CMathFunctions::SinList(Rmin, Rmax, listR, N_r + 1, log_factorR);
    else
        CMathFunctions::ExpList(Rmin, Rmax, listR, N_r + 1, log_factorR);

    if(log_factorPh == 0)
    {
        listPhi[0] = 0;
        listPhi[N_ph] = PIx2;
        for(uint i_ph = 1; i_ph < N_ph; i_ph++)
            bin_reader.read((char*) &listPhi[i_ph], 8);
    }
    else
        CMathFunctions::LinearList(0, PIx2, listPhi, N_ph + 1);

    if(log_factorZ == 0)
    {
        // Custom radial cell list consists of all cell borders,
        // but the inner and outer radial ring!
        for(uint i_r = 0; i_r < N_r; i_r++)
        {
            listZ[i_r] = new double[N_z + 1];
            listZ[i_r][0] = -Zmax;
            listZ[i_r][N_z] = Zmax;
        }
        for(uint i_z = 1; i_z < N_z; i_z++)
        {
            double z;
            bin_reader.read((char*) &z, 8);
            for(uint i_r = 0; i_r < N_r; i_r++)
                listZ[i_r][i_z] = z * conv_length_in_SI;
        }
    }
    else if(log_factorZ == -1.0)
    {
        // Save that additional empty cells are required
        additional_cells = true;

        // Add an empty cell at both ends
        N_z += 2;

        // Variable width of each cylindrical ring in z-direction
        for(uint i_r = 0; i_r < N_r; i_r++)
        {
            listZ[i_r] = new double[N_z + 1];
            
            double dz = 0;
            bin_reader.read((char*) &dz, 8);
            double tmp_zmax = N_z * dz * conv_length_in_SI / 2.;

            listZ[i_r][1] = -tmp_zmax;
            listZ[i_r][N_z - 1] = tmp_zmax;

            for(uint i_z = 2; i_z <= N_z - 2; i_z++)
                listZ[i_r][i_z] = -tmp_zmax + i_z * dz;                
        }
        for(uint i_r = 0; i_r < N_r; i_r++)
        {
            listZ[i_r][0] = -Zmax;
            listZ[i_r][N_z] = Zmax;
        }
    }
    else if(log_factorZ == 1.0)
        for(uint i_r = 0; i_r < N_r; i_r++)
        {
            listZ[i_r] = new double[N_z + 1];
            CMathFunctions::SinList(-Zmax, Zmax, listZ[i_r], N_z + 1, log_factorZ);
        }
    else if(log_factorZ > 1.0)
        for(uint i_r = 0; i_r < N_r; i_r++)
        {
            listZ[i_r] = new double[N_z + 1];
            CMathFunctions::ExpListSym(-Zmax, Zmax, listZ[i_r], N_z + 1, log_factorZ);
        }
    else
        for(uint i_r = 0; i_r < N_r; i_r++)
        {
            listZ[i_r] = new double[N_z + 1];
            CMathFunctions::LinearList(-Zmax, Zmax, listZ[i_r], N_z + 1);
        } 

    if(Rmin <= 0)
    {
        cout << "ERROR: Inner radius (Rmin = " << Rmin << " must be larger than zero!" << endl;
        return false;
    }

    if(Rmax <= 0)
    {
        cout << "ERROR: Outer radius (Rmax = " << Rmax << " must be larger than zero!" << endl;
        return false;
    }

    if(Rmax <= Rmin)
    {
        cout << "ERROR: Outer radius (Rmax = " << Rmax << ") must be larger than inner radius (Rmin = " << Rmin << ")!" << endl;
        return false;
    }

    if(Zmax <= 0)
    {
        cout << "ERROR: Maximum vertical extent must be larger than zero (Zmax = " << Zmax << ")!" << endl;
        return false;
    }

    if(N_r < 3)
    {
        cout << "ERROR: Nr.of cells in R direction has to be larger than 2!" << endl;
        return false;
    }

    max_cells = N_r * N_ph * N_z + N_z;
    line_counter = -1;

    max_len = max(2 * Rmax, 2 * Zmax);
    min_len = listR[1] - listR[0];

    total_volume = PI * Rmax * Rmax * (2 * Zmax);

    while(!bin_reader.eof())
    {
        line_counter++;

        if(r_counter < N_r)
        {
            double dr = listR[r_counter + 1] - listR[r_counter];

            if(dr == 0)
            {
                cout << "ERROR: No step size in r-direction of cylindrical grid!" << endl;
                return false;
            }

            if(dr < min_len)
                min_len = dr;
        }

        if(r_counter == 0)
        {
            if(z_counter < N_z)
            {
                double dz = listZ[0][z_counter + 1] - listZ[0][z_counter];

                if(dz == 0)
                {
                    cout << "ERROR: No step size in z-direction of cylindrical grid!" << endl;
                    return false;
                }

                if(dz < min_len)
                    min_len = dz;
            }

            if(ph_counter < N_ph)
            {
                double dph = listPhi[ph_counter + 1] - listPhi[ph_counter];

                if(dph == 0)
                {
                    cout << "ERROR: No step size in phi-direction of cylindrical grid!" << endl;
                    return false;
                }

                double d = 2 * Rmin * dph / PIx2;

                if(d < min_len)
                    min_len = d;
            }
        }

        // Calculate percentage of total progress per source
        uint percentage = uint(100.0 * double(line_counter) / double(max_cells));

        // Show only new percentage number if it changed
        if((percentage - last_percentage) > PERCENTAGE_STEP)
        {
            char_counter++;
            cout << "-> Loading cylindrical grid file : "
                    << percentage << " [%]      \r" << flush;
            last_percentage = percentage;
        }

        if(z_counter == N_z && b_center == false)
        {
            ph_counter++;
            z_counter = 0;
        }

        if(ph_counter == N_ph && b_center == false)
        {
            r_counter++;
            ph_counter = 0;
        }

        if(r_counter == N_r && b_center == false)
        {
            b_center = true;
            z_counter = 0;
        }

        if(z_counter == N_z && b_center == true)
            break;

        cell_cyl * tmp_cell = 0;

        if(b_center == true)
        {
            center_cells[z_counter] = new cell_cyl;
            center_cells[z_counter]->setRID(MAX_UINT);
            center_cells[z_counter]->setPhID(MAX_UINT);
            center_cells[z_counter]->setZID(z_counter);
            center_cells[z_counter]->resize(data_len + tmp_data_offset);
            tmp_cell = center_cells[z_counter];
            tmp_cell->setID(line_counter);
        }
        else
        {
            grid_cells[r_counter][ph_counter][z_counter] = new cell_cyl;
            grid_cells[r_counter][ph_counter][z_counter]->setRID(r_counter);
            grid_cells[r_counter][ph_counter][z_counter]->setPhID(ph_counter);
            grid_cells[r_counter][ph_counter][z_counter]->setZID(z_counter);
            grid_cells[r_counter][ph_counter][z_counter]->resize(data_len + tmp_data_offset);
            tmp_cell = grid_cells[r_counter][ph_counter][z_counter];
        }

        for(uint i = 0; i < data_offset; i++)
        {
            double tmp_data1 = 0;
            if(additional_cells == false || (z_counter < N_z - 1 && z_counter > 0))
                bin_reader.read((char*) &tmp_data1, 8); 
            tmp_cell->setData(i, tmp_data1);
        }

        updateVelocity(tmp_cell, param);

        if(uint(tmp_cell->getData(data_pos_id)) < 0 ||
            uint(tmp_cell->getData(data_pos_id)) > param.getMaxDustComponentChoice())
        {
            cout << "ERROR: Dust ID in grid exceeds maximum number of dust choices available!" << endl;
            return false;
        }

        tmp_cell->setID(line_counter);

        updateDataRange(tmp_cell);

        double tmp_vol = getVolume(tmp_cell);
        total_gas_mass += getGasMassDensity(tmp_cell) * tmp_vol;
        cell_volume += tmp_vol;
        z_counter++;
    }

    bin_reader.close();

    if(max_cells != uint(line_counter))
    {
        cout << "ERROR: Nr. of read in cells do not match the maximal number of expected cells!" << endl;
        return false;
    }

    data_len += tmp_data_offset;
    data_offset += tmp_data_offset;

    cout << CLR_LINE;
    //cout << "Loading cylindrical grid file        : done" << endl;

    /*for(uint i_ph = 0; i_ph < N_ph; i_ph++)
    {
    cout << i_ph << "\t" << listPhi[i_ph] << "\t" << listPhi[i_ph + 1] << "\t" << grid_cells[4][i_ph][2]->getData(2) << endl;
    }

    photon_package pp;
    pp.setPosition(Vector3D(0.01, 2, 0));

    bool res = positionPhotonInGrid(&pp);
    double t = getGasTemperature(&pp);*/


    return true;
}

bool CGridCylindrical::writeGNUPlotFiles(string path, parameter & param)
{
    nrOfGnuPoints = param.getNrOfGnuPoints();
    nrOfGnuVectors = param.getNrOfGnuVectors();
    maxGridLines = param.getmaxGridLines();

    if(nrOfGnuPoints + nrOfGnuVectors == 0)
        return true;

    if(max_cells == 0)
    {
        cout << "Error: Can't plot spherical grid to gnuplot file to:" << endl;
        cout << path;
        cout << "Not enough tree cells available!" << endl;
        return false;
    }

    plt_gas_dens = (data_pos_gd_list.size() > 0); // 1
    plt_dust_dens = false; //param.getPlot(plIDnd) && (data_pos_dd_list.size() > 0); // 2
    plt_gas_temp =  (data_pos_tg != MAX_UINT); // 3
    plt_dust_temp = (!data_pos_dt_list.empty()); // 4
    plt_rat = (data_pos_aalg != MAX_UINT); // 5
    plt_delta = (data_pos_tg != MAX_UINT) && (data_pos_mx != MAX_UINT); // && (data_pos_td != MAX_UINT); // 6
    plt_larm = (data_pos_tg != MAX_UINT) && (data_pos_mx != MAX_UINT); // && (data_pos_td != MAX_UINT); // 7
    plt_mach = (data_pos_vx != MAX_UINT) && (data_pos_tg != MAX_UINT); // 8

    plt_mag = (data_pos_mx != MAX_UINT); // 0
    plt_vel = (data_pos_vx != MAX_UINT); // 1

    if(nrOfGnuPoints <= 1)
    {
        nrOfGnuPoints = max_cells / 10;

        plt_gas_dens = false;
        plt_dust_dens = false;
        plt_gas_temp = false;
        plt_dust_temp = false;
        plt_rat = false;
        plt_delta = false;
        plt_larm = false;
        plt_mach = false;
    }
    else
        nrOfGnuPoints = max_cells / nrOfGnuPoints;

    if(nrOfGnuVectors <= 1)
    {
        nrOfGnuVectors = max_cells / 10;
        plt_mag = false;
        plt_vel = false;
    }
    else
        nrOfGnuVectors = max_cells / nrOfGnuVectors;

    if(nrOfGnuPoints == 0)
        nrOfGnuPoints = 1;

    if(nrOfGnuVectors == 0)
        nrOfGnuVectors = 1;

    stringstream point_header, vec_header, basic_grid_l0, basic_grid_l1;

    string grid_filename = path + "grid_geometry.plt";
    string dens_gas_filename = path + "grid_gas_density.plt";
    string dens_dust_filename = path + "grid_dust_density.plt";
    string temp_gas_filename = path + "grid_gas_temp.plt";
    string temp_dust_filename = path + "grid_dust_temp.plt";
    string rat_filename = path + "grid_RAT.plt";
    string delta_filename = path + "grid_data.dat";
    string larm_filename = path + "grid_mag.plt";
    string mach_filename = path + "grid_vel.plt";
    string mag_filename = path + "grid_mag.plt";
    string vel_filename = path + "grid_vel.plt";

    ofstream point_fields[9];
    ofstream vec_fields[2];

    point_fields[0].open(grid_filename.c_str());

    if(point_fields[0].fail())
    {
        cout << "ERROR: Can't write to " << grid_filename << endl;
        return false;
    }

    if(plt_gas_dens)
    {
        point_fields[1].open(dens_gas_filename.c_str());

        if(point_fields[1].fail())
        {
            cout << "ERROR: Can't write to " << dens_gas_filename << endl;
            return false;
        }
    }

    if(plt_dust_dens)
    {
        point_fields[2].open(dens_dust_filename.c_str());

        if(point_fields[2].fail())
        {
            cout << "ERROR: Can't write to " << dens_dust_filename << endl;
            return false;
        }
    }

    if(plt_gas_temp)
    {
        point_fields[3].open(temp_gas_filename.c_str());

        if(point_fields[3].fail())
        {
            cout << "ERROR: Can't write to " << temp_gas_filename
                    << endl;
            return false;
        }
    }

    if(plt_dust_temp)
    {
        point_fields[4].open(temp_dust_filename.c_str());

        if(point_fields[4].fail())
        {
            cout << "ERROR: Can't write to " << temp_dust_filename
                    << endl;
            return false;
        }
    }

    if(plt_rat)
    {
        point_fields[5].open(rat_filename.c_str());

        if(point_fields[5].fail())
        {
            cout << "ERROR: Can't write to " << rat_filename << endl;
            return false;
        }
    }

    if(plt_delta)
    {
        point_fields[6].open(delta_filename.c_str());

        if(point_fields[6].fail())
        {
            cout << "ERROR: Can't write to " << delta_filename << endl;
            return false;
        }
    }

    if(plt_larm)
    {
        point_fields[7].open(larm_filename.c_str());

        if(point_fields[7].fail())
        {
            cout << "ERROR: Can't write to " << larm_filename << endl;
            return false;
        }
    }

    if(plt_mach)
    {
        point_fields[8].open(mach_filename.c_str());

        if(point_fields[8].fail())
        {
            cout << "ERROR: Can't write to " << mach_filename << endl;
            return false;
        }
    }

    if(plt_mag)
    {
        vec_fields[0].open(mag_filename.c_str());

        if(vec_fields[0].fail())
        {
            cout << "ERROR: Can't write to " << mag_filename << endl;
            return false;
        }
    }

    if(plt_vel)
    {
        vec_fields[1].open(vel_filename.c_str());

        if(vec_fields[1].fail())
        {
            cout << "ERROR: Can't write to " << vel_filename << endl;
            return false;
        }
    }

    line_counter = 0;
    char_counter = 0;

    point_header.str("");
    point_header << "reset" << endl;
    point_header << "#set terminal postscript" << endl;
    point_header << "#set output \'filename.plt\'" << endl;
    point_header << "set ticslevel 0" << endl;
    point_header << "set size ratio -1" << endl;
    point_header << "set view 45,45" << endl;

    point_header << "set xlabel \'x[m]\'" << endl;
    point_header << "set ylabel \'y[m]\'" << endl;
    point_header << "set zlabel \'z[m]\'" << endl;

    point_header << "set xrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    point_header << "set yrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    point_header << "set zrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;

    point_header << "set style arrow 1 nohead lt 2 lc rgb 0x0000ff lw 0.5" << endl;
    point_header << "set style arrow 2 nohead lt 2 lc rgb 0x0054dd lw 0.5" << endl;

    point_header << "set style line 1 pt 1 ps variable lt palette" << endl;

    point_header << "set grid" << endl;
    point_header << "set nokey" << endl;

    //0 spherical grid
    point_fields[0] << point_header.str();
    point_fields[0] << "set title \'3D spherical grid geometry\' font \'Arial,12\'" << endl;
    point_fields[0] << "set style arrow 3 nohead ls 1 lw 0.5 lc rgb 0x550066" << endl;
    point_fields[0] << "splot '-' with vectors as 3,'-' with vectors as 2,'-' with vectors as 1" << endl;

    //1 gas density
    point_fields[1] << point_header.str();
    point_fields[1] << "set title \'3D gas number density distribution (min: " << min_gas_dens
            << "[m^-3]; max: " << max_gas_dens << "[m^-3])\' font \'Arial,12\'" << endl;
    point_fields[1] << "set cblabel \'gas density[m^-3]\'" << endl;
    point_fields[1] << "set palette defined (0 0.5 0 0, 1 0 0 1, 2 0 1 1, 3 1 1 0)" << endl;

    if(min_gas_dens == 0 && max_gas_dens == 0)
    {
        min_gas_dens = 0.1;
        max_gas_dens = 1;
    }

    if(min_gas_dens == 0)
    {
        min_gas_dens = 0.001 * max_gas_dens;
    }

    if(min_gas_dens / max_gas_dens > 0.9)
        min_gas_dens = 0.9 * max_gas_dens;

    point_fields[1] << "set cbrange[" << log10(min_gas_dens) << ":"
            << log10(max_gas_dens) << "]" << endl;
    point_fields[1] << "set format cb \'%.02g\'" << endl;

    point_fields[1] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    //2 dust density
    point_fields[2] << point_header.str();
    point_fields[2] << "set title \'3D gas number density distribution (min: " << min_dust_dens
            << "[m^-3]; max: " << max_dust_dens << "[m^-3])\' font \'Arial,12\'" << endl;
    point_fields[2] << "set cblabel \'gas density[m^-3]\'" << endl;
    point_fields[2] << "set palette defined (0 0.5 0 0, 1 0 0 1, 2 0 1 1)" << endl;

    if(min_dust_dens == 0 && max_dust_dens == 0)
    {
        min_dust_dens = 0.1;
        max_dust_dens = 1;
    }

    if(min_dust_dens == 0)
    {
        min_dust_dens = 0.001 * max_dust_dens;
    }

    if(min_dust_dens / max_dust_dens > 0.9)
        min_dust_dens = 0.9 * max_dust_dens;

    point_fields[2] << "set cbrange[" << log10(min_dust_dens) << ":"
            << log10(max_dust_dens) << "]" << endl;
    point_fields[2] << "set format cb \'%.02g\'" << endl;

    point_fields[2] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    //3 gas_temp
    point_fields[3] << point_header.str();
    point_fields[3] << "set palette defined (0 0.05 0 0, 0.4 1 0 0, 0.7 1 1 0, 1 1 1 0.5)" << endl;

    point_fields[3] << "set title \'3D gas temperature distribution (min: "
            << min_gas_temp << "[K]; max: " << max_gas_temp
            << "[K])\' font \'Arial,12\'" << endl;
    point_fields[3] << "set cblabel \'temperature [K]\'" << endl;

    if(min_gas_temp == 0 && max_gas_temp == 0)
    {
        min_gas_temp = 0.1;
        max_gas_temp = 1;
    }

    if(min_gas_temp / max_gas_temp > 0.90)
        min_gas_temp = 0.9 * max_gas_temp;

    point_fields[3] << "set cbrange[" << float(min_gas_temp) << ":"
            << float(max_gas_temp) << "]" << endl;
    point_fields[3] << "set format cb \'%.03g\'" << endl;

    point_fields[3] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    //4 dust temp
    point_fields[4] << point_header.str();
    point_fields[4]
            << "set palette defined (0 0.05 0 0, 0.4 1 0 0, 0.7 1 1 0, 1 1 1 0.5)" << endl;

    point_fields[4] << "set title \'3D dust temperature distribution (min: "
            << min_dust_temp << "[K]; max: " << max_dust_temp
            << "[K])\' font \'Arial,12\'" << endl;
    point_fields[4] << "set cblabel \'temperature [K]\'" << endl;

    if(min_dust_temp == 0 && max_dust_temp == 0)
    {
        min_dust_temp = 0.1;
        max_dust_temp = 1;
    }

    if(min_dust_temp / max_dust_temp > 0.9)
        min_dust_temp = 0.9 * max_dust_temp;

    point_fields[4] << "set cbrange[" << float(min_dust_temp) << ":"
            << float(max_dust_temp) << "]" << endl;
    point_fields[4] << "set format cb \'%.03g\'" << endl;

    point_fields[4] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    //5 rat
    point_fields[5] << point_header.str();
    point_fields[5] << "set palette defined (0 0.05 0 0, 0.4 1 0 0, 0.7 1 1 0, 1 1 1 0.5)" << endl;

    point_fields[5] << "set title \'3D aligned radii distribution (min ID: "
            << aalg_min << "; max ID: " << aalg_max
            << ")\' font \'Arial,12\'" << endl;


    point_fields[5] << "set cblabel \'aligned radius ID\'" << endl;

    if(aalg_min == aalg_max)
        aalg_max = 1.01 * aalg_min;

    point_fields[5] << "set cbrange[" << aalg_min << ":" << aalg_max << "]" << endl;

    point_fields[5] << "set format cb \'%.03g\'" << endl;

    point_fields[5] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    vec_header.str("");
    vec_header << "reset" << endl;
    vec_header << "#set terminal postscript" << endl;
    vec_header << "#set output \'\'" << endl;
    vec_header << "set ticslevel 0" << endl;
    vec_header << "set size ratio -1" << endl;
    vec_header << "set view 45,45" << endl;

    vec_header << "set xlabel \'x[m]\'" << endl;
    vec_header << "set ylabel \'y[m]\'" << endl;
    vec_header << "set zlabel \'z[m]\'" << endl;

    vec_header << "set xrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set yrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set zrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;

    vec_header << "set style arrow 1 nohead ls 1 lw 1 lc rgb 0x0000cc" << endl;
    vec_header << "set style arrow 2 nohead ls 1 lw 1 lc rgb 0x5500dd" << endl;
    vec_header << "set style arrow 3 ls 1 lw 1 lc palette" << endl;

    vec_header << "set grid" << endl;
    vec_header << "set nokey" << endl;

    //0 mag
    vec_fields[0] << vec_header.str();
    vec_fields[0] << "set palette defined (0 1 0 0, 0.5 0.0 0.9 0,  0.75 0.0 0.9 1, 0.9 0 0.1 0.9)" << endl;

    if(min_mag == 0 && max_mag == 0)
    {
        min_mag = 1e-45;
        max_mag = 2e-45;
    }

    vec_fields[0]
            << "set title \'3D mag. field distribution (min:" << log10(min_mag) << " log10([T]); max:" << log10(max_mag) << " log10([T])  \' font \'Arial,12\'" << endl;

    if(min_mag / max_mag > 0.9)
        min_mag = 0.9 * max_mag;

    vec_fields[0] << "set cbrange[" << log10(min_mag) << ":" << log10(max_mag) << "]" << endl;
    vec_fields[0] << "set format cb \'%.02g\'" << endl;
    vec_fields[0] << "splot  \'-\' with vectors as 3, \'-\' with vectors as 2, \'-\' with vectors as 1" << endl;

    //1 vel
    vec_fields[1] << vec_header.str();
    vec_fields[1] << "set palette defined (0 1 0 0, 0.5 0.0 0.9 0,  0.75 0.0 0.9 1, 0.9 0 0.1 0.9)" << endl;

    if(min_vel == 0 && max_vel == 0)
    {
        min_vel = 1e-45;
        max_vel = 1e-45;
    }

    vec_fields[1]
            << "set title \'3D vel. field directions (min:" << log10(min_vel) << " log10(m/s); max:" << log10(max_vel) << " log10(m/s)\' font \'Arial,12\'" << endl;

    if(min_vel / max_vel > 0.9)
        min_vel = 0.9 * max_vel;

    vec_fields[1] << "set cbrange[" << float(log10(min_vel)) << ":" << float(log10(max_vel)) << "]" << endl;
    vec_fields[1] << "set format cb \'%.03g\'" << endl;
    vec_fields[1] << "splot  \'-\' with vectors as 3, \'-\' with vectors as 2, \'-\' with vectors as 1" << endl;


    /*for(uint i_r = 0; i_r < N_r; i_r++)
    {
        for(uint i_ph = 0; i_ph < N_ph; i_ph++)
        {
            for(uint i_z = 0; i_z < N_z; i_z++)
            {
                cell_cyl * tmp_cell_pos = grid_cells[i_r][i_ph][i_z];

                Vector3D pos = getCenter(tmp_cell_pos);
                double size = 0.5 + 4 * double(i_r) / double(N_r);

                if(plt_gas_dens)
                {
                    double dens = getGasDensity(tmp_cell_pos);

                    if(dens > 0)
                        point_fields[0] << pos.X() << " " << pos.Y() << " " << pos.Z() << " "
                        << float(size) << " " << log10(dens) << endl;
                }
            }
        }
    }*/

    for(uint pos = 0; pos < 8; pos++)
    {
        //point_fields[pos] << "\ne" << endl;
    }

    for(uint pos = 0; pos < 2; pos++)
    {
        vec_fields[pos] << "\ne" << endl;
    }

    stringstream buffer;
    buffer.str("");

    for(uint i_r = 0; i_r <= N_r; i_r++)
    {
        for(uint i_z = 0; i_z <= N_z; i_z++)
        {
            double r = listR[i_r];
            double z;
            if(i_r == N_r)
                z = listZ[N_r - 1][i_z];
            else
                z = listZ[i_r][i_z];            

            double Nstep = double(int(15 + 25.0 * double(i_r) / double(N_r) + 0.5));
            double tmp_dph = PIx2 / Nstep;

            for(uint i_ph = 0; i_ph < Nstep; i_ph++)
            {
                Vector3D p1(r, tmp_dph*i_ph, z);
                Vector3D p2(r, tmp_dph * (i_ph + 1), z);

                p1.cyl2cart();
                p2.cyl2cart();
                Vector3D dist = p2 - p1;

                buffer << float(p1.X()) << " " << float(p1.Y()) << " "
                        << float(p1.Z()) << " " << float(dist.X()) << " " << float(dist.Y()) << " "
                        << float(dist.Z()) << endl;
            }
        }
    }/**/

    for(uint pos = 0; pos < 8; pos++)
    {
        point_fields[pos] << buffer.str() << "\ne" << endl;
    }

    for(uint pos = 0; pos < 2; pos++)
    {
        vec_fields[pos] << buffer.str() << "\ne" << endl;
    }

    buffer.str("");

    double dPh = PIx2 / double(N_ph);

    for(uint i_r = 0; i_r <= N_r; i_r++)
    {
        for(uint i_ph = 0; i_ph <= N_ph; i_ph++)
        {
            double r = listR[i_r];
            double ph = i_ph*dPh;

            Vector3D p1(r, ph, -Zmax);
            Vector3D p2(r, ph, Zmax);

            p1.cyl2cart();
            p2.cyl2cart();
            Vector3D dist = p2 - p1;

            buffer << float(p1.X()) << " " << float(p1.Y()) << " "
                    << float(p1.Z()) << " " << float(dist.X()) << " " << float(dist.Y()) << " "
                    << float(dist.Z()) << endl;
        }

    }

    for(uint pos = 0; pos < 8; pos++)
    {
        point_fields[pos] << buffer.str() << "\ne" << endl;
    }

    for(uint pos = 0; pos < 2; pos++)
    {
        vec_fields[pos] << "\ne" << endl;
    }

    for(uint i_z = 0; i_z <= N_z; i_z++)
    {
        for(uint i_ph = 0; i_ph <= N_ph; i_ph++)
        {
            double zmin = listZ[0][i_z];
            double zmax = listZ[N_r - 1][i_z];
            double ph = i_ph*dPh;

            Vector3D p1(Rmin, ph, zmin);
            Vector3D p2(Rmax, ph, zmax);

            p1.cyl2cart();
            p2.cyl2cart();
            Vector3D dist = p2 - p1;

            buffer << float(p1.X()) << " " << float(p1.Y()) << " "
                    << float(p1.Z()) << " " << float(dist.X()) << " " << float(dist.Y()) << " "
                    << float(dist.Z()) << endl;
        }

    }

    for(uint pos = 0; pos < 8; pos++)
    {
        point_fields[pos] << buffer.str();
    }


    for(uint pos = 0; pos < 8; pos++)
        point_fields[pos].close();

    for(uint pos = 0; pos < 2; pos++)
        vec_fields[pos].close();

    //cout << "- Writing of gnuplot files      : done" << endl;
    return true;
}

bool CGridCylindrical::saveBinaryGridFile(string filename, ushort id, ushort data_size)
{
    if(max_cells == 0)
    {
        cout << "Error: Cannot save cylindrical grid file to:" << endl;
        cout << filename;
        cout << "Not enough cells available!" << endl;
        return false;
    }

    ofstream bin_writer(filename.c_str(), ios::out | ios::binary);

    if(bin_writer.fail())
    {
        cout << "ERROR: Cannot open cylindrical grid file:" << endl;
        cout << filename;
        return false;
    }

    bin_writer.write((char*) &id, 2);
    bin_writer.write((char*) &data_size, 2);

    if(dataID == GRID_ID_CYL)
    {
        for(uint i = 0; i < data_offset; i++)
        {
            ushort tmp_ids = data_ids[i];
            bin_writer.write((char*) &tmp_ids, 2);
        }
    }
    else
    {
        cout << "Error: Cannot save cylindrical grid file to:" << endl;
        cout << filename;
        cout << "A cylindrical grid requires an ID of " << GRID_ID_CYL << "!" << endl;
        return false;
    }

    bin_writer.write((char*) &Rmin, 8);
    bin_writer.write((char*) &Rmax, 8);
    bin_writer.write((char*) &Zmax, 8);
    bin_writer.write((char*) &N_r, 2);
    bin_writer.write((char*) &N_ph, 2);
    bin_writer.write((char*) &N_z, 2);
    bin_writer.write((char*) &log_factorR, 8);
    bin_writer.write((char*) &log_factorZ, 8);
    if(log_factorR == 0)
        for(uint i_r = 1; i_r < N_r; i_r++)
            bin_writer.write((char*) &listR[i_r], 8);

    for(uint i_r = 0; i_r < N_r; i_r++)
    {
        cout << "-> Writing cylindrical grid file: " <<
                float(100.0 * double(i_r) / double(N_r)) << "      \r" << flush;

        for(uint i_ph = 0; i_ph < N_ph; i_ph++)
        {
            for(uint i_z = 0; i_z < N_z; i_z++)
            {
                for(uint i = 0; i < data_offset; i++)
                {
                    double tmp_data = grid_cells[i_r][i_ph][i_z]->getData(i);
                    bin_writer.write((char*) &tmp_data, 8);
                }
            }
        }
    }

    for(uint i_z = 0; i_z < N_z; i_z++)
    {
        for(uint i = 0; i < data_offset; i++)
        {
            double tmp_data = center_cells[i_z]->getData(i);
            bin_writer.write((char*) &tmp_data, 8);
        }
    }

    bin_writer.close();

    cout << CLR_LINE;
    cout << "- Writing cylindrical grid file : done" << endl;

    return true;
}

/*bool CGridCylindrical::createArtificialGrid(string path)
{
    resetGridValues();

    line_counter = 1;
    char_counter = 0;

    dataID = 6;
    data_offset = 9;
    data_len = 0;
    max_data = 9;
    data_pos_gd_list.push_back(0);
    data_pos_dt_list.push_back(1);
    data_pos_tg = 2;
    data_pos_mx = 3;
    data_pos_my = 4;
    data_pos_mz = 5;
    data_pos_vx = 6;
    data_pos_vy = 7;
    data_pos_vz = 8;

    min_len = 1e30;
    max_len = 2 * Rmax;

    string filename = path;

    dataID = 40;
    Rmin = 1;
    Rmax = 10;
    N_r = 20;
    N_ph = 20;
    N_z = 20;
    log_factorR = 0;
    log_factorZ = 0;

    ofstream bin_writer(filename.c_str(), ios::out | ios::binary);

    if(bin_writer.fail())
    {
        cout << "ERROR: Can not open file:" << endl;
        cout << filename;
        return false;
    }

    bin_writer.write((char*) &dataID, 2);
    bin_writer.write((char*) &data_offset, 2);

    ushort tmp_ids;

    tmp_ids = GRIDgas_dens;
    bin_writer.write((char*) &tmp_ids, 2);

    tmp_ids = GRIDdust_temp;
    bin_writer.write((char*) &tmp_ids, 2);

    tmp_ids = GRIDgas_temp;
    bin_writer.write((char*) &tmp_ids, 2);

    tmp_ids = GRIDmx;
    bin_writer.write((char*) &tmp_ids, 2);

    tmp_ids = GRIDmy;
    bin_writer.write((char*) &tmp_ids, 2);

    tmp_ids = GRIDmz;
    bin_writer.write((char*) &tmp_ids, 2);

    tmp_ids = GRIDvx;
    bin_writer.write((char*) &tmp_ids, 2);

    tmp_ids = GRIDvy;
    bin_writer.write((char*) &tmp_ids, 2);

    tmp_ids = GRIDvz;
    bin_writer.write((char*) &tmp_ids, 2);

    bin_writer.write((char*) &Rmin, 8);
    bin_writer.write((char*) &Rmax, 8);
    bin_writer.write((char*) &N_r, 2);
    bin_writer.write((char*) &N_ph, 2);
    bin_writer.write((char*) &N_z, 2);
    bin_writer.write((char*) &log_factorR, 8);
    bin_writer.write((char*) &log_factorZ, 8);

    for(uint i_r = 0; i_r < N_r; i_r++)
    {
        if(i_r % 50 == 0)
            cout << "-> Creating tree: " << 100.0 * float(i_r) / float(N_r) << " %           \r" << flush;

        for(uint i_ph = 0; i_ph < N_ph; i_ph++)
        {
            for(uint i_z = 0; i_z < N_z; i_z++)
            {
                max_cells++;
                double tmp_data = (i_r + 1)*(i_z + 1)*(i_ph + 1);

                bin_writer.write((char*) &tmp_data, 8);

                tmp_data = 2 * N_ph - i_ph;
                bin_writer.write((char*) &tmp_data, 8);

                tmp_data = 1;
                bin_writer.write((char*) &tmp_data, 8);

                for(uint i = 3; i < 9; i++)
                {
                    tmp_data = double(10);
                    bin_writer.write((char*) &tmp_data, 8);
                }

            }
        }
    }

    for(uint i_z = 0; i_z < N_z; i_z++)
    {
        max_cells++;
        double tmp_data = (i_z + 10);
        bin_writer.write((char*) &tmp_data, 8);

        tmp_data = 1;
        bin_writer.write((char*) &tmp_data, 8);

        tmp_data = 1;
        bin_writer.write((char*) &tmp_data, 8);

        for(uint i = 3; i < 9; i++)
        {
            tmp_data = double(10);
            bin_writer.write((char*) &tmp_data, 8);
        }
    }

    bin_writer.close();

    cout << "min: " << min_gas_dens << "  max_dens: " << max_gas_dens << endl;
    cout << "Creating artificial grid                  : done" << endl;
    cout << "Max cells: " << max_cells << endl;

    return true;
}*/

void CGridCylindrical::printParameter()
{
    if(max_cells == 0)
        cout << "Error: No cylindrical grid parameter available!" << endl;
    else
    {
        cout << CLR_LINE;
        cout << "Cylindrical grid parameter (ID: " << getDataID()
                << ", data len: " << getDataOffset()
                << ", Nr: " << N_r
                << ", Nph: " << N_ph
                << ", Nz: " << N_z
                << ")" << endl;
        cout << SEP_LINE;

        cout << "- Nr. of cylindrical cells      : " << max_cells << endl;

        printPhysicalParameter();
        cout << SEP_LINE;
    }
}

bool CGridCylindrical::positionPhotonInGrid(photon_package * pp)
{
    Vector3D ca_pos = pp->getPosition();

    if(!isInside(ca_pos))
        return false;

    cell_cyl * tmp_cell = 0;
    Vector3D sp_pos = ca_pos.getCylindricalCoord();

    if(ca_pos.X() * ca_pos.X() + ca_pos.Y() * ca_pos.Y() < Rmin * Rmin)
    {
        uint i_Z = CMathFunctions::biListIndexSearch(sp_pos.Z(), listZ[0], N_z + 1);
        cell_cyl * tmp_cell = center_cells[i_Z];
        pp->setPositionCell(tmp_cell);
        return true;
    }

    uint i_R = CMathFunctions::biListIndexSearch(sp_pos.R(), listR, N_r + 1);
    uint i_Phi = CMathFunctions::biListIndexSearch(sp_pos.Phi(), listPhi, N_ph + 1);
    uint i_Z = CMathFunctions::biListIndexSearch(sp_pos.Z(), listZ[i_R], N_z + 1);

    if(i_Z == MAX_UINT)
        return false;

    tmp_cell = grid_cells[i_R][i_Phi][i_Z];

    pp->setPositionCell(tmp_cell);

    return true;
}

bool CGridCylindrical::goToNextCellBorder(photon_package * pp)
{
    Vector3D v_n, v_a;
    cell_cyl * tmp_cell_pos = (cell_cyl*) pp->getPositionCell();
    Vector3D p = pp->getPosition();
    Vector3D d = pp->getDirection();

    bool hit = false;
    double num = 1, den = 1;
    double A = 1, B = 0, C = 0, dscr = 0;
    double path_length = 0, length = 2e300, min_length = 1e300;

    if(tmp_cell_pos->getRID() == MAX_UINT)
    {
        double z1 = listZ[0][tmp_cell_pos->getZID()];
        double z2 = listZ[0][tmp_cell_pos->getZID() + 1];

        for(int i_side = 1; i_side <= 3; i_side++)
        {
            switch(i_side)
            {
                case 1:
                    B = 2 * (p.X() * d.X() + p.Y() * d.Y());
                    C = p.X() * p.X() + p.Y() * p.Y() - Rmin*Rmin;
                    A = d.X() * d.X() + d.Y() * d.Y();
                    dscr = B * B - 4 * A*C;

                    if(dscr >= 0)
                    {
                        dscr = sqrt(dscr);
                        double l1 = (-B + dscr) / (2 * A);
                        double l2 = (-B - dscr) / (2 * A);

                        if(l1 <= 0)
                            l1 = 2e300;

                        if(l2 <= 0)
                            l2 = 2e300;

                        length = min(l1, l2);
                    }

                    break;

                case 2:
                    v_n = Vector3D(0, 0, -1);
                    v_a = Vector3D(0, 0, z1);

                    num = v_n * (p - v_a);
                    den = v_n * d;

                    if(den != 0)
                    {
                        length = -num / den;

                        if(length <= 0)
                            length = 2e300;
                    }
                    else
                        length = 2e300;

                    break;

                case 3:

                    v_n = Vector3D(0, 0, 1);
                    v_a = Vector3D(0, 0, z2);

                    num = v_n * (p - v_a);
                    den = v_n * d;

                    if(den != 0)
                    {
                        length = -num / den;

                        if(length <= 0)
                            length = 2e300;
                    }
                    else
                        length = 2e300;

                    break;
            }

            if(length < min_length)
            {
                min_length = length;
                hit = true;
            }
        }
    }
    else
    {
        double r1 = listR[tmp_cell_pos->getRID()];
        double r2 = listR[tmp_cell_pos->getRID() + 1];
        double ph1 = listPhi[tmp_cell_pos->getPhID()];
        double ph2 = listPhi[tmp_cell_pos->getPhID() + 1];
        double z1 = listZ[tmp_cell_pos->getRID()][tmp_cell_pos->getZID()];
        double z2 = listZ[tmp_cell_pos->getRID()][tmp_cell_pos->getZID() + 1];

        double sin_ph1 = sin(ph1);
        double sin_ph2 = sin(ph2);

        double cos_ph1 = cos(ph1);
        double cos_ph2 = cos(ph2);

        double max_sides = 6;
        if(N_ph == 1)
            max_sides = 4;

        for(int i_side = 1; i_side <= max_sides; i_side++)
        {
            switch(i_side)
            {
                case 1:
                    B = 2 * (p.X() * d.X() + p.Y() * d.Y());
                    C = p.X() * p.X() + p.Y() * p.Y() - r1*r1;
                    A = d.X() * d.X() + d.Y() * d.Y();
                    dscr = B * B - 4 * A*C;

                    if(dscr >= 0)
                    {
                        dscr = sqrt(dscr);
                        double l1 = (-B + dscr) / (2 * A);
                        double l2 = (-B - dscr) / (2 * A);

                        if(l1 <= 0)
                            l1 = 2e300;

                        if(l2 <= 0)
                            l2 = 2e300;

                        length = min(l1, l2);
                    }
                    break;

                case 2:
                    B = 2 * (p.X() * d.X() + p.Y() * d.Y());
                    C = p.X() * p.X() + p.Y() * p.Y() - r2*r2;
                    A = d.X() * d.X() + d.Y() * d.Y();
                    dscr = B * B - 4 * A*C;

                    if(dscr >= 0)
                    {
                        dscr = sqrt(dscr);
                        double l1 = (-B + dscr) / (2 * A);
                        double l2 = (-B - dscr) / (2 * A);

                        if(l1 * l2 < 0)
                        {
                            if(l1 > 0)
                                length = l1;

                            if(l2 > 0)
                                length = l2;
                        }
                        else
                            length = 2e300;
                    }
                    break;

                case 3:
                    v_n = Vector3D(0, 0, -1);
                    v_a = Vector3D(0, 0, z1);

                    num = v_n * (p - v_a);
                    den = v_n * d;

                    if(den != 0)
                    {
                        length = -num / den;

                        if(length <= 0)
                            length = 2e300;
                    }
                    else
                        length = 2e300;
                    break;

                case 4:
                    v_n = Vector3D(0, 0, 1);
                    v_a = Vector3D(0, 0, z2);

                    num = v_n * (p - v_a);
                    den = v_n * d;

                    if(den != 0)
                    {
                        length = -num / den;

                        if(length <= 0)
                            length = 2e300;
                    }
                    else
                        length = 2e300;
                    break;

                case 5:
                    v_n = -Vector3D(-sin_ph1, cos_ph1, 0);
                    v_a = r1 * Vector3D(cos_ph1, sin_ph1, z1);

                    num = v_n * (p - v_a);
                    den = v_n * d;

                    if(den != 0)
                    {
                        length = -num / den;

                        if(length <= 0)
                            length = 2e300;
                    }
                    else
                        length = 2e300;
                    break;

                case 6:
                    v_n = Vector3D(-sin_ph2, cos_ph2, 0);
                    v_a = r1 * Vector3D(cos_ph2, sin_ph2, z1);

                    num = v_n * (p - v_a);
                    den = v_n * d;

                    if(den != 0)
                    {
                        length = -num / den;

                        if(length <= 0)
                            length = 2e300;
                    }
                    else
                        length = 2e300;
                    break;
            }

            if(length < min_length)
            {
                min_length = length;
                hit = true;
            }
        }
    }

    /*if(min_length == 0)
    {
        cout << "ERROR: Min.lenght is zero!" << endl;
        min_length = min_len;
    }

    path_length = 1.0001 * min_length;

    if(min_length < 1e-3 * min_len)
        path_length = 1e-3 * min_len;

    pp->setPosition(p + d * path_length);
    pp->setTmpPathLength(path_length);
    */
    if(hit)
        path_length = min_length + MIN_LEN_STEP * min_len;
    else
    {
        cout << "Error: Wrong cell border!                           " << endl;
        return false;
    }

    pp->setPosition(p + d * path_length);
    pp->setTmpPathLength(path_length);

    return hit;
}

bool CGridCylindrical::updateShortestDistance(photon_package * pp)
{
    Vector3D tmp_pos;
    //    double min_dist, tmp_dist[6];

    //    double loc_x_min, loc_x_max, loc_y_min, loc_y_max, loc_z_min, loc_z_max;
    bool found = false;

    cell_oc * tmp_cell_pos = (cell_oc*) pp->getPositionCell();

    tmp_pos = pp->getPosition();

    /*loc_x_min = tmp_cell_pos->x_min;
    loc_y_min = tmp_cell_pos->y_min;
    loc_z_min = tmp_cell_pos->z_min;

    loc_x_max = tmp_cell_pos->x_max;
    loc_y_max = tmp_cell_pos->y_max;
    loc_z_max = tmp_cell_pos->z_max;

    min_dist = 1E200;

    tmp_dist[0] = abs(loc_x_max - tmp_pos.X());
    tmp_dist[1] = abs(tmp_pos.X() - loc_x_min);

    tmp_dist[2] = abs(loc_y_max - tmp_pos.Y());
    tmp_dist[3] = abs(tmp_pos.Y() - loc_y_min);

    tmp_dist[4] = abs(loc_z_max - tmp_pos.Z());
    tmp_dist[5] = abs(tmp_pos.Z() - loc_z_min);

    for(int i = 0; i < 6; i++)
    {
            if(min_dist > tmp_dist[i])
            {
                    min_dist = tmp_dist[i];
                    found = true;
            }
    }

    pp->setShortestDistance(min_dist);*/
    return found;
}

bool CGridCylindrical::findStartingPoint(photon_package * pp)
{
    uint try_counter = 0;

    double A, B, C, dscr, min_length = 2e300;
    Vector3D v_n, v_a, new_pos;
    Vector3D p = pp->getPosition();
    Vector3D d = pp->getDirection();
    double num, den, length;

    if(isInside(p))
        return true;

    for(int i_side = 1; i_side <= 3; i_side++)
    {
        switch(i_side)
        {
            case 1:
                B = 2 * (p.X() * d.X() + p.Y() * d.Y());
                C = p.X() * p.X() + p.Y() * p.Y() - Rmax*Rmax;
                A = d.X() * d.X() + d.Y() * d.Y();
                dscr = B * B - 4 * A*C;

                if(dscr >= 0 && A > 0)
                {
                    dscr = sqrt(dscr);
                    double l1 = (-B + dscr) / (2 * A);
                    double l2 = (-B - dscr) / (2 * A);

                    if(l1 <= 0)
                        l1 = 2e300;

                    if(l2 <= 0)
                        l2 = 2e300;

                    min_length = min(l1, l2);
                }
                break;

            case 2:
                v_n = Vector3D(0, 0, -1);
                v_a = Vector3D(0, 0, -Zmax);

                num = v_n * (p - v_a);
                den = v_n * d;

                if(den != 0)
                {
                    length = -num / den;

                    if(length > 0)
                        min_length = length;

                }
                break;

            case 3:
                v_n = Vector3D(0, 0, 1);
                v_a = Vector3D(0, 0, Zmax);

                num = v_n * (p - v_a);
                den = v_n * d;

                if(den != 0)
                {
                    length = -num / den;

                    if(length > 0)
                        min_length = length;
                }
                break;
        }

        min_length *= 1.00001;
        new_pos = p + min_length*d;

        if(isInside(new_pos))
        {
            pp->setPosition(new_pos);
            return positionPhotonInGrid(pp);
        }
    }

    return false;
}
