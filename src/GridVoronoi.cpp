/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include <vector>
#include "GridVoronoi.hpp"
#include "CellOcTree.hpp"

bool CGridVoronoi::loadGridFromBinaryFile(parameters & param, uint _data_len)
{
    ushort tmpID, tmpOffset;
    string filename = param.getPathGrid();

    line_counter = 0;
    char_counter = 0;

    min_nrOfNeigbors = uint(1e6);
    max_nrOfNeigbors = 0;

    ifstream bin_reader(filename.c_str(), ios::in | ios::binary);

    if(bin_reader.fail())
    {
        cout << ERROR_LINE << "Cannot load binary Voronoi grid file: \n";
        cout << filename << "\n\n";
        return false;
    }

    resetGridValues();

    turbulent_velocity = param.getTurbulentVelocity();

    max_cells = 0;

    line_counter = 1;
    char_counter = 0;
    float last_percentage = 0;

    bin_reader.read((char *)&tmpID, 2);
    bin_reader.read((char *)&tmpOffset, 2);

    dataID = tmpID;
    data_offset = (uint)tmpOffset;
    data_len = _data_len + data_offset;

    if(dataID == GRID_ID_VOR)
    {
        data_ids.resize(data_offset);

        for(uint i = 0; i < data_offset; i++)
        {
            ushort tmp_ids = 0;
            bin_reader.read((char *)&tmp_ids, 2);
            data_ids[i] = tmp_ids;
        }

        if(!setDataPositionsVariable())
            return false;

        double tmp_val;
        bin_reader.read((char *)&tmp_val, 8);
        bin_reader.read((char *)&max_len, 8);
        max_cells = ulong(tmp_val);
    }
    else
    {
        cout << ERROR_LINE << "A Voronoi grid requires an ID of \"" << GRID_ID_VOR << "\"!               \n";
        return false;
    }

    if(max_cells < 4)
    {
        cout << ERROR_LINE << "A minimum amount of at least four Voronoi cells is required!   "
                "  \n";
        return false;
    }

    uint tmp_data_offset = validateDataPositions(param);
    if(tmp_data_offset == uint(-1))
        return false;

    max_len *= conv_length_in_SI;
    min_len = max_len / double(max_cells);

    total_volume = max_len * max_len * max_len;

    cell_list = new cell_basic *[max_cells];
    vector<h_list> tmp_h_list;
    stree = new search_tree();
    stree->initTree(max_len);

    cout << CLR_LINE;

    line_counter = 0;

    while(!bin_reader.eof())
    {
        if(line_counter == int(max_cells))
            break;

        // Calculate percentage of total progress per source
        float percentage = 100.0 * double(line_counter) / double(max_cells);

        // Show only new percentage number if it changed
        if((percentage - last_percentage) > PERCENTAGE_STEP)
        {
            char_counter++;
            cout << "-> Loading Voronoi grid file: " << percentage << " [%]      \r" << flush;
            last_percentage = percentage;
        }

        cell_vo * tmp_cell = new cell_vo;
        tmp_cell->resize(data_len + tmp_data_offset);
        tmp_cell->setID(line_counter);

        float tmpCX = 0, tmpCY = 0, tmpCZ = 0;
        double tmp_vol = 0;
        bin_reader.read((char *)&tmpCX, 4);
        bin_reader.read((char *)&tmpCY, 4);
        bin_reader.read((char *)&tmpCZ, 4);
        bin_reader.read((char *)&tmp_vol, 8);

        if(tmp_vol == 0)
        {
            cout << ERROR_LINE << "A Voronoi cell requires a non-zero volume                  "
                    "   \n";
            return false;
        }

        tmpCX *= float(conv_length_in_SI);
        tmpCY *= float(conv_length_in_SI);
        tmpCZ *= float(conv_length_in_SI);
        tmp_vol *= conv_length_in_SI * conv_length_in_SI * conv_length_in_SI;

        tmp_cell->setCenter(tmpCX, tmpCY, tmpCZ);
        tmp_cell->setVolume(tmp_vol);

        for(uint i = 0; i < data_offset; i++)
        {
            float tmp_data = 0;
            bin_reader.read((char *)&tmp_data, 4);
            tmp_cell->setData(i, tmp_data);
        }

        int nr_neighbors = 0;
        int tmp_n = 0;

        bin_reader.read((char *)&nr_neighbors, 4);

        if(nr_neighbors == 0)
        {
            cout << WARNING_LINE << "Voronoi cell nr. " << line_counter + 1
                 << " without neighbors!                     \n";
            cout << "           Check your grid for identical cell positions!            "
                    "         \n";
            cout << "           Possible deadlock in any RT simulation!                  "
                    "   \n";
        }

        if(nr_neighbors < 0)
        {
            nr_neighbors *= -1;
            tmp_h_list.push_back(h_list(tmpCX, tmpCY, tmpCZ, line_counter));
        }

        if(max_nrOfNeigbors < uint(nr_neighbors))
            max_nrOfNeigbors = uint(nr_neighbors);

        if(min_nrOfNeigbors > uint(nr_neighbors))
            min_nrOfNeigbors = uint(nr_neighbors);

        tmp_cell->initNeighbors(ushort(nr_neighbors));

        for(uint i = 0; i < uint(nr_neighbors); i++)
        {
            bin_reader.read((char *)&tmp_n, 4);
            tmp_cell->setNeighbor(i, tmp_n);
        }

        cell_list[line_counter] = tmp_cell;

        if(!stree->addCell(tmp_cell))
        {
            cout << ERROR_LINE << "Failed attempt to add Voronoi cell to the search tree!     "
                    "                \n";
            cout << "       Voronoi cell center nr. " << line_counter + 1
                 << " outside of grid boundaries!                     \n";

            return false;
        }

        updateVelocity(tmp_cell, param);

        if(uint(tmp_cell->getData(data_pos_id)) < 0 ||
           uint(tmp_cell->getData(data_pos_id)) > param.getMaxDustComponentChoice())
        {
            cout << ERROR_LINE << "Dust ID in grid exceeds maximum number of dust choices "
                    "available!   \n";
            return false;
        }

        updateDataRange(tmp_cell);

        cell_volume += tmp_vol;
        total_gas_mass += getGasMassDensity(*tmp_cell) * tmp_vol;

        line_counter++;
    }

    bin_reader.close();

    if(max_cells != uint(line_counter))
    {
        cout << ERROR_LINE << "Amount of read in Voronoi cells (" << uint(line_counter)
             << ")\n does not match the maximal number (" << max_cells << ") of expected cells!  \n";
        return false;
    }

    cout << CLR_LINE;

    hull_size = uint(tmp_h_list.size());
    hull_list = new h_list[hull_size];

    for(uint i = 0; i < hull_size; i++)
    {
        hull_list[i] = tmp_h_list[i];
    }

    data_len += tmp_data_offset;
    data_offset += tmp_data_offset;

    // cout << CLR_LINE;
    // cout << "- Loading Voronoi grid file            : done       \n";

    return true;
}

// plot Voronoi as plot file
bool CGridVoronoi::writePlotFiles(string path, parameters & param)
{
    nrOfPlotPoints = param.getNrOfPlotPoints();
    nrOfPlotVectors = param.getNrOfPlotVectors();
    maxPlotLines = param.getMaxPlotLines();

    if(nrOfPlotPoints + nrOfPlotVectors == 0)
        return true;

    if(max_cells == 0)
    {
        cout << ERROR_LINE << "Cannot plot Voronoi grid to:\n";
        cout << path;
        cout << "Not enough Voronoi cells available! \n";
        return false;
    }

    plt_gas_dens = (size_gd_list > 0);  // 1
    plt_mol_dens = (nrOfDensRatios>0);
    plt_dust_dens = false;                       // param.getPlot(plIDnd) && (!data_pos_dd_list.empty()); // 2
    plt_gas_temp = (data_pos_tg != MAX_UINT);    // 3
    plt_dust_temp = (!data_pos_dt_list.empty()); // 4
    plt_rat = (!data_pos_aalg_list.empty());     // 5
    //plt_delta = false;                           // param.getPlot(plIDdelta) && (data_pos_tg != uint(-1)) &&
                       // (data_pos_mx != uint(-1)) && (data_pos_td != uint(-1)); // 6
    plt_larm = false; // param.getPlot(plIDlarm) && (data_pos_tg != uint(-1)) &&
                      // (data_pos_mx != uint(-1)) && (data_pos_td != uint(-1)); // 7
    plt_mach = false; // param.getPlot(plIDmach) && (data_pos_vx != uint(-1)) &&
                      // (data_pos_tg != uint(-1)); // 8

    plt_mag = (data_pos_mx != uint(-1)); // 0
    plt_vel = (data_pos_vx != uint(-1)); // 1

    if(nrOfPlotPoints <= 1)
    {
        nrOfPlotPoints = max_cells / 10;

        plt_gas_dens = false;
        plt_dust_dens = false;
        plt_gas_temp = false;
        plt_dust_temp = false;
        plt_rat = false;
        //plt_delta = false;
        plt_larm = false;
        plt_mach = false;
    }
    else
        nrOfPlotPoints = max_cells / nrOfPlotPoints;

    if(nrOfPlotVectors <= 1)
    {
        nrOfPlotVectors = max_cells / 10;
        plt_mag = false;
        plt_vel = false;
    }
    else
        nrOfPlotVectors = max_cells / nrOfPlotVectors;

    if(maxPlotLines <= 1)
    {
        maxPlotLines = max_cells / 10;
    }
    else
        maxPlotLines = max_cells / maxPlotLines;

    if(nrOfPlotPoints == 0)
        nrOfPlotPoints = 1;

    if(nrOfPlotVectors == 0)
        nrOfPlotVectors = 1;

    stringstream point_header, vec_header, basic_grid, tetra_lines;

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
        cout << ERROR_LINE << "Cannot write to:\n " << grid_filename << "          \n";
        return false;
    }

    if(plt_gas_dens)
    {
        point_fields[1].open(dens_gas_filename.c_str());

        if(point_fields[1].fail())
        {
            cout << ERROR_LINE << "Cannot write to:\n " << dens_gas_filename << "          \n";
            return false;
        }
    }

    if(plt_dust_dens)
    {
        point_fields[2].open(dens_dust_filename.c_str());

        if(point_fields[2].fail())
        {
            cout << ERROR_LINE << "Cannot write to:\n " << dens_dust_filename << "          \n";
            return false;
        }
    }

    if(plt_gas_temp)
    {
        point_fields[3].open(temp_gas_filename.c_str());

        if(point_fields[3].fail())
        {
            cout << ERROR_LINE << "Cannot write to:\n " << temp_gas_filename << "          \n";
            return false;
        }
    }

    if(plt_dust_temp)
    {
        point_fields[4].open(temp_dust_filename.c_str());

        if(point_fields[4].fail())
        {
            cout << ERROR_LINE << "Cannot write to:\n " << temp_dust_filename << "          \n";
            return false;
        }
    }

    if(plt_rat)
    {
        point_fields[5].open(rat_filename.c_str());

        if(point_fields[5].fail())
        {
            cout << ERROR_LINE << "Cannot write to:\n " << rat_filename << "          \n";
            return false;
        }
    }

    /*if(plt_delta)
    {
        point_fields[6].open(delta_filename.c_str());

        if(point_fields[6].fail())
        {
            cout << ERROR_LINE << "Cannot write to:\n " << delta_filename << "          \n";
            return false;
        }
    }*/

    if(plt_larm)
    {
        point_fields[7].open(larm_filename.c_str());

        if(point_fields[7].fail())
        {
            cout << ERROR_LINE << "Cannot write to:\n " << larm_filename << "          \n";
            return false;
        }
    }

    if(plt_mach)
    {
        point_fields[8].open(mach_filename.c_str());

        if(point_fields[8].fail())
        {
            cout << ERROR_LINE << "Cannot write to:\n " << mach_filename << "          \n";
            return false;
        }
    }

    if(plt_mag)
    {
        vec_fields[0].open(mag_filename.c_str());

        if(vec_fields[0].fail())
        {
            cout << ERROR_LINE << "Cannot write to:\n " << mag_filename << "          \n";
            return false;
        }
    }

    if(plt_vel)
    {
        vec_fields[1].open(vel_filename.c_str());

        if(vec_fields[1].fail())
        {
            cout << ERROR_LINE << "Cannot write to:\n " << vel_filename << "          \n";
            return false;
        }
    }

    cout << "-> Creating plot files  ...        \r" << flush;

    line_counter = 0;
    char_counter = 0;

    // Grid boundaries
    basic_grid.str("");
    double xm = -0.5 * max_len;

    basic_grid.str("");

    basic_grid << xm << " " << xm << " " << xm << " " << max_len << " " << 0 << " " << 0 << " \n";
    basic_grid << xm << " " << xm << " " << xm << " " << 0 << " " << max_len << " " << 0 << " \n";
    basic_grid << xm << " " << xm << " " << xm << " " << 0 << " " << 0 << " " << max_len << " \n";

    basic_grid << xm + max_len << " " << xm + max_len << " " << xm + max_len << " " << -max_len << " " << 0
               << " " << 0 << " \n";
    basic_grid << xm + max_len << " " << xm + max_len << " " << xm + max_len << " " << 0 << " " << -max_len
               << " " << 0 << " \n";
    basic_grid << xm + max_len << " " << xm + max_len << " " << xm + max_len << " " << 0 << " " << 0 << " "
               << -max_len << " \n";

    basic_grid << xm + max_len << " " << xm + max_len << " " << xm << " " << -max_len << " " << 0 << " " << 0
               << " \n";
    basic_grid << xm + max_len << " " << xm + max_len << " " << xm << " " << 0 << " " << -max_len << " " << 0
               << " \n";

    basic_grid << xm + max_len << " " << xm + max_len << " " << xm << " " << -max_len << " " << 0 << " " << 0
               << " \n";
    basic_grid << xm + max_len << " " << xm + max_len << " " << xm << " " << 0 << " " << -max_len << " " << 0
               << " \n";

    basic_grid << xm << " " << xm + max_len << " " << xm + max_len << " " << 0 << " " << -max_len << " " << 0
               << " \n";
    basic_grid << xm << " " << xm + max_len << " " << xm + max_len << " " << 0 << " " << 0 << " " << -max_len
               << " \n";

    basic_grid << xm + max_len << " " << xm << " " << xm + max_len << " " << -max_len << " " << 0 << " " << 0
               << " \n";
    basic_grid << xm + max_len << " " << xm << " " << xm + max_len << " " << 0 << " " << 0 << " " << -max_len
               << " \n";

    cout << CLR_LINE;

    point_header.str("");
    point_header << "reset\n";
    point_header << "#set terminal postscript\n";
    point_header << "#set output \'filename.plt\'\n";
    point_header << "set ticslevel 0\n";
    point_header << "set size ratio -1\n";
    point_header << "set view 45,45\n";

    point_header << "set xlabel \'x[m]\'\n";
    point_header << "set ylabel \'y[m]\'\n";
    point_header << "set zlabel \'z[m]\'\n";

    point_header << "set xrange[" << -0.51 * max_len << ":" << 0.51 * max_len << "]\n";
    point_header << "set yrange[" << -0.51 * max_len << ":" << 0.51 * max_len << "]\n";
    point_header << "set zrange[" << -0.51 * max_len << ":" << 0.51 * max_len << "]\n";

    point_header << "set style arrow 1 nohead lt 2 lc rgb 0x0000ff lw 0.5\n";
    point_header << "set style arrow 2 nohead lt 2 lc rgb 0x0054dd lw 0.5\n";

    point_header << "set style line 1 pt 1 ps variable lt palette\n";

    point_header << "set grid\n";
    point_header << "set nokey\n";

    // 0 voronoi grid
    point_fields[0] << point_header.str();
    point_fields[0] << "set style line 2 pt 1 ps 2 lc rgb 0x550066\n";
    point_fields[0] << "set style arrow 3 nohead lt 2 lc rgb 0xff54dd lw 0.5\n";
    point_fields[0] << "set style line 2 pt 7 ps 0.25 lc rgb 0xff0000\n";
    point_fields[0] << "set title \'3D voronoi grid geometry\' font \'Arial,12\'\n";
    point_fields[0] << "splot '-' with vectors as 3,'-' with vectors as 2,'-' with "
                       "vectors as 1, '-' w p ls 2\n";

    // 1 gas density
    point_fields[1] << point_header.str();
    point_fields[1] << "set title \'3D gas number density distribution (min: " << min_gas_dens
                    << "[m^-3]; max: " << max_gas_dens << "[m^-3])\' font \'Arial,12\'\n";
    point_fields[1] << "set cblabel \'gas density[m^-3]\'\n";
    point_fields[1] << "set palette defined (0 0.5 0 0, 1 0 0 1, 2 0 1 1, 3 1 1 0)\n";

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

    point_fields[1] << "set cbrange[" << log10(min_gas_dens) << ":" << log10(max_gas_dens) << "]\n";
    point_fields[1] << "set format cb \'%.02g\'\n";

    point_fields[1] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1\n";

    // 2 dust density
    point_fields[2] << point_header.str();
    point_fields[2] << "set title \'3D gas number density distribution (min: " << min_dust_dens
                    << "[m^-3]; max: " << max_dust_dens << "[m^-3])\' font \'Arial,12\'\n";
    point_fields[2] << "set cblabel \'gas density[m^-3]\'\n";
    point_fields[2] << "set palette defined (0 0.5 0 0, 1 0 0 1, 2 0 1 1)\n";

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

    point_fields[2] << "set cbrange[" << log10(min_dust_dens) << ":" << log10(max_dust_dens) << "]\n";
    point_fields[2] << "set format cb \'%.02g\'\n";

    point_fields[2] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1\n";

    // 3 gas_temp
    point_fields[3] << point_header.str();
    point_fields[3] << "set palette defined (0 0.05 0 0, 0.4 1 0 0, 0.7 1 1 0, 1 1 1 0.5)\n";

    point_fields[3] << "set title \'3D gas temperature distribution (min: " << min_gas_temp
                    << "[K]; max: " << max_gas_temp << "[K])\' font \'Arial,12\'\n";
    point_fields[3] << "set cblabel \'temperature [K]\'\n";

    if(min_gas_temp == 0 && max_gas_temp == 0)
    {
        min_gas_temp = 0.1;
        max_gas_temp = 1;
    }

    if(min_gas_temp / max_gas_temp > 0.90)
        min_gas_temp = 0.9 * max_gas_temp;

    point_fields[3] << "set cbrange[" << float(min_gas_temp) << ":" << float(max_gas_temp) << "]\n";
    point_fields[3] << "set format cb \'%.03g\'\n";

    point_fields[3] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1\n";

    // 4 dust temp
    point_fields[4] << point_header.str();
    point_fields[4] << "set palette defined (0 0.05 0 0, 0.4 1 0 0, 0.7 1 1 0, 1 1 1 0.5)\n";

    point_fields[4] << "set title \'3D dust temperature distribution (min: " << min_dust_temp
                    << "[K]; max: " << max_dust_temp << "[K])\' font \'Arial,12\'\n";
    point_fields[4] << "set cblabel \'temperature [K]\'\n";

    if(min_dust_temp == 0 && max_dust_temp == 0)
    {
        min_dust_temp = 0.1;
        max_dust_temp = 1;
    }

    if(min_dust_temp / max_dust_temp > 0.9)
        min_dust_temp = 0.9 * max_dust_temp;

    point_fields[4] << "set cbrange[" << float(min_dust_temp) << ":" << float(max_dust_temp) << "]\n";
    point_fields[4] << "set format cb \'%.03g\'\n";

    point_fields[4] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1\n";

    // 5 rat
    point_fields[5] << point_header.str();
    point_fields[5] << "set palette defined (0 0.05 0 0, 0.4 1 0 0, 0.7 1 1 0, 1 1 1 0.5)\n";

    point_fields[5] << "set title \'3D aligned radii distribution (min ID: " << aalg_min
                    << "; max ID: " << aalg_max << ")\' font \'Arial,12\'\n";

    point_fields[5] << "set cblabel \'aligned radius ID\'\n";

    if(aalg_min == aalg_max)
        aalg_max = 1.01 * aalg_min;

    point_fields[5] << "set cbrange[" << aalg_min << ":" << aalg_max << "]\n";

    point_fields[5] << "set format cb \'%.03g\'\n";

    point_fields[5] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1\n";

    vec_header.str("");
    vec_header << "reset\n";
    vec_header << "#set terminal postscript\n";
    vec_header << "#set output \'\'\n";
    vec_header << "set ticslevel 0\n";
    vec_header << "set size ratio -1\n";
    vec_header << "set view 45,45\n";

    vec_header << "set xlabel \'x[m]\'\n";
    vec_header << "set ylabel \'y[m]\'\n";
    vec_header << "set zlabel \'z[m]\'\n";

    vec_header << "set xrange[" << -0.51 * max_len << ":" << 0.51 * max_len << "]\n";
    vec_header << "set yrange[" << -0.51 * max_len << ":" << 0.51 * max_len << "]\n";
    vec_header << "set zrange[" << -0.51 * max_len << ":" << 0.51 * max_len << "]\n";

    vec_header << "set style arrow 1 nohead ls 1 lw 1 lc rgb 0x0000cc\n";
    vec_header << "set style arrow 2 nohead ls 1 lw 1 lc rgb 0x5500dd\n";
    vec_header << "set style arrow 3 ls 1 lw 1 lc palette\n";

    vec_header << "set grid\n";
    vec_header << "set nokey\n";

    // 0 mag
    vec_fields[0] << vec_header.str();
    vec_fields[0] << "set palette defined (0 1 0 0, 0.5 0.0 0.9 0,  0.75 0.0 0.9 1, 0.9 "
                     "0 0.1 0.9)\n";

    if(min_mag == 0 && max_mag == 0)
    {
        min_mag = 1e-45;
        max_mag = 2e-45;
    }

    vec_fields[0] << "set title \'3D mag. field distribution (min:" << log10(min_mag)
                  << " log10([T]); max:" << log10(max_mag) << " log10([T])  \' font \'Arial,12\'\n";

    if(min_mag / max_mag > 0.9)
        min_mag = 0.9 * max_mag;

    vec_fields[0] << "set cbrange[" << log10(min_mag) << ":" << log10(max_mag) << "]\n";
    vec_fields[0] << "set format cb \'%.02g\'\n";
    vec_fields[0] << "splot  \'-\' with vectors as 3, \'-\' with vectors as 2, \'-\' "
                     "with vectors as 1\n";

    // 1 vel
    vec_fields[1] << vec_header.str();
    vec_fields[1] << "set palette defined (0 1 0 0, 0.5 0.0 0.9 0,  0.75 0.0 0.9 1, 0.9 "
                     "0 0.1 0.9)\n";

    if(min_vel == 0 && max_vel == 0)
    {
        min_vel = 1e-45;
        max_vel = 1e-45;
    }

    vec_fields[1] << "set title \'3D vel. field directions (min:" << log10(min_vel)
                  << " log10(m/s); max:" << log10(max_vel) << " log10(m/s)\' font \'Arial,12\'\n";

    if(min_vel / max_vel > 0.9)
        min_vel = 0.9 * max_vel;

    vec_fields[1] << "set cbrange[" << float(log10(min_vel)) << ":" << float(log10(max_vel)) << "]\n";
    vec_fields[1] << "set format cb \'%.03g\'\n";
    vec_fields[1] << "splot  \'-\' with vectors as 3, \'-\' with vectors as 2, \'-\' "
                     "with vectors as 1\n";

    line_counter = 0;

    point_fields[0] << basic_grid.str() << "\ne\n";
    point_fields[0] << tetra_lines.str() << "\ne\n";
    point_fields[0] << basic_grid.str() << "\ne\n";

    for(uint i = 0; i < max_cells; i++)
    {
        // if(line_counter % 200 == 0)
        // {
        //     char_counter++;
        //     cout << "-> Creating plot files : " << float(100.0 * double(line_counter) / double(max_cells))
        //          << "      \r" << flush;
        // }

        // Delaunay lines
        if(line_counter % maxPlotLines == 0)
            addPlotLines(i, tetra_lines);

        const cell_vo * tmp_cell_pos = (const cell_vo *)cell_list[i];
        Vector3D c = getCenter(*tmp_cell_pos);

        line_counter++;

        double p_size = 1.0;
        double v_size = 1e-2;

        if(line_counter % nrOfPlotPoints == 0)
        {
            point_fields[0] << c.X() << " " << c.Y() << " " << c.Z() << "\n";

            if(plt_gas_dens)
            {
                double dens = getGasDensity(*tmp_cell_pos);

                if(dens > 0)
                    point_fields[1] << c.X() << " " << c.Y() << " " << c.Z() << " " << float(p_size) << " "
                                    << log10(dens) << "  \n";
            }

            if(plt_gas_temp)
            {
                double Tg = getGasTemperature(*tmp_cell_pos);
                point_fields[3] << c.X() << " " << c.Y() << " " << c.Z() << " " << float(p_size) << " " << Tg
                                << "  \n";
            }

            if(plt_gas_temp)
            {
                double Td = getDustTemperature(*tmp_cell_pos);
                point_fields[4] << c.X() << " " << c.Y() << " " << c.Z() << " " << float(p_size) << " " << Td
                                << "  \n";
            }

            if(plt_rat)
            {
                double a_alg = getAlignedRadius(*tmp_cell_pos, 0);
                point_fields[5] << c.X() << " " << c.Y() << " " << c.Z() << " " << float(p_size) << " "
                                << a_alg << "  \n";
            }
        }

        if(line_counter % nrOfPlotVectors == 0)
        {
            if(plt_mag)
            {
                double mx = getMagField(*tmp_cell_pos).X();
                double my = getMagField(*tmp_cell_pos).Y();
                double mz = getMagField(*tmp_cell_pos).Z();

                double b_len = sqrt(mx * mx + my * my + mz * mz);

                if(b_len != 0)
                {
                    mx = max_len * v_size * mx / b_len;
                    my = max_len * v_size * my / b_len;
                    mz = max_len * v_size * mz / b_len;

                    vec_fields[0] << float(c.X() - mx) << " " << float(c.Y() - my) << " " << float(c.Z() - mz)
                                  << " " << float(2.0 * mx) << " " << float(2.0 * my) << " "
                                  << float(2.0 * mz) << " " << float(log10(b_len)) << "\n";
                }
            }

            if(plt_vel)
            {
                double vx = getVelocityField(*tmp_cell_pos).X();
                double vy = getVelocityField(*tmp_cell_pos).Y();
                double vz = getVelocityField(*tmp_cell_pos).Z();

                double v_len = sqrt(vx * vx + vy * vy + vz * vz);

                if(v_len != 0)
                {
                    vx = max_len * v_size * vx / v_len;
                    vy = max_len * v_size * vy / v_len;
                    vz = max_len * v_size * vz / v_len;

                    vec_fields[1] << float(c.X() - vx) << " " << float(c.Y() - vy) << " " << float(c.Z() - vz)
                                  << " " << float(2.0 * vx) << " " << float(2.0 * vy) << " "
                                  << float(2.0 * vz) << " " << float(log10(v_len)) << "\n";
                }
            }
        }
    }

    for(uint pos = 1; pos < 6; pos++)
    {
        point_fields[pos] << "\ne\n" << basic_grid.str() << "\ne\n";
    }

    for(uint pos = 0; pos < 2; pos++)
    {
        vec_fields[pos] << "\ne\n" << basic_grid.str() << "\ne\n";
    }

    for(uint pos = 0; pos < 8; pos++)
        point_fields[pos].close();

    for(uint pos = 0; pos < 2; pos++)
        vec_fields[pos].close();

    cout << "- Writing of plot files                : done       \n";
    return true;
}

// saves grid in the POLARIS Voronoi grid file format
bool CGridVoronoi::saveBinaryGridFile(string filename, ushort id, ushort data_size)
{
    if(max_cells == 0)
    {
        cout << ERROR_LINE << "Cannot save Voronoi grid file to:\n";
        cout << filename;
        cout << "Not enough cells available! \n";
        return false;
    }

    ofstream bin_writer(filename.c_str(), ios::out | ios::binary);

    if(bin_writer.fail())
    {
        cout << ERROR_LINE << "Cannot open Voronoi grid file: \n";
        cout << filename;
        return false;
    }

    bin_writer.write((char *)&id, 2);
    bin_writer.write((char *)&data_size, 2);

    if(dataID == GRID_ID_VOR)
    {
        for(uint i = 0; i < data_offset; i++)
        {
            ushort tmp_ids = data_ids[i];
            bin_writer.write((char *)&tmp_ids, 2);
        }
    }
    else
    {
        cout << ERROR_LINE << "Cannot save Voronoi grid file to:\n";
        cout << filename;
        cout << "A Voronoi grid requires an ID of " << GRID_ID_VOR << "!               \n";
        return false;
    }

    double tmp_val = double(max_cells);
    bin_writer.write((char *)&tmp_val, 8);
    bin_writer.write((char *)&max_len, 8);

    line_counter = 0;

    for(ulong c_i = 0; c_i < max_cells; c_i++)
    {
        line_counter++;
        // if(line_counter % 100 == 0)
        // {
        //     char_counter++;
        //     cout << "-> Writing binary Voronoi grid file: "
        //          << float(100.0 * double(line_counter) / double(max_cells)) << "      \r" << flush;
        // }

        cell_vo * tmp_cell = (cell_vo *)cell_list[c_i];

        float tmpCX = (float)tmp_cell->getX();
        float tmpCY = (float)tmp_cell->getY();
        float tmpCZ = (float)tmp_cell->getZ();

        double tmp_vol = tmp_cell->getVolume();

        bin_writer.write((char *)&tmpCX, 4);
        bin_writer.write((char *)&tmpCY, 4);
        bin_writer.write((char *)&tmpCZ, 4);
        bin_writer.write((char *)&tmp_vol, 8);

        for(uint i = 0; i < data_offset; i++)
        {
            float tmp_data = float(tmp_cell->getData(i));
            bin_writer.write((char *)&tmp_data, 4);
        }

        int nr_neighbors = (int)tmp_cell->getNrOfNeighbors();
        uint id = tmp_cell->getUniqueID();
        int tmp_n = 0;

        if(isHullPoint(id))
            nr_neighbors *= -1;

        bin_writer.write((char *)&nr_neighbors, 4);

        for(uint i = 0; i < uint(abs(nr_neighbors)); i++)
        {
            tmp_n = tmp_cell->getNeighborID(i);
            bin_writer.write((char *)&tmp_n, 4);
        }
    }

    bin_writer.close();

    cout << CLR_LINE;
    cout << "- Writing Voronoi grid file            : done     \n";
    return true;
}

// function for debug purposes only
bool CGridVoronoi::createArtificialGrid(string path)
{
    resetGridValues();

    line_counter = 1;
    char_counter = 0;

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

    /*double tmpX[] = {-0.588, 0.79, 0.944, 0.114, 0.946, -0.293, -0.578, -0.98, 0.906,
    -0.387, -0.008, -0.987, 0.144, -0.619, -0.509, -0.159, 0.823, -0.611, 0.35, 0.814,
    0.948, 0.219, 0.802, -0.366, -0.601, -0.693, 0.411, -0.972, 0.413, 0.547, -0.531,
    -0.409, 0.703, 0.748, 0.841, 0.926, -0.288, -0.379, 0.675, 0.777, 0.541, 0.959, 0.042,
    0.938, 0.024, -0.206, -0.913, -0.272, 0.105, -0.431, -0.067, -0.857, 0.421, -0.435,
    -0.789, -0.242, -0.488, -0.115, 0.904, -0.018, 0.667, -0.575, 0.852, 0.646, 0.558,
    0.428, 0.18, 0.528, 0.091, 0.545, -0.335, 0.993, 0.486, -0.814, -0.968, -0.953, 0.024,
    0.052, -0.142, 0.548, 0.692, -0.779, 0.28, -0.71, 0.817, -0.93, 0.811, 0.549, -0.83,
    0.579, -0.244, 0.37, 0.845, -0.968, 0.579, -0.958, -0.15, 0.181, -0.435, -0.909,
    0.307, -0.126}; double tmpY[] = {0.138, 0.311, -0.751, -0.067, -0.857, 0.194, 0.515,
    -0.234, -0.25, -0.212, 0.2, -0.01, -0.067, -0.953, 0.994, -0.3, 0.123, 0.144, 0.049,
    0.662, -0.347, -0.333, 0.811, -0.809, -0.317, 0.525, 0.347, 0.737, 0.087, -0.1,
    -0.153, -0.934, 0.868, -0.499, -0.1, 0.208, -0.382, -0.354, -0.941, 0.606, -0.925,
    -0.334, -0.631, -0.883, 0.958, 0.884, -0.405, 0.216, -0.11, -0.333, -0.793, -0.835,
    0.208, 0.145, -0.831, -0.037, 0.05, 0.779, -0.116, 0.453, -0.614, -0.47, 0.951,
    -0.484, -0.272, 0.486, 0.999, -0.117, -0.213, 0.77, 0.436, 0.862, 0.713, 0.359,
    -0.083, -0.656, 0.525, 0.052, -0.857, 0.688, -0.647, -0.654, 0.832, -0.982, -0.949,
    -0.902, 0.122, -0.826, 0.032, -0.42, -0.729, -0.508, 0.773, -0.937, 0.88, 0.951,
    0.762, 0.311, 0.481, 0.059, 0.815, 0.184}; double tmpZ[] = {0.022, 0.131, 0.237,
    -0.89, -0.65, -0.019, -0.433, -0.039, -0.537, -0.867, 0.272, -0.417, -0.168, 0.766,
    0.909, 0.103, -0.382, 0.385, -0.866, -0.078, -0.878, 0.112, -0.956, -0.908, -0.915,
    -0.214, 0.91, 0.053, -0.898, 0.764, -0.583, 0.084, -0.105, 0.087, 0.77, 0.294, -0.36,
    -0.623, 0.464, -0.637, 0.211, 0.277, 0.145, 0.646, -0.454, 0.755, -0.628, 0.019,
    0.078, -0.316, -0.241, 0.426, 0.795, 0.913, 0.87, 0.753, -0.653, 0.564, -0.547, 0.12,
    0.328, 0.493, -0.299, 0.846, -0.986, -0.023, -0.973, 0.807, -0.31, -0.857, -0.471,
    0.059, 0.01, 0.527, -0.746, 0.74, 0.738, -0.379, 0.408, 0.178, -0.387, -0.525, -0.813,
    0.547, -0.537, -0.193, -0.669, 0.384, 0.527, 0.215, 0.568, 0.914, 0.746, 0.764, 0.815,
    0.913, 0.035, 0.209, -0.759, -0.178, 0.094, 0.587}; max_cells = 102; */

    /*double tmpX[] = {-0.14, 0.91, 0.87, -0.02, 0.16, 0.98, -0.62, -0.42, 0.59, -0.29,
    0.89, -0.29, -0.59, -0.03, 0.06, 0.00, 0.73, 0.31, 0.93, -0.84, 0.76, 0.17, -0.06,
    0.44, -0.36, 0.05, -0.56, -0.85, -0.61}; double tmpY[] = {0.25, 0.30, 0.77, 0.80,
    -0.21, 0.62, 0.12, -0.62, 0.19, 0.83, 0.49, -0.97, 0.95, 0.15, -0.43, -0.79, -0.57,
    -0.02, -0.61, 0.88, -0.24, -0.05, -0.21, -0.20, 0.07, 0.85, -0.21, -0.53, -0.11};
    double tmpZ[] = {0.53, 0.28, 0.67, 0.70, -0.39, 0.28, -0.52, 0.62, 0.73, 0.70, -0.55,
    0.00, -0.23, 0.56, -0.88, -0.61, -0.11, 0.00, -0.58, 0.81, 0.16, 0.35, 0.92, 0.36,
    0.78, 0.89, 0.49, -0.15, 0.15}; max_cells = 29; */

    /*double tmpX[] = {-0.55, +0.55,  0   , 0  , 0  , 0};
    double tmpY[] = {+0.55, +0.55, -0.55, 0  , 0  , 0};
    double tmpZ[] = {-0.55, -0.55, -0.55,-0.8, 0.8, 0};
    max_cells = 6;*/

    /*double tmpX[] = {-0.55, +0.55, 0, 0, 0, -0.85, 0};
    double tmpY[] = {+0.55, +0.55, -0.55, 0, 0, -0.85, 0};
    double tmpZ[] = {-0.55, -0.55, -0.55, 0.8, -0.8, -0.85, 0};
    max_cells = 7; */

    double tmpX[] = { -0.25, +0.25, 0, 0, 0 };
    double tmpY[] = { +0.25, +0.25, -0.25, 0, 0 };
    double tmpZ[] = { -0.25, -0.25, -0.25, 0.25, 0.9 };
    max_cells = 5;

    /*double tmpX[] = { 0.0, 0.0, 0.0, 0.15, -0.15};
    double tmpY[] = {-0.5, 0.5, 0.0, 0.0, 0.0};
    double tmpZ[] = {-0.8,-0.8, 0.8, -0.1, -0.1};
    max_cells = 5; */

    /*double tmpZ[] = {0.0, 0.0, 0.0, 0.15, -0.15};
    double tmpY[] = {-0.5, 0.5, 0.0, 0.0, 0.0};
    double tmpX[] = {-0.8, -0.8, 0.8, -0.1, -0.1};
    max_cells = 5; */

    /* double tmpX[] = {0.0, 0.5, 0.0, 0.0, 0.5};
    double tmpY[] = {0.0, 0.0, 0.5, 0.0, 0.5};
    double tmpZ[] = {0.0, 0.0, 0.0, 0.5, 0.5};
    max_cells = 5;*/

    /*double tmpX[] = {-1, 0, 1};
    double tmpY[] = {-1, 0, 1};
    double tmpZ[] = {-1, 0, 1};
    max_cells = 3;*/

    min_len = 1e30;
    max_len = 2;

    // max_cells/=10;

    string filename = path;

    dataID = GRID_ID_VOR;

    ofstream bin_writer(filename.c_str(), ios::out | ios::binary);

    if(bin_writer.fail())
    {
        cout << ERROR_LINE << "Cannot open file: \n";
        cout << filename;
        return false;
    }

    bin_writer.write((char *)&dataID, 2);
    bin_writer.write((char *)&data_offset, 2);

    ushort tmp_ids;

    tmp_ids = GRIDgas_dens;
    bin_writer.write((char *)&tmp_ids, 2);

    tmp_ids = GRIDdust_temp;
    bin_writer.write((char *)&tmp_ids, 2);

    tmp_ids = GRIDgas_temp;
    bin_writer.write((char *)&tmp_ids, 2);

    tmp_ids = GRIDmx;
    bin_writer.write((char *)&tmp_ids, 2);

    tmp_ids = GRIDmy;
    bin_writer.write((char *)&tmp_ids, 2);

    tmp_ids = GRIDmz;
    bin_writer.write((char *)&tmp_ids, 2);

    tmp_ids = GRIDvx;
    bin_writer.write((char *)&tmp_ids, 2);

    tmp_ids = GRIDvy;
    bin_writer.write((char *)&tmp_ids, 2);

    tmp_ids = GRIDvz;
    bin_writer.write((char *)&tmp_ids, 2);

    bin_writer.write((char *)&max_cells, 4);

    bin_writer.write((char *)&max_len, 8);

    for(uint i = 0; i < max_cells; i++)
    {
        // if(i % 50 == 0)
        //     cout << "-> Creating grid: " << 100.0 * float(i) / float(max_cells) << " [%]           \r"
        //          << flush;

        double tmpC;

        tmpC = tmpX[i];
        bin_writer.write((char *)&tmpC, 8);

        tmpC = tmpY[i];
        bin_writer.write((char *)&tmpC, 8);

        tmpC = tmpZ[i];
        bin_writer.write((char *)&tmpC, 8);

        double tmp_data = i;

        bin_writer.write((char *)&tmp_data, 8);

        tmp_data = 10;
        bin_writer.write((char *)&tmp_data, 8);

        tmp_data = 20;
        bin_writer.write((char *)&tmp_data, 8);

        for(uint i = 3; i < 9; i++)
        {
            tmp_data = double(10);
            bin_writer.write((char *)&tmp_data, 8);
        }
    }

    bin_writer.close();

    cout << "min: " << min_gas_dens << "  max_dens: " << max_gas_dens << "   \n";
    cout << "Creating artificial Voronoi grid                  : done    \n";
    cout << "Max cells: " << max_cells << "   \n";

    return true;
}

void CGridVoronoi::printParameters()
{
    if(max_cells == 0)
    {
        cout << ERROR_LINE << "No Voronoi cells available! \n";
    }
    else
    {
        cout << CLR_LINE;
        cout << SEP_LINE;
        cout << "Voronoi grid parameters (ID: " << getDataID() << ", data len.:  " << getDataOffset()
             << ")         \n";
        cout << SEP_LINE;

        if(stree != 0)
        {
            cout << "- Number of search tree depth      : " << stree->getMaxLevel() << " \n";
            cout << "- Number of total tree nodes       : " << stree->getMaxNodes() << " \n";
        }

        cout << "- Number of Voronoi cells          : " << max_cells << " \n";
        cout << "- Number of convex hull points     : " << hull_size << " \n";
        cout << "- Number of neighbors    (min,max) : [" << min_nrOfNeigbors << ", " << max_nrOfNeigbors
             << "]   \n";
        cout << SEP_LINE;

        printPhysicalParameters();
        cout << SEP_LINE;
    }
}

// brute force Voronoi position finder for debug purposes only
bool CGridVoronoi::positionPhotonInGridTest(photon_package * pp)
{
    Vector3D pos = pp->getPosition();
    Vector3D tmp_pos;
    double min_l = -0.5 * max_len;
    double max_l = 0.5 * max_len;

    if(pos.X() < min_l)
        return false;
    if(pos.Y() < min_l)
        return false;
    if(pos.Z() < min_l)
        return false;

    if(pos.X() > max_l)
        return false;
    if(pos.Y() > max_l)
        return false;
    if(pos.Z() > max_l)
        return false;

    double min_radius = 1e300;

    bool found = false;

    for(uint i = 0; i < max_cells; i++)
    {
        cell_vo * tmp_cell = ((cell_vo *)cell_list[i]);
        Vector3D dist = tmp_cell->getCenter() - pos;
        double len = dist.sq_length();

        if(len < min_radius)
        {
            found = true;
            min_radius = len;
            pp->setPositionCell(tmp_cell);
        }
    }

    return found;
}

bool CGridVoronoi::positionPhotonInGrid(photon_package * pp)
{
    Vector3D pos = pp->getPosition();
    double min_l = -0.5 * max_len;
    double max_l = 0.5 * max_len;

    if(pos.X() < min_l)
        return false;
    if(pos.Y() < min_l)
        return false;
    if(pos.Z() < min_l)
        return false;

    if(pos.X() > max_l)
        return false;
    if(pos.Y() > max_l)
        return false;
    if(pos.Z() > max_l)
        return false;

    cell_vo * cell = stree->findClosestCell(pos, cell_list);

    if(cell == 0)
    {
        cout << ERROR_LINE << "Photon package cannot be positioned in Voronoi grid!           "
                "             \n";
        return false;
    }

    pp->setPositionCell((cell_basic *)cell);
    return true;
}

bool CGridVoronoi::goToNextCellBorder(photon_package * pp)
{
    bool hit = false;

    double path_length = 2e300;

    double min_l = -0.5 * max_len;
    double max_l = 0.5 * max_len;

    Vector3D pos = pp->getPosition();
    Vector3D dir = pp->getDirection();

    // length_eps is the minimum step width to ensure that
    // the photon 1) moves and 2) enters the cell !numerically!
    double length_eps_1, length_eps_2;

    cell_vo * center_cell = (cell_vo *)pp->getPositionCell();

    if(center_cell == 0)
        return false;

    uint n_size = center_cell->getNrOfNeighbors();
    Vector3D c_pos = center_cell->getCenter();

    //uint counter = 0;

    if(n_size == 0)
    {
        double volume = center_cell->getVolume();
        double ref_length = pow(3.0 * volume / PIx4, 1.0 / 3.0);
    
        path_length = ref_length;
        hit = true;
    }

    Vector3D v_n, v_a;
    double length, num, den, cell_distance;

    for(uint i = 0; i < n_size; i++)
    {
        if(isNeigboringVoroCell(center_cell, i))
        {
            uint id = center_cell->getNeighborID(i);

            cell_vo * n_cell = ((cell_vo *)cell_list[id]);
            Vector3D n_pos = n_cell->getCenter();

            // v_n is normal vector on the cell border pointing
            // towards next cell
            v_n = n_pos - c_pos;
            cell_distance = v_n.length();
            // v_n should be unit vector
            v_n.normalize();

            // den = cos of angle between cell border normal and photon direction
            den = v_n * dir;

            // den must be positive as long as v_n points outwards
            if(den > 0)
            {
                // v_a is a point on the cell border and
                // on the line between the two center points
                v_a = c_pos + 0.5 * cell_distance * v_n;

                // geometrically, abs(num) is the shortest distance from current
                // position to the cell border (perp to border, ie. parallel to v_n)
                // if num is 0 -> photon is on the border
                num = v_n * (pos - v_a);

                // distance num to border is enlarged to ensure that at least one
                // component of the photon position changes parallel to v_n after step
                // sign(num) is necessary to ensure that abs(num) gets larger
                length_eps_1 = abs(pos * v_n) * MIN_LEN_STEP * EPS_DOUBLE;
                num += Vector3D::sign(num) * length_eps_1;

                length = -num / den;

                if(length > 0 && length < path_length)
                {
                    hit = true;
                    length_eps_2 = abs( (pos + dir * length) * v_n ) / den * MIN_LEN_STEP*EPS_DOUBLE;
                    path_length = length + length_eps_2;
                }
            }
        }
    }

    for(int i_side = 0; i_side < 6; i_side++)
    {
        // v_n points outside of current cell
        v_n = 0;
        // v_a is a point on the cell border
        v_a = 0;
        switch(i_side)
        {
            case 0:
                v_n.setZ(-1);
                v_a.setZ(min_l);
                break;
            case 1:
                v_n.setZ(1);
                v_a.setZ(max_l);
                break;
            case 2:
                v_n.setY(-1);
                v_a.setY(min_l);
                break;
            case 3:
                v_n.setY(1);
                v_a.setY(max_l);
                break;
            case 4:
                v_n.setX(-1);
                v_a.setX(min_l);
                break;
            case 5:
                v_n.setX(1);
                v_a.setX(max_l);
                break;
        }
        // den = cos of angle between cell border normal and photon direction
        den = v_n * dir;

        // den must be positive as long as v_n points outwards
        if(den > 0)
        {
            // geometrically, abs(num) is the shortest distance from current
            // position to the cell border (perp to border, ie. parallel to v_n)
            // if num is 0 -> photon is on the border
            num = v_n * (pos - v_a);

            // distance num to border is enlarged to ensure that at least one
            // component of the photon position changes parallel to v_n after step
            // sign(num) is necessary to ensure that abs(num) gets larger
            length_eps_1 = abs(pos * v_n) * MIN_LEN_STEP * EPS_DOUBLE;
            num += Vector3D::sign(num) * length_eps_1;

            length = -num / den;

            if(length > 0 && length < path_length)
            {
                hit = true;
                length_eps_2 = abs( (pos + dir * length) * v_n ) / den * MIN_LEN_STEP*EPS_DOUBLE;
                path_length = length + length_eps_2;
            }
        }
    }

    pp->setPosition(pos + dir * path_length);

    if(pos == pp->getPosition())
    {
        cout << ERROR_LINE << "Could not transfer photon to the next cell border!   " << endl;
        return false;
    }

    pp->setTmpPathLength(path_length);

    return hit;
}

bool CGridVoronoi::updateShortestDistance(photon_package * pp)
{
    Vector3D tmp_pos;
    double min_dist, tmp_dist[6];

    double loc_x_min, loc_x_max, loc_y_min, loc_y_max, loc_z_min, loc_z_max;
    bool found = false;

    cell_oc * tmp_cell_pos = (cell_oc *)pp->getPositionCell();

    tmp_pos = pp->getPosition();

    loc_x_min = tmp_cell_pos->getXmin();
    loc_y_min = tmp_cell_pos->getYmin();
    loc_z_min = tmp_cell_pos->getZmin();

    loc_x_max = tmp_cell_pos->getXmax();
    loc_y_max = tmp_cell_pos->getYmax();
    loc_z_max = tmp_cell_pos->getZmax();

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

    // pp->setShortestDistance(min_dist);
    return found;
}

bool CGridVoronoi::findStartingPoint(photon_package * pp)
{
    bool hit = false;

    double path_length = 0;

    double min_l = -0.5 * max_len;
    double max_l = 0.5 * max_len;

    Vector3D dir = pp->getDirection();
    Vector3D pos = pp->getPosition();

    if(isInside(pos))
        return true;

    Vector3D v_n, v_a;
    double length, num, den;

    // length_eps is the minimum step width to ensure that
    // the photon 1) moves and 2) enters the cell !numerically!
    double length_eps_1, length_eps_2;

    for(int i_side = 0; i_side < 6; i_side++)
    {
        // v_n points inside the neighboring cells
        // photon is outside of the grid
        // -> v_n has different sign compared to goToNextCellBorder
        v_n = 0;
        // v_a is a point on the cell border
        v_a = 0;
        switch(i_side)
        {
            case 0:
                v_n.setZ(1);
                v_a.setZ(min_l);
                break;
            case 1:
                v_n.setZ(-1);
                v_a.setZ(max_l);
                break;
            case 2:
                v_n.setY(1);
                v_a.setY(min_l);
                break;
            case 3:
                v_n.setY(-1);
                v_a.setY(max_l);
                break;
            case 4:
                v_n.setX(1);
                v_a.setX(min_l);
                break;
            case 5:
                v_n.setX(-1);
                v_a.setX(max_l);
                break;
        }
        // den = cos of angle between cell border normal and photon direction
        den = v_n * dir;

        // den is positive (negative) if v_n points away from (towards) photon
        if(den != 0)
        {
            // geometrically, abs(num) is the shortest distance from current
            // position to the cell border (perp to border, ie. parallel to v_n)
            // if num is 0 -> photon is on the border
            num = v_n * (pos - v_a);

            // distance num to border is enlarged to ensure that at least one
            // component of the photon position changes parallel to v_n after step
            // sign(num) is necessary to ensure that abs(num) gets larger
            length_eps_1 = abs(pos * v_n) * MIN_LEN_STEP * EPS_DOUBLE;
            num += Vector3D::sign(num) * length_eps_1;

            length = -num / den;

            if(length > 0 && isInside(pos + dir * length))
            {
                hit = true;
                length_eps_2 = abs( (pos + dir * length) * v_n ) / den * MIN_LEN_STEP*EPS_DOUBLE;
                path_length = length + length_eps_2;
                break;
            }
        }
    }

    pp->setPosition(pos + dir * path_length);
    pp->setTmpPathLength(0);

    return positionPhotonInGrid(pp);
}

bool CGridVoronoi::isInside(const Vector3D & pos) const
{
    double l_min = -0.5 * max_len;
    double l_max = 0.5 * max_len;

    if(pos.X() < l_min)
        return false;
    if(pos.Y() < l_min)
        return false;
    if(pos.Z() < l_min)
        return false;

    if(pos.X() > l_max)
        return false;
    if(pos.Y() > l_max)
        return false;
    if(pos.Z() > l_max)
        return false;

    return true;
}

Vector3D CGridVoronoi::getCenter(const cell_basic & cell) const
{
    return ((const cell_vo *)&cell)->getCenter();
}

Vector3D CGridVoronoi::getCenter(uint id) const
{
    const cell_vo * cell = ((const cell_vo *)cell_list[id]);
    return cell->getCenter();
}

bool CGridVoronoi::next(photon_package * pp)
{
    if(!positionPhotonInGrid(pp))
        return false;

    if(!goToNextCellBorder(pp))
        return false;

    return true;
}

void CGridVoronoi::getLengths(uint bins, double & step_xy, double & off_xy)
{
    step_xy = 2 * max_len / double(bins);
    off_xy = step_xy / 2.0;
}

bool CGridVoronoi::createCellList()
{
    // cout << CLR_LINE;
    // cout << "- Creating of cell list                : done          \n" << flush;
    return true;
}

double CGridVoronoi::getVolume(const cell_basic & cell) const
{
    const cell_vo * cell_pos = (const cell_vo *)&cell;
    return cell_pos->getVolume();
}

bool CGridVoronoi::saveBinaryGridFile(string filename)
{
    return saveBinaryGridFile(filename, GRID_ID_VOR, data_offset);
}

bool CGridVoronoi::loadGridFromBinaryFile(parameters & param)
{
    return loadGridFromBinaryFile(param, 0);
}

// final cleanup
void CGridVoronoi::clear()
{
    line_counter = 0;
    char_counter = 0;

    if(cell_list != 0)
    {
        delete[] cell_list;
        cell_list = 0;
    }

    if(hull_list != 0)
    {
        delete[] hull_list;
        hull_list = 0;
    }

    cout << CLR_LINE;
    cout << "Final cleanup                                : done     \n";
}

double CGridVoronoi::abs_min(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max)
{
    double res = 1e300;

    if(res > abs(x_min))
        res = abs(x_min);

    if(res > abs(x_max))
        res = abs(x_max);

    if(res > abs(y_min))
        res = abs(y_min);

    if(res > abs(y_max))
        res = abs(y_max);

    if(res > abs(z_min))
        res = abs(z_min);

    if(res > abs(z_max))
        res = abs(z_max);

    return res;
}

bool CGridVoronoi::isHullPoint(uint id)
{
    uint N = hull_size;
    uint min = 0, max = N - 1;

    if(id < hull_list[min].id)
        return false;

    if(id > hull_list[max].id)
        return false;

    while(max - min > 1)
    {
        uint i = min + (max - min) / 2;
        if(hull_list[i].id > id)
            max = i;
        else
            min = i;
    }

    if(id == hull_list[min].id)
        return true;

    uint upper = min + 1;

    if(upper > N)
        upper = N;

    if(id == hull_list[upper].id)
        return true;

    uint lower = min - 1;

    if(lower == uint(-1))
        lower = 0;

    if(id == hull_list[lower].id)
        return true;

    return false;
}

void CGridVoronoi::addPlotLines(uint cID, stringstream & str)
{
    cell_vo * tmp_cell = (cell_vo *)cell_list[cID];
    uint nr_neighbors = tmp_cell->getNrOfNeighbors();

    Vector3D p1 = tmp_cell->getCenter();

    for(uint i = 0; i < nr_neighbors; i++)
    {
        if(isNeigboringVoroCell(tmp_cell, i))
        {
            Vector3D p2 = getNeighborCenter(tmp_cell, i);
            Vector3D tmp_len = p2 - p1;
            str << p1.X() << " " << p1.Y() << " " << p1.Z() << " " << tmp_len.X() << " " << tmp_len.Y()
                << " " << tmp_len.Z() << "\n";
        }
    }
}

bool CGridVoronoi::isNeigboringVoroCell(cell_vo * cell, uint nID)
{
    int id = cell->getNeighborID(nID);

    if(id > int(max_cells))
        return false;

    return id > -1;
}

Vector3D CGridVoronoi::getNeighborCenter(cell_vo * cell, uint nID)
{
    int id = cell->getNeighborID(nID);
    cell_vo * n_cell = ((cell_vo *)cell_list[id]);
    return n_cell->getCenter();
}

bool CGridVoronoi::search_tree::tree_node::nodeIntersection(Vector3D point, double _sq_distance)
{
    double X = point.X();
    double Y = point.Y();
    double Z = point.Z();

    // is point within node?
    if(isPointInNode(point))
        return true;

    // edge intersection
    for(uint ix = 0; ix <= 1; ix++)
        for(uint iy = 0; iy <= 1; iy++)
            for(uint iz = 0; iz <= 1; iz++)
            {
                double X1 = x_min + double(ix) * length;
                double Y1 = y_min + double(iy) * length;
                double Z1 = z_min + double(iz) * length;

                double sq_distance =
                    (X - X1) * (X - X1) + (Y - Y1) * (Y - Y1) + (Z - Z1) * (Z - Z1);

                if(sq_distance <= _sq_distance)
                    return true;
            }

    // surface intersection
    double distance = sqrt(_sq_distance);
    for(int ix = -1; ix <= 1; ix++)
        for(int iy = -1; iy <= 1; iy++)
            for(int iz = -1; iz <= 1; iz++)
            {
                if(ix + iy + iz == 0)
                    continue;

                double x = point.X() + double(ix) * distance;
                double y = point.Y() + double(iy) * distance;
                double z = point.Z() + double(iz) * distance;

                if(isPointInNode(Vector3D(x, y, z)))
                    return true;
            } /**/

    return false;
}

bool CGridVoronoi::search_tree::tree_node::isPointInNode(Vector3D point)
{
    if(point.X() < getXMin() || point.Y() < getYMin() || point.Z() < getZMin())
        return false;

    if(point.X() > getXMax() || point.Y() > getYMax() || point.Z() > getZMax())
        return false;

    return true;
}

void CGridVoronoi::search_tree::tree_node::increment()
{
    size++;
}

void CGridVoronoi::search_tree::tree_node::add_cell(cell_vo * cell)
{
    list_element * new_element = new list_element();
    new_element->cell = cell;

    if(first == 0)
    {
        first = new_element;
        last = new_element;
    }
    else
    {
        last->next = new_element;
        last = new_element;
    }

    size++;
}

void CGridVoronoi::search_tree::tree_node::clear()
{
    if(first == 0)
        return;

    list_element * pos = first;

    while(pos != 0)
    {
        list_element * tmp_element = pos;
        pos = pos->next;

        // only for debugging
        // final cell cleanup happens in Voronoi grid itself
        // delete tmp_element->cell;
        // tmp_element->cell=0;

        delete tmp_element;
        tmp_element = 0;
    }

    first = 0;
    last = 0;
    size = 0;
}

void CGridVoronoi::search_tree::tree_node::push_down()
{
    list_element * pos = first;

    while(pos != 0)
    {
        list_element * tmp_element = pos;
        pos = pos->next;

        Vector3D center = tmp_element->cell->getCenter();
        double X = center.X();
        double Y = center.Y();
        double Z = center.Z();

        double xmid = this->getXCenter();
        double ymid = this->getYCenter();
        double zmid = this->getZCenter();

        tree_node * leaf;

        if(Z < zmid) // z 0 1 2 3
        {
            if(Y < ymid) // y 0 1
            {
                if(X < xmid) // x 0
                    leaf = this->getLeaf(0);
                else
                    // x 1
                    leaf = this->getLeaf(1);
            }
            else // y 2 3
            {
                if(X < xmid) // x 2
                    leaf = this->getLeaf(2);
                else // x 3
                    leaf = this->getLeaf(3);
            }
        }
        else // z 4 5 6 7
        {
            if(Y < ymid) // y 4 5
            {
                if(X < xmid) // x 4
                    leaf = this->getLeaf(4);
                else // x 5
                    leaf = this->getLeaf(5);
            }
            else // y 6 7
            {
                if(X < xmid) // x 6
                    leaf = this->getLeaf(6);
                else // x 7
                    leaf = this->getLeaf(7);
            }
        }

        leaf->add_cell(tmp_element->cell);

        delete tmp_element;
        tmp_element = 0;
    }

    first = 0;
    last = 0;
}

cell_vo * CGridVoronoi::search_tree::tree_node::findClosestCell(Vector3D point, double & _min_distance, tree_node *& f_node)
{
    double min_distance = 1e200;
    cell_vo * res = 0;

    // search in highest level
    if(this->getLeafs() != 0)
    {
        for(uint i = 0; i < 8; i++)
        {
            double tmp_distance = 0;
            cell_vo * cell = leafs[i].findClosestCell(point, tmp_distance, f_node);

            if(tmp_distance < min_distance)
            {
                min_distance = tmp_distance;
                res = cell;
            }
        }
    }
    else // search in level plus all sub-levels within search radius
    {
        list_element * pos = first;

        while(pos != 0)
        {
            cell_vo * cell = pos->cell;
            Vector3D center = cell->getCenter();

            double X = point.X();
            double Y = point.Y();
            double Z = point.Z();

            double X1 = center.X();
            double Y1 = center.Y();
            double Z1 = center.Z();

            double sq_distance = (X - X1) * (X - X1) + (Y - Y1) * (Y - Y1) + (Z - Z1) * (Z - Z1);

            if(sq_distance < min_distance)
            {
                min_distance = sq_distance;
                res = cell;

                f_node = this;
            }

            pos = pos->next;
        }
    }

    _min_distance = min_distance;
    return res;
}

uint CGridVoronoi::search_tree::tree_node::get_size()
{
    return size;
}

bool CGridVoronoi::search_tree::tree_node::is_emty()
{
    return size == 0;
}

double CGridVoronoi::search_tree::tree_node::getXMin()
{
    return x_min;
}

double CGridVoronoi::search_tree::tree_node::getYMin()
{
    return y_min;
}

double CGridVoronoi::search_tree::tree_node::getZMin()
{
    return z_min;
}

double CGridVoronoi::search_tree::tree_node::getXMax()
{
    return x_min + length;
}

double CGridVoronoi::search_tree::tree_node::getYMax()
{
    return y_min + length;
}

double CGridVoronoi::search_tree::tree_node::getZMax()
{
    return z_min + length;
}

double CGridVoronoi::search_tree::tree_node::getXCenter()
{
    return x_min + 0.5 * length;
}

double CGridVoronoi::search_tree::tree_node::getYCenter()
{
    return y_min + 0.5 * length;
}

double CGridVoronoi::search_tree::tree_node::getZCenter()
{
    return z_min + 0.5 * length;
}

void CGridVoronoi::search_tree::tree_node::setXMin(double x)
{
    x_min = x;
}

void CGridVoronoi::search_tree::tree_node::setYMin(double y)
{
    y_min = y;
}

void CGridVoronoi::search_tree::tree_node::setZMin(double y)
{
    z_min = y;
}

double CGridVoronoi::search_tree::tree_node::getLength() const
{
    return length;
}

void CGridVoronoi::search_tree::tree_node::setLength(double l)
{
    length = l;
}

void CGridVoronoi::search_tree::tree_node::setLevel(uint l)
{
    level = l;
}

CGridVoronoi::search_tree::list_element * CGridVoronoi::search_tree::tree_node::get_first()
{
    return first;
}

CGridVoronoi::search_tree::list_element * CGridVoronoi::search_tree::tree_node::get_last()
{
    return last;
}

uint CGridVoronoi::search_tree::tree_node::getLevel()
{
    return level;
}

CGridVoronoi::search_tree::tree_node * CGridVoronoi::search_tree::tree_node::getLeafs()
{
    return leafs;
}

CGridVoronoi::search_tree::tree_node * CGridVoronoi::search_tree::tree_node::getLeaf(uint index)
{
    return &leafs[index];
}

CGridVoronoi::search_tree::tree_node * CGridVoronoi::search_tree::tree_node::getBranch()
{
    return branch;
}

void CGridVoronoi::search_tree::tree_node::setLeafs(tree_node * l)
{
    leafs = l;
}

void CGridVoronoi::search_tree::tree_node::setBranch(tree_node * b)
{
    branch = b;
}

bool CGridVoronoi::search_tree::addCell(cell_vo * cell)
{
    Vector3D center = cell->getCenter();

    if(center.X() < root->getXMin() || center.Y() < root->getYMin() || center.Z() < root->getZMin())
        return false;

    if(center.X() > root->getXMax() || center.Y() > root->getYMax() || center.Z() > root->getZMax())
        return false;

    return addCell(root, cell);
}

void CGridVoronoi::search_tree::initTree(double _side_length)
{
    side_length = _side_length;

    root = new tree_node();
    root->setLength(side_length);
    root->setXMin(-0.5 * side_length);
    root->setYMin(-0.5 * side_length);
    root->setZMin(-0.5 * side_length);

    max_nodes = 1;

    createLeafNodes(root);
}

uint CGridVoronoi::search_tree::getMaxLevel()
{
    return max_level;
}

uint CGridVoronoi::search_tree::getMaxNodes()
{
    return max_nodes;
}

cell_vo * CGridVoronoi::search_tree::findClosestCell(Vector3D point, cell_basic ** cell_list)
{
    // find current node that contains the point
    tree_node * p_node = findMatchingNode(point, MAX_LEVEL);
    double min_distance = 0;

    // no node found
    if(p_node == 0)
        return 0;

    // when node has not enough cells go to next lower level
    while(p_node->is_emty())
    {
        p_node = p_node->getBranch();

        if(p_node->getBranch() == 0)
            break;
    }

    // node that contains the closest cell
    // not necessarily identical with p_node
    tree_node * f_node = 0;

    // find closest cell in current node
    cell_vo * cell = p_node->findClosestCell(point, min_distance, f_node);

    // find possible closer cells in neighboring nodes
    double len = min_distance * 1.00001;
    double n_distance = 1e200;
    cell_vo * n_cell = checkNeighboringNodes(root, p_node, f_node, point, len, n_distance);

    if(n_distance < min_distance)
    {
        cell = n_cell;
        min_distance = n_distance;
    }

    return cell;
}

cell_vo * CGridVoronoi::search_tree::checkNeighboringNodes(tree_node * node,
                                                           tree_node * f_node,
                                                           tree_node * p_node,
                                                           Vector3D point,
                                                           double distance,
                                                           double & min_distance)
{
    cell_vo * res_cell = 0;
    if(node->getLeafs() == 0)
    {
        // do not search f_node again
        if(node != f_node)
        {
            // do not search p_node again
            if(node != p_node)
            {
                tree_node * dummy;
                res_cell = node->findClosestCell(point, min_distance, dummy);
            }
        }
    }
    else
    {
        for(uint i = 0; i < 8; i++)
        {
            if(node->getLeaf(i)->nodeIntersection(point, distance))
            {
                double tmp_distance = 1e200;
                cell_vo * tmp_cell = checkNeighboringNodes(
                    node->getLeaf(i), f_node, p_node, point, distance, tmp_distance);

                if(tmp_distance < min_distance)
                {
                    min_distance = tmp_distance;
                    res_cell = tmp_cell;
                }
            }
        }
    }

    return res_cell;
}

CGridVoronoi::search_tree::tree_node * CGridVoronoi::search_tree::findMatchingNode(Vector3D point)
{
    return findMatchingNode(point, MAX_LEVEL);
}

CGridVoronoi::search_tree::tree_node * CGridVoronoi::search_tree::findMatchingNode(Vector3D point, uint _level)
{
    if(point.X() < root->getXMin() || point.Y() < root->getYMin() || point.Z() < root->getZMin())
        return 0;

    if(point.X() > root->getXMax() || point.Y() > root->getYMax() || point.Z() > root->getZMax())
        return 0;

    return goNextLevelDown(root, point, _level);
}

bool CGridVoronoi::search_tree::addCell(tree_node * node, cell_vo * cell)
{
    Vector3D center = cell->getCenter();
    double X = center.X();
    double Y = center.Y();
    double Z = center.Z();

    if(node->getLeafs() == 0)
    {
        node->add_cell(cell);

        // check for maximal cells in node
        if(node->get_size() >= MAX_CELLS)
        {
            // check for maximal level
            if(node->getLevel() < MAX_LEVEL)
            {
                // do refinement
                createLeafNodes(node);
                node->push_down();
            }
            return true;
        }

        return true;
    }

    node->increment();

    double xmid = node->getXCenter();
    double ymid = node->getYCenter();
    double zmid = node->getZCenter();

    if(Z < zmid) // z 0 1 2 3
    {
        if(Y < ymid) // y 0 1
        {
            if(X < xmid) // x 0
                node = node->getLeaf(0);
            else
                // x 1
                node = node->getLeaf(1);
        }
        else // y 2 3
        {
            if(X < xmid) // x 2
                node = node->getLeaf(2);
            else // x 3
                node = node->getLeaf(3);
        }
    }
    else // z 4 5 6 7
    {
        if(Y < ymid) // y 4 5
        {
            if(X < xmid) // x 4
                node = node->getLeaf(4);
            else // x 5
                node = node->getLeaf(5);
        }
        else // y 6 7
        {
            if(X < xmid) // x 6
                node = node->getLeaf(6);
            else // x 7
                node = node->getLeaf(7);
        }
    }

    return addCell(node, cell);
}

CGridVoronoi::search_tree::tree_node * CGridVoronoi::search_tree::goNextLevelDown(tree_node * node, Vector3D point, uint _level)
{
    if(node->getLeafs() == 0)
        return node;

    if(node->getLevel() == _level)
        return node;

    double xmid = node->getXCenter();
    double ymid = node->getYCenter();
    double zmid = node->getZCenter();

    if(point.Z() < zmid) // z 0 1 2 3
    {
        if(point.Y() < ymid) // y 0 1
        {
            if(point.X() < xmid) // x 0
                node = node->getLeaf(0);
            else
                // x 1
                node = node->getLeaf(1);
        }
        else // y 2 3
        {
            if(point.X() < xmid) // x 2
                node = node->getLeaf(2);
            else // x 3
                node = node->getLeaf(3);
        }
    }
    else // z 4 5 6 7
    {
        if(point.Y() < ymid) // y 4 5
        {
            if(point.X() < xmid) // x 4
                node = node->getLeaf(4);
            else // x 5
                node = node->getLeaf(5);
        }
        else // y 6 7
        {
            if(point.X() < xmid) // x 6
                node = node->getLeaf(6);
            else // x 7
                node = node->getLeaf(7);
        }
    }

    return goNextLevelDown(node, point, _level);
}

void CGridVoronoi::search_tree::clear(tree_node * node)
{
    tree_node * leafs = node->getLeafs();

    node->clear();

    if(leafs == 0)
        return;

    for(uint i = 0; i < 8; i++)
        clear(&leafs[i]);

    delete[] leafs;
    node->setLeafs(0);
}

// create eight sub-nodes for node pointer
bool CGridVoronoi::search_tree::createLeafNodes(tree_node * node)
{
    uint next_level = 1 + node->getLevel();

    if(MAX_LEVEL < next_level)
        return false;

    max_level = next_level;
    max_nodes += 4;

    double ox = node->getXMin();
    double oy = node->getYMin();
    double oz = node->getZMin();

    double tmp_length = 0.5 * node->getLength();
    tree_node * leafs = new tree_node[8];

    leafs[0].setBranch(node);
    leafs[0].setXMin(ox);
    leafs[0].setYMin(oy);
    leafs[0].setZMin(oz);

    leafs[0].setLength(tmp_length);
    leafs[0].setLevel(next_level);

    leafs[1].setBranch(node);
    leafs[1].setXMin(ox + tmp_length);
    leafs[1].setYMin(oy);
    leafs[1].setZMin(oz);

    leafs[1].setLength(tmp_length);
    leafs[1].setLevel(next_level);

    leafs[2].setBranch(node);
    leafs[2].setXMin(ox);
    leafs[2].setYMin(oy + tmp_length);
    leafs[2].setZMin(oz);

    leafs[2].setLength(tmp_length);
    leafs[2].setLevel(next_level);

    leafs[3].setBranch(node);
    leafs[3].setXMin(ox + tmp_length);
    leafs[3].setYMin(oy + tmp_length);
    leafs[3].setZMin(oz);

    leafs[3].setLength(tmp_length);
    leafs[3].setLevel(next_level);

    leafs[4].setBranch(node);
    leafs[4].setXMin(ox);
    leafs[4].setYMin(oy);
    leafs[4].setZMin(oz + tmp_length);

    leafs[4].setLength(tmp_length);
    leafs[4].setLevel(next_level);

    leafs[5].setBranch(node);
    leafs[5].setXMin(ox + tmp_length);
    leafs[5].setYMin(oy);
    leafs[5].setZMin(oz + tmp_length);

    leafs[5].setLength(tmp_length);
    leafs[5].setLevel(next_level);

    leafs[6].setBranch(node);
    leafs[6].setXMin(ox);
    leafs[6].setYMin(oy + tmp_length);
    leafs[6].setZMin(oz + tmp_length);

    leafs[6].setLength(tmp_length);
    leafs[6].setLevel(next_level);

    leafs[7].setBranch(node);
    leafs[7].setXMin(ox + tmp_length);
    leafs[7].setYMin(oy + tmp_length);
    leafs[7].setZMin(oz + tmp_length);

    leafs[7].setLength(tmp_length);
    leafs[7].setLevel(next_level);

    node->setLeafs(leafs);

    return true;
}
