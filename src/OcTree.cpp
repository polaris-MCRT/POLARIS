#include "typedefs.h"
#include "OcTree.h"
#include "MathFunctions.h"
#include "CommandParser.h"
#include <limits>

void CGridOcTree::plotNextDataVector(ofstream * file_streams, cell_oc * cell, uint level)
{
    if(cell->getChildren() == 0)
    {
        double mmx, mmy, mmz, blen = 0;
        double vvx, vvy, vvz, vlen = 0;
        double len = 0.5 * cell->getLength();

        line_counter++;
        if(line_counter % 15000 == 0)
        {
            char_counter++;
            cout << "-> Writing Gnuplot vectors: "
                    << ru[(unsigned int) char_counter % 4] << "           \r";
        }

        if(cell->getLevel() > level)
        {

            rec_counter++;

            if(rec_counter % nrOfGnuVectors == 0)
            {
                if(plt_mag)
                {
                    mmx = cell->getData(data_pos_mx);
                    mmy = cell->getData(data_pos_my);
                    mmz = cell->getData(data_pos_mz);

                    blen = sqrt(mmx * mmx + mmy * mmy + mmz * mmz);

                    if(blen > 0)
                    {
                        if(max_level > 3)
                        {
                            mmx = (len / 0.4) * (mmx / blen);
                            mmy = (len / 0.4) * (mmy / blen);
                            mmz = (len / 0.4) * (mmz / blen);
                        }
                        else
                        {
                            mmx = (len / 1.8) * (mmx / blen);
                            mmy = (len / 1.8) * (mmy / blen);
                            mmz = (len / 1.8) * (mmz / blen);
                        }
                    }
                }

                if(plt_vel)
                {
                    vvx = cell->getData(data_pos_vx);
                    vvy = cell->getData(data_pos_vy);
                    vvz = cell->getData(data_pos_vz);

                    vlen = sqrt(vvx * vvx + vvy * vvy + vvz * vvz);

                    if(vlen > 0)
                    {
                        if(max_level > 3)
                        {
                            vvx = (len / 0.4) * (vvx / vlen);
                            vvy = (len / 0.4) * (vvy / vlen);
                            vvz = (len / 0.4) * (vvz / vlen);
                        }
                        else
                        {
                            vvx = (len / 1.8) * (vvx / vlen);
                            vvy = (len / 1.8) * (vvy / vlen);
                            vvz = (len / 1.8) * (vvz / vlen);
                        }
                    }
                }

                if(blen != 0)
                    file_streams[0]
                        << float(cell->getXmin() + len - mmx) << " "
                    << float(cell->getYmin() + len - mmy) << " "
                    << float(cell->getZmin() + len - mmz) << " "
                    << float(2.0 * mmx) << " " << float(2.0 * mmy) << " "
                    << float(2.0 * mmz) << " "
                    << float(log10(blen)) << endl;

                if(vlen != 0)
                    file_streams[1]
                        << float(cell->getXmin() + len - vvx) << " "
                    << float(cell->getYmin() + len - vvy) << " "
                    << float(cell->getZmin() + len - vvz) << " "
                    << float(2.0 * vvx) << " " << float(2.0 * vvy)
                    << " " << float(2.0 * vvz) << " "
                    << float(vlen) << endl;

                //file_streams[0].flush();
                //file_streams[1].flush();
            }
        }
    }
    else
    {
        for(int i = 0; i < 8; i++)
            plotNextDataVector(file_streams, cell->getChild(i), cell->getLevel());
    }
}

void CGridOcTree::plotNextDataPoint(ofstream * file_streams, cell_oc * cell, uint level)
{
    if(cell->getChildren() == 0)
    {
        Vector3D c = getCenter(cell);
        float x = float(c.X());
        float y = float(c.Y());
        float z = float(c.Z());

        line_counter++;
        if(line_counter % 15000 == 0)
        {
            char_counter++;
            cout << "-> Writing Gnuplot point files: "
                    << ru[(unsigned int) char_counter % 4] << "           \r";
        }

        if(cell->getLevel() > level)
        {
            rec_counter++;

            if(rec_counter % nrOfGnuPoints == 0)
            {
                if(plt_gas_dens)
                    file_streams[1] << x << " " << y << " " << z << " "
                        << float(max_level + 1.5 - cell->getLevel()) / 0.9 << " "
                    << float(log10(getGasDensity(cell))) << endl;

                if(plt_dust_dens)
                    file_streams[2] << x << " " << y << " " << z << " "
                        << float(max_level + 1.5 - cell->getLevel()) / 0.9 << " "
                    << float(log10(getDustDensity(cell))) << endl;

                if(plt_gas_temp)
                    file_streams[3] << x << " " << y << " " << z << " "
                        << float(max_level + 1.5 - cell->getLevel()) / 0.9 << " "
                    << float(log10(getGasTemperature(cell))) << endl;

                if(plt_dust_temp)
                    file_streams[4] << x << " " << y << " " << z << " "
                        << float(max_level + 1.5 - cell->getLevel()) / 0.9 << " "
                    << float(log10(getDustTemperature(cell))) << endl;

                if(plt_rat)
                    file_streams[5] << x << " " << y << " " << z << " "
                        << float(max_level + 1.5 - cell->getLevel()) / 0.9 << " "
                    << getAlignedRadius(cell) << endl;

                /*if(plt_delta)
                {
                        double field = getMagField(cell).length();
                        double Td = getDustTemperature(cell);
                        double Tg = getGasTemperature(cell);
                        double dens = getGasDensity(cell);
                        double delta = CMathFunctions::calc_delta(field, Td, Tg, dens);

                        file_streams[6] << x << " " << y << " " << z << " "
                                << float(max_level + 1.5 - cell->getLevel()) / 0.9 << " "
                                << delta << endl;
                }

                if(plt_larm)
                {
                        double gas_temp = getGasTemperature(cell);
                        double field = getMagField(cell).length();
                        double Td = getDustTemperature(cell);
                        double Tg = getGasTemperature(cell);
                        double dens = getGasDensity(cell);

                        double a_limit = CMathFunctions::calc_eff_limit(field, Td, Tg, dens, 0.5);

                        file_streams[7] << x << " " << y << " " << z << " "
                                << float(max_level + 1.5 - cell->getLevel()) / 0.9 << " "
                                << a_limit << endl;
                }

                if(plt_mach)
                {
                        Vector3D vel_field = getVelocityField(cell);
                        double gas_temp = getGasTemperature(cell);
                        double mach = CMathFunctions::calc_mach(vel_field.length(), gas_temp, mu);

                        file_streams[8] << x << " " << y << " " << z << " "
                                << float(max_level + 1.5 - cell->getLevel()) / 0.9 << " "
                                << mach << endl;
                }*/
            }
        }
    }
    else
    {
        for(int i = 0; i < 8; i++)
            plotNextDataPoint(file_streams, cell->getChild(i), cell->getLevel());
    }
}

void CGridOcTree::plotNextGridCell(ofstream * grid_streams, cell_oc * cell, uint level)
{
    if(cell->getLevel() <= maxGridLines)
    {
        float len_max = float(cell->getLength());
        float len_min = float(0.5 * cell->getLength());

        float xm = float(cell->getXmin());
        float ym = float(cell->getYmin());
        float zm = float(cell->getZmin());

        line_counter++;
        if(line_counter % 15000 == 0)
        {
            char_counter++;
            cout << " -> Writing Gnuplot grid file: " << ru[(unsigned int) char_counter % 4] << "           \r";
        }

        stringstream buffer;
        buffer.str("");

        //center lines
        if(cell->getLevel() != max_level && cell->getLevel() > 1)
        {
            buffer << xm << " " << ym + len_min << " " << zm + len_min << " "
                    << len_max << " " << 0 << " " << 0 << endl;

            buffer << xm + len_min << " " << ym + len_min << " " << zm << " "
                    << 0 << " " << 0 << " " << len_max << endl;

            buffer << xm + len_min << " " << ym << " " << zm + len_min << " "
                    << 0 << " " << len_max << " " << 0 << endl;
        }

        //front face
        buffer << xm + len_min << " " << ym << " " << zm << " "
                << 0 << " " << len_max << " " << 0 << endl;

        buffer << xm << " " << ym + len_min << " " << zm << " "
                << len_max << " " << 0 << " " << 0 << endl;

        //back face
        buffer << xm + len_min << " " << ym << " " << zm + len_max << " "
                << 0 << " " << len_max << " " << 0 << endl;

        buffer << xm << " " << ym + len_min << " " << zm + len_max << " "
                << len_max << " " << 0 << " " << 0 << endl;


        //left face
        buffer << xm << " " << ym << " " << zm + len_min << " "
                << 0 << " " << len_max << " " << 0 << endl;

        buffer << xm << " " << ym + len_min << " " << zm << " "
                << 0 << " " << 0 << " " << len_max << endl;

        //right face
        buffer << xm + len_max << " " << ym << " " << zm + len_min << " "
                << 0 << " " << len_max << " " << 0 << endl;

        buffer << xm + len_max << " " << ym + len_min << " " << zm << " "
                << 0 << " " << 0 << " " << len_max << endl;

        //bottom face
        buffer << xm + len_min << " " << ym << " " << zm << " "
                << 0 << " " << 0 << " " << len_max << endl;

        buffer << xm << " " << ym << " " << zm + len_min << " "
                << len_max << " " << 0 << " " << 0 << endl;

        //top face
        buffer << xm + len_min << " " << ym + len_max << " " << zm << " "
                << 0 << " " << 0 << " " << len_max << endl;

        buffer << xm << " " << ym + len_max << " " << zm + len_min << " "
                << len_max << " " << 0 << " " << 0 << endl;

        grid_streams[0] << buffer.str();
    }

    if(cell->getChildren() != 0)
    {
        for(int i = 0; i < 8; i++)
            plotNextGridCell(grid_streams, cell->getChild(i), cell->getLevel());
    }
}

bool CGridOcTree::reduceBinrayFile(string in_filename, string out_filename, uint tr_level)
{
    parameters param;
    param.setCommand(CMD_TEMP);
    param.setPathGrid(in_filename);

    if(!loadGridFromBinrayFile(param))
        return false;

    reduceLevelOfBinrayFile(cell_oc_root, tr_level);

    if(!saveBinaryGridFile(out_filename))
        return false;

    return true;
}

bool CGridOcTree::reduceLevelOfBinrayFile(cell_oc * cell, uint tr_level)
{
    line_counter++;
    if(line_counter % 1000 == 0)
    {
        char_counter++;
        cout << "-> reducing tree: " << ru[(unsigned int) char_counter % 4] << "           \r";
    }

    if(cell->getChildren() == 0)
    {
        if(cell->getLevel() > tr_level)
            return true;

        return false;
    }
    else
    {
        bool comb = true;
        for(int i = 0; i < 8; i++)
            comb &= reduceLevelOfBinrayFile(&cell->getChildren()[i], tr_level);

        if(comb)
        {
            double dens = 0, Tg = 0, Td = 0;
            double mx = 0, my = 0, mz = 0;
            double vx = 0, vy = 0, vz = 0;

            for(int i = 0; i < 8; i++)
            {
                dens += cell->getChildren()[i].getData(data_pos_gd_list[0]);
                Tg += cell->getChildren()[i].getData(data_pos_tg);
                Td += cell->getChildren()[i].getData(data_pos_dt_list[0]);
                mx += cell->getChildren()[i].getData(data_pos_mx);
                my += cell->getChildren()[i].getData(data_pos_my);
                mz += cell->getChildren()[i].getData(data_pos_mz);

                vx += cell->getChildren()[i].getData(data_pos_vx);
                vy += cell->getChildren()[i].getData(data_pos_vy);
                vz += cell->getChildren()[i].getData(data_pos_vz);
            }

            dens /= 8;
            Tg /= 8;
            Td /= 8;
            mx /= 8;
            my /= 8;
            mz /= 8;
            vx /= 8;
            vy /= 8;
            vz /= 8;

            delete[] cell->getChildren();
            cell->setChildren(0);
            max_data = 9;
            cell->resize(max_data);
            cell->setData(data_pos_gd_list[0], dens);
            cell->setData(data_pos_tg, Tg);
            cell->setData(data_pos_dt_list[0], Td);
            cell->setData(data_pos_mx, mx);
            cell->setData(data_pos_my, my);
            cell->setData(data_pos_mz, mz);
            cell->setData(data_pos_vx, vx);
            cell->setData(data_pos_vy, vy);
            cell->setData(data_pos_vz, vz);

            if(cell->getLevel() > tr_level)
                return true;

            return false;
        }
    }

    return false;
}

bool CGridOcTree::loadGridFromBinrayFile(parameters & param, uint _data_len)
{
    double cube_length;
    int cube_pos;

    ushort tmpID, tmpOffset;
    ushort isleaf, level;

    string filename = param.getPathGrid();
    float tmp_data;

    line_counter = 0;
    char_counter = 0;

    turbulent_velocity = param.getTurbulentVelocity();

    ifstream bin_reader(filename.c_str(), ios::in | ios::binary);

    if(bin_reader.fail())
    {
        cout << "\nERROR: Cannot load octree grid file:" << endl;
        cout << filename << "\n" << endl;
        return false;
    }

    cell_oc_root = new cell_oc();
    cell_oc_pos = cell_oc_root;

    resetGridValues();

    line_counter = 1;
    char_counter = 0;
    cube_pos = -1;

    bin_reader.read((char*) &tmpID, 2);
    bin_reader.read((char*) &tmpOffset, 2);

    dataID = tmpID;
    data_offset = (uint) tmpOffset;
    data_len = _data_len + data_offset;

    data_ids.resize(data_offset);

    if(dataID != GRID_ID_OCT)
    {
        if(!createCompatibleTree())
            return false;

        double min_len;
        double max_len;

        bin_reader.read((char*) &min_len, 8);
        bin_reader.read((char*) &max_len, 8);

        bin_reader.read((char*) &min_len, 8);
        bin_reader.read((char*) &max_len, 8);

        bin_reader.read((char*) &min_len, 8);
        bin_reader.read((char*) &max_len, 8);

        cube_length = max_len - min_len;
    }
    else
    {
        for(uint i = 0; i < data_offset; i++)
        {
            ushort tmp_ids = 0;
            bin_reader.read((char*) &tmp_ids, 2);
            data_ids[i] = tmp_ids;
        }

        if(!setDataPositionsVariable())
            return false;

        bin_reader.read((char*) &cube_length, 8);
    }

    uint tmp_data_offset = validateDataPositions(param);
    if(tmp_data_offset == MAX_UINT)
        return false;

    cube_length *= conv_length_in_SI;

    if(cube_length <= 0)
    {
        cout << "\nERROR: Octree cube length must be larger than zero!" << endl;
        return false;
    }

    max_len = cube_length;

    total_volume = max_len * max_len*max_len;

    cell_oc_root->setXmin(-0.5 * max_len);
    cell_oc_root->setYmin(-0.5 * max_len);
    cell_oc_root->setZmin(-0.5 * max_len);

    cell_oc_root->setLength(max_len);
    cell_oc_root->setLevel(0);

    while(!bin_reader.eof())
    {
        line_counter++;

        if(line_counter % 5000 == 0)
        {
            char_counter++;
            cout << "-> Loading octree grid file: "
                    << ru[(unsigned int) char_counter % 4] << "           \r";
        }

        bin_reader.read((char*) &isleaf, 2);
        bin_reader.read((char*) &level, 2);

        if(isleaf == 1)
        {
            cube_pos++;

            cell_oc_pos->getChildren()[cube_pos].resize(data_len + tmp_data_offset);
            cell_oc_pos->getChildren()[cube_pos].setLevel((uchar) level);
            cell_oc_pos->getChildren()[cube_pos].setID(cube_pos);

            for(uint i = 0; i < data_offset; i++)
            {
                bin_reader.read((char*) &tmp_data, 4);
                cell_oc_pos->getChildren()[cube_pos].setData(i, tmp_data);
            }

            updateVelocity(cell_oc_pos, param);

            if(uint(cell_oc_pos->getData(data_pos_id)) < 0 ||
                uint(cell_oc_pos->getData(data_pos_id)) > param.getMaxDustComponentChoice())
            {
                cout << "\nERROR: Dust ID in grid exceeds maximum number of dust choices available! " << endl;
                return false;
            }

            //assignOpiateID(&cell_oc_pos->getChildren()[cube_pos]);
            updateDataRange(&cell_oc_pos->getChildren()[cube_pos]);

            if(cube_pos > 7)
            {
                cout << "Error in octree grid file:" << endl;
                cout << filename;
                cout << "Data set nr.: " << line_counter << " level: " << level
                        << " cell: " << cube_pos + 1;
                cout << "\nMore then 8 low level boxes!";
                return false;
            }

            if(cell_oc_pos->getLevel() + 1 != level)
            {
                cout << "Error in octree grid file:" << endl;
                cout << filename;
                cout << "Data set nr.: " << line_counter << " level: " << level
                        << " cell: " << cube_pos + 1;
                cout << "\nWrong number of low level boxes!";
                return false;
            }

            double tmp_vol = getVolume(cell_oc_pos->getChild(cube_pos));
            total_gas_mass += getGasMassDensity(cell_oc_pos->getChild(cube_pos)) * tmp_vol;
            cell_volume += tmp_vol;

            if(level > max_level)
                max_level = level;

            max_cells++;

            if(cube_pos == 7)
            {
                bool is_closed = false;

                if(cell_oc_pos->getLevel() == 0)
                    is_closed = true;
                else
                {
                    do
                    {
                        cube_pos = cell_oc_pos->getID();

                        if(cell_oc_pos->getLevel() == 0)
                        {
                            is_closed = true;
                            break;
                        }

                        cell_oc_pos = cell_oc_pos->getParent();

                    }
                    while(cube_pos == 7);
                }

                if(is_closed == true)
                    break;
            }
        }
        else
        {
            if(level > 0)
            {
                cube_pos++;

                if(cube_pos > 7)
                {
                    cout << "Error in octree grid file:" << endl;
                    cout << line_counter;
                    cout << "\nMore then 8 low level boxes!";
                    return false;
                }

                cell_oc_pos->getChildren()[cube_pos].setChildren(new cell_oc[8]);
                cell_oc_pos = &cell_oc_pos->getChildren()[cube_pos];

                createBoundingCell();

                cell_oc_pos->setLevel((uchar) level);
                cell_oc_pos->setID((uint) cube_pos);
                cube_pos = -1;
            }
            else
            {
                cell_oc_pos->setChildren(new cell_oc[8]);

                createBoundingCell();
                cell_oc_pos->setLevel((uchar) level);
            }
        }
    }

    //delete[] data;
    bin_reader.close();

    if(max_cells == 0)
    {
        cout << "\nERROR: No cells in octree grid file:" << endl;
        cout << filename << "\n" << endl;
        return false;
    }

    data_offset += tmp_data_offset;
    data_len += tmp_data_offset;

    cout << CLR_LINE;
    cout << "- Loading octree grid file             : done" << endl;

    return true;
}

bool CGridOcTree::writeGNUPlotFiles(string path, parameters & param)
{
    nrOfGnuPoints = param.getNrOfGnuPoints();
    nrOfGnuVectors = param.getNrOfGnuVectors();
    maxGridLines = param.getmaxGridLines();

    if(nrOfGnuPoints + nrOfGnuVectors == 0)
        return true;

    if(cell_oc_root == 0)
    {
        cout << "\nERROR: Cannot plot octree to:" << endl;
        cout << path;
        cout << "No tree loaded!" << endl;
        return false;
    }

    if(max_level < maxGridLines)
    {
        cout << "\nWARNING: Number of max. grid level is higher than max. tree level!" << endl;
        maxGridLines = uint(max_level);
        return false;
    }

    if(cell_oc_root->getChildren() == 0)
    {
        cout << "\nERROR: Cannot plot octree to:" << endl;
        cout << path;
        cout << "Wrong amount of tree level!" << endl;
        return false;
    }

    if(max_cells == 0)
    {
        cout << "\nERROR: Cannot plot octree to:" << endl;
        cout << path;
        cout << "Not enough tree cells available! " << endl;
        return false;
    }

    plt_gas_dens = (data_pos_gd_list.size() > 0); // 1
    plt_dust_dens = false; //param.getPlot(plIDnd) && (data_pos_dd_list.size() > 0); // 2
    plt_gas_temp = (data_pos_tg != MAX_UINT); // 3
    plt_dust_temp = (!data_pos_dt_list.empty()); // 4
    plt_rat = (data_pos_aalg != MAX_UINT); // 5
    plt_delta = (data_pos_tg != MAX_UINT) && (data_pos_mx != MAX_UINT);// && (data_pos_td != MAX_UINT); // 6
    plt_larm = (data_pos_tg != MAX_UINT) && (data_pos_mx != MAX_UINT);// && (data_pos_td != MAX_UINT); // 7
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
        cout << "\nERROR: Cannot write to:\n " << grid_filename << endl;
        return false;
    }

    if(plt_gas_dens)
    {
        point_fields[1].open(dens_gas_filename.c_str());

        if(point_fields[1].fail())
        {
            cout << "\nERROR: Cannot write to:\n " << dens_gas_filename << endl;
            return false;
        }
    }

    if(plt_dust_dens)
    {
        point_fields[2].open(dens_dust_filename.c_str());

        if(point_fields[2].fail())
        {
            cout << "\nERROR: Cannot write to:\n " << dens_dust_filename << endl;
            return false;
        }
    }

    if(plt_gas_temp)
    {
        point_fields[3].open(temp_gas_filename.c_str());

        if(point_fields[3].fail())
        {
            cout << "\nERROR: Cannot write to:\n " << temp_gas_filename
                    << endl;
            return false;
        }
    }

    if(plt_dust_temp)
    {
        point_fields[4].open(temp_dust_filename.c_str());

        if(point_fields[4].fail())
        {
            cout << "\nERROR: Cannot write to:\n " << temp_dust_filename
                    << endl;
            return false;
        }
    }

    if(plt_rat)
    {
        point_fields[5].open(rat_filename.c_str());

        if(point_fields[5].fail())
        {
            cout << "\nERROR: Cannot write to:\n " << rat_filename << endl;
            return false;
        }
    }

    if(plt_delta)
    {
        point_fields[6].open(delta_filename.c_str());

        if(point_fields[6].fail())
        {
            cout << "\nERROR: Cannot write to:\n " << delta_filename << endl;
            return false;
        }
    }

    if(plt_larm)
    {
        point_fields[7].open(larm_filename.c_str());

        if(point_fields[7].fail())
        {
            cout << "\nERROR: Cannot write to:\n " << larm_filename << endl;
            return false;
        }
    }

    if(plt_mach)
    {
        point_fields[8].open(mach_filename.c_str());

        if(point_fields[8].fail())
        {
            cout << "\nERROR: Cannot write to:\n " << mach_filename << endl;
            return false;
        }
    }

    if(plt_mag)
    {
        vec_fields[0].open(mag_filename.c_str());

        if(vec_fields[0].fail())
        {
            cout << "\nERROR: Cannot write to:\n " << mag_filename << endl;
            return false;
        }
    }

    if(plt_vel)
    {
        vec_fields[1].open(vel_filename.c_str());

        if(vec_fields[1].fail())
        {
            cout << "\nERROR: Cannot write to:\n " << vel_filename << endl;
            return false;
        }
    }

    line_counter = 0;
    char_counter = 0;

    double xm = cell_oc_root->getXmin();
    double ym = cell_oc_root->getYmin();
    double zm = cell_oc_root->getZmin();
    double len_max = cell_oc_root->getLength();
    double len_min = 0.5 * cell_oc_root->getLength();

    basic_grid_l0.str("");

    basic_grid_l0 << xm << " " << ym << " " << zm << " "
            << len_max << " " << 0 << " " << 0 << endl;
    basic_grid_l0 << xm << " " << ym << " " << zm << " "
            << 0 << " " << len_max << " " << 0 << endl;
    basic_grid_l0 << xm << " " << ym << " " << zm << " "
            << 0 << " " << 0 << " " << len_max << endl;

    basic_grid_l0 << xm + len_max << " " << ym + len_max << " " << zm + len_max << " "
            << -len_max << " " << 0 << " " << 0 << endl;
    basic_grid_l0 << xm + len_max << " " << ym + len_max << " " << zm + len_max << " "
            << 0 << " " << -len_max << " " << 0 << endl;
    basic_grid_l0 << xm + len_max << " " << ym + len_max << " " << zm + len_max << " "
            << 0 << " " << 0 << " " << -len_max << endl;


    basic_grid_l0 << xm + len_max << " " << ym + len_max << " " << zm << " "
            << -len_max << " " << 0 << " " << 0 << endl;
    basic_grid_l0 << xm + len_max << " " << ym + len_max << " " << zm << " "
            << 0 << " " << -len_max << " " << 0 << endl;


    basic_grid_l0 << xm + len_max << " " << ym + len_max << " " << zm << " "
            << -len_max << " " << 0 << " " << 0 << endl;
    basic_grid_l0 << xm + len_max << " " << ym + len_max << " " << zm << " "
            << 0 << " " << -len_max << " " << 0 << endl;


    basic_grid_l0 << xm << " " << ym + len_max << " " << zm + len_max << " "
            << 0 << " " << -len_max << " " << 0 << endl;
    basic_grid_l0 << xm << " " << ym + len_max << " " << zm + len_max << " "
            << 0 << " " << 0 << " " << -len_max << endl;


    basic_grid_l0 << xm + len_max << " " << ym << " " << zm + len_max << " "
            << -len_max << " " << 0 << " " << 0 << endl;
    basic_grid_l0 << xm + len_max << " " << ym << " " << zm + len_max << " "
            << 0 << " " << 0 << " " << -len_max << endl;












    //center lines

    basic_grid_l1 << xm << " " << ym + len_min << " " << zm + len_min << " "
            << len_max << " " << 0 << " " << 0 << endl;

    basic_grid_l1 << xm + len_min << " " << ym + len_min << " " << zm << " "
            << 0 << " " << 0 << " " << len_max << endl;

    basic_grid_l1 << xm + len_min << " " << ym << " " << zm + len_min << " "
            << 0 << " " << len_max << " " << 0 << endl;


    //front face
    basic_grid_l1 << xm + len_min << " " << ym << " " << zm << " "
            << 0 << " " << len_max << " " << 0 << endl;

    basic_grid_l1 << xm << " " << ym + len_min << " " << zm << " "
            << len_max << " " << 0 << " " << 0 << endl;

    //back face
    basic_grid_l1 << xm + len_min << " " << ym << " " << zm + len_max << " "
            << 0 << " " << len_max << " " << 0 << endl;

    basic_grid_l1 << xm << " " << ym + len_min << " " << zm + len_max << " "
            << len_max << " " << 0 << " " << 0 << endl;


    //left face
    basic_grid_l1 << xm << " " << ym << " " << zm + len_min << " "
            << 0 << " " << len_max << " " << 0 << endl;

    basic_grid_l1 << xm << " " << ym + len_min << " " << zm << " "
            << 0 << " " << 0 << " " << len_max << endl;

    //right face
    basic_grid_l1 << xm + len_max << " " << ym << " " << zm + len_min << " "
            << 0 << " " << len_max << " " << 0 << endl;

    basic_grid_l1 << xm + len_max << " " << ym + len_min << " " << zm << " "
            << 0 << " " << 0 << " " << len_max << endl;

    //bottom face
    basic_grid_l1 << xm + len_min << " " << ym << " " << zm << " "
            << 0 << " " << 0 << " " << len_max << endl;

    basic_grid_l1 << xm << " " << ym << " " << zm + len_min << " "
            << len_max << " " << 0 << " " << 0 << endl;

    //top face
    basic_grid_l1 << xm + len_min << " " << ym + len_max << " " << zm << " "
            << 0 << " " << 0 << " " << len_max << endl;

    basic_grid_l1 << xm << " " << ym + len_max << " " << zm + len_min << " "
            << len_max << " " << 0 << " " << 0 << endl;

    cout << CLR_LINE;

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

    point_header << "set xrange[" << 1.01 * cell_oc_root->getXmin() << ":"
            << 1.01 * cell_oc_root->getXmax() << "]" << endl;
    point_header << "set yrange[" << 1.01 * cell_oc_root->getYmin() << ":"
            << 1.01 * cell_oc_root->getYmax() << "]" << endl;
    point_header << "set zrange[" << 1.01 * cell_oc_root->getZmin() << ":"
            << 1.01 * cell_oc_root->getZmax() << "]" << endl;

    point_header << "set style arrow 1 nohead ls 1 lw 1 lc rgb 0x0000cc" << endl;
    point_header << "set style arrow 2 nohead ls 1 lw 1 lc rgb 0xaa0000" << endl;
    point_header << "set style line 1 pt 1 ps variable lt palette" << endl;

    point_header << "set grid" << endl;
    point_header << "set nokey" << endl;

    //0 octre grid
    point_fields[0] << point_header.str();
    point_fields[0] << "set title \'3D octree grid geometry\' font \'Arial,12\'" << endl;
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

    if(min_gas_dens / max_gas_dens > 0.9)
        min_gas_dens = 0.9 * max_gas_dens;

    point_fields[1] << "set cbrange[" << log10(min_gas_dens) << ":"
            << log10(max_gas_dens) << "]" << endl;
    point_fields[1] << "set format cb \'%.02g\'" << endl;

    point_fields[1] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;


    //2 dust density
    point_fields[2] << point_header.str();
    point_fields[2] << "set title \'3D dust number density distribution (min: " << min_dust_dens
            << "[m^-3]; max: " << max_dust_dens << "[m^-3])\' font \'Arial,12\'" << endl;
    point_fields[2] << "set cblabel \'gas density[m^-3]\'" << endl;
    point_fields[2] << "set palette defined (0 0.5 0 0, 1 0 0 1, 2 0 1 1)" << endl;

    if(min_dust_dens == 0 && max_dust_dens == 0)
    {
        min_dust_dens = 0.1;
        max_dust_dens = 1;
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

    vec_header << "set xlabel \'x\'" << endl;
    vec_header << "set ylabel \'y\'" << endl;
    vec_header << "set zlabel \'z\'" << endl;

    vec_header << "set xrange[" << 1.01 * cell_oc_root->getXmin() << ":"
            << 1.01 * cell_oc_root->getXmax() << "]" << endl;
    vec_header << "set yrange[" << 1.01 * cell_oc_root->getYmin() << ":"
            << 1.01 * cell_oc_root->getYmax() << "]" << endl;
    vec_header << "set zrange[" << 1.01 * cell_oc_root->getZmin() << ":"
            << 1.01 * cell_oc_root->getZmax() << "]" << endl;

    vec_header << "set style arrow 1 nohead ls 1 lw 1 lc rgb 0x0000cc" << endl;
    vec_header << "set style arrow 2 nohead ls 1 lw 1 lc rgb 0xbb00cc" << endl;
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

    line_counter = 0;
    rec_counter = 0;
    plotNextDataPoint(point_fields, cell_oc_root, 0);

    line_counter = 0;
    rec_counter = 0;
    plotNextDataVector(vec_fields, cell_oc_root, 0);

    line_counter = 0;
    rec_counter = 0;
    plotNextGridCell(&point_fields[0], cell_oc_root, 0);

    for(uint pos = 0; pos < 7; pos++)
    {
        point_fields[pos]
                << "\n e \n" << basic_grid_l1.str()
                << "\n e \n" << basic_grid_l0.str() << "\n e" << endl;
    }

    for(uint pos = 0; pos < 2; pos++)
    {
        vec_fields[pos] << "\n e \n" << basic_grid_l1.str()
                << "\n e \n" << basic_grid_l0.str() << "\n e" << endl;
    }

    for(uint pos = 0; pos < 8; pos++)
        point_fields[pos].close();

    for(uint pos = 0; pos < 2; pos++)
        vec_fields[pos].close();

    cout << CLR_LINE;
    cout << "- Writing of Gnuplot files             : done" << endl;
    return true;
}

bool CGridOcTree::saveBinaryGridFile(string filename, ushort id, ushort data_size)
{
    if(cell_oc_root == 0)
    {
        cout << "\nERROR: Cannot save octree grid file to:" << endl;
        cout << filename;
        cout << "No tree loaded!" << endl;
        return false;
    }

    if(cell_oc_root->getChildren() == 0)
    {
        cout << "\nERROR: Cannot save octree grid file to:" << endl;
        cout << filename;
        cout << "Octree has to be refined at least to level 1!" << endl;
        return false;
    }

    if(max_cells == 0)
    {
        cout << "\nERROR: Cannot save octree grid file to:" << endl;
        cout << filename;
        cout << "No cells available! " << endl;
        return false;
    }

    double x_min, y_min, z_min, cube_length;
    line_counter = 0;
    char_counter = 0;

    ofstream bin_writer(filename.c_str(), ios::out | ios::binary);

    if(bin_writer.fail())
    {
        cout << "\nERROR: Cannot open octree grid file:" << endl;
        cout << filename;
        return false;
    }

    x_min = y_min = z_min = cell_oc_root->getXmin();
    cube_length = cell_oc_root->getLength();

    bin_writer.write((char*) &id, 2);
    bin_writer.write((char*) &data_size, 2);

    if(dataID == GRID_ID_OCT)
    {
        for(uint i = 0; i < data_offset; i++)
        {
            ushort tmp_ids = data_ids[i];
            bin_writer.write((char*) &tmp_ids, 2);
        }
    }
    else
    {
        cout << "\nERROR: Cannot save octree grid file to:" << endl;
        cout << filename;
        cout << "An octree grid requires an ID of " << GRID_ID_OCT << "!" << endl;
        return false;
    }

    bin_writer.write((char*) &cube_length, 8);

    nextBinaryDataCell(bin_writer, cell_oc_root, data_size);

    bin_writer.close();
    cout << CLR_LINE;
    cout << "- Writing octree grid file             : done" << endl;

    return true;
}

void CGridOcTree::nextBinaryDataCell(ofstream & file_stream, cell_oc * cell, uint data_size)
{
    ushort isleaf, level;
    float data;

    if(cell->getChildren() == 0)
    {
        line_counter++;
        if(line_counter % 15000 == 0)
        {
            char_counter++;
            cout << "-> Writing binary octree grid file: "
                    << ru[(unsigned int) char_counter % 4] << "           \r" << flush;
        }

        isleaf = (ushort) 1;
        level = (ushort) cell->getLevel();

        file_stream.write((char*) &isleaf, 2);
        file_stream.write((char*) &level, 2);

        for(uint pos = 0; pos < data_size; pos++)
        {
            data = (float) cell->getData(pos);
            file_stream.write((char*) &data, 4);
        }
    }
    else
    {
        isleaf = (ushort) 0;
        level = (ushort) cell->getLevel();

        file_stream.write((char*) &isleaf, 2);
        file_stream.write((char*) &level, 2);

        for(int i = 0; i < 8; i++)
            nextBinaryDataCell(file_stream, &cell->getChildren()[i], data_size);
    }
}

void CGridOcTree::printParameters()
{
    if(max_cells == 0)
        cout << "\nERROR: No octree grid parameters available! " << endl;
    else
    {
        ulong tmp_cells = ulong(pow(double(8), double(max_level)));
        cout << CLR_LINE;
        cout << "OcTree parameters (ID: " << getDataID() << "; data len.: " << getDataOffset() << "; level: " << max_level << ")" << endl;
        cout << SEP_LINE;
        cout << "- Number of OcTree cells        : " << max_cells << "(data), " << tmp_cells << " (max)" << endl;

        printPhysicalParameters();
        cout << SEP_LINE;
    }
}

bool CGridOcTree::createArtificialGrid(string path)
{
    resetGridValues();

    cell_oc_root = new cell_oc();
    cell_oc_pos = cell_oc_root;
    double cube_length = 8 * con_pc;
    max_level = 5;

    dataID = GRID_ID_OCT;
    data_offset = 9;
    max_data = 9;

    data_pos_dt_list.push_back(0);
    data_pos_dt_list.push_back(1);
    data_pos_tg = 2;
    data_pos_mx = 3;
    data_pos_my = 4;
    data_pos_mz = 5;
    data_pos_vx = 6;
    data_pos_vy = 7;
    data_pos_vz = 8;

    data_ids.resize(data_offset);

    data_ids[0] = GRIDgas_dens;
    data_ids[1] = GRIDdust_temp;
    data_ids[2] = GRIDgas_temp;
    data_ids[3] = GRIDmx;
    data_ids[4] = GRIDmy;
    data_ids[5] = GRIDmz;
    data_ids[6] = GRIDvx;
    data_ids[7] = GRIDvy;
    data_ids[8] = GRIDvz;

    cell_oc_root->setXmin(-cube_length / 2.0);
    cell_oc_root->setYmin(-cube_length / 2.0);
    cell_oc_root->setZmin(-cube_length / 2.0);
    cell_oc_root->setLength(cube_length);

    max_cells = ulong(pow(8.0, double(max_level)));

    line_counter = 1;
    char_counter = 0;

    createNextLevel(cell_oc_root);

    cout << "min: " << min_gas_dens << "  max_dens: " << max_gas_dens << endl;
    cout << "Creating artificial tree                    : done" << endl;
    cout << "Max cells: " << max_cells << endl;
    return true;
}

void CGridOcTree::createNextLevel(cell_oc * cell)
{
    Vector3D p = getCenter(cell);
    double r = p.length();
    uint tmp_level = uint(max_level);

    if(cell->getLevel() >= tmp_level)
    {
        double Tg = 10;
        double Td = 20;

        double dens = 1e37 / (r * r + con_pc);

        Vector3D mag(0, 0, 10);
        Vector3D vel(0, 0, 10);

        line_counter++;
        if(line_counter % 15000 == 0)
        {
            char_counter++;
            cout << " - Creating artificial tree: "
                    << ru[(unsigned int) char_counter % 4] << "           \r";
        }

        cell->resize(max_data);
        cell->setData(data_pos_gd_list[0], dens);
        cell->setData(data_pos_tg, Tg);
        cell->setData(data_pos_dt_list[0], Td);

        cell->setData(data_pos_mx, mag.X());
        cell->setData(data_pos_my, mag.Y());
        cell->setData(data_pos_mz, mag.Z());

        cell->setData(data_pos_vx, (vel.X()));
        cell->setData(data_pos_vy, (vel.Y()));
        cell->setData(data_pos_vz, (vel.Z()));
    }
    else
    {
        cell->setChildren(new cell_oc[8]);
        bool found = false;
        cell_oc_pos = cell;
        createBoundingCell();

        double avg_dens = 0;
        double avg_tg = 0;
        double avg_td = 0;
        double avg_mx = 0;
        double avg_my = 0;
        double avg_mz = 0;
        double avg_vx = 0;
        double avg_vy = 0;
        double avg_vz = 0;

        for(uint i = 0; i < 8; i++)
        {
            cell->getChildren()[i].setLevel(cell->getLevel() + 1);
            createNextLevel(&cell->getChildren()[i]);
        }

        bool nl = true;
        double factor = 1e-58;

        for(uint i = 0; i < 8; i++)
        {
            if(cell->getChildren()[i].getChildren() != 0)
                nl = false;
        }

        if(nl)
        {
            nl = true;
            for(unsigned int i = 0; i < 8; i++)
            {
                if(cell->getChildren()[i].getData(data_pos_gd_list[0]) > factor)
                {
                    nl = false;
                    break;
                }

                avg_dens += cell->getChildren()[i].getData(data_pos_gd_list[0]);
                avg_tg += cell->getChildren()[i].getData(data_pos_tg);
                avg_td += cell->getChildren()[i].getData(data_pos_dt_list[0]);
                avg_mx += cell->getChildren()[i].getData(data_pos_mx);
                avg_my += cell->getChildren()[i].getData(data_pos_my);
                avg_mz += cell->getChildren()[i].getData(data_pos_mz);

                avg_vx += cell->getChildren()[i].getData(data_pos_vx);
                avg_vy += cell->getChildren()[i].getData(data_pos_vy);
                avg_vz += cell->getChildren()[i].getData(data_pos_vz);
            }

            if(nl)
            {
                cell->resize(max_data);
                cell->setData(data_pos_gd_list[0], avg_dens / 8.0);
                cell->setData(data_pos_tg, avg_tg / 8.0);
                cell->setData(data_pos_dt_list[0], avg_td / 8.0);
                cell->setData(data_pos_mx, (avg_mx / 8.0));
                cell->setData(data_pos_my, (avg_my / 8.0));
                cell->setData(data_pos_mz, (avg_mz / 8.0));

                cell->setData(data_pos_vx, (avg_vx / 8.0));
                cell->setData(data_pos_vy, (avg_vy / 8.0));
                cell->setData(data_pos_vz, (avg_vz / 8.0));

                delete[] cell->getChildren();
                cell->setChildren(0);
                max_cells -= 7;
            }
        }
    }
}

void CGridOcTree::createBoundingCell()
{
    double ox, oy, oz, length;

    ox = cell_oc_pos->getXmin();
    oy = cell_oc_pos->getYmin();
    oz = cell_oc_pos->getZmin();

    length = 0.5 * cell_oc_pos->getLength();

    if(length < min_len)
        min_len = length;

    cell_oc_pos->getChildren()[0].setParent(cell_oc_pos);
    cell_oc_pos->getChildren()[0].setXmin(ox);
    cell_oc_pos->getChildren()[0].setYmin(oy);
    cell_oc_pos->getChildren()[0].setZmin(oz);

    cell_oc_pos->getChildren()[0].setLength(length);


    cell_oc_pos->getChildren()[1].setParent(cell_oc_pos);
    cell_oc_pos->getChildren()[1].setXmin(ox + length);
    cell_oc_pos->getChildren()[1].setYmin(oy);
    cell_oc_pos->getChildren()[1].setZmin(oz);

    cell_oc_pos->getChildren()[1].setLength(length);

    cell_oc_pos->getChildren()[2].setParent(cell_oc_pos);
    cell_oc_pos->getChildren()[2].setXmin(ox);
    cell_oc_pos->getChildren()[2].setYmin(oy + length);
    cell_oc_pos->getChildren()[2].setZmin(oz);

    cell_oc_pos->getChildren()[2].setLength(length);

    cell_oc_pos->getChildren()[3].setParent(cell_oc_pos);
    cell_oc_pos->getChildren()[3].setXmin(ox + length);
    cell_oc_pos->getChildren()[3].setYmin(oy + length);
    cell_oc_pos->getChildren()[3].setZmin(oz);

    cell_oc_pos->getChildren()[3].setLength(length);

    cell_oc_pos->getChildren()[4].setParent(cell_oc_pos);
    cell_oc_pos->getChildren()[4].setXmin(ox);
    cell_oc_pos->getChildren()[4].setYmin(oy);
    cell_oc_pos->getChildren()[4].setZmin(oz + length);

    cell_oc_pos->getChildren()[4].setLength(length);

    cell_oc_pos->getChildren()[5].setParent(cell_oc_pos);
    cell_oc_pos->getChildren()[5].setXmin(ox + length);
    cell_oc_pos->getChildren()[5].setYmin(oy);
    cell_oc_pos->getChildren()[5].setZmin(oz + length);

    cell_oc_pos->getChildren()[5].setLength(length);

    cell_oc_pos->getChildren()[6].setParent(cell_oc_pos);
    cell_oc_pos->getChildren()[6].setXmin(ox);
    cell_oc_pos->getChildren()[6].setYmin(oy + length);
    cell_oc_pos->getChildren()[6].setZmin(oz + length);

    cell_oc_pos->getChildren()[6].setLength(length);

    cell_oc_pos->getChildren()[7].setParent(cell_oc_pos);
    cell_oc_pos->getChildren()[7].setXmin(ox + length);
    cell_oc_pos->getChildren()[7].setYmin(oy + length);
    cell_oc_pos->getChildren()[7].setZmin(oz + length);

    cell_oc_pos->getChildren()[7].setLength(length);
}

bool CGridOcTree::goToNextCellBorder(photon_package * pp)
{
    Vector3D tmp_pos_xyz, tmp_dir_xyz;
    Vector3D v_n, v_a, v_S, v_ds;
    double d_ls, min_dl, num, den;

    double loc_x_min, loc_x_max, loc_y_min, loc_y_max, loc_z_min, loc_z_max;
    double loc_dx, loc_dy, loc_dz;
    bool hit = false;

    double path_length = 0;
    cell_oc * tmp_cell_pos = (cell_oc*) pp->getPositionCell();

    tmp_pos_xyz = pp->getPosition();
    tmp_dir_xyz = pp->getDirection();

    loc_x_min = tmp_cell_pos->getXmin();
    loc_y_min = tmp_cell_pos->getYmin();
    loc_z_min = tmp_cell_pos->getZmin();

    loc_x_max = tmp_cell_pos->getXmax();
    loc_y_max = tmp_cell_pos->getYmax();
    loc_z_max = tmp_cell_pos->getZmax();

    loc_dx = loc_x_max - loc_x_min;
    loc_dy = loc_y_max - loc_y_min;
    loc_dz = loc_z_max - loc_z_min;

    min_dl = 1E200;
    d_ls = -1;

    for(int i_side = 1; i_side <= 6; i_side++)
    {
        v_n = 0;
        v_a = 0;

        switch(i_side)
        {
            case 1:
                v_n.setZ(-loc_dx * loc_dy);
                v_a.setZ(loc_z_min);
                break;
            case 2:
                v_n.setZ(loc_dx * loc_dy);
                v_a.setZ(loc_z_max);
                break;
            case 3:
                v_n.setY(-loc_dx * loc_dz);
                v_a.setY(loc_y_min);
                break;
            case 4:
                v_n.setY(loc_dx * loc_dz);
                v_a.setY(loc_y_max);
                break;
            case 5:
                v_n.setX(-loc_dy * loc_dz);
                v_a.setX(loc_x_min);
                break;
            case 6:
                v_n.setX(loc_dy * loc_dz);
                v_a.setX(loc_x_max);
                break;
        }

        num = v_n * (tmp_pos_xyz - v_a);
        den = v_n * (tmp_dir_xyz);

        if(den != 0)
        {
            d_ls = -num / den;

            if(d_ls >= 0 && d_ls < min_dl)
            {
                min_dl = d_ls;
                hit = true;
            }
        }
    }

    if(hit == true)
    {
        path_length = min_dl;
        path_length *= 1.00001;

        if(path_length < 1e-3 * loc_dx)
        {
            path_length = 1e-3 * loc_dx;
        }

        pp->setPosition(tmp_pos_xyz + tmp_dir_xyz * path_length);
        pp->setTmpPathLength(path_length);
    }

    return hit;
}

bool CGridOcTree::updateShortestDistance(photon_package * pp)
{
    Vector3D tmp_pos;
    double min_dist, tmp_dist[6];

    double loc_x_min, loc_x_max, loc_y_min, loc_y_max, loc_z_min, loc_z_max;
    bool found = false;

    cell_oc * tmp_cell_pos = (cell_oc*) pp->getPositionCell();

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

    pp->setShortestDistance(min_dist);
    return found;
}

/*bool CGridOcTree::updateShortestDistance(photon_package * pp)
 {
 Vector3D tmp_pos_xyz, tmp_dir_xyz;
 Vector3D v_n, v_a, v_S, v_ds;
 double len_a, len_n, min_dist, tmp_dist;

 double loc_x_min, loc_x_max, loc_y_min, loc_y_max, loc_z_min, loc_z_max;
 double loc_dx, loc_dy, loc_dz;
 bool hit = false;

 double path_length = 0;
 cell_oc * tmp_cell_pos = (cell_oc*)pp->getPositionCell();

 tmp_pos_xyz = pp->getPosition();
 tmp_dir_xyz = pp->getDirection();

 loc_x_min = tmp_cell_pos->getXmin();
 loc_y_min = tmp_cell_pos->getYmin();
 loc_z_min = tmp_cell_pos->getZmin();

 loc_x_max = tmp_cell_pos->getXmax();
 loc_y_max = tmp_cell_pos->getYmax();
 loc_z_max = tmp_cell_pos->getZmax();

 loc_dx = loc_x_max - loc_x_min;
 loc_dy = loc_y_max - loc_y_min;
 loc_dz = loc_z_max - loc_z_min;

 min_dist = 1E50;

 for(int i_side = 1; i_side <= 6; i_side++)
 {
 v_n = 0;
 v_a = 0;

 switch(i_side)
 {
 case 1:
 v_n.setZ(-loc_dx*loc_dy);
 v_a.setZ(loc_z_min);
 break;
 case 2:
 v_n.setZ(loc_dx*loc_dy);
 v_a.setZ(loc_z_max);
 break;
 case 3:
 v_n.setY(-loc_dx*loc_dz);
 v_a.setY(loc_y_min);
 break;
 case 4:
 v_n.setY(loc_dx*loc_dz);
 v_a.setY(loc_y_max);
 break;
 case 5:
 v_n.setX(-loc_dy*loc_dz);
 v_a.setX(loc_x_min);
 break;
 case 6:
 v_n.setX(loc_dy*loc_dz);
 v_a.setX(loc_x_max);
 break;
 }

 len_n = v_n.length();
 //den = v_n * tmp_dir_xyz;
 len_a = (v_n*v_a);

 //num = den - len_a;
 tmp_dist = abs((v_n*tmp_pos_xyz - len_a) / len_n);

 if(min_dist>tmp_dist)
 min_dist = tmp_dist;
 }

 pp->setShortestDistance(min_dist);
 return hit;
 }*/

bool CGridOcTree::nextLowLevelCell()
{
    uint cube_pos;

    if(cell_oc_pos == cell_oc_root)
    {
        cube_pos = 0;
    }

    if(cell_oc_pos->getChildren() != 0)
    {
        while(cell_oc_pos->getChildren() != 0)
            cell_oc_pos = &cell_oc_pos->getChildren()[0];

        return true;
    }

    cube_pos = cell_oc_pos->getID() + 1;

    if(cube_pos > 7)
    {
        do
        {
            cell_oc_pos = cell_oc_pos->getParent();
            if(cell_oc_pos == 0)
                return false;

            if(cell_oc_pos == cell_oc_root)
                return false;

            cube_pos = cell_oc_pos->getID();
        }
        while(cube_pos == 7);

        if(cell_oc_pos == 0)
            return false;

        if(cell_oc_pos == cell_oc_root)
            return false;

        cube_pos = cell_oc_pos->getID() + 1;
    }

#pragma warning(suppress: 6011)
    cell_oc_pos = &cell_oc_pos->getParent()->getChildren()[cube_pos];

    if(cell_oc_pos->getChildren() != 0)
        return nextLowLevelCell();
    else
        return true;

    return false;
}

bool CGridOcTree::nextLowLevelCell(cell_basic * cell)
{
    uint cube_pos;
    cell_oc * extern_cell = (cell_oc*) cell;

    if(extern_cell == 0)
    {
        extern_cell = cell_oc_root;
    }

    if(extern_cell == cell_oc_root)
    {
        cube_pos = 0;
    }

    if(extern_cell->getChildren() != 0)
    {
        while(extern_cell->getChildren() != 0)
            extern_cell = &extern_cell->getChildren()[0];

        return true;
    }

    cube_pos = extern_cell->getID() + 1;

    if(cube_pos > 7)
    {
        do
        {
            extern_cell = extern_cell->getParent();
            if(extern_cell == 0)
                return false;

            if(extern_cell == cell_oc_root)
                return false;

            cube_pos = extern_cell->getID();
        }
        while(cube_pos == 7);

        if(extern_cell == 0)
            return false;

        if(extern_cell == cell_oc_root)
            return false;

        cube_pos = extern_cell->getID() + 1;
    }

#pragma warning(suppress: 6011)
    extern_cell = &extern_cell->getParent()->getChildren()[cube_pos];

    if(extern_cell->getChildren() != 0)
        return nextLowLevelCell(cell);
    else
        return true;

    return false;
}

bool CGridOcTree::findStartingPoint(photon_package * pp)
{
    Vector3D pos = pp->getPosition();
    if(isInside(pos))
        return true;

    positionPhotonInGrid(pp);
    uint try_counter = 0;
    bool res = false;

    while(goToNextCellBorder(pp))
    {
        try_counter++;
        if(try_counter > 6)
            return false;

        if(findMatchingCell(pp))
        {
            res = true;
            break;
        }
    }

    return res;
}

void CGridOcTree::clear(cell_oc * cell)
{
    if(cell->getChildren() != 0)
    {
        line_counter++;
        if(line_counter % 15000 == 0)
        {
            char_counter++;
            cout << " -> Final cleanup: "
                    << ru[(unsigned int) char_counter % 4] << "              \r";
        }
        
//#pragma omp parallel for schedule(dynamic)
        for(unsigned int i = 0; i < 8; i++)
            clear(&cell->getChildren()[i]);

        delete[] cell->getChildren();
        cell->setChildren(0);
    }

    CGridOcTree();
}

bool CGridOcTree::initiateTreeFromFile(uint _nx, uint _max_level, double _fa,
        double _length, string str_dens, string str_temp,
        string str_magx, string str_magy, string str_magz)
{
    ifstream reader_dens(str_dens.c_str());
    ifstream reader_temp(str_temp.c_str());
    ifstream reader_magx(str_magx.c_str());
    ifstream reader_magy(str_magy.c_str());
    ifstream reader_magz(str_magz.c_str());


    double dens, temp, magx, magy, magz;

    unsigned int per_counter = 0;

    max_value = 0;

    nx = _nx;
    factor = _fa;


    dataID = GRID_ID_OCT;
    data_offset = 9;
    max_data = 9;

    data_pos_dt_list.push_back(0);
    data_pos_dt_list.push_back(1);
    data_pos_tg = 2;
    data_pos_mx = 3;
    data_pos_my = 4;
    data_pos_mz = 5;
    data_pos_vx = 6;
    data_pos_vy = 7;
    data_pos_vz = 8;

    data_ids.resize(data_offset);

    data_ids[0] = GRIDgas_dens;
    data_ids[1] = GRIDdust_temp;
    data_ids[2] = GRIDgas_temp;
    data_ids[3] = GRIDmx;
    data_ids[4] = GRIDmy;
    data_ids[5] = GRIDmz;
    data_ids[6] = GRIDvx;
    data_ids[7] = GRIDvy;
    data_ids[8] = GRIDvz;


    if(reader_dens.fail())
    {
        cout << "todo :ERROR" << endl;
        return false;
    }

    datdens.resize(nx);
    dattemp.resize(nx);
    datmx.resize(nx);
    datmy.resize(nx);
    datmz.resize(nx);

    max_gas_dens = -1e30;
    min_gas_dens = 1e30;

    max_delta = -1e30;
    min_delta = 1e30;

    max_mag = -1e30;
    min_mag = 1e30;

    max_gas_temp = -1e30;
    min_gas_temp = 1e30;

    max_dust_temp = -1e30;
    min_dust_temp = 1e30;

    max_delta = -1e30;
    min_delta = 1e30;

    max_level = _max_level;

    min_len = 1e30;
    max_len = _length;
    min_len = max_len / (pow(2, max_level));

    max_cells = (uint) pow(double(8.0), double(max_level));

    for(uint i = 0; i < nx; ++i)
    {
        datdens[i].resize(nx);
        dattemp[i].resize(nx);
        datmx[i].resize(nx);
        datmy[i].resize(nx);
        datmz[i].resize(nx);

        for(uint j = 0; j < nx; ++j)
        {
            datdens[i][j].resize(nx);
            dattemp[i][j].resize(nx);
            datmx[i][j].resize(nx);
            datmy[i][j].resize(nx);
            datmz[i][j].resize(nx);

            per_counter++;

            if(j % 10 == 0)
                cout << " -> Reading input data: "
                    << 100.0 * float(per_counter) / float(nx * nx)
                << " [%]                \r";

            for(uint k = 0; k < nx; k++)
            {
                //string line;
                //getline(reader_dens, line);


                reader_dens >> dens;
                reader_temp >> temp;
                reader_magx>> magx;
                reader_magy >> magy;
                reader_magz >> magz;

                //temp = 0;
                //magx = 0;
                //magy = 0;
                //magz = 1;

                datdens[i][j][k] = dens;
                dattemp[i][j][k] = temp;
                datmx[i][j][k] = magx;
                datmy[i][j][k] = magy;
                datmz[i][j][k] = magz;
            }
        }
    }

    double xyz_min = -0.5 * max_len;


    f_min = 1e30;
    f_max = -1e30;

    treelevel_counter = 0;
    tagged_cells = 0;
    cell_oc_root = new cell_oc();
    createTree(cell_oc_root, 0, 0, 0, max_len, 0);
    //(cell_oc * parent, double _x_min, double _y_min, double _z_min, double _length, uint _level)
    reader_dens.close();

    cell_oc_root->setLength(4.7305E+19);
    cout << CLR_LINE;
    cout << "min , max 1   " << f_min << "\t" << f_max << endl;
    return true;
}/**/

bool CGridOcTree::createTree(cell_oc * parent, double _x_min, double _y_min, double _z_min, double _length, uint _level)
{
    //double dx, dy, dz,
    double field, delta;
    float max_cells;

    //dx = (x_max - x_min) / 2.0;
    //dy = (y_max - y_min) / 2.0;
    //dz = (z_max - z_min) / 2.0;

    //double qx = x_max - dx - 0.5;
    //double qy = y_max - dy - 0.5;
    //double qz = z_max - dy - 0.5;

    //double len = sqrt(qx * qx + qy * qy + qz * qz);

    //if(len>0.9)
    //cout << len << endl;

    parent->setXmin(_x_min);
    parent->setYmin(_y_min);
    parent->setZmin(_z_min);
    parent->setLevel(_level);
    parent->setLength(_length);

    double scale_factor = 1.0; //double(nx) / (_length);

    double px, py, pz;
    double dens, gas_temp, dust_temp, mx, my, mz, vx, vy, vz;
    int p_x, p_y, p_z;

    //calcMeanValue(parent);
    max_cells = (float) pow(8.0, max_level);
    //
    if(parent->getLevel() == max_level)
    {
        treelevel_counter++;
        Vector3D center = getCenter((cell_basic*) parent);

        px = center.X();
        py = center.Y();
        pz = center.Z();

        if(treelevel_counter % 10000 == 0)
        {
            cout << "Creating tree: "
                    << float(100.0 * treelevel_counter) / max_cells
                    << " [%]                          \r";
        }

        p_x = int(px * scale_factor);
        p_y = int(py * scale_factor);
        p_z = int(pz * scale_factor);

        //cout << px << "\t" << p_x << endl;

        if(p_x > (int) nx - 1)
            p_x = (int) nx - 1;
        if(p_y > (int) nx - 1)
            p_y = (int) nx - 1;
        if(p_z > (int) nx - 1)
            p_z = (int) nx - 1;

        if(p_x < 0)
            p_x = 0;
        if(p_y < 0)
            p_y = 0;
        if(p_z < 0)
            p_z = 0;

        //cout << "1" << endl;

        dens = 1 * datdens[p_x][p_y][p_z];
        gas_temp = dattemp[p_x][p_y][p_z];
        dust_temp = 10;
        vx = mx = datmx[p_x][p_y][p_z];
        vy = my = datmy[p_x][p_y][p_z];
        vz = mz = datmz[p_x][p_y][p_z];

        //if(dens<1e-25) dens=1e-25;
        px = center.X() - 64;
        py = center.Y() - 64;
        pz = center.Z() - 64;

        double lx = sqrt(px * px + py * py + pz * pz) / 128.0 / 0.85926;

        //if(lx>1)
        {
            //double f=1000/1.1*(lx-1.9)+1;
            double mx = pow((lx - 1), 6);
            double px = pow((lx + 1), 6);
            double f = (exp(-0.9 * px) + exp(-0.9 * mx) - 0.81319) / 0.186806; //= exp(-9.21*xx*xx);

            if(lx > 0.6)
            {
                double max = 100;
                double d = (max - 1) / (1 - 0.6)*(lx - 0.6) + 1;
                f /= d;

                if(d < f_min)
                    f_min = d;
                if(d > f_max)
                    f_max = d;
            }


            dens *= f;
            mx *= f;
            my *= f;
            mz *= f;


        }



        //if(dens<1e-24) dens = 1e-24;

        field = sqrt(mx * mx + my * my + mz * mz);

        if(field < 1e-12)
        {
            mx *= 1e-12 / field;
            my *= 1e-12 / field;
            mz *= 1e-12 / field;
        }

        Vector3D vel;
        dust_temp = 2.5 * (30 + log10(dens));

        max_data = 9;
        parent->resize(max_data);
        parent->setData(data_pos_gd_list[0], dens);
        parent->setData(data_pos_tg, gas_temp);
        parent->setData(data_pos_dt_list[0], dust_temp);
        parent->setData(data_pos_mx, mx);
        parent->setData(data_pos_my, my);
        parent->setData(data_pos_mz, mz);
        parent->setData(data_pos_vx, vx);
        parent->setData(data_pos_vy, vy);
        parent->setData(data_pos_vz, vz);

        vel.set(1.5 * my, -mz, 2 * mx);
        vel.normalize();
        vel *= 1000 * 60000 / abs(log10(dens));

        if(vel.length() < 1)
            vel.set(1, 1, 1);

        parent->setData(data_pos_vx, (vel.X()));
        parent->setData(data_pos_vy, (vel.Y()));
        parent->setData(data_pos_vz, (vel.Z()));

        meanBdir += Vector3D(mx, my, mz);
        field = sqrt(mx * mx + my * my + mz * mz);

        delta = CMathFunctions::calc_delta(field, dust_temp, gas_temp, dens);

        if(dens > max_gas_dens)
            max_gas_dens = dens;
        if(dens < min_gas_dens)
            min_gas_dens = dens;

        if(field > max_mag)
            max_mag = field;
        if(field < min_mag)
            min_mag = field;

        if(gas_temp > max_gas_temp)
            max_gas_temp = gas_temp;
        if(gas_temp < min_gas_temp)
            min_gas_temp = gas_temp;

        if(dust_temp > max_dust_temp)
            max_dust_temp = dust_temp;
        if(dust_temp < min_dust_temp)
            min_dust_temp = dust_temp;

        if(delta > max_delta)
            max_delta = delta;
        if(delta < min_delta)
            min_delta = delta;

        return true;
    }

    parent->setChildren(new cell_oc[8]);
    uint level = parent->getLevel() + 1;
    double length = 0.5 * parent->getLength();

    createTree(parent->getChild(0), _x_min, _y_min, _z_min, length, level);
    createTree(parent->getChild(1), _x_min + length, _y_min, _z_min, length, level);
    createTree(parent->getChild(2), _x_min, _y_min + length, _z_min, length, level);
    createTree(parent->getChild(3), _x_min + length, _y_min + length, _z_min, length, level);

    createTree(parent->getChild(4), _x_min, _y_min, _z_min + length, length, level);
    createTree(parent->getChild(5), _x_min + length, _y_min, _z_min + length, length, level);
    createTree(parent->getChild(6), _x_min, _y_min + length, _z_min + length, length, level);
    createTree(parent->getChild(7), _x_min + length, _y_min + length, _z_min + length, length, level);

    for(int i = 0; i < 8; i++)
        parent->getChild(i)->setParent(parent);

    if(parent->getChild(0)->getChildren() == 0 && parent->getChild(1)->getChildren() == 0
            && parent->getChild(2)->getChildren() == 0
            && parent->getChild(3)->getChildren() == 0
            && parent->getChild(4)->getChildren() == 0
            && parent->getChild(5)->getChildren() == 0
            && parent->getChild(6)->getChildren() == 0
            && parent->getChild(7)->getChildren() == 0)
    {
        double limit = 2e-23;

        //if(len<0.51)
        //   limit=0;

        if(parent->getChild(0)->getData(data_pos_gd_list[0]) < limit
                && parent->getChild(1)->getData(data_pos_gd_list[0]) < limit
                && parent->getChild(2)->getData(data_pos_gd_list[0]) < limit
                && parent->getChild(3)->getData(data_pos_gd_list[0]) < limit
                && parent->getChild(4)->getData(data_pos_gd_list[0]) < limit
                && parent->getChild(5)->getData(data_pos_gd_list[0]) < limit
                && parent->getChild(6)->getData(data_pos_gd_list[0]) < limit
                && parent->getChild(7)->getData(data_pos_gd_list[0]) < limit)
        {

            dens = gas_temp = dust_temp = mx = my = mz = vx = vy = vz = 0;

            for(unsigned int i = 0; i < 8; i++)
            {
                dens += parent->getChild(i)->getData(data_pos_gd_list[0]);
                gas_temp += parent->getChild(i)->getData(data_pos_tg);
                dust_temp += parent->getChild(i)->getData(data_pos_dt_list[0]);
                mx += parent->getChild(i)->getData(data_pos_mx);
                my += parent->getChild(i)->getData(data_pos_my);
                mz += parent->getChild(i)->getData(data_pos_mz);

                vx += parent->getChild(i)->getData(data_pos_vx);
                vy += parent->getChild(i)->getData(data_pos_vy);
                vz += parent->getChild(i)->getData(data_pos_vz);
            }

            dens /= 8.0;
            gas_temp /= 8.0;
            dust_temp /= 8.0;
            mx /= 8.0;
            my /= 8.0;
            mz /= 8.0;
            vx /= 8.0;
            vy /= 8.0;
            vz /= 8.0;

            /*if(level < 5)
             {
             dens /= 5;
             mx /= 5;
             my /= 5;
             mz /= 5;
             }*/

            /*if(level < 5)
             {
             dens /=10;
             mx /= 10;
             my /= 10;
             mz /= 10;
             dust_temp/=2;

             if(dust_temp<3)
             dust_temp=3;

             }*/

            Vector3D vel(vx, vy, vz);
            if(vel.length() < 1)
                vel.normalize();

            //if(dust_temp < 1) dust_temp = 0.1;

            //if(dens<1e-20) dens=1e-20;
            parent->resize(max_data);
            parent->setData(data_pos_gd_list[0], dens);
            parent->setData(data_pos_tg, gas_temp);
            parent->setData(data_pos_dt_list[0], dust_temp);
            parent->setData(data_pos_mx, (mx));
            parent->setData(data_pos_my, (my));
            parent->setData(data_pos_mz, (mz));

            parent->setData(data_pos_vx, (vx));
            parent->setData(data_pos_vy, (vy));
            parent->setData(data_pos_vz, (vz));

            meanBdir += Vector3D(mx, my, mz);
            field = sqrt(mx * mx + my * my + mz * mz);

            delta = CMathFunctions::calc_delta(field, dust_temp, gas_temp,
                    dens);

            if(dens > max_gas_dens)
                max_gas_dens = dens;
            if(dens < min_gas_dens)
                min_gas_dens = dens;

            if(field > max_mag)
                max_mag = field;
            if(field < min_mag)
                min_mag = field;

            if(gas_temp > max_gas_temp)
                max_gas_temp = gas_temp;
            if(gas_temp < min_gas_temp)
                min_gas_temp = gas_temp;

            if(dust_temp > max_dust_temp)
                max_dust_temp = dust_temp;
            if(dust_temp < min_dust_temp)
                min_dust_temp = dust_temp;

            if(delta > max_delta)
                max_delta = delta;
            if(delta < min_delta)
                min_delta = delta;

            delete[] parent->getChildren();
            parent->setChildren(0);
            max_cells -= 7;
        }

    }

    return true;
}
