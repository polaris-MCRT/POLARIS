/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "GridOcTree.hpp"
#include "CommandParser.hpp"
#include "MathFunctions.hpp"
#include "Parameters.hpp"
#include "Typedefs.hpp"

bool CGridOcTree::reduceBinaryFile(string in_filename, string out_filename, uint tr_level)
{
    parameters param;
    param.setCommand(CMD_TEMP);
    param.setPathGrid(in_filename);

    if(!loadGridFromBinaryFile(param))
        return false;

    createCellList();

    printParameters();

    cell_oc_root=&cell_oc_root->getChildren()[6];


    ulong max_cells = getMaxDataCells();

    cout << CLR_LINE;


    #pragma omp parallel for schedule(dynamic)
    for(long c_i = 0; c_i < long(max_cells); c_i++)
    {
        cell_basic * cell = getCellFromIndex(c_i);
        double dens0 = getGasDensity(*cell, 0);
        double dens1 = getGasDensity(*cell, 1);


        // if(c_i%5000==0)
        //     cout << "-> " << float(100*c_i)/float(max_cells) << "                       \r" << flush;

        double dens = dens0+dens1;


        double length = ((cell_oc *)cell)->getLength();

        double cx = ((cell_oc *)cell)->getXmin()+0.5*length;
        double cy = ((cell_oc *)cell)->getYmin()+0.5*length;
        double cz = ((cell_oc *)cell)->getZmin()+0.5*length;

        vector<Vector3D> vlist;

        length=0.51*length;

        Vector3D center=Vector3D(cx,cy,cz);

        vlist.push_back(Vector3D(cx+length,cy,cz));
        vlist.push_back(Vector3D(cx-length,cy,cz));

        vlist.push_back(Vector3D(cx,cy+length,cz));
        vlist.push_back(Vector3D(cx,cy-length,cz));

        vlist.push_back(Vector3D(cx,cy,cz+length));
        vlist.push_back(Vector3D(cx,cy,cz-length));

        length=1.732*length;

        vlist.push_back(Vector3D(cx+length,cy+length,cz+length));
        vlist.push_back(Vector3D(cx+length,cy-length,cz+length));
        vlist.push_back(Vector3D(cx-length,cy-length,cz+length));
        vlist.push_back(Vector3D(cx-length,cy+length,cz+length));

        vlist.push_back(Vector3D(cx+length,cy+length,cz-length));
        vlist.push_back(Vector3D(cx+length,cy-length,cz-length));
        vlist.push_back(Vector3D(cx-length,cy-length,cz-length));
        vlist.push_back(Vector3D(cx-length,cy+length,cz-length));

        photon_package pp;

        pp.setPosition(center);

        if(positionPhotonInGrid(&pp))
        {
            //cell_basic * cell = pp.getPositionCell()
            double tg= 0.8*getGasTemperature(pp);

            for(uint g=0;g<vlist.size();g++)
            {
                pp.setPosition(vlist[g]);

                if(positionPhotonInGrid(&pp))
                {
                    tg+=0.8/14.0*getGasTemperature(pp);
                }
            }

            setGasTemperature(cell, tg);
        }


        if(dens*254098.7886>1e13)
        {
            setGasDensity(cell, 0, 0.0*dens);
            setGasDensity(cell, 1, 1.0*dens);
        }
        else
        {
            setGasDensity(cell, 0, 1.0*dens);
            setGasDensity(cell, 1, 0.0*dens);
        }
    }

    //reduceLevelOfBinaryFile(cell_oc_root, tr_level);

    if(!saveBinaryGridFile(out_filename))
        return false;

    return true;
}

bool CGridOcTree::reduceLevelOfBinaryFile(cell_oc * cell, uint tr_level)
{
    line_counter++;
    if(line_counter % 1000 == 0)
    {
        char_counter++;
        cout << "-> reducing tree: " << ru[(unsigned int)char_counter % 4] << "           \r";
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
            comb &= reduceLevelOfBinaryFile(&cell->getChildren()[i], tr_level);

        if(comb)
        {
            double * tmp_data=new double[data_len];

            for(uint j=0;j<data_len;j++)
            {
                tmp_data[j]=0;
            }

            for(int i = 0; i < 8; i++)
            {
                for(uint j =0;j<data_len;j++)
                {
                    tmp_data[j]+=cell->getChildren()[i].getData(j)/8.0;
                }
            }

            delete[] cell->getChildren();
            cell->setChildren(0);
            cell->resize(data_len);

            for(uint j=0;j<data_len;j++)
            {
                cell->setData(j, tmp_data[j]); ;
            }

            if(cell->getLevel() > tr_level)
                return true;

            return false;
        }
    }

    return false;
}

bool CGridOcTree::loadGridFromBinaryFile(parameters & param, uint _data_len)
{
    double cube_length;
    int cube_pos;

    ushort tmpID, tmpOffset;
    ushort isleaf, level;

    string filename = param.getPathGrid();
    float tmp_data;

    line_counter = 0;
    char_counter = 0;

    ifstream bin_reader(filename.c_str(), ios::in | ios::binary);

    if(bin_reader.fail())
    {
        cout << ERROR_LINE << "Cannot load octree grid file:" << endl;
        cout << filename << "\n" << endl;
        return false;
    }

    cell_oc_root = new cell_oc();
    cell_oc_pos = cell_oc_root;

//    resetGridValues();

    turbulent_velocity = param.getTurbulentVelocity();

    line_counter = 1;
    char_counter = 0;
    cube_pos = -1;
    float last_percentage = 0;

    bin_reader.read((char *)&tmpID, 2);
    bin_reader.read((char *)&tmpOffset, 2);

    dataID = tmpID;
    data_offset = (uint)tmpOffset;
    data_len = _data_len + data_offset;

    data_ids.resize(data_offset);

    if(dataID != GRID_ID_OCT)
    {
        if(!createCompatibleTree())
            return false;

        double min_len;
        double max_len;

        bin_reader.read((char *)&min_len, 8);
        bin_reader.read((char *)&max_len, 8);

        bin_reader.read((char *)&min_len, 8);
        bin_reader.read((char *)&max_len, 8);

        bin_reader.read((char *)&min_len, 8);
        bin_reader.read((char *)&max_len, 8);

        cube_length = max_len - min_len;
    }
    else
    {
        for(uint i = 0; i < data_offset; i++)
        {
            ushort tmp_ids = 0;
            bin_reader.read((char *)&tmp_ids, 2);
            data_ids[i] = tmp_ids;
        }

        if(!setDataPositionsVariable())
            return false;

        bin_reader.read((char *)&cube_length, 8);
    }

    uint tmp_data_offset = validateDataPositions(param);
    if(tmp_data_offset == MAX_UINT)
        return false;

    cube_length *= conv_length_in_SI;

    if(cube_length <= 0)
    {
        cout << ERROR_LINE << "Octree cube length must be larger than zero!" << endl;
        return false;
    }

    max_len = cube_length;

    total_volume = max_len * max_len * max_len;

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
             cout << "-> Loading octree grid file: " << ru[(unsigned int)char_counter % 4] << "           \r";
        }

        /*// Calculate percentage of total progress per source
        float percentage = 100.0 * double(line_counter) / double(max_cells);

        // Show only new percentage number if it changed
        if((percentage - last_percentage) > PERCENTAGE_STEP)
        {
            char_counter++;
            cout << "-> Loading octree grid file: " << percentage << " [%]      \r" << flush;
            last_percentage = percentage;
        }*/

        bin_reader.read((char *)&isleaf, 2);
        bin_reader.read((char *)&level, 2);

        if(isleaf == 1)
        {
            cube_pos++;

            cell_oc_pos->getChildren()[cube_pos].resize(data_len + tmp_data_offset);
            cell_oc_pos->getChildren()[cube_pos].setLevel((uchar)level);
            cell_oc_pos->getChildren()[cube_pos].setID(cube_pos);

            for(uint i = 0; i < data_offset; i++)
            {
                bin_reader.read((char *)&tmp_data, 4);
                cell_oc_pos->getChildren()[cube_pos].setData(i, tmp_data);
            }

            updateVelocity(cell_oc_pos, param);

            if(uint(cell_oc_pos->getData(data_pos_id)) < 0 ||
               uint(cell_oc_pos->getData(data_pos_id)) > param.getMaxDustComponentChoice())
            {
                cout << ERROR_LINE << "Dust ID in grid exceeds maximum number of dust choices "
                        "available! "
                     << endl;
                return false;
            }

            // assignOpiateID(&cell_oc_pos->getChildren()[cube_pos]);
            updateDataRange(&cell_oc_pos->getChildren()[cube_pos]);

            if(cube_pos > 7)
            {
                cout << "Error in octree grid file:" << endl;
                cout << filename;
                cout << "Data set nr.: " << line_counter << " level: " << level << " cell: " << cube_pos + 1;
                cout << "\nMore then 8 low level boxes!";
                return false;
            }

            if(cell_oc_pos->getLevel() + 1 != level)
            {
                cout << "Error in octree grid file:" << endl;
                cout << filename;
                cout << "Data set nr.: " << line_counter << " level: " << level << " cell: " << cube_pos + 1;
                cout << "\nWrong number of low level boxes!";
                return false;
            }

            double tmp_vol = getVolume(*cell_oc_pos->getChild(cube_pos));
            total_gas_mass += getGasMassDensity(*cell_oc_pos->getChild(cube_pos)) * tmp_vol;
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

                    } while(cube_pos == 7);
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

                cell_oc_pos->setLevel((uchar)level);
                cell_oc_pos->setID((uint)cube_pos);
                cube_pos = -1;
            }
            else
            {
                cell_oc_pos->setChildren(new cell_oc[8]);

                createBoundingCell();
                cell_oc_pos->setLevel((uchar)level);
            }
        }
    }

    // delete[] data;
    bin_reader.close();

    if(max_cells == 0)
    {
        cout << ERROR_LINE << "No cells in octree grid file:" << endl;
        cout << filename << "\n" << endl;
        return false;
    }

    data_offset += tmp_data_offset;
    data_len += tmp_data_offset;

    // cout << CLR_LINE;
    // cout << "- Loading octree grid file             : done" << endl;

    return true;
}


bool CGridOcTree::saveBinaryGridFile(string filename, ushort id, ushort data_size)
{
    if(cell_oc_root == 0)
    {
        cout << ERROR_LINE << "Cannot save octree grid file to:" << endl;
        cout << filename;
        cout << "No tree loaded!" << endl;
        return false;
    }

    if(cell_oc_root->getChildren() == 0)
    {
        cout << ERROR_LINE << "Cannot save octree grid file to:" << endl;
        cout << filename;
        cout << "Octree has to be refined at least to level 1!" << endl;
        return false;
    }

    if(max_cells == 0)
    {
        cout << ERROR_LINE << "Cannot save octree grid file to:" << endl;
        cout << filename;
        cout << "No cells available! " << endl;
        return false;
    }

    line_counter = 0;
    char_counter = 0;

    ofstream bin_writer(filename.c_str(), ios::out | ios::binary);

    if(bin_writer.fail())
    {
        cout << ERROR_LINE << "Cannot open octree grid file:" << endl;
        cout << filename;
        return false;
    }

    // double x_min, y_min, z_min;
    // x_min = y_min = z_min = cell_oc_root->getXmin();
    double cube_length = cell_oc_root->getLength();

    bin_writer.write((char *)&id, 2);
    bin_writer.write((char *)&data_size, 2);

    if(dataID == GRID_ID_OCT)
    {
        for(uint i = 0; i < data_offset; i++)
        {
            ushort tmp_ids = data_ids[i];
            bin_writer.write((char *)&tmp_ids, 2);
        }
    }
    else
    {
        cout << ERROR_LINE << "Cannot save octree grid file to:" << endl;
        cout << filename;
        cout << "An octree grid requires an ID of " << GRID_ID_OCT << "!" << endl;
        return false;
    }

    bin_writer.write((char *)&cube_length, 8);

    nextBinaryDataCell(bin_writer, cell_oc_root, data_size);

    bin_writer.close();
    cout << CLR_LINE;
    cout << "- Writing octree grid file      : done" << endl;

    return true;
}

void CGridOcTree::nextBinaryDataCell(ofstream & file_stream, cell_oc * cell, uint data_size)
{
    ushort isleaf, level;
    float data;

    if(cell->getChildren() == 0)
    {
        line_counter++;
        // if(line_counter % 15000 == 0)
        // {
        //     char_counter++;
        //     cout << "-> Writing binary octree grid file: " << ru[(unsigned int)char_counter % 4]
        //          << "           \r" << flush;
        // }

        isleaf = (ushort)1;
        level = (ushort)cell->getLevel();//-1;

        file_stream.write((char *)&isleaf, 2);
        file_stream.write((char *)&level, 2);

        for(uint pos = 0; pos < data_size; pos++)
        {
            data = (float)cell->getData(pos);
            file_stream.write((char *)&data, 4);
        }
    }
    else
    {
        isleaf = (ushort)0;
        level = (ushort)cell->getLevel();//-1;

        file_stream.write((char *)&isleaf, 2);
        file_stream.write((char *)&level, 2);

        for(int i = 0; i < 8; i++)
            nextBinaryDataCell(file_stream, &cell->getChildren()[i], data_size);
    }
}

void CGridOcTree::printParameters()
{
    if(max_cells == 0)
        cout << ERROR_LINE << "No octree grid parameters available! " << endl;
    else
    {
        ulong tmp_cells = ulong(pow(double(8), double(max_level)));
        cout << CLR_LINE;
        cout << "OcTree parameters (ID: " << getDataID() << "; data len.: " << getDataOffset()
             << "; level: " << max_level << ")" << endl;
        cout << SEP_LINE;
        cout << "- Number of OcTree cells        : " << max_cells << "(data), " << tmp_cells << " (max)"
             << endl;

        printPhysicalParameters();
        cout << SEP_LINE;
    }
}

bool CGridOcTree::createArtificialGrid(string path)
{
//    resetGridValues();

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
    Vector3D p = getCenter(*cell);
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
            cout << " - Creating artificial tree: " << ru[(unsigned int)char_counter % 4] << "           \r";
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
    cell_oc * tmp_cell = (cell_oc *)pp->getPositionCell();

    bool hit = false;
    double path_length = 1e300;

    Vector3D pos = pp->getPosition();
    Vector3D dir = pp->getDirection();

    double loc_x_min = tmp_cell->getXmin();
    double loc_y_min = tmp_cell->getYmin();
    double loc_z_min = tmp_cell->getZmin();

    double loc_x_max = tmp_cell->getXmax();
    double loc_y_max = tmp_cell->getYmax();
    double loc_z_max = tmp_cell->getZmax();

    Vector3D v_n, v_a;
    double num, den, length;
    // length_eps is the minimum step width to ensure that
    // the photon 1) moves and 2) enters the cell !numerically!
    double length_eps_1, length_eps_2;

    for(uint i_side = 0; i_side < 6; i_side++)
    {
        // v_n is normal vector on the cell border pointing
        // towards next cell
        v_n = 0;
        // v_a is a point on the cell border
        v_a = 0;
        switch(i_side)
        {
            case 0:
                v_n.setZ(-1);
                v_a.setZ(loc_z_min);
                break;
            case 1:
                v_n.setZ(1);
                v_a.setZ(loc_z_max);
                break;
            case 2:
                v_n.setY(-1);
                v_a.setY(loc_y_min);
                break;
            case 3:
                v_n.setY(1);
                v_a.setY(loc_y_max);
                break;
            case 4:
                v_n.setX(-1);
                v_a.setX(loc_x_min);
                break;
            case 5:
                v_n.setX(1);
                v_a.setX(loc_x_max);
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
            // if num > 0 -> border is behind the photon
            num = v_n * (pos - v_a);

            // distance num to border is enlarged to ensure that at least one
            // component of the photon position changes parallel to v_n after step
            // sign(num) is necessary to ensure that abs(num) gets larger
            length_eps_1 = abs(pos * v_n) * MIN_LEN_STEP * EPS_DOUBLE;
            num += Vector3D::sign(num) * length_eps_1;

            length = -num / den;
            
            //position is exactly at the cell wall
            if(abs(num)<=EPS_DOUBLE)
            {
                const double eps_len = (loc_z_max - loc_z_min) * EPS_DOUBLE;
                
                if (eps_len < path_length)
                {
                    hit = true;
                    path_length = eps_len;
                    continue;
                }
            }

            if(length > 0 && length < path_length)
            {
                hit = true;
                length_eps_2 = abs( (pos + dir * length) * v_n ) / den * MIN_LEN_STEP*EPS_DOUBLE;
                path_length = length + length_eps_2;
            }
        }
    }

    if(!hit)
    {
        cout << ERROR_LINE << "Wrong cell border!                                   " << endl;
        return false;
    }

    pp->setPosition(pos + dir * path_length);

    if(pos == pp->getPosition())
    {
        cout << ERROR_LINE << "Could not transfer photon to the next cell border!   " << endl;
        return false;
    }

    pp->setTmpPathLength(path_length);

    return true;
}

bool CGridOcTree::updateShortestDistance(photon_package * pp)
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
        } while(cube_pos == 7);

        if(cell_oc_pos == 0)
            return false;

        if(cell_oc_pos == cell_oc_root)
            return false;

        cube_pos = cell_oc_pos->getID() + 1;
    }

    #pragma warning(suppress : 6011)
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
    cell_oc * extern_cell = (cell_oc *)cell;

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
        } while(cube_pos == 7);

        if(extern_cell == 0)
            return false;

        if(extern_cell == cell_oc_root)
            return false;

        cube_pos = extern_cell->getID() + 1;
    }

    #pragma warning(suppress : 6011)
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
        return positionPhotonInGrid(pp);

    bool hit = false;
    double path_length = 1e300;

    Vector3D dir = pp->getDirection();

    double loc_x_min = cell_oc_root->getXmin();
    double loc_y_min = cell_oc_root->getYmin();
    double loc_z_min = cell_oc_root->getZmin();

    double loc_x_max = cell_oc_root->getXmax();
    double loc_y_max = cell_oc_root->getYmax();
    double loc_z_max = cell_oc_root->getZmax();

    Vector3D v_n, v_a;
    double num, den, length;
    // length_eps is the minimum step width to ensure that
    // the photon 1) moves and 2) enters the cell !numerically!
    double length_eps_1, length_eps_2;

    for(uint i_side = 0; i_side < 6; i_side++)
    {
        // v_n is normal vector on the cell border pointing inside that cell
        // signs are switched compared to goToNextCellBorder
        // because we are outside the cell
        v_n = 0;
        // v_a is a point on the cell border
        v_a = 0;
        switch(i_side)
        {
            case 0:
                v_n.setZ(1);
                v_a.setZ(loc_z_min);
                break;
            case 1:
                v_n.setZ(-1);
                v_a.setZ(loc_z_max);
                break;
            case 2:
                v_n.setY(1);
                v_a.setY(loc_y_min);
                break;
            case 3:
                v_n.setY(-1);
                v_a.setY(loc_y_max);
                break;
            case 4:
                v_n.setX(1);
                v_a.setX(loc_x_min);
                break;
            case 5:
                v_n.setX(-1);
                v_a.setX(loc_x_max);
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

            if(length > 0 && length < path_length)
            {
                if(isInside(pos + dir * length))
                {
                    hit = true;
                    length_eps_2 = abs( (pos + dir * length) * v_n ) / den * MIN_LEN_STEP*EPS_DOUBLE;
                    path_length = length + length_eps_2;
                }
            }
        }
    }

    if(!hit)
    {
        cout << ERROR_LINE << "Wrong cell border!                                   " << endl;
        return false;
    }

    pp->setPosition(pos + dir * path_length);

    return positionPhotonInGrid(pp);
}

void CGridOcTree::clear(cell_oc * cell)
{
    if(cell->getChildren() != 0)
    {
        line_counter++;
        if(line_counter % 15000 == 0)
        {
            char_counter++;
            cout << " -> Final cleanup: " << ru[(unsigned int)char_counter % 4] << "              \r";
        }

        //#pragma omp parallel for schedule(dynamic)
        for(int i = 0; i < 8; i++)
            clear(&cell->getChildren()[i]);

        delete[] cell->getChildren();
        cell->setChildren(0);
    }

    CGridOcTree();
}

bool CGridOcTree::initiateTreeFromFile(uint _nx,
                                       uint _max_level,
                                       double _fa,
                                       double _length,
                                       string str_dens,
                                       string str_temp,
                                       string str_magx,
                                       string str_magy,
                                       string str_magz)
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

    max_mag = -1e30;
    min_mag = 1e30;

    max_gas_temp = -1e30;
    min_gas_temp = 1e30;

    max_dust_temp = -1e30;
    min_dust_temp = 1e30;


    max_level = _max_level;

    min_len = 1e30;
    max_len = _length;
    min_len = max_len / (pow(2, max_level));

    max_cells = (uint)pow(double(8.0), double(max_level));

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

            // if(j % 10 == 0)
            //     cout << " -> Reading input data: " << 100.0 * float(per_counter) / float(nx * nx)
            //          << " [%]                \r";

            for(uint k = 0; k < nx; k++)
            {
                // string line;
                // getline(reader_dens, line);

                reader_dens >> dens;
                reader_temp >> temp;
                reader_magx >> magx;
                reader_magy >> magy;
                reader_magz >> magz;

                // temp = 0;
                // magx = 0;
                // magy = 0;
                // magz = 1;

                datdens[i][j][k] = dens;
                dattemp[i][j][k] = temp;
                datmx[i][j][k] = magx;
                datmy[i][j][k] = magy;
                datmz[i][j][k] = magz;
            }
        }
    }

    f_min = 1e30;
    f_max = -1e30;

    treelevel_counter = 0;
    tagged_cells = 0;
    cell_oc_root = new cell_oc();
    createTree(cell_oc_root, 0, 0, 0, max_len, 0);
    //(cell_oc * parent, double _x_min, double _y_min, double _z_min, double _length, uint
    //_level)
    reader_dens.close();

    cell_oc_root->setLength(4.7305E+19);
    cout << CLR_LINE;
    cout << "min , max 1   " << f_min << "\t" << f_max << endl;
    return true;
} /**/

bool CGridOcTree::createTree(cell_oc * parent,
                             double _x_min,
                             double _y_min,
                             double _z_min,
                             double _length,
                             uint _level)
{
    // double dx, dy, dz,
    double field, delta;
    float max_cells;

    // dx = (x_max - x_min) / 2.0;
    // dy = (y_max - y_min) / 2.0;
    // dz = (z_max - z_min) / 2.0;

    // double qx = x_max - dx - 0.5;
    // double qy = y_max - dy - 0.5;
    // double qz = z_max - dy - 0.5;

    // double len = sqrt(qx * qx + qy * qy + qz * qz);

    // if(len>0.9)
    // cout << len << endl;

    parent->setXmin(_x_min);
    parent->setYmin(_y_min);
    parent->setZmin(_z_min);
    parent->setLevel(_level);
    parent->setLength(_length);

    double scale_factor = 1.0; // double(nx) / (_length);

    double px, py, pz;
    double dens, gas_temp, dust_temp, mx, my, mz, vx, vy, vz;
    int p_x, p_y, p_z;

    // calcMeanValue(parent);
    max_cells = (float)pow(8.0, max_level);
    //
    if(parent->getLevel() == max_level)
    {
        treelevel_counter++;
        Vector3D center = getCenter((const cell_basic &)parent);

        px = center.X();
        py = center.Y();
        pz = center.Z();

        if(treelevel_counter % 10000 == 0)
        {
            cout << "Creating tree: " << float(100.0 * treelevel_counter) / max_cells
                 << " [%]                          \r";
        }

        p_x = int(px * scale_factor);
        p_y = int(py * scale_factor);
        p_z = int(pz * scale_factor);

        // cout << px << "\t" << p_x << endl;

        if(p_x > (int)nx - 1)
            p_x = (int)nx - 1;
        if(p_y > (int)nx - 1)
            p_y = (int)nx - 1;
        if(p_z > (int)nx - 1)
            p_z = (int)nx - 1;

        if(p_x < 0)
            p_x = 0;
        if(p_y < 0)
            p_y = 0;
        if(p_z < 0)
            p_z = 0;

        // cout << "1" << endl;

        dens = 1 * datdens[p_x][p_y][p_z];
        gas_temp = dattemp[p_x][p_y][p_z];
        dust_temp = 10;
        vx = mx = datmx[p_x][p_y][p_z];
        vy = my = datmy[p_x][p_y][p_z];
        vz = mz = datmz[p_x][p_y][p_z];

        // if(dens<1e-25) dens=1e-25;
        px = center.X() - 64;
        py = center.Y() - 64;
        pz = center.Z() - 64;

        double lx = sqrt(px * px + py * py + pz * pz) / 128.0 / 0.85926;

        // if(lx>1)
        {
            // double f=1000/1.1*(lx-1.9)+1;
            double mx = pow((lx - 1), 6);
            double px = pow((lx + 1), 6);
            double f = (exp(-0.9 * px) + exp(-0.9 * mx) - 0.81319) / 0.186806; //= exp(-9.21*xx*xx);

            if(lx > 0.6)
            {
                double max = 100;
                double d = (max - 1) / (1 - 0.6) * (lx - 0.6) + 1;
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

        // if(dens<1e-24) dens = 1e-24;

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
    createTree(parent->getChild(7),
               _x_min + length,
               _y_min + length,
               _z_min + length,
               length,
               level);

    for(int i = 0; i < 8; i++)
        parent->getChild(i)->setParent(parent);

    if(parent->getChild(0)->getChildren() == 0 && parent->getChild(1)->getChildren() == 0 &&
       parent->getChild(2)->getChildren() == 0 && parent->getChild(3)->getChildren() == 0 &&
       parent->getChild(4)->getChildren() == 0 && parent->getChild(5)->getChildren() == 0 &&
       parent->getChild(6)->getChildren() == 0 && parent->getChild(7)->getChildren() == 0)
    {
        double limit = 2e-23;

        // if(len<0.51)
        //   limit=0;

        if(parent->getChild(0)->getData(data_pos_gd_list[0]) < limit &&
           parent->getChild(1)->getData(data_pos_gd_list[0]) < limit &&
           parent->getChild(2)->getData(data_pos_gd_list[0]) < limit &&
           parent->getChild(3)->getData(data_pos_gd_list[0]) < limit &&
           parent->getChild(4)->getData(data_pos_gd_list[0]) < limit &&
           parent->getChild(5)->getData(data_pos_gd_list[0]) < limit &&
           parent->getChild(6)->getData(data_pos_gd_list[0]) < limit &&
           parent->getChild(7)->getData(data_pos_gd_list[0]) < limit)
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

            // if(dust_temp < 1) dust_temp = 0.1;

            // if(dens<1e-20) dens=1e-20;
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

            delete[] parent->getChildren();
            parent->setChildren(0);
            max_cells -= 7;
        }
    }

    return true;
}

void CGridOcTree::goToRoot()
{
    cell_oc_pos = cell_oc_root;
}

Vector3D CGridOcTree::getCenter(const cell_basic & cell) const
{
    Vector3D center;
    const cell_oc * tmp_cell = (const cell_oc *)&cell;

    center.setX(tmp_cell->getXmin() + 0.5 * tmp_cell->getLength());
    center.setY(tmp_cell->getYmin() + 0.5 * tmp_cell->getLength());
    center.setZ(tmp_cell->getZmin() + 0.5 * tmp_cell->getLength());

    return center;
}

Vector3D CGridOcTree::getMidplaneCenter(cell_basic * cell)
{
    Vector3D center;
    cell_oc * tmp_cell = (cell_oc *)cell;

    center.setX(tmp_cell->getXmin() + 0.5 * tmp_cell->getLength());
    center.setY(tmp_cell->getYmin() + 0.5 * tmp_cell->getLength());
    center.setZ(0);

    return center;
}

bool CGridOcTree::createCellList()
{
    if(max_cells == 0)
    {
        cout << ERROR_LINE << "OcTree grid contains no cells!" << endl;
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
        // if(pos_counter >= 2097152)
        //     return true;

        cell_list[pos_counter] = (cell_basic *)cell_oc_pos;
        cell_oc_pos->setUniqueID(pos_counter);

        pos_counter++;
    }

    // cout << CLR_LINE;
    // cout << "- Creating cell list                   : done          \n" << flush;
    return true;
}

bool CGridOcTree::findMatchingCell(photon_package * pp)
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

bool CGridOcTree::next(photon_package * pp)
{
    if(!findMatchingCell(pp))
        return false;

    if(!goToNextCellBorder(pp))
        return false;

    return true;
}

void CGridOcTree::getLengths(uint bins, double & step_xy, double & off_xy)
{
    step_xy = (cell_oc_root->getLength()) / double(bins);

    off_xy = step_xy / 2.0;
}

double CGridOcTree::getVolume(const cell_basic & cell) const
{
    const cell_oc * cell_pos = (const cell_oc *)&cell;

    double volume = cell_pos->getLength();
    volume = volume * volume * volume;

    return volume;
}

bool CGridOcTree::positionPhotonInGrid(photon_package * pp)
{
    pp->setPositionCell(cell_oc_root);
    return findMatchingCell(pp);
}

const cell_oc * CGridOcTree::getTopLevelCell() const
{
    return cell_oc_root;
}

const cell_oc * CGridOcTree::getCurrentCell() const
{
    return cell_oc_pos;
}

bool CGridOcTree::saveBinaryGridFile(string filename)
{
    return saveBinaryGridFile(filename, GRID_ID_OCT, data_offset);
}

bool CGridOcTree::loadGridFromBinaryFile(parameters & param)
{
    return loadGridFromBinaryFile(param, 0);
}

void CGridOcTree::clear()
{
    line_counter = 0;
    char_counter = 0;
    clear(cell_oc_root);
    cell_oc_root = 0;
    cell_oc_pos = 0;
    cout << "Final cleanup                                : done" << endl;
}

void CGridOcTree::goNextLevelDown(photon_package * pp)
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

void CGridOcTree::goNextLevelUp(photon_package * pp)
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

bool CGridOcTree::isInside(const Vector3D & pos, const cell_basic & _cell) const
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

bool CGridOcTree::isInside(const Vector3D & pos) const
{
    if(pos.X() < cell_oc_root->getXmin() || pos.Y() < cell_oc_root->getYmin() ||
        pos.Z() < cell_oc_root->getZmin())
        return false;

    if(pos.X() > cell_oc_root->getXmax() || pos.Y() > cell_oc_root->getYmax() ||
        pos.Z() > cell_oc_root->getZmax())
        return false;

    return true;
}

// bool CGridOcTree::isInside(photon_package * pp, Vector3D & pos)
// {
//     cell_oc * cell = (cell_oc *)pp->getPositionCell();
//     if(pos.X() < cell->getXmin() || pos.Y() < cell->getYmin() || pos.Z() < cell->getZmin())
//         return false;

//     if(pos.X() > cell->getXmax() || pos.Y() > cell->getYmax() || pos.Z() > cell->getZmax())
//         return false;

//     return true;
// }

void CGridOcTree::setRndPositionInCell(photon_package * pp, CRandomGenerator * rand_gen)
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
