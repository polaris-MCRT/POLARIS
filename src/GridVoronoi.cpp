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

//    resetGridValues();

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
    double ref_length=0.01 * center_cell->getRefLength();

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
    
    if(pos == pp->getPosition())
    {
        path_length += ref_length;
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
