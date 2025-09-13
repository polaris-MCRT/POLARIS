/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "GridSpherical.hpp"
#include "CommandParser.hpp"
#include "CellOcTree.hpp"
#include "MathFunctions.hpp"
#include "Typedefs.hpp"
#include "Parameters.hpp"

bool CGridSpherical::loadGridFromBinaryFile(parameters & param, uint _data_len)
{
    ushort tmpID, tmpOffset;
    string filename = param.getPathGrid();

    uint r_counter = 0;
    uint ph_counter = 0;
    uint th_counter = 0;

    line_counter = 0;
    char_counter = 0;

    ifstream bin_reader(filename.c_str(), ios::in | ios::binary);

    if(bin_reader.fail())
    {
        cout << ERROR_LINE << "Cannot write to:\n Cannot load binary spherical grid file:" << endl;
        cout << filename << "\n" << endl;
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

    if(dataID == GRID_ID_SPH)
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
    }
    else
    {
        cout << ERROR_LINE << "Cannot write to:\n A spherical grid requires an ID of \"" << GRID_ID_SPH << "\"!"
             << endl;
        return false;
    }

    uint tmp_data_offset = validateDataPositions(param);
    if(tmp_data_offset == MAX_UINT)
        return false;

    bin_reader.read((char *)&Rmin, 8);
    bin_reader.read((char *)&Rmax, 8);
    bin_reader.read((char *)&N_r, 2);
    bin_reader.read((char *)&N_ph, 2);
    bin_reader.read((char *)&N_th, 2);
    bin_reader.read((char *)&log_factorR, 8);
    bin_reader.read((char *)&log_factorPh, 8);
    bin_reader.read((char *)&log_factorTh, 8);

    // Convert borders with conversion factors
    min_len = Rmin;
    Rmin *= conv_length_in_SI;
    Rmax *= conv_length_in_SI;

    total_volume = PIx4 * Rmax * Rmax * Rmax / 3.0;

    listR = new double[N_r + 1];
    listPh = new double[N_ph + 1];
    listTh = new double[N_th + 1];

    // --------------------------------------
    // ---------- Radial-direction ----------
    // --------------------------------------

    // Init radial cell border
    listR = new double[N_r + 1];
    if(log_factorR == 0)
    {
        // Allow user defined radius list, if log_factorR is zero

        // The global borders are already in the grid
        listR[0] = Rmin;
        listR[N_r] = Rmax;

        // Set the cell borders
        for(uint i_r = 1; i_r < N_r; i_r++)
        {
            // Read radial cell border position
            bin_reader.read((char *)&listR[i_r], 8);

            // Update radial position with conversion factors
            listR[i_r] *= conv_length_in_SI;
        }
    }
    else if(log_factorR == 1.0)
    {
        // Sinus shaped list, which emphasizes the middle rings
        CMathFunctions::SinList(Rmin, Rmax, listR, N_r + 1, log_factorR);
    }
    else if(log_factorR > 1.0)
    {
        // Exponentially increasing width of the cells in radial direction
        CMathFunctions::ExpList(Rmin, Rmax, listR, N_r + 1, log_factorR);
    }
    else
    {
        // Linear width of the cells in radial direction
        CMathFunctions::LinearList(Rmin, Rmax, listR, N_r + 1);
    }

    // -----------------------------------
    // ---------- Phi-direction ----------
    // -----------------------------------

    // Init phi cell border
    listPh = new double[N_ph + 1];
    if(log_factorPh == 0)
    {
        // Allow user defined phi list, if log_factorPh is zero

        // The global borders are already in the grid
        listPh[0] = 0;
        listPh[N_ph] = PIx2;

        // Set the cell borders
        for(uint i_ph = 1; i_ph < N_ph; i_ph++)
            bin_reader.read((char *)&listPh[i_ph], 8);
    }
    else
    {
        // Linear width of the cells in phi direction
        CMathFunctions::LinearList(0, PIx2, listPh, N_ph + 1);
    }

    // -------------------------------------
    // ---------- Theta-direction ----------
    // -------------------------------------
    if(log_factorTh == 0)
    {
        // The global borders are already in the grid
        listTh[0] = 0;
        listTh[N_th] = PI;

        // Read cell border in theta direction
        for(uint i_th = 1; i_th < N_th; i_th++)
            bin_reader.read((char *)&listTh[i_th], 8);
    }
    else if(log_factorTh == 1.0)
    {
        // Sinus shaped list, which emphasizes the midplane
        CMathFunctions::SinList(0, PI, listTh, N_th + 1, log_factorTh);
    }
    else if(log_factorTh > 1.0)
    {
        // Exponentially increasing width of the cells in z-direction (symmetrically)
        CMathFunctions::ExpListSym(0, PI, listTh, N_th + 1, log_factorTh);
    }
    else
    {
        // Linear width of the cells in theta direction
        CMathFunctions::LinearList(0, PI, listTh, N_th + 1);
    }

    // -----------------------------------------
    // ---------- Check of the limits ----------
    // -----------------------------------------
    if(Rmin <= 0)
    {
        cout << ERROR_LINE << "Cannot write to:\n Inner radius (Rmin = " << Rmin << ") must be larger than zero!"
             << endl;
        return false;
    }

    if(Rmax <= 0)
    {
        cout << ERROR_LINE << "Cannot write to:\n Outer radius (Rmax = " << Rmax << ") must be larger than zero!"
             << endl;
        return false;
    }

    if(Rmax <= Rmin)
    {
        cout << ERROR_LINE << "Cannot write to:\n Outer radius (Rmax = " << Rmax
             << ") must be larger than inner radius (Rmin = " << Rmin << ")!" << endl;
        return false;
    }

    // Init grid cells
    grid_cells = new cell_sp ***[N_r];

    for(uint i_r = 0; i_r < N_r; i_r++)
    {
        grid_cells[i_r] = new cell_sp **[N_ph];

        // cout << "Allocating memory for spherical grid cells: " << float(100.0 * double(i_r) / double(N_r))
        //      << "      \r" << flush;

        for(uint i_ph = 0; i_ph < N_ph; i_ph++)
        {
            grid_cells[i_r][i_ph] = new cell_sp *[N_th];

            for(uint i_th = 0; i_th < N_th; i_th++)
            {
                grid_cells[i_r][i_ph][i_th] = 0;
            }
        }
    }

    // Clear user output
    cout << CLR_LINE;

    max_cells = N_r * N_ph * N_th + 1;
    line_counter = -1;

    while(!bin_reader.eof())
    {
        line_counter++;

        if(r_counter < N_r)
        {
            double dr = listR[r_counter + 1] - listR[r_counter];

            if(dr == 0)
            {
                cout << ERROR_LINE << "No step size in r-direction of spherical grid!" << endl;
                return false;
            }

            if(dr < min_len)
                min_len = dr;
        }

        if(r_counter == 0)
        {
            double d;

            if(ph_counter < N_ph)
            {
                double dph = listPh[ph_counter + 1] - listPh[ph_counter];

                if(dph == 0)
                {
                    cout << ERROR_LINE << "No step size in phi-direction of spherical grid!" << endl;
                    return false;
                }

                if(N_ph > 2)
                {
                    d = 2 * Rmin * tan(dph / 2);
                    if(d < min_len)
                        min_len = d;
                }
            }

            if(th_counter < N_th)
            {
                double dth = listTh[th_counter + 1] - listTh[th_counter];

                if(dth == 0)
                {
                    cout << ERROR_LINE << "No step size in theta-direction of spherical grid!" << endl;
                    return false;
                }

                d = 2 * Rmin * tan(dth / 2);
                if(d < min_len)
                    min_len = d;
            }
        }

        // Calculate percentage of total progress per source
        float percentage = 100.0 * double(line_counter) / double(max_cells);

        // Show only new percentage number if it changed
        if((percentage - last_percentage) > PERCENTAGE_STEP)
        {
            char_counter++;
            cout << "-> Loading spherical grid file: " << percentage << " [%]      \r" << flush;
            last_percentage = percentage;
        }

        if(th_counter == N_th)
        {
            ph_counter++;
            th_counter = 0;
        }

        if(ph_counter == N_ph)
        {
            r_counter++;
            ph_counter = 0;
        }

        if(r_counter > N_r)
            break;

        cell_sp * tmp_cell = 0;

        if(r_counter == N_r)
        {
            center_cell = new cell_sp;
            center_cell->setRID(MAX_UINT);
            center_cell->setPhID(MAX_UINT);
            center_cell->setThID(MAX_UINT);
            center_cell->resize(data_len + tmp_data_offset);
            tmp_cell = center_cell;
            tmp_cell->setID(line_counter);
            r_counter++;
        }
        else
        {
            grid_cells[r_counter][ph_counter][th_counter] = new cell_sp;
            grid_cells[r_counter][ph_counter][th_counter]->setRID(r_counter);
            grid_cells[r_counter][ph_counter][th_counter]->setPhID(ph_counter);
            grid_cells[r_counter][ph_counter][th_counter]->setThID(th_counter);
            grid_cells[r_counter][ph_counter][th_counter]->resize(data_len + tmp_data_offset);
            tmp_cell = grid_cells[r_counter][ph_counter][th_counter];
            tmp_cell->setID(line_counter);
        }

        for(uint i = 0; i < data_offset; i++)
        {
            double tmp_data1 = 0;
            bin_reader.read((char *)&tmp_data1, 8);
            // cout << tmp_data1 << " ";
            tmp_cell->setData(i, tmp_data1);
        }

        updateVelocity(tmp_cell, param);

        if(uint(tmp_cell->getData(data_pos_id)) < 0 ||
           uint(tmp_cell->getData(data_pos_id)) > param.getMaxDustComponentChoice())
        {
            cout << ERROR_LINE << "Cannot write to:\n Dust ID in grid exceeds maximum number "
                    "of dust choices available! "
                 << endl;
            return false;
        }

        updateDataRange(tmp_cell);

        double tmp_vol = getVolume(*tmp_cell);
        total_gas_mass += getGasMassDensity(*tmp_cell) * tmp_vol;
        cell_volume += tmp_vol;
        th_counter++;
    }

    bin_reader.close();

    if(max_cells != uint(line_counter))
    {
        cout << ERROR_LINE << "Cannot write to:\n Number of read in cells do not match the "
                "maximal number of expected cells!"
             << endl;
        return false;
    }

    data_len += tmp_data_offset;
    data_offset += tmp_data_offset;

    max_len = 2 * Rmax;
    // min_len = listR[1] - listR[0];

    // cout << CLR_LINE;
    // cout << "- Loading spherical grid file          : done" << endl;

    return true;
}



bool CGridSpherical::saveBinaryGridFile(string filename, ushort id, ushort data_size)
{
    if(max_cells == 0)
    {
        cout << ERROR_LINE << "Cannot save spherical grid file to:" << endl;
        cout << filename;
        cout << "Not enough cells available! " << endl;
        return false;
    }

    ofstream bin_writer(filename.c_str(), ios::out | ios::binary);

    if(bin_writer.fail())
    {
        cout << ERROR_LINE << "Cannot write to:\n Cannot open spherical grid file:" << endl;
        cout << filename;
        return false;
    }

    bin_writer.write((char *)&id, 2);
    bin_writer.write((char *)&data_size, 2);

    if(dataID == GRID_ID_SPH)
    {
        for(uint i = 0; i < data_offset; i++)
        {
            ushort tmp_ids = data_ids[i];
            bin_writer.write((char *)&tmp_ids, 2);
        }
    }
    else
    {
        cout << ERROR_LINE << "Cannot save spherical grid file to:" << endl;
        cout << filename;
        cout << "A spherical grid requires an ID of " << GRID_ID_SPH << "!" << endl;
        return false;
    }

    bin_writer.write((char *)&Rmin, 8);
    bin_writer.write((char *)&Rmax, 8);
    bin_writer.write((char *)&N_r, 2);
    bin_writer.write((char *)&N_ph, 2);
    bin_writer.write((char *)&N_th, 2);
    bin_writer.write((char *)&log_factorR, 8);
    bin_writer.write((char *)&log_factorPh, 8);
    bin_writer.write((char *)&log_factorTh, 8);
    if(log_factorR == 0)
        for(uint i_r = 1; i_r < N_r; i_r++)
            bin_writer.write((char *)&listR[i_r], 8);
    if(log_factorPh == 0)
        for(uint i_ph = 1; i_ph < N_ph; i_ph++)
            bin_writer.write((char *)&listPh[i_ph], 8);
    if(log_factorTh == 0)
        for(uint i_th = 1; i_th < N_th; i_th++)
            bin_writer.write((char *)&listTh[i_th], 8);

    for(uint i_r = 0; i_r < N_r; i_r++)
    {
        // cout << "-> Writing binary spherical grid file: " << float(100.0 * double(i_r) / double(N_r))
        //      << "      \r" << flush;

        for(uint i_ph = 0; i_ph < N_ph; i_ph++)
        {
            for(uint i_th = 0; i_th < N_th; i_th++)
            {
                for(uint i = 0; i < data_offset; i++)
                {
                    double tmp_data = grid_cells[i_r][i_ph][i_th]->getData(i);
                    bin_writer.write((char *)&tmp_data, 8);
                }
            }
        }
    }

    for(uint i = 0; i < data_offset; i++)
    {
        double tmp_data = center_cell->getData(i);
        bin_writer.write((char *)&tmp_data, 8);
    }

    bin_writer.close();

    cout << CLR_LINE;
    cout << "- Writing spherical grid file   : done" << endl;

    return true;
}

bool CGridSpherical::createArtificialGrid(string path)
{
//    resetGridValues();

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

    min_len = 1e30;
    max_len = 2 * Rmax;

    string filename = path;

    dataID = GRID_ID_SPH;
    Rmin = 1; //* con_pc;
    Rmax = 4; //*con_pc;
    N_r = 3;
    N_ph = 4;
    N_th = 4;
    log_factorR = 0;
    log_factorPh = 0;
    log_factorTh = 0;

    ofstream bin_writer(filename.c_str(), ios::out | ios::binary);

    if(bin_writer.fail())
    {
        cout << ERROR_LINE << "Cannot write to:\n" << endl;
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

    bin_writer.write((char *)&Rmin, 8);
    bin_writer.write((char *)&Rmax, 8);
    bin_writer.write((char *)&N_r, 2);
    bin_writer.write((char *)&N_ph, 2);
    bin_writer.write((char *)&N_th, 2);
    bin_writer.write((char *)&log_factorR, 8);
    bin_writer.write((char *)&log_factorPh, 8);
    bin_writer.write((char *)&log_factorTh, 8);

    for(uint i_r = 0; i_r < N_r; i_r++)
    {
        // if(i_r % 50 == 0)
        //     cout << "-> Creating tree: " << 100.0 * float(i_r) / float(N_r) << " [%]           \r" << flush;

        for(uint i_ph = 0; i_ph < N_ph; i_ph++)
        {
            for(uint i_th = 0; i_th < N_th; i_th++)
            {
                max_cells++;

                double tmp_data = 1.0 / double(i_r * i_r + 1);

                bin_writer.write((char *)&tmp_data, 8);

                tmp_data = 1;
                bin_writer.write((char *)&tmp_data, 8);

                tmp_data = i_ph + 10;
                bin_writer.write((char *)&tmp_data, 8);

                for(uint i = 3; i < 9; i++)
                {
                    tmp_data = double(10);
                    bin_writer.write((char *)&tmp_data, 8);
                }
            }
        }
    }

    double tmp_data = 1e-5;
    bin_writer.write((char *)&tmp_data, 8);

    tmp_data = 1;
    bin_writer.write((char *)&tmp_data, 8);

    tmp_data = 2;
    bin_writer.write((char *)&tmp_data, 8);

    for(uint i = 3; i < 9; i++)
    {
        tmp_data = double(1);
        bin_writer.write((char *)&tmp_data, 8);
    }

    bin_writer.close();

    cout << "min: " << min_gas_dens << "  max_dens: " << max_gas_dens << endl;
    cout << "Creating artificial sphere                  : done" << endl;
    cout << "Max cells: " << max_cells << endl;

    return true;
}

void CGridSpherical::printParameters()
{
    if(max_cells == 0)
        cout << ERROR_LINE << "No tree parameters available! " << endl;
    else
    {
        cout << CLR_LINE;
        cout << "Spherical grid parameters (ID: " << getDataID() << ", data len.: " << getDataOffset()
             << ", Nr: " << N_r << ", Nph: " << N_ph << ", Nth: " << N_th << ")" << endl;
        cout << SEP_LINE;

        cout << "- Number of spherical cells     : " << max_cells << endl;
        printPhysicalParameters();
        cout << SEP_LINE;
    }
}

bool CGridSpherical::positionPhotonInGrid(photon_package * pp)
{
    uint dirID = pp->getDirectionID();
    uint rID = MAX_UINT, thID = MAX_UINT, phID = MAX_UINT;
    if(dirID < 6 && pp->getPositionCell() != 0)
    {
        cell_sp * tmp_cell = (cell_sp *)pp->getPositionCell();

        rID = tmp_cell->getRID();
        thID = tmp_cell->getThID();
        phID = tmp_cell->getPhID();

        // Update index of next cell
        switch(dirID)
        {
            case 0:
                rID--;
                break;

            case 1:
                rID++;
                break;

            case 2:
                thID--;
                break;

            case 3:
                thID++;
                break;

            case 4:
                if(phID == 0)
                    phID += N_ph;
                phID--;
                break;

            case 5:
                phID++;
                if(phID >= N_ph)
                    phID -= N_ph;
                break;

            default:
                return false;
                break;
        }
    }

    Vector3D pos = pp->getPosition();
    double sp_r = pos.length();

    if(sp_r < Rmin)
    {
        pp->setPositionCell(center_cell);
        return true;
    }

    uint i_r = 0, i_ph = 0, i_th = 0;

    i_r = CMathFunctions::biListIndexSearch(sp_r, listR, N_r + 1);
    if(i_r == MAX_UINT)
        return false;

    if(N_ph > 1)
    {
        double tmp_phi = pos.getPhiCoord();
        i_ph = CMathFunctions::biListIndexSearch(tmp_phi, listPh, N_ph + 1);
        if(i_ph == MAX_UINT)
            return false;
    }

    double tmp_theta = acos( pos.Z() / sp_r );

    i_th = CMathFunctions::biListIndexSearch(tmp_theta, listTh, N_th + 1);
    if(i_th == MAX_UINT)
        return false;

    pp->setPositionCell(grid_cells[i_r][i_ph][i_th]);

    return true;
}

bool CGridSpherical::goToNextCellBorder(photon_package * pp)
{
    cell_sp * tmp_cell = (cell_sp *)pp->getPositionCell();
    Vector3D p = pp->getPosition();
    Vector3D d = pp->getDirection();

    bool hit = false;
    double path_length = 1e300;
    uint dirID = MAX_UINT;

    uint rID = tmp_cell->getRID();

    if(rID == MAX_UINT)
    {
        double r2 = Rmin * (1 + MIN_LEN_STEP*EPS_DOUBLE);

        double B = p * d;
        double C = p.sq_length() - r2 * r2;
        // dscr is always >=0, we are inside the inner cell
        double dscr = B * B - C;

        dscr = sqrt(dscr);
        // "-"-solution is not needed for outer cells; only the "+"-solution can be correct
        double length = -B + dscr;

        if(length > 0 && length < path_length)
        {
            path_length = length;
            hit = true;
            dirID = 1;
        }
    }
    else
    {
        // --- Radial cell borders ---

        double r1 = listR[rID] * (1 - MIN_LEN_STEP*EPS_DOUBLE);
        double r2 = listR[rID + 1] * (1 + MIN_LEN_STEP*EPS_DOUBLE);

        double p_sq = p.sq_length();
        double B = p * d;
        double B_sq = pow(B, 2);

        double C1 = p_sq - r1 * r1;
        double C2 = p_sq - r2 * r2;

        double dscr1 = B_sq - C1;
        // dscr2 is always >= 0
        double dscr2 = B_sq - C2;

        if(dscr1 > 0)
        {
            dscr1 = sqrt(dscr1);
            // "+"-solution is not needed for inner cells; only the "-"-solution can be correct
            double length = -B - dscr1;

            if(length > 0 && length < path_length)
            {
                path_length = length;
                hit = true;
                dirID = 0;
            }
        }

        dscr2 = sqrt(dscr2);
        // "-"-solution is not needed for outer cells; only the "+"-solution can be correct
        double length = -B + dscr2;

        if(length != 0 && length < path_length)
        {
            path_length = length;
            hit = true;
            dirID = 1;
        }

        // --- Theta cell borders ---
        if(N_th > 1)
        {
            uint thID = tmp_cell->getThID();

            double th1 = listTh[thID] * (1 - MIN_LEN_STEP*EPS_DOUBLE);
            double cos_th1 = cos(th1);

            double cos_th1_sq = cos_th1 * cos_th1;
            double A1 = cos_th1_sq - d.Z() * d.Z();
            double B1 = cos_th1_sq * (d.X() * p.X() + d.Y() * p.Y()) - d.Z() * p.Z() * (1 - cos_th1_sq);
            double C3 = cos_th1_sq * (p.X() * p.X() + p.Y() * p.Y()) - p.Z() * p.Z() * (1 - cos_th1_sq);

            double dscr3 = B1 * B1 - A1 * C3;

            // dscr < 0 should not happen, but might if d.Z = 1, p = p.Z, and th1 = 0 or PI
            if(dscr3 >= 0)
            {
                dscr3 = sqrt(dscr3);

                double length[2];
                length[0] = (-B1 + dscr3) / A1;
                length[1] = (-B1 - dscr3) / A1;

                for(uint i=0; i<2; i++)
                    if(length[i] > 0 && length[i] < path_length)
                    {
                        path_length = length[i];
                        hit = true;
                        dirID = 2;
                    }
            }

            double th2 = listTh[thID + 1] * (1 + MIN_LEN_STEP*EPS_DOUBLE);
            double cos_th2 = cos(th2);

            double cos_th2_sq = cos_th2 * cos_th2;
            double A2 = cos_th2_sq - d.Z() * d.Z();
            double B2 = cos_th2_sq * (d.X() * p.X() + d.Y() * p.Y()) - d.Z() * p.Z() * (1 - cos_th2_sq);
            double C4 = cos_th2_sq * (p.X() * p.X() + p.Y() * p.Y()) - p.Z() * p.Z() * (1 - cos_th2_sq);

            double dscr4 = B2 * B2 - A2 * C4;

            // dscr < 0 should not happen, but might if d.Z = -1, p = p.Z, and th1 = 0 or PI
            if(dscr4 >= 0)
            {
                dscr4 = sqrt(dscr4);

                double length[2];
                length[0] = (-B2 + dscr4) / A2;
                length[1] = (-B2 - dscr4) / A2;

                for(uint i=0; i<2; i++)
                    if(length[i] > 0 && length[i] < path_length)
                    {
                        path_length = length[i];
                        hit = true;
                        dirID = 3;
                    }
            }
        }

        // --- Phi cell borders ---
        if(N_ph > 1)
        {
            uint phID = tmp_cell->getPhID();

            double r = sqrt(p.sq_length());
            double rho = sqrt(p.X() * p.X() + p.Y() * p.Y());

            double sin_th = rho / r;
            double cos_th = p.Z() / r;

            double ph1 = listPh[phID] * (1 - MIN_LEN_STEP*EPS_DOUBLE) - MIN_LEN_STEP*EPS_DOUBLE;
            double ph2 = listPh[phID + 1] * (1 + MIN_LEN_STEP*EPS_DOUBLE) + MIN_LEN_STEP*EPS_DOUBLE;

            double sin_ph1 = sin(ph1);
            double sin_ph2 = sin(ph2);

            double cos_ph1 = cos(ph1);
            double cos_ph2 = cos(ph2);

            Vector3D v_n1 = -Vector3D(-sin_ph1, cos_ph1, 0);
            Vector3D v_a1 = r * Vector3D(sin_th * cos_ph1, sin_th * sin_ph1, cos_th);

            double den1 = v_n1 * d;
            if(den1 != 0)
            {
                double num = v_n1 * (p - v_a1);
                double length = -num / den1;

                if(length > 0 && length < path_length)
                {
                    path_length = length;
                    hit = true;
                    dirID = 4;
                }
            }

            Vector3D v_n2 = Vector3D(-sin_ph2, cos_ph2, 0);
            Vector3D v_a2 = r * Vector3D(sin_th * cos_ph2, sin_th * sin_ph2, cos_th);

            double den2 = v_n2 * d;
            if(den2 != 0)
            {
                double num = v_n2 * (p - v_a2);
                double length = -num / den2;

                if(length > 0 && length < path_length)
                {
                    path_length = length;
                    hit = true;
                    dirID = 5;
                }
            }
        }
    }

    if(!hit)
    {
        cout << ERROR_LINE << "Wrong cell border!                                   " << endl;
        return false;
    }

    pp->setPosition(p + d * path_length);

    if(p == pp->getPosition())
    {
        cout << ERROR_LINE << "Could not transfer photon to the next cell border!   " << endl;
        return false;
    }

    pp->setTmpPathLength(path_length);
    pp->setDirectionID(dirID);
    return true;
}

bool CGridSpherical::updateShortestDistance(photon_package * pp)
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

bool CGridSpherical::findStartingPoint(photon_package * pp)
{
    Vector3D p = pp->getPosition();
    Vector3D d = pp->getDirection();

    if(isInside(p))
        return positionPhotonInGrid(pp);

    double path_length = 1e300;
    bool hit = false;

    double r2 = Rmax * (1 - MIN_LEN_STEP*EPS_DOUBLE);

    double B = p * d;
    // C is positive, we are outside of the cell
    double C = p.sq_length() - r2 * r2;
    double dscr = B * B - C;

    if(dscr > 0)
    {
        dscr = sqrt(dscr);
        // "+"-solution is not needed for inner cells; only the "-"-solution can be correct
        double length = -B - dscr;

        if(length > 0 && length < path_length)
        {
            path_length = length;
            hit = true;
        }
    }

    if(!hit)
        return false;

    pp->setPosition(p + d * path_length);
    pp->setDirectionID(MAX_UINT);
    return positionPhotonInGrid(pp);
}

Vector3D CGridSpherical::getCenter(const cell_basic & cell) const
{
    Vector3D center;
    const cell_sp * tmp_cell = (const cell_sp *)&cell;

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

void CGridSpherical::setRndPositionInCell(photon_package * pp, CRandomGenerator * rand_gen)
{
    Vector3D pos;
    cell_sp * tmp_cell = (cell_sp *)pp->getPositionCell();
    double r1, r2, ph1, ph2, th1, th2;

    double rnd_r = rand_gen->getRND();
    double rnd_ph = rand_gen->getRND();
    double rnd_th = rand_gen->getRND();

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

bool CGridSpherical::next(photon_package * pp)
{
    if(!positionPhotonInGrid(pp))
        return false;

    if(!goToNextCellBorder(pp))
        return false;

    return true;
}

/*
void CGridSpherical::getBoundingPoints(Vector3D & p_min, Vector3D & p_max)
{
    p_min.set(cell_oc_root->x_min, cell_oc_root->y_min,
            cell_oc_root->z_min);
    p_max.set(cell_oc_root->x_max, cell_oc_root->y_max,
            cell_oc_root->z_max);
}

void CGridSpherical::getBoundingPoints(cell_basic * cell, Vector3D & p_min,
        Vector3D & p_max)
{
    cell_oc * curr_cell = (cell_oc*) cell;
    p_min.set(curr_cell->x_min, curr_cell->y_min, curr_cell->z_min);
    p_max.set(curr_cell->x_max, curr_cell->y_max, curr_cell->z_max);
}
*/

void CGridSpherical::getLengths(uint bins, double & step_xy, double & off_xy)
{
    step_xy = 2 * Rmax / double(bins);
    off_xy = step_xy / 2.0;
}

bool CGridSpherical::createCellList()
{
    if(max_cells == 0)
    {
        cout << ERROR_LINE << "Spherical grid contains no cells!" << endl;
        cout << "       Cell list cannot be created!" << endl;
        return false;
    }

    cell_list = new cell_basic *[max_cells];
    ulong pos_counter = 0;

    // cout << CLR_LINE;
    // cout << "-> Creating cell list    : 0.0 [%]           \r" << flush;

    for(uint i_r = 0; i_r < N_r; i_r++)
    {
        // cout << "-> Creating cell list     : " << 100.0 * float(i_r) / float(N_r) << " %        \r"
        //      << flush;

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

double CGridSpherical::getVolume(const cell_basic & cell) const
{
    const cell_sp * cell_pos = (const cell_sp *)&cell;

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

Vector3D CGridSpherical::rotateToCenter(const photon_package & pp, Vector3D dir, bool inv, bool phi_only) const
{
    const cell_sp * cell_pos = (const cell_sp *)pp.getPositionCell();
    Vector3D pos = pp.getPosition().getSphericalCoord();

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

bool CGridSpherical::saveBinaryGridFile(string filename)
{
    return saveBinaryGridFile(filename, GRID_ID_SPH, data_offset);
}

bool CGridSpherical::loadGridFromBinaryFile(parameters & param)
{
    return loadGridFromBinaryFile(param, 0);
}

void CGridSpherical::clear()
{
    line_counter = 0;
    char_counter = 0;
    cout << "Final cleanup                                : done" << endl;
}

bool CGridSpherical::getPolarRTGridParameter(double max_len,
                                             double pixel_width,
                                             uint max_subpixel_lvl,
                                             dlist & _listR,
                                             uint & N_polar_r,
                                             uint *& N_polar_ph)
{
    return CGridBasic::getPolarRTGridParameterWorker(max_len,
                                                     pixel_width,
                                                     max_subpixel_lvl,
                                                     _listR,
                                                     N_polar_r,
                                                     N_polar_ph,
                                                     N_r,
                                                     listR);
}

bool CGridSpherical::isInside(const Vector3D & pos) const
{
    if(pos.sq_length() > Rmax * Rmax)
        return false;

    return true;
}
