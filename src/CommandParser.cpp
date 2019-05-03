#include "CommandParser.h"
#include <omp.h>

dlist CCommandParser::parseDataString(string data)
{
    string::size_type start, stop, pos;
    string col, val;
    uint ccount;
    dlist res;
    double tmp_val = -1;

    while(data.find("\"") != string::npos)
    {
        pos = data.find("\"");
        data.erase(pos, 1);
    }

    while(data.find("(") != string::npos)
    {
        start = data.find_first_of("(");
        stop = data.find_first_of(")");
        ccount = 0;
        tmp_val = -1;

        if(start == string::npos)
        {
            cout << "\nERROR: Missing \"(\" in line " << line_counter << "!\nColor bar set to default."
                 << endl;
            res.clear();
            return res;
        }

        if(stop == string::npos)
        {
            cout << "\nERROR: Missing \")\" in line " << line_counter << "!\nColor bar set to default.";
            res.clear();
            return res;
        }

        if(stop - start == 1)
        {
            cout << "\nERROR: No color values between (...) in line " << line_counter
                 << "!\nColor bar set to default." << endl;
            res.clear();
            return res;
        }

        col = data.substr(start + 1, stop - start - 1);
        data.erase(start, stop - start + 1);

        while(col.find(":") != string::npos)
        {
            ccount++;
            stop = col.find_first_of(":");
            val = col.substr(0, stop);
            col.erase(0, stop + 1);
            res.push_back(atof(val.c_str()));
        }

        if(ccount != 3 || col.size() == 0)
        {
            cout << "\nERROR: False amount of values between (...)  in line " << line_counter
                 << "!\nColor bar set to default." << endl;
            res.clear();
            return res;
        }

        res.push_back(atof(col.c_str()));
    }

    if(res.size() < 8)
    {
        cout << "\nERROR: At least two colors are required in line " << line_counter
             << "!\nColor bar set to default." << endl;
        res.clear();
        return res;
    }

    if(res.size() % 4 != 0)
    {
        cout << "\nERROR: False number of values between (...)  in line " << line_counter
             << "!\nColor bar set to default." << endl;
        res.clear();
        return res;
    }

    if(res[0] != 0)
    {
        cout << "\nERROR: First position in line " << line_counter
             << " has to be 0!\nColor bar set to default." << endl;
        res.clear();
        return res;
    }

    if(res[res.size() - 4] != 1)
    {
        cout << "\nERROR: Last position in line " << line_counter
             << " has to be 1!\nColor bar set to default." << endl;
        res.clear();
        return res;
    }

    for(uint i = 0; i < res.size(); i++)
    {
        if((i) % 4 == 0 || i == 0)
        {
            if(res[i] < 0 || res[i] > 1)
            {
                cout << "\nERROR: Position in color bar  in line " << line_counter
                     << " has to be between 0 and 1!\nColor bar set to default." << endl;
                res.clear();
                return res;
            }

            if(res[i] <= tmp_val)
            {
                cout << "\nERROR: Position values in line " << line_counter
                     << " have to be in ascending order!\nColor bar set to default." << endl;
                res.clear();
                return res;
            }

            tmp_val = res[i];
        }
        else
        {
            if(res[i] < 0 || res[i] > 255)
            {
                cout << "\nERROR: Color values in line " << line_counter
                     << " have to be between 0 and 255!\nColor bar set to default." << endl;
                res.clear();
                return res;
            }
        }
    }

    return res;
}

dlist CCommandParser::parseValues(string & str)
{
    uint pos;
    dlist values;
    string v;

    formatLine(str);

    if(str.size() == 0)
        return values;

    while(str.find(" ") != string::npos)
    {
        pos = uint(str.find(" "));
        v = str.substr(0, pos);
        str.erase(0, pos + 1);
        values.push_back(atof(v.c_str()));
    }

    values.push_back(atof(str.c_str()));

    return values;
}

void CCommandParser::formatLine(string & line)
{
    string::size_type pos = 0;
    string tmp_str = seperateString(line);

    if(line.size() == 0)
        return;

    if(line.find(">") != string::npos)
    {
        pos = line.find(">");
        line.replace(pos, 1, "> ");
    }

    if(line.find("=") != string::npos)
    {
        pos = line.find("=");
        line.replace(pos, 1, " = ");
    }

    while(line.find(";") != string::npos)
    {
        pos = line.find(";");
        line.replace(pos, 1, " ");
    }

    while(line.find("?") != string::npos)
    {
        pos = line.find("?");
        line.replace(pos, 1, " ");
    }

    while(line.find("*") != string::npos)
    {
        pos = line.find("*");
        line.replace(pos, 1, " ");
    }

    while(line.find('\t') != string::npos)
    {
        pos = line.find('\t');
        line.replace(pos, 1, " ");
    }

    while(line.find(" \r\n") != string::npos)
    {
        pos = line.find(" \r\n");
        line.replace(pos, 3, " ");
    }

    while(line.find(" \r") != string::npos)
    {
        pos = line.find(" \r");
        line.replace(pos, 2, " ");
    }

    while(line.find(" \n") != string::npos)
    {
        pos = line.find(" \n");
        line.replace(pos, 2, " ");
    }

    while(line.find("\r\n") != string::npos)
    {
        pos = line.find("\r\n");
        line.replace(pos, 2, " ");
    }

    while(line.find("\r") != string::npos)
    {
        pos = line.find("\r");
        line.replace(pos, 1, " ");
    }

    while(line.find("\n") != string::npos)
    {
        pos = line.find("\n");
        line.replace(pos, 1, " ");
    }
    //***

    while(line.find("  ") != string::npos)
    {
        pos = line.find("  ");
        line.replace(pos, 2, " ");
    }

    while(line.find(",") != string::npos)
    {
        pos = line.find(",");
        line.replace(pos, 1, ".");
    }

    if(line == " ")
        line = "";

    if(line.size() > 0)
    {
        while(line.c_str()[line.size() - 1] == ' ')
        {
            pos = line.find_last_of(' ');
            line.erase(pos, 1);
        }
    }

    while(line.c_str()[0] == ' ')
    {
        // pos=line.find_first_of(' ');
        line.erase(0, 1);
    }

    if(line.find_first_of("#") != string::npos)
    {
        pos = line.find("#");
        if(pos == 0)
            tmp_str = "";

        line.erase(pos, line.length() - pos);
    }

    if(line.find_first_of("!") != string::npos)
    {
        pos = line.find("!");
        if(pos == 0)
            tmp_str = "";

        line.erase(pos, line.length() - pos);
    }

    if(tmp_str.size() != 0)
        line += " \"" + tmp_str + "\"";
}

bool CCommandParser::parse()
{
    tag = 1000;
    int task_id = 0;
    uint task_pos = -1;
    bool apply = true;
    ifstream reader(cmd_filename.c_str());
    string str_common = "";
    vector<string> tasks;

    string line, command;
    string::size_type pos = 0;

    line_counter = 0;

    if(reader.fail())
    {
        cout << "\nERROR: Cannot open file:\n" << cmd_filename << endl;
        return false;
    }

    while(getline(reader, line))
    {
        line_counter++;

        if(line.compare("") == 0)
            continue;
        if(line.c_str()[0] == '#')
            continue;

        formatLine(line);
        pos = line.find(">");

        if(pos != (uint)(-1))
        {
            command = line.substr(0, pos + 1);
            line.erase(0, pos + 2);
        }
        else
        {
            command = line;
            line = "";
        }

        if(command.compare("") == 0)
            continue;
        if(command.compare(" ") == 0)
            continue;

        if(command.c_str()[0] != '<')
        {
            cout << "\nWARNING: Unknown command " << command << " in line " << line_counter << "!" << endl;
            continue;
        }

        if(command.compare("<task>") == 0)
        {
            if(task_id == -id_tsk)
            {
                cout << "\nERROR: Closing tag \"</task>\" is missing in line " << line_counter << "!" << endl;
                return false;
            }

            if(task_id == -id_cmn)
            {
                cout << "\nERROR: Closing tag \"</common>\" is missing in line " << line_counter << "!"
                     << endl;
                return false;
            }

            if(line.size() == 0)
                apply = true;
            else
            {
                apply = atob(atoi(line.c_str()));
            }

            command = "";
            line = "";
            task_id = -id_tsk;

            if(apply)
            {
                task_pos++;
                tasks.push_back("");
            }
        }

        if(command.compare("</task>") == 0)
        {
            if(task_id != -id_tsk)
            {
                cout << "\nERROR: Opening tag \"<task>\" is missing in line " << line_counter << "!" << endl;
                return false;
            }

            command = "";
            line = "";
            task_id = id_tsk;
        }

        if(command.compare("<common>") == 0)
        {
            if(task_id == -id_tsk)
            {
                cout << "\nERROR: Closing tag \"</task>\" is missing in line " << line_counter << "!" << endl;
                return false;
            }

            if(task_id == -id_cmn)
            {
                cout << "\nERROR: Closing tag \"</common>\" is missing in line " << line_counter << "!"
                     << endl;
                return false;
            }

            command = "";
            line = "";
            task_id = -id_cmn;
        }

        if(command.compare("</common>") == 0)
        {
            if(task_id != -id_cmn)
            {
                cout << "\nERROR: Opening tag \"<common>\" is missing in line " << line_counter << "!"
                     << endl;
                return false;
            }

            command = "";
            line = "";
            task_id = id_cmn;
        }

        switch(task_id)
        {
            case -id_tsk:
                if(apply)
                    tasks[task_pos] += (command + " " + line + "?");
                break;

            case -id_cmn:
                str_common += (command + " " + line + "?");
                break;
        };
    }

    param_list.resize(tasks.size());
    tasks.insert(tasks.begin(), str_common);

    parameters common_param;

    for(uint i = 0; i < tasks.size(); i++)
    {
        line_counter = 0;
        int i_pos = i - 1;
        string str_task = tasks[i];
        string::size_type pos1 = 0, pos2 = 0;
        string ret = "";
        int len = -1;
        parameters * param;

        if(i == 0)
            param = &common_param;
        else
            param = &param_list[i_pos];

        param->setTaskID(i);

        while(str_task.find("?") != string::npos)
        {
            pos1 = str_task.find("?");
            pos2 = str_task.find("?", pos1 + 1);

            len = int(pos2 - pos1 + 1);

            if(len < 0)
                break;

            line = str_task.substr(pos1, len);
            str_task.erase(pos1, len - 1);

            line_counter++;

            if(line.compare("") == 0)
                continue;
            if(line.c_str()[0] == '#')
                continue;

            formatLine(line);
            pos = line.find(">");

            if(pos != (uint)(-1))
            {
                command = line.substr(0, pos + 1);
                line.erase(0, pos + 2);
            }
            else
            {
                command = line;
                line = "";
            }

            if(command.compare("") == 0)
                continue;
            if(command.compare(" ") == 0)
                continue;

            if(command.c_str()[0] != '<')
            {
                cout << "\nWARNING: Unknown command " << command << " in line " << line_counter << "!"
                     << endl;
                continue;
            }

            if(!parseLine(param, command, line, i))
            {
                cout << "\nWARNING: Unknown command " << command << " in line " << line_counter << "!"
                     << endl;
                continue;
            }
        }

        uint start = param->getStart();
        uint stop = param->getStop();
        uint map_max = 1;

        switch(param->getCommand())
        {
            case CMD_DUST_EMISSION:
                if(i == 0)
                    break;

                map_max = uint(param->getDustRayDetectors().size() / NR_OF_RAY_DET);

                if(start == MAX_UINT)
                    start = 0;
                else
                    start--;

                if(stop == MAX_UINT)
                    stop = map_max - 1;
                else
                    stop--;

                if(stop > map_max - 1)
                {
                    stop = map_max - 1;
                    cout << "\nWARNING: <stop> value larger than number of raytracing "
                            "detectors!"
                         << endl;
                    cout << " Value set to " << stop + 1 << "." << endl;
                }

                if(start > map_max - 1)
                {
                    start = 0;
                    cout << "\nWARNING: <start> value larger than number of raytracing "
                            "detectors!"
                         << endl;
                    cout << " Value set to 1." << endl;
                }

                if(map_max == 0)
                {
                    start = 0;
                    stop = 0;
                }
                break;

            case CMD_SYNCHROTRON:
                if(i == 0)
                    break;

                map_max = uint(param->getSyncRayDetectors().size() / NR_OF_RAY_DET);

                if(start == MAX_UINT)
                    start = 0;
                else
                    start--;

                if(stop == MAX_UINT)
                    stop = map_max - 1;
                else
                    stop--;

                if(stop > map_max - 1)
                {
                    stop = map_max - 1;
                    cout << "\nWARNING: <stop> value larger than number of raytracing "
                            "detectors!"
                         << endl;
                    cout << " Value set to " << stop + 1 << "." << endl;
                }

                if(start > map_max - 1)
                {
                    start = 0;
                    cout << "\nWARNING: <start> value larger than number of raytracing "
                            "detectors!"
                         << endl;
                    cout << " Value set to 1." << endl;
                }

                if(map_max == 0)
                {
                    start = 0;
                    stop = 0;
                }
                break;

            case CMD_LINE_EMISSION:
                if(i == 0)
                    break;

                map_max = param->getNrOfGasSpecies();

                if(start == MAX_UINT)
                    start = 0;
                else
                    start--;

                if(stop == MAX_UINT)
                    stop = map_max - 1;
                else
                    stop--;

                if(stop > map_max - 1)
                {
                    stop = map_max - 1;
                    cout << "\nWARNING: <stop> value larger than number of gas species!" << endl;
                    cout << " Value set to " << stop + 1 << "." << endl;
                }

                if(start > map_max - 1)
                {
                    start = 0;
                    cout << "\nWARNING: <start> value larger than number of gas species!" << endl;
                    cout << " Value set to 1." << endl;
                }
                break;
        }

        param->setStart(start);
        param->setStop(stop);

        if(i == 0)
        {
            for(uint j = 0; j < param_list.size(); j++)
                param_list[j] = common_param;
        }
        else
            param_list[i_pos] = *param;
    }

    return true;
}

string CCommandParser::seperateString(string & str)
{
    string::size_type pos1 = 0, pos2 = 0;
    string ret = "";
    int len = -1;

    if(str.find_first_of("\"") != string::npos)
    {
        pos1 = str.find("\"");
        pos2 = str.find("\"", pos1 + 1);

        len = int(pos2 - pos1 + 1);

        if(len < 0)
            return ret;

        ret = str.substr(pos1, len);
        str.erase(pos1, len);
    }

    while(ret.find("\"") != string::npos)
    {
        pos1 = ret.find("\"");
        ret.erase(pos1, 1);
    }

    return ret;
}

bool CCommandParser::parseLine(parameters * param, string cmd, string data, uint id)
{
    if(cmd.compare("<cmd>") == 0)
    {
        if(data.compare("CMD_DUST_EMISSION") == 0)
        {
            param->setCommand(CMD_DUST_EMISSION);
            return true;
        }

        if(data.compare("CMD_OPIATE") == 0)
        {
            param->setCommand(CMD_OPIATE);
            return true;
        }

        if(data.compare("CMD_FORCE") == 0)
        {
            param->setCommand(CMD_FORCE);
            return true;
        }

        if(data.compare("CMD_LINE_EMISSION") == 0)
        {
            param->setCommand(CMD_LINE_EMISSION);
            return true;
        }

        if(data.compare("CMD_TEMP") == 0)
        {
            param->setCommand(CMD_TEMP);
            return true;
        }

        if(data.compare("CMD_TEMP_RAT") == 0)
        {
            param->setCommand(CMD_TEMP_RAT);
            return true;
        }

        /*if(data.compare("CMD_PLOT") == 0)
        {
            param->setCommand(CMD_PLOT);
            return true;
        }*/

        if(data.compare("CMD_RAT") == 0)
        {
            param->setCommand(CMD_RAT);
            return true;
        }

        if(data.compare("CMD_DUST_SCATTERING") == 0)
        {
            param->setCommand(CMD_DUST_SCATTERING);
            return true;
        }

        if(data.compare("CMD_SYNCHROTRON") == 0)
        {
            param->setCommand(CMD_SYNCHROTRON);
            return true;
        }

        cout << "\nERROR: Command cannot be recognized!" << endl;
        return false;
    }

    if(cmd.compare("<delta0>") == 0)
    {
        double val = atof(data.c_str());
        param->setDelta0(val);
        return true;
    }

    if(cmd.compare("<larm_f>") == 0)
    {
        double val = atof(data.c_str());
        param->setLarmF(val);
        return true;
    }

    if(cmd.compare("<plot_list>") == 0)
    {
        dlist ids = parseValues(data);

        if(ids.empty())
        {
            cout << "\nWARNING: The list of plot IDs is empty! " << endl;
            cout << "         Only integers and no text is allowed!    " << endl;
            // return false;
        }

        for(uint i = 0; i < ids.size(); i++)
        {
            uint id = ids[i];

            if(id >= minGRID && id <= maxGRID)
                param->addToPlotList(id);
            else
            {
                cout << "\nWARNING: Unknown grid ID in line " << line_counter << "!    " << endl;
                cout << "         A plot ID of " << id
                     << " is not a valid POLARIS grid ID (see manual, Table 3.3)!     " << endl;
                // return false;
            }
        }

        return true;
    }

    if(cmd.compare("<plot_list>") == 0)
    {
        dlist ids = parseValues(data);

        if(ids.empty())
        {
            cout << "\nWARNING: The list of plot IDs is empty! " << endl;
            cout << "         Only integers and no text is allowed!    " << endl;
            // return false;
        }

        for(uint i = 0; i < ids.size(); i++)
        {
            uint id = ids[i];

            if(id >= minGRID && id <= maxGRID)
                param->addToPlotList(id);
            else
            {
                cout << "\nWARNING: Unknown grid ID in line " << line_counter << "!    " << endl;
                cout << "         A plot ID of " << id
                     << " is not a valid POLARIS grid ID (see manual, Table 3.3)!     " << endl;
                // return false;
            }
        }

        return true;
    }

    if(cmd.compare("<phase_function>") == 0)
    {
        if(data.compare("PH_HG") == 0)
        {
            param->setPhaseFunctionID(PH_HG);
            return true;
        }

        if(data.compare("PH_MIE") == 0)
        {
            param->setPhaseFunctionID(PH_MIE);
            return true;
        }

        if(data.compare("PH_ISO") == 0)
        {
            param->setPhaseFunctionID(PH_ISO);
            return true;
        }

        cout << "\nERROR: Phase function cannot be recognized!" << endl;
        return false;
    }

    if(cmd.compare("<star_mass>") == 0)
    {
        dlist mass = parseValues(data);

        for(uint i = 0; i < mass.size(); i++)
        {
            param->addStarMass(mass[i] * M_sun);
        }
        return true;
    }

    if(cmd.compare("<detector_line nr_pixel = vel_channels = >") == 0)
    {
        string str1 = seperateString(data);
        string str2 = seperateString(data);
        formatLine(str1);
        dlist nr_of_pixel = parseValues(str1);
        formatLine(str2);
        dlist nr_of_channels = parseValues(str2);

        formatLine(data);
        dlist values = parseValues(data);

        while(values[3] < 0)
            values[3] += 360;
        while(values[4] < 0)
            values[4] += 360;

        if(values.size() == NR_OF_LINE_DET - 11)
        {
            // Only gas_species_id, transition_id, source_id, max_velocity, rot_angle_1
            // and rot_angle_2 Set distance to 1
            values.push_back(1.0);
            // Set sidelength in x-direction to cube sidelength
            values.push_back(-1);
            // Set sidelength in y-direction to cube sidelength
            values.push_back(-1);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_LINE_DET - 10)
        {
            // As above, but with distance to observer
            // Set sidelength in x-direction of cube sidelength
            values.push_back(-1);
            // Set sidelength in y-direction of cube sidelength
            values.push_back(-1);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_LINE_DET - 9)
        {
            // As above, but with one sidelength for both directions of cube sidelength
            // Set given sidelength also for y-direction of cube sidelength
            values.push_back(values[NR_OF_LINE_DET - 10]);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_LINE_DET - 8)
        {
            // As above, but with two sidelengths for x- and y-directions of cube
            // sidelength Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_LINE_DET - 6)
        {
            // As above, but with x- and y-shift of the detector map
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
        }

        values.push_back(DET_PLANE);
        if(!checkPixel(values, nr_of_pixel))
            return false;
        if(!checkVelChannels(values, nr_of_channels))
            return false;

        // HINT: +1 because the gas species id is not saved in the gas_species list!
        if(values.size() != (NR_OF_LINE_DET + 1))
        {
            cout << "\nERROR: Number of parameters in line detector could not be "
                    "recognized!"
                 << endl;
            return false;
        }

        param->addLineRayDetector(values);
        param->updateDetectorAngles(values[4], values[5]);
        param->updateObserverDistance(values[6]);
        param->updateMapSidelength(values[7], values[8]);
        param->updateRayGridShift(values[9], values[10]);
        param->updateDetectorPixel(uint(values[NR_OF_LINE_DET - 2]), uint(values[NR_OF_LINE_DET - 1]));

        return true;
    }

    if(cmd.compare("<detector_line_healpix nr_sides = vel_channels = >") == 0)
    {
        string str1 = seperateString(data);
        string str2 = seperateString(data);
        formatLine(str1);
        dlist nr_of_sides = parseValues(str1);
        formatLine(str2);
        dlist nr_of_channels = parseValues(str2);

        formatLine(data);
        dlist values = parseValues(data);

        if(values.size() == NR_OF_LINE_DET - 10)
        {
            // Only gas_species_id, transition_id, source_id, max_velocity, position X, Y,
            // and Z Set galactic coordinate l (Longitude) to [-180°, 180°]
            values.push_back(-180);
            values.push_back(180);
            // Set galactic coordinate b (Latitude) to [-90°, 90°]
            values.push_back(-90);
            values.push_back(90);
            // Set velocity of the observer to (0, 0, 0) [m/s]
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_LINE_DET - 7)
        {
            // As above, but with galactic coordinates
            // Set velocity of the observer to (0, 0, 0) [m/s]
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
        }

        values.push_back(DET_SPHER);
        if(!checkPixel(values, nr_of_sides, true))
            return false;
        if(!checkVelChannels(values, nr_of_channels))
            return false;

        if(values.size() != (NR_OF_LINE_DET + 1))
        {
            cout << "\nERROR: Number of parameters in healpix line detector could not be "
                    "recognized!"
                 << endl;
            return false;
        }

        param->addLineRayDetector(values);
        param->updateDetectorPixel(uint(nr_of_sides[0]), 0);

        double distance = sqrt(values[4] * values[4] + values[5] * values[5] + values[6] * values[6]);
        param->updateObserverDistance(distance);

        // Showing full sphere coverage
        param->updateDetectorAngles(-90, -180);
        param->updateDetectorAngles(90, 180);

        return true;
    }

    if(cmd.compare("<detector_line_polar nr_pixel = vel_channels = >") == 0)
    {
        string str1 = seperateString(data);
        string str2 = seperateString(data);
        formatLine(str1);
        dlist nr_of_pixel = parseValues(str1);
        formatLine(str2);
        dlist nr_of_channels = parseValues(str2);

        formatLine(data);
        dlist values = parseValues(data);

        while(values[3] < 0)
            values[3] += 360;
        while(values[4] < 0)
            values[4] += 360;

        if(values.size() == NR_OF_LINE_DET - 11)
        {
            // Only gas_species_id, transition_id, source_id, max_velocity, rot_angle_1
            // and rot_angle_2 Set distance to 1
            values.push_back(1.0);
            // Set sidelength in x-direction to cube sidelength
            values.push_back(-1);
            // Set sidelength in y-direction to cube sidelength
            values.push_back(-1);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_LINE_DET - 10)
        {
            // As above, but with distance to observer
            // Set sidelength in x-direction of cube sidelength
            values.push_back(-1);
            // Set sidelength in y-direction of cube sidelength
            values.push_back(-1);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_LINE_DET - 9)
        {
            // As above, but with one sidelength for both directions of cube sidelength
            // Set given sidelength also for y-direction of cube sidelength
            values.push_back(values[NR_OF_LINE_DET - 10]);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
        }

        values.push_back(DET_POLAR);
        if(!checkPixel(values, nr_of_pixel))
            return false;
        if(!checkVelChannels(values, nr_of_channels))
            return false;

        if(values.size() != (NR_OF_LINE_DET + 1))
        {
            cout << "\nERROR: Number of parameters in polar line detector could not be "
                    "recognized!"
                 << endl;
            return false;
        }

        param->addLineRayDetector(values);
        param->updateDetectorAngles(values[4], values[5]);
        param->updateObserverDistance(values[6]);
        param->updateMapSidelength(values[7], values[8]);
        param->updateDetectorPixel(uint(values[NR_OF_LINE_DET - 2]), uint(values[NR_OF_LINE_DET - 1]));

        return true;
    }

    if(cmd.compare("<detector_line_slice nr_pixel = vel_channels = >") == 0)
    {
        string str1 = seperateString(data);
        string str2 = seperateString(data);
        formatLine(str1);
        dlist nr_of_pixel = parseValues(str1);
        formatLine(str2);
        dlist nr_of_channels = parseValues(str2);

        formatLine(data);
        dlist values = parseValues(data);

        while(values[3] < 0)
            values[3] += 360;
        while(values[4] < 0)
            values[4] += 360;

        if(values.size() == NR_OF_LINE_DET - 11)
        {
            // Only gas_species_id, transition_id, source_id, max_velocity, rot_angle_1
            // and rot_angle_2 Set distance to 1
            values.push_back(1.0);
            // Set sidelength in x-direction to cube sidelength
            values.push_back(-1);
            // Set sidelength in y-direction to cube sidelength
            values.push_back(-1);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_LINE_DET - 10)
        {
            // As above, but with distance to observer
            // Set sidelength in x-direction of cube sidelength
            values.push_back(-1);
            // Set sidelength in y-direction of cube sidelength
            values.push_back(-1);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_LINE_DET - 9)
        {
            // As above, but with one sidelength for both directions of cube sidelength
            // Set given sidelength also for y-direction of cube sidelength
            values.push_back(values[NR_OF_LINE_DET - 10]);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_LINE_DET - 8)
        {
            // As above, but with two sidelengths for x- and y-directions of cube
            // sidelength Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_LINE_DET - 6)
        {
            // As above, but with x- and y-shift of the detector map
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
            values.push_back(0.0);
        }

        values.push_back(DET_SLICE);
        if(!checkPixel(values, nr_of_pixel))
            return false;
        if(!checkVelChannels(values, nr_of_channels))
            return false;

        // HINT: +1 because the gas species id is not saved in the gas_species list!
        if(values.size() != (NR_OF_LINE_DET + 1))
        {
            cout << "\nERROR: Number of parameters in line detector could not be "
                    "recognized!"
                 << endl;
            return false;
        }

        param->addLineRayDetector(values);
        param->updateDetectorAngles(values[4], values[5]);
        param->updateObserverDistance(values[6]);
        param->updateMapSidelength(values[7], values[8]);
        param->updateRayGridShift(values[9], values[10]);
        param->updateDetectorPixel(uint(values[NR_OF_LINE_DET - 2]), uint(values[NR_OF_LINE_DET - 1]));

        return true;
    }

    if(cmd.compare("<detector_dust nr_pixel = >") == 0)
    {
        string str = seperateString(data);
        formatLine(str);
        dlist nr_of_pixel = parseValues(str);

        formatLine(data);
        dlist values = parseValues(data);

        while(values[4] < 0)
            values[4] += 360;
        while(values[5] < 0)
            values[5] += 360;

        if(values.size() == NR_OF_RAY_DET - 8)
        {
            // Only wl_min, wl_max, wl_skip, source_id, rot_angle_1 and rot_angle_2
            // Set distance to 1
            values.push_back(1.0);
            // Set sidelength in x-direction to cube sidelength
            values.push_back(-1);
            // Set sidelength in y-direction to cube sidelength
            values.push_back(-1);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_RAY_DET - 7)
        {
            // As above, but with distance to observer
            // Set sidelength in x-direction of cube sidelength
            values.push_back(-1);
            // Set sidelength in y-direction of cube sidelength
            values.push_back(-1);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_RAY_DET - 6)
        {
            // As above, but with one sidelength for both directions of cube sidelength
            // Set given sidelength also for y-direction of cube sidelength
            values.push_back(values[NR_OF_RAY_DET - 7]);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_RAY_DET - 5)
        {
            // As above, but with two sidelengths for x- and y-directions of cube
            // sidelength Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }

        values.push_back(DET_PLANE);
        if(!checkPixel(values, nr_of_pixel))
            return false;

        if(values.size() != NR_OF_RAY_DET)
        {
            cout << "\nERROR: Number of parameters in raytracing detector could not be "
                    "recognized!"
                 << endl;
            return false;
        }

        if(values[2] < 1)
        {
            cout << "\nERROR: Number of wavelengths needs to be at least 1!" << endl;
            return false;
        }
        else if(values[2] > 1 && values[0] == values[1])
        {
            cout << "\nERROR: Minimum and maximum wavelength cannot be the "
                 << "same if the number of wavelengths is larger than 1!" << endl;
            return false;
        }

        param->addDustRayDetector(values);
        param->updateDetectorAngles(values[4], values[5]);
        param->updateObserverDistance(values[6]);
        param->updateMapSidelength(values[7], values[8]);
        param->updateRayGridShift(values[9], values[10]);
        param->updateDetectorPixel(uint(values[NR_OF_RAY_DET - 2]), uint(values[NR_OF_RAY_DET - 1]));

        return true;
    }

    if(cmd.compare("<detector_dust_healpix nr_sides = >") == 0)
    {
        string str = seperateString(data);
        formatLine(str);
        dlist nr_of_sides = parseValues(str);

        if(nr_of_sides.size() != 1)
        {
            cout << "\nERROR: Number of sides in healpix raytracing detector could not "
                    "be recognized!"
                 << endl;
            return false;
        }

        if(nr_of_sides[0] <= 0)
        {
            cout << "\nERROR: Number of sides in healpix raytracing detector could not "
                    "be recognized!"
                 << endl;
            return false;
        }

        if(!CMathFunctions::isPowerOfTwo(int(nr_of_sides[0])))
        {
            cout << "\nERROR: Number of sides must be a power of two!" << endl;
            return false;
        }

        formatLine(data);
        dlist values = parseValues(data);

        if(values.size() == NR_OF_RAY_DET - 7)
        {
            // Only wl_min, wl_max, wl_skip, source_id, position X, Y, and Z
            // Set galactic coordinate l (Longitude) to [-180°, 180°]
            values.push_back(-180);
            values.push_back(180);
            // Set galactic coordinate b (Latitude) to [-90°, 90°]
            values.push_back(-90);
            values.push_back(90);
        }

        values.push_back(DET_SPHER);
        if(!checkPixel(values, nr_of_sides, true))
            return false;

        if(values.size() != NR_OF_RAY_DET)
        {
            cout << "\nERROR: Number of parameters in healpix raytracing detector could "
                    "not be recognized!"
                 << endl;
            return false;
        }

        if(values[2] < 1)
        {
            cout << "\nERROR: Number of wavelengths needs to be at least 1!" << endl;
            return false;
        }
        else if(values[2] > 1 && values[0] == values[1])
        {
            cout << "\nERROR: Minimum and maximum wavelength cannot be the "
                 << "same if the number of wavelengths is larger than 1!" << endl;
            return false;
        }

        param->addDustRayDetector(values);
        param->updateDetectorPixel(uint(nr_of_sides[0]), 0);

        double distance = sqrt(values[3] * values[3] + values[4] * values[4] + values[5] * values[5]);
        param->updateObserverDistance(distance);

        // Showing full sphere coverage
        param->updateDetectorAngles(0, 0);
        param->updateDetectorAngles(180, 360);

        return true;
    }

    if(cmd.compare("<detector_dust_polar nr_pixel = >") == 0)
    {
        string str = seperateString(data);
        formatLine(str);
        dlist nr_of_pixel = parseValues(str);

        formatLine(data);
        dlist values = parseValues(data);

        while(values[4] < 0)
            values[4] += 360;
        while(values[5] < 0)
            values[5] += 360;

        if(values.size() == NR_OF_RAY_DET - 8)
        {
            // Only wl_min, wl_max, wl_skip, source_id, rot_angle_1 and rot_angle_2
            // Set distance to 1
            values.push_back(1.0);
            // Set sidelength in x-direction to cube sidelength
            values.push_back(-1);
            // Set sidelength in y-direction to cube sidelength
            values.push_back(-1);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_RAY_DET - 7)
        {
            // As above, but with distance to observer
            // Set sidelength in x-direction of cube sidelength
            values.push_back(-1);
            // Set sidelength in y-direction of cube sidelength
            values.push_back(-1);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_RAY_DET - 6)
        {
            // As above, but with one sidelength for both directions of cube sidelength
            // Set given sidelength also for y-direction of cube sidelength
            values.push_back(values[NR_OF_RAY_DET - 7]);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_RAY_DET - 5)
        {
            // As above, but with two sidelengths for x- and y-directions of cube
            // sidelength Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }

        values.push_back(DET_POLAR);
        if(!checkPixel(values, nr_of_pixel))
            return false;

        if(values.size() != NR_OF_RAY_DET)
        {
            cout << "\nERROR: Number of parameters in polar raytracing detector could "
                    "not be recognized!"
                 << endl;
            return false;
        }

        if(values[2] < 1)
        {
            cout << "\nERROR: Number of wavelengths needs to be at least 1!" << endl;
            return false;
        }
        else if(values[2] > 1 && values[0] == values[1])
        {
            cout << "\nERROR: Minimum and maximum wavelength cannot be the "
                 << "same if the number of wavelengths is larger than 1!" << endl;
            return false;
        }

        param->addDustRayDetector(values);
        param->updateDetectorAngles(values[4], values[5]);
        param->updateObserverDistance(values[6]);
        param->updateMapSidelength(values[7], values[8]);
        param->updateDetectorPixel(uint(values[NR_OF_RAY_DET - 2]), uint(values[NR_OF_RAY_DET - 1]));

        return true;
    }

    if(cmd.compare("<detector_dust_slice nr_pixel = >") == 0)
    {
        string str = seperateString(data);
        formatLine(str);
        dlist nr_of_pixel = parseValues(str);

        formatLine(data);
        dlist values = parseValues(data);

        while(values[4] < 0)
            values[4] += 360;
        while(values[5] < 0)
            values[5] += 360;

        if(values.size() == NR_OF_RAY_DET - 8)
        {
            // Only wl_min, wl_max, wl_skip, source_id, rot_angle_1 and rot_angle_2
            // Set distance to 1
            values.push_back(1.0);
            // Set sidelength in x-direction to cube sidelength
            values.push_back(-1);
            // Set sidelength in y-direction to cube sidelength
            values.push_back(-1);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_RAY_DET - 7)
        {
            // As above, but with distance to observer
            // Set sidelength in x-direction of cube sidelength
            values.push_back(-1);
            // Set sidelength in y-direction of cube sidelength
            values.push_back(-1);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_RAY_DET - 6)
        {
            // As above, but with one sidelength for both directions of cube sidelength
            // Set given sidelength also for y-direction of cube sidelength
            values.push_back(values[NR_OF_RAY_DET - 7]);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_RAY_DET - 5)
        {
            // As above, but with two sidelengths for x- and y-directions of cube
            // sidelength Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }

        values.push_back(DET_SLICE);
        if(!checkPixel(values, nr_of_pixel))
            return false;

        if(values.size() != NR_OF_RAY_DET)
        {
            cout << "\nERROR: Number of parameters in raytracing detector could not be "
                    "recognized!"
                 << endl;
            return false;
        }

        if(values[2] < 1)
        {
            cout << "\nERROR: Number of wavelengths needs to be at least 1!" << endl;
            return false;
        }
        else if(values[2] > 1 && values[0] == values[1])
        {
            cout << "\nERROR: Minimum and maximum wavelength cannot be the "
                 << "same if the number of wavelengths is larger than 1!" << endl;
            return false;
        }

        param->addDustRayDetector(values);
        param->updateDetectorAngles(values[4], values[5]);
        param->updateObserverDistance(values[6]);
        param->updateMapSidelength(values[7], values[8]);
        param->updateRayGridShift(values[9], values[10]);
        param->updateDetectorPixel(uint(values[NR_OF_RAY_DET - 2]), uint(values[NR_OF_RAY_DET - 1]));

        return true;
    }

    if(cmd.compare("<detector_dust_mc nr_pixel = >") == 0)
    {
        string str = seperateString(data);
        formatLine(str);
        dlist nr_of_pixel = parseValues(str);

        formatLine(data);
        dlist values = parseValues(data);

        while(values[3] < 0)
            values[3] += 360;
        while(values[4] < 0)
            values[4] += 360;

        if(values.size() == NR_OF_MC_DET - 4)
        {
            // Only wl_min, wl_max, wl_skip, source_id, rot_angle_1 and rot_angle_2
            // Set sidelength in x-direction of cube sidelength
            values.push_back(-1);
            // Set sidelength in y-direction of cube sidelength
            values.push_back(-1);
        }
        else if(values.size() == NR_OF_MC_DET - 3)
        {
            // As above, but with one sidelength for both directions of cube sidelength
            // Set given sidelength also for y-direction of cube sidelength
            values.push_back(values[NR_OF_MC_DET - 4]);
        }

        if(!checkPixel(values, nr_of_pixel))
            return false;

        if(values.size() != NR_OF_MC_DET)
        {
            cout << "\nERROR: Number of pixel in Monte-Carlo detector could not be "
                    "recognized!"
                 << endl;
            return false;
        }

        if(values[2] < 1)
        {
            cout << "\nERROR: Number of wavelengths needs to be at least 1!" << endl;
            return false;
        }
        else if(values[2] > 1 && values[0] == values[1])
        {
            cout << "\nERROR: Minimum and maximum wavelength cannot be the "
                 << "same if the number of wavelengths is larger than 1!" << endl;
            return false;
        }

        param->addDustMCDetector(values);
        param->updateDetectorAngles(values[3], values[4]);
        param->updateObserverDistance(values[5]);
        param->updateMapSidelength(values[6], values[7]);
        param->updateDetectorPixel(uint(values[NR_OF_MC_DET - 2]), uint(values[NR_OF_MC_DET - 1]));

        return true;
    }

    if(cmd.compare("<detector_sync nr_pixel = >") == 0)
    {
        string str = seperateString(data);
        formatLine(str);
        dlist nr_of_pixel = parseValues(str);

        formatLine(data);
        dlist values = parseValues(data);

        while(values[4] < 0)
            values[4] += 360;
        while(values[5] < 0)
            values[5] += 360;

        if(values.size() == NR_OF_RAY_DET - 8)
        {
            // Only wl_min, wl_max, wl_skip, source_id, rot_angle_1 and rot_angle_2
            // Set distance to 1
            values.push_back(1.0);
            // Set sidelength in x-direction to cube sidelength
            values.push_back(-1);
            // Set sidelength in y-direction to cube sidelength
            values.push_back(-1);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_RAY_DET - 7)
        {
            // As above, but with distance to observer
            // Set sidelength in x-direction of cube sidelength
            values.push_back(-1);
            // Set sidelength in y-direction of cube sidelength
            values.push_back(-1);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_RAY_DET - 6)
        {
            // As above, but with one sidelength for both directions of cube sidelength
            // Set given sidelength also for y-direction of cube sidelength
            values.push_back(values[NR_OF_RAY_DET - 7]);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_RAY_DET - 5)
        {
            // As above, but with two sidelengths for x- and y-directions of cube
            // sidelength Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }

        values.push_back(DET_PLANE);
        if(!checkPixel(values, nr_of_pixel))
            return false;

        if(values.size() != NR_OF_RAY_DET)
        {
            cout << "\nERROR: Number of parameters in raytracing detector could not be "
                    "recognized!"
                 << endl;
            return false;
        }

        if(values[2] < 1)
        {
            cout << "\nERROR: Number of wavelengths needs to be at least 1!" << endl;
            return false;
        }
        else if(values[2] > 1 && values[0] == values[1])
        {
            cout << "\nERROR: Minimum and maximum wavelength cannot be the "
                 << "same if the number of wavelengths is larger than 1!" << endl;
            return false;
        }

        param->addSyncRayDetector(values);
        param->updateDetectorAngles(values[4], values[5]);
        param->updateObserverDistance(values[6]);
        param->updateMapSidelength(values[7], values[8]);
        param->updateRayGridShift(values[9], values[10]);
        param->updateDetectorPixel(uint(values[NR_OF_RAY_DET - 2]), uint(values[NR_OF_RAY_DET - 1]));

        return true;
    }

    if(cmd.compare("<detector_sync_slice nr_pixel = >") == 0)
    {
        string str = seperateString(data);
        formatLine(str);
        dlist nr_of_pixel = parseValues(str);

        formatLine(data);
        dlist values = parseValues(data);

        while(values[4] < 0)
            values[4] += 360;
        while(values[5] < 0)
            values[5] += 360;

        if(values.size() == NR_OF_RAY_DET - 8)
        {
            // Only wl_min, wl_max, wl_skip, source_id, rot_angle_1 and rot_angle_2
            // Set distance to 1
            values.push_back(1.0);
            // Set sidelength in x-direction to cube sidelength
            values.push_back(-1);
            // Set sidelength in y-direction to cube sidelength
            values.push_back(-1);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_RAY_DET - 7)
        {
            // As above, but with distance to observer
            // Set sidelength in x-direction of cube sidelength
            values.push_back(-1);
            // Set sidelength in y-direction of cube sidelength
            values.push_back(-1);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_RAY_DET - 6)
        {
            // As above, but with one sidelength for both directions of cube sidelength
            // Set given sidelength also for y-direction of cube sidelength
            values.push_back(values[NR_OF_RAY_DET - 7]);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_RAY_DET - 5)
        {
            // As above, but with two sidelengths for x- and y-directions of cube
            // sidelength Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }

        values.push_back(DET_SLICE);
        if(!checkPixel(values, nr_of_pixel))
            return false;

        if(values.size() != NR_OF_RAY_DET)
        {
            cout << "\nERROR: Number of parameters in raytracing detector could not be "
                    "recognized!"
                 << endl;
            return false;
        }

        if(values[2] < 1)
        {
            cout << "\nERROR: Number of wavelengths needs to be at least 1!" << endl;
            return false;
        }
        else if(values[2] > 1 && values[0] == values[1])
        {
            cout << "\nERROR: Minimum and maximum wavelength cannot be the "
                 << "same if the number of wavelengths is larger than 1!" << endl;
            return false;
        }

        param->addSyncRayDetector(values);
        param->updateDetectorAngles(values[4], values[5]);
        param->updateObserverDistance(values[6]);
        param->updateMapSidelength(values[7], values[8]);
        param->updateRayGridShift(values[9], values[10]);
        param->updateDetectorPixel(uint(values[NR_OF_RAY_DET - 2]), uint(values[NR_OF_RAY_DET - 1]));

        return true;
    }

    if(cmd.compare("<detector_sync_healpix nr_sides = >") == 0)
    {
        string str = seperateString(data);
        formatLine(str);
        dlist nr_of_sides = parseValues(str);

        if(nr_of_sides.size() != 1)
        {
            cout << "\nERROR: Number of sides in healpix raytracing detector could not "
                    "be recognized!"
                 << endl;
            return false;
        }

        if(nr_of_sides[0] <= 0)
        {
            cout << "\nERROR: Number of sides in healpix raytracing detector could not "
                    "be recognized!"
                 << endl;
            return false;
        }

        if(!CMathFunctions::isPowerOfTwo(int(nr_of_sides[0])))
        {
            cout << "\nERROR: Number of sides must be a power of two!" << endl;
            return false;
        }

        formatLine(data);
        dlist values = parseValues(data);

        if(values.size() == NR_OF_RAY_DET - 7)
        {
            // Only wl_min, wl_max, wl_skip, source_id, position X, Y, and Z
            // Set galactic coordinate l (Longitude) to [-180°, 180°]
            values.push_back(-180);
            values.push_back(180);
            // Set galactic coordinate b (Latitude) to [-90°, 90°]
            values.push_back(-90);
            values.push_back(90);
        }

        values.push_back(DET_SPHER);
        if(!checkPixel(values, nr_of_sides, true))
            return false;

        if(values.size() != NR_OF_RAY_DET)
        {
            cout << "\nERROR: Number of parameters in healpix raytracing detector could "
                    "not be recognized!"
                 << endl;
            return false;
        }

        if(values[2] < 1)
        {
            cout << "\nERROR: Number of wavelengths needs to be at least 1!" << endl;
            return false;
        }
        else if(values[2] > 1 && values[0] == values[1])
        {
            cout << "\nERROR: Minimum and maximum wavelength cannot be the "
                 << "same if the number of wavelengths is larger than 1!" << endl;
            return false;
        }

        param->addSyncRayDetector(values);
        param->updateDetectorPixel(uint(12 * nr_of_sides[0] * nr_of_sides[0]), 0);

        double distance = sqrt(values[3] * values[3] + values[4] * values[4] + values[5] * values[5]);
        param->updateObserverDistance(distance);

        // Showing full sphere coverage
        param->updateDetectorAngles(0, 0);
        param->updateDetectorAngles(180, 360);

        return true;
    }

    if(cmd.compare("<detector_dust_polar nr_pixel = >") == 0)
    {
        string str = seperateString(data);
        formatLine(str);
        dlist nr_of_pixel = parseValues(str);

        formatLine(data);
        dlist values = parseValues(data);

        while(values[4] < 0)
            values[4] += 360;
        while(values[5] < 0)
            values[5] += 360;

        if(values.size() == NR_OF_RAY_DET - 8)
        {
            // Only wl_min, wl_max, wl_skip, source_id, rot_angle_1 and rot_angle_2
            // Set distance to 1
            values.push_back(1.0);
            // Set sidelength in x-direction to cube sidelength
            values.push_back(-1);
            // Set sidelength in y-direction to cube sidelength
            values.push_back(-1);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_RAY_DET - 7)
        {
            // As above, but with distance to observer
            // Set sidelength in x-direction of cube sidelength
            values.push_back(-1);
            // Set sidelength in y-direction of cube sidelength
            values.push_back(-1);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_RAY_DET - 6)
        {
            // As above, but with one sidelength for both directions of cube sidelength
            // Set given sidelength also for y-direction of cube sidelength
            values.push_back(values[NR_OF_RAY_DET - 7]);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_RAY_DET - 5)
        {
            // As above, but with two sidelengths for x- and y-directions of cube
            // sidelength Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }

        values.push_back(DET_POLAR);
        if(!checkPixel(values, nr_of_pixel))
            return false;

        if(values.size() != NR_OF_RAY_DET)
        {
            cout << "\nERROR: Number of parameters in polar raytracing detector could "
                    "not be recognized!"
                 << endl;
            return false;
        }

        if(values[2] < 1)
        {
            cout << "\nERROR: Number of wavelengths needs to be at least 1!" << endl;
            return false;
        }
        else if(values[2] > 1 && values[0] == values[1])
        {
            cout << "\nERROR: Minimum and maximum wavelength cannot be the "
                 << "same if the number of wavelengths is larger than 1!" << endl;
            return false;
        }

        param->addSyncRayDetector(values);
        param->updateDetectorAngles(values[4], values[5]);
        param->updateObserverDistance(values[6]);
        param->updateMapSidelength(values[7], values[8]);
        param->updateDetectorPixel(uint(values[NR_OF_RAY_DET - 2]), uint(values[NR_OF_RAY_DET - 1]));

        return true;
    }

    if(cmd.compare("<detector_dust_slice nr_pixel = >") == 0)
    {
        string str = seperateString(data);
        formatLine(str);
        dlist nr_of_pixel = parseValues(str);

        formatLine(data);
        dlist values = parseValues(data);

        while(values[4] < 0)
            values[4] += 360;
        while(values[5] < 0)
            values[5] += 360;

        if(values.size() == NR_OF_RAY_DET - 8)
        {
            // Only wl_min, wl_max, wl_skip, source_id, rot_angle_1 and rot_angle_2
            // Set distance to 1
            values.push_back(1.0);
            // Set sidelength in x-direction to cube sidelength
            values.push_back(-1);
            // Set sidelength in y-direction to cube sidelength
            values.push_back(-1);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_RAY_DET - 7)
        {
            // As above, but with distance to observer
            // Set sidelength in x-direction of cube sidelength
            values.push_back(-1);
            // Set sidelength in y-direction of cube sidelength
            values.push_back(-1);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_RAY_DET - 6)
        {
            // As above, but with one sidelength for both directions of cube sidelength
            // Set given sidelength also for y-direction of cube sidelength
            values.push_back(values[NR_OF_RAY_DET - 7]);
            // Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }
        else if(values.size() == NR_OF_RAY_DET - 5)
        {
            // As above, but with two sidelengths for x- and y-directions of cube
            // sidelength Do not use the other values
            values.push_back(0.0);
            values.push_back(0.0);
        }

        values.push_back(DET_SLICE);
        if(!checkPixel(values, nr_of_pixel))
            return false;

        if(values.size() != NR_OF_RAY_DET)
        {
            cout << "\nERROR: Number of parameters in raytracing detector could not be "
                    "recognized!"
                 << endl;
            return false;
        }

        if(values[2] < 1)
        {
            cout << "\nERROR: Number of wavelengths needs to be at least 1!" << endl;
            return false;
        }
        else if(values[2] > 1 && values[0] == values[1])
        {
            cout << "\nERROR: Minimum and maximum wavelength cannot be the "
                 << "same if the number of wavelengths is larger than 1!" << endl;
            return false;
        }

        param->addSyncRayDetector(values);
        param->updateDetectorAngles(values[4], values[5]);
        param->updateObserverDistance(values[6]);
        param->updateMapSidelength(values[7], values[8]);
        param->updateRayGridShift(values[9], values[10]);
        param->updateDetectorPixel(uint(values[NR_OF_RAY_DET - 2]), uint(values[NR_OF_RAY_DET - 1]));

        return true;
    }

    if(cmd.compare("<gas_species>") == 0)
    {
        string::size_type pos = 0;

        pos = data.find("POP_LTE");
        if(pos != string::npos)
        {
            string str_id = " ";
            str_id[0] = char(POP_LTE + 48);
            data.replace(pos, 7, str_id);
        }

        pos = data.find("POP_FEP");
        if(pos != string::npos)
        {
            string str_id = " ";
            str_id[0] = char(POP_FEP + 48);
            data.replace(pos, 7, str_id);
        }

        pos = data.find("POP_LVG");
        if(pos != string::npos)
        {
            string str_id = " ";
            str_id[0] = char(POP_LVG + 48);
            data.replace(pos, 7, str_id);
        }

        string str1 = seperateString(data);
        string str2 = seperateString(data);

        string gas_species_path, zeeman_path;

        if(str2.size() == 0)
            gas_species_path = str1;
        else
        {
            gas_species_path = str1;
            zeeman_path = str2;
        }

        dlist values = parseValues(data);
        if(values.size() == 2)
            param->addGasSpecies(gas_species_path, zeeman_path, values);
        else
        {
            cout << "\nERROR: False amount of parameters for gas species line transfer "
                    "in line "
                 << line_counter << "!" << endl;
            return false;
        }

        return true;
    }

    if(cmd.compare("<source_star nr_photons = >") == 0)
    {
        string str = seperateString(data);
        ullong nr_of_photons = ullong(atof(str.c_str()));
        string ps_path = seperateString(data);

        if(nr_of_photons <= 0)
        {
            cout << "\nERROR: Number of star photons could not be recognized!" << endl;
            return false;
        }

        dlist values = parseValues(data);

        if(ps_path.size() != 0)
        {
            if(values.size() == NR_OF_POINT_SOURCES - 5)
            {
                values.push_back(0);
                values.push_back(0);
                values.push_back(0);
                values.push_back(0);
            }
            else
            {
                cout << "\nERROR: False amount of parameters for source star in line " << line_counter << "!"
                     << endl;
                return false;
            }
        }
        else
        {
            if(values.size() == NR_OF_POINT_SOURCES - 3)
            {
                values.push_back(0);
                values.push_back(0);
            }
            else if(values.size() != NR_OF_POINT_SOURCES - 1)
            {
                cout << "\nERROR: False amount of parameters for source star in line " << line_counter << "!"
                     << endl;
                return false;
            }
        }

        double P_l = sqrt(pow(values[5], 2) + pow(values[6], 2));
        if(P_l > 1.0)
        {
            cout << "\nERROR: Chosen polarization of source star is larger than 1!" << endl;
            return false;
        }
        else if(P_l < 0)
        {
            cout << "\nHINT: Chosen polarization of source star is less than 0 (now set "
                    "to 0)!"
                 << endl;
            values[5] = 0;
            values[6] = 0;
        }

        values.push_back(double(nr_of_photons));
        param->addPointSource(values, ps_path);
        return true;
    }

    if(cmd.compare("<source_starfield nr_photons = >") == 0)
    {
        string str = seperateString(data);
        ullong nr_of_photons = ullong(atof(str.c_str()));
        string ps_path = seperateString(data);

        if(nr_of_photons <= 0)
        {
            cout << "\nERROR: Number of star photons could not be recognized!" << endl;
            return false;
        }

        dlist values = parseValues(data);

        if(ps_path.size() != 0)
        {
            if(values.size() == NR_OF_DIFF_SOURCES - 5)
            {
                values.push_back(0);
                values.push_back(0);
                values.push_back(0);
                values.push_back(0);
            }
            else
            {
                cout << "\nERROR: False amount of parameters for source starfield in line " << line_counter
                     << "!" << endl;
                return false;
            }
        }
        else
        {
            if(values.size() == NR_OF_DIFF_SOURCES - 3)
            {
                values.push_back(0);
                values.push_back(0);
            }
            else if(values.size() != NR_OF_DIFF_SOURCES - 1)
            {
                cout << "\nERROR: False amount of parameters for source starfield in line " << line_counter
                     << "!" << endl;
                return false;
            }
        }

        double P_l = sqrt(pow(values[6], 2) + pow(values[7], 2));
        if(P_l > 1.0)
        {
            cout << "\nERROR: Chosen polarization of source star is larger than 1!" << endl;
            return false;
        }
        else if(P_l < 0)
        {
            cout << "\nHINT: Chosen polarization of source star is less than 0 (now set "
                    "to 0)!"
                 << endl;
            values[6] = 0;
            values[7] = 0;
        }

        values.push_back(double(nr_of_photons));
        param->addDiffuseSource(values, ps_path);
        return true;
    }

    if(cmd.compare("<source_background nr_photons = >") == 0)
    {
        string str = seperateString(data);
        ullong nr_of_photons = ullong(atof(str.c_str()));

        if(nr_of_photons <= 0)
        {
            cout << "\nERROR: Number of point source photons could not be recognized!" << endl;
            return false;
        }

        string path = seperateString(data);
        dlist values = parseValues(data);
        values.push_back(nr_of_photons);

        if(values.size() == NR_OF_BG_SOURCES)
            param->addBackgroundSource(values);
        else
        {
            if(path.size() > 3 && values.size() == NR_OF_BG_SOURCES - 5)
                param->addBackgroundSource(path, values);
            else
            {
                if(path.size() > 3 && values.size() < NR_OF_BG_SOURCES - 5)
                    param->addBackgroundSource(path);
                else
                {
                    if(values.size() == NR_OF_BG_SOURCES - 2)
                    {
                        values.push_back(0);
                        values.push_back(0);

                        param->addBackgroundSource(values);
                    }
                    else
                    {
                        cout << "\nERROR: Cannot detect path for background parameters "
                                "file in line: ";
                        cout << line_counter << endl;
                        return false;
                    }
                }
            }
        }
        return true;
    }

    if(cmd.compare("<source_laser nr_photons = >") == 0)
    {
        string str = seperateString(data);
        ullong nr_of_photons = ullong(atof(str.c_str()));

        if(nr_of_photons <= 0)
        {
            cout << "\nERROR: Number of laser photons could not be recognized!" << endl;
            return false;
        }

        dlist values = parseValues(data);

        if(values.size() == NR_OF_LASER_SOURCES - 3)
        {
            values.push_back(0);
            values.push_back(0);
        }
        else if(values.size() != NR_OF_LASER_SOURCES - 1)
        {
            cout << "\nERROR: False amount of parameters for source laser in line " << line_counter << "!"
                 << endl;
            return false;
        }

        double P_l = sqrt(pow(values[9], 2) + pow(values[10], 2));
        if(P_l > 1.0)
        {
            cout << "\nERROR: Chosen polarization of source laser is larger than 1!" << endl;
            return false;
        }

        values.push_back(double(nr_of_photons));
        param->addLaserSource(values);
        return true;
    }

    if(cmd.compare("<axis1>") == 0)
    {
        dlist values = parseValues(data);

        if(values.size() != 3)
        {
            cout << "\nERROR: Values for first axis are not a vector!" << endl;
            return false;
        }

        param->setAxis1(values[0], values[1], values[2]);

        return true;
    }

    if(cmd.compare("<axis2>") == 0)
    {
        dlist values = parseValues(data);

        if(values.size() != 3)
        {
            cout << "\nERROR: Values for second axis are not a vector!" << endl;
            return false;
        }

        param->setAxis2(values[0], values[1], values[2]);
        return true;
    }

    if(cmd.compare("<align>") == 0)
    {
        if(data.compare("ALIG_INTERNAL") == 0)
        {
            param->addAlignmentMechanism(ALIG_INTERNAL);
            return true;
        }

        if(data.compare("ALIG_PA") == 0)
        {
            param->addAlignmentMechanism(ALIG_PA);
            return true;
        }

        if(data.compare("ALIG_IDG") == 0)
        {
            param->addAlignmentMechanism(ALIG_IDG);
            return true;
        }

        if(data.compare("ALIG_RAT") == 0)
        {
            param->addAlignmentMechanism(ALIG_RAT);
            return true;
        }

        if(data.compare("ALIG_GOLD") == 0)
        {
            param->addAlignmentMechanism(ALIG_GOLD);
            return true;
        }
    }

    if(cmd.compare("<mu>") == 0)
    {
        param->setMu(atof(data.c_str()));
        return true;
    }

    if(cmd.compare("<xy_min>") == 0)
    {
        param->setXYMin(atof(data.c_str()));
        param->setAutoScale(false);
        return true;
    }

    if(cmd.compare("<xy_max>") == 0)
    {
        param->setXYMax(atof(data.c_str()));
        param->setAutoScale(false);
        return true;
    }

    if(cmd.compare("<xy_steps>") == 0)
    {
        param->setXYSteps(uint(atof(data.c_str())));
        param->setAutoScale(false);
        return true;
    }

    if(cmd.compare("<xy_bins>") == 0)
    {
        param->setXYBins(uint(atof(data.c_str())));
        return true;
    }

    if(cmd.compare("<xy_label>") == 0)
    {
        string label = seperateString(data);
        param->setXYLabel(label);
        param->setAutoScale(false);
        return true;
    }

    if(cmd.compare("<path_grid>") == 0)
    {
        string path = seperateString(data);
        param->setPathGrid(path);
        return true;
    }

    if(cmd.compare("<path_grid_cgs>") == 0)
    {
        string path = seperateString(data);
        param->setPathGrid(path);

        param->updateSIConvDH(1e3);
        param->updateSIConvLength(1e-2);
        param->updateSIConvBField(1e-4);
        param->updateSIConvVField(1e-2);

        return true;
    }

    if(cmd.compare("<sub_dust>") == 0)
    {
        param->setSublimate(atob(atoi(data.c_str())));
        return true;
    }

    if(cmd.compare("<vel_maps>") == 0)
    {
        bool b = atob(atoi(data.c_str()));
        param->setVelMaps(b);
        return true;
    }

    if(cmd.compare("<max_subpixel_lvl>") == 0)
    {
        param->setMaxSubpixelLvl(int(atof(data.c_str())));
        return true;
    }

    if(cmd.compare("<path_input>") == 0)
    {
        string path = seperateString(data);
        param->setPathInput(path);
        return true;
    }

    if(cmd.compare("<dust_component>") == 0 || cmd.compare("<dust_component id = >") == 0)
    {
        string size_keyword, path;

        uint dust_component_choice = 0;
        if(cmd.compare("<dust_component id = >") == 0)
        {
            string str1 = seperateString(data);
            string str2 = seperateString(data);
            string str3 = seperateString(data);
            if(str3.size() != 0)
            {
                size_keyword = str1;
                dust_component_choice = uint(atof(str2.c_str()));
                path = str3;
            }
            else
            {
                size_keyword = "";
                dust_component_choice = uint(atof(str1.c_str()));
                path = str2;
            }

            if(dust_component_choice < 0)
            {
                cout << "Error ID " << dust_component_choice << " from <dust_component id = > is not valid!"
                     << endl;
                return false;
            }
        }
        else
        {
            path = seperateString(data);
            size_keyword = seperateString(data);
        }
        param->AddDustComponentChoice(dust_component_choice);

        uint nr_size_parameter = 0;
        if(size_keyword.size() >= 3)
        {
            if(size_keyword.find("plaw") != std::string::npos)
            {
                nr_size_parameter++;
                if(size_keyword.find("-ed") != std::string::npos)
                    nr_size_parameter += 3;
                if(size_keyword.find("-cv") != std::string::npos)
                    nr_size_parameter += 3;
            }
            else if(size_keyword.find("logn") != std::string::npos)
                nr_size_parameter += 2;
            else if(size_keyword.find("zda") != std::string::npos)
                nr_size_parameter += 14;
            else
            {
                cout << "\nERROR: Dust size distribution keyword not known!" << endl;
                return false;
            }
        }
        else
            size_keyword = "plaw";
        dlist fr = parseValues(data);

        dlist size_parameter;
        if(path.size() > 3 && nr_size_parameter > 0 && fr.size() == nr_size_parameter + 4)
        {
            for(uint i = 0; i < nr_size_parameter; i++)
                size_parameter.push_back(fr[i + 4]);
            for(uint i = nr_size_parameter; i < NR_OF_SIZE_DIST_PARAM; i++)
                size_parameter.push_back(0);
            param->addDustComponent(path, size_keyword, fr[0], fr[1], fr[2], fr[3], size_parameter);
        }
        else if(path.size() > 3 && nr_size_parameter == 0)
        {
            size_parameter.push_back(-3.5);
            for(uint i = 1; i < NR_OF_SIZE_DIST_PARAM; i++)
                size_parameter.push_back(0);
            if(fr.size() == 4)
            {
                size_parameter[0] = fr[1];
                param->addDustComponent(path, size_keyword, fr[0], 0, fr[2], fr[3], size_parameter);
            }
            else if(fr.size() == 2)
            {
                size_parameter[0] = fr[1];
                param->addDustComponent(path, size_keyword, fr[0], 0, 0, 0, size_parameter);
            }
            else if(fr.size() == 1)
                param->addDustComponent(path, size_keyword, fr[0], 0, 0, 0, size_parameter);
            else if(fr.size() == 0)
                param->addDustComponent(path, size_keyword, 1.0, 0, 0, 0, size_parameter);
        }
        else
        {
            cout << "\nWARNING: False parameters set for dust component in line " << line_counter << "!"
                 << endl;
            return false;
        }

        return true;
    }

    if(cmd.compare("<path_out>") == 0)
    {
        string path = seperateString(data);
        param->setPathOutput(path);
        return true;
    }

    if(cmd.compare("<opiata_param_path>") == 0)
    {
        string path = seperateString(data);
        param->setOpiateParamPath(path);
        return true;
    }

    if(cmd.compare("<opiate_data_path>") == 0)
    {
        string path = seperateString(data);
        param->setOpiateDataPath(path);
        return true;
    }

    if(cmd.compare("<nr_gnu_points>") == 0)
    {
        param->setNrOfGnuPoints(uint(atof(data.c_str())));
        return true;
    }

    if(cmd.compare("<nr_gnu_vectors>") == 0)
    {
        param->setNrOfGnuVectors(uint(atof(data.c_str())));
        return true;
    }

    if(cmd.compare("<f_highJ>") == 0)
    {
        param->setFhighJ(atof(data.c_str()));
        return true;
    }

    if(cmd.compare("<f_c>") == 0)
    {
        param->setFcorr(atof(data.c_str()));
        return true;
    }

    if(cmd.compare("<adj_tgas>") == 0)
    {
        param->setAdjTgas(atof(data.c_str()));
        return true;
    }

    if(cmd.compare("<max_lines>") == 0)
    {
        param->setmaxGridLines(uint(atof(data.c_str())));
        return true;
    }

    if(cmd.compare("<start>") == 0)
    {
        param->setStart(uint(atof(data.c_str())));
        return true;
    }

    if(cmd.compare("<stop>") == 0)
    {
        param->setStop(uint(atof(data.c_str())));
        return true;
    }

    if(cmd.compare("<conv_dens>") == 0)
    {
        double value = atof(data.c_str());

        if(value < 0)
        {
            cout << "\nWARNING: Negative conversion factors are no longer supported!" << endl;
            cout << "         Grid must always contain number densities" << endl;
            value = -value;
        }

        param->updateSIConvDH(atof(data.c_str()));
        return true;
    }

    if(cmd.compare("<conv_len>") == 0)
    {
        param->updateSIConvLength(atof(data.c_str()));
        return true;
    }

    if(cmd.compare("<conv_mag>") == 0)
    {
        param->updateSIConvBField(atof(data.c_str()));
        return true;
    }

    if(cmd.compare("<conv_vel>") == 0)
    {
        double conv = atof(data.c_str());

        if(conv < 0)
        {
            cout << "\nWARNING: Negative conversion factor are no longer allowed!" << endl;
            cout << "         The grid can only contain number densities.                "
                    "    \n"
                 << endl;
            conv = abs(conv);
        }

        param->updateSIConvVField(conv);
        return true;
    }

    if(cmd.compare("<mass_fraction>") == 0)
    {
        if(atof(data.c_str()) == 0)
        {
            // Use the dust fractions as dust-to-gas mass ratios
            param->setDustMassFraction(1.0);
            param->setIndividualDustMassFractions(true);
        }
        else
            param->setDustMassFraction(atof(data.c_str()));
        return true;
    }

    if(cmd.compare("<mrw>") == 0)
    {
        if(atob(atoi(data.c_str())))
            param->setMRW(true);

        cout << "\nWARNING: MRW currently unavailable! " << endl;

        return true;
    }

    if(cmd.compare("<pda>") == 0)
    {
        if(atob(atoi(data.c_str())))
            param->setPDA(true);

        cout << "\nWARNING: PDA currently unavailable! " << endl;

        return true;
    }

    if(cmd.compare("<dust_offset>") == 0)
    {
        if(atob(atoi(data.c_str())))
            param->setDustOffset(true);

        return true;
    }

    if(cmd.compare("<dust_offset min_gas_density = >") == 0)
    {
        if(atob(atoi(data.c_str())))
        {
            string str = seperateString(data);
            formatLine(str);
            dlist min_gas_density = parseValues(str);
            param->setDustOffset(min_gas_density[0]);
        }
        return true;
    }

    if(cmd.compare("<dust_gas_coupling>") == 0)
    {
        if(atob(atoi(data.c_str())))
            param->setDustGasCoupling(true);

        return true;
    }

    if(cmd.compare("<dust_gas_coupling min_gas_density = >") == 0)
    {
        if(atob(atoi(data.c_str())))
        {
            string str = seperateString(data);
            formatLine(str);
            dlist min_gas_density = parseValues(str);
            param->setDustGasCoupling(min_gas_density[0]);
        }
        return true;
    }

    if(cmd.compare("<radiation_field>") == 0)
    {
        if(atob(atoi(data.c_str())))
            param->setSaveRadiationField(true);

        return true;
    }

    if(cmd.compare("<rt_scattering>") == 0)
    {
        if(atob(atoi(data.c_str())))
            param->setScatteringToRay(true);
        else
            param->setScatteringToRay(false);

        return true;
    }

    if(cmd.compare("<full_dust_temp>") == 0)
    {
        if(atob(atoi(data.c_str())))
            param->setFullDustTemp(true);

        return true;
    }

    if(cmd.compare("<stochastic_heating>") == 0)
    {
        double max_size = atof(data.c_str());

        if(max_size >= 0)
        {
            param->setStochasticHeatingMaxSize(max_size);
            return true;
        }
        else
        {
            cout << "\nERROR: For stochastic heating, a non-negative dust grain size "
                    "limit needs to be chosen!"
                 << endl;
            return false;
        }
    }

    if(cmd.compare("<source_dust nr_photons = >") == 0)
    {
        string str = seperateString(data);
        ullong nr_of_photons = ullong(atof(str.c_str()));

        if(nr_of_photons > 0)
            param->setNrOfDustPhotons(long(nr_of_photons));
        else
        {
            cout << "\nERROR: Number of dust photons could not be recognized!" << endl;
            return false;
        }

        return true;
    }

    if(cmd.compare("<source_isrf nr_photons = >") == 0)
    {
        string str = seperateString(data);
        ullong nr_of_photons = ullong(atof(str.c_str()));

        if(nr_of_photons > 0)
            param->setNrOfISRFPhotons(long(nr_of_photons));
        else
        {
            cout << "\nERROR: Number of ISRF photons could not be recognized!" << endl;
            return false;
        }

        str = seperateString(data);
        dlist values = parseValues(data);

        if(str.size() > 0)
            param->setISRF(str);
        else if(values.size() == 1 && values[0] > 0)
            param->setISRF("", values[0]);
        else if(values.size() == 2 && values[0] > 0 && values[1] >= 1)
            param->setISRF("", values[0], values[1]);
        else
        {
            cout << "\nERROR: Path of ISRF SED path could not be recognized!" << endl;
            return false;
        }

        return true;
    }

    if(cmd.compare("<foreground_extinction>") == 0)
    {
        dlist values = parseValues(data);

        if(values.size() == 1)
            param->setForegroundExtinction(values[0]);
        else if(values.size() == 2)
            param->setForegroundExtinction(values[0], values[1]);
        else if(values.size() == 3)
            param->setForegroundExtinction(values[0], values[1], values[2]);
        else
        {
            cout << "\nERROR: Foreground extinction command could not be recognized!" << endl;
            return false;
        }

        return true;
    }

    if(cmd.compare("<enfsca>") == 0)
    {
        if(atob(atoi(data.c_str())))
            param->setEnfScattering(true);
        else
            param->setEnfScattering(false);

        return true;
    }

    if(cmd.compare("<peel_off>") == 0)
    {
        if(atob(atoi(data.c_str())))
            param->setPeelOff(true);
        else
            param->setPeelOff(false);

        return true;
    }

    if(cmd.compare("<acceptance_angle>") == 0)
    {
        double angle = atof(data.c_str());
        if(angle > 0)
            param->setAcceptanceAngle(angle);
        return true;
    }

    if(cmd.compare("<nr_threads>") == 0)
    {
        int max_t = omp_get_max_threads();
        int tr = int(atof(data.c_str()));

        if(tr == -1)
            tr = max_t;

        if(tr > max_t)
        {
            tr = max_t;
            cout << "\nWARNING: Max. nr. of threads is:  " << max_t << endl;
        }

        if(tr <= 0)
        {
            tr = 1;
            cout << "\nWARNING: Max. nr. of threads is set to: " << 1 << endl;
        }

        param->setNrOfThreads(tr);
        return true;
    }

    if(cmd.compare("<vel_is_speed_of_sound>") == 0)
    {
        if(atob(atoi(data.c_str())))
            param->setIsSpeedOfSound(true);

        return true;
    }

    if(cmd.compare("<amira_inp_points>") == 0)
    {
        uint points = uint(atof(data.c_str()));

        // if(points % 2 != 0)
        //    points++;

        param->setInpAMIRAPoints(points);
        return true;
    }

    if(cmd.compare("<amira_out_points>") == 0)
    {
        uint points = uint(atof(data.c_str()));

        // if(points % 2 != 0)
        //    points++;

        param->setOutAMIRAPoints(points);
        return true;
    }

    if(cmd.compare("<plot_inp_midplanes>") == 0)
    {
        bool plot = atob(atoi(data.c_str()));
        param->setInpMidPlot(plot);
        return true;
    }

    if(cmd.compare("<plot_out_midplanes>") == 0)
    {
        bool plot = atob(atoi(data.c_str()));
        param->setOutMidPlot(plot);
        return true;
    }

    if(cmd.compare("<write_3d_midplanes>") == 0)
    {
        dlist values = parseValues(data);

        uint plane, nr_of_slices = 0;
        double z_min = 0, z_max = 0;
        if(values.size() >= 1)
            plane = uint(values[0]);
        if(values.size() >= 2)
            nr_of_slices = uint(values[1]);
        if(values.size() == 4)
        {
            if(values[2] > values[3])
            {
                cout << "\nERROR: z_min is larger than z_max for 3D midplane creation!" << endl;
                return false;
            }
            z_min = values[2];
            z_max = values[3];
        }
        if(values.size() > 4 || values.size() == 0)
        {
            cout << "\nERROR: Number of parameters for 3D midplane files in line " << line_counter
                 << "could not be recognized!" << endl;
            return false;
        }

        if(plane > 3 || plane < 0)
        {
            cout << "\nERROR: Wrong plane for 3D midplane files in line " << line_counter << "!" << endl;
            return false;
        }

        param->set3dMidplane(plane, nr_of_slices, z_min, z_max);
        return true;
    }

    if(cmd.compare("<write_inp_midplanes>") == 0)
    {
        uint points = uint(atof(data.c_str()));

        // if(points % 2 != 0)
        //    points++;

        param->setInpMidDataPoints(points);
        return true;
    }

    if(cmd.compare("<write_out_midplanes>") == 0)
    {
        uint points = uint(atof(data.c_str()));

        param->setOutMidDataPoints(points);
        return true;
    }

    if(cmd.compare("<write_radiation_field>") == 0)
    {
        uint val = atoi(data.c_str());

        if(val > 3)
        {
            cout << "\nWARNING: Command \"<write_radiation_field>\" accepts only paramers between 0 to 3!"
                 << endl;
            param->setWriteRadiationField(0);
            return true;
        }

        param->setWriteRadiationField(val);

        return true;
    }

    if(cmd.compare("<write_full_radiation_field>") == 0)
    {
        cout << "\nWARNING: Command <write_full_radiation_field> is no longer available!" << endl;

        return true;
    }

    if(cmd.compare("<write_g_zero>") == 0)
    {
        if(atob(atoi(data.c_str())))
            param->setWriteGZero(true);
        else
            param->setWriteGZero(false);

        return true;
    }

    if(cmd.compare("<midplane_zoom>") == 0)
    {
        param->setMidplaneZoom(uint(atof(data.c_str())));
        return true;
    }

    if(cmd.compare("<kepler_star_mass>") == 0)
    {
        param->setKeplerStarMass(atof(data.c_str()));
        return true;
    }

    if(cmd.compare("<turbulent_velocity>") == 0)
    {
        param->setTurbulentVelocity(atof(data.c_str()));
        return true;
    }

    return false;
}

bool CCommandParser::checkPixel(dlist & values, dlist nr_of_pixel, bool nsides_as_pixel)
{
    if(nsides_as_pixel)
    {
        if(nr_of_pixel.size() != 1)
        {
            cout << "\nERROR: Number of sides in healpix "
                 << "line detector could not be recognized!" << endl;
            return false;
        }
        if(nr_of_pixel[0] <= 0)
        {
            cout << "\nERROR: Number of sides in healpix line detector could not be "
                    "recognized!"
                 << endl;
            return false;
        }
        if(!CMathFunctions::isPowerOfTwo(int(nr_of_pixel[0])))
        {
            cout << "\nERROR: Number of sides must be a power of two!" << endl;
            return false;
        }
        values.push_back(uint(nr_of_pixel[0]));
        values.push_back(uint(nr_of_pixel[0]));
    }
    else
    {
        if(nr_of_pixel.size() == 1)
        {
            // Pixel per axis
            values.push_back(uint(nr_of_pixel[0]));
            values.push_back(uint(nr_of_pixel[0]));
        }
        else if(nr_of_pixel.size() == 2 && !nsides_as_pixel)
        {
            // Individual pixel per axi
            values.push_back(uint(nr_of_pixel[0]));
            values.push_back(uint(nr_of_pixel[1]));
        }
        else
        {
            cout << "\nERROR: Number of pixel in Raytracing detector could not be "
                    "recognized!"
                 << endl;
            return false;
        }
        for(int i = 0; i < nr_of_pixel.size(); i++)
        {
            if(nr_of_pixel[i] <= 0)
            {
                cout << "\nERROR: Number of pixel in Raytracing detector could not be "
                        "recognized!"
                     << endl;
                return false;
            }
        }
    }
    return true;
}

bool CCommandParser::checkVelChannels(dlist & values, dlist nr_of_channels)
{
    if(nr_of_channels.size() == 1)
    {
        if(nr_of_channels[0] <= 0)
        {
            cout << "\nERROR: Number of velocity channels in Raytracing detector could "
                    "not be recognized!"
                 << endl;
            return false;
        }
        values.push_back(uint(nr_of_channels[0]));
    }
    else
    {
        cout << "\nERROR: Number of velocity channels in Raytracing detector could not "
                "be recognized!"
             << endl;
        return false;
    }
    return true;
}
