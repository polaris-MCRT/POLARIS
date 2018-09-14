/****************************************************************************/
/* Code:       POLARIS_v2                                                   */
/*                                                                          */
/* Author:     -Stefan Reissl                                               */
/*              reissl@uni-heidelberg.de                                    */
/*                                                                          */
/*              University of Heidelberg,                                   */
/*              Institute of Theoretical Astrophysics (ITA),                */
/*              Albert-Ueberle-Str. 2, 69120 Heidelberg, Germany            */
/*                                                                          */
/*                                                                          */
/* COPYRIGHT:                                                               */
/* The code is free of charge for any scientific purpose.                   */
/* This software is provided in the hope that it will be useful but         */
/* without any warranty of ability or fitness of a particular purpose.      */
/* We also reject any responsibility for incorrect results that may be      */
/* produced with this code.                                                 */
/* Any publication that makes use of the POLARIS software                   */
/* (completely or in part) must mention the name of the POLARIS code and    */
/* cite the paper Reissl et al. 2016                                        */
/* http://esoads.eso.org/abs/2016A%26A...593A..87R                          */
/*                                                                          */
/* History:   29.11.2016                                                    */
/****************************************************************************/

#include "OPIATE.h"

/*void OPIATE::formatLine(string &line)
{
    string::size_type pos = 0;

    if(line.find_first_of("#") != string::npos)
    {
        pos = line.find("#");
        line.erase(pos, line.length() - pos - 1);
    }

    if(line.find_first_of("!") != string::npos)
    {
        pos = line.find("!");
        line.erase(pos, line.length() - pos - 1);
    }

    if(line.size() < 3)
    {
        line = "";
        return;
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

    while(line.find("  ") != string::npos)
    {
        pos = line.find("  ");
        line.replace(pos, 2, " ");
    }

    if(line == " ")
        line = "";

    while(line.find(",") != string::npos)
    {
        pos = line.find(",");
        line.replace(pos, 1, ".");
    }

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
        line.erase(0, 1);
    }
}

bool OPIATE::loadDataFile(string filename)
{
    uint line_counter = 0;
    uint char_counter = 0;
    string line;

    ifstream reader(filename.c_str());

    cout << CLR_LINE1;

    if(reader.fail())
    {
        cout << "ERROR: Cannot load OPIAE data file:" << endl;
        cout << filename << "\n" << endl;
        return false;
    }

    while(getline(reader, line))
    {
        if(line_counter % 1000 == 0)
        {
            char_counter++;
            cout << " -> Loading OPIATE data file : "
                    << ru[(unsigned int) char_counter % 4] << "           \r" << flush;
        }

        formatLine(line);

        if(line.size() == 0)
            continue;

        dlist values = parseValues(line);

        if(line_counter == 0)
        {
            if(values.size() < 2)
            {
                cout << "ERROR in line: " << line_counter << endl;
                cout << "Opiate data files rquire more than one column! \n" << endl;
                return false;
            }

            max_entries = uint(values.size());
            column_length = max_entries - 1;
        }
        else
        {
            if(values.size() != max_entries)
            {
                cout << "ERROR in line: " << line_counter << endl;
                cout << " Wrong ammount of columns. \n" << endl;
                return false;
            }
        }

        opiate_data * data = new opiate_data(values);

        insertInList(data);

        line_counter++;
    }

    reader.close();
    cout << CLR_LINE1;
    cout << " Loading OPIATE data file              : done" << endl;
    return true;
}

bool OPIATE::loadUniqueParamFile(string filename)
{
    uint line_counter = 0;
    uint char_counter = 0;
    string line;

    ifstream reader(filename.c_str());

    if(reader.fail())
    {
        cout << "ERROR: Cannot load uinique parameter file:" << endl;
        cout << filename << "\n" << endl;
        return false;
    }

    while(getline(reader, line))
    {
        line_counter++;

        if(line_counter % 1000 == 0)
        {
            char_counter++;
            cout << "-> Loading uinique parameter file : "
                    << ru[(unsigned int) char_counter % 4] << "           \r" << flush;
        }

        formatLine(line);

        if(line.size() == 0)
            continue;

        dlist values = parseValues(line);

        if(values.size() != 9)
        {
            cout << "ERROR: in line: " << line_counter << endl;
            cout << " Wrong ammount of columns. \n" << endl;
            return false;
        }

        uint UniqID = uint(values[0]);
        double dx = values[1];
        double dens = values[2];
        double temp = values[3];
        double flge = values[4];
        double fluv = values[5];
        double flih = values[6];
        double fli2 = values[7];
        uint N = uint(values[8]);

        insertInIDList(UniqID, dx, dens, temp, flge, fluv, flih, fli2, N);
    }

    reader.close();

    cout << CLR_LINE1;
    cout << " Loading OPIATE uinique parameter file : done" << endl;
    return true;
}*/

