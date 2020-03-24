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
/* History:   05.12.2016                                                    */
/****************************************************************************/

#include <CCfits/CCfits>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <utility>
#include <vector>

#include "CommandParser.h"
#include "Grid.h"
#include "Parameters.h"

using namespace std;
using namespace CCfits;

#ifndef OPIATE_H
#define OPIATE_H

class COpiateDataBase
{
  public:
    COpiateDataBase()
    {
        current_index = MAX_UINT;
        list_freq = 0;
        list_IDs = 0;
        list_weight = 0;
        max_species = 0;
        max_ids = 0;

        database_counter = 0;

        max_velocity = 1;
        velocity_channel = 0;

        path_emi="";
        path_abs="";
    }

    ~COpiateDataBase()
    {
        cout << CLR_LINE;
        cout << "Final OPIATE data cleanup ...   \r" << flush;

        if(list_freq != 0)
        {
            delete[] list_freq;
            list_freq = 0;
        }

        if(list_IDs != 0)
        {
            delete[] list_IDs;
            list_IDs = 0;
        }

        if(list_weight != 0)
        {
            delete[] list_weight;
            list_weight = 0;
        }

        if(velocity_channel != 0)
        {
            delete[] velocity_channel;
            velocity_channel = 0;
        }

        cout << CLR_LINE;
    }

    double getGaussA(double temp_gas, double v_turb)
    {
        double v_th = sqrt(2.0 * con_kB * temp_gas / (list_weight[current_index] * 1e-3 / con_Na));
        double gauss_a = 1.0 / sqrt(pow(v_th, 2) + pow(v_turb, 2));
        return gauss_a;
    }

    void calcLineBroadening(CGridBasic * grid)
    {
        long max_cells = grid->getMaxDataCells();

        cout << CLR_LINE;
        cout << "-> Calculating line broadening for each grid cell ...     \r" << flush;

    #pragma omp parallel for schedule(dynamic)
        for(long i_cell = 0; i_cell < long(max_cells); i_cell++)
        {
            cell_basic * cell = grid->getCellFromIndex(i_cell);

            // Get necessary quantities from the current cell
            double temp_gas = grid->getGasTemperature(*cell);
            double turbulent_velocity = grid->getTurbulentVelocity(cell);

            // Set gauss_a for each transition only once
            grid->setGaussA(cell, getGaussA(temp_gas, turbulent_velocity));
        }

        cout << CLR_LINE;
    }

    uint getMaxSpecies()
    {
        return max_species;
    }

    double getFrequency(uint pos)
    {
        return list_freq[pos];
    }


    bool readOpiateDataBase(parameters & param);

    void printParameters(parameters & param, CGridBasic * grid);

    bool findIndexByName(string name)
    {
        for(uint i = 0; i < list_names.size(); i++)
        {
            if(name.compare(list_names[i]) == 0)
            {
                current_index = i;
                return true;
            }
        }

        cout << CLR_LINE;
        cout << "ERROR: A species by the name of \"" << name
             << "\" is not listed in the loaded OPIATE databases!              \n";
        return false;
    }

    bool readEmissivityData(string filename)
    {
        path_emi=filename;
        return readFitsData(filename, mat_emissivity);
    }

    bool readAbsorptionData(string filename)
    {
        path_abs=filename;
        return readFitsData(filename, mat_absorption);
    }

    bool readFitsData(string filename, Matrix2D & mat);

    double getFrequency()
    {
        return list_freq[current_index];
    }

    double getMolecularWeight()
    {
        return list_weight[current_index];
    }

    bool initVelChannels(uint nr_of_channels, double max_vel)
    {
        if(velocity_channel != 0)
        {
            delete[] velocity_channel;
            velocity_channel = 0;
        }

        nr_of_velocity_channels = nr_of_channels;
        max_velocity = max_vel;

        velocity_channel = new double[nr_of_velocity_channels];

        if(nr_of_velocity_channels > 1)
        {
            for(uint i = 0; i < nr_of_velocity_channels; i++)
            {
                velocity_channel[i] =
                    2 * (float)i / ((float)nr_of_velocity_channels - 1) * max_velocity - max_velocity;
            };
        }
        else if(nr_of_velocity_channels == 1)
        {
            velocity_channel[0] = 0;
        }
        else
        {
            cout << CLR_LINE;
            cout << "ERROR: Number of velocity channels is not larger than zero!                 \n";
            return false;
        }

        return true;
    }

    double getVelocityChannel(uint vch)
    {
        return velocity_channel[vch];
    }

    uint getNrOfVelChannels()
    {
        return nr_of_velocity_channels;
    }

    double getMaxVelocity()
    {
        return max_velocity;
    }

    uint biListIndexDataSearch(uint id)
    {
        uint N = max_ids;
        uint min = 0, max = N - 1;

        // cout << list_IDs[0] << "\n";
        // cout << list_IDs[N-1] << "\n";

        if(N == 0)
            return MAX_UINT;

        if(id < list_IDs[0] || id > list_IDs[N - 1])
            return MAX_UINT;

        while(max - min > 1)
        {
            uint i = min + (max - min) / 2;
            if(list_IDs[i] >= id)
                max = i;
            else
                min = i;
        }

        if(list_IDs[min] == id)
            return min;

        uint lower = min - 1;
        uint upper = min + 1;

        if(lower == MAX_UINT)
            lower = 0;

        if(upper >= N)
            upper = N - 1;

        // check neighboring indices just to be sure
        if(list_IDs[lower] == id)
            return lower;

        if(list_IDs[upper] == id)
            return upper;

        return MAX_UINT;
    } /**/

  private:
    // void formatLine(string &line);

    /*dlist parseValues(string & str)
    {
        uint pos;
        dlist values;
        string v;

        if(str.size() == 0)
            return values;

        formatLine(str);

        while(str.find(" ") != string::npos)
        {
            pos = uint(str.find(" "));
            v = str.substr(0, pos);
            str.erase(0, pos + 1);
            values.push_back(atof(v.c_str()));
        }

        values.push_back(atof(str.c_str()));

        return values;
    }/**/

    Matrix2D mat_emissivity;
    Matrix2D mat_absorption;

    uint current_index;
    uint max_species;
    uint max_ids;

    string path_emi;
    string path_abs;

    double * list_freq;
    double * list_IDs;
    double * list_weight;
    strlist list_names;

    double * velocity_channel;
    uint nr_of_velocity_channels;
    double max_velocity;

    uint database_counter;
};

#endif
