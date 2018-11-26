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

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <stdlib.h> 
#include <algorithm>
#include <fstream>

using namespace std;

#ifndef OPIATE_H
#define OPIATE_H

typedef unsigned int uint;
#define MAX_UINT uint(-1)
#define CLR_LINE1 "                                                          \r"

class opiate_link
{
public:

    opiate_link()
    {
        value = 1e-300;
        list = 0;
    }

    opiate_link(double val)
    {
        value = val;
        list = 0;
    }

    double value;
    vector<opiate_link> * list;
};

typedef vector<opiate_link> linked_list;

typedef vector<double> dlist;

class opiate_data
{
public:

    opiate_data()
    {
        id = MAX_UINT;
        data = 0;
    }

    opiate_data(dlist values)
    {
        id = uint(values[0]);
        uint tmp_len = uint(values.size());
        data = new double[tmp_len - 1];

        for(uint i = 1; i < tmp_len; i++)
            data[i - 1] = values[i];
    }

    ~opiate_data()
    {
        if(data != 0)
            delete[] data;

        data = 0;
    }

    uint id;
    double * data;
};

class opiate_ids
{
public:

    opiate_ids()
    {
        N = MAX_UINT;
        uiniqD = MAX_UINT;

        dx = 1e-300;
        dens = 1e-300;
        temp = 1e-300;
        flge = 1e-300;
        fluv = 1e-300;
        flih = 1e-300;
        fli2 = 1e-300;
    }
    //UniqID	dx	dens	temp	flge	fluv	flih	fli2	N

    opiate_ids(uint _N, uint _uiniqD, double _dx, double _dens,
            double _temp, double _flge, double _fluv, double _flih, double _fli2)
    {
        N = _N;
        uiniqD = _uiniqD;
        dx = _dx;
        dens = _dens;
        temp = _temp;
        flge = _flge;
        fluv = _fluv;
        flih = _flih;
        fli2 = _fli2;
    }

    opiate_ids(double _dx, double _dens,
            double _temp, double _flge, double _fluv, double _flih, double _fli2)
    {
        N = MAX_UINT;
        uiniqD = MAX_UINT;
        dx = _dx;
        dens = _dens;
        temp = _temp;
        flge = _flge;
        fluv = _fluv;
        flih = _flih;
        fli2 = _fli2;
    }

    ~opiate_ids()
    {
    };

    uint N;
    uint uiniqD;
    double dx;
    double dens;
    double temp;
    double flge;
    double fluv;
    double flih;
    double fli2;
};

typedef vector<opiate_data*> opiate_data_list;

class COpiate
{
public:

    COpiate()
    {
        ru[0] = '|';
        ru[1] = '/';
        ru[2] = '-';
        ru[3] = '\\';

        data_len = 0;
        data_rows = 0;
        max_entries = 0;
        column_length = 0;
        velocity_channel = 0;

        nr_of_velocity_channels = 0;
        max_velocity = 0;
        
        abundance = 0;
        mol_weight = 0;
    }

    ~COpiate()
    {
        cout << CLR_LINE1;
        cout << "Final OPIATE cleanup ...   \r" << flush;
        if(data_list.size() != 0)
        {
            for(uint i = 0; i < data_list.size(); i++)
            {
                opiate_data *tmp_data = data_list[i];

                delete tmp_data;
                tmp_data = 0;
            }

            data_list.clear();
        }

        if(velocity_channel != 0)
        {
            delete[] velocity_channel;
            velocity_channel = 0;
        }


        if(id_list.size() != 0)
        {

            clearList(&id_list);
        }

        cout << CLR_LINE1;
    }
    
    
    void setAbundance(double a)
    {
        abundance = a;
    }
    
    void setMolWeight(double w)
    {
        mol_weight=w;
    }
    
    double getAbundance()
    {
        return abundance;
    }
    
    double getMolWeight()
    {
        return mol_weight;
    }

    linked_list * insertInNextLevelIDList(linked_list * list, double val)
    {
        linked_list & tmp_list = *list;
        uint N = uint(tmp_list.size());
        uint min = 0, max = N - 1;
        uint lower, upper;
        linked_list * res_list = 0;

        if(N == 0)
        {
            opiate_link link;
            link.value = val;
            tmp_list.push_back(link);

            if(tmp_list[0].list == 0)
                tmp_list[0].list = new linked_list;

            return tmp_list[0].list;
        }

        if(val < tmp_list[0].value)
        {
            linked_list::iterator it = tmp_list.begin();
            opiate_link link;
            link.value = val;
            link.list = new linked_list;

            tmp_list.insert(it, link);
            return link.list;
        }

        if(val > tmp_list[max].value)
        {
            opiate_link link;
            link.value = val;
            link.list = new linked_list;

            tmp_list.push_back(link);

            return link.list;
        }

        while(max - min > 1)
        {
            uint i = min + (max - min) / 2;
            if(tmp_list[i].value >= val)
                max = i;
            else
                min = i;
        }

        if(tmp_list[min].value == val)
            return tmp_list[min].list;

        lower = min - 1;
        upper = min + 1;

        if(lower == MAX_UINT)
            lower = 0;

        if(upper > N)
            upper = N;

        if(tmp_list[lower].value == val)
            return tmp_list[lower].list;

        if(tmp_list[upper].value == val)
            return tmp_list[upper].list;

        opiate_link link;
        linked_list::iterator it = tmp_list.begin();
        link.value = val;
        link.list = new linked_list;

        tmp_list.insert(it + min + 1, link);

        return link.list;
    }

    opiate_ids findID(double dx, double dens, double temp, double flge, double fluv, double flih, double fli2)
    {
        opiate_ids res;

        linked_list * next_level;

        next_level = findNextLevelList(&id_list, fli2);

        if(next_level == 0)
            return res;

        res.fli2 = fli2;

        next_level = findNextLevelList(next_level, flih);

        if(next_level == 0)
            return res;

        res.flih = flih;

        next_level = findNextLevelList(next_level, fluv);

        if(next_level == 0)
            return res;

        res.fluv = fluv;

        next_level = findNextLevelList(next_level, flge);

        if(next_level == 0)
            return res;

        res.flge = flge;

        next_level = findNextLevelList(next_level, temp);

        if(next_level == 0)
            return res;

        res.temp = temp;

        next_level = findNextLevelList(next_level, dens);

        if(next_level == 0)
            return res;

        res.dens = dens;

        next_level = findNextLevelList(next_level, dx);

        if(next_level == 0)
            return res;

        res.dx = dx;

        if(next_level->size() < 2)
        {
            cout << "\nERROR: Linked list is inconsistent!  " << endl;
            return res;
        }


        res.uiniqD = uint((*next_level)[0].value);
        res.N = uint((*next_level)[1].value);

        return res;
    }

    bool loadUniqueParamFile(string filename);
    bool loadDataFile(string filename);

    double getData(uint id, uint pos)
    {

        if(column_length == 0)
        {
            cout << "\nERROR: No data file loaded yet!" << endl;
            return -1;
        }

        uint index = biListIndexDataSearch(id);

        if(index == MAX_UINT)
        {
            cout << "\nERROR: ID out of boundaries!" << endl;
            return -1;
        }

        if(pos >= column_length)
        {
            cout << "\nERROR: Position out of boundaries!" << endl;
            return -1;
        }

        opiate_data *data = data_list[index];
        double res = data->data[pos];

        return res;
    }

    bool initVelChannels(uint nr_of_channels, double max_vel)
    {
        nr_of_velocity_channels = nr_of_channels;
        max_velocity = max_vel;
		
        velocity_channel = new double[nr_of_velocity_channels];

        if(nr_of_velocity_channels > 1)
        {
            for(uint i = 0; i < nr_of_velocity_channels; i++)
            {
                velocity_channel[i] = 2 * (float) i
                        / ((float) nr_of_velocity_channels - 1) * max_velocity
                        - max_velocity;
            };
        }
        else if(nr_of_velocity_channels == 1)
        {
            velocity_channel[0] = 0;
        }
        else
        {
            cout
                    << "Number of velocity channels is not larger than zero!          " << endl;
            return false;
        }

		return true;
    }
    
    
    double getVelocityChannel(uint vch)
    { 
        return velocity_channel[vch];
    }
    
    uint getNrOfVelChannels(uint vch)
    { 
        return nr_of_velocity_channels;
    }
    
    double getMaxVelocity()
    {
        return max_velocity;
    }

	linked_list id_list; //top level linked list fli2 flih fluv flge temp dens dx

private:
    unsigned char ru[4];
    uint data_len;
    uint data_rows;
    uint max_entries;
    uint column_length;
    opiate_data_list data_list;
    
    double * velocity_channel;
    uint nr_of_velocity_channels;
    double max_velocity;
    double abundance;
    double mol_weight;

    void formatLine(string &line);

    uint insertInList(opiate_data * data)
    {
        uint N = uint(data_list.size());
        uint id = data->id;

        if(N == 0)
        {
            data_list.push_back(data);
            return 0;
        }

        if(id < data_list[0]->id)
        {
            opiate_data_list::iterator it = data_list.begin();
            data_list.insert(it, data);
            return 0;
        }

        if(id > data_list[N - 1]->id)
        {
            data_list.push_back(data);
            return uint(data_list.size() - 1);
        }

        opiate_data_list::iterator it = data_list.begin();
        uint index = biListIndexDataSearch(id);

        if(index == MAX_UINT)
            return MAX_UINT;

        if(data_list[index + 1]->id == id)
            index++;

        data_list.insert(it + index + 1, data);

        return index + 1;
    }

    uint biListIndexDataSearch(uint id)
    {
        uint N = uint(data_list.size());
        uint min = 0, max = N - 1;

        if(id < data_list[0]->id || id > data_list[max]->id)
            return MAX_UINT;

        while(max - min > 1)
        {
            uint i = min + (max - min) / 2;
            if(data_list[i]->id >= id)
                max = i;
            else
                min = i;
        }

        if(data_list[min]->id == id)
            return min;

        uint lower = min - 1;
        uint upper = min + 1;

        if(lower == MAX_UINT)
            lower = 0;

        if(upper >= data_list.size())
            upper = uint(data_list.size()) - 1;

        if(data_list[lower]->id == id)
            return lower;

        if(data_list[upper]->id == id)
            return upper;

        return min;
    }

    dlist parseValues(string & str)
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
    }

    void clearList(linked_list * list)
    {
        linked_list & tmp_list = *list;

        for(uint i = 0; i < tmp_list.size(); i++)
        {
            if(tmp_list[i].list == 0)
                continue;

            clearList(tmp_list[i].list);
            tmp_list[i].list->clear();
            delete tmp_list[i].list;
            tmp_list[i].list = 0;
        }

        list->clear();
    }

    linked_list * findNextLevelList(linked_list * list, double val)
    {
        linked_list & tmp_list = *list;
        uint N = uint(tmp_list.size());
        uint min = 0, max = N - 1;
        uint lower, upper;

        if(N == 0)
            return 0;

        if(val < tmp_list[0].value)
            return 0;

        if(val > tmp_list[max].value)
            return 0;

        while(max - min > 1)
        {
            uint i = min + (max - min) / 2;
            if(tmp_list[i].value >= val)
                max = i;
            else
                min = i;
        }

        if(tmp_list[min].value == val)
            return tmp_list[min].list;

        lower = min - 1;
        upper = min + 1;

        if(lower == MAX_UINT)
            lower = 0;

        if(upper > N)
            upper = N;

        if(tmp_list[upper].value == val)
            return tmp_list[upper].list;

        if(tmp_list[lower].value == val)
            return tmp_list[lower].list;

        return 0;
    }

    //fli2 flih fluv flge temp dens dx

    void insertInIDList(uint uniqID, double dx, double dens, double temp, double flge, double fluv, double flih, double fli2, uint N)
    {
        linked_list * next_level;

        next_level = insertInNextLevelIDList(&id_list, fli2);
        next_level = insertInNextLevelIDList(next_level, flih);
        next_level = insertInNextLevelIDList(next_level, fluv);
        next_level = insertInNextLevelIDList(next_level, flge);
        next_level = insertInNextLevelIDList(next_level, temp);
        next_level = insertInNextLevelIDList(next_level, dens);
        next_level = insertInNextLevelIDList(next_level, dx);
        next_level->push_back(opiate_link(double(uniqID)));
        next_level->push_back(opiate_link(double(N)));
    }


};

#endif

