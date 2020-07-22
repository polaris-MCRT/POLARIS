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
        
        has_emi_data=false;
        has_abs_data=false;
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

        if(entries.size()>0)
        {
            for(uint i =0;i<entries.size();i++)
            {
                COpiateEntry * tmp=entries[i];
                delete tmp;
                tmp=0;
                entries[i]=0;
            }
            
            entries.clear();
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
    
    double getProjCellVelocityInterp(const Vector3D & tmp_pos,
                                     const Vector3D & dir_map_xyz,
                                     const VelFieldInterp & vel_field_interp)
    {
        double cell_velocity = 0;

        // Get the velocity in the photon direction of the current position
        if(vel_field_interp.vel_field.size() > 0 && !vel_field_interp.zero_vel_field)
        {
            // Get velocity from grid cell with interpolation
            Vector3D rel_pos = tmp_pos - vel_field_interp.start_pos;
            cell_velocity = vel_field_interp.vel_field.getValue(rel_pos.length());
        }
        return cell_velocity;
    }
    
    uint getIndex(uint op_id) const
    {
        uint N = max_ids;
        uint min = 0, max = N - 1;

        if(op_id < list_IDs[0] || op_id > list_IDs[max])
            return MAX_UINT;

        while(max - min > 1)
        {
            uint i = min + (max - min) / 2;
            if(list_IDs[i] >= op_id)
                max = i;
            else
                min = i;
        }

        if(list_IDs[min] == op_id)
            return min;

        uint lower = min - 1;
        uint upper = min + 1;

        if(lower == MAX_UINT)
            lower = 0;

        if(upper >= max_ids)
            upper = max_ids - 1;

        if(list_IDs[lower] == op_id)
            return lower;

        if(list_IDs[upper] == op_id)
            return upper;

        return min;
    }
    
    
    void getMatrices(CGridBasic * grid,
                                 const photon_package * pp,
                                 uint i_spec,
                                 uint i_trans,
                                 double velocity,
                                 const LineBroadening & line_broadening,
                                 const MagFieldInfo & mag_field_info,
                                 StokesVector * line_emissivity,
                                 Matrix2D * line_absorption_matrix) const
    {
        double emission = 0.0;
        double absorption = 0.0;
        
        uint op_id=grid->getOpiateID(pp);
        uint index = getIndex(op_id);
        
        // Reset absorption matrix and emissivity
        line_absorption_matrix->resize(4, 4);
        *line_emissivity = 0;
              
        if(index!=MAX_UINT)
        {
            if(has_abs_data==true)
                absorption = mat_absorption(index,current_index)*line_broadening.gauss_a;

            if(has_emi_data==true)
                emission = mat_emissivity(index,current_index)*line_broadening.gauss_a;
        }

        // Calculate the line matrix from rotation polarization matrix and line shape
        // getGaussLineMatrix(grid, pp, velocity, line_absorption_matrix);
        double line_amplitude = exp(-(pow(velocity, 2) * pow(line_broadening.gauss_a, 2))) / PIsq;
        
        // Only diagonal without polarization rotation matrix elements
        for(uint i = 0; i < 4; i++)
        {
            line_absorption_matrix->setValue(i, i, line_amplitude);
        }

        // Calculate the Emissivity of the gas particles in the current cell
        *line_emissivity = *line_absorption_matrix * StokesVector(emission, 0, 0, 0);

        // Calculate the line matrix from rotation polarization matrix and line shape
        *line_absorption_matrix *= absorption;
    }
    
    uint getMaxSpecies()
    {
        return max_species;
    }
    
    double getFrequency(uint pos)
    {
        return list_freq[pos];
    }
    
    
    bool readDataBase(string filename);
    
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
    
    string getCurrentName()
    {
        return list_names[current_index];
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

    double getCurrentFrequency()
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

  private:
    class COpiateEntry
    {
        public:
            COpiateEntry()
            {
                freq=0;
                weight=0;
                size=0;
                name="Empty";
                em=0;
                ex=0;
            }

            COpiateEntry(uint _size)
            {
                freq=0;
                weight=0;
                size=_size;
                name="Empty";
                em=new double[size];
                ex=new double[size];

                for(uint i=0;i<size;i++)
                {
                    em[i]=0;
                    ex[i]=0;
                }               
            }

            ~COpiateEntry()
            {
                if(em != 0)
                {
                    delete[] em;
                    em = 0;
                }

                if(ex != 0)
                {
                    delete[] ex;
                    ex = 0;
                }
            }

            void setData(uint pos, double _em, double _ex)
            {
                em[pos]=_em;
                ex[pos]=_ex;
            }

            double freq;
            double weight;
            uint size;
            string name;
            double * em;
            double * ex;
    };    

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
    
    bool has_emi_data;
    bool has_abs_data;
    
    vector < COpiateEntry * > entries;
};

#endif
