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

#include <stdio.h>
#include <string.h>
#include <memory>
#include <valarray>

#include "CCfits/CCfits.h"
#include "CCfits/FITS.h"
#include "CCfits/FITSUtilT.h"
#include "CCfits/HDU.h"
#include "CCfits/PHDU.h"
#include "CCfits/PHDUT.h"
#include "Parameters.h"

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
}*/

void COpiateDataBase::printParameters(parameters & param, CGridBasic * grid)
{
    cout << CLR_LINE;
    cout << "OPIATE parameters                             " << endl;
    cout << SEP_LINE;

    if(database_counter==0)
    {
        cout << "\nERROR: No OPIATE database available!               \n" ;
        return;
    }

    if(database_counter==1)
    {
        cout << "- Emission   database:\n   " << path_emi << endl;
        cout << "- Absorption database: none                 \n";
    }

    if(database_counter==2)
    {
        cout << "- Emission   database:\n   " << path_emi << endl;
        cout << "- Absorption database:\n   " << path_abs << endl;
    }

    cout << "- Velocity field     : ";

    if(grid->isVelocityFieldAvailable())
        cout << "taken from grid" << endl;
    else
        cout << "none" << endl;

    cout << "- Turbulent Velocity : ";
    if(param.getTurbulentVelocity() > 0)
        cout << param.getTurbulentVelocity() << " [m/s]" << endl;
    else if(grid->isTurbulentVelocityAvailable())
        cout << "taken from grid" << endl;
    else
        cout << "none" << endl;

    cout << "\n- Available species:    " << endl;
    for(uint i = 0; i<max_species; i++)
    {
        cout << "  - Species nr. " << i +1 << ", label: " <<  list_names[i];

        if(list_weight[i]>0)
            cout << ", weight: " << list_weight[i];

        double freq=list_freq[i];

        if(freq/1.0e9>1.0)
        {
            cout << ", freq.: " << freq/1.0e9 << " GHz" << endl;
        }
        else if(freq/1.0e6>1.0)
        {
            cout << ", freq.: " << freq/1.0e6 << " MHz" << endl;
        }
        else if(freq/1.0e3>1.0)
        {
            cout << ",freq.: " << freq/1.0e6 << " kHz" << endl;
        }
        else
        {
            cout << "   freq. : " << freq/1.0e6 << " Hz" << endl;
        }
    }

    cout << "\n- Selected species :    " << endl;

    for(uint i=0;i<param.getNrOfOPIATESpecies();i++)
    {
        cout << "  - Detector nr. " << i+1 << ", label: " << param.getOpiateSpec(i) << endl;
    }

    cout << SEP_LINE;
}


bool COpiateDataBase::readOpiateDataBase(parameters & param)
{
    string path_emi=param.getOpiatePathEmission();
    string path_abs=param.getOpiatePathAbsorption();


    if(path_emi.size()>0)
    {
        if(!readEmissivityData(path_emi))
            return false;
    }
    else
    {
        cout << CLR_LINE;
        cout << "\nERROR: A path to an OPIATE emissivity database is required!               \n" ;
        return false;
    }

    if(path_abs.size()>0)
    {
        if(!readAbsorptionData(path_abs))
            return false;
    }

    return true;
}

bool COpiateDataBase::readFitsData(string filename, Matrix2D & mat)
{
        //
        //cout << "Reading OPIATE fits data from:\n       " << filename <<  "               \n" << flush;
        // auto_ptr<FITS> pInfile(0);
        unique_ptr<FITS> pInfile;

        cout << CLR_LINE;
        cout << "-> Reading fits data...              \r" << flush;

        try
        {
            pInfile.reset(new FITS(filename.c_str(),Read,true));
        }
        catch(CCfits::FITS::CantOpen)
        {
            cout << CLR_LINE;
            cout << "\nERROR: Cannot open OPIATE file:\n" << filename << "   \n" ;
            cout << "         Check path and file format!                   \n" ;
            return false;
        }

        PHDU& image = pInfile->pHDU();
        valarray<double>  contents;

        image.readAllKeys();
        image.read(contents);

        long max_row=image.axis(1);
        long max_col=image.axis(0);

        if(database_counter==0)
        {
            max_species=max_col-1;
            max_ids=max_row;

            list_freq=new double[max_species];
            list_weight=new double[max_species];
            list_IDs=new double[max_row];
        }

        if(database_counter==1)
        {
            if(max_species!=max_col-1)
            {
                cout << CLR_LINE;
                cout << "\nERROR: Ammount of species in differen OPIATE files do not match!            \n" ;
                cout << "         Identical order and values are required in all files!           \n" ;
                return false;
            }

            if(max_ids!=max_row)
            {
                cout << CLR_LINE;
                cout << "\nERROR: Ammount of unique OPIATE files do not match!            \n" ;
                cout << "         Identical order and values are required in all files!           \n" ;
                return false;
            }
        }

        if(database_counter>1)
        {
            cout << CLR_LINE;
            cout << "\nERROR: Only two databases are currently supported!              \n" ;
            return false;
        }

        cout << CLR_LINE;
        cout << "-> Scanning for keys ...              \r" << flush;

        for(uint i=0; i<max_species;i++)
        {
            char str_tmp[1024];
            char str_end[1024];

            //copy for WINDOWS has to be adjusted here

            strcpy(str_tmp, "%03d");
            sprintf(str_end, str_tmp, i+1);

            string key_name="SNA_";
            string key_weight="SWE_";
            string key_freq="SFR_";

            key_name+=str_end;
            key_weight+=str_end;
            key_freq+=str_end;

            string s_name;
            double s_weight;
            double s_frequ;

            try
            {
                image.readKey(key_name, s_name);
            }
            catch(CCfits::HDU::NoSuchKeyword)
            {
                cout << CLR_LINE;
                cout << "\nERROR: Keyword \""<< key_name << "\" is required in file:\n       " << filename <<  "               \n" ;
                return false;
            }

            try
            {
                image.readKey(key_weight, s_weight);
            }
            catch(CCfits::HDU::NoSuchKeyword)
            {
                cout << CLR_LINE;
                cout << "\nERROR: Keyword \""<< key_weight << "\" is required in file:\n       " << filename <<  "               \n" ;
                return false;
            }

            try
            {
                image.readKey(key_freq, s_frequ);
            }
            catch(CCfits::HDU::NoSuchKeyword)
            {
                cout << CLR_LINE;
                cout << "\nERROR: Keyword \""<< key_freq << "\" is required in file:\n       " << filename <<  "               \n" ;
                return false;
            }

            //cout << "Data:" << key_name << ":\t\"" << s_name << "\"\n" << flush;
            //cout << "Data:" << key_weight << ":\t\"" << s_weight << "\"\n" << flush;
            //cout << "Data:" << key_freq << ":\t\"" << s_frequ << "\"\n" << flush;

            if(database_counter==0)
            {
                list_names.push_back(s_name);
                list_freq[i]=s_frequ;
                list_weight[i]=s_weight;
            }
            else
            {
                if(list_freq[i]!=s_frequ)
                {
                    cout << CLR_LINE;
                    cout << "\nERROR: Frequencies of different OPIATE files do not match!            \n" ;
                    cout << "         Identical order and values are required in all files!           \n" ;
                    return false;
                }

                if(list_names[i].compare(s_name)!=0)
                {
                    cout << CLR_LINE;
                    cout << "\nERROR: Names of different OPIATE files do not match!            \n" ;
                    cout << "         Identical order and values are required in all files!           \n" ;
                    return false;
                }

                if(list_weight[i]!=s_weight)
                {
                    cout << CLR_LINE;
                    cout << "\nERROR: Weigths of different OPIATE files do not match!            \n" ;
                    cout << "         Identical order and values are required in all files!           \n" ;
                    return false;
                }


            }
        }

        Matrix2D tmp_mat;

        //cout << "row: " << max_row << "\n";
        //cout << "col: " << max_col << "\n";

        mat.resize(max_row,max_col-1);
        tmp_mat.resize(max_row,max_col);
        //cout << contents.size() << "\n";

        cout << CLR_LINE;
        cout << "Creating matrix ...                \r" << flush;
        for (long i = 0; i < max_col*max_row; i++)
        {
            tmp_mat.set(i,contents[i]);
        }

        //cout << "\n";

        //tmp_mat.printMatrix();

        cout << CLR_LINE;
        cout << "-> Creating ID table ...              \r" << flush;
        for (long i = 0; i < max_row; i++)
        {
            if(database_counter==0)
            {
                list_IDs[i]=tmp_mat(i,0);
            }
            else
            {
                if(list_IDs[i]!=tmp_mat(i,0))
                {
                    cout << "\nERROR: Unique OPIATE IDs of different OPIATE files do not match!            \n" ;
                    return false;
                }
            }

            cout << i << "\t" <<list_IDs[i] << "\n";

            if(i>0)
            {
                if(list_IDs[i-1]>list_IDs[i])
                {
                    cout << CLR_LINE;
                    cout << "\nERROR: OPIATE IDs are not in ascending order!               \n" ;
                    cout << "         Check first column of your input fits file!                   \n" ;
                    return false;
                }

                if(list_IDs[i-1]==list_IDs[i])
                {
                    cout << CLR_LINE;
                    cout << "\nERROR: Identical OPIATE IDs detected!               \n" ;
                    cout << "         Check first column of your input fits file!                   \n" ;
                    return false;
                }
            }
        }

        cout << "\nMatrix:\n";

        for (long i = 0; i < max_row; i++)
        {
            for (long j = 1; j < max_col; j++)
            {
                mat(i,j-1) = tmp_mat(i,j);
            }
        }

        //cout << "\n";

        mat.printMatrix();


    database_counter++;
    return true;
}
