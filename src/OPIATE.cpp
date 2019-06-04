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
}*/

bool OPIATE::readFitsData(string filename)
{
        cout << CLR_LINE; 
        cout << "Reading OPIATE fits data from:\n       " << filename <<  "               \n" ;
        auto_ptr<FITS> pInfile(0);
      
        try
        {
            pInfile.reset(new FITS(filename.c_str(),Read,true));
        }
        catch(CCfits::FITS::CantOpen)
        {
            cout << "ERROR: Cannot open OPIATE file!               \n" ;
            cout << "       Check path and file name!                   \n" ;
            return false;
        }
        
        PHDU& image = pInfile->pHDU();
        valarray<double>  contents;
        
        image.readAllKeys();
        image.read(contents);
                
        long max_row=image.axis(1);
        long max_col=image.axis(0);
        max_species=max_col-1;
        max_ids=max_row;
        
        list_freq=new double[max_species];
        list_weight=new double[max_species];
        list_IDs=new double[max_row];
        
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
                cout << "ERROR: Keyword \""<< key_name << "\" is required in file:\n       " << filename <<  "               \n" ;
                return false;
            }
            
            try
            {
                image.readKey(key_weight, s_weight);
            }
            catch(CCfits::HDU::NoSuchKeyword)
            {
                cout << "ERROR: Keyword \""<< key_weight << "\" is required in file:\n       " << filename <<  "               \n" ;
                return false;
            }
            
            try
            {
                image.readKey(key_freq, s_frequ);
            }
            catch(CCfits::HDU::NoSuchKeyword)
            {
                cout << "ERROR: Keyword \""<< key_freq << "\" is required in file:\n       " << filename <<  "               \n" ;
                return false;
            }
            
            cout << "Data:" << key_name << ":\t\"" << s_name << "\"\n" << flush; 
            cout << "Data:" << key_weight << ":\t\"" << s_weight << "\"\n" << flush; 
            cout << "Data:" << key_freq << ":\t\"" << s_frequ << "\"\n" << flush; 
            
            list_names.push_back(s_name);
            list_freq[i]=s_frequ;
            list_weight[i]=s_weight;
        }
        
        Matrix2D tmp_mat;
        Matrix2D mat;
        
        cout << "row: " << max_row << "\n";
        cout << "col: " << max_col << "\n";
        
        mat.resize(max_row,max_col-1);
        tmp_mat.resize(max_row,max_col);
        cout << contents.size() << "\n";
        
        for (long i = 0; i < max_col*max_row; i++)
        {
            tmp_mat.set(i,contents[i]);
        }
        
        cout << "\n";
        
        tmp_mat.printMatrix();
        
        cout << "\nIDs:\n";
        
        
        for (long i = 0; i < max_row; i++)
        {
            list_IDs[i]=tmp_mat(i,0);
            cout << i << "\t" <<list_IDs[i] << "\n";
            
            if(i>0)
            {
                if(list_IDs[i-1]>list_IDs[i])
                {
                    cout << "ERROR: OPIATE IDs are not in ascending order!               \n" ;
                    cout << "       Check first column of your input fits file!                   \n" ;
                    return false;
                }
                
                if(list_IDs[i-1]==list_IDs[i])
                {
                    cout << "ERROR: Identical OPIATE IDs detected!               \n" ;
                    cout << "       Check first column of your input fits file!                   \n" ;
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
        
        cout << "\n";
        
        mat.printMatrix();
        
        

    return true;
}


