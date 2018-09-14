#pragma once
#include "typedefs.h"
#include "chelper.h"

#ifndef COMMANDPARSER
#define COMMANDPARSER

#define id_tsk 1
#define id_cmn 2

typedef vector<parameter> parameter_list;
//int counter=0;

class CCommandParser
{
public:

    CCommandParser(void)
    {
        cmd_filename = "";
        tag = 0;
        line_counter = 0;
    }

    CCommandParser(string filename)
    {
        CCommandParser();
        cmd_filename = filename;
        tag = 0;
        line_counter = 0;
    }

    CCommandParser(const char * filename)
    {
        CCommandParser();
        cmd_filename = filename;
        tag = 0;
        line_counter = 0;
    }

    ~CCommandParser(void)
    {
    }

    parameter_list getParameterList()
    {
        return param_list;
    }


    bool checkPixel(dlist & values, dlist nr_of_pixel, bool nsides_as_pixel=false);
    bool checkVelChannels(dlist & values, dlist nr_of_channels);

    dlist parseValues(string & str);
    void formatLine(string &line);
    bool parse();
    bool parseLine(parameter * param, string cmd, string data, uint id);
    dlist parseDataString(string data);
    string seperateString(string & str);

private:
    int tag;
    string cmd_filename;
    //parameter param;
    uint line_counter;
    parameter_list param_list;

    bool atob(int b)
    {
        if(b == 1)
            return true;
        return false;
    }
};

#endif
