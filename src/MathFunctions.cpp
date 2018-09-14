#include "MathFunctions.h"

bool CMathFunctions::writeSEDStatistics(string path, bool fin, Vector3D axis1,
        Vector3D axis2)
{
    string data_path = path + "statistiscs.dat";

    ofstream data_writer(data_path.c_str());

    cout << CLR_LINE;
    cout << "-> Writing statistics                               \r";

    if(data_writer.fail())
    {
        cout << "ERROR: Can't write to " << data_path << endl;
        return false;
    }

    for(uint i = 0; i < nr_ofSeq; i++)
    {
        stat_data(i, 6) = stat_data(i, 7);

        if(stat_data(i, 5) != 0)
        {
            stat_data(i, 7) /= stat_data(i, 5);
            stat_data(i, 10) /= stat_data(i, 5);
            stat_data(i, 13) /= stat_data(i, 5);
            stat_data(i, 16) /= stat_data(i, 5);
            stat_data(i, 19) /= stat_data(i, 5);
            stat_data(i, 22) /= stat_data(i, 5);
            stat_data(i, 25) /= stat_data(i, 5);
        }

        if(stat_data(i, 8) == 1e200)
            stat_data(i, 8) = 0;
        if(stat_data(i, 9) == -1e200)
            stat_data(i, 9) = 0;

        if(stat_data(i, 11) == 1e200)
            stat_data(i, 11) = 0;
        if(stat_data(i, 12) == -1e200)
            stat_data(i, 12) = 0;

        if(stat_data(i, 14) == 1e200)
            stat_data(i, 14) = 0;
        if(stat_data(i, 15) == -1e200)
            stat_data(i, 15) = 0;

        if(stat_data(i, 17) == 1e200)
            stat_data(i, 17) = 0;
        if(stat_data(i, 18) == -1e200)
            stat_data(i, 18) = 0;

        if(stat_data(i, 20) == 1e200)
            stat_data(i, 20) = 0;
        if(stat_data(i, 21) == -1e200)
            stat_data(i, 21) = 0;

        if(stat_data(i, 23) == 1e200)
            stat_data(i, 23) = 0;
        if(stat_data(i, 24) == -1e200)
            stat_data(i, 24) = 0;

        if(stat_data(i, 26) == 1e200)
            stat_data(i, 26) = 0;
        if(stat_data(i, 27) == -1e200)
            stat_data(i, 27) = 0;

        //stat_data(i, 3) *= 180.0 / PI;
        //stat_data(i, 4) *= 180.0 / PI;
    }

    data_writer << "#rot. axis: n1 = (" << axis1.X() << "," << axis1.Y() << ","
            << axis1.Z() << "); ";
    data_writer << "n2 = (" << axis2.X() << ", " << axis2.Y() << ", "
            << axis2.Z() << ")" << endl;

    if(fin)
        data_writer
            << "#dust ray det. ID\tsource ID\twave [m]\tangle1 [deg]\tangle2 [deg]\tvalid pixel\tI_SED [Jy]\tI_mean [Jy]\tI_min  [Jy]\tI_max  [Jy]\tQ_mean [Jy]\tQ_min [Jy]\tQ_max [Jy]\tU_mean [Jy]\tU_min [Jy]\tU_max [Jy]\tV_mean [Jy]\tV_min [Jy]\tV_max [Jy]\ttau_mean\ttau_min\ttau_max\tPl_mean\tPl_min\tPl_max\tPc_mean\tPc_min\tPc_max" << endl;
    else
        data_writer
            << "#wave. ID\tdet. ID\twave [m]\tangle1 [deg]\tangle2 [deg]\tvalid pixel\tI_SED [Jy]\tI_mean [Jy]\tI_min  [Jy]\tI_max  [Jy]\tQ_mean [Jy]\tQ_min [Jy]\tQ_max [Jy]\tU_mean [Jy]\tU_min [Jy]\tU_max [Jy]\tV_mean [Jy]\tV_min [Jy]\tV_max [Jy]\ttau_mean\ttau_min\ttau_max\tPl_mean\tPl_min\tPl_max\tPc_mean\tPc_min\tPc_max" << endl;

    for(uint i = 0; i < nr_ofSeq; i++)
    {
        cout << "-> Writing statistics file pos. " << i + 1 << " of "
                << nr_ofSeq << "     \r";

        for(uint j = 0; j < 28; j++)
        {
            if(j < 2 || (j > 2 && j < 6))
                data_writer << uint(stat_data(i, j) + 0.5) << "\t";
            else
                data_writer << stat_data(i, j) << "\t";
        }

        data_writer << endl;
    }

    data_writer.close();
    cout << CLR_LINE;

    if(fin)
        cout << "- Writing statistics file  : done" << endl;
    return true;
}

void CMathFunctions::initSEDStatistics(uint _nr_ofSeq)
{
    nr_ofSeq = _nr_ofSeq;
    stat_data.resize(_nr_ofSeq, 28);

    for(uint i = 0; i < nr_ofSeq; i++)
    {
        stat_data(i, 8) = 1e200;
        stat_data(i, 11) = 1e200;
        stat_data(i, 14) = 1e200;
        stat_data(i, 17) = 1e200;
        stat_data(i, 20) = 1e200;
        stat_data(i, 23) = 1e200;
        stat_data(i, 26) = 1e200;

        stat_data(i, 9) = -1e200;
        stat_data(i, 12) = -1e200;
        stat_data(i, 15) = -1e200;
        stat_data(i, 18) = -1e200;
        stat_data(i, 21) = -1e200;
        stat_data(i, 24) = -1e200;
        stat_data(i, 27) = -1e200;
    }
}
