#pragma once
#include "typedefs.h"
#include "Vector.h"
#include "Stokes.h"
#include "chelper.h"
#include "GasSpecies.h"
#include <CCfits/CCfits>
#include <cmath>

#ifndef DETECTOR
#define DETECTOR

class CDetector
{
public:

    CDetector()
    {
        detector_id = -1;

        rot_angle1 = 0;
        rot_angle2 = 0;
        ex.set(1, 0, 0);
        ey.set(0, 1, 0);
        ez.set(0, 0, 1);

        obs_pos.set(0, 0, 0);
        obs_vel.set(0, 0, 0);
        l_min = 0;
        l_max = 0;
        b_min = 0;
        b_max = 0;

        axis1.set(1, 0, 0);
        axis2.set(0, 1, 0);
        path = "";
        ID = 0;
        bins_x = 1;
        bins_y = 1;
        max_cells = 1;
        distance = 1;
        i_trans = 0;

        sidelength_x = 1;
        sidelength_y = 1;

        map_shift_x = 0;
        map_shift_y = 0;

        matrixI = 0;
        matrixQ = 0;
        matrixU = 0;
        matrixV = 0;
        matrixT = 0;
        matrixS = 0;

        sedI = 0;
        sedQ = 0;
        sedU = 0;
        sedV = 0;
        sedT = 0;

        lam_min = 0;
        lam_max = 0;

        nr_of_spectral_bins = 1;
        channel_width = 0;

        cos_acceptance_angle = 0;
    }

    // Detector for dust and synchrotron
    // Plane detector
    CDetector(uint _detector_id, string _path, uint _bins_x, uint _bins_y, uint _id, double _sidelength_x,
        double _sidelength_y, double _map_shift_x, double _map_shift_y, double _distance,
        double _l_min, double _l_max, uint _nr_of_spectral_bins, uint nr_extra = 1)
    {
        detector_id = _detector_id;

        rot_angle1 = 0;
        rot_angle2 = 0;
        ex.set(1, 0, 0);
        ey.set(0, 1, 0);
        ez.set(0, 0, 1);
        path = _path;
        path += "data";
        path += SEP;

        ID = _id;
        bins_x = _bins_x;
        bins_y = _bins_y;
        max_cells = bins_x * bins_y;

        sidelength_x = _sidelength_x;
        sidelength_y = _sidelength_y;

        map_shift_x = _map_shift_x;
        map_shift_y = _map_shift_y;

        channel_width = 0;
        i_trans = 0;
        cos_acceptance_angle = 0;

        distance = _distance;

        lam_min = _l_min;
        lam_max = _l_max;
        nr_of_spectral_bins = _nr_of_spectral_bins;

        wavelength_list_det.resize(nr_of_spectral_bins);
        CMathFunctions::LogList(lam_min, lam_max, wavelength_list_det, 10);

        sedI = new double[nr_extra * nr_of_spectral_bins];
        sedQ = new double[nr_extra * nr_of_spectral_bins];
        sedU = new double[nr_extra * nr_of_spectral_bins];
        sedV = new double[nr_extra * nr_of_spectral_bins];
        sedT = new double[nr_extra * nr_of_spectral_bins];

        matrixI = new Matrix2D[nr_extra * nr_of_spectral_bins];
        matrixQ = new Matrix2D[nr_extra * nr_of_spectral_bins];
        matrixU = new Matrix2D[nr_extra * nr_of_spectral_bins];
        matrixV = new Matrix2D[nr_extra * nr_of_spectral_bins];
        matrixT = new Matrix2D[nr_extra * nr_of_spectral_bins];
        matrixS = new Matrix2D[nr_extra * nr_of_spectral_bins];

        for(uint w = 0; w < nr_of_spectral_bins; w++)
        {
            for(uint i_extra = 0; i_extra < nr_extra; i_extra++)
            {
                matrixI[w + i_extra * nr_of_spectral_bins].resize(bins_x, bins_y);
                matrixQ[w + i_extra * nr_of_spectral_bins].resize(bins_x, bins_y);
                matrixU[w + i_extra * nr_of_spectral_bins].resize(bins_x, bins_y);
                matrixV[w + i_extra * nr_of_spectral_bins].resize(bins_x, bins_y);
                matrixT[w + i_extra * nr_of_spectral_bins].resize(bins_x, bins_y);
                matrixS[w + i_extra * nr_of_spectral_bins].resize(bins_x, bins_y);

                sedI[w + i_extra * nr_of_spectral_bins] = 0;
                sedQ[w + i_extra * nr_of_spectral_bins] = 0;
                sedU[w + i_extra * nr_of_spectral_bins] = 0;
                sedV[w + i_extra * nr_of_spectral_bins] = 0;
                sedT[w + i_extra * nr_of_spectral_bins] = 0;
            }
        }
    }

    // Spherical detector
    CDetector(string _path, uint _bins, uint _id, Vector3D obs_pos, double _sidelength,
        double _l_min, double _l_max, uint _nr_of_spectral_bins, uint nr_extra = 1)
    {
        detector_id = DET_SPHER;

        rot_angle1 = 0;
        rot_angle2 = 0;
        ex.set(1, 0, 0);
        ey.set(0, 1, 0);
        ez.set(0, 0, 1);
        path = _path;
        path += "data";
        path += SEP;

        ID = _id;
        bins_x = _bins;
        bins_y = 1;
        max_cells = bins_x * bins_y;

        sidelength_x = _sidelength;
        sidelength_y = _sidelength;

        map_shift_x = 0;
        map_shift_y = 0;

        pos = obs_pos;

        channel_width = 0;
        i_trans = 0;
        cos_acceptance_angle = 0;

        lam_min = _l_min;
        lam_max = _l_max;
        nr_of_spectral_bins = _nr_of_spectral_bins;

        wavelength_list_det.resize(nr_of_spectral_bins);
        CMathFunctions::LogList(lam_min, lam_max, wavelength_list_det, 10);

        sedI = new double[nr_extra * nr_of_spectral_bins];
        sedQ = new double[nr_extra * nr_of_spectral_bins];
        sedU = new double[nr_extra * nr_of_spectral_bins];
        sedV = new double[nr_extra * nr_of_spectral_bins];
        sedT = new double[nr_extra * nr_of_spectral_bins];

        matrixI = new Matrix2D[nr_extra * nr_of_spectral_bins];
        matrixQ = new Matrix2D[nr_extra * nr_of_spectral_bins];
        matrixU = new Matrix2D[nr_extra * nr_of_spectral_bins];
        matrixV = new Matrix2D[nr_extra * nr_of_spectral_bins];
        matrixT = new Matrix2D[nr_extra * nr_of_spectral_bins];
        matrixS = new Matrix2D[nr_extra * nr_of_spectral_bins];

        for(uint w = 0; w < nr_of_spectral_bins; w++)
        {
            for(uint i_extra = 0; i_extra < nr_extra; i_extra++)
            {
                matrixI[w + i_extra * nr_of_spectral_bins].resize(bins_x, bins_y);
                matrixQ[w + i_extra * nr_of_spectral_bins].resize(bins_x, bins_y);
                matrixU[w + i_extra * nr_of_spectral_bins].resize(bins_x, bins_y);
                matrixV[w + i_extra * nr_of_spectral_bins].resize(bins_x, bins_y);
                matrixT[w + i_extra * nr_of_spectral_bins].resize(bins_x, bins_y);
                matrixS[w + i_extra * nr_of_spectral_bins].resize(bins_x, bins_y);

                sedI[w + i_extra * nr_of_spectral_bins] = 0;
                sedQ[w + i_extra * nr_of_spectral_bins] = 0;
                sedU[w + i_extra * nr_of_spectral_bins] = 0;
                sedV[w + i_extra * nr_of_spectral_bins] = 0;
                sedT[w + i_extra * nr_of_spectral_bins] = 0;
            }
        }
        
    }
    // End detector dust and synchrotron

    //Detector for line
    // Plane detector
    CDetector(uint _detector_id, string _path, uint _bins_x, uint _bins_y, uint _id, double _sidelength_x,
        double _sidelength_y, double _map_shift_x, double _map_shift_y, double _distance, uint _i_trans,
        uint _nr_of_spectral_bins, double _max_velocity)
    {
        detector_id = _detector_id;

        rot_angle1 = 0;
        rot_angle2 = 0;
        ex.set(1, 0, 0);
        ey.set(0, 1, 0);
        ez.set(0, 0, 1);

        obs_pos.set(0, 0, 0);
        obs_vel.set(0, 0, 0);
        l_min = 0;
        l_max = 0;
        b_min = 0;
        b_max = 0;

        path = _path;
        path += "data";
        path += SEP;

        ID = _id;
        bins_x = _bins_x;
        bins_y = _bins_y;
        max_cells = bins_x * bins_y;
        distance = _distance;

        calcVelocityChannels(_nr_of_spectral_bins, _max_velocity);

        sidelength_x = _sidelength_x;
        sidelength_y = _sidelength_y;

        map_shift_x = _map_shift_x;
        map_shift_y = _map_shift_y;

        wavelength_list_det.resize(nr_of_spectral_bins);

        i_trans = _i_trans;

        cos_acceptance_angle = 0;

        sedI = new double[nr_of_spectral_bins];
        sedQ = new double[nr_of_spectral_bins];
        sedU = new double[nr_of_spectral_bins];
        sedV = new double[nr_of_spectral_bins];
        sedT = new double[nr_of_spectral_bins];

        matrixI = new Matrix2D[nr_of_spectral_bins];
        matrixQ = new Matrix2D[nr_of_spectral_bins];
        matrixU = new Matrix2D[nr_of_spectral_bins];
        matrixV = new Matrix2D[nr_of_spectral_bins];
        matrixT = new Matrix2D[nr_of_spectral_bins];
        matrixS = new Matrix2D[nr_of_spectral_bins];

        for(uint vch = 0; vch < nr_of_spectral_bins; vch++)
        {
            matrixI[vch].resize(bins_x, bins_y);
            matrixQ[vch].resize(bins_x, bins_y);
            matrixU[vch].resize(bins_x, bins_y);
            matrixV[vch].resize(bins_x, bins_y);
            matrixT[vch].resize(bins_x, bins_y);
            matrixS[vch].resize(bins_x, bins_y);

            sedI[vch] = 0;
            sedQ[vch] = 0;
            sedU[vch] = 0;
            sedV[vch] = 0;
            sedT[vch] = 0;
        }

        lam_min = 0;
        lam_max = 0;
    }

    // Spherical detector
    CDetector(string _path, uint _bins, uint _id, Vector3D obs_pos, double _sidelength, uint _i_trans,
        uint _nr_of_spectral_bins, double _max_velocity)
    {
        detector_id = DET_SPHER;

        rot_angle1 = 0;
        rot_angle2 = 0;
        ex.set(1, 0, 0);
        ey.set(0, 1, 0);
        ez.set(0, 0, 1);

        obs_pos.set(0, 0, 0);
        obs_vel.set(0, 0, 0);
        l_min = 0;
        l_max = 0;
        b_min = 0;
        b_max = 0;

        path = _path;
        path += "data";
        path += SEP;

        ID = _id;
        bins_x = _bins;
        bins_y = 1;
        max_cells = bins_x * bins_y;
        distance = 1;

        calcVelocityChannels(_nr_of_spectral_bins, _max_velocity);

        sidelength_x = _sidelength;
        sidelength_y = _sidelength;

        map_shift_x = 0;
        map_shift_y = 0;

        pos = obs_pos;

        wavelength_list_det.resize(nr_of_spectral_bins);

        i_trans = _i_trans;

        cos_acceptance_angle = 0;

        sedI = new double[nr_of_spectral_bins];
        sedQ = new double[nr_of_spectral_bins];
        sedU = new double[nr_of_spectral_bins];
        sedV = new double[nr_of_spectral_bins];
        sedT = new double[nr_of_spectral_bins];

        matrixI = new Matrix2D[nr_of_spectral_bins];
        matrixQ = new Matrix2D[nr_of_spectral_bins];
        matrixU = new Matrix2D[nr_of_spectral_bins];
        matrixV = new Matrix2D[nr_of_spectral_bins];
        matrixT = new Matrix2D[nr_of_spectral_bins];
        matrixS = new Matrix2D[nr_of_spectral_bins];

        for(uint vch = 0; vch < nr_of_spectral_bins; vch++)
        {
            matrixI[vch].resize(bins_x, bins_y);
            matrixQ[vch].resize(bins_x, bins_y);
            matrixU[vch].resize(bins_x, bins_y);
            matrixV[vch].resize(bins_x, bins_y);
            matrixT[vch].resize(bins_x, bins_y);
            matrixS[vch].resize(bins_x, bins_y);

            sedI[vch] = 0;
            sedQ[vch] = 0;
            sedU[vch] = 0;
            sedV[vch] = 0;
            sedT[vch] = 0;
        }

        lam_min = 0;
        lam_max = 0;
    }
    // End detector line

    // Dust scattering detector
    void init(string _path, uint _bins_x, uint _bins_y, double _sidelength_x, double _sidelength_y,
        double _distance, double _l_min, double _l_max, uint _nr_of_spectral_bins)
    {
        rot_angle1 = 0;
        rot_angle2 = 0;
        ex.set(1, 0, 0);
        ey.set(0, 1, 0);
        ez.set(0, 0, 1);

        obs_pos.set(0, 0, 0);
        obs_vel.set(0, 0, 0);
        l_min = 0;
        l_max = 0;
        b_min = 0;
        b_max = 0;

        path = _path;

        ID = 0;
        bins_x = _bins_x;
        bins_y = _bins_y;
        max_cells = bins_x * bins_y;

        sidelength_x = _sidelength_x;
        sidelength_y = _sidelength_y;

        map_shift_x = 0;
        map_shift_y = 0;

        channel_width = 0;
        i_trans = 0,
        cos_acceptance_angle = 0;

        distance = _distance;

        lam_min = _l_min;
        lam_max = _l_max;
        nr_of_spectral_bins = _nr_of_spectral_bins;

        wavelength_list_det.resize(nr_of_spectral_bins);
        CMathFunctions::LogList(lam_min, lam_max, wavelength_list_det, 10);

        sedI = new double[nr_of_spectral_bins];
        sedQ = new double[nr_of_spectral_bins];
        sedU = new double[nr_of_spectral_bins];
        sedV = new double[nr_of_spectral_bins];
        sedT = new double[nr_of_spectral_bins];

        matrixI = new Matrix2D[nr_of_spectral_bins];
        matrixQ = new Matrix2D[nr_of_spectral_bins];
        matrixU = new Matrix2D[nr_of_spectral_bins];
        matrixV = new Matrix2D[nr_of_spectral_bins];
        matrixT = new Matrix2D[nr_of_spectral_bins];
        matrixS = new Matrix2D[nr_of_spectral_bins];

        for(uint w = 0; w < nr_of_spectral_bins; w++)
        {
            matrixI[w].resize(bins_x, bins_y);
            matrixQ[w].resize(bins_x, bins_y);
            matrixU[w].resize(bins_x, bins_y);
            matrixV[w].resize(bins_x, bins_y);
            matrixT[w].resize(bins_x, bins_y);
            matrixS[w].resize(bins_x, bins_y);

            sedI[w] = 0;
            sedQ[w] = 0;
            sedU[w] = 0;
            sedV[w] = 0;
            sedT[w] = 0;
        }
    }

    ~CDetector()
    {
        if(matrixI != 0)
            delete[] matrixI;
        if(matrixQ != 0)
            delete[] matrixQ;
        if(matrixU != 0)
            delete[] matrixU;
        if(matrixV != 0)
            delete[] matrixV;
        if(matrixT != 0)
            delete[] matrixT;
        if(matrixS != 0)
            delete[] matrixS;

        if(sedI != 0)
            delete[] sedI;
        if(sedQ != 0)
            delete[] sedQ;
        if(sedU != 0)
            delete[] sedU;
        if(sedV != 0)
            delete[] sedV;
        if(sedT != 0)
            delete[] sedT;
    }

    Vector3D getDirection()
    {
        return ez;
    }

    double getDistance()
    {
        return distance;
    }

    double getAcceptanceAngle()
    {
        return cos_acceptance_angle;
    }

    Vector3D getEX()
    {
        return ex;
    }

    Vector3D getEY()
    {
        return ey;
    }

    Vector3D getEZ()
    {
        return ez;
    }

    void setOrientation(Vector3D n1, Vector3D n2, double _rot_angle1, double _rot_angle2)
    {
        rot_angle1 = _rot_angle1;
        rot_angle2 = _rot_angle2;

        axis1 = n1;
        axis2 = n2;

        ex.set(1, 0, 0);
        ey.set(0, 1, 0);
        ez.set(0, 0, 1);

        double cos_a = cos(rot_angle1);
        double sin_a = sin(rot_angle1);

        ex.rot(n1, cos_a, sin_a);
        ey.rot(n1, cos_a, sin_a);
        ez.rot(n1, cos_a, sin_a);

        cos_a = cos(rot_angle2);
        sin_a = sin(rot_angle2);

        ex.rot(n2, cos_a, sin_a);
        ey.rot(n2, cos_a, sin_a);
        ez.rot(n2, cos_a, sin_a);

        ex.normalize();
        ey.normalize();
        ez.normalize();
    }

    void setObsPosition(Vector3D _obs_pos, Vector3D _obs_vel,
            double _l_min, double _l_max, double _b_min, double _b_max)
    {
        obs_pos = _obs_pos;
        obs_vel = _obs_vel;

        l_min = _l_min;
        l_max = _l_max;

        b_min = _b_min;
        b_max = _b_max;
    }

    void setAcceptanceAngle(double acceptance_angle)
    {
        cos_acceptance_angle = cos(acceptance_angle);
    }

    void addToSedDetector(StokesVector st, uint i_spectral)
    {
        if(sedI != 0)
        {
#pragma omp atomic update
            sedI[i_spectral] += st.I();
#pragma omp atomic update
            sedQ[i_spectral] += st.Q();
#pragma omp atomic update
            sedU[i_spectral] += st.U();
#pragma omp atomic update
            sedV[i_spectral] += st.V();
#pragma omp atomic update
            sedT[i_spectral] += st.T() / max_cells;
        }
    }

    void addToRaytracingSedDetector(photon_package * pp, uint spectral_offset=0)
    {
        uint i_spectral = pp->getWavelengthID();
        StokesVector st = pp->getMultiStokesVector(i_spectral);  
        addToSedDetector(st, i_spectral + spectral_offset);
    }

    void addToRaytracingDetector(photon_package * pp, uint pos_id, uint spectral_offset)
    {
        uint i_spectral = pp->getWavelengthID();
        StokesVector st = pp->getMultiStokesVector(i_spectral);        

        matrixI[i_spectral + spectral_offset].addValue(pos_id, st.I());
        matrixQ[i_spectral + spectral_offset].addValue(pos_id, st.Q());
        matrixU[i_spectral + spectral_offset].addValue(pos_id, st.U());
        matrixV[i_spectral + spectral_offset].addValue(pos_id, st.V());
        matrixT[i_spectral + spectral_offset].addValue(pos_id, st.T());
        matrixS[i_spectral + spectral_offset].addValue(pos_id, st.Sp());
    }

    void addToRaytracingDetector(photon_package * pp, uint spectral_offset)
    {
        uint i_spectral = pp->getWavelengthID();
        StokesVector st = pp->getMultiStokesVector(i_spectral);
        Vector3D pos = pp->getPosition();

        uint x = uint((pos.X() + 0.5 * sidelength_x - map_shift_x) / sidelength_x * double(bins_x));

        if(x < 0 || x >= int(bins_x))
            return;

        uint y = uint((pos.Y() + 0.5 * sidelength_y - map_shift_y) / sidelength_y * double(bins_y));

        if(y < 0 || y >= int(bins_y))
            return;

        matrixI[i_spectral + spectral_offset].addValue(x, y, st.I());
        matrixQ[i_spectral + spectral_offset].addValue(x, y, st.Q());
        matrixU[i_spectral + spectral_offset].addValue(x, y, st.U());
        matrixV[i_spectral + spectral_offset].addValue(x, y, st.V());
        matrixT[i_spectral + spectral_offset].addValue(x, y, st.T());
        matrixS[i_spectral + spectral_offset].addValue(x, y, st.Sp());
    }

    void addToMonteCarloDetector(photon_package * pp, uint i_spectral, uint radiation_type)
    {
        Vector3D pos = pp->getPosition();
        Vector3D dir = pp->getDirection();

        if(dir.length() == 0)
            return;

        double lx = pos * ex;
        double ly = pos * ey;

        int x = int((lx + 0.5 * sidelength_x) / sidelength_x * double(bins_x));

        if(x < 0 || x >= int(bins_x))
            return;

        int y = int((ly + 0.5 * sidelength_y) / sidelength_y * double(bins_y));

        if(y < 0 || y >= int(bins_y))
            return;

        StokesVector st = pp->getStokesVector();

        matrixI[i_spectral].addValue(x, y, st.I());
        matrixQ[i_spectral].addValue(x, y, st.Q());
        matrixU[i_spectral].addValue(x, y, st.U());
        matrixV[i_spectral].addValue(x, y, st.V());

        // Add to SED
        if(sedI != 0)
        {
#pragma omp atomic update
            sedI[i_spectral] += st.I();
#pragma omp atomic update
            sedQ[i_spectral] += st.Q();
#pragma omp atomic update
            sedU[i_spectral] += st.U();
#pragma omp atomic update
            sedV[i_spectral] += st.V();
        }

        if(radiation_type == DIRECT_STAR)
        {
            matrixT[i_spectral].addValue(x, y, st.I());
        }
        else
        {
            matrixS[i_spectral].addValue(x, y, st.I());
            if(sedI != 0)
            {
#pragma omp atomic update
                sedT[i_spectral] += st.I();
            }
        }
    }

    bool writeMap(uint nr, uint results_type)
    {
        cout << CLR_LINE << flush;
        cout << " -> Writing map ...  \r" << flush;
        auto_ptr<CCfits::FITS> pFits(0);
        //unique_ptr<CCfits::FITS> pFits;

        try
        {
            char str_tmp[1024];
            char str_end[1024];
            nr++;

#ifdef WINDOWS
            strcpy_s(str_tmp, "polaris_detector_nr%04d");
            sprintf_s(str_end, str_tmp, nr);
#else
            strcpy(str_tmp, "polaris_detector_nr%04d");
            sprintf(str_end, str_tmp, nr);
#endif

            string path_out = path + str_end + ".fits";
            remove(path_out.c_str());

            long naxis = 4;
            long naxes[4] = {uint(bins_x), uint(bins_y), nr_of_spectral_bins, 6};

            pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
        }
        catch(CCfits::FITS::CantCreate){ return false; }

        long nelements = uint(bins_x) * uint(bins_y);
        if(max_cells != nelements)
        {
            cout << "\nWARNING: Max cells are not equal to bins x bins!" << endl;
            return false;
        }

        vector<long> fpixel(4);
        fpixel[0] = 1;
        fpixel[1] = 1;

        for(uint i_spectral = 0; i_spectral < nr_of_spectral_bins; i_spectral++)
        {
            fpixel[2] = i_spectral + 1;

            valarray<double> array_I(nelements);
            valarray<double> array_Q(nelements);
            valarray<double> array_U(nelements);
            valarray<double> array_V(nelements);
            valarray<double> array_T(nelements);
            valarray<double> array_S(nelements);

            uint i = 0;
            for(uint i_y = 0; i_y < bins_y; i_y++)
                for(uint i_x = 0; i_x < bins_x; i_x++)
                {
                    array_I[i] = matrixI[i_spectral](i_x, i_y);
                    array_Q[i] = matrixQ[i_spectral](i_x, i_y);
                    array_U[i] = matrixU[i_spectral](i_x, i_y);
                    array_V[i] = matrixV[i_spectral](i_x, i_y);
                    array_T[i] = matrixT[i_spectral](i_x, i_y);
                    array_S[i] = matrixS[i_spectral](i_x, i_y);
                    i++;
                }

            fpixel[3] = 1;
            pFits->pHDU().write(fpixel, nelements, array_I);
            fpixel[3] = 2;
            pFits->pHDU().write(fpixel, nelements, array_Q);
            fpixel[3] = 3;
            pFits->pHDU().write(fpixel, nelements, array_U);
            fpixel[3] = 4;
            pFits->pHDU().write(fpixel, nelements, array_V);
            fpixel[3] = 5;
            pFits->pHDU().write(fpixel, nelements, array_T);
            fpixel[3] = 6;
            pFits->pHDU().write(fpixel, nelements, array_S);
        }

        double bin_width_x, first_pix_val_x, deg_per_pix_x, first_pix_val_deg_x;
        calcCoordinateParameters(sidelength_x, bins_x, map_shift_x, distance, bin_width_x, first_pix_val_x, deg_per_pix_x, first_pix_val_deg_x);
        double bin_width_y, first_pix_val_y, deg_per_pix_y, first_pix_val_deg_y;
        calcCoordinateParameters(sidelength_y, bins_y, map_shift_y, distance, bin_width_y, first_pix_val_y, deg_per_pix_y, first_pix_val_deg_y);

        // Grid
        pFits->pHDU().addKey("CTYPE1", "PARAM", "type of unit 1");
        pFits->pHDU().addKey("CRVAL1", first_pix_val_x, "value of axis 1");
        pFits->pHDU().addKey("CRPIX1", 1, "pixel where CRVAL1 is defined ");
        pFits->pHDU().addKey("CDELT1", bin_width_x, "delta of axis 1");
        pFits->pHDU().addKey("CUNIT1", "m", "unit of axis 1");

        // Alternatively as arc seconds grid
        pFits->pHDU().addKey("CTYPE1A", "RA---TAN", "type of unit 1");
        pFits->pHDU().addKey("CRVAL1A", first_pix_val_deg_x, "value of axis 1");
        pFits->pHDU().addKey("CRPIX1A", bins_x, "pixel where CRVAL1A is defined ");
        pFits->pHDU().addKey("CDELT1A", -deg_per_pix_x, "delta of axis 1");
        pFits->pHDU().addKey("CUNIT1A", "degree", "unit of axis 1");

        // Alternatively as AU grid
        pFits->pHDU().addKey("CTYPE1B", "PARAM", "type of unit 1");
        pFits->pHDU().addKey("CRVAL1B", first_pix_val_x / con_AU, "value of axis 1");
        pFits->pHDU().addKey("CRPIX1B", 1, "pixel where CRVAL1 is defined ");
        pFits->pHDU().addKey("CDELT1B", bin_width_x / con_AU, "delta of axis 1");
        pFits->pHDU().addKey("CUNIT1B", "AU", "unit of axis 1");

        // Alternatively as pc grid
        pFits->pHDU().addKey("CTYPE1C", "PARAM", "type of unit 1");
        pFits->pHDU().addKey("CRVAL1C", first_pix_val_x / con_pc, "value of axis 1");
        pFits->pHDU().addKey("CRPIX1C", 1, "pixel where CRVAL1 is defined ");
        pFits->pHDU().addKey("CDELT1C", bin_width_x / con_pc, "delta of axis 1");
        pFits->pHDU().addKey("CUNIT1C", "pc", "unit of axis 1");

        // Grid
        pFits->pHDU().addKey("CTYPE2", "PARAM", "type of unit 2");
        pFits->pHDU().addKey("CRVAL2", first_pix_val_y, "value of axis 2");
        pFits->pHDU().addKey("CRPIX2", 1, "pixel where CRVAL2 is defined ");
        pFits->pHDU().addKey("CDELT2", bin_width_y, "delta of axis 2");
        pFits->pHDU().addKey("CUNIT2", "m", "unit of axis 2");

        // Alternatively as arc seconds grid
        pFits->pHDU().addKey("CTYPE2A", "DEC--TAN", "type of unit 2");
        pFits->pHDU().addKey("CRVAL2A", first_pix_val_deg_y, "value of axis 2");
        pFits->pHDU().addKey("CRPIX2A", 1, "pixel where CRVAL2A is defined ");
        pFits->pHDU().addKey("CDELT2A", deg_per_pix_y, "delta of axis 2");
        pFits->pHDU().addKey("CUNIT2A", "degree", "unit of axis 2");

        // Alternatively as AU grid
        pFits->pHDU().addKey("CTYPE2B", "PARAM", "type of unit 2");
        pFits->pHDU().addKey("CRVAL2B", first_pix_val_y / con_AU, "value of axis 2");
        pFits->pHDU().addKey("CRPIX2B", 1, "pixel where CRVAL2 is defined ");
        pFits->pHDU().addKey("CDELT2B", bin_width_y / con_AU, "delta of axis 2");
        pFits->pHDU().addKey("CUNIT2B", "AU", "unit of axis 2");

        // Alternatively as pc grid
        pFits->pHDU().addKey("CTYPE2C", "PARAM", "type of unit 2");
        pFits->pHDU().addKey("CRVAL2C", first_pix_val_y / con_pc, "value of axis 2");
        pFits->pHDU().addKey("CRPIX2C", 1, "pixel where CRVAL2 is defined ");
        pFits->pHDU().addKey("CDELT2C", bin_width_y / con_pc, "delta of axis 2");
        pFits->pHDU().addKey("CUNIT2C", "pc", "unit of axis 2");

        // Wavelength IDs
        pFits->pHDU().addKey("CTYPE3", "PARAM", "type of unit 3");
        pFits->pHDU().addKey("CRVAL3", 1, "value of axis 3");
        pFits->pHDU().addKey("CRPIX3", 1, "pixel where CRVAL3 is defined ");
        pFits->pHDU().addKey("CDELT3", 1, "delta of axis 3");
        pFits->pHDU().addKey("CUNIT3", "Wavelength index", "unit of axis 3");

        // Simulated quantities
        pFits->pHDU().addKey("CTYPE4", "PARAM", "type of unit 4");
        pFits->pHDU().addKey("CRVAL4", 1, "value of axis 4");
        pFits->pHDU().addKey("CRPIX4", 1, "pixel where CRVAL4 is defined ");
        pFits->pHDU().addKey("CDELT4", 1, "delta of axis 4");
        if(results_type == RESULTS_RAY)
        {
            pFits->pHDU().addKey("CUNIT4",
                "I, Q, U, V [Jy/px], optical depth, column density [m^-2]", "unit of axis 4");
            pFits->pHDU().addKey("ETYPE", "thermal emission", "type of emission");
        }
        else if(results_type == RESULTS_FULL)
        {
            pFits->pHDU().addKey("CUNIT4",
                "I, Q, U, V [Jy/px], optical depth, column density [m^-2]", "unit of axis 4");
            pFits->pHDU().addKey("ETYPE", "thermal emission (+scattering)", "type of emission");
        }
        else
        {
            pFits->pHDU().addKey("CUNIT4", "I, Q, U, V, I_direct, I_scat [Jy/px]", "unit of axis 4");
            pFits->pHDU().addKey("ETYPE", "scattered emission / direct stellar emission", "type of emission");
        }
        pFits->pHDU().addKey("ID", nr, "detector ID");

        for(uint i_spectral = 0; i_spectral < nr_of_spectral_bins; i_spectral++)
        {
            char str_1[1024];
            char str_2[1024];
#ifdef WINDOWS
            sprintf_s(str_1, "WAVELENGTH%i", i_spectral + 1);
            sprintf_s(str_2, "value of %i. wavelength", i_spectral + 1);
#else
            sprintf(str_1, "WAVELENGTH%i", i_spectral + 1);
            sprintf(str_2, "value of %i. wavelength", i_spectral + 1);
#endif
            pFits->pHDU().addKey(str_1, wavelength_list_det[i_spectral], str_2);
        }
        pFits->pHDU().addKey("DISTANCE", distance, "distance to object");
        pFits->pHDU().addKey("RAXIS1X", axis1.X(), "rotation axes 1 (x component)");
        pFits->pHDU().addKey("RAXIS1Y", axis1.Y(), "rotation axes 1 (y component)");
        pFits->pHDU().addKey("RAXIS1Z", axis1.Z(), "rotation axes 1 (z component)");
        pFits->pHDU().addKey("RANGLE1", float(180.0 * rot_angle1 / PI), "rotation angle 1 [deg]");
        pFits->pHDU().addKey("RAXIS2X", axis2.X(), "rotation axes 2 (x component)");
        pFits->pHDU().addKey("RAXIS2Y", axis2.Y(), "rotation axes 2 (y component)");
        pFits->pHDU().addKey("RAXIS2Z", axis2.Z(), "rotation axes 2 (z component)");
        pFits->pHDU().addKey("RANGLE2", float(180.0 * rot_angle2 / PI), "rotation angle 2 [deg]");
        pFits->pHDU().addKey("DETGRID", getDetectorGridDescription(), "description of the detector grid");

        cout << CLR_LINE << flush;
        return true;
    }

    bool writeSed(uint nr, uint results_type)
    {
        cout << CLR_LINE << flush;
        cout << " -> Writing line spectrum ...  \r" << flush;
        auto_ptr<CCfits::FITS> pFits(0);
        //unique_ptr<CCfits::FITS> pFits;

        try
        {
            char str_tmp[1024];
            char str_end[1024];
            nr++;

#ifdef WINDOWS
            strcpy_s(str_tmp, "polaris_detector_nr%04d_sed");
            sprintf_s(str_end, str_tmp, nr);
#else
            strcpy(str_tmp, "polaris_detector_nr%04d_sed");
            sprintf(str_end, str_tmp, nr);
#endif

            string path_out = path + str_end + ".fits";
            remove(path_out.c_str());

            long naxis = 3;
            long naxes[3] = {nr_of_spectral_bins, 1, 5};
            pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
        }
        catch(CCfits::FITS::CantCreate){ return false; }

        long nelements = uint(nr_of_spectral_bins);

        vector<long> fpixel(3);
        fpixel[0] = 1;
        fpixel[1] = 1;

        std::valarray<double> array_I(nelements);
        std::valarray<double> array_Q(nelements);
        std::valarray<double> array_U(nelements);
        std::valarray<double> array_V(nelements);
        std::valarray<double> array_T(nelements);
        for(uint i_spectral = 0; i_spectral < nr_of_spectral_bins; i_spectral++)
        {
            array_I[i_spectral] = sedI[i_spectral];
            array_Q[i_spectral] = sedQ[i_spectral];
            array_U[i_spectral] = sedU[i_spectral];
            array_V[i_spectral] = sedV[i_spectral];
            array_T[i_spectral] = sedT[i_spectral];

            fpixel[2] = 1;
            pFits->pHDU().write(fpixel, nelements, array_I);
            fpixel[2] = 2;
            pFits->pHDU().write(fpixel, nelements, array_Q);
            fpixel[2] = 3;
            pFits->pHDU().write(fpixel, nelements, array_U);
            fpixel[2] = 4;
            pFits->pHDU().write(fpixel, nelements, array_V);
            fpixel[2] = 5;
            pFits->pHDU().write(fpixel, nelements, array_T);
        }

        pFits->pHDU().addKey("CTYPE1", "PARAM", "type of unit 1");
        pFits->pHDU().addKey("CRVAL1", 1, "value of axis 1");
        pFits->pHDU().addKey("CRPIX1", 1, "pixel where CRVAL1 is defined ");
        pFits->pHDU().addKey("CDELT1", 1, "delta of axis 1");
        pFits->pHDU().addKey("CUNIT1", "Wavelength index", "unit of axis 1");
        pFits->pHDU().addKey("CTYPE2", "PARAM", "type of unit 2");
        pFits->pHDU().addKey("CRVAL2", 1, "value of axis 2");
        pFits->pHDU().addKey("CRPIX2", 1, "pixel where CRVAL2 is defined ");
        pFits->pHDU().addKey("CDELT2", 1, "delta of axis 2");
        pFits->pHDU().addKey("CUNIT2", "None", "unit of axis 2");
        if(results_type == RESULTS_RAY)
        {
            pFits->pHDU().addKey("CUNIT3", "I, Q, U, V [Jy], optical depth", "unit of axis 3");
            pFits->pHDU().addKey("ETYPE", "thermal emission", "type of emission");
        }
        else if(results_type == RESULTS_FULL)
        {
            pFits->pHDU().addKey("CUNIT3", "I, Q, U, V [Jy], optical depth", "unit of axis 3");
            pFits->pHDU().addKey("ETYPE", "thermal emission (+scattering)", "type of emission");
        }
        else
        {
            pFits->pHDU().addKey("CUNIT3", "I, Q, U, V, I_scat [Jy/px]", "unit of axis 3");
            pFits->pHDU().addKey("ETYPE", "scattered emission / direct stellar emission", "type of emission");
        }
        pFits->pHDU().addKey("ID", nr, "detector ID");
        for(uint i_spectral = 0; i_spectral < nr_of_spectral_bins; i_spectral++)
        {
            char str_1[1024];
            char str_2[1024];
#ifdef WINDOWS
            sprintf_s(str_1, "WAVELENGTH%i", i_spectral + 1);
            sprintf_s(str_2, "value of %i. wavelength", i_spectral + 1);
#else
            sprintf(str_1, "WAVELENGTH%i", i_spectral + 1);
            sprintf(str_2, "value of %i. wavelength", i_spectral + 1);
#endif
            pFits->pHDU().addKey(str_1, wavelength_list_det[i_spectral], str_2);
        }
        if(detector_id == DET_SPHER)
        {
            pFits->pHDU().addKey("OBS_POSITION_X", obs_pos.X(), "x-axis position of observer");
            pFits->pHDU().addKey("OBS_POSITION_Y", obs_pos.Y(), "y-axis position of observer");
            pFits->pHDU().addKey("OBS_POSITION_Z", obs_pos.Z(), "z-axis position of observer");
            pFits->pHDU().addKey("OBS_VELOCITY_X", obs_vel.X(), "velocity of observer in x direction");
            pFits->pHDU().addKey("OBS_VELOCITY_Y", obs_vel.Y(), "velocity of observer in y direction");
            pFits->pHDU().addKey("OBS_VELOCITY_Z", obs_vel.Z(), "velocity of observer in z direction");
            pFits->pHDU().addKey("LONGITUDE_MIN", l_min, "minimum considered galactic longitude");
            pFits->pHDU().addKey("LONGITUDE_MAX", l_max, "maximum considered galactic longitude");
            pFits->pHDU().addKey("LATITUDE_MIN", b_min, "minimum considered galactic latitude");
            pFits->pHDU().addKey("LATITUDE_MAX", b_max, "maximum considered galactic latitude");
        }
        else
        {
            pFits->pHDU().addKey("DISTANCE", distance, "distance to object");
            pFits->pHDU().addKey("RAXIS1X", axis1.X(), "rotation axes 1 (x component)");
            pFits->pHDU().addKey("RAXIS1Y", axis1.Y(), "rotation axes 1 (y component)");
            pFits->pHDU().addKey("RAXIS1Z", axis1.Z(), "rotation axes 1 (z component)");
            pFits->pHDU().addKey("RANGLE1", float(180.0 * rot_angle1 / PI), "rotation angle 1 [deg]");
            pFits->pHDU().addKey("RAXIS2X", axis2.X(), "rotation axes 2 (x component)");
            pFits->pHDU().addKey("RAXIS2Y", axis2.Y(), "rotation axes 2 (y component)");
            pFits->pHDU().addKey("RAXIS2Z", axis2.Z(), "rotation axes 2 (z component)");
            pFits->pHDU().addKey("RANGLE2", float(180.0 * rot_angle2 / PI), "rotation angle 2 [deg]");
        }
        pFits->pHDU().addKey("DETGRID", getDetectorGridDescription(), "description of the detector grid");

        cout << CLR_LINE << flush;
       
        return true;
    }

    bool writeHealMaps(uint nr, uint results_type)
    {
        cout << CLR_LINE << flush;
        cout << " -> Writing healpix map ...  \r" << flush;
        auto_ptr<CCfits::FITS> pFits(0);
        //unique_ptr<CCfits::FITS> pFits;

        try
        {
            char str_tmp[1024];
            char str_end[1024];
            nr++;

#ifdef WINDOWS
            strcpy_s(str_tmp, "polaris_detector_nr%04d");
            sprintf_s(str_end, str_tmp, nr);
#else
            strcpy(str_tmp, "polaris_detector_nr%04d");
            sprintf(str_end, str_tmp, nr);
#endif

            string path_out = path + str_end + ".fits";
            remove(path_out.c_str());

            long naxis = 1;
            long naxes[1] = {0};
            pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
        }
        catch(CCfits::FITS::CantCreate){  return false; }

        long nelements = uint(bins_x) * uint(bins_y);
        if(max_cells != nelements)
        {
            cout << "\nWARNING: Max cells are not equal to bins x bins!" << endl;
            return false;
        }

        // Init columns
        string newName("HEALPIX_EXTENSION");
        uint nr_of_quantities = 5;
        vector<string> colName(nr_of_quantities * nr_of_spectral_bins + 1, "");
        vector<string> colForm(nr_of_quantities * nr_of_spectral_bins + 1, "");
        vector<string> colUnit(nr_of_quantities * nr_of_spectral_bins + 1, "");

        for(uint i_spectral = 0; i_spectral < nr_of_spectral_bins; i_spectral++)
        {
            char str_1[1024];
#ifdef WINDOWS
            sprintf_s(str_1, "I_STOKES (WAVELENGTH = %e [m])", wavelength_list_det[i_spectral]);
#else
            sprintf(str_1, "I_STOKES (WAVELENGTH = %e [m])", wavelength_list_det[i_spectral]);
#endif
            colName[nr_of_quantities * i_spectral + 0] = str_1;
#ifdef WINDOWS
            sprintf_s(str_1, "Q_STOKES (WAVELENGTH = %e [m])", wavelength_list_det[i_spectral]);
#else
            sprintf(str_1, "Q_STOKES (WAVELENGTH = %e [m])", wavelength_list_det[i_spectral]);
#endif
            colName[nr_of_quantities * i_spectral + 1] = str_1;
#ifdef WINDOWS
            sprintf_s(str_1, "U_STOKES (WAVELENGTH = %e [m])", wavelength_list_det[i_spectral]);
#else
            sprintf(str_1, "U_STOKES (WAVELENGTH = %e [m])", wavelength_list_det[i_spectral]);
#endif
            colName[nr_of_quantities * i_spectral + 2] =  str_1;
#ifdef WINDOWS
            sprintf_s(str_1, "V_STOKES (WAVELENGTH = %e [m])", wavelength_list_det[i_spectral]);
#else
            sprintf(str_1, "V_STOKES (WAVELENGTH = %e [m])", wavelength_list_det[i_spectral]);
#endif
            colName[nr_of_quantities * i_spectral + 3] = str_1;
#ifdef WINDOWS
            sprintf_s(str_1, "OPTICAL_DEPTH (WAVELENGTH = %e [m])", wavelength_list_det[i_spectral]);
#else
            sprintf(str_1, "OPTICAL_DEPTH (WAVELENGTH = %e [m])", wavelength_list_det[i_spectral]);
#endif
            colName[nr_of_quantities * i_spectral + 4] = str_1;

            colForm[nr_of_quantities * i_spectral + 0] = "D";
            colForm[nr_of_quantities * i_spectral + 1] = "D";
            colForm[nr_of_quantities * i_spectral + 2] = "D";
            colForm[nr_of_quantities * i_spectral + 3] = "D";
            colForm[nr_of_quantities * i_spectral + 4] = "D";

            colUnit[nr_of_quantities * i_spectral + 0] = "Jy/px";
            colUnit[nr_of_quantities * i_spectral + 1] = "Jy/px";
            colUnit[nr_of_quantities * i_spectral + 2] = "Jy/px";
            colUnit[nr_of_quantities * i_spectral + 3] = "Jy/px";
            colUnit[nr_of_quantities * i_spectral + 4] = "";

        }
        colName[nr_of_quantities * nr_of_spectral_bins] = "COLUMN_DENSITY";
        colForm[nr_of_quantities * nr_of_spectral_bins] = "D";
        colUnit[nr_of_quantities * nr_of_spectral_bins] = "m^-2";

        CCfits::Table* newTable = pFits->addTable(newName, nelements, colName, colForm, colUnit);

        for(uint i_spectral = 0; i_spectral < nr_of_spectral_bins; i_spectral++)
        {
            valarray<double> array_I(nelements);
            valarray<double> array_Q(nelements);
            valarray<double> array_U(nelements);
            valarray<double> array_V(nelements);
            valarray<double> array_T(nelements);

            uint i = 0;
            for(uint i_y = 0; i_y < bins_y; i_y++)
                for(uint i_x = 0; i_x < bins_x; i_x++)
                {
                    array_I[i] = matrixI[i_spectral](i_x, i_y);
                    array_Q[i] = matrixQ[i_spectral](i_x, i_y);
                    array_U[i] = matrixU[i_spectral](i_x, i_y);
                    array_V[i] = matrixV[i_spectral](i_x, i_y);
                    array_T[i] = matrixT[i_spectral](i_x, i_y);
                    i++;
                }

            newTable->column(colName[nr_of_quantities * i_spectral + 0]).write(array_I, 1);
            newTable->column(colName[nr_of_quantities * i_spectral + 1]).write(array_Q, 1);
            newTable->column(colName[nr_of_quantities * i_spectral + 2]).write(array_U, 1);
            newTable->column(colName[nr_of_quantities * i_spectral + 3]).write(array_V, 1);
            newTable->column(colName[nr_of_quantities * i_spectral + 4]).write(array_T, 1);
        }

        // Column density
        valarray<double> array_C(nelements);
        uint i = 0;
        for(uint i_y = 0; i_y < bins_y; i_y++)
            for(uint i_x = 0; i_x < bins_x; i_x++)
            {
                array_C[i] = matrixS[0](i_x, i_y);
                i++;
            }
        newTable->column(colName[nr_of_quantities * nr_of_spectral_bins]).write(array_C, 1);

        // Healpix parameter
        newTable->addKey("PIXTYPE", "HEALPIX", "Pixel algorigthm");
        newTable->addKey("ORDERING", "RING", "Ordering scheme");
        newTable->addKey("INDXSCHM", "IMPLICIT", "Indexing scheme");
        newTable->addKey("NSIDE", uint(sqrt(bins_x / 12.0)), "Resolution Parameter");
        newTable->addKey("FIRSTPIX", 0, "First pixel (0 based)");
        newTable->addKey("LASTPIX", max_cells - 1, "Last pixel (0 based)");
        newTable->addKey("CROTA2", 0, "Rotation Angle (Degrees)");

        if(results_type == RESULTS_RAY)
            newTable->addKey("ETYPE", "thermal emission", "type of emission");
        else if(results_type == RESULTS_FULL)
            newTable->addKey("ETYPE", "thermal emission (+scattering)", "type of emission");
        else
            newTable->addKey("ETYPE", "scattered emission / direct stellar emission", "type of emission");
        newTable->addKey("ID", nr, "detector ID");
        newTable->addKey("OBS_POSITION_X", obs_pos.X(), "x-axis position of observer");
        newTable->addKey("OBS_POSITION_Y", obs_pos.Y(), "y-axis position of observer");
        newTable->addKey("OBS_POSITION_Z", obs_pos.Z(), "z-axis position of observer");
        newTable->addKey("OBS_VELOCITY_X", obs_vel.X(), "velocity of observer in x direction");
        newTable->addKey("OBS_VELOCITY_Y", obs_vel.Y(), "velocity of observer in y direction");
        newTable->addKey("OBS_VELOCITY_Z", obs_vel.Z(), "velocity of observer in z direction");
        newTable->addKey("LONGITUDE_MIN", l_min, "minimum considered galactic longitude");
        newTable->addKey("LONGITUDE_MAX", l_max, "maximum considered galactic longitude");
        newTable->addKey("LATITUDE_MIN", b_min, "minimum considered galactic latitude");
        newTable->addKey("LATITUDE_MAX", b_max, "maximum considered galactic latitude");
        cout << CLR_LINE << flush;
        return true;
    }

    bool writeSyncMap(uint nr)
    {
        cout << CLR_LINE << flush;
        cout << " -> Writing synchrotron plane detector ...  \r" << flush;
        auto_ptr<CCfits::FITS> pFits(0);
        //unique_ptr<CCfits::FITS> pFits;

        try
        {
            char str_tmp[1024];
            char str_end[1024];
            nr++;

#ifdef WINDOWS
            strcpy_s(str_tmp, "polaris_detector_nr%04d");
            sprintf_s(str_end, str_tmp, nr);
#else
            strcpy(str_tmp, "polaris_detector_nr%04d");
            sprintf(str_end, str_tmp, nr);
#endif

            string path_out = path + str_end + ".fits";
            remove(path_out.c_str());

            long naxis = 4;
            long naxes[4] = {uint(bins_x), uint(bins_y), nr_of_spectral_bins, 6 * 3};
            pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
        }
        catch(CCfits::FITS::CantCreate){ return false; }

        long nelements = uint(bins_x) * uint(bins_y);
        if(max_cells != nelements)
        {
            cout << "\nWARNING: Max cells are not equal to bins x bins!" << endl;
            return false;
        }

        vector<long> fpixel(4);
        fpixel[0] = 1;
        fpixel[1] = 1;

        for(uint i_electrons = 0; i_electrons < 3; i_electrons++)
        {
            for(uint i_spectral = 0; i_spectral < nr_of_spectral_bins; i_spectral++)
            {
                fpixel[2] = i_spectral + 1;

                valarray<double> array_I(nelements);
                valarray<double> array_Q(nelements);
                valarray<double> array_U(nelements);
                valarray<double> array_V(nelements);
                valarray<double> array_T(nelements);
                valarray<double> array_S(nelements);

                uint i = 0;
                for(uint i_y = 0; i_y < bins_y; i_y++)
                    for(uint i_x = 0; i_x < bins_x; i_x++)
                    {
                        array_I[i] = matrixI[i_spectral + i_electrons * nr_of_spectral_bins](i_x, i_y);
                        array_Q[i] = matrixQ[i_spectral + i_electrons * nr_of_spectral_bins](i_x, i_y);
                        array_U[i] = matrixU[i_spectral + i_electrons * nr_of_spectral_bins](i_x, i_y);
                        array_V[i] = matrixV[i_spectral + i_electrons * nr_of_spectral_bins](i_x, i_y);
                        array_T[i] = matrixT[i_spectral + i_electrons * nr_of_spectral_bins](i_x, i_y);
                        array_S[i] = matrixS[i_spectral + i_electrons * nr_of_spectral_bins](i_x, i_y);
                        i++;
                    }

                fpixel[3] = 6 * i_electrons + 1;
                pFits->pHDU().write(fpixel, nelements, array_I);
                fpixel[3] = 6 * i_electrons + 2;
                pFits->pHDU().write(fpixel, nelements, array_Q);
                fpixel[3] = 6 * i_electrons + 3;
                pFits->pHDU().write(fpixel, nelements, array_U);
                fpixel[3] = 6 * i_electrons + 4;
                pFits->pHDU().write(fpixel, nelements, array_V);
                fpixel[3] = 6 * i_electrons + 5;
                pFits->pHDU().write(fpixel, nelements, array_T);
                fpixel[3] = 6 * i_electrons + 6;
                pFits->pHDU().write(fpixel, nelements, array_S);
            }
        }

        double bin_width_x, first_pix_val_x, deg_per_pix_x, first_pix_val_deg_x;
        calcCoordinateParameters(sidelength_x, bins_x, map_shift_x, distance, bin_width_x, first_pix_val_x, deg_per_pix_x, first_pix_val_deg_x);
        double bin_width_y, first_pix_val_y, deg_per_pix_y, first_pix_val_deg_y;
        calcCoordinateParameters(sidelength_y, bins_y, map_shift_y, distance, bin_width_y, first_pix_val_y, deg_per_pix_y, first_pix_val_deg_y);

        // Grid
        pFits->pHDU().addKey("CTYPE1", "PARAM", "type of unit 1");
        pFits->pHDU().addKey("CRVAL1", first_pix_val_x, "value of axis 1");
        pFits->pHDU().addKey("CRPIX1", 1, "pixel where CRVAL1 is defined ");
        pFits->pHDU().addKey("CDELT1", bin_width_x, "delta of axis 1");
        pFits->pHDU().addKey("CUNIT1", "m", "unit of axis 1");

        // Alternatively as arc seconds grid
        pFits->pHDU().addKey("CTYPE1A", "RA---TAN", "type of unit 1");
        pFits->pHDU().addKey("CRVAL1A", first_pix_val_deg_x, "value of axis 1");
        pFits->pHDU().addKey("CRPIX1A", bins_x, "pixel where CRVAL1A is defined ");
        pFits->pHDU().addKey("CDELT1A", -deg_per_pix_x, "delta of axis 1");
        pFits->pHDU().addKey("CUNIT1A", "degree", "unit of axis 1");

        // Alternatively as AU grid
        pFits->pHDU().addKey("CTYPE1B", "PARAM", "type of unit 1");
        pFits->pHDU().addKey("CRVAL1B", first_pix_val_x / con_AU, "value of axis 1");
        pFits->pHDU().addKey("CRPIX1B", 1, "pixel where CRVAL1 is defined ");
        pFits->pHDU().addKey("CDELT1B", bin_width_x / con_AU, "delta of axis 1");
        pFits->pHDU().addKey("CUNIT1B", "AU", "unit of axis 1");

        // Alternatively as pc grid
        pFits->pHDU().addKey("CTYPE1C", "PARAM", "type of unit 1");
        pFits->pHDU().addKey("CRVAL1C", first_pix_val_x / con_pc, "value of axis 1");
        pFits->pHDU().addKey("CRPIX1C", 1, "pixel where CRVAL1 is defined ");
        pFits->pHDU().addKey("CDELT1C", bin_width_x / con_pc, "delta of axis 1");
        pFits->pHDU().addKey("CUNIT1C", "pc", "unit of axis 1");

        // Grid
        pFits->pHDU().addKey("CTYPE2", "PARAM", "type of unit 2");
        pFits->pHDU().addKey("CRVAL2", first_pix_val_y, "value of axis 2");
        pFits->pHDU().addKey("CRPIX2", 1, "pixel where CRVAL2 is defined ");
        pFits->pHDU().addKey("CDELT2", bin_width_y, "delta of axis 2");
        pFits->pHDU().addKey("CUNIT2", "m", "unit of axis 2");

        // Alternatively as arc seconds grid
        pFits->pHDU().addKey("CTYPE2A", "DEC--TAN", "type of unit 2");
        pFits->pHDU().addKey("CRVAL2A", first_pix_val_deg_y, "value of axis 2");
        pFits->pHDU().addKey("CRPIX2A", 1, "pixel where CRVAL2A is defined ");
        pFits->pHDU().addKey("CDELT2A", deg_per_pix_y, "delta of axis 2");
        pFits->pHDU().addKey("CUNIT2A", "degree", "unit of axis 2");

        // Alternatively as AU grid
        pFits->pHDU().addKey("CTYPE2B", "PARAM", "type of unit 2");
        pFits->pHDU().addKey("CRVAL2B", first_pix_val_y / con_AU, "value of axis 2");
        pFits->pHDU().addKey("CRPIX2B", 1, "pixel where CRVAL2 is defined ");
        pFits->pHDU().addKey("CDELT2B", bin_width_y / con_AU, "delta of axis 2");
        pFits->pHDU().addKey("CUNIT2B", "AU", "unit of axis 2");

        // Alternatively as pc grid
        pFits->pHDU().addKey("CTYPE2C", "PARAM", "type of unit 2");
        pFits->pHDU().addKey("CRVAL2C", first_pix_val_y / con_pc, "value of axis 2");
        pFits->pHDU().addKey("CRPIX2C", 1, "pixel where CRVAL2 is defined ");
        pFits->pHDU().addKey("CDELT2C", bin_width_y / con_pc, "delta of axis 2");
        pFits->pHDU().addKey("CUNIT2C", "pc", "unit of axis 2");

        // Wavelength IDs
        pFits->pHDU().addKey("CTYPE3", "PARAM", "type of unit 3");
        pFits->pHDU().addKey("CRVAL3", 1, "value of axis 3");
        pFits->pHDU().addKey("CRPIX3", 1, "pixel where CRVAL3 is defined ");
        pFits->pHDU().addKey("CDELT3", 1, "delta of axis 3");
        pFits->pHDU().addKey("CUNIT3", "Wavelength index", "unit of axis 3");

        // Simulated quantities
        pFits->pHDU().addKey("CTYPE4", "PARAM", "type of unit 4");
        pFits->pHDU().addKey("CRVAL4", 1, "value of axis 4");
        pFits->pHDU().addKey("CRPIX4", 1, "pixel where CRVAL4 is defined ");
        pFits->pHDU().addKey("CDELT4", 1, "delta of axis 4");

        pFits->pHDU().addKey("CUNIT4", "I, Q, U, V, [Jy/px], ...", "unit of axis 4");
        pFits->pHDU().addKey("ETYPE", "synchotron emission", "type of emission");

        pFits->pHDU().addKey("ID", nr, "detector ID");

        for(uint i_spectral = 0; i_spectral < nr_of_spectral_bins; i_spectral++)
        {
            char str_1[1024];
            char str_2[1024];
#ifdef WINDOWS
            sprintf_s(str_1, "WAVELENGTH%i", i_spectral + 1);
            sprintf_s(str_2, "value of %i. wavelength", i_spectral + 1);
#else
            sprintf(str_1, "WAVELENGTH%i", i_spectral + 1);
            sprintf(str_2, "value of %i. wavelength", i_spectral + 1);
#endif
            pFits->pHDU().addKey(str_1, wavelength_list_det[i_spectral], str_2);
        }
        pFits->pHDU().addKey("DISTANCE", distance, "distance to object");
        pFits->pHDU().addKey("RAXIS1X", axis1.X(), "rotation axes 1 (x component)");
        pFits->pHDU().addKey("RAXIS1Y", axis1.Y(), "rotation axes 1 (y component)");
        pFits->pHDU().addKey("RAXIS1Z", axis1.Z(), "rotation axes 1 (z component)");
        pFits->pHDU().addKey("RANGLE1", float(180.0 * rot_angle1 / PI), "rotation angle 1 [deg]");
        pFits->pHDU().addKey("RAXIS2X", axis2.X(), "rotation axes 2 (x component)");
        pFits->pHDU().addKey("RAXIS2Y", axis2.Y(), "rotation axes 2 (y component)");
        pFits->pHDU().addKey("RAXIS2Z", axis2.Z(), "rotation axes 2 (z component)");
        pFits->pHDU().addKey("RANGLE2", float(180.0 * rot_angle2 / PI), "rotation angle 2 [deg]");
        pFits->pHDU().addKey("DETGRID", getDetectorGridDescription(), "description of the detector grid");

        return true;
    }

    bool writeSyncHealMap(uint nr)
    {
        cout << CLR_LINE << flush;
        cout << " -> Writing synchrotron healpix detector ...  \r" << flush;
        auto_ptr<CCfits::FITS> pFits(0);
        //unique_ptr<CCfits::FITS> pFits;

        try
        {
            char str_tmp[1024];
            char str_end[1024];
            nr++;

#ifdef WINDOWS
            strcpy_s(str_tmp, "polaris_detector_nr%04d");
            sprintf_s(str_end, str_tmp, nr);
#else
            strcpy(str_tmp, "polaris_detector_nr%04d");
            sprintf(str_end, str_tmp, nr);
#endif

            string path_out = path + str_end + ".fits";
            remove(path_out.c_str());

            long naxis = 1;
            long naxes[1] = {0};
            pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
        }
        catch(CCfits::FITS::CantCreate){  return false; }

        long nelements = uint(bins_x) * uint(bins_y);
        if(max_cells != nelements)
        {
            cout << "\nWARNING: Max cells are not equal to bins x bins!" << endl;
            return false;
        }

        // Init columns
        string newName("HEALPIX_EXTENSION");
        uint nr_of_quantities = 5;
        vector<string> colName(3 * nr_of_quantities * nr_of_spectral_bins + 3, "");
        vector<string> colForm(3 * nr_of_quantities * nr_of_spectral_bins + 3, "");
        vector<string> colUnit(3 * nr_of_quantities * nr_of_spectral_bins + 3, "");

        strlist e_description;
        e_description.push_back("THERMAL ELECTRONS");
        e_description.push_back("COSMIC RAYS");
        e_description.push_back("TOTAL");

        for(uint i_electrons = 0; i_electrons < 3; i_electrons++)
        {
            for(uint i_spectral = 0; i_spectral < nr_of_spectral_bins; i_spectral++)
            {
                uint id = i_electrons * nr_of_quantities * nr_of_spectral_bins + nr_of_quantities * i_spectral;
                char str_1[1024];
    #ifdef WINDOWS
                sprintf_s(str_1, "I_STOKES (WAVELENGTH = %e [m], %s)", wavelength_list_det[i_spectral],
                    e_description[i_electrons].c_str());
    #else
                sprintf(str_1, "I_STOKES (WAVELENGTH = %e [m], %s)", wavelength_list_det[i_spectral],
                    e_description[i_electrons].c_str());
    #endif
                colName[id + 0] = str_1;
    #ifdef WINDOWS
                sprintf_s(str_1, "Q_STOKES (WAVELENGTH = %e [m], %s)", wavelength_list_det[i_spectral],
                    e_description[i_electrons].c_str());
    #else
                sprintf(str_1, "Q_STOKES (WAVELENGTH = %e [m], %s)", wavelength_list_det[i_spectral],
                    e_description[i_electrons].c_str());
    #endif
                colName[id + 1] = str_1;
    #ifdef WINDOWS
                sprintf_s(str_1, "U_STOKES (WAVELENGTH = %e [m], %s)", wavelength_list_det[i_spectral],
                    e_description[i_electrons].c_str());
    #else
                sprintf(str_1, "U_STOKES (WAVELENGTH = %e [m], %s)", wavelength_list_det[i_spectral],
                    e_description[i_electrons].c_str());
    #endif
                colName[id + 2] =  str_1;
    #ifdef WINDOWS
                sprintf_s(str_1, "V_STOKES (WAVELENGTH = %e [m], %s)", wavelength_list_det[i_spectral],
                    e_description[i_electrons].c_str());
    #else
                sprintf(str_1, "V_STOKES (WAVELENGTH = %e [m], %s)", wavelength_list_det[i_spectral],
                    e_description[i_electrons].c_str());
    #endif
                colName[id + 3] = str_1;
    #ifdef WINDOWS
                sprintf_s(str_1, "FARADY ROTATION (WAVELENGTH = %e [m], %s)", wavelength_list_det[i_spectral],
                    e_description[i_electrons].c_str());
    #else
                sprintf(str_1, "FARADY ROTATION (WAVELENGTH = %e [m], %s)", wavelength_list_det[i_spectral],
                    e_description[i_electrons].c_str());
    #endif
                colName[id + 4] = str_1;

                colForm[id + 0] = "D";
                colForm[id + 1] = "D";
                colForm[id + 2] = "D";
                colForm[id + 3] = "D";
                colForm[id + 4] = "D";

                colUnit[id + 0] = "Jy/px";
                colUnit[id + 1] = "Jy/px";
                colUnit[id + 2] = "Jy/px";
                colUnit[id + 3] = "Jy/px";
                colUnit[id + 4] = "rad/m^2";
            }
            char str_1[1024];
#ifdef WINDOWS
            sprintf_s(str_1, "COLUMN_DENSITY (%s)", e_description[i_electrons].c_str());
#else
            sprintf(str_1, "COLUMN_DENSITY (%s)", e_description[i_electrons].c_str());
#endif
            colName[3 * nr_of_quantities * nr_of_spectral_bins + i_electrons] = str_1;
            colForm[3 * nr_of_quantities * nr_of_spectral_bins + i_electrons] = "D";
            colUnit[3 * nr_of_quantities * nr_of_spectral_bins + i_electrons] = "m^-2";
        }

        CCfits::Table* newTable = pFits->addTable(newName, nelements, colName, colForm, colUnit);

        for(uint i_electrons = 0; i_electrons < 3; i_electrons++)
        {
            for(uint i_spectral = 0; i_spectral < nr_of_spectral_bins; i_spectral++)
            {
                valarray<double> array_I(nelements);
                valarray<double> array_Q(nelements);
                valarray<double> array_U(nelements);
                valarray<double> array_V(nelements);
                valarray<double> array_T(nelements);

                uint i = 0;
                for(uint i_y = 0; i_y < bins_y; i_y++)
                    for(uint i_x = 0; i_x < bins_x; i_x++)
                    {
                        array_I[i] = matrixI[i_spectral + i_electrons * nr_of_spectral_bins](i_x, i_y);
                        array_Q[i] = matrixQ[i_spectral + i_electrons * nr_of_spectral_bins](i_x, i_y);
                        array_U[i] = matrixU[i_spectral + i_electrons * nr_of_spectral_bins](i_x, i_y);
                        array_V[i] = matrixV[i_spectral + i_electrons * nr_of_spectral_bins](i_x, i_y);
                        array_T[i] = matrixT[i_spectral + i_electrons * nr_of_spectral_bins](i_x, i_y);
                        i++;
                    }
                uint id = i_electrons * nr_of_quantities * nr_of_spectral_bins + nr_of_quantities * i_spectral;
                newTable->column(colName[id + 0]).write(array_I, 1);
                newTable->column(colName[id + 1]).write(array_Q, 1);
                newTable->column(colName[id + 2]).write(array_U, 1);
                newTable->column(colName[id + 3]).write(array_V, 1);
                newTable->column(colName[id + 4]).write(array_T, 1);
            }

            // Column density
            valarray<double> array_C(nelements);
            uint i = 0;
            for(uint i_y = 0; i_y < bins_y; i_y++)
                for(uint i_x = 0; i_x < bins_x; i_x++)
                {
                    array_C[i] = matrixS[0 + i_electrons * nr_of_spectral_bins](i_x, i_y);
                    i++;
                }
            newTable->column(colName[3 * nr_of_quantities * nr_of_spectral_bins + i_electrons]).write(array_C, 1);
        }

        // Healpix parameter
        newTable->addKey("PIXTYPE", "HEALPIX", "Pixel algorigthm");
        newTable->addKey("ORDERING", "RING", "Ordering scheme");
        newTable->addKey("INDXSCHM", "IMPLICIT", "Indexing scheme");
        newTable->addKey("NSIDE", uint(sqrt(bins_x / 12.0)), "Resolution Parameter");
        newTable->addKey("FIRSTPIX", 0, "First pixel (0 based)");
        newTable->addKey("LASTPIX", max_cells - 1, "Last pixel (0 based)");
        newTable->addKey("CROTA2", 0, "Rotation Angle (Degrees)");

        newTable->addKey("ETYPE", "synchotron emission", "type of emission");
        newTable->addKey("ID", nr, "detector ID");
        newTable->addKey("OBS_POSITION_X", obs_pos.X(), "x-axis position of observer");
        newTable->addKey("OBS_POSITION_Y", obs_pos.Y(), "y-axis position of observer");
        newTable->addKey("OBS_POSITION_Z", obs_pos.Z(), "z-axis position of observer");
        newTable->addKey("OBS_VELOCITY_X", obs_vel.X(), "velocity of observer in x direction");
        newTable->addKey("OBS_VELOCITY_Y", obs_vel.Y(), "velocity of observer in y direction");
        newTable->addKey("OBS_VELOCITY_Z", obs_vel.Z(), "velocity of observer in z direction");
        newTable->addKey("LONGITUDE_MIN", l_min, "minimum considered galactic longitude");
        newTable->addKey("LONGITUDE_MAX", l_max, "maximum considered galactic longitude");
        newTable->addKey("LATITUDE_MIN", b_min, "minimum considered galactic latitude");
        newTable->addKey("LATITUDE_MAX", b_max, "maximum considered galactic latitude");

        cout << CLR_LINE << flush;

        return true;
    }

    bool writeLineSpectrum(CGasMixture * gas, uint i_species, uint i_line)
    {
        cout << CLR_LINE << flush;
        cout << " -> Writing line spectrum ...  \r" << flush;
        auto_ptr<CCfits::FITS> pFits(0);
        //unique_ptr<CCfits::FITS> pFits;

        try
        {
            char str_tmp[1024];
            char str_end[1024];

#ifdef WINDOWS
            strcpy_s(str_tmp, "line_spectrum_species_%04d_line_%04d");
            sprintf_s(str_end, str_tmp, i_species + 1, i_line + 1);
#else
            strcpy(str_tmp, "line_spectrum_species_%04d_line_%04d");
            sprintf(str_end, str_tmp, i_species + 1, i_line + 1);
#endif

            string path_out = path + str_end + ".fits";
            remove(path_out.c_str());

            long naxis = 3;
            long naxes[3] = {nr_of_spectral_bins, 1, 5};
            pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
        }
        catch(CCfits::FITS::CantCreate)
        {
            return false;
        }

        long nelements = uint(nr_of_spectral_bins);

        vector<long> fpixel(3);
        fpixel[0] = 1;
        fpixel[1] = 1;

        std::valarray<double> array_I(nelements);
        std::valarray<double> array_Q(nelements);
        std::valarray<double> array_U(nelements);
        std::valarray<double> array_V(nelements);
        std::valarray<double> array_T(nelements);
        for(uint vch = 0; vch < nr_of_spectral_bins; vch++)
        {
            array_I[vch] = sedI[vch];
            array_Q[vch] = sedQ[vch];
            array_U[vch] = sedU[vch];
            array_V[vch] = sedV[vch];
            array_T[vch] = sedT[vch];

            fpixel[2] = 1;
            pFits->pHDU().write(fpixel, nelements, array_I);
            fpixel[2] = 2;
            pFits->pHDU().write(fpixel, nelements, array_Q);
            fpixel[2] = 3;
            pFits->pHDU().write(fpixel, nelements, array_U);
            fpixel[2] = 4;
            pFits->pHDU().write(fpixel, nelements, array_V);
            fpixel[2] = 5;
            pFits->pHDU().write(fpixel, nelements, array_T);
        }

        // Frequency
        pFits->pHDU().addKey("CTYPE1", "VELO", "type of unit 1");
        pFits->pHDU().addKey("CRVAL1", -max_velocity, "value of axis 1");
        pFits->pHDU().addKey("CRPIX1", 1, "pixel where CRVAL1 is defined ");
        pFits->pHDU().addKey("CDELT1", channel_width, "delta of axis 1");
        pFits->pHDU().addKey("CUNIT1", "m/s", "unit of axis 1");

        // Simulated quantities
        pFits->pHDU().addKey("CTYPE2", "PARAM", "type of unit 2");
        pFits->pHDU().addKey("CRVAL2", 1, "value of axis 2");
        pFits->pHDU().addKey("CRPIX2", 1, "pixel where CRVAL2 is defined ");
        pFits->pHDU().addKey("CDELT2", 1, "delta of axis 2");
        pFits->pHDU().addKey("CUNIT2", "I, Q, U, V [Jy/px], optical depth", "unit of axis 2");

        pFits->pHDU().addKey("GAS_SPECIES", gas->getGasSpeciesName(i_species), "name of the observed gas species");
        pFits->pHDU().addKey("TRANS", i_trans + 1, "transition index number (see leiden database)");
        pFits->pHDU().addKey("LEVEL_UPPER", gas->getUpperTransition(i_species, i_trans),
            "upper energy level index number (see leiden database)");
        pFits->pHDU().addKey("LEVEL_LOWER", gas->getLowerTransition(i_species, i_trans),
            "lower energy level index number (see leiden database)");
        pFits->pHDU().addKey("FREQ", gas->getTransitionFrequency(i_species, i_trans),
            "frequency of the simulated transition");
        pFits->pHDU().addKey("CHANNELS", nr_of_spectral_bins, "number of velocity channels");
        pFits->pHDU().addKey("MAXVEL", max_velocity, "velocity of the velocity channels (-maxvel to maxvel)");
        pFits->pHDU().addKey("ZEEMAN", gas->getZeemanSplitting(i_species),
            "is zeeman splitting in the simulations considered (1=yes/0=no)");
        pFits->pHDU().addKey("DISTANCE", distance, "distance to object");
        pFits->pHDU().addKey("RAXIS1X", axis1.X(), "rotation axes 1 (x component)");
        pFits->pHDU().addKey("RAXIS1Y", axis1.Y(), "rotation axes 1 (y component)");
        pFits->pHDU().addKey("RAXIS1Z", axis1.Z(), "rotation axes 1 (z component)");
        pFits->pHDU().addKey("RANGLE1", float(180.0 * rot_angle1 / PI), "rotation angle 1 [deg]");
        pFits->pHDU().addKey("RAXIS2X", axis2.X(), "rotation axes 2 (x component)");
        pFits->pHDU().addKey("RAXIS2Y", axis2.Y(), "rotation axes 2 (y component)");
        pFits->pHDU().addKey("RAXIS2Z", axis2.Z(), "rotation axes 2 (z component)");
        pFits->pHDU().addKey("RANGLE2", float(180.0 * rot_angle2 / PI), "rotation angle 2 [deg]");

        cout << CLR_LINE << flush;
        return true;
    }

    bool writeVelChannelMaps(CGasMixture * gas, uint i_species, uint i_line)
    {
        cout << CLR_LINE << flush;
        cout << " -> Writing velocity channel maps: 0%     \r" << flush;

        for(uint vch = 0; vch < nr_of_spectral_bins; vch++)
        {
            long naxis = 3;
            long naxes[3] = {uint(bins_x), uint(bins_y), 5};

            std::auto_ptr<CCfits::FITS> pFits(0);
            //std::unique_ptr<CCfits::FITS> pFits;

            try
            {
                char str_tmp[1024];
                char str_end[1024];

    #ifdef WINDOWS
                strcpy_s(str_tmp, "vel_channel_maps_species_%04d_line_%04d_vel_%04d");
                sprintf_s(str_end, str_tmp, i_species + 1, i_line + 1,vch+1);
    #else
                strcpy(str_tmp, "vel_channel_maps_species_%04d_line_%04d_vel_%04d");
                sprintf(str_end, str_tmp, i_species + 1, i_line + 1,vch+1);
    #endif

                string path_out = path + str_end + ".fits";
                remove(path_out.c_str());
                pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
            }
            catch(CCfits::FITS::CantCreate)
            {
                return false;
            }

            long nelements = uint(bins_x) * uint(bins_y);
            if(max_cells != nelements)
            {
                cout << "\nWARNING: Max cells are not equal to bins x bins!" << endl;
                return false;
            }

            vector<long> fpixel(3);
            fpixel[0] = 1;
            fpixel[1] = 1;
            fpixel[2] = 1;


            cout << " -> Writing velocity channel maps: " << int(100.0 * vch / (nr_of_spectral_bins - 1))
                << "%     \r" << flush;

            std::valarray<double> array_I(nelements);
            std::valarray<double> array_Q(nelements);
            std::valarray<double> array_U(nelements);
            std::valarray<double> array_V(nelements);
            std::valarray<double> array_T(nelements);

            uint i = 0;
            for(uint i_y = 0; i_y < bins_y; i_y++)
                for(uint i_x = 0; i_x < bins_x; i_x++)
                {
                    array_I[i] = matrixI[vch](i_x, i_y);
                    array_Q[i] = matrixQ[vch](i_x, i_y);
                    array_U[i] = matrixU[vch](i_x, i_y);
                    array_V[i] = matrixV[vch](i_x, i_y);
                    array_T[i] = matrixT[vch](i_x, i_y);
                    i++;
                }
            fpixel[2] = 1;
            pFits->pHDU().write(fpixel, nelements, array_I);
            fpixel[2] = 2;
            pFits->pHDU().write(fpixel, nelements, array_Q);
            fpixel[2] = 3;
            pFits->pHDU().write(fpixel, nelements, array_U);
            fpixel[2] = 4;
            pFits->pHDU().write(fpixel, nelements, array_V);
            fpixel[2] = 5;
            pFits->pHDU().write(fpixel, nelements, array_T);

            double bin_width_x, first_pix_val_x, deg_per_pix_x, first_pix_val_deg_x;
            calcCoordinateParameters(sidelength_x, bins_x, map_shift_x, distance, bin_width_x, first_pix_val_x, deg_per_pix_x, first_pix_val_deg_x);
            double bin_width_y, first_pix_val_y, deg_per_pix_y, first_pix_val_deg_y;
            calcCoordinateParameters(sidelength_y, bins_y, map_shift_y, distance, bin_width_y, first_pix_val_y, deg_per_pix_y, first_pix_val_deg_y);

            // Grid
            pFits->pHDU().addKey("CTYPE1", "PARAM", "type of unit 1");
            pFits->pHDU().addKey("CRVAL1", first_pix_val_x, "value of axis 1");
            pFits->pHDU().addKey("CRPIX1", 1, "pixel where CRVAL1 is defined ");
            pFits->pHDU().addKey("CDELT1", bin_width_x, "delta of axis 1");
            pFits->pHDU().addKey("CUNIT1", "m", "unit of axis 1");

            // Alternatively as arc seconds grid
            pFits->pHDU().addKey("CTYPE1A", "RA---TAN", "type of unit 1");
            pFits->pHDU().addKey("CRVAL1A", first_pix_val_deg_x, "value of axis 1");
            pFits->pHDU().addKey("CRPIX1A", bins_x, "pixel where CRVAL1A is defined ");
            pFits->pHDU().addKey("CDELT1A", -deg_per_pix_x, "delta of axis 1");
            pFits->pHDU().addKey("CUNIT1A", "degree", "unit of axis 1");

            // Alternatively as AU grid
            pFits->pHDU().addKey("CTYPE1B", "PARAM", "type of unit 1");
            pFits->pHDU().addKey("CRVAL1B", first_pix_val_x / con_AU, "value of axis 1");
            pFits->pHDU().addKey("CRPIX1B", 1, "pixel where CRVAL1 is defined ");
            pFits->pHDU().addKey("CDELT1B", bin_width_x / con_AU, "delta of axis 1");
            pFits->pHDU().addKey("CUNIT1B", "AU", "unit of axis 1");

            // Alternatively as pc grid
            pFits->pHDU().addKey("CTYPE1C", "PARAM", "type of unit 1");
            pFits->pHDU().addKey("CRVAL1C", first_pix_val_x / con_pc, "value of axis 1");
            pFits->pHDU().addKey("CRPIX1C", 1, "pixel where CRVAL1 is defined ");
            pFits->pHDU().addKey("CDELT1C", bin_width_x / con_pc, "delta of axis 1");
            pFits->pHDU().addKey("CUNIT1C", "pc", "unit of axis 1");

            // Grid
            pFits->pHDU().addKey("CTYPE2", "PARAM", "type of unit 2");
            pFits->pHDU().addKey("CRVAL2", first_pix_val_y, "value of axis 2");
            pFits->pHDU().addKey("CRPIX2", 1, "pixel where CRVAL2 is defined ");
            pFits->pHDU().addKey("CDELT2", bin_width_y, "delta of axis 2");
            pFits->pHDU().addKey("CUNIT2", "m", "unit of axis 2");

            // Alternatively as arc seconds grid
            pFits->pHDU().addKey("CTYPE2A", "DEC--TAN", "type of unit 2");
            pFits->pHDU().addKey("CRVAL2A", first_pix_val_deg_y, "value of axis 2");
            pFits->pHDU().addKey("CRPIX2A", 1, "pixel where CRVAL2A is defined ");
            pFits->pHDU().addKey("CDELT2A", deg_per_pix_y, "delta of axis 2");
            pFits->pHDU().addKey("CUNIT2A", "degree", "unit of axis 2");

            // Alternatively as AU grid
            pFits->pHDU().addKey("CTYPE2B", "PARAM", "type of unit 2");
            pFits->pHDU().addKey("CRVAL2B", first_pix_val_y / con_AU, "value of axis 2");
            pFits->pHDU().addKey("CRPIX2B", 1, "pixel where CRVAL2 is defined ");
            pFits->pHDU().addKey("CDELT2B", bin_width_y / con_AU, "delta of axis 2");
            pFits->pHDU().addKey("CUNIT2B", "AU", "unit of axis 2");

            // Alternatively as pc grid
            pFits->pHDU().addKey("CTYPE2C", "PARAM", "type of unit 2");
            pFits->pHDU().addKey("CRVAL2C", first_pix_val_y / con_pc, "value of axis 2");
            pFits->pHDU().addKey("CRPIX2C", 1, "pixel where CRVAL2 is defined ");
            pFits->pHDU().addKey("CDELT2C", bin_width_y / con_pc, "delta of axis 2");
            pFits->pHDU().addKey("CUNIT2C", "pc", "unit of axis 2");

            // Simulated quantities
            pFits->pHDU().addKey("CTYPE3", "PARAM", "type of unit 3");
            pFits->pHDU().addKey("CRVAL3", 1, "value of axis 3");
            pFits->pHDU().addKey("CRPIX3", 1, "pixel where CRVAL3 is defined ");
            pFits->pHDU().addKey("CDELT3", 1, "delta of axis 3");
            pFits->pHDU().addKey("CUNIT3", "I, Q, U, V [Jy/px], optical depth, column density [m^-2]", "unit of axis 3");

            pFits->pHDU().addKey("GAS_SPECIES", gas->getGasSpeciesName(i_species), "name of the observed gas species");
            pFits->pHDU().addKey("TRANS", i_trans + 1, "transition index number (see leiden database)");
            pFits->pHDU().addKey("LEVEL_UPPER", gas->getUpperTransition(i_species, i_trans),
                "upper energy level index number (see leiden database)");
            pFits->pHDU().addKey("LEVEL_LOWER", gas->getLowerTransition(i_species, i_trans),
                "lower energy level index number (see leiden database)");
            pFits->pHDU().addKey("FREQ", gas->getTransitionFrequency(i_species, i_trans),
                "frequency of the simulated transition");
            pFits->pHDU().addKey("VCH", vch, "current velocity channel");
            pFits->pHDU().addKey("CHANNELS", nr_of_spectral_bins, "number of velocity channels");
            pFits->pHDU().addKey("MAXVEL", max_velocity, "velocity of the velocity channels (-maxvel to maxvel)");
            pFits->pHDU().addKey("ZEEMAN", gas->getZeemanSplitting(i_species),
                "is zeeman splitting in the simulations considered (1=yes/0=no)");
            pFits->pHDU().addKey("DISTANCE", distance, "distance to object");
            pFits->pHDU().addKey("RAXIS1X", axis1.X(), "rotation axes 1 (x component)");
            pFits->pHDU().addKey("RAXIS1Y", axis1.Y(), "rotation axes 1 (y component)");
            pFits->pHDU().addKey("RAXIS1Z", axis1.Z(), "rotation axes 1 (z component)");
            pFits->pHDU().addKey("RANGLE1", float(180.0 * rot_angle1 / PI), "rotation angle 1 [deg]");
            pFits->pHDU().addKey("RAXIS2X", axis2.X(), "rotation axes 2 (x component)");
            pFits->pHDU().addKey("RAXIS2Y", axis2.Y(), "rotation axes 2 (y component)");
            pFits->pHDU().addKey("RAXIS2Z", axis2.Z(), "rotation axes 2 (z component)");
            pFits->pHDU().addKey("RANGLE2", float(180.0 * rot_angle2 / PI), "rotation angle 2 [deg]");
        }

        // Writing extra fits file for Zeeman data or column density
        uint nr_of_quantities = 1;
        if(gas->getZeemanSplitting(i_species))
            nr_of_quantities += 4;

        long naxis = 3;
        long naxes[3] = {uint(bins_x), uint(bins_y), nr_of_quantities};

        std::auto_ptr<CCfits::FITS> pFits(0);
        //std::unique_ptr<CCfits::FITS> pFits;

        try
        {
            char str_tmp[1024];
            char str_end[1024];

#ifdef WINDOWS
            strcpy_s(str_tmp, "vel_channel_maps_species_%04d_line_%04d_extra");
            sprintf_s(str_end, str_tmp, i_species + 1, i_line + 1);
#else
            strcpy(str_tmp, "vel_channel_maps_species_%04d_line_%04d_extra");
            sprintf(str_end, str_tmp, i_species + 1, i_line + 1);
#endif

            string path_out = path + str_end + ".fits";
            remove(path_out.c_str());
            pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
        }
        catch(CCfits::FITS::CantCreate)
        {
            return false;
        }

        long nelements = uint(bins_x) * uint(bins_y);
        if(max_cells != nelements)
        {
            cout << "\nWARNING: Max cells are not equal to bins x bins!" << endl;
            return false;
        }

        vector<long> fpixel(3);
        fpixel[0] = 1;
        fpixel[1] = 1;
        fpixel[2] = 1;

        std::valarray<double> array_S(nelements);

        for(uint i_extra = 0; i_extra < nr_of_quantities; i_extra++)
        {
            uint i = 0;
            for(uint i_y = 0; i_y < bins_y; i_y++)
                for(uint i_x = 0; i_x < bins_x; i_x++)
                {
                    if(i_extra == 0)
                    {
                        //column density of the total gas
                        if(gas->getZeemanSplitting(i_species))
                            array_S[i] = matrixS[5](i_x, i_y);
                        else
                            array_S[i] = matrixS[0](i_x, i_y);
                    }
                    else if(i_extra == 1)
                    {
                        //intensity weighted LOS magnetic field
                        if(matrixS[2](i_x, i_y) > 0)
                            array_S[i] = matrixS[0](i_x, i_y) / matrixS[2](i_x, i_y);
                        else
                            array_S[i] = 0;
                    }
                    else if(i_extra == 2)
                    {
                        //intensity weighted total magnetic field
                        if(matrixS[2](i_x, i_y) > 0)
                            array_S[i] = matrixS[1](i_x, i_y) / matrixS[2](i_x, i_y);
                        else
                            array_S[i] = 0;
                    }
                    else if(i_extra == 3)
                    {
                        //density weighted LOS magnetic field
                        if(matrixS[5](i_x, i_y) > 0)
                            array_S[i] = matrixS[3](i_x, i_y) / matrixS[5](i_x, i_y);
                        else
                           array_S[i] = 0;
                    }
                    else if(i_extra == 4)
                    {
                        //density weighted magnetic field
                        if(matrixS[5](i_x, i_y) > 0)
                            array_S[i] = matrixS[4](i_x, i_y) / matrixS[5](i_x, i_y);
                        else
                            array_S[i] = 0;
                    }
                    i++;
                }
            fpixel[2] = i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_S);
        }

        double bin_width_x, first_pix_val_x, deg_per_pix_x, first_pix_val_deg_x;
        calcCoordinateParameters(sidelength_x, bins_x, map_shift_x, distance, bin_width_x, first_pix_val_x, deg_per_pix_x, first_pix_val_deg_x);
        double bin_width_y, first_pix_val_y, deg_per_pix_y, first_pix_val_deg_y;
        calcCoordinateParameters(sidelength_y, bins_y, map_shift_y, distance, bin_width_y, first_pix_val_y, deg_per_pix_y, first_pix_val_deg_y);

        // Grid
        pFits->pHDU().addKey("CTYPE1", "PARAM", "type of unit 1");
        pFits->pHDU().addKey("CRVAL1", first_pix_val_x, "value of axis 1");
        pFits->pHDU().addKey("CRPIX1", 1, "pixel where CRVAL1 is defined ");
        pFits->pHDU().addKey("CDELT1", bin_width_x, "delta of axis 1");
        pFits->pHDU().addKey("CUNIT1", "m", "unit of axis 1");

        // Alternatively as arc seconds grid
        pFits->pHDU().addKey("CTYPE1A", "RA---TAN", "type of unit 1");
        pFits->pHDU().addKey("CRVAL1A", first_pix_val_deg_x, "value of axis 1");
        pFits->pHDU().addKey("CRPIX1A", bins_x, "pixel where CRVAL1A is defined ");
        pFits->pHDU().addKey("CDELT1A", -deg_per_pix_x, "delta of axis 1");
        pFits->pHDU().addKey("CUNIT1A", "degree", "unit of axis 1");

        // Alternatively as AU grid
        pFits->pHDU().addKey("CTYPE1B", "PARAM", "type of unit 1");
        pFits->pHDU().addKey("CRVAL1B", first_pix_val_x / con_AU, "value of axis 1");
        pFits->pHDU().addKey("CRPIX1B", 1, "pixel where CRVAL1 is defined ");
        pFits->pHDU().addKey("CDELT1B", bin_width_x / con_AU, "delta of axis 1");
        pFits->pHDU().addKey("CUNIT1B", "AU", "unit of axis 1");

        // Alternatively as pc grid
        pFits->pHDU().addKey("CTYPE1C", "PARAM", "type of unit 1");
        pFits->pHDU().addKey("CRVAL1C", first_pix_val_x / con_pc, "value of axis 1");
        pFits->pHDU().addKey("CRPIX1C", 1, "pixel where CRVAL1 is defined ");
        pFits->pHDU().addKey("CDELT1C", bin_width_x / con_pc, "delta of axis 1");
        pFits->pHDU().addKey("CUNIT1C", "pc", "unit of axis 1");

        // Grid
        pFits->pHDU().addKey("CTYPE2", "PARAM", "type of unit 2");
        pFits->pHDU().addKey("CRVAL2", first_pix_val_y, "value of axis 2");
        pFits->pHDU().addKey("CRPIX2", 1, "pixel where CRVAL2 is defined ");
        pFits->pHDU().addKey("CDELT2", bin_width_y, "delta of axis 2");
        pFits->pHDU().addKey("CUNIT2", "m", "unit of axis 2");

        // Alternatively as arc seconds grid
        pFits->pHDU().addKey("CTYPE2A", "DEC--TAN", "type of unit 2");
        pFits->pHDU().addKey("CRVAL2A", first_pix_val_deg_y, "value of axis 2");
        pFits->pHDU().addKey("CRPIX2A", 1, "pixel where CRVAL2A is defined ");
        pFits->pHDU().addKey("CDELT2A", deg_per_pix_y, "delta of axis 2");
        pFits->pHDU().addKey("CUNIT2A", "degree", "unit of axis 2");

        // Alternatively as AU grid
        pFits->pHDU().addKey("CTYPE2B", "PARAM", "type of unit 2");
        pFits->pHDU().addKey("CRVAL2B", first_pix_val_y / con_AU, "value of axis 2");
        pFits->pHDU().addKey("CRPIX2B", 1, "pixel where CRVAL2 is defined ");
        pFits->pHDU().addKey("CDELT2B", bin_width_y / con_AU, "delta of axis 2");
        pFits->pHDU().addKey("CUNIT2B", "AU", "unit of axis 2");

        // Alternatively as pc grid
        pFits->pHDU().addKey("CTYPE2C", "PARAM", "type of unit 2");
        pFits->pHDU().addKey("CRVAL2C", first_pix_val_y / con_pc, "value of axis 2");
        pFits->pHDU().addKey("CRPIX2C", 1, "pixel where CRVAL2 is defined ");
        pFits->pHDU().addKey("CDELT2C", bin_width_y / con_pc, "delta of axis 2");
        pFits->pHDU().addKey("CUNIT2C", "pc", "unit of axis 2");

        // Simulated quantities
        pFits->pHDU().addKey("CTYPE3", "PARAM", "type of unit 3");
        pFits->pHDU().addKey("CRVAL3", 1, "value of axis 3");
        pFits->pHDU().addKey("CRPIX3", 1, "pixel where CRVAL3 is defined ");
        pFits->pHDU().addKey("CDELT3", 1, "delta of axis 3");
        if(nr_of_quantities == 5)
            pFits->pHDU().addKey("CUNIT3", "Column density [m^-2], intensity weighted B_LOS, B, "
                "dens. weighted B_LOS, B", "unit of axis 3");
        else
            pFits->pHDU().addKey("CUNIT3", "Column density [m^-2]", "unit of axis 3");

        pFits->pHDU().addKey("GAS_SPECIES", gas->getGasSpeciesName(i_species), "name of the observed gas species");
        pFits->pHDU().addKey("TRANS", i_trans + 1, "transition index number (see leiden database)");
        pFits->pHDU().addKey("LEVEL_UPPER", gas->getUpperTransition(i_species, i_trans),
            "upper energy level index number (see leiden database)");
        pFits->pHDU().addKey("LEVEL_LOWER", gas->getLowerTransition(i_species, i_trans),
            "lower energy level index number (see leiden database)");
        pFits->pHDU().addKey("FREQ", gas->getTransitionFrequency(i_species, i_trans),
            "frequency of the simulated transition");
        pFits->pHDU().addKey("CHANNELS", nr_of_spectral_bins, "number of velocity channels");
        pFits->pHDU().addKey("MAXVEL", max_velocity, "velocity of the velocity channels (-maxvel to maxvel)");
        pFits->pHDU().addKey("ZEEMAN", gas->getZeemanSplitting(i_species),
            "is zeeman splitting in the simulations considered (1=yes/0=no)");
        pFits->pHDU().addKey("DISTANCE", distance, "distance to object");
        pFits->pHDU().addKey("RAXIS1X", axis1.X(), "rotation axes 1 (x component)");
        pFits->pHDU().addKey("RAXIS1Y", axis1.Y(), "rotation axes 1 (y component)");
        pFits->pHDU().addKey("RAXIS1Z", axis1.Z(), "rotation axes 1 (z component)");
        pFits->pHDU().addKey("RANGLE1", float(180.0 * rot_angle1 / PI), "rotation angle 1 [deg]");
        pFits->pHDU().addKey("RAXIS2X", axis2.X(), "rotation axes 2 (x component)");
        pFits->pHDU().addKey("RAXIS2Y", axis2.Y(), "rotation axes 2 (y component)");
        pFits->pHDU().addKey("RAXIS2Z", axis2.Z(), "rotation axes 2 (z component)");
        pFits->pHDU().addKey("RANGLE2", float(180.0 * rot_angle2 / PI), "rotation angle 2 [deg]");

        return true;
    }

    bool writeIntChannelMaps(CGasMixture * gas, uint i_species, uint i_line)
    {
        cout << CLR_LINE << flush;
        cout << " -> Writing integrated velocity channel maps: 0%  \r" << flush;

        long naxis = 3;
        long naxes[3] = {uint(bins_x), uint(bins_y), 6};

        std::auto_ptr<CCfits::FITS> pFits(0);
        //std::unique_ptr<CCfits::FITS> pFits;

        try
        {
            char str_tmp[1024];
            char str_end[1024];

#ifdef WINDOWS
            strcpy_s(str_tmp, "int_channel_map_species_%04d_line_%04d");
            sprintf_s(str_end, str_tmp, i_species + 1, i_line + 1);
#else
            strcpy(str_tmp, "int_channel_map_species_%04d_line_%04d");
            sprintf(str_end, str_tmp, i_species + 1, i_line + 1);
#endif

            string path_out = path + str_end + ".fits";
            remove(path_out.c_str());
            pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
        }
        catch(CCfits::FITS::CantCreate)
        {
            return false;
        }

        long nelements = uint(bins_x) * uint(bins_y);
        if(max_cells != nelements)
        {
            cout << "\nWARNING: Max cells are not equal to bins x bins!" << endl;
            return false;
        }

        vector<long> fpixel(3);
        fpixel[0] = 1;
        fpixel[1] = 1;
        fpixel[2] = 1;

        std::valarray<double> array_I(nelements);
        std::valarray<double> array_Q(nelements);
        std::valarray<double> array_U(nelements);
        std::valarray<double> array_V(nelements);
        std::valarray<double> array_T(nelements);
        std::valarray<double> array_S(nelements);

        uint i = 0;
        for(uint i_y = 0; i_y < bins_y; i_y++)
            for(uint i_x = 0; i_x < bins_x; i_x++)
            {
                for(uint vch = 0; vch < nr_of_spectral_bins; vch++)
                {
                    array_I[i] += matrixI[vch](i_x, i_y) * (2 * max_velocity / nr_of_spectral_bins) * 1.0e-3;
                    array_Q[i] += matrixQ[vch](i_x, i_y) * (2 * max_velocity / nr_of_spectral_bins) * 1.0e-3;
                    array_U[i] += matrixU[vch](i_x, i_y) * (2 * max_velocity / nr_of_spectral_bins) * 1.0e-3;
                    array_V[i] += matrixV[vch](i_x, i_y) * (2 * max_velocity / nr_of_spectral_bins) * 1.0e-3;
                }
                array_T[i] = matrixT[int(nr_of_spectral_bins / 2.0)](i_x, i_y);

                if(matrixS[2](i_x, i_y) == 0)
                {
                    // LOS magnetic field strength
                    array_S[i] = 0;
                }
                else
                    if(matrixS[2](i_x, i_y) > 0)
                        array_S[i] = matrixS[0](i_x, i_y) / matrixS[2](i_x, i_y);
                    else
                        array_S[i] = 0;
                i++;
            }

        fpixel[2] = 1;
        pFits->pHDU().write(fpixel, nelements, array_I);
        fpixel[2] = 2;
        pFits->pHDU().write(fpixel, nelements, array_Q);
        fpixel[2] = 3;
        pFits->pHDU().write(fpixel, nelements, array_U);
        fpixel[2] = 4;
        pFits->pHDU().write(fpixel, nelements, array_V);
        fpixel[2] = 5;
        pFits->pHDU().write(fpixel, nelements, array_T);
        fpixel[2] = 6;
        pFits->pHDU().write(fpixel, nelements, array_S);

        double bin_width_x, first_pix_val_x, deg_per_pix_x, first_pix_val_deg_x;
        calcCoordinateParameters(sidelength_x, bins_x, map_shift_x, distance, bin_width_x, first_pix_val_x, deg_per_pix_x, first_pix_val_deg_x);
        double bin_width_y, first_pix_val_y, deg_per_pix_y, first_pix_val_deg_y;
        calcCoordinateParameters(sidelength_y, bins_y, map_shift_y, distance, bin_width_y, first_pix_val_y, deg_per_pix_y, first_pix_val_deg_y);

        // Grid
        pFits->pHDU().addKey("CTYPE1", "PARAM", "type of unit 1");
        pFits->pHDU().addKey("CRVAL1", first_pix_val_x, "value of axis 1");
        pFits->pHDU().addKey("CRPIX1", 1, "pixel where CRVAL1 is defined ");
        pFits->pHDU().addKey("CDELT1", bin_width_x, "delta of axis 1");
        pFits->pHDU().addKey("CUNIT1", "m", "unit of axis 1");

        // Alternatively as arc seconds grid
        pFits->pHDU().addKey("CTYPE1A", "RA---TAN", "type of unit 1");
        pFits->pHDU().addKey("CRVAL1A", first_pix_val_deg_x, "value of axis 1");
        pFits->pHDU().addKey("CRPIX1A", bins_x, "pixel where CRVAL1A is defined ");
        pFits->pHDU().addKey("CDELT1A", -deg_per_pix_x, "delta of axis 1");
        pFits->pHDU().addKey("CUNIT1A", "degree", "unit of axis 1");

        // Alternatively as AU grid
        pFits->pHDU().addKey("CTYPE1B", "PARAM", "type of unit 1");
        pFits->pHDU().addKey("CRVAL1B", first_pix_val_x / con_AU, "value of axis 1");
        pFits->pHDU().addKey("CRPIX1B", 1, "pixel where CRVAL1 is defined ");
        pFits->pHDU().addKey("CDELT1B", bin_width_x / con_AU, "delta of axis 1");
        pFits->pHDU().addKey("CUNIT1B", "AU", "unit of axis 1");

        // Alternatively as pc grid
        pFits->pHDU().addKey("CTYPE1C", "PARAM", "type of unit 1");
        pFits->pHDU().addKey("CRVAL1C", first_pix_val_x / con_pc, "value of axis 1");
        pFits->pHDU().addKey("CRPIX1C", 1, "pixel where CRVAL1 is defined ");
        pFits->pHDU().addKey("CDELT1C", bin_width_x / con_pc, "delta of axis 1");
        pFits->pHDU().addKey("CUNIT1C", "pc", "unit of axis 1");

        // Grid
        pFits->pHDU().addKey("CTYPE2", "PARAM", "type of unit 2");
        pFits->pHDU().addKey("CRVAL2", first_pix_val_y, "value of axis 2");
        pFits->pHDU().addKey("CRPIX2", 1, "pixel where CRVAL2 is defined ");
        pFits->pHDU().addKey("CDELT2", bin_width_y, "delta of axis 2");
        pFits->pHDU().addKey("CUNIT2", "m", "unit of axis 2");

        // Alternatively as arc seconds grid
        pFits->pHDU().addKey("CTYPE2A", "DEC--TAN", "type of unit 2");
        pFits->pHDU().addKey("CRVAL2A", first_pix_val_deg_y, "value of axis 2");
        pFits->pHDU().addKey("CRPIX2A", 1, "pixel where CRVAL2A is defined ");
        pFits->pHDU().addKey("CDELT2A", deg_per_pix_y, "delta of axis 2");
        pFits->pHDU().addKey("CUNIT2A", "degree", "unit of axis 2");

        // Alternatively as AU grid
        pFits->pHDU().addKey("CTYPE2B", "PARAM", "type of unit 2");
        pFits->pHDU().addKey("CRVAL2B", first_pix_val_y / con_AU, "value of axis 2");
        pFits->pHDU().addKey("CRPIX2B", 1, "pixel where CRVAL2 is defined ");
        pFits->pHDU().addKey("CDELT2B", bin_width_y / con_AU, "delta of axis 2");
        pFits->pHDU().addKey("CUNIT2B", "AU", "unit of axis 2");

        // Alternatively as pc grid
        pFits->pHDU().addKey("CTYPE2C", "PARAM", "type of unit 2");
        pFits->pHDU().addKey("CRVAL2C", first_pix_val_y / con_pc, "value of axis 2");
        pFits->pHDU().addKey("CRPIX2C", 1, "pixel where CRVAL2 is defined ");
        pFits->pHDU().addKey("CDELT2C", bin_width_y / con_pc, "delta of axis 2");
        pFits->pHDU().addKey("CUNIT2C", "pc", "unit of axis 2");

        // Simulated quantities
        pFits->pHDU().addKey("CTYPE3", "PARAM", "type of unit 3");
        pFits->pHDU().addKey("CRVAL3", 1, "value of axis 3");
        pFits->pHDU().addKey("CRPIX3", 1, "pixel where CRVAL3 is defined ");
        pFits->pHDU().addKey("CDELT3", 1, "delta of axis 3");
        pFits->pHDU().addKey("CUNIT3", "I, Q, U, V [Jy/px], optical depth, column density [m^-2]", "unit of axis 3");

        pFits->pHDU().addKey("GAS_SPECIES", gas->getGasSpeciesName(i_species), "name of the observed gas species");
        pFits->pHDU().addKey("TRANS", i_trans + 1, "transition index number (see leiden database)");
        pFits->pHDU().addKey("LEVEL_UPPER", gas->getUpperTransition(i_species, i_trans),
            "upper energy level index number (see leiden database)");
        pFits->pHDU().addKey("LEVEL_LOWER", gas->getLowerTransition(i_species, i_trans),
            "lower energy level index number (see leiden database)");
        pFits->pHDU().addKey("FREQ", gas->getTransitionFrequency(i_species, i_trans),
            "frequency of the simulated transition");
        pFits->pHDU().addKey("CHANNELS", nr_of_spectral_bins, "number of velocity channels");
        pFits->pHDU().addKey("MAXVEL", max_velocity, "velocity of the velocity channels (-maxvel to maxvel)");
        pFits->pHDU().addKey("ZEEMAN", gas->getZeemanSplitting(i_species),
            "is zeeman splitting in the simulations considered (1=yes/0=no)");
        pFits->pHDU().addKey("DISTANCE", distance, "distance to object");
        pFits->pHDU().addKey("RAXIS1X", axis1.X(), "rotation axes 1 (x component)");
        pFits->pHDU().addKey("RAXIS1Y", axis1.Y(), "rotation axes 1 (y component)");
        pFits->pHDU().addKey("RAXIS1Z", axis1.Z(), "rotation axes 1 (z component)");
        pFits->pHDU().addKey("RANGLE1", float(180.0 * rot_angle1 / PI), "rotation angle 1 [deg]");
        pFits->pHDU().addKey("RAXIS2X", axis2.X(), "rotation axes 2 (x component)");
        pFits->pHDU().addKey("RAXIS2Y", axis2.Y(), "rotation axes 2 (y component)");
        pFits->pHDU().addKey("RAXIS2Z", axis2.Z(), "rotation axes 2 (z component)");
        pFits->pHDU().addKey("RANGLE2", float(180.0 * rot_angle2 / PI), "rotation angle 2 [deg]");

        return true;
    }

    bool writeVelChannelHealMaps(CGasMixture * gas, uint i_species, uint i_line)
    {
        cout << CLR_LINE << flush;
        cout << " -> Writing velocity channel maps: 0%     \r" << flush;

        for(uint vch = 0; vch < nr_of_spectral_bins; vch++)
        {
            long naxis = 1;
            long naxes[1] = {0};

            std::auto_ptr<CCfits::FITS> pFits(0);
            //std::unique_ptr<CCfits::FITS> pFits;

            try
            {
                char str_tmp[1024];
                char str_end[1024];

    #ifdef WINDOWS
                strcpy_s(str_tmp, "vel_channel_maps_species_%04d_line_%04d_vel_%04d");
                sprintf_s(str_end, str_tmp, i_species + 1, i_line + 1, vch+1);
    #else
                strcpy(str_tmp, "vel_channel_maps_species_%04d_line_%04d_vel_%04d");
                sprintf(str_end, str_tmp, i_species + 1, i_line + 1, vch+1);
    #endif

                string path_out = path + str_end + ".fits";
                remove(path_out.c_str());
                pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
            }
            catch(CCfits::FITS::CantCreate){ return false; }

            long nelements = uint(bins_x) * uint(bins_y);
            if(max_cells != nelements)
            {
                cout << "\nWARNING: Max cells are not equal to bins x bins!" << endl;
                return false;
            }

            // Init columns
            string newName("HEALPIX_EXTENSION");
            uint nr_of_quantities = 5;
            vector<string> colName(nr_of_quantities, "");
            vector<string> colForm(nr_of_quantities, "");
            vector<string> colUnit(nr_of_quantities, "");

            colName[0] = "I_STOKES";
            colName[1] = "Q_STOKES";
            colName[2] = "U_STOKES";
            colName[3] = "V_STOKES";
            colName[4] = "TAU";

            colForm[0] = "D";
            colForm[1] = "D";
            colForm[2] = "D";
            colForm[3] = "D";
            colForm[4] = "D";

            colUnit[0] = "Jy/px";
            colUnit[1] = "Jy/px";
            colUnit[2] = "Jy/px";
            colUnit[3] = "Jy/px";
            colUnit[4] = "";

            CCfits::Table* newTable = pFits->addTable(newName, nelements, colName, colForm, colUnit);

            cout << " -> Writing velocity channel maps: " << int(100.0 * vch / (nr_of_spectral_bins - 1))
                << "%     \r" << flush;

            valarray<double> array_I(nelements);
            valarray<double> array_Q(nelements);
            valarray<double> array_U(nelements);
            valarray<double> array_V(nelements);
            valarray<double> array_T(nelements);

            uint i = 0;
            for(uint i_y = 0; i_y < bins_y; i_y++)
                for(uint i_x = 0; i_x < bins_x; i_x++)
                {
                    array_I[i] = matrixI[vch](i_x, i_y);
                    array_Q[i] = matrixQ[vch](i_x, i_y);
                    array_U[i] = matrixU[vch](i_x, i_y);
                    array_V[i] = matrixV[vch](i_x, i_y);
                    array_T[i] = matrixT[vch](i_x, i_y);
                    i++;
                }

            newTable->column(colName[0]).write(array_I, 1);
            newTable->column(colName[1]).write(array_Q, 1);
            newTable->column(colName[2]).write(array_U, 1);
            newTable->column(colName[3]).write(array_V, 1);
            newTable->column(colName[4]).write(array_T, 1);

            // Healpix parameter
            newTable->addKey("PIXTYPE", "HEALPIX", "Pixel algorigthm");
            newTable->addKey("ORDERING", "RING", "Ordering scheme");
            newTable->addKey("INDXSCHM", "IMPLICIT", "Indexing scheme");
            newTable->addKey("NSIDE", uint(sqrt(bins_x / 12.0)), "Resolution Parameter");
            newTable->addKey("FIRSTPIX", 0, "First pixel (0 based)");
            newTable->addKey("LASTPIX", max_cells - 1, "Last pixel (0 based)");
            newTable->addKey("CROTA2", 0, "Rotation Angle (Degrees)");

            newTable->addKey("GAS_SPECIES", gas->getGasSpeciesName(i_species), "name of the observed gas species");
            newTable->addKey("TRANS", i_trans + 1, "transition index number (see leiden database)");
            newTable->addKey("LEVEL_UPPER", gas->getUpperTransition(i_species, i_trans),
                "upper energy level index number (see leiden database)");
            newTable->addKey("LEVEL_LOWER", gas->getLowerTransition(i_species, i_trans),
                "lower energy level index number (see leiden database)");
            newTable->addKey("FREQ", gas->getTransitionFrequency(i_species, i_trans),
                "frequency of the simulated transition");
            newTable->addKey("VCH", vch, "current velocity channel");
            newTable->addKey("CHANNELS", nr_of_spectral_bins, "number of velocity channels");
            newTable->addKey("MAXVEL", max_velocity, "velocity of the velocity channels (-maxvel to maxvel)");
            newTable->addKey("ZEEMAN", gas->getZeemanSplitting(i_species),
                "is zeeman splitting in the simulations considered (1=yes/0=no)");
            newTable->addKey("OBS_POSITION_X", obs_pos.X(), "x-axis position of observer");
            newTable->addKey("OBS_POSITION_Y", obs_pos.Y(), "y-axis position of observer");
            newTable->addKey("OBS_POSITION_Z", obs_pos.Z(), "z-axis position of observer");
            newTable->addKey("OBS_VELOCITY_X", obs_vel.X(), "velocity of observer in x direction");
            newTable->addKey("OBS_VELOCITY_Y", obs_vel.Y(), "velocity of observer in y direction");
            newTable->addKey("OBS_VELOCITY_Z", obs_vel.Z(), "velocity of observer in z direction");
            newTable->addKey("LONGITUDE_MIN", l_min, "minimum considered galactic longitude");
            newTable->addKey("LONGITUDE_MAX", l_max, "maximum considered galactic longitude");
            newTable->addKey("LATITUDE_MIN", b_min, "minimum considered galactic latitude");
            newTable->addKey("LATITUDE_MAX", b_max, "maximum considered galactic latitude");
        }

        long naxis = 1;
        long naxes[1] = {0};

        std::auto_ptr<CCfits::FITS> pFits(0);
        //std::unique_ptr<CCfits::FITS> pFits;

        try
        {
            char str_tmp[1024];
            char str_end[1024];

#ifdef WINDOWS
            strcpy_s(str_tmp, "vel_channel_maps_species_%04d_line_%04d_extra");
            sprintf_s(str_end, str_tmp, i_species + 1, i_line + 1);
#else
            strcpy(str_tmp, "vel_channel_maps_species_%04d_line_%04d_extra");
            sprintf(str_end, str_tmp, i_species + 1, i_line + 1);
#endif

            string path_out = path + str_end + ".fits";
            remove(path_out.c_str());
            pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
        }

        catch(CCfits::FITS::CantCreate){ return false; }

        long nelements = uint(bins_x) * uint(bins_y);
        if(max_cells != nelements)
        {
            cout << "\nWARNING: Max cells are not equal to bins x bins!" << endl;
            return false;
        }

        // Init columns
        string newName("HEALPIX_EXTENSION");
        uint nr_of_quantities = 1;
        if(gas->getZeemanSplitting(i_species))
            nr_of_quantities += 4;
        vector<string> colName(nr_of_quantities, "");
        vector<string> colForm(nr_of_quantities, "");
        vector<string> colUnit(nr_of_quantities, "");

        colName[0] = "COLUMN_DENSITY";
        colForm[0] = "D";
        colUnit[0] = "m^-2";

        if(gas->getZeemanSplitting(i_species))
        {
            colName[1] = "INTENSITY WEIGHTED LOS MAGNETIC FIELD";
            colName[2] = "INTENSITY WEIGHTED TOTAL MAGNETIC FIELD";
            colName[3] = "DENSITY WEIGHTED LOS MAGNETIC FIELD";
            colName[4] = "DENSITY WEIGHTED TOTAL MAGNETIC FIELD";

            colForm[1] = "D";
            colForm[2] = "D";
            colForm[3] = "D";
            colForm[4] = "D";

            colUnit[1] = "T";
            colUnit[2] = "T";
            colUnit[3] = "T";
            colUnit[4] = "T";
        }

        CCfits::Table* newTable = pFits->addTable(newName, nelements, colName, colForm, colUnit);

        // Column density / extra data
        for (int i_extra = 0; i_extra < nr_of_quantities; i_extra++)
        {
            valarray<double> array_S(nelements);
            uint i = 0;
            for(uint i_y = 0; i_y < bins_y; i_y++)
                for(uint i_x = 0; i_x < bins_x; i_x++)
                {
                    if(i_extra == 0)
                    {
                        //column density of the total gas
                        if(gas->getZeemanSplitting(i_species))
                            array_S[i] = matrixS[5](i_x, i_y);
                        else
                            array_S[i] = matrixS[0](i_x, i_y);
                    }
                    else if(i_extra == 1)
                    {
                        //intensity weighted LOS magnetic field
                        if(matrixS[2](i_x, i_y) > 0)
                            array_S[i] = matrixS[0](i_x, i_y) / matrixS[2](i_x, i_y);
                        else
                            array_S[i] = 0;
                    }
                    else if(i_extra == 2)
                    {
                        //intensity weighted total magnetic field
                        if(matrixS[2](i_x, i_y) > 0)
                            array_S[i] = matrixS[1](i_x, i_y) / matrixS[2](i_x, i_y);
                        else
                            array_S[i] = 0;
                    }
                    else if(i_extra == 3)
                    {
                        //density weighted LOS magnetic field
                        if(matrixS[5](i_x, i_y) > 0)
                            array_S[i] = matrixS[3](i_x, i_y) / matrixS[5](i_x, i_y);
                        else
                            array_S[i] = 0;
                    }
                    else if(i_extra == 4)
                    {
                        //density weighted magnetic field
                        if(matrixS[5](i_x, i_y) > 0)
                            array_S[i] = matrixS[4](i_x, i_y) / matrixS[5](i_x, i_y);
                        else
                            array_S[i] = 0;
                    }
                    i++;
                }
            newTable->column(colName[i_extra]).write(array_S, 1);
        }

        // Healpix parameter
        newTable->addKey("PIXTYPE", "HEALPIX", "Pixel algorigthm");
        newTable->addKey("ORDERING", "RING", "Ordering scheme");
        newTable->addKey("INDXSCHM", "IMPLICIT", "Indexing scheme");
        newTable->addKey("NSIDE", uint(sqrt(bins_x / 12.0)), "Resolution Parameter");
        newTable->addKey("FIRSTPIX", 0, "First pixel (0 based)");
        newTable->addKey("LASTPIX", max_cells - 1, "Last pixel (0 based)");
        newTable->addKey("CROTA2", 0, "Rotation Angle (Degrees)");

        newTable->addKey("GAS_SPECIES", gas->getGasSpeciesName(i_species), "name of the observed gas species");
        newTable->addKey("TRANS", i_trans + 1, "transition index number (see leiden database)");
        newTable->addKey("LEVEL_UPPER", gas->getUpperTransition(i_species, i_trans),
            "upper energy level index number (see leiden database)");
        newTable->addKey("LEVEL_LOWER", gas->getLowerTransition(i_species, i_trans),
            "lower energy level index number (see leiden database)");
        newTable->addKey("FREQ", gas->getTransitionFrequency(i_species, i_trans),
            "frequency of the simulated transition");
        newTable->addKey("CHANNELS", nr_of_spectral_bins, "number of velocity channels");
        newTable->addKey("MAXVEL", max_velocity, "velocity of the velocity channels (-maxvel to maxvel)");
        newTable->addKey("ZEEMAN", gas->getZeemanSplitting(i_species),
            "is zeeman splitting in the simulations considered (1=yes/0=no)");
        newTable->addKey("OBS_POSITION_X", obs_pos.X(), "x-axis position of observer");
        newTable->addKey("OBS_POSITION_Y", obs_pos.Y(), "y-axis position of observer");
        newTable->addKey("OBS_POSITION_Z", obs_pos.Z(), "z-axis position of observer");
        newTable->addKey("OBS_VELOCITY_X", obs_vel.X(), "velocity of observer in x direction");
        newTable->addKey("OBS_VELOCITY_Y", obs_vel.Y(), "velocity of observer in y direction");
        newTable->addKey("OBS_VELOCITY_Z", obs_vel.Z(), "velocity of observer in z direction");
        newTable->addKey("LONGITUDE_MIN", l_min, "minimum considered galactic longitude");
        newTable->addKey("LONGITUDE_MAX", l_max, "maximum considered galactic longitude");
        newTable->addKey("LATITUDE_MIN", b_min, "minimum considered galactic latitude");
        newTable->addKey("LATITUDE_MAX", b_max, "maximum considered galactic latitude");

        return true;
    }

    bool writeIntVelChannelHealMaps(CGasMixture * gas, uint i_species, uint i_line)
    {
        cout << CLR_LINE << flush;
        cout << " -> Writing integrated velocity channel maps: 0%     \r" << flush;

        long naxis = 1;
        long naxes[1] = {0};

        std::auto_ptr<CCfits::FITS> pFits(0);
        //std::unique_ptr<CCfits::FITS> pFits;

        try
        {
            char str_tmp[1024];
            char str_end[1024];

#ifdef WINDOWS
            strcpy_s(str_tmp, "int_channel_map_species_%04d_line_%04d");
            sprintf_s(str_end, str_tmp, i_species + 1, i_line + 1);
#else
            strcpy(str_tmp, "int_channel_map_species_%04d_line_%04d");
            sprintf(str_end, str_tmp, i_species + 1, i_line + 1);
#endif

            string path_out = path + str_end + ".fits";
            remove(path_out.c_str());
            pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
        }
        catch(CCfits::FITS::CantCreate)
        {
            return false;
        }

        long nelements = uint(bins_x) * uint(bins_y);
        if(max_cells != nelements)
        {
            cout << "\nWARNING: Max cells are not equal to bins x bins!" << endl;
            return false;
        }

        string newName("HEALPIX_EXTENSION");
        vector<string> colName(6, "");
        vector<string> colForm(6, "");
        vector<string> colUnit(6, "");


        colName[0] = "I_STOKES";
        colName[1] = "Q_STOKES";
        colName[2] = "U_STOKES";
        colName[3] = "V_STOKES";
        colName[4] = "TAU";
        colName[5] = "COLUMN_DENSITY";

        colForm[0] = "D";
        colForm[1] = "D";
        colForm[2] = "D";
        colForm[3] = "D";
        colForm[4] = "D";
        colForm[5] = "D";

        colUnit[0] = "Jy/px";
        colUnit[1] = "Jy/px";
        colUnit[2] = "Jy/px";
        colUnit[3] = "Jy/px";
        colUnit[4] = "";
        colUnit[5] = "m^-2";

        CCfits::Table* newTable = pFits->addTable(newName, nelements, colName, colForm, colUnit);

        valarray<double> array_I(nelements);
        valarray<double> array_Q(nelements);
        valarray<double> array_U(nelements);
        valarray<double> array_V(nelements);
        valarray<double> array_T(nelements);
        valarray<double> array_S(nelements);

        uint i = 0;
        for(uint i_y = 0; i_y < bins_y; i_y++)
            for(uint i_x = 0; i_x < bins_x; i_x++)
            {
                for(uint vch = 0; vch < nr_of_spectral_bins; vch++)
                {
                    array_I[i] += matrixI[vch](i_x, i_y) * (2 * max_velocity / nr_of_spectral_bins) * 1.0e-3;
                    array_Q[i] += matrixQ[vch](i_x, i_y) * (2 * max_velocity / nr_of_spectral_bins) * 1.0e-3;
                    array_U[i] += matrixU[vch](i_x, i_y) * (2 * max_velocity / nr_of_spectral_bins) * 1.0e-3;
                    array_V[i] += matrixV[vch](i_x, i_y) * (2 * max_velocity / nr_of_spectral_bins) * 1.0e-3;
                }
                array_T[i] = matrixT[int(nr_of_spectral_bins / 2.0)](i_x, i_y);

                if(matrixS[2](i_x, i_y) == 0.0)
                {
                    // LOS magnetic field strength
                    array_S[i] = 0.0;
                }
                else
                    array_S[i] = matrixS[0](i_x, i_y) / matrixS[2](i_x, i_y);
                i++;
            }

        newTable->column(colName[0]).write(array_I, 1);
        newTable->column(colName[1]).write(array_Q, 1);
        newTable->column(colName[2]).write(array_U, 1);
        newTable->column(colName[3]).write(array_V, 1);
        newTable->column(colName[4]).write(array_T, 1);
        newTable->column(colName[5]).write(array_S, 1);

        // Healpix parameter
        newTable->addKey("PIXTYPE", "HEALPIX", "Pixel algorigthm");
        newTable->addKey("ORDERING", "RING", "Ordering scheme");
        newTable->addKey("INDXSCHM", "IMPLICIT", "Indexing scheme");
        newTable->addKey("NSIDE", uint(sqrt(bins_x / 12.0)), "Resolution Parameter");
        newTable->addKey("FIRSTPIX", 0, "First pixel (0 based)");
        newTable->addKey("LASTPIX", max_cells - 1, "Last pixel (0 based)");
        newTable->addKey("CROTA2", 0, "Rotation Angle (Degrees)");

        newTable->addKey("GAS_SPECIES", gas->getGasSpeciesName(i_species), "name of the observed gas species");
        newTable->addKey("TRANS", i_trans + 1, "transition index number (see leiden database)");
        newTable->addKey("LEVEL_UPPER", gas->getUpperTransition(i_species, i_trans),
            "upper energy level index number (see leiden database)");
        newTable->addKey("LEVEL_LOWER", gas->getLowerTransition(i_species, i_trans),
            "lower energy level index number (see leiden database)");
        newTable->addKey("FREQ", gas->getTransitionFrequency(i_species, i_trans),
            "frequency of the simulated transition");
        newTable->addKey("CHANNELS", nr_of_spectral_bins, "number of velocity channels");
        newTable->addKey("MAXVEL", max_velocity, "velocity of the velocity channels (-maxvel to maxvel)");
        newTable->addKey("ZEEMAN", gas->getZeemanSplitting(i_species),
            "is zeeman splitting in the simulations considered (1=yes/0=no)");
        newTable->addKey("OBS_POSITION_X", obs_pos.X(), "x-axis position of observer");
        newTable->addKey("OBS_POSITION_Y", obs_pos.Y(), "y-axis position of observer");
        newTable->addKey("OBS_POSITION_Z", obs_pos.Z(), "z-axis position of observer");
        newTable->addKey("OBS_VELOCITY_X", obs_vel.X(), "velocity of observer in x direction");
        newTable->addKey("OBS_VELOCITY_Y", obs_vel.Y(), "velocity of observer in y direction");
        newTable->addKey("OBS_VELOCITY_Z", obs_vel.Z(), "velocity of observer in z direction");
        newTable->addKey("LONGITUDE_MIN", l_min, "minimum considered galactic longitude");
        newTable->addKey("LONGITUDE_MAX", l_max, "maximum considered galactic longitude");
        newTable->addKey("LATITUDE_MIN", b_min, "minimum considered galactic latitude");
        newTable->addKey("LATITUDE_MAX", b_max, "maximum considered galactic latitude");

        return true;
    }

    void calcCoordinateParameters(double sidelength, double bins, double map_shift, double distance,
        double & bin_width, double & first_pix_val, double & deg_per_pix, double & first_pix_val_deg)
    {
        bin_width = sidelength / bins;

        first_pix_val = -sidelength / 2.0 + map_shift;
        first_pix_val += (bin_width / 2.0);

        deg_per_pix = bin_width / con_AU;
        deg_per_pix /= (distance / con_pc);
        deg_per_pix /= 3600.0;

        first_pix_val_deg = -sidelength / (2.0 * con_AU) + map_shift / con_AU;
        first_pix_val_deg /= (distance / con_pc);
        first_pix_val_deg /= 3600.0;
        first_pix_val_deg += (deg_per_pix / 2.0);
    }

    void calcVelocityChannels(uint _nr_of_spectral_bins, double _max_velocity)
    {
        nr_of_spectral_bins = _nr_of_spectral_bins;
        max_velocity = _max_velocity;

        if(nr_of_spectral_bins > 1)
            channel_width = (2.0 * max_velocity) / (nr_of_spectral_bins - 1);
        else
            channel_width = (2.0 * max_velocity);

        velocity_channel.resize(nr_of_spectral_bins);
        if(nr_of_spectral_bins > 1)
            for(uint i = 0; i < nr_of_spectral_bins; i++)
                velocity_channel[i] = 2 * (float) i / ((float) nr_of_spectral_bins - 1) * max_velocity - max_velocity;
        else if(nr_of_spectral_bins == 1)
            velocity_channel[0] = 0;
    }

    string getPath()
    {
        return path;
    }

    double getLamMin()
    {
        return lam_min;
    }

    double getLamMax()
    {
        return lam_max;
    }

    uint getDetectorWavelengthID(double wavelength)
    {
        dlist::iterator it = find(wavelength_list_det.begin(), wavelength_list_det.end(), wavelength);
        if(it != wavelength_list_det.end())
            return it - wavelength_list_det.begin();
        else
            return MAX_UINT;
    }

    double getWavelength(uint i_wave)
    {
        return wavelength_list_det[i_wave];
    }

    uint getNrOfSpectralBins()
    {
        return nr_of_spectral_bins;
    }

    double getChannelWidth()
    {
        return channel_width;
    }

    double getVelocityChannel(uint vch)
    {
        return velocity_channel[vch];
    }

    string getDetectorGridDescription()
    {
        switch(detector_id)
        {
            case DET_PLANE:
                return "Plane / Cartesian background grid";
                break;

            case DET_SPHER:
                return "Spherical / Healpix background grid (all-sky map)";
                break;

            case DET_POLAR:
                return  "Polar background grid";
                break;

            case DET_SLICE:
                return "Cartesian grid as a slice through the grid (y-direction along LOS)";
                break;

            default:
                return "Photon packages launched from emission sources (Monte-Carlo)";
                break;
        }
    }

private:
    double cos_acceptance_angle;
    double rot_angle1, rot_angle2, distance;
    double l_min, l_max, b_min, b_max;
    double lam_min, lam_max;
    double sidelength_x, sidelength_y;
    double map_shift_x, map_shift_y;
    double channel_width, max_velocity;
    Vector3D ex, ey, ez;
    Vector3D obs_pos, obs_vel;
    string path;
    uint ID, detector_id;
    uint bins_x, bins_y;
    uint max_cells;
    uint nr_of_spectral_bins;
    uint nr_velocity_channels;
    uint i_trans;
    Matrix2D *matrixI, *matrixQ, *matrixU, *matrixV, *matrixT, *matrixS;
    double *sedI, *sedQ, *sedU, *sedV, *sedT;
    dlist wavelength_list_det;
    dlist velocity_channel;
    Vector3D axis1, axis2, pos;
};

#endif
