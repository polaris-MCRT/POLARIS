/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CDETECTOR_H
#define CDETECTOR_H

#include "GasSpecies.hpp"
#include "GasMixture.hpp"
#include "Stokes.hpp"
#include "Typedefs.hpp"
#include "Vector3D.hpp"
#include "OPIATE.hpp"
#include <CCfits/CCfits>

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
        rad_bubble=0;

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
        matrixS1 = 0;
        matrixS2 = 0;
        matrixS3 = 0;
        matrixS4 = 0;
        
        special_param = 1;

        sedI = 0;
        sedQ = 0;
        sedU = 0;
        sedV = 0;
        sedT = 0;
        sedS = 0;

        lam_min = 0;
        lam_max = 0;

        nr_extra = 1;

        nr_spectral_bins = 1;
        channel_width = 0;

        cos_acceptance_angle = 0;
        alignment = ALIG_RND;
        processing_method = 0;

        w1_I = 0;
        w1_Q = 0;
        w1_U = 0;
        w1_V = 0;
        w1_PI = 0;

        w2_I = 0;
        w2_Q = 0;
        w2_U = 0;
        w2_V = 0;
        w2_PI = 0;

        w3_I = 0;
        w3_Q = 0;
        w3_U = 0;
        w4_V = 0;
        w3_PI = 0;

        w4_I = 0;
        w4_Q = 0;
        w4_U = 0;
        w4_V = 0;
        w4_PI = 0;

        N_photon = 0;
    }

    // Detector for dust, synchrotron, and free free
    // Plane detector
    CDetector(uint _detector_id,
              string _path,
              uint _bins_x,
              uint _bins_y,
              uint _id,
              double _sidelength_x,
              double _sidelength_y,
              double _map_shift_x,
              double _map_shift_y,
              double _distance,
              double _l_min,
              double _l_max,
              uint _nr_spectral_bins,
              uint _nr_extra,
              uint _special_param,
              uint _alignment = ALIG_RND)
    {
        detector_id = _detector_id;

        rot_angle1 = 0;
        rot_angle2 = 0;
        rad_bubble = 0;
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
        alignment = _alignment;

        distance = _distance;

        lam_min = _l_min;
        lam_max = _l_max;
        nr_spectral_bins = _nr_spectral_bins;
        nr_extra = _nr_extra;
        
        special_param = _special_param;

        wavelength_list_det.resize(nr_spectral_bins);
        if(USE_LOG_SPACING)
            {CMathFunctions::LogList(lam_min, lam_max, wavelength_list_det, 10);}
        else
            {CMathFunctions::LinearList(lam_min, lam_max, wavelength_list_det);}

        sedI = new double[nr_extra * nr_spectral_bins];
        sedQ = new double[nr_extra * nr_spectral_bins];
        sedU = new double[nr_extra * nr_spectral_bins];
        sedV = new double[nr_extra * nr_spectral_bins];
        sedT = new double[nr_extra * nr_spectral_bins];
        sedS = new double[nr_extra * nr_spectral_bins];

        matrixI = new Matrix2D[nr_extra * nr_spectral_bins];
        matrixQ = new Matrix2D[nr_extra * nr_spectral_bins];
        matrixU = new Matrix2D[nr_extra * nr_spectral_bins];
        matrixV = new Matrix2D[nr_extra * nr_spectral_bins];
        
        matrixT = new Matrix2D[nr_extra * nr_spectral_bins];
        
        matrixS1 = new Matrix2D[nr_extra * nr_spectral_bins];
        
        matrixS2 = 0;
        matrixS3 = 0;
        matrixS4 = 0;
        
        if(special_param==3)
        {
            matrixS2 = new Matrix2D[nr_extra * nr_spectral_bins];
            matrixS3 = new Matrix2D[nr_extra * nr_spectral_bins];
        }
        

        
        N_photon = new Matrix2D[nr_extra * nr_spectral_bins];

        for(uint i_extra = 0; i_extra < nr_extra; i_extra++)
        {
            for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
            {
                matrixI[i_spectral + i_extra * nr_spectral_bins].resize(bins_x, bins_y);
                matrixQ[i_spectral + i_extra * nr_spectral_bins].resize(bins_x, bins_y);
                matrixU[i_spectral + i_extra * nr_spectral_bins].resize(bins_x, bins_y);
                matrixV[i_spectral + i_extra * nr_spectral_bins].resize(bins_x, bins_y);
                matrixT[i_spectral + i_extra * nr_spectral_bins].resize(bins_x, bins_y);
                
                if(matrixS1!=0)
                    matrixS1[i_spectral + i_extra * nr_spectral_bins].resize(bins_x, bins_y);
                
                if(matrixS2!=0)
                    matrixS2[i_spectral + i_extra * nr_spectral_bins].resize(bins_x, bins_y);
                
                if(matrixS3!=0)
                    matrixS3[i_spectral + i_extra * nr_spectral_bins].resize(bins_x, bins_y);
                
                if(matrixS4!=0)
                    matrixS4[i_spectral + i_extra * nr_spectral_bins].resize(bins_x, bins_y);
                
                N_photon[i_spectral + i_extra * nr_spectral_bins].resize(bins_x, bins_y);

                sedI[i_spectral + i_extra * nr_spectral_bins] = 0;
                sedQ[i_spectral + i_extra * nr_spectral_bins] = 0;
                sedU[i_spectral + i_extra * nr_spectral_bins] = 0;
                sedV[i_spectral + i_extra * nr_spectral_bins] = 0;
                sedT[i_spectral + i_extra * nr_spectral_bins] = 0;
                sedS[i_spectral + i_extra * nr_spectral_bins] = 0;
            }
        }
    }

    // Spherical dust detector
    CDetector(string _path,
              uint _bins,
              uint _id,
              Vector3D obs_pos,
              double _sidelength,
              double _l_min,
              double _l_max,
              double _rad_bubble,
              uint _nr_of_spectral_bins,
              uint _nr_extra,
              uint _special_param,
              uint _alignment,
              bool compact)
    {
        detector_id = DET_SPHER;

        rot_angle1 = 0;
        rot_angle2 = 0;
        rad_bubble = _rad_bubble;
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
        alignment = _alignment;
        
        special_param = _special_param;

        lam_min = _l_min;
        lam_max = _l_max;
        nr_spectral_bins = _nr_of_spectral_bins;
        nr_extra = _nr_extra;

        wavelength_list_det.resize(nr_spectral_bins);
        if(USE_LOG_SPACING)
            {CMathFunctions::LogList(lam_min, lam_max, wavelength_list_det, 10);}
        else
            {CMathFunctions::LinearList(lam_min, lam_max, wavelength_list_det);}

        sedI = new double[nr_extra * nr_spectral_bins];
        sedQ = new double[nr_extra * nr_spectral_bins];
        sedU = new double[nr_extra * nr_spectral_bins];
        sedV = new double[nr_extra * nr_spectral_bins];
        sedT = new double[nr_extra * nr_spectral_bins];
        sedS = new double[nr_extra * nr_spectral_bins];

        matrixI = new Matrix2D[nr_extra * nr_spectral_bins];
        
        matrixS2 = 0;
        matrixS3 = 0;
        matrixS4 = 0;
        
        if(compact)
        {
            matrixQ = 0;
            matrixU = 0;
            matrixV = 0;
            matrixT = 0;
            
            matrixS1 = 0;

            N_photon = 0;
        }
        else
        {
            matrixT = new Matrix2D[nr_extra * nr_spectral_bins];
            matrixS1 = new Matrix2D[nr_extra * nr_spectral_bins];
            
            if(special_param==3)
            {
                matrixS2 = new Matrix2D[nr_extra * nr_spectral_bins];
                matrixS3 = new Matrix2D[nr_extra * nr_spectral_bins];
            }
            
            N_photon = new Matrix2D[nr_extra * nr_spectral_bins];
            
            if(alignment!=ALIG_RND)
            {
                matrixQ = new Matrix2D[nr_extra * nr_spectral_bins];
                matrixU = new Matrix2D[nr_extra * nr_spectral_bins];
                matrixV = new Matrix2D[nr_extra * nr_spectral_bins];
            }
            else
            {
                matrixQ = 0;
                matrixU = 0;
                matrixV = 0;
            }
        }
        
        for(uint i_extra = 0; i_extra < nr_extra; i_extra++)
        {
            for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
            {
                
                matrixI[i_spectral + i_extra * nr_spectral_bins].resize(bins_x, bins_y);
                
                if(matrixQ!=0)
                    matrixQ[i_spectral + i_extra * nr_spectral_bins].resize(bins_x, bins_y);
                
                if(matrixU!=0)
                    matrixU[i_spectral + i_extra * nr_spectral_bins].resize(bins_x, bins_y);
                
                if(matrixV!=0)
                    matrixV[i_spectral + i_extra * nr_spectral_bins].resize(bins_x, bins_y);
                
                if(matrixT!=0)
                    matrixT[i_spectral + i_extra * nr_spectral_bins].resize(bins_x, bins_y);
                
                if(matrixS1!=0)
                    matrixS1[i_spectral + i_extra * nr_spectral_bins].resize(bins_x, bins_y);
                
                if(matrixS2!=0)
                    matrixS2[i_spectral + i_extra * nr_spectral_bins].resize(bins_x, bins_y);
                
                if(matrixS3!=0)
                    matrixS3[i_spectral + i_extra * nr_spectral_bins].resize(bins_x, bins_y);
                
                if(matrixS4!=0)
                    matrixS4[i_spectral + i_extra * nr_spectral_bins].resize(bins_x, bins_y);
                
                if(N_photon!=0)
                    N_photon[i_spectral + i_extra * nr_spectral_bins].resize(bins_x, bins_y);

                sedI[i_spectral + i_extra * nr_spectral_bins] = 0;
                sedQ[i_spectral + i_extra * nr_spectral_bins] = 0;
                sedU[i_spectral + i_extra * nr_spectral_bins] = 0;
                sedV[i_spectral + i_extra * nr_spectral_bins] = 0;
                sedT[i_spectral + i_extra * nr_spectral_bins] = 0;
                sedS[i_spectral + i_extra * nr_spectral_bins] = 0;
            }
        }
    }
    // End detector dust, synchrotron, and free-free

    // Detector for line
    // Plane detector
    CDetector(uint _detector_id,
              string _path,
              uint _bins_x,
              uint _bins_y,
              uint _id,
              double _sidelength_x,
              double _sidelength_y,
              double _map_shift_x,
              double _map_shift_y,
              double _distance,
              uint _i_trans,
              uint _nr_spectral_bins,
              double _min_velocity,
              double _max_velocity,
              bool hasZeeman)
    {
        detector_id = _detector_id;

        rot_angle1 = 0;
        rot_angle2 = 0;
        rad_bubble = 0;
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

        nr_extra = 1;
        special_param=1;

        calcVelocityChannels(_nr_spectral_bins, _min_velocity,  _max_velocity);

        sidelength_x = _sidelength_x;
        sidelength_y = _sidelength_y;

        map_shift_x = _map_shift_x;
        map_shift_y = _map_shift_y;

        wavelength_list_det.resize(nr_spectral_bins);

        i_trans = _i_trans;

        cos_acceptance_angle = 0;

        sedI = new double[nr_spectral_bins];
        sedQ = new double[nr_spectral_bins];
        sedU = new double[nr_spectral_bins];
        sedV = new double[nr_spectral_bins];
        sedT = new double[nr_spectral_bins];
        sedS = new double[nr_spectral_bins];

        matrixI = new Matrix2D[nr_spectral_bins];
        matrixQ = new Matrix2D[nr_spectral_bins];
        matrixU = new Matrix2D[nr_spectral_bins];
        matrixV = new Matrix2D[nr_spectral_bins];
        matrixT = new Matrix2D[nr_spectral_bins];
        matrixS1 = new Matrix2D[nr_spectral_bins];
        N_photon = new Matrix2D[nr_spectral_bins];
        
        matrixS2 = 0;
        matrixS3 = 0;
        matrixS4 = 0;

        for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
        {
            matrixI[i_spectral].resize(bins_x, bins_y);
            matrixQ[i_spectral].resize(bins_x, bins_y);
            matrixU[i_spectral].resize(bins_x, bins_y);
            matrixV[i_spectral].resize(bins_x, bins_y);
            matrixT[i_spectral].resize(bins_x, bins_y);
            matrixS1[i_spectral].resize(bins_x, bins_y);
            N_photon[i_spectral].resize(bins_x, bins_y);

            sedI[i_spectral] = 0;
            sedQ[i_spectral] = 0;
            sedU[i_spectral] = 0;
            sedV[i_spectral] = 0;
            sedT[i_spectral] = 0;
            sedS[i_spectral] = 0;
        }

        lam_min = 0;
        lam_max = 0;
    }

    // Spherical line detector
    CDetector(string _path,
              uint _bins,
              uint _id,
              Vector3D obs_pos,
              double _sidelength,
              uint _i_trans,
              uint _nr_spectral_bins,
              double _min_velocity,
              double _max_velocity,
              bool hasZeeman)
    {
        detector_id = DET_SPHER;

        rot_angle1 = 0;
        rot_angle2 = 0;
        rad_bubble = 0; //tbd: implementation for line detectors!
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
        
        special_param = 1;

        calcVelocityChannels(_nr_spectral_bins, _min_velocity, _max_velocity);

        sidelength_x = _sidelength;
        sidelength_y = _sidelength;

        map_shift_x = 0;
        map_shift_y = 0;

        pos = obs_pos;

        wavelength_list_det.resize(nr_spectral_bins);

        i_trans = _i_trans;

        cos_acceptance_angle = 0;

        sedI = new double[nr_spectral_bins];
        sedQ = new double[nr_spectral_bins];
        sedU = new double[nr_spectral_bins];
        sedV = new double[nr_spectral_bins];
        sedT = new double[nr_spectral_bins];
        sedS = new double[nr_spectral_bins];

        matrixI = new Matrix2D[nr_spectral_bins];
        
        if(hasZeeman)
        {
            matrixQ = new Matrix2D[nr_spectral_bins];
            matrixU = new Matrix2D[nr_spectral_bins];
            matrixV = new Matrix2D[nr_spectral_bins]; 
        }
        else
        {
            matrixQ = 0;
            matrixU = 0;
            matrixV = 0; 
        }
        
        matrixT = new Matrix2D[nr_spectral_bins];
        matrixS1 = new Matrix2D[nr_spectral_bins];
        
        matrixS2 = 0;
        matrixS3 = 0;
        matrixS4 = 0;
        
        N_photon = 0;

        for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
        {
            matrixI[i_spectral].resize(bins_x, bins_y);
            
            if(hasZeeman)
            {
                matrixQ[i_spectral].resize(bins_x, bins_y);
                matrixU[i_spectral].resize(bins_x, bins_y);
                matrixV[i_spectral].resize(bins_x, bins_y);
            }
            
            matrixT[i_spectral].resize(bins_x, bins_y);
            
            matrixS1[i_spectral].resize(bins_x, bins_y);
            
            //N_photon[i_spectral].resize(bins_x, bins_y);

            sedI[i_spectral] = 0;
            sedQ[i_spectral] = 0;
            sedU[i_spectral] = 0;
            sedV[i_spectral] = 0;
            sedT[i_spectral] = 0;
            sedS[i_spectral] = 0;
        }

        lam_min = 0;
        lam_max = 0;
    }
    // End detector line

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
            
        if(matrixS1 != 0)
            delete[] matrixS1;
        if(matrixS2 != 0)
            delete[] matrixS2;
        if(matrixS3 != 0)
            delete[] matrixS3;
        if(matrixS4 != 0)
            delete[] matrixS4;

        if(w1_I != 0)
            delete[] w1_I;
        if(w1_Q != 0)
            delete[] w1_Q;
        if(w1_U != 0)
            delete[] w1_U;
        if(w1_V != 0)
            delete[] w1_V;
        if(w1_PI != 0)
            delete[] w1_PI;

        if(w2_I != 0)
            delete[] w2_I;
        if(w2_Q != 0)
            delete[] w2_Q;
        if(w2_U != 0)
            delete[] w2_U;
        if(w2_V != 0)
            delete[] w2_V;
        if(w2_PI != 0)
            delete[] w2_PI;

        if(w3_I != 0)
            delete[] w3_I;
        if(w3_Q != 0)
            delete[] w3_Q;
        if(w3_U != 0)
            delete[] w3_U;
        if(w3_V != 0)
            delete[] w3_V;
        if(w3_PI != 0)
            delete[] w3_PI;

        if(w4_I != 0)
            delete[] w4_I;
        if(w4_Q != 0)
            delete[] w4_Q;
        if(w4_U != 0)
            delete[] w4_U;
        if(w4_V != 0)
            delete[] w4_V;
        if(w4_PI != 0)
            delete[] w4_PI;

        if(N_photon != 0)
            delete[] N_photon;

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
        if(sedS != 0)
            delete[] sedS;
    }

    void init(string _path,
              uint _bins_x,
              uint _bins_y,
              double _sidelength_x,
              double _sidelength_y,
              double _map_shift_x,
              double _map_shift_y,
              double _distance,
              double _l_min,
              double _l_max,
              uint _nr_spectral_bins);

    Vector3D getDirection();

    double getDistance();

    double getAcceptanceAngle();

    Vector3D getEX();

    Vector3D getEY();

    Vector3D getEZ();

    void setOrientation(Vector3D n1, Vector3D n2, double _rot_angle1, double _rot_angle2);

    void setObsPosition(Vector3D _obs_pos,
                        Vector3D _obs_vel,
                        double _l_min,
                        double _l_max,
                        double _b_min,
                        double _b_max);

    void setAcceptanceAngle(double acceptance_angle);

    void setProcessingMethod(uint val);

    void addToRaytracingSedDetector(const photon_package & pp);

    void addToRaytracingDetector(const photon_package & pp, uint pos_id = MAX_UINT);

    void addToMonteCarloDetector(const photon_package & pp, uint i_det_spectral, uint radiation_type);

    double calc_R(uint i_spectral, uint x, uint y, uint quantity);

    double calc_VOV(uint i_spectral, uint x, uint y, uint quantity);

    bool writeMap(uint nr, uint results_type);

    bool writeMapStats(uint nr, uint results_type);

    bool writeSed(uint nr, uint results_type);

    bool writeHealMaps(uint nr, uint results_type);
    bool writeHealMapsTiny(uint nr, uint results_type);
    
    bool writeFreeFreeMap(uint nr);
    bool writeFreeFreeHealMap(uint nr);

    bool writeSyncMap(uint nr);

    bool writeSyncHealMap(uint nr);

    bool writeLineSpectrum(CGasMixture * gas, uint i_species, uint i_line);

    bool writeOPIATESpectrum(COpiateDataBase *op, uint det_id);

    bool writeOPIATEVelChannelMaps(COpiateDataBase * op, uint det_id);

    bool writeOPIATEIntChannelMaps(COpiateDataBase * op, uint det_id);

    bool writeVelChannelMaps(CGasMixture * gas, uint i_species, uint i_line);

    bool writeIntChannelMaps(CGasMixture * gas, uint i_species, uint i_line);

    bool writeOPIATEVelChannelHealMaps(COpiateDataBase * op, uint det_id);

    bool writeVelChannelHealMaps(CGasMixture * gas, uint i_species, uint i_line);
    
    bool writeVelChannelHealMapsTiny(CGasMixture * gas, uint i_species, uint i_line);

    bool writeOPIATEIntVelChannelHealMaps(COpiateDataBase * op, uint det_id);

    bool writeIntVelChannelHealMaps(CGasMixture * gas, uint i_species, uint i_line);

    void calcCoordinateParameters(double sidelength,
                                  double bins,
                                  double map_shift,
                                  double distance,
                                  double & bin_width,
                                  double & first_pix_val,
                                  double & deg_per_pix,
                                  double & first_pix_val_deg);

    void calcVelocityChannels(uint _nr_spectral_bins, double _min_velocity, double _max_velocity);

    string getPath();

    double getLamMin();

    double getLamMax();

    uint getDetectorWavelengthID(double wavelength);

    double getWavelength(uint i_wave);

    uint getNrOfSpectralBins();

    double getChannelWidth();

    double getVelocityChannel(uint i_spectral);

    string getDetectorGridDescription();

    string getAlignmentDescription();

private:
    double cos_acceptance_angle;
    double rot_angle1, rot_angle2, distance;
    double rad_bubble;
    double l_min, l_max, b_min, b_max;
    double lam_min, lam_max;
    double sidelength_x, sidelength_y;
    double map_shift_x, map_shift_y;
    double channel_width, min_velocity, max_velocity1;
    Vector3D ex, ey, ez;
    Vector3D obs_pos, obs_vel;
    string path;
    uint ID, detector_id;
    uint bins_x, bins_y;
    uint special_param;
    uint max_cells;
    uint nr_spectral_bins;
    uint nr_velocity_channels;
    uint i_trans;
    uint nr_extra;
    uint alignment;
    uint processing_method;
    
    //Stokes components
    Matrix2D *matrixI, *matrixQ, *matrixU, *matrixV;
    
    //special parameters
    Matrix2D *matrixT, *matrixS1, *matrixS2, *matrixS3, *matrixS4;
    
    Matrix2D *w1_I, *w2_I, *w3_I, *w4_I;
    Matrix2D *w1_Q, *w2_Q, *w3_Q, *w4_Q;
    Matrix2D *w1_U, *w2_U, *w3_U, *w4_U;
    Matrix2D *w1_V, *w2_V, *w3_V, *w4_V;
    Matrix2D *w1_PI, *w2_PI, *w3_PI, *w4_PI;
    Matrix2D *N_photon;
    double *sedI, *sedQ, *sedU, *sedV, *sedT, *sedS;
    dlist wavelength_list_det;
    dlist velocity_channel;
    Vector3D axis1, axis2, pos;
};

#endif /* CDETECTOR_H */
