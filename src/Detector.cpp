/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "Detector.hpp"

// Dust scattering detector
void CDetector::init(string _path,
                     uint _bins_x,
                     uint _bins_y,
                     double _sidelength_x,
                     double _sidelength_y,
                     double _map_shift_x,
                     double _map_shift_y,
                     double _distance,
                     double _l_min,
                     double _l_max,
                     uint _nr_spectral_bins)
{
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

    ID = 0;
    bins_x = _bins_x;
    bins_y = _bins_y;
    max_cells = bins_x * bins_y;

    sidelength_x = _sidelength_x;
    sidelength_y = _sidelength_y;

    map_shift_x = _map_shift_x;
    map_shift_y = _map_shift_y;

    channel_width = 0;
    i_trans = 0, cos_acceptance_angle = 0;

    distance = _distance;

    lam_min = _l_min;
    lam_max = _l_max;
    nr_spectral_bins = _nr_spectral_bins;

    wavelength_list_det.resize(nr_spectral_bins);
    if(USE_LOG_SPACING)
        {CMathFunctions::LogList(lam_min, lam_max, wavelength_list_det, 10);}
    else
        {CMathFunctions::LinearList(lam_min, lam_max, wavelength_list_det);}

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

    w1_I = new Matrix2D[nr_spectral_bins];
    w1_Q = new Matrix2D[nr_spectral_bins];
    w1_U = new Matrix2D[nr_spectral_bins];
    w1_V = new Matrix2D[nr_spectral_bins];
    w1_PI = new Matrix2D[nr_spectral_bins];

    w2_I = new Matrix2D[nr_spectral_bins];
    w2_Q = new Matrix2D[nr_spectral_bins];
    w2_U = new Matrix2D[nr_spectral_bins];
    w2_V = new Matrix2D[nr_spectral_bins];
    w2_PI = new Matrix2D[nr_spectral_bins];

    w3_I = new Matrix2D[nr_spectral_bins];
    w3_Q = new Matrix2D[nr_spectral_bins];
    w3_U = new Matrix2D[nr_spectral_bins];
    w3_V = new Matrix2D[nr_spectral_bins];
    w3_PI = new Matrix2D[nr_spectral_bins];

    w4_I = new Matrix2D[nr_spectral_bins];
    w4_Q = new Matrix2D[nr_spectral_bins];
    w4_U = new Matrix2D[nr_spectral_bins];
    w4_V = new Matrix2D[nr_spectral_bins];
    w4_PI = new Matrix2D[nr_spectral_bins];

    N_photon = new Matrix2D[nr_spectral_bins];

    for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
    {
        matrixI[i_spectral].resize(bins_x, bins_y);
        matrixQ[i_spectral].resize(bins_x, bins_y);
        matrixU[i_spectral].resize(bins_x, bins_y);
        matrixV[i_spectral].resize(bins_x, bins_y);
        matrixT[i_spectral].resize(bins_x, bins_y);
        matrixS1[i_spectral].resize(bins_x, bins_y);

        sedI[i_spectral] = 0;
        sedQ[i_spectral] = 0;
        sedU[i_spectral] = 0;
        sedV[i_spectral] = 0;
        sedT[i_spectral] = 0;
        sedS[i_spectral] = 0;

        w1_I[i_spectral].resize(bins_x, bins_y);
        w1_Q[i_spectral].resize(bins_x, bins_y);
        w1_U[i_spectral].resize(bins_x, bins_y);
        w1_V[i_spectral].resize(bins_x, bins_y);
        w1_PI[i_spectral].resize(bins_x, bins_y);

        w2_I[i_spectral].resize(bins_x, bins_y);
        w2_Q[i_spectral].resize(bins_x, bins_y);
        w2_U[i_spectral].resize(bins_x, bins_y);
        w2_V[i_spectral].resize(bins_x, bins_y);
        w2_PI[i_spectral].resize(bins_x, bins_y);

        w3_I[i_spectral].resize(bins_x, bins_y);
        w3_Q[i_spectral].resize(bins_x, bins_y);
        w3_U[i_spectral].resize(bins_x, bins_y);
        w3_V[i_spectral].resize(bins_x, bins_y);
        w3_PI[i_spectral].resize(bins_x, bins_y);

        w4_I[i_spectral].resize(bins_x, bins_y);
        w4_Q[i_spectral].resize(bins_x, bins_y);
        w4_U[i_spectral].resize(bins_x, bins_y);
        w4_V[i_spectral].resize(bins_x, bins_y);
        w4_PI[i_spectral].resize(bins_x, bins_y);

        N_photon[i_spectral].resize(bins_x, bins_y);
    }
}

Vector3D CDetector::getDirection()
{
    return ez;
}

double CDetector::getDistance()
{
    return distance;
}

double CDetector::getAcceptanceAngle()
{
    return cos_acceptance_angle;
}

StokesVector CDetector::getStokesVector(int64_t i_pix, int64_t pos)
{
    StokesVector st;
    
    double I = 0;
    double Q = 0;
    double U = 0;
    double V = 0;
    double T = 0;
    double S1 = 0;
    
    if(matrixI!=0)
        I = matrixI[pos](i_pix, 0);
    
    if(matrixQ!=0)
        Q = matrixQ[pos](i_pix, 0);
    
    if(matrixU!=0)
        U = matrixU[pos](i_pix, 0);
    
    if(matrixV!=0)
        V = matrixV[pos](i_pix, 0);
    
    if(matrixT!=0)
        T = matrixT[pos](i_pix, 0);

    if(matrixS1!=0)
        S1 = matrixS1[pos](i_pix, 0);
    
    st.set(I,Q,U,V,T);
    st.setSp1(S1);
    
    return st;
}

Vector3D CDetector::getEX()
{
    return ex;
}

Vector3D CDetector::getEY()
{
    return ey;
}

Vector3D CDetector::getEZ()
{
    return ez;
}

void CDetector::setOrientation(Vector3D n1, Vector3D n2, double _rot_angle1, double _rot_angle2)
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

void CDetector::setObsPosition(Vector3D _obs_pos,
                    Vector3D _obs_vel,
                    double _l_min,
                    double _l_max,
                    double _b_min,
                    double _b_max)
{
    obs_pos = _obs_pos;
    obs_vel = _obs_vel;

    l_min = _l_min;
    l_max = _l_max;

    b_min = _b_min;
    b_max = _b_max;
}

void CDetector::setAcceptanceAngle(double acceptance_angle)
{
    cos_acceptance_angle = cos(acceptance_angle);
}

void CDetector::setProcessingMethod(uint val)
{
    processing_method = val;
}

void CDetector::addToRaytracingSedDetector(const photon_package & pp)
{
    StokesVector st = pp.getStokesVector();
    uint i_spectral = pp.getSpectralID();
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

void CDetector::addToRaytracingDetector(const photon_package & pp, uint pos_id)
{
    // Get Stokes vector from photon package
    StokesVector st = pp.getStokesVector();

    // Get spectral index from photon package
    uint i_spectral = pp.getSpectralID();

    if(pos_id == MAX_UINT)
    {
        Vector3D pos = pp.getPosition();

        double float_x = (pos.X() + 0.5 * sidelength_x - map_shift_x) / sidelength_x * double(bins_x);
        if (float_x < 0.) // integer casting leads to zeros from negative and positive side which results in too much flux in the lowest row
            return; 

        uint x = uint(float_x);
        if(x >= bins_x) 
            return;

        double float_y = (pos.Y() + 0.5 * sidelength_y - map_shift_y) / sidelength_y * double(bins_y);
        if (float_y < 0.) // integer casting leads to zeros from negative and positive side which results in too much flux in the lowest column
            return;

        uint y = uint(float_y);
        if(y >= bins_y)
            return;

        if(matrixI!=0)
            matrixI[i_spectral].addValue(x, y, st.I());
        
        if(matrixQ!=0)
            matrixQ[i_spectral].addValue(x, y, st.Q());
        
        if(matrixU!=0)
            matrixU[i_spectral].addValue(x, y, st.U());
        
        if(matrixV!=0)
            matrixV[i_spectral].addValue(x, y, st.V());
        
        if(matrixT!=0)
            matrixT[i_spectral].addValue(x, y, st.T());

        if(N_photon!=0)
            N_photon[i_spectral].addValue(x, y, 1);
        
        if(matrixS1!=0)
            matrixS1[i_spectral].addValue(x, y, st.Sp1());
        
        if(matrixS2!=0)
            matrixS2[i_spectral].addValue(x, y, st.Sp2());
        
        if(matrixS3!=0)
            matrixS3[i_spectral].addValue(x, y, st.Sp3());
        
        if(matrixS4!=0)
            matrixS4[i_spectral].addValue(x, y, st.Sp4());
        
    }
    else if(pos_id == MAX_UINT-1)
    {
        Vector3D pos = pp.getPosition();
        
        int x = pos.X();
        int y = pos.Y();

        if(x < 0)
            return; 
        
        if(x >= bins_x) 
            return;

        if (y < 0)
            return;

        if(y >= bins_y)
            return;

        if(matrixI!=0)
            matrixI[i_spectral].addValue(x, y, st.I());
        
        if(matrixQ!=0)
            matrixQ[i_spectral].addValue(x, y, st.Q());
        
        if(matrixU!=0)
            matrixU[i_spectral].addValue(x, y, st.U());
        
        if(matrixV!=0)
            matrixV[i_spectral].addValue(x, y, st.V());
        
        if(matrixT!=0)
            matrixT[i_spectral].addValue(x, y, st.T());

        if(N_photon!=0)
            N_photon[i_spectral].addValue(x, y, 1);
        
        if(matrixS1!=0)
            matrixS1[i_spectral].addValue(x, y, st.Sp1());
        
        if(matrixS2!=0)
            matrixS2[i_spectral].addValue(x, y, st.Sp2());
        
        if(matrixS3!=0)
            matrixS3[i_spectral].addValue(x, y, st.Sp3());
        
        if(matrixS4!=0)
            matrixS4[i_spectral].addValue(x, y, st.Sp4());
        
    }
    else
    {
        if(matrixI!=0)
            matrixI[i_spectral].addValue(pos_id, st.I());
        
        if(matrixQ!=0)
            matrixQ[i_spectral].addValue(pos_id, st.Q());
            
        if(matrixU!=0)    
            matrixU[i_spectral].addValue(pos_id, st.U());
            
        if(matrixV!=0)    
            matrixV[i_spectral].addValue(pos_id, st.V());
        
        if(matrixT!=0)
            matrixT[i_spectral].addValue(pos_id, st.T());
        
        if(matrixS1!=0)
            matrixS1[i_spectral].addValue(pos_id, st.Sp1());
        
        if(matrixS2!=0)
            matrixS2[i_spectral].addValue(pos_id, st.Sp2());
        
        if(matrixS3!=0)
            matrixS3[i_spectral].addValue(pos_id, st.Sp3());
        
        if(matrixS4!=0)
            matrixS4[i_spectral].addValue(pos_id, st.Sp4());
    }
}

void CDetector::addToMonteCarloDetector(const photon_package & pp, uint i_det_spectral, uint radiation_type)
{
    Vector3D pos = pp.getPosition();
    Vector3D dir = pp.getDirection();

    if(dir.length() == 0)
        return;

    double lx = pos * ex;
    double ly = pos * ey;

    uint x = uint((lx + 0.5 * sidelength_x - map_shift_x) / sidelength_x * double(bins_x));

    if(x < 0 || x >= bins_x)
        return;

    uint y = uint((ly + 0.5 * sidelength_y - map_shift_y) / sidelength_y * double(bins_y));

    if(y < 0 || y >= bins_y)
        return;

    StokesVector st = pp.getStokesVector();

    matrixI[i_det_spectral].addValue(x, y, st.I());
    matrixQ[i_det_spectral].addValue(x, y, st.Q());
    matrixU[i_det_spectral].addValue(x, y, st.U());
    matrixV[i_det_spectral].addValue(x, y, st.V());

    if(radiation_type == SCATTERED_DUST)
    {
        w1_I[i_det_spectral].addValue(x, y, st.I());
        w1_Q[i_det_spectral].addValue(x, y, st.Q());
        w1_U[i_det_spectral].addValue(x, y, st.U());
        w1_V[i_det_spectral].addValue(x, y, st.V());
        w1_PI[i_det_spectral].addValue(x, y, st.tPol());

        w2_I[i_det_spectral].addValue(x, y, pow(st.I(),2));
        w2_Q[i_det_spectral].addValue(x, y, pow(st.Q(),2));
        w2_U[i_det_spectral].addValue(x, y, pow(st.U(),2));
        w2_V[i_det_spectral].addValue(x, y, pow(st.V(),2));
        w2_PI[i_det_spectral].addValue(x, y, pow(st.tPol(),2));

        w3_I[i_det_spectral].addValue(x, y, pow(st.I(),3));
        w3_Q[i_det_spectral].addValue(x, y, pow(st.Q(),3));
        w3_U[i_det_spectral].addValue(x, y, pow(st.U(),3));
        w3_V[i_det_spectral].addValue(x, y, pow(st.V(),3));
        w3_PI[i_det_spectral].addValue(x, y, pow(st.tPol(),3));

        w4_I[i_det_spectral].addValue(x, y, pow(st.I(),4));
        w4_Q[i_det_spectral].addValue(x, y, pow(st.Q(),4));
        w4_U[i_det_spectral].addValue(x, y, pow(st.U(),4));
        w4_V[i_det_spectral].addValue(x, y, pow(st.V(),4));
        w4_PI[i_det_spectral].addValue(x, y, pow(st.tPol(),4));

        N_photon[i_det_spectral].addValue(x, y, 1);
    }

    // Add to SED
    if(sedI != 0)
    {
        #pragma omp atomic update
        sedI[i_det_spectral] += st.I();
        #pragma omp atomic update
        sedQ[i_det_spectral] += st.Q();
        #pragma omp atomic update
        sedU[i_det_spectral] += st.U();
        #pragma omp atomic update
        sedV[i_det_spectral] += st.V();
    }

    if(radiation_type == DIRECT_STAR)
    {
        matrixT[i_det_spectral].addValue(x, y, st.I());
        if(sedI != 0)
        {
            #pragma omp atomic update
            sedT[i_det_spectral] += st.I();
        }
    }
    else
    {
        matrixS1[i_det_spectral].addValue(x, y, st.I());
        if(sedI != 0)
        {
            #pragma omp atomic update
            sedS[i_det_spectral] += st.I();
        }
    }
}

double CDetector::calc_R(uint i_spectral, uint x, uint y, uint quantity)
{
    double w1, w2, N;
    double R = 1;

    N = N_photon[i_spectral](x,y);
    if(N >= 2)
    {
        switch(quantity)
        {
            case 0:
            {
                w1 = w1_I[i_spectral](x,y);
                w2 = w2_I[i_spectral](x,y);
                break;
            }
            case 1:
            {
                w1 = w1_Q[i_spectral](x,y);
                w2 = w2_Q[i_spectral](x,y);
                break;
            }
            case 2:
            {
                w1 = w1_U[i_spectral](x,y);
                w2 = w2_U[i_spectral](x,y);
                break;
            }
            case 3:
            {
                w1 = w1_V[i_spectral](x,y);
                w2 = w2_V[i_spectral](x,y);
                break;
            }
            case 4:
            {
                w1 = w1_PI[i_spectral](x,y);
                w2 = w2_PI[i_spectral](x,y);
                break;
            }
        }

        if(w1 != 0)
        {
            R = w2 / pow(w1,2) - 1/N;
            R = sqrt(R);
        }
    }
    return R;
}

double CDetector::calc_VOV(uint i_spectral, uint x, uint y, uint quantity)
{
    double w1, w2, w3, w4, N;
    double VOV = 1;

    N = N_photon[i_spectral](x,y);
    if(N >= 2)
    {
        switch(quantity)
        {
            case 0:
            {
                w1 = w1_I[i_spectral](x,y);
                w2 = w2_I[i_spectral](x,y);
                w3 = w3_I[i_spectral](x,y);
                w4 = w4_I[i_spectral](x,y);
                break;
            }
            case 1:
            {
                w1 = w1_Q[i_spectral](x,y);
                w2 = w2_Q[i_spectral](x,y);
                w3 = w3_Q[i_spectral](x,y);
                w4 = w4_Q[i_spectral](x,y);
                break;
            }
            case 2:
            {
                w1 = w1_U[i_spectral](x,y);
                w2 = w2_U[i_spectral](x,y);
                w3 = w3_U[i_spectral](x,y);
                w4 = w4_U[i_spectral](x,y);
                break;
            }
            case 3:
            {
                w1 = w1_V[i_spectral](x,y);
                w2 = w2_V[i_spectral](x,y);
                w3 = w3_V[i_spectral](x,y);
                w4 = w4_V[i_spectral](x,y);
                break;
            }
            case 4:
            {
                w1 = w1_PI[i_spectral](x,y);
                w2 = w2_PI[i_spectral](x,y);
                w3 = w3_PI[i_spectral](x,y);
                w4 = w4_PI[i_spectral](x,y);
                break;
            }
        }

        VOV = w4 + w1 / N * ( -3*w1 * pow(w1 / N, 2) + 6*w1*w2 / N - 4*w3 );
        VOV /= pow(w2 - pow(w1,2) / N, 2);
        VOV -= 1 / N;
    }
    return VOV;
}

bool CDetector::writeProjMap(uint nr, uint results_type)
{
    cout << CLR_LINE << flush;
    cout << " -> Writing map ...  \r" << flush;
    // auto_ptr<CCfits::FITS> pFits(0);
    unique_ptr<CCfits::FITS> pFits;
    
    uint nr_quantities=0;
    
    if (matrixI != 0)
        nr_quantities++;

    if (matrixQ != 0)
        nr_quantities++;

    if (matrixU != 0)
        nr_quantities++;

    if (matrixV != 0)
        nr_quantities++;

    if (matrixT != 0)
        nr_quantities++;


    if (matrixS1 != 0)
        nr_quantities++;

    if (matrixS2 != 0)
        nr_quantities++;

    if (matrixS3 != 0)
        nr_quantities++;

    try
    {
        char str_tmp[1024];
        char str_end[1024];
        nr++;

#ifdef WINDOWS
        strcpy_s(str_tmp, "polaris_detector_nr%04d");
        sprintf_s(str_end, str_tmp, nr);
#else
        strcpy(str_tmp, "proj_polaris_detector_nr%04d");
        sprintf(str_end, str_tmp, nr);
#endif

        string path_out = path + str_end + FITS_COMPRESS_EXT;
        remove(path_out.c_str());

        long naxis = 4;
        long naxes[4] = { bins_x, bins_y, nr_spectral_bins, nr_quantities * nr_extra };

        pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
    }
    catch(CCfits::FITS::CantCreate)
    {
        return false;
    }

    long nelements = bins_x * bins_y;
    if(max_cells != nelements)
    {
        cout << WARNING_LINE << "Max cells are not equal to bins x bins!" << endl;
        return false;
    }

    vector<long> fpixel(4);
    fpixel[0] = 1;
    fpixel[1] = 1;

    for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
    {
        fpixel[2] = i_spectral + 1;
        for(uint i_extra = 0; i_extra < nr_extra; i_extra++)
        {
            /*valarray<double> array_I(nelements);
            valarray<double> array_Q(nelements);
            valarray<double> array_U(nelements);
            valarray<double> array_V(nelements);
            valarray<double> array_T(nelements);
            valarray<double> array_S1(nelements);
            valarray<double> array_S2(nelements);
            valarray<double> array_S3(nelements);*/
            
            valarray<double> array_I = {};
            valarray<double> array_Q = {};
            valarray<double> array_U = {};
            valarray<double> array_V = {};
            valarray<double> array_T = {};
            
            valarray<double> array_S1 = {};
            valarray<double> array_S2 = {};
            valarray<double> array_S3 = {};

            if (matrixI != 0)
                array_I = valarray<double>(double(0), long(nelements));

            if (matrixQ != 0)
                array_Q = valarray<double>(double(0), long(nelements));

            if (matrixU != 0)
                array_U = valarray<double>(double(0), long(nelements));

            if (matrixV != 0)
                array_V = valarray<double>(double(0), long(nelements));

            if (matrixT != 0)
                array_T = valarray<double>(double(0), long(nelements));
            
            
            if (matrixS1 != 0)
                array_S1 = valarray<double>(double(0), long(nelements));

            if (matrixS2 != 0)
                array_S2 = valarray<double>(double(0), long(nelements));

            if (matrixS3 != 0)
                array_S3 = valarray<double>(double(0), long(nelements));            
            
            uint i = 0;
            for(uint i_y = 0; i_y < bins_y; i_y++)
            {
                for(uint i_x = 0; i_x < bins_x; i_x++)
                {
                    if(matrixI != 0)
                        array_I[i] = matrixI[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    
                    if(matrixQ != 0)
                        array_Q[i] = matrixQ[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    
                    if(matrixU != 0)
                        array_U[i] = matrixU[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    
                    if(matrixV != 0)
                        array_V[i] = matrixV[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    
                    if(matrixT != 0)
                        array_T[i] = matrixT[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    
                    if(matrixS1 != 0)
                        array_S1[i] = matrixS1[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    
                    if(special_param==3)
                    {
                        array_S2[i] = matrixS2[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                        array_S3[i] = matrixS3[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    }
                    
                    double N_ph = 1;
                    
                    if(N_photon != 0)
                        N_ph = N_photon[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);

                    if((results_type == RESULTS_RAY || results_type == RESULTS_FULL) &&
                        detector_id == DET_POLAR && processing_method == NEAREST && N_ph > 0)
                    {
                        if(matrixT != 0)
                            array_T[i] /= N_ph;
                        
                        if(matrixS1!=0)
                            array_S1[i] /= N_ph;
                    }
                    
                    i++;
                }
            }
            
            uint off_counter = 0;
            
            if(matrixI != 0)
            {
                fpixel[3] = off_counter * nr_extra + i_extra + 1;
                pFits->pHDU().write(fpixel, nelements, array_I);
                off_counter++;
            }
            
            if(matrixQ != 0)
            {
                fpixel[3] = off_counter * nr_extra + i_extra + 1;
                pFits->pHDU().write(fpixel, nelements, array_Q);
                off_counter++;
            }
            
            if(matrixU != 0)
            {
                fpixel[3] = off_counter * nr_extra + i_extra + 1;
                pFits->pHDU().write(fpixel, nelements, array_U);
                off_counter++;
            }
            
            if(matrixV != 0)
            {
                fpixel[3] = off_counter * nr_extra + i_extra + 1;
                pFits->pHDU().write(fpixel, nelements, array_V);
                off_counter++;
            }
            
            if(matrixT != 0)
            {
                fpixel[3] = off_counter * nr_extra + i_extra + 1;
                pFits->pHDU().write(fpixel, nelements, array_T);
                off_counter++;
            }
            
            if(matrixS1 != 0)
            {
                fpixel[3] = off_counter * nr_extra + i_extra + 1;
                pFits->pHDU().write(fpixel, nelements, array_S1);
                off_counter++;
            }
            
            if(special_param==3)
            {
                if(matrixS2 != 0)
                {
                    fpixel[3] = off_counter * nr_extra + i_extra + 1;
                    pFits->pHDU().write(fpixel, nelements, array_S2);
                    off_counter++;
                }

                if(matrixS3 != 0)
                {
                    fpixel[3] = off_counter * nr_extra + i_extra + 1;
                    pFits->pHDU().write(fpixel, nelements, array_S3);
                    off_counter++;
                }
            }
        }
    }

    double bin_width_x, first_pix_val_x, deg_per_pix_x, first_pix_val_deg_x;
    
    calcCoordinateParameters(sidelength_x,
                                bins_x,
                                map_shift_x,
                                distance,
                                bin_width_x,
                                first_pix_val_x,
                                deg_per_pix_x,
                                first_pix_val_deg_x);
    
    double bin_width_y, first_pix_val_y, deg_per_pix_y, first_pix_val_deg_y;
    
    calcCoordinateParameters(sidelength_y,
                                bins_y,
                                map_shift_y,
                                distance,
                                bin_width_y,
                                first_pix_val_y,
                                deg_per_pix_y,
                                first_pix_val_deg_y);

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
        pFits->pHDU().addKey(
            "CUNIT4", "I, Q, U, V [Jy/px], optical depth, column density [m^-2]", "unit of axis 4");
        pFits->pHDU().addKey("ETYPE", "thermal emission", "type of emission");
    }
    else if(results_type == RESULTS_FULL)
    {
        pFits->pHDU().addKey(
            "CUNIT4", "I, Q, U, V [Jy/px], optical depth, column density [m^-2]", "unit of axis 4");
        pFits->pHDU().addKey("ETYPE", "thermal emission (+scattering)", "type of emission");
    }
    else if(results_type == RESULTS_MC)
    {
        pFits->pHDU().addKey("CUNIT4", "I, Q, U, V, I_direct, I_scat [Jy/px]", "unit of axis 4");
        pFits->pHDU().addKey("ETYPE", "scattered emission / direct stellar emission", "type of emission");
    }
    else if(nr_extra > 1)
    {
        pFits->pHDU().addKey(
            "CUNIT4", "I, Q, U, V [Jy/px], optical depth, column density [m^-2]", "unit of axis 4");
        pFits->pHDU().addKey(
            "ETYPE", "Each Stokes entry for FULL / SCAT / THERMAL / STOCHASTIC", "type of emission");
    }
    pFits->pHDU().addKey("ID", nr, "detector ID");

    for(uint i_extra = 0; i_extra < nr_extra; i_extra++)
    {
        for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
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
    string alignment_descr = getAlignmentDescription();
    if(alignment_descr != "")
        pFits->pHDU().addKey("ALIGNMENT", alignment_descr, "alignment method of dust grains");

    cout << CLR_LINE << flush;
    return true;
}

bool CDetector::writeMap(uint nr, uint results_type)
{
    cout << CLR_LINE << flush;
    cout << " -> Writing map ...  \r" << flush;
    // auto_ptr<CCfits::FITS> pFits(0);
    unique_ptr<CCfits::FITS> pFits;
    
    uint nr_quantities=0;
    
    if (matrixI != 0)
        nr_quantities++;

    if (matrixQ != 0)
        nr_quantities++;

    if (matrixU != 0)
        nr_quantities++;

    if (matrixV != 0)
        nr_quantities++;

    if (matrixT != 0)
        nr_quantities++;


    if (matrixS1 != 0)
        nr_quantities++;

    if (matrixS2 != 0)
        nr_quantities++;

    if (matrixS3 != 0)
        nr_quantities++;

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

        string path_out = path + str_end + FITS_COMPRESS_EXT;
        remove(path_out.c_str());

        long naxis = 4;
        long naxes[4] = { bins_x, bins_y, nr_spectral_bins, nr_quantities * nr_extra };

        pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
    }
    catch(CCfits::FITS::CantCreate)
    {
        return false;
    }

    long nelements = bins_x * bins_y;
    if(max_cells != nelements)
    {
        cout << WARNING_LINE << "Max cells are not equal to bins x bins!" << endl;
        return false;
    }

    vector<long> fpixel(4);
    fpixel[0] = 1;
    fpixel[1] = 1;

    for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
    {
        fpixel[2] = i_spectral + 1;
        for(uint i_extra = 0; i_extra < nr_extra; i_extra++)
        {
            /*valarray<double> array_I(nelements);
            valarray<double> array_Q(nelements);
            valarray<double> array_U(nelements);
            valarray<double> array_V(nelements);
            valarray<double> array_T(nelements);
            valarray<double> array_S1(nelements);
            valarray<double> array_S2(nelements);
            valarray<double> array_S3(nelements);*/
            
            valarray<double> array_I = {};
            valarray<double> array_Q = {};
            valarray<double> array_U = {};
            valarray<double> array_V = {};
            valarray<double> array_T = {};
            
            valarray<double> array_S1 = {};
            valarray<double> array_S2 = {};
            valarray<double> array_S3 = {};

            if (matrixI != 0)
                array_I = valarray<double>(double(0), long(nelements));

            if (matrixQ != 0)
                array_Q = valarray<double>(double(0), long(nelements));

            if (matrixU != 0)
                array_U = valarray<double>(double(0), long(nelements));

            if (matrixV != 0)
                array_V = valarray<double>(double(0), long(nelements));

            if (matrixT != 0)
                array_T = valarray<double>(double(0), long(nelements));
            
            
            if (matrixS1 != 0)
                array_S1 = valarray<double>(double(0), long(nelements));

            if (matrixS2 != 0)
                array_S2 = valarray<double>(double(0), long(nelements));

            if (matrixS3 != 0)
                array_S3 = valarray<double>(double(0), long(nelements));            
            
            uint i = 0;
            for(uint i_y = 0; i_y < bins_y; i_y++)
            {
                for(uint i_x = 0; i_x < bins_x; i_x++)
                {
                    if(matrixI != 0)
                        array_I[i] = matrixI[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    
                    if(matrixQ != 0)
                        array_Q[i] = matrixQ[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    
                    if(matrixU != 0)
                        array_U[i] = matrixU[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    
                    if(matrixV != 0)
                        array_V[i] = matrixV[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    
                    if(matrixT != 0)
                        array_T[i] = matrixT[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    
                    if(matrixS1 != 0)
                        array_S1[i] = matrixS1[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    
                    if(special_param==3)
                    {
                        array_S2[i] = matrixS2[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                        array_S3[i] = matrixS3[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    }
                    
                    double N_ph = 1;
                    
                    if(N_photon != 0)
                        N_ph = N_photon[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);

                    if((results_type == RESULTS_RAY || results_type == RESULTS_FULL) &&
                        detector_id == DET_POLAR && processing_method == NEAREST && N_ph > 0)
                    {
                        if(matrixT != 0)
                            array_T[i] /= N_ph;
                        
                        if(matrixS1!=0)
                            array_S1[i] /= N_ph;
                    }
                    
                    i++;
                }
            }
            
            uint off_counter = 0;
            
            if(matrixI != 0)
            {
                fpixel[3] = off_counter * nr_extra + i_extra + 1;
                pFits->pHDU().write(fpixel, nelements, array_I);
                off_counter++;
            }
            
            if(matrixQ != 0)
            {
                fpixel[3] = off_counter * nr_extra + i_extra + 1;
                pFits->pHDU().write(fpixel, nelements, array_Q);
                off_counter++;
            }
            
            if(matrixU != 0)
            {
                fpixel[3] = off_counter * nr_extra + i_extra + 1;
                pFits->pHDU().write(fpixel, nelements, array_U);
                off_counter++;
            }
            
            if(matrixV != 0)
            {
                fpixel[3] = off_counter * nr_extra + i_extra + 1;
                pFits->pHDU().write(fpixel, nelements, array_V);
                off_counter++;
            }
            
            if(matrixT != 0)
            {
                fpixel[3] = off_counter * nr_extra + i_extra + 1;
                pFits->pHDU().write(fpixel, nelements, array_T);
                off_counter++;
            }
            
            if(matrixS1 != 0)
            {
                fpixel[3] = off_counter * nr_extra + i_extra + 1;
                pFits->pHDU().write(fpixel, nelements, array_S1);
                off_counter++;
            }
            
            if(special_param==3)
            {
                if(matrixS2 != 0)
                {
                    fpixel[3] = off_counter * nr_extra + i_extra + 1;
                    pFits->pHDU().write(fpixel, nelements, array_S2);
                    off_counter++;
                }

                if(matrixS3 != 0)
                {
                    fpixel[3] = off_counter * nr_extra + i_extra + 1;
                    pFits->pHDU().write(fpixel, nelements, array_S3);
                    off_counter++;
                }
            }
        }
    }

    double bin_width_x, first_pix_val_x, deg_per_pix_x, first_pix_val_deg_x;
    
    calcCoordinateParameters(sidelength_x,
                                bins_x,
                                map_shift_x,
                                distance,
                                bin_width_x,
                                first_pix_val_x,
                                deg_per_pix_x,
                                first_pix_val_deg_x);
    
    double bin_width_y, first_pix_val_y, deg_per_pix_y, first_pix_val_deg_y;
    
    calcCoordinateParameters(sidelength_y,
                                bins_y,
                                map_shift_y,
                                distance,
                                bin_width_y,
                                first_pix_val_y,
                                deg_per_pix_y,
                                first_pix_val_deg_y);

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
        pFits->pHDU().addKey(
            "CUNIT4", "I, Q, U, V [Jy/px], optical depth, column density [m^-2]", "unit of axis 4");
        pFits->pHDU().addKey("ETYPE", "thermal emission", "type of emission");
    }
    else if(results_type == RESULTS_FULL)
    {
        pFits->pHDU().addKey(
            "CUNIT4", "I, Q, U, V [Jy/px], optical depth, column density [m^-2]", "unit of axis 4");
        pFits->pHDU().addKey("ETYPE", "thermal emission (+scattering)", "type of emission");
    }
    else if(results_type == RESULTS_MC)
    {
        pFits->pHDU().addKey("CUNIT4", "I, Q, U, V, I_direct, I_scat [Jy/px]", "unit of axis 4");
        pFits->pHDU().addKey("ETYPE", "scattered emission / direct stellar emission", "type of emission");
    }
    else if(nr_extra > 1)
    {
        pFits->pHDU().addKey(
            "CUNIT4", "I, Q, U, V [Jy/px], optical depth, column density [m^-2]", "unit of axis 4");
        pFits->pHDU().addKey(
            "ETYPE", "Each Stokes entry for FULL / SCAT / THERMAL / STOCHASTIC", "type of emission");
    }
    pFits->pHDU().addKey("ID", nr, "detector ID");

    for(uint i_extra = 0; i_extra < nr_extra; i_extra++)
    {
        for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
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
    string alignment_descr = getAlignmentDescription();
    if(alignment_descr != "")
        pFits->pHDU().addKey("ALIGNMENT", alignment_descr, "alignment method of dust grains");

    cout << CLR_LINE << flush;
    return true;
}

bool CDetector::writeMapStats(uint nr, uint results_type)
{
    cout << CLR_LINE << flush;
    cout << " -> Writing map with statistics ...  \r" << flush;
    // auto_ptr<CCfits::FITS> pFits(0);
    unique_ptr<CCfits::FITS> pFits;

    try
    {
        char str_tmp[1024];
        char str_end[1024];
        nr++;

#ifdef WINDOWS
        strcpy_s(str_tmp, "polaris_detector_nr%04d_stats");
        sprintf_s(str_end, str_tmp, nr);
#else
        strcpy(str_tmp, "polaris_detector_nr%04d_stats");
        sprintf(str_end, str_tmp, nr);
#endif

        string path_out = path + str_end + FITS_COMPRESS_EXT;
        remove(path_out.c_str());

        long naxis = 4;
        long naxes[4] = { bins_x, bins_y, nr_spectral_bins, 11 * nr_extra };

        pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
    }
    catch(CCfits::FITS::CantCreate)
    {
        return false;
    }

    long nelements = bins_x * bins_y;
    if(max_cells != nelements)
    {
        cout << WARNING_LINE << "Max cells are not equal to bins x bins!" << endl;
        return false;
    }

    vector<long> fpixel(4);
    fpixel[0] = 1;
    fpixel[1] = 1;

    for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
    {
        fpixel[2] = i_spectral + 1;
        for(uint i_extra = 0; i_extra < nr_extra; i_extra++)
        {
            valarray<double> array_N_photon(nelements);
            valarray<double> array_R_I(nelements);
            valarray<double> array_VOV_I(nelements);
            valarray<double> array_R_Q(nelements);
            valarray<double> array_VOV_Q(nelements);
            valarray<double> array_R_U(nelements);
            valarray<double> array_VOV_U(nelements);
            valarray<double> array_R_V(nelements);
            valarray<double> array_VOV_V(nelements);
            valarray<double> array_R_PI(nelements);
            valarray<double> array_VOV_PI(nelements);

            uint i = 0;
            for(uint i_y = 0; i_y < bins_y; i_y++)
                for(uint i_x = 0; i_x < bins_x; i_x++)
                {
                    array_N_photon[i] = N_photon[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    array_R_I[i] = calc_R(i_spectral + i_extra * nr_spectral_bins, i_x, i_y, 0);
                    array_VOV_I[i] = calc_VOV(i_spectral + i_extra * nr_spectral_bins, i_x, i_y, 0);
                    array_R_Q[i] = calc_R(i_spectral + i_extra * nr_spectral_bins, i_x, i_y, 1);
                    array_VOV_Q[i] = calc_VOV(i_spectral + i_extra * nr_spectral_bins, i_x, i_y, 1);
                    array_R_U[i] = calc_R(i_spectral + i_extra * nr_spectral_bins, i_x, i_y, 2);
                    array_VOV_U[i] = calc_VOV(i_spectral + i_extra * nr_spectral_bins, i_x, i_y, 2);
                    array_R_V[i] = calc_R(i_spectral + i_extra * nr_spectral_bins, i_x, i_y, 3);
                    array_VOV_V[i] = calc_VOV(i_spectral + i_extra * nr_spectral_bins, i_x, i_y, 3);
                    array_R_PI[i] = calc_R(i_spectral + i_extra * nr_spectral_bins, i_x, i_y, 4);
                    array_VOV_PI[i] = calc_VOV(i_spectral + i_extra * nr_spectral_bins, i_x, i_y, 4);
                    i++;
                }

            fpixel[3] = 0 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_N_photon);
            fpixel[3] = 1 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_R_I);
            fpixel[3] = 2 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_VOV_I);
            fpixel[3] = 3 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_R_Q);
            fpixel[3] = 4 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_VOV_Q);
            fpixel[3] = 5 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_R_U);
            fpixel[3] = 6 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_VOV_U);
            fpixel[3] = 7 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_R_V);
            fpixel[3] = 8 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_VOV_V);
            fpixel[3] = 9 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_R_PI);
            fpixel[3] = 10 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_VOV_PI);
        }
    }

    double bin_width_x, first_pix_val_x, deg_per_pix_x, first_pix_val_deg_x;
    calcCoordinateParameters(sidelength_x,
                                bins_x,
                                map_shift_x,
                                distance,
                                bin_width_x,
                                first_pix_val_x,
                                deg_per_pix_x,
                                first_pix_val_deg_x);
    double bin_width_y, first_pix_val_y, deg_per_pix_y, first_pix_val_deg_y;
    calcCoordinateParameters(sidelength_y,
                                bins_y,
                                map_shift_y,
                                distance,
                                bin_width_y,
                                first_pix_val_y,
                                deg_per_pix_y,
                                first_pix_val_deg_y);

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
    pFits->pHDU().addKey("CUNIT4", "dimensionless", "unit of axis 4");
    pFits->pHDU().addKey("ETYPE", "N, I(R,VOV), Q(R,VOV), U(R,VOV), V(R,VOV), PI(R,VOV)", "quantities");
    pFits->pHDU().writeComment("N: Number of photon packages,");
    pFits->pHDU().writeComment("R: relative error, VOV: variance of the variance");
    pFits->pHDU().writeComment("(see Camps & Baes 2018, ApJ 861:80, Eqs. 14 and 15)");

    pFits->pHDU().addKey("ID", nr, "detector ID");

    for(uint i_extra = 0; i_extra < nr_extra; i_extra++)
    {
        for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
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
    string alignment_descr = getAlignmentDescription();
    if(alignment_descr != "")
        pFits->pHDU().addKey("ALIGNMENT", alignment_descr, "alignment method of dust grains");

    cout << CLR_LINE << flush;
    return true;
}

bool CDetector::writeSed(uint nr, uint results_type)
{
    cout << CLR_LINE << flush;
    cout << " -> Writing SED ...  \r" << flush;
    // auto_ptr<CCfits::FITS> pFits(0);
    unique_ptr<CCfits::FITS> pFits;

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

        string path_out = path + str_end + FITS_COMPRESS_EXT;
        remove(path_out.c_str());

        long naxis = 3;
        long naxes[3] = { nr_spectral_bins, nr_extra, 5 };
        pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
    }
    catch(CCfits::FITS::CantCreate)
    {
        return false;
    }

    long nelements = uint(nr_spectral_bins);

    vector<long> fpixel(3);
    fpixel[0] = 1;

    std::valarray<double> array_I(nelements);
    std::valarray<double> array_Q(nelements);
    std::valarray<double> array_U(nelements);
    std::valarray<double> array_V(nelements);
    std::valarray<double> array_T(nelements);
    std::valarray<double> array_S(nelements);
    for(uint i_extra = 0; i_extra < nr_extra; i_extra++)
    {
        fpixel[1] = i_extra + 1;

        for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
        {
            array_I[i_spectral] = sedI[i_spectral + i_extra * nr_spectral_bins];
            array_Q[i_spectral] = sedQ[i_spectral + i_extra * nr_spectral_bins];
            array_U[i_spectral] = sedU[i_spectral + i_extra * nr_spectral_bins];
            array_V[i_spectral] = sedV[i_spectral + i_extra * nr_spectral_bins];
            array_T[i_spectral] = sedT[i_spectral + i_extra * nr_spectral_bins];
            if(results_type == RESULTS_MC)
                array_S[i_spectral] = sedS[i_spectral + i_extra * nr_spectral_bins];

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

            if(results_type == RESULTS_MC)
            {
                fpixel[2] = 6;
                pFits->pHDU().write(fpixel, nelements, array_S);
            }
        }
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
        pFits->pHDU().addKey("CUNIT3", "I, Q, U, V, I_direct, I_scat [Jy]", "unit of axis 3");
        pFits->pHDU().addKey("ETYPE", "scattered emission / direct stellar emission", "type of emission");
    }
    pFits->pHDU().addKey("ID", nr, "detector ID");
    for(uint i_extra = 0; i_extra < nr_extra; i_extra++)
    {
        for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
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
    string alignment_descr = getAlignmentDescription();
    if(alignment_descr != "")
        pFits->pHDU().addKey("ALIGNMENT", alignment_descr, "alignment method of dust grains");

    cout << CLR_LINE << flush;

    return true;
}

bool CDetector::writeDustHealMaps(uint nr, ilist64 & heal_index, long npix, uint results_type)
{
    cout << CLR_LINE;
    cout << " -> Writing healpix map ...  \r" << flush;
    // auto_ptr<CCfits::FITS> pFits(0);
    unique_ptr<CCfits::FITS> pFits;

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

        string path_out = path + str_end + FITS_COMPRESS_EXT;
        remove(path_out.c_str());

        long naxis = 1;
        long naxes[1] = { 0 };
        pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
    }
    catch(CCfits::FITS::CantCreate)
    {
        return false;
    }

    long nelements = bins_x * bins_y;
    if(max_cells != nelements)
    {
        cout << ERROR_LINE << "Max. cells are not equal to bins x bins!" << endl;
        return false;
    }
    
    if(npix < nelements)
    {
        cout << ERROR_LINE << "Nr. of pixel is smaller than number of ray-traced elements!" << endl;
        return false;
    }
    
    for (uint i_x = 0; i_x < bins_x; ++i_x)
    {
        int64_t idx = heal_index[i_x];
        if (idx < 0 || idx >= npix)
        {
            cout << ERROR_LINE << "heal_index[" << i_x << "] = "
                      << idx << " out of range 0.." << (npix - 1) << endl << flush;
            
            return false;
        }
    }

    // Init columns
    string newName("HEALPIX_EXTENSION");
    uint nr_of_quantities = 0;
    uint off=0;
    
    if(matrixI!=0)
        nr_of_quantities++;
    
    if(matrixQ!=0)
        nr_of_quantities++;
    
    if(matrixU!=0)
        nr_of_quantities++;
    
    if(matrixV!=0)
        nr_of_quantities++;
        
    if(matrixT!=0)
        nr_of_quantities++;
    
    if(matrixS1!=0)
        off=1;
    
    if(nr_of_quantities==0)
    {
        cout << ERROR_LINE << "No quantities selected to write!" << endl;
        return false;
    }
    
    vector<string> colName(nr_of_quantities * nr_spectral_bins + off, "");
    vector<string> colForm(nr_of_quantities * nr_spectral_bins + off, "");
    vector<string> colUnit(nr_of_quantities * nr_spectral_bins + off, "");

    for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
    {
        char str_1[1024];
        uint off_counter=0;
        
        if(matrixI!=0)
        {
            sprintf(str_1, "I_STOKES (WAVELENGTH = %e [m])", wavelength_list_det[i_spectral]);
            
            colName[nr_of_quantities * i_spectral + off_counter] = str_1;
            colForm[nr_of_quantities * i_spectral + off_counter] = "D";
            colUnit[nr_of_quantities * i_spectral + off_counter] = "Jy/px";
            
            off_counter++;
        }
        
        if(matrixQ!=0)
        {
            sprintf(str_1, "Q_STOKES (WAVELENGTH = %e [m])", wavelength_list_det[i_spectral]);
            
            colName[nr_of_quantities * i_spectral + off_counter] = str_1;
            colForm[nr_of_quantities * i_spectral + off_counter] = "D";
            colUnit[nr_of_quantities * i_spectral + off_counter] = "Jy/px";
            
            off_counter++;
        }
        
        if(matrixU!=0)
        {
            sprintf(str_1, "U_STOKES (WAVELENGTH = %e [m])", wavelength_list_det[i_spectral]);
            
            colName[nr_of_quantities * i_spectral + off_counter] = str_1;
            colForm[nr_of_quantities * i_spectral + off_counter] = "D";
            colUnit[nr_of_quantities * i_spectral + off_counter] = "Jy/px";
            
            off_counter++;
        }
        
        if(matrixV!=0)
        {
            sprintf(str_1, "V_STOKES (WAVELENGTH = %e [m])", wavelength_list_det[i_spectral]);
            
            colName[nr_of_quantities * i_spectral + off_counter] = str_1;
            colForm[nr_of_quantities * i_spectral + off_counter] = "D";
            colUnit[nr_of_quantities * i_spectral + off_counter] = "Jy/px";
            
            off_counter++;
        }
        
        if(matrixT!=0)
        {
            sprintf(str_1, "OPTICAL_DEPTH (WAVELENGTH = %e [m])", wavelength_list_det[i_spectral]);
            
            colName[nr_of_quantities * i_spectral + off_counter] = str_1;
            colForm[nr_of_quantities * i_spectral + off_counter] = "D";
            colUnit[nr_of_quantities * i_spectral + off_counter] = "";
            
            off_counter++;
        }
    }
    
    if(off==1)
    {
        colName[nr_of_quantities * nr_spectral_bins] = "COLUMN_DENSITY";
        colForm[nr_of_quantities * nr_spectral_bins] = "D";
        colUnit[nr_of_quantities * nr_spectral_bins] = "m^-2";
    }

    CCfits::Table * newTable = pFits->addTable(newName, npix, colName, colForm, colUnit);
    
    for (uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
    {
        valarray<double> array_I = {};
        valarray<double> array_Q = {};
        valarray<double> array_U = {};
        valarray<double> array_V = {};
        valarray<double> array_T = {};

        if (matrixI != 0)
            array_I = valarray<double>(double(0), long(npix));

        if (matrixQ != 0)
            array_Q = valarray<double>(double(0), long(npix));

        if (matrixU != 0)
            array_U = valarray<double>(double(0), long(npix));

        if (matrixV != 0)
            array_V = valarray<double>(double(0), long(npix));

        if (matrixT != 0)
            array_T = valarray<double>(double(0), long(npix));

        uint i_y = 0;
        uint off_counter = 0;

        for (uint i_x = 0; i_x < bins_x; i_x++)
        {
            long index = (long)heal_index[i_x];

            if (index >= npix)
                continue;

            if (matrixI != 0)
                array_I[index] = matrixI[i_spectral](i_x, i_y);

            if (matrixQ != 0)
                array_Q[index] = matrixQ[i_spectral](i_x, i_y);

            if (matrixU != 0)
                array_U[index] = matrixU[i_spectral](i_x, i_y);

            if (matrixV != 0)
                array_V[index] = matrixV[i_spectral](i_x, i_y);

            if (matrixT != 0)
                array_T[index] = matrixT[i_spectral](i_x, i_y);
        }


        if (matrixI != 0)
        {
            newTable->column(colName[nr_of_quantities * i_spectral + off_counter]).write(array_I, 1);
            off_counter++;
        }

        if (matrixQ != 0)
        {
            newTable->column(colName[nr_of_quantities * i_spectral + off_counter]).write(array_Q, 1);
            off_counter++;
        }

        if (matrixU != 0)
        {
            newTable->column(colName[nr_of_quantities * i_spectral + off_counter]).write(array_U, 1);
            off_counter++;
        }

        if (matrixV != 0)
        {
            newTable->column(colName[nr_of_quantities * i_spectral + off_counter]).write(array_V, 1);
            off_counter++;
        }

        if (matrixT != 0)
        {
            newTable->column(colName[nr_of_quantities * i_spectral + off_counter]).write(array_T, 1);
            off_counter++;
        }
    }

    // Column density
    if(off==1)
    {
        valarray<double> array_C(double(0), long(npix));

        uint i_y = 0;
        for (uint i_x = 0; i_x < bins_x; i_x++)
        {
            long index = (long)heal_index[i_x];
            if (index >= npix)
                continue;

            array_C[index] = matrixS1[0](i_x, i_y);
        }

        newTable->column(colName[nr_of_quantities * nr_spectral_bins]).write(array_C, 1);

    }
    
    // Healpix parameter
    newTable->addKey("PIXTYPE", "HEALPIX", "Pixel algorigthm");
    newTable->addKey("ORDERING", "RING", "Ordering scheme");
    newTable->addKey("INDXSCHM", "IMPLICIT", "Indexing scheme");
    newTable->addKey("NSIDE", uint(sqrt(npix / 12.0)), "Resolution Parameter");
    newTable->addKey("FIRSTPIX", 0, "First pixel (0 based)");
    newTable->addKey("LASTPIX", npix - 1, "Last pixel (0 based)");
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
    newTable->addKey("BUBBLE_RADIUS", rad_bubble, "bubble radius [m]");

    string alignment_descr = getAlignmentDescription();
    if(alignment_descr != "")
        pFits->pHDU().addKey("ALIGNMENT", alignment_descr, "alignment method of dust grains");

    cout << CLR_LINE << flush;
    return true;
}

bool CDetector::writeDustHealMapsTiny(uint nr, ilist64 & heal_index, long npix, uint results_type)
{
    cout << CLR_LINE << flush;
    cout << " -> Writing healpix min. map ...  \r" << flush;
    // auto_ptr<CCfits::FITS> pFits(0);
    unique_ptr<CCfits::FITS> pFits;

    try
    {
        char str_tmp[1024];
        char str_end[1024];
        nr++;

#ifdef WINDOWS
        strcpy_s(str_tmp, "min_polaris_detector_nr%04d");
        sprintf_s(str_end, str_tmp, nr);
#else
        strcpy(str_tmp, "min_polaris_detector_nr%04d");
        sprintf(str_end, str_tmp, nr);
#endif

        string path_out = path + str_end + FITS_COMPRESS_EXT;
        remove(path_out.c_str());

        long naxis = 1;
        long naxes[1] = { 0 };
        pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
    }
    catch(CCfits::FITS::CantCreate)
    {
        return false;
    }

    long nelements = bins_x * bins_y;
    if(max_cells != nelements)
    {
        cout << WARNING_LINE << "Max cells are not equal to bins x bins!" << endl;
        return false;
    }
    
    for (uint i_x = 0; i_x < bins_x; ++i_x)
    {
        int64_t index = heal_index[i_x];
        if (index < 0 || index >= npix)
        {
            cout << ERROR_LINE << "heal_index[" << i_x << "] = "
                      << index << " out of range 0.." << (npix - 1) << endl << flush;
            
            return false;
        }
    }

    // Init columns
    string newName("DATA_EXTENSION");
    uint nr_of_quantities = 0;
    uint off=1;
    
    if(matrixI!=0)
        nr_of_quantities++;
    
    if(matrixQ!=0)
        nr_of_quantities++;
    
    if(matrixU!=0)
        nr_of_quantities++;
    
    if(matrixV!=0)
        nr_of_quantities++;
        
    if(matrixT!=0)
        nr_of_quantities++;
    
    if(matrixS1!=0)
        off=2;
    
    if(nr_of_quantities==0)
    {
        cout << ERROR_LINE << "No quantities selected to write!" << endl;
        return false;
    }
    
    vector<string> colName(nr_of_quantities * nr_spectral_bins + off, "");
    vector<string> colForm(nr_of_quantities * nr_spectral_bins + off, "");
    vector<string> colUnit(nr_of_quantities * nr_spectral_bins + off, "");

    for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
    {
        char str_1[1024];
        uint off_counter=0;
        
        if(matrixI!=0)
        {
            sprintf(str_1, "I_STOKES (WAVELENGTH = %e [m])", wavelength_list_det[i_spectral]);
            
            colName[nr_of_quantities * i_spectral + off_counter] = str_1;
            colForm[nr_of_quantities * i_spectral + off_counter] = "D";
            colUnit[nr_of_quantities * i_spectral + off_counter] = "Jy/px";
            
            off_counter++;
        }
        
        if(matrixQ!=0)
        {
            sprintf(str_1, "Q_STOKES (WAVELENGTH = %e [m])", wavelength_list_det[i_spectral]);
            
            colName[nr_of_quantities * i_spectral + off_counter] = str_1;
            colForm[nr_of_quantities * i_spectral + off_counter] = "D";
            colUnit[nr_of_quantities * i_spectral + off_counter] = "Jy/px";
            
            off_counter++;
        }
        
        if(matrixU!=0)
        {
            sprintf(str_1, "U_STOKES (WAVELENGTH = %e [m])", wavelength_list_det[i_spectral]);
            
            colName[nr_of_quantities * i_spectral + off_counter] = str_1;
            colForm[nr_of_quantities * i_spectral + off_counter] = "D";
            colUnit[nr_of_quantities * i_spectral + off_counter] = "Jy/px";
            
            off_counter++;
        }
        
        if(matrixV!=0)
        {
            sprintf(str_1, "V_STOKES (WAVELENGTH = %e [m])", wavelength_list_det[i_spectral]);
            
            colName[nr_of_quantities * i_spectral + off_counter] = str_1;
            colForm[nr_of_quantities * i_spectral + off_counter] = "D";
            colUnit[nr_of_quantities * i_spectral + off_counter] = "Jy/px";
            
            off_counter++;
        }
        
        if(matrixT!=0)
        {
            sprintf(str_1, "OPTICAL_DEPTH (WAVELENGTH = %e [m])", wavelength_list_det[i_spectral]);
            
            colName[nr_of_quantities * i_spectral + off_counter] = str_1;
            colForm[nr_of_quantities * i_spectral + off_counter] = "D";
            colUnit[nr_of_quantities * i_spectral + off_counter] = "";
            
            off_counter++;
        }
    }
    
    colName[nr_of_quantities * nr_spectral_bins] = "HEAL_INDEX";
    colForm[nr_of_quantities * nr_spectral_bins] = "J";
    colUnit[nr_of_quantities * nr_spectral_bins] = "";
    
    if(off==2)
    {
        colName[nr_of_quantities * nr_spectral_bins+1] = "COLUMN_DENSITY";
        colForm[nr_of_quantities * nr_spectral_bins+1] = "D";
        colUnit[nr_of_quantities * nr_spectral_bins+1] = "m^-2";
    }
        
    CCfits::Table * newTable = pFits->addTable(newName, nelements, colName, colForm, colUnit);
    
    for (uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
    {
        valarray<double> array_I = {};
        valarray<double> array_Q = {};
        valarray<double> array_U = {};
        valarray<double> array_V = {};
        valarray<double> array_T = {};

        if (matrixI != 0)
            array_I = valarray<double>(double(0), long(nelements));

        if (matrixQ != 0)
            array_Q = valarray<double>(double(0), long(nelements));

        if (matrixU != 0)
            array_U = valarray<double>(double(0), long(nelements));

        if (matrixV != 0)
            array_V = valarray<double>(double(0), long(nelements));

        if (matrixT != 0)
            array_T = valarray<double>(double(0), long(nelements));

        uint i_y = 0;
        uint off_counter = 0;

        for (uint i_x = 0; i_x < bins_x; i_x++)
        {
            if (matrixI != 0)
                array_I[i_x] = matrixI[i_spectral](i_x, i_y);

            if (matrixQ != 0)
                array_Q[i_x] = matrixQ[i_spectral](i_x, i_y);

            if (matrixU != 0)
                array_U[i_x] = matrixU[i_spectral](i_x, i_y);

            if (matrixV != 0)
                array_V[i_x] = matrixV[i_spectral](i_x, i_y);

            if (matrixT != 0)
                array_T[i_x] = matrixT[i_spectral](i_x, i_y);
        }

        if (matrixI != 0)
        {
            newTable->column(colName[nr_of_quantities * i_spectral + off_counter]).write(array_I, 1);
            off_counter++;
        }

        if (matrixQ != 0)
        {
            newTable->column(colName[nr_of_quantities * i_spectral + off_counter]).write(array_Q, 1);
            off_counter++;
        }

        if (matrixU != 0)
        {
            newTable->column(colName[nr_of_quantities * i_spectral + off_counter]).write(array_U, 1);
            off_counter++;
        }

        if (matrixV != 0)
        {
            newTable->column(colName[nr_of_quantities * i_spectral + off_counter]).write(array_V, 1);
            off_counter++;
        }

        if (matrixT != 0)
        {
            newTable->column(colName[nr_of_quantities * i_spectral + off_counter]).write(array_T, 1);
            off_counter++;
        }
    }

    // Column density
    valarray<long> array_index(nelements);
    valarray<double> array_C = {};

    if (matrixS1 != 0)
        array_C = valarray<double>(double(0), long(npix));
    
    long i_y=0;
    for(long i_x = 0; i_x < bins_x; i_x++)
    {
        array_index[i_x]=(long)heal_index[i_x];
        
        if(matrixS1!=0)
            array_C[i_x] = matrixS1[0](i_x, i_y);
    }
        
    newTable->column(colName[nr_of_quantities * nr_spectral_bins + 0]).write(array_index, 1);
    
    if(matrixS1!=0)
        newTable->column(colName[nr_of_quantities * nr_spectral_bins + 1]).write(array_C, 1);

    // Healpix parameter
    newTable->addKey("PIXTYPE", "HEALPIX", "Pixel algorigthm");
    newTable->addKey("ORDERING", "RING", "Ordering scheme");
    newTable->addKey("INDXSCHM", "EXPLICIT", "Indexing scheme");
    newTable->addKey("NSIDE", uint(sqrt(npix / 12.0)), "Resolution Parameter");
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
    newTable->addKey("BUBBLE_RADIUS", rad_bubble, "bubble radius [m]");
    
    string alignment_descr = getAlignmentDescription();
    if(alignment_descr != "")
        pFits->pHDU().addKey("ALIGNMENT", alignment_descr, "alignment method of dust grains");

    cout << CLR_LINE << flush;
    return true;
}

// ame maps

bool CDetector::writeDustAMEMap(uint nr)
{
    cout << CLR_LINE << flush;
    cout << " -> Writing ame map ...  \r" << flush;
    // auto_ptr<CCfits::FITS> pFits(0);
    unique_ptr<CCfits::FITS> pFits;

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

        string path_out = path + str_end + FITS_COMPRESS_EXT;
        remove(path_out.c_str());

        long naxis = 4;
        long naxes[4] = { bins_x, bins_y, nr_spectral_bins, 3 * nr_extra };

        pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
    }
    catch(CCfits::FITS::CantCreate)
    {
        return false;
    }

    long nelements = bins_x * bins_y;
    if(max_cells != nelements)
    {
        cout << WARNING_LINE << "Max cells are not equal to bins x bins!" << endl;
        return false;
    }

    vector<long> fpixel(4);
    fpixel[0] = 1;
    fpixel[1] = 1;

    for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
    {
        fpixel[2] = i_spectral + 1;
        for(uint i_extra = 0; i_extra < nr_extra; i_extra++)
        {
            valarray<double> array_I(nelements);
            /*valarray<double> array_Q(nelements);
            valarray<double> array_U(nelements);
            valarray<double> array_V(nelements);
            valarray<double> array_T(nelements);*/
            valarray<double> array_S1(nelements);
            valarray<double> array_S2(nelements);
            //valarray<double> array_S3(nelements);

            uint i = 0;
            for(uint i_y = 0; i_y < bins_y; i_y++)
            {
                for(uint i_x = 0; i_x < bins_x; i_x++)
                {
                    array_I[i] = matrixI[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    /*array_Q[i] = matrixQ[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    array_U[i] = matrixU[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    array_V[i] = matrixV[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    array_T[i] = matrixT[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);*/
                    
                    array_S1[i] = matrixS1[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    array_S2[i] = matrixS2[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    
                    double N_ph = N_photon[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);

                    if(N_ph > 0)
                    {
                        //array_T[i] /= N_ph;
                        array_S1[i] /= N_ph;
                        array_S2[i] /= N_ph;
                    }
                    i++;
                }
            }    
               
            fpixel[3] = 0 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_I);
            
            /*fpixel[3] = 1 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_Q);
            
            fpixel[3] = 2 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_U);
            
            fpixel[3] = 3 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_V);
            
            fpixel[3] = 4 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_T);*/
            
            fpixel[3] = 1 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_S1);
            
            fpixel[3] = 2 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_S2);
        }
    }

    double bin_width_x, first_pix_val_x, deg_per_pix_x, first_pix_val_deg_x;
    
    calcCoordinateParameters(sidelength_x,
                                bins_x,
                                map_shift_x,
                                distance,
                                bin_width_x,
                                first_pix_val_x,
                                deg_per_pix_x,
                                first_pix_val_deg_x);
    
    double bin_width_y, first_pix_val_y, deg_per_pix_y, first_pix_val_deg_y;
    
    calcCoordinateParameters(sidelength_y,
                                bins_y,
                                map_shift_y,
                                distance,
                                bin_width_y,
                                first_pix_val_y,
                                deg_per_pix_y,
                                first_pix_val_deg_y);

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
    
    pFits->pHDU().addKey(
        "CUNIT4", "I [Jy/px], ND, NH [m^-2]", "unit of axis 4");
    pFits->pHDU().addKey("ETYPE", "thermal emission", "type of emission");
    
    pFits->pHDU().addKey("ID", nr, "detector ID");

    for(uint i_extra = 0; i_extra < nr_extra; i_extra++)
    {
        for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
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
    
    /*string alignment_descr = getAlignmentDescription();
    if(alignment_descr != "")
        pFits->pHDU().addKey("ALIGNMENT", alignment_descr, "alignment method of dust grains");*/

    cout << CLR_LINE << flush;
    return true;
}

bool CDetector::writeDustAMEHealMap(uint nr)
{
    cout << CLR_LINE << flush;
    cout << " -> Writing ame map ...  \r" << flush;
    // auto_ptr<CCfits::FITS> pFits(0);
    unique_ptr<CCfits::FITS> pFits;

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

        string path_out = path + str_end + FITS_COMPRESS_EXT;
        remove(path_out.c_str());

        long naxis = 4;
        long naxes[4] = { bins_x, bins_y, nr_spectral_bins, 6 * nr_extra };

        pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
    }
    catch(CCfits::FITS::CantCreate)
    {
        return false;
    }

    long nelements = bins_x * bins_y;
    if(max_cells != nelements)
    {
        cout << WARNING_LINE << "Max cells are not equal to bins x bins!" << endl;
        return false;
    }

    vector<long> fpixel(4);
    fpixel[0] = 1;
    fpixel[1] = 1;

    for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
    {
        fpixel[2] = i_spectral + 1;
        for(uint i_extra = 0; i_extra < nr_extra; i_extra++)
        {
            valarray<double> array_I(nelements);
            /*valarray<double> array_Q(nelements);
            valarray<double> array_U(nelements);
            valarray<double> array_V(nelements);
            valarray<double> array_T(nelements);*/
            valarray<double> array_S1(nelements);
            valarray<double> array_S2(nelements);
            //valarray<double> array_S3(nelements);

            uint i = 0;
            for(uint i_y = 0; i_y < bins_y; i_y++)
                for(uint i_x = 0; i_x < bins_x; i_x++)
                {
                    array_I[i] = matrixI[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                   /* array_Q[i] = matrixQ[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    array_U[i] = matrixU[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    array_V[i] = matrixV[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    array_T[i] = matrixT[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);*/
                    
                    array_S1[i] = matrixS1[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    array_S2[i] = matrixS2[i_spectral + i_extra * nr_spectral_bins](i_x, i_y);
                    
                    i++;
                }

            fpixel[3] = 0 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_I);
            
            /*fpixel[3] = 1 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_Q);
            
            fpixel[3] = 2 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_U);
            
            fpixel[3] = 3 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_V);
            
            fpixel[3] = 4 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_T);*/
            
            fpixel[3] = 1 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_S1);
            
            fpixel[3] = 2 * nr_extra + i_extra + 1;
            pFits->pHDU().write(fpixel, nelements, array_S2);
            
        }
    }

    double bin_width_x, first_pix_val_x, deg_per_pix_x, first_pix_val_deg_x;
    
    calcCoordinateParameters(sidelength_x,
                                bins_x,
                                map_shift_x,
                                distance,
                                bin_width_x,
                                first_pix_val_x,
                                deg_per_pix_x,
                                first_pix_val_deg_x);
    
    double bin_width_y, first_pix_val_y, deg_per_pix_y, first_pix_val_deg_y;
    
    calcCoordinateParameters(sidelength_y,
                                bins_y,
                                map_shift_y,
                                distance,
                                bin_width_y,
                                first_pix_val_y,
                                deg_per_pix_y,
                                first_pix_val_deg_y);

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
    
    pFits->pHDU().addKey(
            "CUNIT4", "I [Jy/px], ND, NH [m^-2]", "unit of axis 4");
        pFits->pHDU().addKey("ETYPE", "thermal emission", "type of emission");
    
    pFits->pHDU().addKey("ID", nr, "detector ID");

    for(uint i_extra = 0; i_extra < nr_extra; i_extra++)
    {
        for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
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
    string alignment_descr = getAlignmentDescription();
    
    //if(alignment_descr != "")
    //    pFits->pHDU().addKey("ALIGNMENT", alignment_descr, "alignment method of dust grains");

    cout << CLR_LINE << flush;
    return true;
}

// free-free maps

bool CDetector::writeFreeFreeMap(uint nr)
{
    cout << CLR_LINE << flush;
    cout << " -> Writing free-free plane detector ...  \r" << flush;
    // auto_ptr<CCfits::FITS> pFits(0);
    unique_ptr<CCfits::FITS> pFits;

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

        string path_out = path + str_end + FITS_COMPRESS_EXT;
        remove(path_out.c_str());

        long naxis = 4;
        long naxes[4] = { bins_x, bins_y, nr_spectral_bins, 5 * nr_extra };
        pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
    }
    catch(CCfits::FITS::CantCreate)
    {
        return false;
    }

    long nelements = bins_x * bins_y;
    if(max_cells != nelements)
    {
        cout << WARNING_LINE << "Max cells are not equal to bins x bins!" << endl;
        return false;
    }

    vector<long> fpixel(4);
    fpixel[0] = 1;
    fpixel[1] = 1;

    for(uint i_electrons = 0; i_electrons < 1; i_electrons++)
    {
        for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
        {
            fpixel[2] = i_spectral + 1;

            valarray<double> array_I(nelements);
            //valarray<double> array_Q(nelements);
            //valarray<double> array_U(nelements);
            //valarray<double> array_V(nelements);
            valarray<double> array_T(nelements);
            valarray<double> array_NH(nelements);
            valarray<double> array_NE(nelements);
            valarray<double> array_NI(nelements);

            uint i = 0;
            for(uint i_y = 0; i_y < bins_y; i_y++)
                for(uint i_x = 0; i_x < bins_x; i_x++)
                {
                    array_I[i] = matrixI[i_spectral + i_electrons * nr_spectral_bins](i_x, i_y);
                    
                    /*array_Q[i] = matrixQ[i_spectral + i_electrons * nr_spectral_bins](i_x, i_y);
                    array_U[i] = matrixU[i_spectral + i_electrons * nr_spectral_bins](i_x, i_y);
                    array_V[i] = matrixV[i_spectral + i_electrons * nr_spectral_bins](i_x, i_y);*/
                    array_T[i] = matrixT[i_spectral + i_electrons * nr_spectral_bins](i_x, i_y);
                    
                    array_NH[i] = matrixS1[i_spectral + i_electrons * nr_spectral_bins](i_x, i_y);
                    array_NE[i] = matrixS2[i_spectral + i_electrons * nr_spectral_bins](i_x, i_y);
                    array_NI[i] = matrixS3[i_spectral + i_electrons * nr_spectral_bins](i_x, i_y);
                    i++;
                }

            fpixel[3] = 6 * i_electrons + 1;
            pFits->pHDU().write(fpixel, nelements, array_I);
            
            /*fpixel[3] = 6 * i_electrons + 2;
            pFits->pHDU().write(fpixel, nelements, array_Q);
            fpixel[3] = 6 * i_electrons + 3;
            pFits->pHDU().write(fpixel, nelements, array_U);
            fpixel[3] = 6 * i_electrons + 4;
            pFits->pHDU().write(fpixel, nelements, array_V);*/
            
            fpixel[3] = 6 * i_electrons + 2;
            pFits->pHDU().write(fpixel, nelements, array_T);
            fpixel[3] = 6 * i_electrons + 3;
            pFits->pHDU().write(fpixel, nelements, array_NH);
            fpixel[3] = 6 * i_electrons + 4;
            pFits->pHDU().write(fpixel, nelements, array_NE);
            fpixel[3] = 6 * i_electrons + 5;
            pFits->pHDU().write(fpixel, nelements, array_NI);
        }
    }

    double bin_width_x, first_pix_val_x, deg_per_pix_x, first_pix_val_deg_x;
    calcCoordinateParameters(sidelength_x,
                                bins_x,
                                map_shift_x,
                                distance,
                                bin_width_x,
                                first_pix_val_x,
                                deg_per_pix_x,
                                first_pix_val_deg_x);
    double bin_width_y, first_pix_val_y, deg_per_pix_y, first_pix_val_deg_y;
    calcCoordinateParameters(sidelength_y,
                                bins_y,
                                map_shift_y,
                                distance,
                                bin_width_y,
                                first_pix_val_y,
                                deg_per_pix_y,
                                first_pix_val_deg_y);

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

    pFits->pHDU().addKey("CUNIT4",
                            "I, Q, U, V [Jy/px], tau, NH, NE, NI [m^-2]",
                            "unit of axis 4");
    
    pFits->pHDU().addKey("CUNIT4",
                            "I [Jy/px], tau [1], NH, NE, NI [m^-2]",
                            "unit of axis 4");
        
    pFits->pHDU().addKey("ETYPE", "free-free emission", "type of emission");

    pFits->pHDU().addKey("ID", nr, "detector ID");
    pFits->pHDU().writeComment("Jy/px is Jy per surface area of that pixel");

    for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
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


bool CDetector::writeFreeFreeHealMap(uint nr)
{
    cout << CLR_LINE << flush;
    cout << " -> Writing free-free healpix detector ...  \r" << flush;
    // auto_ptr<CCfits::FITS> pFits(0);
    unique_ptr<CCfits::FITS> pFits;

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

        string path_out = path + str_end + FITS_COMPRESS_EXT;
        remove(path_out.c_str());

        long naxis = 1;
        long naxes[1] = { 0 };
        pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
    }
    catch(CCfits::FITS::CantCreate)
    {
        return false;
    }

    long nelements = bins_x * bins_y;
    if(max_cells != nelements)
    {
        cout << WARNING_LINE << "Max cells are not equal to bins x bins!" << endl;
        return false;
    }

    // Init columns
    string newName("HEALPIX_EXTENSION");
    uint nr_of_quantities = 5;
    vector<string> colName(nr_of_quantities * nr_spectral_bins, "");
    vector<string> colForm(nr_of_quantities * nr_spectral_bins, "");
    vector<string> colUnit(nr_of_quantities * nr_spectral_bins, "");

    for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
    {
        uint id = nr_of_quantities * i_spectral;
        char str_1[1024];
#ifdef WINDOWS
        sprintf_s(str_1,
                    "I_STOKES (WAVELENGTH = %e [m])",
                    wavelength_list_det[i_spectral]);
#else
        sprintf(str_1,
                "I_STOKES (WAVELENGTH = %e [m])",
                wavelength_list_det[i_spectral]);
#endif
        colName[id + 0] = str_1;
#ifdef WINDOWS
        sprintf_s(str_1,
                    "TAU (WAVELENGTH = %e [m])",
                    wavelength_list_det[i_spectral]);
#else
        sprintf(str_1,                
                    "TAU (WAVELENGTH = %e [m])",
                    wavelength_list_det[i_spectral]);
#endif
        colName[id + 1] = str_1;
#ifdef WINDOWS
        sprintf_s(str_1,
                    "N_H (WAVELENGTH = %e [m])",
                    wavelength_list_det[i_spectral]);
#else
        sprintf(str_1,
                    "N_H (WAVELENGTH = %e [m])",
                    wavelength_list_det[i_spectral]);
#endif
        colName[id + 2] = str_1;
#ifdef WINDOWS
        sprintf_s(str_1,
                    "N_ion (WAVELENGTH = %e [m])",
                    wavelength_list_det[i_spectral]);
#else
        sprintf(str_1,
                    "N_ion (WAVELENGTH = %e [m])",
                    wavelength_list_det[i_spectral]);
#endif
        colName[id + 3] = str_1;
#ifdef WINDOWS
        sprintf_s(str_1,
                    "N_el (WAVELENGTH = %e [m])",
                    wavelength_list_det[i_spectral]);
#else
        sprintf(str_1,
                    "N_el (WAVELENGTH = %e [m])",
                    wavelength_list_det[i_spectral]);
#endif
        colName[id + 4] = str_1;
        
        colForm[id + 0] = "D";
        colForm[id + 1] = "D";
        colForm[id + 2] = "D";
        colForm[id + 3] = "D";
        colForm[id + 4] = "D";

        colUnit[id + 0] = "Jy/px";
        colUnit[id + 1] = "";
        colUnit[id + 2] = "m^-2";
        colUnit[id + 3] = "m^-2";
        colUnit[id + 4] = "m^-2";
    }

    CCfits::Table * newTable = pFits->addTable(newName, nelements, colName, colForm, colUnit);

    for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
    {
        valarray<double> array_I(nelements);
        valarray<double> array_T(nelements);
        valarray<double> array_NH(nelements);
        valarray<double> array_NI(nelements);
        valarray<double> array_NE(nelements);

        uint i = 0;
        for(uint i_y = 0; i_y < bins_y; i_y++)
        {
            for(uint i_x = 0; i_x < bins_x; i_x++)
            {
                array_I[i] = matrixI[i_spectral](i_x, i_y);
                array_T[i] = matrixT[i_spectral](i_x, i_y);
                
                array_NH[i] = matrixS1[i_spectral](i_x, i_y);
                array_NE[i] = matrixS2[i_spectral](i_x, i_y);
                array_NI[i] = matrixS3[i_spectral](i_x, i_y);
                
                i++;
            }
        }
            
        uint id = nr_of_quantities * i_spectral;
        newTable->column(colName[id + 0]).write(array_I, 1);
        newTable->column(colName[id + 1]).write(array_T, 1);
        newTable->column(colName[id + 2]).write(array_NH, 1);
        newTable->column(colName[id + 3]).write(array_NE, 1);
        newTable->column(colName[id + 4]).write(array_NI, 1);
    }
        
    // Healpix parameter
    newTable->addKey("PIXTYPE", "HEALPIX", "Pixel algorigthm");
    newTable->addKey("ORDERING", "RING", "Ordering scheme");
    newTable->addKey("INDXSCHM", "IMPLICIT", "Indexing scheme");
    newTable->addKey("NSIDE", uint(sqrt(bins_x / 12.0)), "Resolution Parameter");
    newTable->addKey("FIRSTPIX", 0, "First pixel (0 based)");
    newTable->addKey("LASTPIX", max_cells - 1, "Last pixel (0 based)");
    newTable->addKey("CROTA2", 0, "Rotation Angle (Degrees)");

    newTable->addKey("ETYPE", "free-free emission", "type of emission");
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
    newTable->addKey("BUBBLE_RADIUS", rad_bubble, "bubble radius [m]");

    cout << CLR_LINE << flush;

    return true;
}

// sync maps

bool CDetector::writeSyncMap(uint nr)
{
    cout << CLR_LINE << flush;
    cout << " -> Writing synchrotron plane detector ...  \r" << flush;
    // auto_ptr<CCfits::FITS> pFits(0);
    unique_ptr<CCfits::FITS> pFits;

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

        string path_out = path + str_end + FITS_COMPRESS_EXT;
        remove(path_out.c_str());

        long naxis = 4;
        long naxes[4] = { bins_x, bins_y, nr_spectral_bins, 6 * nr_extra };
        pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
    }
    catch(CCfits::FITS::CantCreate)
    {
        return false;
    }

    long nelements = bins_x * bins_y;
    if(max_cells != nelements)
    {
        cout << WARNING_LINE << "Max cells are not equal to bins x bins!" << endl;
        return false;
    }

    vector<long> fpixel(4);
    fpixel[0] = 1;
    fpixel[1] = 1;

    for(uint i_electrons = 0; i_electrons < 2; i_electrons++)
    {
        for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
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
                    array_I[i] = matrixI[i_spectral + i_electrons * nr_spectral_bins](i_x, i_y);
                    array_Q[i] = matrixQ[i_spectral + i_electrons * nr_spectral_bins](i_x, i_y);
                    array_U[i] = matrixU[i_spectral + i_electrons * nr_spectral_bins](i_x, i_y);
                    array_V[i] = matrixV[i_spectral + i_electrons * nr_spectral_bins](i_x, i_y);
                    array_T[i] = matrixT[i_spectral + i_electrons * nr_spectral_bins](i_x, i_y);
                    array_S[i] = matrixS1[i_spectral + i_electrons * nr_spectral_bins](i_x, i_y);
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
    calcCoordinateParameters(sidelength_x,
                                bins_x,
                                map_shift_x,
                                distance,
                                bin_width_x,
                                first_pix_val_x,
                                deg_per_pix_x,
                                first_pix_val_deg_x);
    double bin_width_y, first_pix_val_y, deg_per_pix_y, first_pix_val_deg_y;
    calcCoordinateParameters(sidelength_y,
                                bins_y,
                                map_shift_y,
                                distance,
                                bin_width_y,
                                first_pix_val_y,
                                deg_per_pix_y,
                                first_pix_val_deg_y);

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

    pFits->pHDU().addKey("CUNIT4",
                            "CR only (I, Q, U, V [Jy/px], 0, Nth [m^-2]), CR+TH (I, Q, U, V [Jy/px], "
                            "lambda^2 RM [rad], Ncr [m^-2])",
                            "unit of axis 4");
    pFits->pHDU().addKey("ETYPE", "synchotron emission", "type of emission");

    pFits->pHDU().addKey("ID", nr, "detector ID");

    for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
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

bool CDetector::writeSyncHealMap(uint nr)
{
    cout << CLR_LINE << flush;
    cout << " -> Writing synchrotron healpix detector ...  \r" << flush;
    // auto_ptr<CCfits::FITS> pFits(0);
    unique_ptr<CCfits::FITS> pFits;

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

        string path_out = path + str_end + FITS_COMPRESS_EXT;
        remove(path_out.c_str());

        long naxis = 1;
        long naxes[1] = { 0 };
        pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
    }
    catch(CCfits::FITS::CantCreate)
    {
        return false;
    }

    long nelements = bins_x * bins_y;
    if(max_cells != nelements)
    {
        cout << WARNING_LINE << "Max cells are not equal to bins x bins!" << endl;
        return false;
    }

    // Init columns
    string newName("HEALPIX_EXTENSION");
    uint nr_of_quantities = 5;
    vector<string> colName(2 * nr_of_quantities * nr_spectral_bins + 2, "");
    vector<string> colForm(2 * nr_of_quantities * nr_spectral_bins + 2, "");
    vector<string> colUnit(2 * nr_of_quantities * nr_spectral_bins + 2, "");

    strlist e_description;
    e_description.push_back("THERMAL_ELECTRONS");
    e_description.push_back("COSMIC_RAYS");

    for(uint i_electrons = 0; i_electrons < 2; i_electrons++)
    {
        for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
        {
            uint id = i_electrons * nr_of_quantities * nr_spectral_bins + nr_of_quantities * i_spectral;
            char str_1[1024];
#ifdef WINDOWS
            sprintf_s(str_1,
                        "I_STOKES (WAVELENGTH = %e [m], %s)",
                        wavelength_list_det[i_spectral],
                        e_description[i_electrons].c_str());
#else
            sprintf(str_1,
                    "I_STOKES (WAVELENGTH = %e [m], %s)",
                    wavelength_list_det[i_spectral],
                    e_description[i_electrons].c_str());
#endif
            colName[id + 0] = str_1;
#ifdef WINDOWS
            sprintf_s(str_1,
                        "Q_STOKES (WAVELENGTH = %e [m], %s)",
                        wavelength_list_det[i_spectral],
                        e_description[i_electrons].c_str());
#else
            sprintf(str_1,
                    "Q_STOKES (WAVELENGTH = %e [m], %s)",
                    wavelength_list_det[i_spectral],
                    e_description[i_electrons].c_str());
#endif
            colName[id + 1] = str_1;
#ifdef WINDOWS
            sprintf_s(str_1,
                        "U_STOKES (WAVELENGTH = %e [m], %s)",
                        wavelength_list_det[i_spectral],
                        e_description[i_electrons].c_str());
#else
            sprintf(str_1,
                    "U_STOKES (WAVELENGTH = %e [m], %s)",
                    wavelength_list_det[i_spectral],
                    e_description[i_electrons].c_str());
#endif
            colName[id + 2] = str_1;
#ifdef WINDOWS
            sprintf_s(str_1,
                        "V_STOKES (WAVELENGTH = %e [m], %s)",
                        wavelength_list_det[i_spectral],
                        e_description[i_electrons].c_str());
#else
            sprintf(str_1,
                    "V_STOKES (WAVELENGTH = %e [m], %s)",
                    wavelength_list_det[i_spectral],
                    e_description[i_electrons].c_str());
#endif
            colName[id + 3] = str_1;
#ifdef WINDOWS
            sprintf_s(str_1,
                        "FARADY ROTATION (WAVELENGTH = %e [m], %s)",
                        wavelength_list_det[i_spectral],
                        e_description[i_electrons].c_str());
#else
            sprintf(str_1,
                    "WAVELENGTH^2*FARADY_ROTATION (WAVELENGTH = %e [m], %s)",
                    wavelength_list_det[i_spectral],
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
            colUnit[id + 4] = "rad";
        }
        char str_1[1024];
#ifdef WINDOWS
        sprintf_s(str_1, "COLUMN_DENSITY (%s)", e_description[i_electrons].c_str());
#else
        sprintf(str_1, "COLUMN_DENSITY (%s)", e_description[i_electrons].c_str());
#endif
        colName[2 * nr_of_quantities * nr_spectral_bins + i_electrons] = str_1;
        colForm[2 * nr_of_quantities * nr_spectral_bins + i_electrons] = "D";
        colUnit[2 * nr_of_quantities * nr_spectral_bins + i_electrons] = "m^-2";
    }

    CCfits::Table * newTable = pFits->addTable(newName, nelements, colName, colForm, colUnit);

    for(uint i_electrons = 0; i_electrons < 2; i_electrons++)
    {
        for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
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
                    array_I[i] = matrixI[i_spectral + i_electrons * nr_spectral_bins](i_x, i_y);
                    array_Q[i] = matrixQ[i_spectral + i_electrons * nr_spectral_bins](i_x, i_y);
                    array_U[i] = matrixU[i_spectral + i_electrons * nr_spectral_bins](i_x, i_y);
                    array_V[i] = matrixV[i_spectral + i_electrons * nr_spectral_bins](i_x, i_y);
                    array_T[i] = matrixT[i_spectral + i_electrons * nr_spectral_bins](i_x, i_y);
                    i++;
                }
            uint id = i_electrons * nr_of_quantities * nr_spectral_bins + nr_of_quantities * i_spectral;
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
                array_C[i] = matrixS1[0 + i_electrons * nr_spectral_bins](i_x, i_y);
                i++;
            }
        newTable->column(colName[2 * nr_of_quantities * nr_spectral_bins + i_electrons])
            .write(array_C, 1);
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
    newTable->addKey("BUBBLE_RADIUS", rad_bubble, "bubble radius [m]");

    cout << CLR_LINE << flush;

    return true;
}

bool CDetector::writeLineSpectrum(CGasMixture * gas, uint i_species, uint i_line)
{
    cout << CLR_LINE << flush;
    cout << " -> Writing line spectrum ...  \r" << flush;
    // auto_ptr<CCfits::FITS> pFits(0);
    unique_ptr<CCfits::FITS> pFits;

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

        string path_out = path + str_end + FITS_COMPRESS_EXT;
        remove(path_out.c_str());

        long naxis = 3;
        long naxes[3] = { nr_spectral_bins, 1, 5 };
        pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
    }
    catch(CCfits::FITS::CantCreate)
    {
        return false;
    }

    long nelements = uint(nr_spectral_bins);

    vector<long> fpixel(3);
    fpixel[0] = 1;
    fpixel[1] = 1;

    std::valarray<double> array_I(nelements);
    std::valarray<double> array_Q(nelements);
    std::valarray<double> array_U(nelements);
    std::valarray<double> array_V(nelements);
    std::valarray<double> array_T(nelements);
    for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
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

    // Frequency
    pFits->pHDU().addKey("CTYPE1", "VELO", "type of unit 1");
    pFits->pHDU().addKey("CRVAL1", min_velocity, "value of axis 1");
    pFits->pHDU().addKey("CRPIX1", 1, "pixel where CRVAL1 is defined ");
    pFits->pHDU().addKey("CDELT1", channel_width, "delta of axis 1");
    pFits->pHDU().addKey("CUNIT1", "m/s", "unit of axis 1");

    // Simulated quantities
    pFits->pHDU().addKey("CTYPE2", "PARAM", "type of unit 2");
    pFits->pHDU().addKey("CRVAL2", 1, "value of axis 2");
    pFits->pHDU().addKey("CRPIX2", 1, "pixel where CRVAL2 is defined ");
    pFits->pHDU().addKey("CDELT2", 1, "delta of axis 2");
    pFits->pHDU().addKey("CUNIT2", "I, Q, U, V [Jy/px], optical depth", "unit of axis 2");

    pFits->pHDU().addKey(
        "GAS_SPECIES", gas->getGasSpeciesName(i_species), "name of the observed gas species");
    pFits->pHDU().addKey("TRANS", i_trans + 1, "transition index number (see leiden database)");
    pFits->pHDU().addKey("LEVEL_UPPER",
                            gas->getUpperEnergyLevel(i_species, i_trans),
                            "upper energy level index number (see leiden database)");
    pFits->pHDU().addKey("LEVEL_LOWER",
                            gas->getLowerEnergyLevel(i_species, i_trans),
                            "lower energy level index number (see leiden database)");
    pFits->pHDU().addKey(
        "FREQ", gas->getTransitionFrequency(i_species, i_trans), "frequency of the simulated transition");
    pFits->pHDU().addKey("CHANNELS", nr_spectral_bins, "number of velocity channels");
    pFits->pHDU().addKey("MINVEL", min_velocity, "velocity of the velocity channels (minvel to maxvel)");
    pFits->pHDU().addKey("MAXVEL", max_velocity1, "velocity of the velocity channels (minvel to maxvel)");
    pFits->pHDU().addKey("ZEEMAN",
                            gas->isTransZeemanSplit(i_species, i_trans),
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

bool CDetector::writeOPIATESpectrum(COpiateDataBase *op, uint det_id)
{
    cout << CLR_LINE << flush;
    cout << " -> Writing line spectrum ...  \r" << flush;
    // auto_ptr<CCfits::FITS> pFits(0);
    unique_ptr<CCfits::FITS> pFits;

    try
    {
        char str_tmp[1024];
        char str_end[1024];

#ifdef WINDOWS
        strcpy_s(str_tmp, "line_spectrum_species_%04d_line_%04d");
        sprintf_s(str_end, str_tmp, i_species + 1, i_line + 1);
#else
        strcpy(str_tmp, "opiate_spectrum_%04d");
        sprintf(str_end, str_tmp, det_id+1);
#endif

        string path_out = path + str_end + FITS_COMPRESS_EXT;
        remove(path_out.c_str());

        long naxis = 3;
        long naxes[3] = { nr_spectral_bins, 1, 5 };
        pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
    }
    catch(CCfits::FITS::CantCreate)
    {
        return false;
    }

    long nelements = uint(nr_spectral_bins);

    vector<long> fpixel(3);
    fpixel[0] = 1;
    fpixel[1] = 1;

    std::valarray<double> array_I(nelements);
    std::valarray<double> array_Q(nelements);
    std::valarray<double> array_U(nelements);
    std::valarray<double> array_V(nelements);
    std::valarray<double> array_T(nelements);
    for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
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

    // Frequency
    pFits->pHDU().addKey("CTYPE1", "VELO", "type of unit 1");
    pFits->pHDU().addKey("CRVAL1", min_velocity, "value of axis 1");
    pFits->pHDU().addKey("CRPIX1", 1, "pixel where CRVAL1 is defined ");
    pFits->pHDU().addKey("CDELT1", channel_width, "delta of axis 1");
    pFits->pHDU().addKey("CUNIT1", "m/s", "unit of axis 1");

    // Simulated quantities
    pFits->pHDU().addKey("CTYPE2", "PARAM", "type of unit 2");
    pFits->pHDU().addKey("CRVAL2", 1, "value of axis 2");
    pFits->pHDU().addKey("CRPIX2", 1, "pixel where CRVAL2 is defined ");
    pFits->pHDU().addKey("CDELT2", 1, "delta of axis 2");
    pFits->pHDU().addKey("CUNIT2", "I, Q, U, V [Jy/px], optical depth", "unit of axis 2");

    pFits->pHDU().addKey("GAS_SPECIES", op->getCurrentName(), "name of the observed gas species");
    pFits->pHDU().addKey("TRANS", i_trans + 1, "transition index number (see leiden database)");
    pFits->pHDU().addKey("FREQ", op->getCurrentFrequency(), "frequency of the simulated transition");
    pFits->pHDU().addKey("CHANNELS", nr_spectral_bins, "number of velocity channels");
    pFits->pHDU().addKey("MINVEL", min_velocity, "velocity of the velocity channels (minvel to maxvel)");
    pFits->pHDU().addKey("MAXVEL", max_velocity1, "velocity of the velocity channels (minvel to maxvel)");
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

bool CDetector::writeOPIATEVelChannelMaps(COpiateDataBase * op, uint det_id)
{
    cout << CLR_LINE << flush;
    cout << " -> Writing velocity channel map(s) ...  \r" << flush;

    for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
    {
        long naxis = 3;
        long naxes[3] = { uint(bins_x), uint(bins_y), 5 };

        // auto_ptr<CCfits::FITS> pFits(0);
        unique_ptr<CCfits::FITS> pFits;

        try
        {
            char str_tmp[1024];
            char str_end[1024];

#ifdef WINDOWS
            strcpy(str_tmp, "opiate_channel_map_%04d");
            sprintf(str_end, str_tmp, det_id);
#else
            strcpy(str_tmp, "opiate_channel_map_%04d_vel_%04d");
            sprintf(str_end, str_tmp, det_id+1, i_spectral+1);
#endif

            string path_out = path + str_end + FITS_COMPRESS_EXT;
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
            cout << WARNING_LINE << "Max cells are not equal to bins x bins!" << endl;
            return false;
        }

        vector<long> fpixel(3);
        fpixel[0] = 1;
        fpixel[1] = 1;
        fpixel[2] = 1;

        // cout << " -> Writing velocity channel map(s): " << int(100.0 * i_spectral / (nr_spectral_bins - 1))
        //      << "%     \r" << flush;

        std::valarray<double> array_I(nelements);
        std::valarray<double> array_Q(nelements);
        std::valarray<double> array_U(nelements);
        std::valarray<double> array_V(nelements);
        std::valarray<double> array_T(nelements);

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
        calcCoordinateParameters(sidelength_x,
                                    bins_x,
                                    map_shift_x,
                                    distance,
                                    bin_width_x,
                                    first_pix_val_x,
                                    deg_per_pix_x,
                                    first_pix_val_deg_x);
        double bin_width_y, first_pix_val_y, deg_per_pix_y, first_pix_val_deg_y;
        calcCoordinateParameters(sidelength_y,
                                    bins_y,
                                    map_shift_y,
                                    distance,
                                    bin_width_y,
                                    first_pix_val_y,
                                    deg_per_pix_y,
                                    first_pix_val_deg_y);

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
        pFits->pHDU().addKey(
            "CUNIT3", "I, Q, U, V [Jy/px], optical depth, column density [m^-2]", "unit of axis 3");

        pFits->pHDU().addKey(
            "GAS_SPECIES", op->getCurrentName(), "name of the observed gas species");
        pFits->pHDU().addKey("TRANS", i_trans + 1, "transition index number (see leiden database)");
        pFits->pHDU().addKey("FREQ", op->getCurrentFrequency(),
                                "frequency of the simulated transition");
        pFits->pHDU().addKey("VCH", i_spectral, "current velocity channel");
        pFits->pHDU().addKey("CHANNELS", nr_spectral_bins, "number of velocity channels");
        pFits->pHDU().addKey("MINVEL", min_velocity, "velocity of the velocity channels (minvel to maxvel)");
        pFits->pHDU().addKey("MAXVEL", max_velocity1, "velocity of the velocity channels (minvel to maxvel)");
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

    long naxis = 3;
    long naxes[3] = { uint(bins_x), uint(bins_y), nr_of_quantities };

    // auto_ptr<CCfits::FITS> pFits(0);
    unique_ptr<CCfits::FITS> pFits;

    try
    {
        char str_tmp[1024];
        char str_end[1024];

#ifdef WINDOWS
        strcpy_s(str_tmp, "channel_maps_species_%04d_line_%04d_extra");
        sprintf_s(str_end, str_tmp, i_species + 1, i_line + 1);
#else
        strcpy(str_tmp, "opiate_channel_map_%04d_extra");
        sprintf(str_end, str_tmp, det_id+1);
#endif

        string path_out = path + str_end + FITS_COMPRESS_EXT;
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
        cout << WARNING_LINE << "Max cells are not equal to bins x bins!" << endl;
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
                    // column density of the total gas
                    array_S[i] = matrixS1[0](i_x, i_y);
                }
                else if(i_extra == 1)
                {
                    // intensity weighted LOS magnetic field
                    if(nr_spectral_bins >= 2 && matrixS1[2](i_x, i_y) > 0)
                        array_S[i] = matrixS1[0](i_x, i_y) / matrixS1[2](i_x, i_y);
                    else
                        array_S[i] = 0;
                }
                else if(i_extra == 2)
                {
                    // intensity weighted total magnetic field
                    if(nr_spectral_bins >= 2 && matrixS1[2](i_x, i_y) > 0)
                        array_S[i] = matrixS1[1](i_x, i_y) / matrixS1[2](i_x, i_y);
                    else
                        array_S[i] = 0;
                }
                else if(i_extra == 3)
                {
                    // density weighted LOS magnetic field
                    if(nr_spectral_bins >= 5 && matrixS1[5](i_x, i_y) > 0)
                        array_S[i] = matrixS1[3](i_x, i_y) / matrixS1[5](i_x, i_y);
                    else
                        array_S[i] = 0;
                }
                else if(i_extra == 4)
                {
                    // density weighted magnetic field
                    if(nr_spectral_bins >= 5 && matrixS1[5](i_x, i_y) > 0)
                        array_S[i] = matrixS1[4](i_x, i_y) / matrixS1[5](i_x, i_y);
                    else
                        array_S[i] = 0;
                }
                i++;
            }
        fpixel[2] = i_extra + 1;
        pFits->pHDU().write(fpixel, nelements, array_S);
    }

    double bin_width_x, first_pix_val_x, deg_per_pix_x, first_pix_val_deg_x;
    calcCoordinateParameters(sidelength_x,
                                bins_x,
                                map_shift_x,
                                distance,
                                bin_width_x,
                                first_pix_val_x,
                                deg_per_pix_x,
                                first_pix_val_deg_x);
    double bin_width_y, first_pix_val_y, deg_per_pix_y, first_pix_val_deg_y;
    calcCoordinateParameters(sidelength_y,
                                bins_y,
                                map_shift_y,
                                distance,
                                bin_width_y,
                                first_pix_val_y,
                                deg_per_pix_y,
                                first_pix_val_deg_y);

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
        pFits->pHDU().addKey("CUNIT3",
                                "Column density [m^-2], intensity weighted B_LOS, B, "
                                "dens. weighted B_LOS, B",
                                "unit of axis 3");
    else
        pFits->pHDU().addKey("CUNIT3", "Column density [m^-2]", "unit of axis 3");

    pFits->pHDU().addKey(
        "GAS_SPECIES", op->getCurrentName(), "name of the observed gas species");
    pFits->pHDU().addKey("TRANS", i_trans + 1, "transition index number (see leiden database)");
    pFits->pHDU().addKey("FREQ", op->getCurrentFrequency(), "frequency of the simulated transition");
    pFits->pHDU().addKey("CHANNELS", nr_spectral_bins, "number of velocity channels");
    pFits->pHDU().addKey("MINVEL", min_velocity, "velocity of the velocity channels (minvel to maxvel)");
    pFits->pHDU().addKey("MAXVEL", max_velocity1, "velocity of the velocity channels (minvel to maxvel)");
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

bool CDetector::writeOPIATEIntChannelMaps(COpiateDataBase * op, uint det_id)
{
    cout << CLR_LINE << flush;
    cout << " -> Writing integrated velocity channel map(s) ...  \r" << flush;

    long naxis = 3;
    long naxes[3] = { uint(bins_x), uint(bins_y), 6 };

    // auto_ptr<CCfits::FITS> pFits(0);
    unique_ptr<CCfits::FITS> pFits;

    try
    {
        char str_tmp[1024];
        char str_end[1024];

#ifdef WINDOWS
        strcpy_s(str_tmp, "int_channel_map_species_%04d_line_%04d");
        sprintf_s(str_end, str_tmp, i_species + 1, i_line + 1);
#else
        strcpy(str_tmp, "opiate_int_map__%04d");
        sprintf(str_end, str_tmp, det_id+1);
#endif

        string path_out = path + str_end + FITS_COMPRESS_EXT;
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
        cout << WARNING_LINE << "Max cells are not equal to bins x bins!" << endl;
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
            for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
            {
                array_I[i] +=
                    matrixI[i_spectral](i_x, i_y) * (channel_width / nr_spectral_bins);// * 1e-3;
                array_Q[i] +=
                    matrixQ[i_spectral](i_x, i_y) * (channel_width / nr_spectral_bins);// * 1e-3;
                array_U[i] +=
                    matrixU[i_spectral](i_x, i_y) * (channel_width / nr_spectral_bins);// * 1e-3;
                array_V[i] +=
                    matrixV[i_spectral](i_x, i_y) * (channel_width / nr_spectral_bins);// * 1e-3;
            }
            array_T[i] = matrixT[int(nr_spectral_bins / 2.0)](i_x, i_y);

            if(nr_spectral_bins >= 2 && matrixS1[2](i_x, i_y) > 0)
            {
                // LOS magnetic field strength
                array_S[i] = matrixS1[0](i_x, i_y) / matrixS1[2](i_x, i_y);
            }
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
    calcCoordinateParameters(sidelength_x,
                                bins_x,
                                map_shift_x,
                                distance,
                                bin_width_x,
                                first_pix_val_x,
                                deg_per_pix_x,
                                first_pix_val_deg_x);
    double bin_width_y, first_pix_val_y, deg_per_pix_y, first_pix_val_deg_y;
    calcCoordinateParameters(sidelength_y,
                                bins_y,
                                map_shift_y,
                                distance,
                                bin_width_y,
                                first_pix_val_y,
                                deg_per_pix_y,
                                first_pix_val_deg_y);

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
    pFits->pHDU().addKey(
        "CUNIT3", "I, Q, U, V [Jy/px], optical depth, column density [m^-2]", "unit of axis 3");

    pFits->pHDU().addKey("GAS_SPECIES", op->getCurrentName(), "name of the observed gas species");
    pFits->pHDU().addKey("TRANS", i_trans + 1, "transition index number (see leiden database)");
    pFits->pHDU().addKey("FREQ", op->getCurrentFrequency(), "frequency of the simulated transition");
    pFits->pHDU().addKey("CHANNELS", nr_spectral_bins, "number of velocity channels");
    pFits->pHDU().addKey("MINVEL", min_velocity, "velocity of the velocity channels (minvel to maxvel)");
    pFits->pHDU().addKey("MAXVEL", max_velocity1, "velocity of the velocity channels (minvel to maxvel)");
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

bool CDetector::writeVelChannelMaps(CGasMixture * gas, uint i_species, uint i_line)
{
    cout << CLR_LINE << flush;
    cout << " -> Writing velocity channel map(s) ...  \r" << flush;

    for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
    {
        long naxis = 3;
        long naxes[3] = { bins_x, bins_y, 5 };

        // auto_ptr<CCfits::FITS> pFits(0);
        unique_ptr<CCfits::FITS> pFits;

        try
        {
            char str_tmp[1024];
            char str_end[1024];

#ifdef WINDOWS
            strcpy_s(str_tmp, "channel_maps_species_%04d_line_%04d_vel_%04d");
            sprintf_s(str_end, str_tmp, i_species + 1, i_line + 1, i_spectral + 1);
#else
            strcpy(str_tmp, "channel_maps_species_%04d_line_%04d_vel_%04d");
            sprintf(str_end, str_tmp, i_species + 1, i_line + 1, i_spectral + 1);
#endif

            string path_out = path + str_end + FITS_COMPRESS_EXT;
            remove(path_out.c_str());
            pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
        }
        catch(CCfits::FITS::CantCreate)
        {
            return false;
        }

        long nelements = bins_x * bins_y;
        if(max_cells != nelements)
        {
            cout << WARNING_LINE << "Max cells are not equal to bins x bins!" << endl;
            return false;
        }

        vector<long> fpixel(3);
        fpixel[0] = 1;
        fpixel[1] = 1;
        fpixel[2] = 1;

        // cout << " -> Writing velocity channel map(s): " << int(100.0 * i_spectral / (nr_spectral_bins - 1))
        //      << "%     \r" << flush;

        std::valarray<double> array_I(nelements);
        std::valarray<double> array_Q(nelements);
        std::valarray<double> array_U(nelements);
        std::valarray<double> array_V(nelements);
        std::valarray<double> array_T(nelements);

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
        calcCoordinateParameters(sidelength_x,
                                    bins_x,
                                    map_shift_x,
                                    distance,
                                    bin_width_x,
                                    first_pix_val_x,
                                    deg_per_pix_x,
                                    first_pix_val_deg_x);
        double bin_width_y, first_pix_val_y, deg_per_pix_y, first_pix_val_deg_y;
        calcCoordinateParameters(sidelength_y,
                                    bins_y,
                                    map_shift_y,
                                    distance,
                                    bin_width_y,
                                    first_pix_val_y,
                                    deg_per_pix_y,
                                    first_pix_val_deg_y);

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
        pFits->pHDU().addKey(
            "CUNIT3", "I, Q, U, V [Jy/px], optical depth", "unit of axis 3");

        pFits->pHDU().addKey(
            "GAS_SPECIES", gas->getGasSpeciesName(i_species), "name of the observed gas species");
        pFits->pHDU().addKey("TRANS", i_trans + 1, "transition index number (see leiden database)");
        pFits->pHDU().addKey("LEVEL_UPPER",
                                gas->getUpperEnergyLevel(i_species, i_trans) + 1,
                                "upper energy level index number (see leiden database)");
        pFits->pHDU().addKey("LEVEL_LOWER",
                                gas->getLowerEnergyLevel(i_species, i_trans) + 1,
                                "lower energy level index number (see leiden database)");
        pFits->pHDU().addKey("FREQ",
                                gas->getTransitionFrequency(i_species, i_trans),
                                "frequency of the simulated transition");
        pFits->pHDU().addKey("VCH", i_spectral, "current velocity channel");
        pFits->pHDU().addKey("CHANNELS", nr_spectral_bins, "number of velocity channels");
        pFits->pHDU().addKey("MINVEL", min_velocity, "velocity of the velocity channels (minvel to maxvel)");
        pFits->pHDU().addKey("MAXVEL", max_velocity1, "velocity of the velocity channels (minvel to maxvel)");
        pFits->pHDU().addKey("ZEEMAN",
                                gas->isTransZeemanSplit(i_species, i_trans),
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
    uint nr_of_quantities = 2;
    bool isZeeman=gas->isTransZeemanSplit(i_species, i_trans);

    if(isZeeman==true)
        nr_of_quantities = 6;

    long naxis = 3;
    long naxes[3] = { bins_x, bins_y, nr_of_quantities };

    // auto_ptr<CCfits::FITS> pFits(0);
    unique_ptr<CCfits::FITS> pFits;

    try
    {
        char str_tmp[1024];
        char str_end[1024];

#ifdef WINDOWS
        strcpy_s(str_tmp, "channel_maps_species_%04d_line_%04d_extra");
        sprintf_s(str_end, str_tmp, i_species + 1, i_line + 1);
#else
        strcpy(str_tmp, "channel_maps_species_%04d_line_%04d_extra");
        sprintf(str_end, str_tmp, i_species + 1, i_line + 1);
#endif

        string path_out = path + str_end + FITS_COMPRESS_EXT;
        remove(path_out.c_str());
        pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
    }
    catch(CCfits::FITS::CantCreate)
    {
        return false;
    }

    long nelements = bins_x * bins_y;
    if(max_cells != nelements)
    {
        cout << WARNING_LINE << "Max cells are not equal to bins x bins!" << endl;
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
                    // column density of the total gas
                    if(gas->isTransZeemanSplit(i_species, i_trans))
                        array_S[i] = matrixS1[5](i_x, i_y);
                    else
                        array_S[i] = matrixS1[0](i_x, i_y);
                }
                else if(i_extra == 1)
                {
                    // intensity weighted LOS magnetic field
                    if(nr_spectral_bins >= 2 && matrixS1[2](i_x, i_y) > 0)
                        array_S[i] = matrixS1[0](i_x, i_y) / matrixS1[2](i_x, i_y);
                    else
                        array_S[i] = 0;
                }
                else if(i_extra == 2)
                {
                    // intensity weighted total magnetic field
                    if(nr_spectral_bins >= 2 && matrixS1[2](i_x, i_y) > 0)
                        array_S[i] = matrixS1[1](i_x, i_y) / matrixS1[2](i_x, i_y);
                    else
                        array_S[i] = 0;
                }
                else if(i_extra == 3)
                {
                    // density weighted LOS magnetic field
                    if(nr_spectral_bins >= 5 && matrixS1[5](i_x, i_y) > 0)
                        array_S[i] = matrixS1[3](i_x, i_y) / matrixS1[5](i_x, i_y);
                    else
                        array_S[i] = 0;
                }
                else if(i_extra == 4)
                {
                    // density weighted magnetic field
                    if(nr_spectral_bins >= 5 && matrixS1[5](i_x, i_y) > 0)
                        array_S[i] = matrixS1[4](i_x, i_y) / matrixS1[5](i_x, i_y);
                    else
                        array_S[i] = 0;
                }
                i++;
            }
        fpixel[2] = i_extra + 1;
        pFits->pHDU().write(fpixel, nelements, array_S);
    }

    double bin_width_x, first_pix_val_x, deg_per_pix_x, first_pix_val_deg_x;
    calcCoordinateParameters(sidelength_x,
                                bins_x,
                                map_shift_x,
                                distance,
                                bin_width_x,
                                first_pix_val_x,
                                deg_per_pix_x,
                                first_pix_val_deg_x);
    double bin_width_y, first_pix_val_y, deg_per_pix_y, first_pix_val_deg_y;
    calcCoordinateParameters(sidelength_y,
                                bins_y,
                                map_shift_y,
                                distance,
                                bin_width_y,
                                first_pix_val_y,
                                deg_per_pix_y,
                                first_pix_val_deg_y);

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
    if(gas->isTransZeemanSplit(i_species, i_trans))
        pFits->pHDU().addKey("CUNIT3",
                                "Column density [m^-2], Column density [m^-2], intensity weighted B_LOS, B, "
                                "dens. weighted B_LOS, B",
                                "unit of axis 3");
    else
        pFits->pHDU().addKey("CUNIT3", "Column density [m^-2], Column density [m^-2]", "unit of axis 3");

    pFits->pHDU().addKey(
        "GAS_SPECIES", gas->getGasSpeciesName(i_species), "name of the observed gas species");
    pFits->pHDU().addKey("TRANS", i_trans + 1, "transition index number (see leiden database)");
    pFits->pHDU().addKey("LEVEL_UPPER",
                            gas->getUpperEnergyLevel(i_species, i_trans),
                            "upper energy level index number (see leiden database)");
    pFits->pHDU().addKey("LEVEL_LOWER",
                            gas->getLowerEnergyLevel(i_species, i_trans),
                            "lower energy level index number (see leiden database)");
    pFits->pHDU().addKey(
        "FREQ", gas->getTransitionFrequency(i_species, i_trans), "frequency of the simulated transition");
    pFits->pHDU().addKey("CHANNELS", nr_spectral_bins, "number of velocity channels");
    pFits->pHDU().addKey("MINVEL", min_velocity, "velocity of the velocity channels (minvel to maxvel)");
    pFits->pHDU().addKey("MAXVEL", max_velocity1, "velocity of the velocity channels (minvel to maxvel)");
    pFits->pHDU().addKey("ZEEMAN",
                            gas->isTransZeemanSplit(i_species, i_trans),
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

bool CDetector::writeIntChannelMaps(CGasMixture * gas, uint i_species, uint i_line)
{
    cout << CLR_LINE << flush;
    cout << " -> Writing integrated velocity channel map(s) ...  \r" << flush;

    long naxis = 3;

    uint nr_of_quantities = 6;
    bool isZeeman=gas->isTransZeemanSplit(i_species, i_trans);

    if(isZeeman==true)
        nr_of_quantities = 10;

    long naxes[3] = { bins_x, bins_y, nr_of_quantities};

    // auto_ptr<CCfits::FITS> pFits(0);
    unique_ptr<CCfits::FITS> pFits;

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

        string path_out = path + str_end + FITS_COMPRESS_EXT;
        remove(path_out.c_str());
        pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
    }
    catch(CCfits::FITS::CantCreate)
    {
        return false;
    }

    long nelements = bins_x * bins_y;
    if(max_cells != nelements)
    {
        cout << WARNING_LINE << "Max cells are not equal to bins x bins!" << endl;
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
    //std::valarray<double> array_T(nelements);

    // column density of the gas species
    std::valarray<double> array_Nx(nelements);

    // column density of the total gas
    std::valarray<double> array_NH(nelements);

    // intensity weighted LOS magnetic field
    std::valarray<double> array_Bpar_I(nelements);

    // intensity weighted total magnetic field
    std::valarray<double> array_Btot_I(nelements);

    // intensity weighted LOS magnetic field
    std::valarray<double> array_Bpar_dens(nelements);

    // intensity weighted total magnetic field
    std::valarray<double> array_Btot_dens(nelements);

    uint i = 0;
    for(uint i_y = 0; i_y < bins_y; i_y++)
        for(uint i_x = 0; i_x < bins_x; i_x++)
        {
            for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
            {
                array_I[i] +=
                    matrixI[i_spectral](i_x, i_y) * (channel_width / nr_spectral_bins);// * 1e-3;
                array_Q[i] +=
                    matrixQ[i_spectral](i_x, i_y) * (channel_width / nr_spectral_bins);// * 1e-3;
                array_U[i] +=
                    matrixU[i_spectral](i_x, i_y) * (channel_width / nr_spectral_bins);// * 1e-3;
                array_V[i] +=
                    matrixV[i_spectral](i_x, i_y) * (channel_width / nr_spectral_bins);// * 1e-3;
            }

            if(isZeeman==true)
            {
                // column density of the total gas
                array_NH[i] = matrixS1[5](i_x, i_y);

                // column density of the gas species
                array_Nx[i] = matrixS1[6](i_x, i_y);

                // intensity weighted LOS magnetic field
                if(nr_spectral_bins >= 2 && matrixS1[2](i_x, i_y) > 0)
                    array_Bpar_I[i] = matrixS1[0](i_x, i_y) / matrixS1[2](i_x, i_y);

                // intensity weighted total magnetic field
                if(nr_spectral_bins >= 2 && matrixS1[2](i_x, i_y) > 0)
                    array_Btot_I[i] = matrixS1[1](i_x, i_y) / matrixS1[2](i_x, i_y);

                // density weighted LOS magnetic field
                if(nr_spectral_bins >= 5 && matrixS1[5](i_x, i_y) > 0)
                    array_Bpar_dens[i] = matrixS1[3](i_x, i_y) / matrixS1[5](i_x, i_y);

                // density weighted magnetic field
                if(nr_spectral_bins >= 5 && matrixS1[5](i_x, i_y) > 0)
                    array_Btot_dens[i] = matrixS1[4](i_x, i_y) / matrixS1[5](i_x, i_y);
            }
            else
            {
                // column density of the total gas
                array_NH[i] = matrixS1[0](i_x, i_y);

                // column density of the gas species
                array_Nx[i] = matrixS1[1](i_x, i_y);
            }

//                array_S[i] = matrixS[0](i_x, i_y);
//                array_T[i] = matrixT[int(nr_spectral_bins / 2.0)](i_x, i_y);
//
//                if(nr_spectral_bins >= 2 && matrixS[2](i_x, i_y) > 0)
//                {
//                    // LOS magnetic field strength
//                    array_S[i] = matrixS[0](i_x, i_y) / matrixS[2](i_x, i_y);
//                }
//                else
//                    array_S[i] = matrixS[0](i_x, i_y);
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

    if(isZeeman)
    {
        fpixel[2] = 5;
        pFits->pHDU().write(fpixel, nelements, array_Nx);

        fpixel[2] = 6;
        pFits->pHDU().write(fpixel, nelements, array_NH);

        fpixel[2] = 7;
        pFits->pHDU().write(fpixel, nelements, array_Bpar_I);

        fpixel[2] = 8;
        pFits->pHDU().write(fpixel, nelements, array_Btot_I);

        fpixel[2] = 9;
        pFits->pHDU().write(fpixel, nelements, array_Bpar_dens);

        fpixel[2] = 10;
        pFits->pHDU().write(fpixel, nelements, array_Btot_dens);
    }
    else
    {
        fpixel[2] = 5;
        pFits->pHDU().write(fpixel, nelements, array_Nx);

        fpixel[2] = 6;
        pFits->pHDU().write(fpixel, nelements, array_NH);
    }




//        fpixel[2] = 5;
//        pFits->pHDU().write(fpixel, nelements, array_S);
//        fpixel[2] = 6;
//        pFits->pHDU().write(fpixel, nelements, array_S);

    double bin_width_x, first_pix_val_x, deg_per_pix_x, first_pix_val_deg_x;
    calcCoordinateParameters(sidelength_x,
                                bins_x,
                                map_shift_x,
                                distance,
                                bin_width_x,
                                first_pix_val_x,
                                deg_per_pix_x,
                                first_pix_val_deg_x);
    double bin_width_y, first_pix_val_y, deg_per_pix_y, first_pix_val_deg_y;
    calcCoordinateParameters(sidelength_y,
                                bins_y,
                                map_shift_y,
                                distance,
                                bin_width_y,
                                first_pix_val_y,
                                deg_per_pix_y,
                                first_pix_val_deg_y);

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


    if(isZeeman)
    {
        pFits->pHDU().addKey("CUNIT3", "I Q U V [Jy/px], Nx NH [m^-3], Bp_I, Bt_I, Bp_d, Bt_d [T]", "unit of axis 3");
    }
    else
    {
        pFits->pHDU().addKey("CUNIT3", "I Q U V [Jy/px], Nx NH [m^-3] ", "unit of axis 3");
    }

    pFits->pHDU().addKey("GAS_SPECIES", gas->getGasSpeciesName(i_species), "name of the observed gas species");
    pFits->pHDU().addKey("TRANS", i_trans + 1, "transition index number (see leiden database)");
    pFits->pHDU().addKey("LEVEL_UPPER",
                            gas->getUpperEnergyLevel(i_species, i_trans),
                            "upper energy level index number (see leiden database)");
    pFits->pHDU().addKey("LEVEL_LOWER",
                            gas->getLowerEnergyLevel(i_species, i_trans),
                            "lower energy level index number (see leiden database)");
    pFits->pHDU().addKey("FREQ", gas->getTransitionFrequency(i_species, i_trans), "frequency of the simulated transition");
    pFits->pHDU().addKey("CHANNELS", nr_spectral_bins, "number of velocity channels");
    pFits->pHDU().addKey("MINVEL", min_velocity, "velocity of the velocity channels (minvel to maxvel)");
    pFits->pHDU().addKey("MAXVEL", max_velocity1, "velocity of the velocity channels (minvel to maxvel)");
    pFits->pHDU().addKey("ZEEMAN",isZeeman,
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

bool CDetector::writeOPIATEVelChannelHealMaps(COpiateDataBase * op, uint det_id)
{
    cout << CLR_LINE << flush;
    cout << " -> Writing velocity channel map(s) ...  \r" << flush;

    for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
    {
        long naxis = 1;
        long naxes[1] = { 0 };

        // auto_ptr<CCfits::FITS> pFits(0);
        unique_ptr<CCfits::FITS> pFits;

        try
        {
            char str_tmp[1024];
            char str_end[1024];

#ifdef WINDOWS
            strcpy_s(str_tmp, "channel_maps_species_%04d_line_%04d_vel_%04d");
            sprintf_s(str_end, str_tmp, i_species + 1, i_line + 1, i_spectral + 1);
#else
            strcpy(str_tmp, "opiate_channel_map_%04d_vel_%04d");
            sprintf(str_end, str_tmp, det_id+1, i_spectral+1);
#endif

            string path_out = path + str_end + FITS_COMPRESS_EXT;
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
            cout << WARNING_LINE << "Max cells are not equal to bins x bins!" << endl;
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

        CCfits::Table * newTable = pFits->addTable(newName, nelements, colName, colForm, colUnit);

        // cout << " -> Writing velocity channel map(s): " << int(100.0 * i_spectral / (nr_spectral_bins - 1))
        //      << "%     \r" << flush;

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

        newTable->addKey("GAS_SPECIES", op->getCurrentName(), "name of the observed gas species");
        newTable->addKey("TRANS", i_trans + 1, "transition index number (see leiden database)");
        newTable->addKey("FREQ", op->getCurrentFrequency(),"frequency of the simulated transition");
        newTable->addKey("VCH", i_spectral, "current velocity channel");
        newTable->addKey("CHANNELS", nr_spectral_bins, "number of velocity channels");
        newTable->addKey("MINVEL", min_velocity, "velocity of the velocity channels (minvel to maxvel)");
        newTable->addKey("MAXVEL", max_velocity1, "velocity of the velocity channels (minvel to maxvel)");
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
    long naxes[1] = { 0 };

    // auto_ptr<CCfits::FITS> pFits(0);
    unique_ptr<CCfits::FITS> pFits;

    try
    {
        char str_tmp[1024];
        char str_end[1024];

#ifdef WINDOWS
        strcpy_s(str_tmp, "channel_maps_species_%04d_line_%04d_extra");
        sprintf_s(str_end, str_tmp, i_species + 1, i_line + 1);
#else
        strcpy(str_tmp, "opiate_map_%04d_extra");
        sprintf(str_end, str_tmp, det_id);
#endif

        string path_out = path + str_end + FITS_COMPRESS_EXT;
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
        cout << WARNING_LINE << "Max cells are not equal to bins x bins!" << endl;
        return false;
    }

    // Init columns
    string newName("HEALPIX_EXTENSION");
    uint nr_of_quantities = 1;

    //if(gas->isTransZeemanSplit(i_species, i_trans))
    //    nr_of_quantities += 4;

    vector<string> colName(nr_of_quantities, "");
    vector<string> colForm(nr_of_quantities, "");
    vector<string> colUnit(nr_of_quantities, "");

    colName[0] = "COLUMN_DENSITY";
    colForm[0] = "D";
    colUnit[0] = "m^-2";

    CCfits::Table * newTable = pFits->addTable(newName, nelements, colName, colForm, colUnit);

    // Column density / extra data
    for(int i_extra = 0; i_extra < nr_of_quantities; i_extra++)
    {
        valarray<double> array_S(nelements);
        uint i = 0;
        for(uint i_y = 0; i_y < bins_y; i_y++)
            for(uint i_x = 0; i_x < bins_x; i_x++)
            {
                if(i_extra == 0)
                {
                    // column density of the total gas
                    array_S[i] = matrixS1[0](i_x, i_y);
                }
                else if(i_extra == 1)
                {
                    // intensity weighted LOS magnetic field
                    if(nr_spectral_bins >= 2 && matrixS1[2](i_x, i_y) > 0)
                        array_S[i] = matrixS1[0](i_x, i_y) / matrixS1[2](i_x, i_y);
                    else
                        array_S[i] = 0;
                }
                else if(i_extra == 2)
                {
                    // intensity weighted total magnetic field
                    if(nr_spectral_bins >= 2 && matrixS1[2](i_x, i_y) > 0)
                        array_S[i] = matrixS1[1](i_x, i_y) / matrixS1[2](i_x, i_y);
                    else
                        array_S[i] = 0;
                }
                else if(i_extra == 3)
                {
                    // density weighted LOS magnetic field
                    if(nr_spectral_bins >= 5 && matrixS1[5](i_x, i_y) > 0)
                        array_S[i] = matrixS1[3](i_x, i_y) / matrixS1[5](i_x, i_y);
                    else
                        array_S[i] = 0;
                }
                else if(i_extra == 4)
                {
                    // density weighted magnetic field
                    if(nr_spectral_bins >= 5 && matrixS1[5](i_x, i_y) > 0)
                        array_S[i] = matrixS1[4](i_x, i_y) / matrixS1[5](i_x, i_y);
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

    newTable->addKey("GAS_SPECIES", op->getCurrentName(), "name of the observed gas species");
    newTable->addKey("TRANS", i_trans + 1, "transition index number (see leiden database)");
    newTable->addKey("FREQ", op->getCurrentFrequency(), "frequency of the simulated transition");
    newTable->addKey("CHANNELS", nr_spectral_bins, "number of velocity channels");
    newTable->addKey("MINVEL", min_velocity, "velocity of the velocity channels (minvel to maxvel)");
    newTable->addKey("MAXVEL", max_velocity1, "velocity of the velocity channels (minvel to maxvel)");
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

bool CDetector::writeVelChannelHealMaps(CGasMixture * gas, uint i_species, uint i_line)
{
    cout << CLR_LINE << flush;
    cout << " -> Writing velocity channel map(s) ...  \r" << flush;
    bool hasZeeman = gas->isTransZeemanSplit(i_species, i_trans);

    for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
    {
        long naxis = 1;
        long naxes[1] = { 0 };

        // auto_ptr<CCfits::FITS> pFits(0);
        unique_ptr<CCfits::FITS> pFits;

        try
        {
            char str_tmp[1024];
            char str_end[1024];

#ifdef WINDOWS
            strcpy_s(str_tmp, "channel_maps_species_%04d_line_%04d_vel_%04d");
            sprintf_s(str_end, str_tmp, i_species + 1, i_line + 1, i_spectral + 1);
#else
            strcpy(str_tmp, "channel_maps_species_%04d_line_%04d_vel_%04d");
            sprintf(str_end, str_tmp, i_species + 1, i_line + 1, i_spectral + 1);
#endif

            string path_out = path + str_end + FITS_COMPRESS_EXT;
            remove(path_out.c_str());
            pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
        }
        catch(CCfits::FITS::CantCreate)
        {
            return false;
        }

        long nelements = bins_x * bins_y;
        if(max_cells != nelements)
        {
            cout << WARNING_LINE << "Max cells are not equal to bins x bins!" << endl;
            return false;
        }

        // Init columns
        string newName("HEALPIX_EXTENSION");
        uint nr_of_quantities = 2;
        
        if(hasZeeman)
            nr_of_quantities=5;
        
        vector<string> colName(nr_of_quantities, "");
        vector<string> colForm(nr_of_quantities, "");
        vector<string> colUnit(nr_of_quantities, "");

        if(hasZeeman)
        {
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
        }
        else
        {
            colName[0] = "I_STOKES";
            colName[1] = "TAU";

            colForm[0] = "D";
            colForm[1] = "D";

            colUnit[0] = "Jy/px";
            colUnit[1] = "";
        }
        
        CCfits::Table * newTable = pFits->addTable(newName, nelements, colName, colForm, colUnit);

        // cout << " -> Writing velocity channel map(s): " << int(100.0 * i_spectral / (nr_spectral_bins - 1))
        //      << "%     \r" << flush;

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
                
                if(hasZeeman)
                {
                    array_Q[i] = matrixQ[i_spectral](i_x, i_y);
                    array_U[i] = matrixU[i_spectral](i_x, i_y);
                    array_V[i] = matrixV[i_spectral](i_x, i_y);
                }
                
                array_T[i] = matrixT[i_spectral](i_x, i_y);
                i++;
            }

        newTable->column(colName[0]).write(array_I, 1);
        
        if(hasZeeman)
        {
            newTable->column(colName[1]).write(array_Q, 1);
            newTable->column(colName[2]).write(array_U, 1);
            newTable->column(colName[3]).write(array_V, 1);
            newTable->column(colName[4]).write(array_T, 1);
        }
        else
        {
            newTable->column(colName[1]).write(array_T, 1);
        }

        // Healpix parameter
        newTable->addKey("PIXTYPE", "HEALPIX", "Pixel algorigthm");
        newTable->addKey("ORDERING", "RING", "Ordering scheme");
        newTable->addKey("INDXSCHM", "IMPLICIT", "Indexing scheme");
        newTable->addKey("NSIDE", uint(sqrt(bins_x / 12.0)), "Resolution Parameter");
        newTable->addKey("FIRSTPIX", 0, "First pixel (0 based)");
        newTable->addKey("LASTPIX", max_cells - 1, "Last pixel (0 based)");
        newTable->addKey("CROTA2", 0, "Rotation Angle (Degrees)");

        newTable->addKey(
            "GAS_SPECIES", gas->getGasSpeciesName(i_species), "name of the observed gas species");
        newTable->addKey("TRANS", i_trans + 1, "transition index number (see leiden database)");
        newTable->addKey("LEVEL_UPPER",
                            gas->getUpperEnergyLevel(i_species, i_trans),
                            "upper energy level index number (see leiden database)");
        newTable->addKey("LEVEL_LOWER",
                            gas->getLowerEnergyLevel(i_species, i_trans),
                            "lower energy level index number (see leiden database)");
        newTable->addKey("FREQ",
                            gas->getTransitionFrequency(i_species, i_trans),
                            "frequency of the simulated transition");
        newTable->addKey("VCH", i_spectral, "current velocity channel");
        newTable->addKey("CHANNELS", nr_spectral_bins, "number of velocity channels");
        newTable->addKey("MINVEL", min_velocity, "velocity of the velocity channels (minvel to maxvel)");
        newTable->addKey("MAXVEL", max_velocity1, "velocity of the velocity channels (minvel to maxvel)");
        newTable->addKey("ZEEMAN",
                            gas->isTransZeemanSplit(i_species, i_trans),
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
    long naxes[1] = { 0 };

    // auto_ptr<CCfits::FITS> pFits(0);
    unique_ptr<CCfits::FITS> pFits;

    try
    {
        char str_tmp[1024];
        char str_end[1024];

#ifdef WINDOWS
        strcpy_s(str_tmp, "channel_maps_species_%04d_line_%04d_extra");
        sprintf_s(str_end, str_tmp, i_species + 1, i_line + 1);
#else
        strcpy(str_tmp, "channel_maps_species_%04d_line_%04d_extra");
        sprintf(str_end, str_tmp, i_species + 1, i_line + 1);
#endif

        string path_out = path + str_end + FITS_COMPRESS_EXT;
        remove(path_out.c_str());
        pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
    }

    catch(CCfits::FITS::CantCreate)
    {
        return false;
    }

    long nelements = bins_x * bins_y;
    if(max_cells != nelements)
    {
        cout << WARNING_LINE << "Max cells are not equal to bins x bins!" << endl;
        return false;
    }

    // Init columns
    string newName("HEALPIX_EXTENSION");
    uint nr_of_quantities = 1;

    if(gas->isTransZeemanSplit(i_species, i_trans))
        nr_of_quantities += 4;

    vector<string> colName(nr_of_quantities, "");
    vector<string> colForm(nr_of_quantities, "");
    vector<string> colUnit(nr_of_quantities, "");

    colName[0] = "COLUMN_DENSITY";
    colForm[0] = "D";
    colUnit[0] = "m^-2";

    if(gas->isTransZeemanSplit(i_species, i_trans))
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

    CCfits::Table * newTable = pFits->addTable(newName, nelements, colName, colForm, colUnit);

    // Column density / extra data
    for(uint i_extra = 0; i_extra < nr_of_quantities; i_extra++)
    {
        valarray<double> array_S(nelements);
        uint i = 0;
        for(uint i_y = 0; i_y < bins_y; i_y++)
            for(uint i_x = 0; i_x < bins_x; i_x++)
            {
                if(i_extra == 0)
                {
                    // column density of the total gas
                    if(gas->isTransZeemanSplit(i_species, i_trans))
                        array_S[i] = matrixS1[5](i_x, i_y);
                    else
                        array_S[i] = matrixS1[0](i_x, i_y);
                }
                else if(i_extra == 1)
                {
                    // intensity weighted LOS magnetic field
                    if(nr_spectral_bins >= 2 && matrixS1[2](i_x, i_y) > 0)
                        array_S[i] = matrixS1[0](i_x, i_y) / matrixS1[2](i_x, i_y);
                    else
                        array_S[i] = 0;
                }
                else if(i_extra == 2)
                {
                    // intensity weighted total magnetic field
                    if(nr_spectral_bins >= 2 && matrixS1[2](i_x, i_y) > 0)
                        array_S[i] = matrixS1[1](i_x, i_y) / matrixS1[2](i_x, i_y);
                    else
                        array_S[i] = 0;
                }
                else if(i_extra == 3)
                {
                    // density weighted LOS magnetic field
                    if(nr_spectral_bins >= 5 && matrixS1[5](i_x, i_y) > 0)
                        array_S[i] = matrixS1[3](i_x, i_y) / matrixS1[5](i_x, i_y);
                    else
                        array_S[i] = 0;
                }
                else if(i_extra == 4)
                {
                    // density weighted magnetic field
                    if(nr_spectral_bins >= 5 && matrixS1[5](i_x, i_y) > 0)
                        array_S[i] = matrixS1[4](i_x, i_y) / matrixS1[5](i_x, i_y);
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

    newTable->addKey(
        "GAS_SPECIES", gas->getGasSpeciesName(i_species), "name of the observed gas species");
    newTable->addKey("TRANS", i_trans + 1, "transition index number (see leiden database)");
    newTable->addKey("LEVEL_UPPER",
                        gas->getUpperEnergyLevel(i_species, i_trans),
                        "upper energy level index number (see leiden database)");
    newTable->addKey("LEVEL_LOWER",
                        gas->getLowerEnergyLevel(i_species, i_trans),
                        "lower energy level index number (see leiden database)");
    newTable->addKey(
        "FREQ", gas->getTransitionFrequency(i_species, i_trans), "frequency of the simulated transition");
    newTable->addKey("CHANNELS", nr_spectral_bins, "number of velocity channels");
    newTable->addKey("MINVEL", min_velocity, "velocity of the velocity channels (minvel to maxvel)");
    newTable->addKey("MAXVEL", max_velocity1, "velocity of the velocity channels (minvel to maxvel)");
    newTable->addKey("ZEEMAN",
                        gas->isTransZeemanSplit(i_species, i_trans),
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

bool CDetector::writeVelChannelHealMapsTiny(CGasMixture * gas, uint i_species, uint i_line)
{
    cout << CLR_LINE << flush;
    cout << " -> Writing channel map(s) ...  \r" << flush;
    bool hasZeeman = gas->isTransZeemanSplit(i_species, i_trans);
    
    llist indices;
    llist bx, by;
    valarray<long> array_index;
    long nelements = 0;

    for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
    {
        long naxis = 1;
        long naxes[1] = { 0 };

        // auto_ptr<CCfits::FITS> pFits(0);
        unique_ptr<CCfits::FITS> pFits;

        try
        {
            char str_tmp[1024];
            char str_end[1024];

#ifdef WINDOWS
            strcpy_s(str_tmp, "min_channel_maps_species_%04d_line_%04d_vel_%04d");
            sprintf_s(str_end, str_tmp, i_species + 1, i_line + 1, i_spectral + 1);
#else
            strcpy(str_tmp, "min_channel_maps_species_%04d_line_%04d_vel_%04d");
            sprintf(str_end, str_tmp, i_species + 1, i_line + 1, i_spectral + 1);
#endif

            string path_out = path + str_end + FITS_COMPRESS_EXT;
            remove(path_out.c_str());
            pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
        }
        catch(CCfits::FITS::CantCreate)
        {
            return false;
        }

        long tmp_nelements = bins_x * bins_y;
        if(max_cells != tmp_nelements)
        {
            cout << CLR_LINE;
            cout << ERROR_LINE << "Max cells are not equal to bins x bins!" << endl;
            return false;
        }

        long index = 0;
        
        nelements=0;
        
        indices.clear();
        bx.clear();
        by.clear();

        for(uint i_y = 0; i_y < bins_y; i_y++)
        {
            for(uint i_x = 0; i_x < bins_x; i_x++)
            {
                if(matrixI[i_spectral](i_x, i_y)>1.0e-100)
                {
                    indices.push_back(index);  
                    bx.push_back(i_x);
                    by.push_back(i_y);
                }

                index++;
            }
        }
        
        nelements=indices.size();        
        
        if(nelements==0)
        {
            cout << CLR_LINE;
            cout << ERROR_LINE << "Healpix map contains no data at all!" << endl;
            return false;
        }        
        
        cout << CLR_LINE;
        cout << " -> Writing channel map " << i_spectral << " of " << nr_spectral_bins << "   \r" << flush;

        // Init columns
        string newName("DATA_EXTENSION");
        uint nr_of_quantities = 2;
        
        if(hasZeeman)
            nr_of_quantities = 5;
        
        vector<string> colName(nr_of_quantities, "");
        vector<string> colForm(nr_of_quantities, "");
        vector<string> colUnit(nr_of_quantities, "");

        if(hasZeeman)
        {
            colName[0] = "HEAL_INDEX";
            colName[1] = "I_STOKES";
            colName[2] = "Q_STOKES";
            colName[3] = "U_STOKES";
            colName[4] = "V_STOKES";

            colForm[0] = "J";
            colForm[1] = "D";
            colForm[2] = "D";
            colForm[3] = "D";
            colForm[4] = "D";

            colUnit[0] = "";
            colUnit[1] = "Jy/px";
            colUnit[2] = "Jy/px";
            colUnit[3] = "Jy/px";
            colUnit[4] = "Jy/px";
        }
        else
        {
            colName[0] = "HEAL_INDEX";
            colName[1] = "I_STOKES";

            colForm[0] = "J";
            colForm[1] = "D";

            colUnit[0] = "";
            colUnit[1] = "Jy/px";
        }
        
        CCfits::Table * newTable = pFits->addTable(newName, nelements, colName, colForm, colUnit);
        
        array_index.resize(nelements,0);
        valarray<double> array_I(nelements);
        valarray<double> array_Q(nelements);
        valarray<double> array_U(nelements);
        valarray<double> array_V(nelements);

        for(uint i = 0; i < nelements; i++)
        {
            long i_x = bx[i];
            long i_y = by[i];
            
            array_index[i]=indices[i];
            array_I[i] = matrixI[i_spectral](i_x, i_y);

            if(hasZeeman)
            {
                array_Q[i] = matrixQ[i_spectral](i_x, i_y);
                array_U[i] = matrixU[i_spectral](i_x, i_y);
                array_V[i] = matrixV[i_spectral](i_x, i_y);
            }
        }

        newTable->column(colName[0]).write(array_index, 1);
        newTable->column(colName[1]).write(array_I, 1);
        
        if(hasZeeman)
        {
            newTable->column(colName[2]).write(array_Q, 1);
            newTable->column(colName[3]).write(array_U, 1);
            newTable->column(colName[4]).write(array_V, 1);
        }

        // Healpix parameter
        newTable->addKey("PIXTYPE", "HEALPIX", "Pixel algorigthm");
        newTable->addKey("ORDERING", "RING", "Ordering scheme");
        newTable->addKey("INDXSCHM", "IMPLICIT", "Indexing scheme");
        newTable->addKey("NSIDE", uint(sqrt(bins_x / 12.0)), "Resolution Parameter");
        newTable->addKey("FIRSTPIX", 0, "First pixel (0 based)");
        newTable->addKey("LASTPIX", max_cells - 1, "Last pixel (0 based)");
        newTable->addKey("CROTA2", 0, "Rotation Angle (Degrees)");

        newTable->addKey(
            "GAS_SPECIES", gas->getGasSpeciesName(i_species), "name of the observed gas species");
        newTable->addKey("TRANS", i_trans + 1, "transition index number (see leiden database)");
        newTable->addKey("LEVEL_UPPER",
                            gas->getUpperEnergyLevel(i_species, i_trans),
                            "upper energy level index number (see leiden database)");
        newTable->addKey("LEVEL_LOWER",
                            gas->getLowerEnergyLevel(i_species, i_trans),
                            "lower energy level index number (see leiden database)");
        newTable->addKey("FREQ",
                            gas->getTransitionFrequency(i_species, i_trans),
                            "frequency of the simulated transition");
        newTable->addKey("VCH", i_spectral, "current velocity channel");
        newTable->addKey("CHANNELS", nr_spectral_bins, "number of velocity channels");
        newTable->addKey("MINVEL", min_velocity, "velocity of the velocity channels (minvel to maxvel)");
        newTable->addKey("MAXVEL", max_velocity1, "velocity of the velocity channels (minvel to maxvel)");
        newTable->addKey("ZEEMAN",
                            gas->isTransZeemanSplit(i_species, i_trans),
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
    long naxes[1] = { 0 };

    // auto_ptr<CCfits::FITS> pFits(0);
    unique_ptr<CCfits::FITS> pFits;

    try
    {
        char str_tmp[1024];
        char str_end[1024];

#ifdef WINDOWS
        strcpy_s(str_tmp, "tiny_channel_maps_species_%04d_line_%04d_extra");
        sprintf_s(str_end, str_tmp, i_species + 1, i_line + 1);
#else
        strcpy(str_tmp, "tiny_channel_maps_species_%04d_line_%04d_extra");
        sprintf(str_end, str_tmp, i_species + 1, i_line + 1);
#endif

        string path_out = path + str_end + FITS_COMPRESS_EXT;
        remove(path_out.c_str());
        pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
    }

    catch(CCfits::FITS::CantCreate)
    {   
        cout << CLR_LINE;
        return false;
    }

    // Init columns
    string newName("DATA_EXTENSION");
    uint nr_of_quantities = 2;

    if(gas->isTransZeemanSplit(i_species, i_trans))
        nr_of_quantities += 4;

    vector<string> colName(nr_of_quantities, "");
    vector<string> colForm(nr_of_quantities, "");
    vector<string> colUnit(nr_of_quantities, "");
    
    colName[0] = "HEAL_INDEX";
    colForm[0] = "J";
    colUnit[0] = "";

    colName[1] = "COLUMN_DENSITY";
    colForm[1] = "D";
    colUnit[1] = "m^-2";

    if(gas->isTransZeemanSplit(i_species, i_trans))
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
    
    
    if(nelements==0)
    {
        cout << CLR_LINE;
        cout << ERROR_LINE << "Healpix map contains no data at all!" << endl;
        return false;
    } 

    CCfits::Table * newTable = pFits->addTable(newName, nelements, colName, colForm, colUnit);
    
    newTable->column(colName[0]).write(array_index, 1);

    // Column density / extra data
    for(uint i_extra = 1; i_extra < nr_of_quantities; i_extra++)
    {
        valarray<double> array_S(nelements);
        uint i = 0;
        
        for(uint i = 0; i < nelements; i++)
        {
            long i_x = bx[i];
            long i_y = by[i];
            
            array_index[i]=indices[i];
            
            if(i_extra == 0+1)
            {
                // column density of the total gas
                if(gas->isTransZeemanSplit(i_species, i_trans))
                    array_S[i] = matrixS1[5](i_x, i_y);
                else
                    array_S[i] = matrixS1[0](i_x, i_y);
            }
            else if(i_extra == 1+1)
            {
                // intensity weighted LOS magnetic field
                if(nr_spectral_bins >= 2 && matrixS1[2](i_x, i_y) > 0)
                    array_S[i] = matrixS1[0](i_x, i_y) / matrixS1[2](i_x, i_y);
                else
                    array_S[i] = 0;
            }
            else if(i_extra == 2+1)
            {
                // intensity weighted total magnetic field
                if(nr_spectral_bins >= 2 && matrixS1[2](i_x, i_y) > 0)
                    array_S[i] = matrixS1[1](i_x, i_y) / matrixS1[2](i_x, i_y);
                else
                    array_S[i] = 0;
            }
            else if(i_extra == 3+1)
            {
                // density weighted LOS magnetic field
                if(nr_spectral_bins >= 5 && matrixS1[5](i_x, i_y) > 0)
                    array_S[i] = matrixS1[3](i_x, i_y) / matrixS1[5](i_x, i_y);
                else
                    array_S[i] = 0;
            }
            else if(i_extra == 4+1)
            {
                // density weighted magnetic field
                if(nr_spectral_bins >= 5 && matrixS1[5](i_x, i_y) > 0)
                    array_S[i] = matrixS1[4](i_x, i_y) / matrixS1[5](i_x, i_y);
                else
                    array_S[i] = 0;
            }
            
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

    newTable->addKey(
        "GAS_SPECIES", gas->getGasSpeciesName(i_species), "name of the observed gas species");
    newTable->addKey("TRANS", i_trans + 1, "transition index number (see leiden database)");
    newTable->addKey("LEVEL_UPPER",
                        gas->getUpperEnergyLevel(i_species, i_trans),
                        "upper energy level index number (see leiden database)");
    newTable->addKey("LEVEL_LOWER",
                        gas->getLowerEnergyLevel(i_species, i_trans),
                        "lower energy level index number (see leiden database)");
    newTable->addKey(
        "FREQ", gas->getTransitionFrequency(i_species, i_trans), "frequency of the simulated transition");
    newTable->addKey("CHANNELS", nr_spectral_bins, "number of velocity channels");
    newTable->addKey("MINVEL", min_velocity, "velocity of the velocity channels (minvel to maxvel)");
    newTable->addKey("MAXVEL", max_velocity1, "velocity of the velocity channels (minvel to maxvel)");
    newTable->addKey("ZEEMAN",
                        gas->isTransZeemanSplit(i_species, i_trans),
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

    cout << CLR_LINE;
    return true;
}

bool CDetector::writeOPIATEIntVelChannelHealMaps(COpiateDataBase * op, uint det_id)
{
    cout << CLR_LINE << flush;
    cout << " -> Writing integrated velocity channel map(s) ...  \r" << flush;

    long naxis = 1;
    long naxes[1] = { 0 };

    // auto_ptr<CCfits::FITS> pFits(0);
    unique_ptr<CCfits::FITS> pFits;

    try
    {
        char str_tmp[1024];
        char str_end[1024];

#ifdef WINDOWS
        strcpy_s(str_tmp, "int_channel_map_species_%04d_line_%04d");
        sprintf_s(str_end, str_tmp, i_species + 1, i_line + 1);
#else
        strcpy(str_tmp, "opiate_int_map__%04d");
        sprintf(str_end, str_tmp, det_id+1);
#endif

        string path_out = path + str_end + FITS_COMPRESS_EXT;
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
        cout << WARNING_LINE << "Max cells are not equal to bins x bins!" << endl;
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

    CCfits::Table * newTable = pFits->addTable(newName, nelements, colName, colForm, colUnit);

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
            for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
            {
                array_I[i] +=
                    matrixI[i_spectral](i_x, i_y) * (channel_width / nr_spectral_bins);// * 1e-3;
                array_Q[i] +=
                    matrixQ[i_spectral](i_x, i_y) * (channel_width / nr_spectral_bins);// * 1e-3;
                array_U[i] +=
                    matrixU[i_spectral](i_x, i_y) * (channel_width / nr_spectral_bins);// * 1e-3;
                array_V[i] +=
                    matrixV[i_spectral](i_x, i_y) * (channel_width / nr_spectral_bins);// * 1e-3;
            }
            array_T[i] = matrixT[int(nr_spectral_bins / 2.0)](i_x, i_y);

            if(nr_spectral_bins >= 2 && matrixS1[2](i_x, i_y) == 0)
            {
                // LOS magnetic field strength
                array_S[i] = 0.0;
            }
            else
                array_S[i] = matrixS1[0](i_x, i_y) / matrixS1[2](i_x, i_y);
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

    newTable->addKey("GAS_SPECIES", op->getCurrentName(), "name of the observed gas species");
    newTable->addKey("TRANS", i_trans + 1, "transition index number (see leiden database)");
    newTable->addKey("FREQ", op->getCurrentFrequency(), "frequency of the simulated transition");
    newTable->addKey("CHANNELS", nr_spectral_bins, "number of velocity channels");
    newTable->addKey("MINVEL", min_velocity, "velocity of the velocity channels (minvel to maxvel)");
    newTable->addKey("MAXVEL", max_velocity1, "velocity of the velocity channels (minvel to maxvel)");
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

bool CDetector::writeIntVelChannelHealMaps(CGasMixture * gas, uint i_species, uint i_line)
{
    cout << CLR_LINE << flush;
    cout << " -> Writing integrated velocity channel map(s) ...  \r" << flush;

    long naxis = 1;
    long naxes[1] = { 0 };
    bool hasZeeman = gas->isTransZeemanSplit(i_species, i_trans);
    
    uint nr_of_quantities = 3;
    
    if(hasZeeman)
        nr_of_quantities = 6;
            

    // auto_ptr<CCfits::FITS> pFits(0);
    unique_ptr<CCfits::FITS> pFits;

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

        string path_out = path + str_end + FITS_COMPRESS_EXT;
        remove(path_out.c_str());
        pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
    }
    catch(CCfits::FITS::CantCreate)
    {
        return false;
    }

    long nelements = bins_x * bins_y;
    if(max_cells != nelements)
    {
        cout << WARNING_LINE << "Max cells are not equal to bins x bins!" << endl;
        return false;
    }

    string newName("HEALPIX_EXTENSION");
    vector<string> colName(nr_of_quantities, "");
    vector<string> colForm(nr_of_quantities, "");
    vector<string> colUnit(nr_of_quantities, "");

    if(hasZeeman)
    {
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
    }
    else
    {
        colName[0] = "I_STOKES";
        colName[1] = "TAU";
        colName[2] = "COLUMN_DENSITY";

        colForm[0] = "D";
        colForm[1] = "D";
        colForm[2] = "D";

        colUnit[0] = "Jy/px";
        colUnit[1] = "";
        colUnit[2] = "m^-2";
    }

    CCfits::Table * newTable = pFits->addTable(newName, nelements, colName, colForm, colUnit);

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
            for(uint i_spectral = 0; i_spectral < nr_spectral_bins; i_spectral++)
            {
                array_I[i] += matrixI[i_spectral](i_x, i_y) * (channel_width / nr_spectral_bins);// * 1e-3;
                
                if(hasZeeman)
                {
                    array_Q[i] += matrixQ[i_spectral](i_x, i_y) * (channel_width / nr_spectral_bins);// * 1e-3;
                    array_U[i] += matrixU[i_spectral](i_x, i_y) * (channel_width / nr_spectral_bins);// * 1e-3;
                    array_V[i] += matrixV[i_spectral](i_x, i_y) * (channel_width / nr_spectral_bins);// * 1e-3;
                }
            }
            
            array_T[i] = matrixT[int(nr_spectral_bins / 2.0)](i_x, i_y);

            if(nr_spectral_bins >= 2 && matrixS1[2](i_x, i_y) == 0)
            {
                // LOS magnetic field strength
                array_S[i] = 0.0;
            }
            else
                array_S[i] = matrixS1[0](i_x, i_y) / matrixS1[2](i_x, i_y);
            
            i++;
        }

    if(hasZeeman)
    {
        newTable->column(colName[0]).write(array_I, 1);
        newTable->column(colName[1]).write(array_Q, 1);
        newTable->column(colName[2]).write(array_U, 1);
        newTable->column(colName[3]).write(array_V, 1);
        newTable->column(colName[4]).write(array_T, 1);
        newTable->column(colName[5]).write(array_S, 1);
    }
    else
    {
        newTable->column(colName[0]).write(array_I, 1);
        newTable->column(colName[1]).write(array_T, 1);
        newTable->column(colName[2]).write(array_S, 1);
    }

    // Healpix parameter
    newTable->addKey("PIXTYPE", "HEALPIX", "Pixel algorigthm");
    newTable->addKey("ORDERING", "RING", "Ordering scheme");
    newTable->addKey("INDXSCHM", "IMPLICIT", "Indexing scheme");
    newTable->addKey("NSIDE", uint(sqrt(bins_x / 12.0)), "Resolution Parameter");
    newTable->addKey("FIRSTPIX", 0, "First pixel (0 based)");
    newTable->addKey("LASTPIX", max_cells - 1, "Last pixel (0 based)");
    newTable->addKey("CROTA2", 0, "Rotation Angle (Degrees)");

    newTable->addKey(
        "GAS_SPECIES", gas->getGasSpeciesName(i_species), "name of the observed gas species");
    newTable->addKey("TRANS", i_trans + 1, "transition index number (see leiden database)");
    newTable->addKey("LEVEL_UPPER",
                        gas->getUpperEnergyLevel(i_species, i_trans),
                        "upper energy level index number (see leiden database)");
    newTable->addKey("LEVEL_LOWER",
                        gas->getLowerEnergyLevel(i_species, i_trans),
                        "lower energy level index number (see leiden database)");
    newTable->addKey(
        "FREQ", gas->getTransitionFrequency(i_species, i_trans), "frequency of the simulated transition");
    newTable->addKey("CHANNELS", nr_spectral_bins, "number of velocity channels");
    newTable->addKey("MINVEL", min_velocity, "velocity of the velocity channels (minvel to maxvel)");
    newTable->addKey("MAXVEL", max_velocity1, "velocity of the velocity channels (minvel to maxvel)");
    newTable->addKey("ZEEMAN",
                        gas->isTransZeemanSplit(i_species, i_trans),
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

void CDetector::calcCoordinateParameters(double sidelength,
                              double bins,
                              double map_shift,
                              double distance,
                              double & bin_width,
                              double & first_pix_val,
                              double & deg_per_pix,
                              double & first_pix_val_deg)
{
    bin_width = sidelength / bins;

    first_pix_val = -sidelength / 2.0 + map_shift;
    first_pix_val += (bin_width / 2.0);

    if(distance>0)
    {
        deg_per_pix = bin_width / con_AU;
        deg_per_pix /= (distance / con_pc);
        deg_per_pix /= 3600.0;

        first_pix_val_deg = -sidelength / (2.0 * con_AU) + map_shift / con_AU;
        first_pix_val_deg /= (distance / con_pc);
        first_pix_val_deg /= 3600.0;
        first_pix_val_deg += (deg_per_pix / 2.0);
    }
}

void CDetector::calcVelocityChannels(uint _nr_spectral_bins, double _min_velocity, double _max_velocity)
{
    nr_spectral_bins = _nr_spectral_bins;
    min_velocity = _min_velocity;
    max_velocity1 = _max_velocity;

    if(nr_spectral_bins > 1)
        channel_width = (_max_velocity-_min_velocity) / (nr_spectral_bins - 1);
    else
        channel_width = (_max_velocity-_min_velocity);

    velocity_channel.resize(nr_spectral_bins);
    
    if(nr_spectral_bins > 1)
    {
        for(uint i = 0; i < nr_spectral_bins; i++)
        {
            velocity_channel[i] = _min_velocity + i* channel_width;
        }
        
        velocity_channel[nr_spectral_bins-1]=_max_velocity;
    }           
    else if(nr_spectral_bins == 1)
    {
        velocity_channel[0] = 0.5*channel_width;
    }
}

/*void CDetector::calcVelocityChannels(uint _nr_spectral_bins, double _min_velocity, double _max_velocity)
{
    nr_spectral_bins = _nr_spectral_bins;
    max_velocity = _max_velocity;

    if(nr_spectral_bins > 1)
        channel_width = (2.0 * max_velocity) / (nr_spectral_bins - 1);
    else
        channel_width = (2.0 * max_velocity);

    velocity_channel.resize(nr_spectral_bins);
    if(nr_spectral_bins > 1)
        for(uint i = 0; i < nr_spectral_bins; i++)
            velocity_channel[i] =
                2 * (float)i / ((float)nr_spectral_bins - 1) * max_velocity - max_velocity;
    else if(nr_spectral_bins == 1)
        velocity_channel[0] = 0;
}
*/

string CDetector::getPath()
{
    return path;
}

double CDetector::getLamMin()
{
    return lam_min;
}

double CDetector::getLamMax()
{
    return lam_max;
}

uint CDetector::getDetectorWavelengthID(double wavelength)
{
    dlist::iterator it = find(wavelength_list_det.begin(), wavelength_list_det.end(), wavelength);
    if(it != wavelength_list_det.end())
        return it - wavelength_list_det.begin();
    else
        return MAX_UINT;
}

double CDetector::getWavelength(uint i_wave)
{
    return wavelength_list_det[i_wave];
}

uint CDetector::getNrOfSpectralBins()
{
    return nr_spectral_bins;
}

double CDetector::getChannelWidth()
{
    return channel_width;
}

double CDetector::getVelocityChannel(uint i_spectral)
{
    return velocity_channel[i_spectral];
}

string CDetector::getDetectorGridDescription()
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
            return "Polar background grid";
            break;

        case DET_SLICE:
            return "Cartesian grid as a slice through the grid (y-direction along "
                    "LOS)";
            break;

        default:
            return "Photon packages launched from emission sources (Monte-Carlo)";
            break;
    }
}

void CDetector::init()
{
    //Stokes components
    matrixI = 0;
    matrixQ = 0;
    matrixU = 0;
    matrixV = 0;

    //special parameters
    matrixT = 0;
    matrixS1 = 0;
    matrixS2 = 0;
    matrixS3 = 0;
    matrixS4 = 0;

    w1_I = 0;
    w2_I = 0;
    w3_I = 0;
    w4_I = 0;
    w1_Q = 0;
    w2_Q = 0;
    w3_Q = 0;
    w4_Q = 0;
    w1_U = 0;
    w2_U = 0;
    w3_U = 0;
    w4_U = 0;
    w1_V = 0;
    w2_V = 0;
    w3_V = 0;
    w4_V = 0;
    w1_PI = 0;
    w2_PI = 0;
    w3_PI = 0;
    w4_PI = 0;
    N_photon = 0;
    sedI = 0;
    sedQ = 0;
    sedU = 0;
    sedV = 0;
    sedT = 0;
    sedS = 0;
}

string CDetector::getAlignmentDescription()
{
    string alignment_string = "";

    if((alignment & ALIG_INTERNAL) == ALIG_INTERNAL)
        alignment_string += "INTERNAL";
    if((alignment & ALIG_PA) == ALIG_PA)
    {
        if(alignment_string != "")
            alignment_string += ", ";
        alignment_string += "PA";
    }
    if((alignment & ALIG_NONPA) == ALIG_NONPA)
    {
        if(alignment_string != "")
            alignment_string += ", ";
        alignment_string += "NONPA";
    }
    if((alignment & ALIG_IDG) == ALIG_IDG)
    {
        if(alignment_string != "")
            alignment_string += ", ";
        alignment_string += "IDG";
    }
    if((alignment & ALIG_RAT) == ALIG_RAT)
    {
        if(alignment_string != "")
            alignment_string += ", ";
        alignment_string += "RAT";
    }
    if((alignment & ALIG_GOLD) == ALIG_GOLD)
    {
        if(alignment_string != "")
            alignment_string += ", ";
        alignment_string += "GOLD";
    }
    if((alignment & ALIG_KRAT) == ALIG_KRAT)
    {
        if(alignment_string != "")
            alignment_string += ", ";
        alignment_string += "KRAT";
    }
    return alignment_string;
}
