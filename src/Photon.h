#include "Cell.h"
#include "MathFunctions.h"
#include "Stokes.h"
#include "Typedefs.h"
#include "Vector.h"

#ifndef PHOTON
#define PHOTON

class photon_package
{
  public:
    photon_package(uint nr_bins = 1)
    {
        // Photon package for dust simulations

        // Init variables
        cell_pos = 0;
        tmp_path = 0;
        i_spectral = 0;
        velocity = 0;
        trans_frequency = 0;
        dirID = MAX_UINT;

        // Set total number of wavelength saved in photon package (usually 1)
        nr_of_spectral_bins = nr_bins;

        // Init pointer arrays
        wavelength = new double[nr_of_spectral_bins];
        wID = new uint[nr_of_spectral_bins];
        multi_stokes = new StokesVector[nr_of_spectral_bins];

        // Set initial values
        for(uint i_spectral = 0; i_spectral < nr_of_spectral_bins; i_spectral++)
        {
            wavelength[i_spectral] = 0;

            // Dust wavelength index has to be MAX_UINT if not set in routine
            wID[i_spectral] = MAX_UINT;
        }
    }

    photon_package(double trans_freq, uint _wID, uint nr_bins = 1)
    {
        // Photon package for spectral line simulations

        // Init variables
        cell_pos = 0;
        tmp_path = 0;
        i_spectral = 0;
        dirID = MAX_UINT;

        // Set number of velocity channels saved in photon package
        nr_of_spectral_bins = nr_bins;

        // Set rest frequency of transition (used for dust properties)
        trans_frequency = trans_freq;

        // Init pointer arrays
        velocity = new double[nr_of_spectral_bins];
        multi_stokes = new StokesVector[nr_of_spectral_bins];

        // Set initial values
        for(uint i_spectral = 0; i_spectral < nr_of_spectral_bins; i_spectral++)
        {
            velocity[i_spectral] = 0;
        }

        // Set wavelength for dust contribution
        wavelength = new double[1];
        wavelength[0] = con_c / getTransFrequency();

        // Set wavelength index to get the right dust properties
        wID = new uint[1];
        wID[0] = _wID;
    }

    ~photon_package()
    {
        if(multi_stokes != 0)
            delete[] multi_stokes;

        if(wID != 0)
            delete[] wID;

        if(wavelength != 0)
            delete[] wavelength;

        if(velocity != 0)
            delete[] velocity;
    }

    // photon_package & operator=(photon_package * rhs)
    // {
    //     nr_of_spectral_bins = rhs->getNrOfSpectralBins();

    //     if(multi_stokes != 0)
    //         delete[] multi_stokes;

    //     if(wID != 0)
    //         delete[] wID;

    //     if(wavelength != 0)
    //         delete[] wavelength;

    //     if(velocity != 0)
    //         delete[] velocity;

    //     trans_frequency = rhs->getTransFrequency();
    //     if(trans_frequency != 0)
    //     {
    //         // Init pointer arrays
    //         velocity = new double[nr_of_spectral_bins];
    //         multi_stokes = new StokesVector[nr_of_spectral_bins];

    //         // Set initial values
    //         multi_stokes = rhs->getMultiStokes();
    //         velocity = rhs->getVelocities();

    //         // Set wavelength for dust contribution
    //         wavelength = new double[1];
    //         wavelength[0] = con_c / rhs->getTransFrequency();

    //         // Set wavelength index to get the right dust properties
    //         wID = new uint[1];
    //         wID[0] = rhs->getDustWavelengthID();
    //     }
    //     else
    //     {
    //         // Init pointer arrays
    //         wavelength = new double[nr_of_spectral_bins];
    //         wID = new uint[nr_of_spectral_bins];
    //         multi_stokes = new StokesVector[nr_of_spectral_bins];

    //         // Set initial values
    //         for(uint i_spectral = 0; i_spectral < nr_of_spectral_bins; i_spectral++)
    //         {
    //             wavelength[i_spectral] = 0;

    //             // Dust wavelength index has to be MAX_UINT if not set in routine
    //             wID[i_spectral] = MAX_UINT;
    //         }
    //     }

    //     delete rhs;
    //     rhs = 0;

    //     return *this;
    // }

    void setStokesVector(StokesVector st, uint i = MAX_UINT)
    {
        if(i == MAX_UINT)
            multi_stokes[i_spectral] = st;
        else
            multi_stokes[i] = st;
    }

    const StokesVector & getStokesVector(uint i = MAX_UINT) const
    {
        if(i == MAX_UINT)
            return multi_stokes[i_spectral];
        else
            return multi_stokes[i];
    }

    StokesVector * getStokesVector(uint i = MAX_UINT)
    {
        if(i == MAX_UINT)
            return &multi_stokes[i_spectral];
        else
            return &multi_stokes[i];
    }

    StokesVector * getMultiStokes()
    {
        return multi_stokes;
    }

    void setSpectralID(uint i)
    {
        // Current index of the photon package (default wavelength or frequency)
        if(i < nr_of_spectral_bins)
            i_spectral = i;
    }

    uint getSpectralID() const
    {
        return i_spectral;
    }

    uint getNrOfSpectralBins()
    {
        return nr_of_spectral_bins;
    }

    uint getDustWavelengthID() const
    {
        if(trans_frequency > 0)
            return wID[0];
        return wID[i_spectral];
    }

    void setWavelength(double wave, uint _wID)
    {
        wavelength[i_spectral] = wave;
        wID[i_spectral] = _wID;
    }

    double getWavelength() const
    {
        if(trans_frequency > 0)
            return wavelength[0];
        return wavelength[i_spectral];
    }

    double getFrequency() const
    {
        if(trans_frequency > 0)
            return getTransFrequency() + getVelocity() / con_c * getTransFrequency();
        return con_c / getWavelength();
    }

    void setVelocity(double vel)
    {
        velocity[i_spectral] = vel;
    }

    double getVelocity() const
    {
        return velocity[i_spectral];
    }

    double * getVelocities() const
    {
        return velocity;
    }

    double getTransFrequency() const
    {
        return trans_frequency;
    }

    double randn(double mu, double sigma)
    {
        double U1, U2, W, mult;
        double X1, X2;
        int call = 0;

        if(call == 1)
        {
            call = !call;
            return (mu + sigma * (double)X2);
        }

        do
        {
            U1 = -1 + getRND() * 2;
            U2 = -1 + getRND() * 2;
            W = pow(U1, 2) + pow(U2, 2);
        } while(W >= 1 || W == 0);

        mult = sqrt((-2 * log(W)) / W);
        X1 = U1 * mult;
        X2 = U2 * mult;

        call = !call;

        double res = mu + sigma * (double)X1;

        if(res < 0)
            return randn(mu, sigma);

        return res;
    }

    void initRandomGenerator(ullong seed)
    {
        rand_gen.setSeed(seed);
    }

    void setD(Matrix2D _mD)
    {
        mD = _mD;
    }

    const Matrix2D & getD() const
    {
        return mD;
    }

    const Vector3D & getPosition() const
    {
        return pos;
    }

    const Vector3D & getDirection() const
    {
        return ez;
    }

    double getRND()
    {
        return rand_gen.getValue();
    }

    double getTmpPathLength() const
    {
        return tmp_path;
    }

    cell_basic * getPositionCell()
    {
        return cell_pos;
    }

    const cell_basic * getPositionCell() const
    {
        return cell_pos;
    }

    void calcRandomDirection()
    {
        ez.rndDir(getRND(), getRND());
    }

    void initCoordSystem()
    {
        double phi = atan3(ez.Y(), -ez.X());
        double theta = acos(ez.Z());
        double cos_phi = cos(phi);
        double sin_phi = sin(phi);
        double cos_theta = ez.Z();
        double sin_theta = sin(theta);

        mD = CMathFunctions::getRotationMatrix(cos_phi, sin_phi, cos_theta, sin_theta);

        ex = mD * Vector3D(1, 0, 0);
        ey = mD * Vector3D(0, 1, 0);
    }

    void updateCoordSystem()
    {
        double phi = atan3(ez.Y(), -ez.X());
        double theta = acos(ez.Z());
        double cos_phi = cos(phi);
        double sin_phi = sin(phi);
        double cos_theta = ez.Z();
        double sin_theta = sin(theta);

        Matrix2D D_help = CMathFunctions::getRotationMatrix(cos_phi, sin_phi, cos_theta, sin_theta);
        mD = mD * D_help;

        ex = mD * Vector3D(1, 0, 0);
        ey = mD * Vector3D(0, 1, 0);
        ez = mD * Vector3D(0, 0, 1);
    }

    void updateCoordSystem(double phi, double theta)
    {
        double cos_phi = cos(phi);
        double sin_phi = sin(phi);
        double cos_theta = cos(theta);
        double sin_theta = sin(theta);

        Matrix2D D_help = CMathFunctions::getRotationMatrix(cos_phi, sin_phi, cos_theta, sin_theta);
        mD = mD * D_help;

        ex = mD * Vector3D(1, 0, 0);
        ey = mD * Vector3D(0, 1, 0);
        ez = mD * Vector3D(0, 0, 1);
    }

    void adjustPosition(Vector3D _pos, double _len)
    {
        tmp_path = _len;
        pos = _pos + tmp_path * ez;

        // Calculate cell by position
        if(tmp_path > 0)
            dirID = MAX_UINT;
    }

    void setPosition(Vector3D val)
    {
        pos = val;
    }

    void addPosition(Vector3D val)
    {
        pos += val;
    }

    void setBackupPosition(Vector3D val)
    {
        backup_pos = val;
    }

    void updateBackupPosition()
    {
        backup_pos = pos;
    }

    void resetPositionToBackupPos()
    {
        pos = backup_pos;
    }

    bool reachedBackupPosition()
    {
        if(ez * (backup_pos - pos) <= 0)
            return true;
        return false;
    }

    bool reachedBackupPosition(Vector3D tmp_pos)
    {
        if(ez * (backup_pos - tmp_pos) <= 0)
            return true;
        return false;
    }

    const Vector3D & getBackupPosition() const
    {
        return backup_pos;
    }

    void setDetectorProjection()
    {
        double tmp_vec = ey * pos;
        pos.setX(ex * pos);
        pos.setY(tmp_vec);
    }

    void setRelativePosition(double tx, double ty, double tz)
    {
        pos = tx * ex + ty * ey + tz * ez;
    }

    void setDirection(Vector3D val)
    {
        ez = val;
    }

    void addDirection(Vector3D val)
    {
        ez += val;
    }

    void multDirection(double val)
    {
        ez *= val;
    }

    const Vector3D & getEX() const
    {
        return ex;
    }

    const Vector3D & getEY() const
    {
        return ey;
    }

    const Vector3D & getEZ() const
    {
        return ez;
    }

    void setEX(Vector3D _e)
    {
        ex = _e;
    }

    void setEY(Vector3D _e)
    {
        ey = _e;
    }

    void setEZ(Vector3D _e)
    {
        ez = _e;
    }

    void setCoordinateSystem(const Vector3D & _ex, const Vector3D & _ey, const Vector3D & _ez)
    {
        ex = _ex;
        ey = _ey;
        ez = _ez;
    }

    void getCoordinateSystem(Vector3D * _ex, Vector3D * _ey, Vector3D * _ez) const
    {
        *_ex = ex;
        *_ey = ey;
        *_ez = ez;
    }

    void setTmpPathLength(double val)
    {
        tmp_path = val;
    }

    void setPositionCell(cell_basic * val)
    {
        cell_pos = val;
    }

    void setDirectionID(uint val)
    {
        dirID = val;
    }

    uint getDirectionID()
    {
        return dirID;
    }

  private:
    CRandomGenerator rand_gen;
    Vector3D pos;
    Vector3D backup_pos;
    Vector3D ex;
    Vector3D ey;
    Vector3D ez;
    Matrix2D mD;
    StokesVector * multi_stokes;

    uint dirID;
    uint i_spectral;
    uint nr_of_spectral_bins;

    double tmp_path;
    double trans_frequency;

    double * velocity;
    double * wavelength;
    uint * wID;

    cell_basic * cell_pos;
};
#endif
