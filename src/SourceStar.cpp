/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "SourceStar.hpp"
#include "CommandParser.hpp"

bool CSourceStar::initSource(uint id, uint max, bool use_energy_density)
{
    // Initial output
    cout << CLR_LINE << flush;
    cout << "-> Initiating source star         \r" << flush;

    // Init variables
    dlist star_emi;
    double tmp_luminosity, diff_luminosity, max_flux = 0;
    ullong kill_counter = 0;

    if(use_energy_density && !is_ext)
    {
        // For using energy density, only the photon number is required
        cout << "- Source (" << id + 1 << " of " << max << ") STAR: " << float(L / L_sun)
             << " [L_sun], photons FOR EACH wavelength: " << float(nr_of_photons) << "      " << endl;
    }
    else
    {
        for(uint w = 0; w < getNrOfWavelength(); w++)
        {
            double pl = CMathFunctions::planck(wavelength_list[w], T); //[W m^-2 m^-1 sr^-1]
            double sp_energy;

            if(is_ext)
                sp_energy = sp_ext.getValue(wavelength_list[w],LINEAR);
            else
                sp_energy =
                    4.0 * PI * PI * (R * R_sun) * (R * R_sun) * pl; //[W m^-1] energy per second an wavelength

            star_emi.push_back(sp_energy);
        }

        tmp_luminosity = CMathFunctions::integ(wavelength_list, star_emi, 0, getNrOfWavelength() - 1);
        L = tmp_luminosity;

        cout << "- Source (" << id + 1 << " of " << max << ") STAR: " << float(L / L_sun)
             << " [L_sun], photons: " << float(nr_of_photons) << endl;

        for(uint w = 0; w < getNrOfWavelength(); w++)
        {
            if(wavelength_list[w] * star_emi[w] > max_flux)
                max_flux = wavelength_list[w] * star_emi[w];
        }

        max_flux *= ACC_SELECT_LEVEL;

        for(uint w = 0; w < getNrOfWavelength(); w++)
            if(wavelength_list[w] * star_emi[w] < max_flux)
            {
                kill_counter++;
                star_emi[w] = 0;
            }

        diff_luminosity = CMathFunctions::integ(wavelength_list, star_emi, 0, getNrOfWavelength() - 1);
        diff_luminosity -= tmp_luminosity;

        cout << "    wavelengths: " << getNrOfWavelength() - kill_counter << " of " << getNrOfWavelength()
             << ", neglected energy: " << float(100.0 * diff_luminosity / tmp_luminosity) << " [%]" << endl;

        double fr;
        lam_pf.resize(getNrOfWavelength());
        for(uint l = 0; l < getNrOfWavelength(); l++)
        {
            fr = CMathFunctions::integ(wavelength_list, star_emi, 0, l) / tmp_luminosity;
            lam_pf.setValue(l, fr, double(l));
        }
    }

    return true;
}

bool CSourceStar::setParameterFromFile(parameters & param, uint p)
{
    dlist values = param.getPointSources();
    string filename = param.getPointSourceString(p / NR_OF_POINT_SOURCES);

    ifstream reader(filename.c_str());
    int line_counter = 0;
    string line;
    CCommandParser ps;

    double w_min = 1e300;
    double w_max = 0;

    is_ext = true;

    pos = Vector3D(values[p], values[p + 1], values[p + 2]);
    R = values[p + 3];
    T = values[p + 4];
    
    r_sub = values[p + 5];

    nr_of_photons = ullong(values[p + NR_OF_POINT_SOURCES - 1]);
    cout << CLR_LINE << flush;
    cout << "-> Loading spectrum for source star...           \r" << flush;

    if(reader.fail())
    {
        cout << ERROR_LINE << "Cannot open spectrum file: \n" << filename << "  \n" << endl;
        return false;
    }

    while(getline(reader, line))
    {
        ps.formatLine(line);

        if(line.size() == 0)
            continue;

        dlist value = ps.parseValues(line);

        if(value.size() != 4 && value.size() != 2)
        {
            cout << ERROR_LINE << "In spectrum file:\n" << filename << endl;
            cout << "Wrong amount of values in line " << line_counter + 1 << "!" << endl;
            return false;
        }

        line_counter++;
        sp_ext.setDynValue(value[0], value[1]);

        if(value.size() == 4)
        {
            sp_ext_q.setDynValue(value[0], value[2]);
            sp_ext_u.setDynValue(value[0], value[3]);
        }
        else
        {
            sp_ext_q.setDynValue(value[0], 0);
            sp_ext_u.setDynValue(value[0], 0);
        }

        if(w_min > value[0])
            w_min = value[0];

        if(w_max < value[0])
            w_max = value[0];
    }

    sp_ext.createDynSpline();
    sp_ext_q.createDynSpline();
    sp_ext_u.createDynSpline();
    reader.close();

    return true;
}

void CSourceStar::createNextRay(photon_package * pp, CRandomGenerator * rand_gen)
{
    // Init variables
    StokesVector tmp_stokes_vector;
    double energy;
    uint wID;

    #if BENCHMARK == TRUST
        // cos(th) should be between -1 and -0.5, thus only a quarter of the sphere
        pp->setRandomDirectionTRUST(rand_gen->getRND(), rand_gen->getRND());
    #else
        // send more photons close to z=0/theta=0, i.e. midplane of a disk
        // theta is drawn from cos(random)^exponent
        // exponent = 1 is isotropic
        // exponent = 3 is ususally sufficient for protoplanetary or debris disks
        // exponent must be an odd int
        uint exponentThetaBias = 1;

        pp->setRandomDirection(rand_gen->getRND(), rand_gen->getRND(), exponentThetaBias);
    #endif


    if(pp->getDustWavelengthID() != MAX_UINT)
    {
        wID = pp->getDustWavelengthID();
        if(is_ext)
        {
            energy = abs(sp_ext.getValue(wavelength_list[wID])) / nr_of_photons;
            double tmp_q = sp_ext_q.getValue(wavelength_list[wID]);
            double tmp_u = sp_ext_u.getValue(wavelength_list[wID]);
            tmp_stokes_vector = energy * StokesVector(1.0, tmp_q, tmp_u, 0);
        }
        else
        {
            double pl = CMathFunctions::planck(wavelength_list[wID], T);
            energy = PIx4 * PI * (R * R_sun) * (R * R_sun) * pl / nr_of_photons;
            tmp_stokes_vector = energy * StokesVector(1.0, q, u, 0);
        }
    }
    else
    {
        energy = L / nr_of_photons;
        wID = lam_pf.getXIndex(rand_gen->getRND());

        if(is_ext)
        {
            double tmp_q = sp_ext_q.getValue(wavelength_list[wID]);
            double tmp_u = sp_ext_u.getValue(wavelength_list[wID]);
            tmp_stokes_vector = energy * StokesVector(1.0, tmp_q, tmp_u, 0);
        }
        else
            tmp_stokes_vector = energy * StokesVector(1.0, q, u, 0);

        // Mol3D uses the upper value of the wavelength interval,
        // used for the selection of the emitting wavelengths from source!
        pp->setWavelength(wavelength_list[wID + 1], wID + 1);
    }

    #if BENCHMARK == TRUST
        // the photon direction was restricted to a quarter of the sphere
        // thus the energy of each photon must be reduced by the same factor
        tmp_stokes_vector /= 4.0;
    #else
        // if the exponent of cos(theta) is not 1, i.e. not isotropic
        // the photon energy must be rescaled accordingly
        if(exponentThetaBias != 1)
        {
            double theta = pp->getDirection().getSphericalCoord().Theta();
            tmp_stokes_vector *= exponentThetaBias * pow(abs(cos(theta)),(double) 1.-1./exponentThetaBias);
        }
    #endif
    
    pp->setPosition(pos);
    
    if(r_sub>0)
    {
        Vector3D ndir = pp->getDirection();
        ndir.normalize();
        
        Vector3D rel_pos = pos + r_sub*ndir;
        pp->setPosition(rel_pos);        
    }   
    
    pp->setStokesVector(tmp_stokes_vector);
    pp->initCoordSystem();
}

void CSourceStar::createDirectRay(photon_package * pp, CRandomGenerator * rand_gen, Vector3D dir_obs)
{
    // Init variables
    StokesVector tmp_stokes_vector;
    double energy;
    uint wID;

    if(pp->getDustWavelengthID() != MAX_UINT)
    {
        wID = pp->getDustWavelengthID();
        if(is_ext)
        {
            energy = sp_ext.getValue(wavelength_list[wID]) / PIx4;
            double tmp_q = sp_ext_q.getValue(wavelength_list[wID]);
            double tmp_u = sp_ext_u.getValue(wavelength_list[wID]);
            tmp_stokes_vector = energy * StokesVector(1.0, tmp_q, tmp_u, 0);
        }
        else
        {
            double pl = CMathFunctions::planck(wavelength_list[wID], T);
            energy = PI * (R * R_sun) * (R * R_sun) * pl;
            tmp_stokes_vector = energy * StokesVector(1.0, q, u, 0);
        }
    }
    else
    {
        energy = L / PIx4;
        wID = lam_pf.getXIndex(rand_gen->getRND());

        if(is_ext)
        {
            double tmp_q = sp_ext_q.getValue(wavelength_list[wID]);
            double tmp_u = sp_ext_u.getValue(wavelength_list[wID]);
            tmp_stokes_vector = energy * StokesVector(1.0, tmp_q, tmp_u, 0);
        }
        else
            tmp_stokes_vector = energy * StokesVector(1.0, q, u, 0);

        // Mol3D uses the upper value of the wavelength interval,
        // used for the selection of the emitting wavelengths from source!
        pp->setWavelength(wavelength_list[wID + 1], wID + 1);
    }

    // Set direction of the photon package to the observer
    if(dir_obs.length() > 0)
    {
        pp->setDirection(dir_obs);
        pp->initCoordSystem();
    }
    
    if(r_sub>0)
    {
        Vector3D ndir = dir_obs.normalized();
        Vector3D rel_pos = pos + r_sub*ndir;
        pp->setPosition(rel_pos);        
    }
    else
        pp->setPosition(pos);
    
    pp->setStokesVector(tmp_stokes_vector);
}

void CSourceStar::setParameter(parameters & param, uint p)
{
    dlist values = param.getPointSources();

    pos = Vector3D(values[p], values[p + 1], values[p + 2]);
    R = values[p + 3];
    T = values[p + 4];
    
    r_sub = values[p + 5];

    q = values[p + 6];
    u = values[p + 7];
    
    nr_of_photons = (ullong)values[p + NR_OF_POINT_SOURCES - 1];

    L = PIx4 * con_sigma * (R * R_sun) * (R * R_sun) * T * T * T * T;
    
    int sub_status = param.getSubStatus();
    
    if((sub_status & SUB_RADIUS) == SUB_RADIUS)
    {
        if(r_sub<=0)
        {
            double T_sub = dust->getMaxSubTemperature();
            
            r_sub = R_SUB * sqrt(L / L_sun) * pow(T_sub / 1500, -2.8); 
        }
    }
}