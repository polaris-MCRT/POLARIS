#include <ext/alloc_traits.h>
#include <stdlib.h>
#include <iostream>
#include <string>

#include "Source.h"
#include "CommandParser.h"
#include "Grid.h"
#include "MathFunctions.h"
#include "Parameters.h"
#include "Photon.h"

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
             << " [L_sun], photons FOR EACH wavelength: " << nr_of_photons << "      " << endl;
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
             << " [L_sun], photons: " << nr_of_photons << endl;

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

    nr_of_photons = ullong(values[p + NR_OF_POINT_SOURCES - 1]);
    cout << CLR_LINE << flush;
    cout << "-> Loading spectrum for source star...           \r" << flush;

    if(reader.fail())
    {
        cout << "\nERROR: Cannot open spectrum file: \n" << filename << "  \n" << endl;
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
            cout << "\nERROR: In spectrum file:\n" << filename << endl;
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

void CSourceStar::createNextRay(photon_package * pp, ullong i_phot)
{
    // Init variables
    StokesVector tmp_stokes_vector;
    double energy;
    uint wID;

    pp->initRandomGenerator(i_phot);
    pp->calcRandomDirection();

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
        wID = lam_pf.getXIndex(pp->getRND());

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

    pp->setPosition(pos);
    pp->setStokesVector(tmp_stokes_vector);
    pp->initCoordSystem();
}

void CSourceStar::createDirectRay(photon_package * pp, Vector3D dir_obs)
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
        wID = lam_pf.getXIndex(pp->getRND());

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

    pp->setPosition(pos);
    pp->setStokesVector(tmp_stokes_vector);
}

bool CSourceAGN::initSource(uint id, uint max, bool use_energy_density)
{
    // Initial output
    cout << CLR_LINE << flush;
    cout << "-> Initiating source star         \r" << flush;

    // Init variables
    dlist star_emi;
    double tmp_luminosity, diff_luminosity, max_flux = 0;
    ullong kill_counter = 0;

    if(use_energy_density)
    {
        // For using energy density, only the photon number is required
        cout << "- Source (" << id + 1 << " of " << max << ") STAR: " << float(L / L_sun)
             << " [L_sun], photons per wavelength: " << nr_of_photons << "      " << endl;
    }
    else
    {
        for(uint w = 0; w < getNrOfWavelength(); w++)
        {
            double pl = CMathFunctions::planck(wavelength_list[w], T); //[W m^-2 m^-1 sr^-1]
            double sp_energy;

            if(is_ext)
                sp_energy = sp_ext.getValue(wavelength_list[w]);
            else
                sp_energy =
                    4.0 * PI * PI * (R * R_sun) * (R * R_sun) * pl; //[W m^-1] energy per second an wavelength

            star_emi.push_back(sp_energy);
        }

        tmp_luminosity = CMathFunctions::integ(wavelength_list, star_emi, 0, getNrOfWavelength() - 1);
        L = tmp_luminosity;

        cout << "- Source (" << id + 1 << " of " << max << ") STAR: " << float(L / L_sun)
             << " [L_sun], photons: " << nr_of_photons << endl;

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

bool CSourceAGN::setParameterFromFile(parameters & param, uint p)
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

    nr_of_photons = ullong(values[p + NR_OF_POINT_SOURCES - 1]);
    cout << CLR_LINE << flush;
    cout << "-> Loading spectrum for source star...           \r" << flush;

    if(reader.fail())
    {
        cout << "\nERROR: Cannot open spectrum file: \n" << filename << "  \n" << endl;
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
            cout << "\nERROR: In spectrum file:\n" << filename << endl;
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

void CSourceAGN::createNextRay(photon_package * pp, ullong i_phot)
{
    // Init variables
    StokesVector tmp_stokes_vector;
    double energy;
    uint wID;

    pp->initRandomGenerator(i_phot);
    pp->calcRandomDirection();

    if(pp->getDustWavelengthID() != MAX_UINT)
    {
        wID = pp->getDustWavelengthID();
        if(is_ext)
        {
            energy = sp_ext.getValue(wavelength_list[wID]) / nr_of_photons;
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
        wID = lam_pf.getXIndex(pp->getRND());

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

    pp->setPosition(pos);
    pp->setStokesVector(tmp_stokes_vector);
    pp->initCoordSystem();
}

void CSourceAGN::createDirectRay(photon_package * pp, Vector3D dir_obs)
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
        wID = lam_pf.getXIndex(pp->getRND());

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

    pp->setPosition(pos);
    pp->setStokesVector(tmp_stokes_vector);
}

bool CSourceStarField::initSource(uint id, uint max, bool use_energy_density)
{
    // Initial output
    cout << CLR_LINE << flush;
    cout << "-> Initiating source star field      \r" << flush;

    if(use_energy_density && !is_ext)
        cout << "- Source (" << id + 1 << " of " << max << ") STARFIELD: " << float(L / L_sun)
             << " [L_sun], photons FOR EACH wavelength: " << nr_of_photons << "      " << endl;
    else
    {
        // Init variables
        dlist star_emi;
        double tmp_luminosity, diff_luminosity, max_flux = 0;
        ullong kill_counter = 0;

        for(uint w = 0; w < getNrOfWavelength(); w++)
        {
            double pl = CMathFunctions::planck(wavelength_list[w], T); //[W m^-2 m^-1 sr^-1]
            double sp_energy;

            if(is_ext)
                sp_energy = sp_ext.getValue(wavelength_list[w]);
            else
                sp_energy = R * pl; //[W m^-1] energy per second and wavelength

            star_emi.push_back(sp_energy);
        }

        tmp_luminosity = CMathFunctions::integ(wavelength_list, star_emi, 0, getNrOfWavelength() - 1);
        L = tmp_luminosity;

        cout << "- Source (" << id + 1 << " of " << max << ") STARFIELD: " << float(L / L_sun)
             << " [L_sun], photons: " << nr_of_photons << endl;

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
             << ", neglected energy: " << float(100.0 * diff_luminosity / tmp_luminosity) << " [%]           "
             << endl;

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

bool CSourceStarField::setParameterFromFile(parameters & param, uint p)
{
    dlist values = param.getDiffuseSources();
    string filename = param.getDiffuseSourceString(p / NR_OF_DIFF_SOURCES);

    ifstream reader(filename.c_str());
    int line_counter = 0;
    string line;
    CCommandParser ps;

    double w_min = 1e300;
    double w_max = 0;

    is_ext = true;

    pos = Vector3D(values[p + 0], values[p + 1], values[p + 2]);
    R = 0;
    T = 0;
    var = values[p + 3];
    nr_of_photons = ullong(values[p + NR_OF_DIFF_SOURCES - 1]);

    cout << CLR_LINE << flush;
    cout << "-> Loading spectrum for source star...           \r" << flush;

    if(reader.fail())
    {
        cout << "\nERROR: Cannot open spectrum file: \n" << filename << "  \n" << endl;
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
            cout << "\nERROR: In spectrum file:\n" << filename << endl;
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

void CSourceStarField::createNextRay(photon_package * pp, ullong i_phot)
{
    StokesVector tmp_stokes_vector;
    double energy;
    uint wID;

    pp->initRandomGenerator(i_phot);
    pp->calcRandomDirection();

    double len = pp->randn(0, var);

    pos = len * pp->getDirection();
    pos += pos;

    if(pp->getDustWavelengthID() != MAX_UINT)
    {
        wID = pp->getDustWavelengthID();
        if(is_ext)
        {
            energy = sp_ext.getValue(wavelength_list[wID]) / nr_of_photons;
            double tmp_q = sp_ext_q.getValue(wavelength_list[wID]);
            double tmp_u = sp_ext_u.getValue(wavelength_list[wID]);
            tmp_stokes_vector = energy * StokesVector(1.0, tmp_q, tmp_u, 0);
        }
        else
        {
            double pl = CMathFunctions::planck(wavelength_list[wID], T);
            energy = 4.0 * PI * PI * (R * R_sun) * (R * R_sun) * pl / nr_of_photons;
            tmp_stokes_vector = energy * StokesVector(1.0, q, u, 0);
        }
    }
    else
    {
        energy = L / nr_of_photons;
        wID = lam_pf.getXIndex(pp->getRND());

        tmp_stokes_vector = energy * StokesVector(1.0, 0, 0, 0);
        // Mol3D uses the upper value of the wavelength interval,
        // used for the selection of the emitting wavelengths from source!
        pp->setWavelength(wavelength_list[wID + 1], wID + 1);
    }

    pp->setPosition(pos);
    pp->setStokesVector(tmp_stokes_vector);
    pp->initCoordSystem();
}

void CSourceStarField::createDirectRay(photon_package * pp, Vector3D dir_obs)
{
    StokesVector tmp_stokes_vector;
    double energy;
    uint wID;

    double len = pp->randn(0, var);

    pos = len * pp->getDirection();
    pos += pos;

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
        wID = lam_pf.getXIndex(pp->getRND());

        tmp_stokes_vector = energy * StokesVector(1.0, 0, 0, 0);
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

    pp->setPosition(pos);
    pp->setStokesVector(tmp_stokes_vector);
}

bool CSourceBackground::initSource(uint id, uint max, bool use_energy_density)
{
    lam_pf = new spline[max_len];
    L = new double[max_len];
    double * star_emi = new double[getNrOfWavelength()];
    sidelength = grid->getMaxLength();
    step_xy = sidelength / double(bins);
    off_xy = step_xy / 2.0;

    cout << CLR_LINE;

    if(constant)
    {
        cout << "Initiating constant background source             \r";

        for(uint w = 0; w < getNrOfWavelength(); w++)
        {
            double pl = 0.0;
            double sp_energy = 0.0;

            if(c_f>=0)
            {
                pl=CMathFunctions::planck(wavelength_list[w], c_temp); //[W m^-2 m^-1]
                sp_energy = abs(c_f) * pl; //[W m^-1] energy per second an wavelength
            }
            else
            {
                sp_energy = abs(c_f);
            }

            star_emi[w] = sp_energy;
        }

        L[0] = CMathFunctions::integ(wavelength_list, star_emi, 0, getNrOfWavelength() - 1);

        double fr, sum = 0;
        lam_pf[0].resize(getNrOfWavelength());

        for(uint w = 0; w < getNrOfWavelength(); w++)
        {
            if(w > 0)
                sum += (wavelength_list[w] - wavelength_list[w - 1]) * star_emi[w - 1] +
                       0.5 * (wavelength_list[w] - wavelength_list[w - 1]) * (star_emi[w] - star_emi[w - 1]);

            fr = sum / L[0];
            lam_pf[0].setValue(w, fr, double(w));
        }

        cout << CLR_LINE;
        if(nr_of_photons==0)
        {
            cout << "Source (" << id + 1 << " of " << max << ") BACKGROUND (const.) initiated \n";
        }
        else
        {
            if(use_energy_density)
                cout << "Source (" << id + 1 << " of " << max << ") BACKGROUND (const.) initiated \n"
                     << "with " << nr_of_photons << " photons per cell and wavelength" << endl;
            else
                cout << "Source (" << id + 1 << " of " << max << ") BACKGROUND (const.) initiated \n"
                     << "with " << nr_of_photons << " photons per cell" << endl;
        }

    }
    else
    {
        cout << CLR_LINE;
        cout << "Initiating variable background source: 0.0 [%]          \r";

        for(uint p = 0; p < max_len; p++)
        {
            double F = f(p);
            double T = temp(p);

            double fr, sum = 0;
            lam_pf[p].resize(getNrOfWavelength());

            for(uint w = 0; w < getNrOfWavelength(); w++)
            {
                double pl = CMathFunctions::planck(wavelength_list[w], T); //[W m^-2 m^-1]
                double sp_energy = F * pl; //[W m^-1] energy per second an wavelength
                star_emi[w] = sp_energy;
            }

            L[p] = CMathFunctions::integ(wavelength_list, star_emi, 0, getNrOfWavelength() - 1);

            for(uint w = 0; w < getNrOfWavelength(); w++)
            {
                if(w > 0)
                    sum +=
                        (wavelength_list[w] - wavelength_list[w - 1]) * star_emi[w - 1] +
                        0.5 * (wavelength_list[w] - wavelength_list[w - 1]) * (star_emi[w] - star_emi[w - 1]);

                fr = sum / L[p];
                lam_pf[p].setValue(w, fr, double(w));
            }

            if(p % 500 == 0)
                cout << "Initiating variable background source: " << 100 * float(p) / float(max_len)
                     << " [%]       \r";
        }

        if(use_energy_density)
            cout << "Source (" << id + 1 << " of " << max << ") BACKGROUND (var.) initiated \n"
                 << "with " << nr_of_photons << " photons per cell and wavelength" << endl;
        else
            cout << "Source (" << id + 1 << " of " << max << ") BACKGROUND (var.) initiated \n"
                 << "with " << nr_of_photons << " photons per cell" << endl;
    }

    // temp.clear();
    // f.clear();

    delete[] star_emi;
    return true;
}

bool CSourceBackground::setParameterFromFile(parameters & param, uint p)
{
    dlist values = param.getDiffuseSources();
    string filename = param.getBackgroundSourceString(p / NR_OF_BG_SOURCES);

    ifstream reader(filename.c_str());
    int line_counter = -2;
    string line;
    CCommandParser ps;

    rot_angle1 = values[p + 5];
    rot_angle2 = values[p + 6];
    nr_of_photons = ullong(values[p + 7]);

    if(reader.fail())
    {
        cout << "\nERROR: Cannot open file:\n" << filename << endl;
        return false;
    }

    while(getline(reader, line))
    {
        ps.formatLine(line);

        if(line.size() == 0)
            continue;

        dlist value = ps.parseValues(line);

        if(value.size() == 0)
            continue;

        line_counter++;

        if(line_counter == -1)
        {
            if(value.size() != 1)
            {
                cout << "\nERROR: Wrong amount of values in:\n " << filename << endl;
                cout << "1 value expected in line " << line_counter + 5 << " !" << endl;
                return false;
            }

            bins = uint(value[0]);
            max_len = bins * bins;

            temp.resize(bins, bins);
            f.resize(bins, bins);
            q.resize(bins, bins);
            u.resize(bins, bins);
            v.resize(bins, bins);
        }
        else
        {
            if(line_counter > (int)max_len)
            {
                cout << "\nERROR: To many background values in : " << filename << endl;
                cout << max_len << " lines expected!" << endl;
                return false;
            }

            if(value.size() == 5)
            {
                f.set(uint(line_counter), value[0]);
                temp.set(uint(line_counter), value[1]);
                q.set(uint(line_counter), value[2]);
                u.set(uint(line_counter), value[3]);
                v.set(uint(line_counter), value[4]);
            }
            else
            {
                cout << "\nERROR: File : " << filename << endl;
                cout << " 5 values in line " << line_counter + 1 << " expected!" << endl;
                return false;
            }
        }
    }

    if(line_counter + 1 < (int)max_len)
    {
        cout << "\nERROR: Not enough background values in : " << filename << endl;
        cout << max_len << " lines expected!" << endl;
        return false;
    }

    constant = false;
    reader.close();

    return true;
}

StokesVector CSourceBackground::getStokesVector(photon_package * pp)
{
    double F, T, I, Q, U, V, pl;
    StokesVector res;
    Vector3D pos = pp->getPosition();

    uint wID = pp->getDustWavelengthID();

    uint x = uint((pos.X() + 0.5 * sidelength) / sidelength * double(bins));
    uint y = uint((pos.Y() + 0.5 * sidelength) / sidelength * double(bins));

    if(constant)
    {
        if(c_f>=0)
        {
            F = c_f;
            T = c_temp;
            Q = c_q;
            U = c_u;
            V = c_v;

            pl = CMathFunctions::planck(wavelength_list[wID], T); //[W m^-2 m^-1]
            I = F * pl;                                           //[W m^-1] energy per second and wavelength
            Q *= I;
            U *= I;
            V *= I;
        }
        else
        {
            I = abs(c_f); //[W m^-1]
            Q = c_q*I;
            U = c_u*I;
            V = c_v*I;
        }
    }
    else
    {
        F = f(x, y);
        T = temp(x, y);
        Q = q(x, y);
        U = u(x, y);
        V = v(x, y);

        pl = CMathFunctions::planck(wavelength_list[wID], T); //[W m^-2 m^-1]
        I = F * pl;                                           //[W m^-1] energy per second and wavelength
        Q *= I;
        U *= I;
        V *= I;
    }

    res.set(I, Q, U, V);

    return res;
}

bool CSourceISRF::initSource(uint id, uint max, bool use_energy_density)
{
    double * star_emi = new double[getNrOfWavelength()];

    cout << CLR_LINE << flush;
    cout << "-> Initiating interstellar radiation field          \r" << flush;

    for(uint w = 0; w < getNrOfWavelength(); w++)
    {
        double pl = sp_ext.getValue(wavelength_list[w]); //[W m^-2 m^-1 sr^-1]
        double sp_energy = pl * PI * PI * 3 * pow(radius * grid->getMaxLength(), 2);    //[W m^-1] energy per second and wavelength
        // if(g_zero > 0)
        //     sp_energy *= PIx2; //[W m^-1]
        star_emi[w] = sp_energy;
    }

    L = CMathFunctions::integ(wavelength_list, star_emi, 0, getNrOfWavelength() - 1);

    if(use_energy_density)
    {
        // For using energy density, only the photon number is required
        cout << "- Source (" << id + 1 << " of " << max << ") ISRF initiated with " << nr_of_photons
             << " photons FOR EACH wavelength"
             << "      " << endl;
    }
    else
    {
        double fr, sum = 0;
        lam_pf.resize(getNrOfWavelength());

        for(uint w = 0; w < getNrOfWavelength(); w++)
        {
            if(w > 0)
                sum += (wavelength_list[w] - wavelength_list[w - 1]) * star_emi[w - 1] +
                       0.5 * (wavelength_list[w] - wavelength_list[w - 1]) * (star_emi[w] - star_emi[w - 1]);

            fr = sum / L;
            lam_pf.setValue(w, fr, double(w));
        }
        cout << "- Source (" << id + 1 << " of " << max << ") ISRF initiated with " << nr_of_photons
             << " photons"
             << "      " << endl;
    }

    cout << "    Radius multiplier: " << radius << ", ";
    if(g_zero > 0)
        cout << "G_0: " << g_zero << " (see Mathis et al. 1983)" << endl;
    else
        cout << "luminosity: " << float(L / L_sun) << " [L_sun]" << endl;

    delete[] star_emi;

    return true;
}

bool CSourceISRF::setParameterFromFile(parameters & param, uint p)
{
    nr_of_photons = param.getNrOfISRFPhotons();
    string filename = param.getISRFPath();
    radius = param.getISRFRadius();

    spline sp_ext_wl;

    ifstream reader(filename.c_str());
    int line_counter = -2;
    string line;
    CCommandParser ps;

    double w_min = 1e300;
    double w_max = 0;

    cout << CLR_LINE << flush;
    cout << "-> Loading spectrum for source ISRF ...           \r" << flush;

    if(reader.fail())
    {
        cout << "\nERROR: Cannot open file: " << filename << endl;
        return false;
    }

    while(getline(reader, line))
    {
        ps.formatLine(line);

        if(line.size() == 0)
            continue;

        dlist value = ps.parseValues(line);

        if(value.size() == 0)
            continue;

        line_counter++;

        if(line_counter == -1)
        {
            if(value.size() != 3)
            {
                cout << "\nERROR: In ISRF file:\n" << filename << endl;
                cout << "Wrong amount of values in line " << line_counter + 1 << "!" << endl;
                return false;
            }

            c_q = value[0];
            c_u = value[1];
            c_v = value[2];
        }
        else if(line_counter >= 0)
        {
            if(value.size() == 2)
            {
                sp_ext_wl.setDynValue(value[0], value[1]);

                if(w_min > value[0])
                    w_min = value[0];

                if(w_max < value[0])
                    w_max = value[0];
            }
            else
            {
                cout << "\nERROR: In ISRF file:\n" << filename << endl;
                cout << "Wrong amount of values in line " << line_counter + 1 << "!" << endl;
                return false;
            }
        }
    }

    sp_ext_wl.createDynSpline();
    reader.close();

    sp_ext.resize(getNrOfWavelength());
    for(uint w = 0; w < getNrOfWavelength(); w++)
    {
        double rad_field = 0;
        if(wavelength_list[w] > w_min && wavelength_list[w] < w_max)
        {
            // Get radiation field from dynamic spline
            // Divide by 4PI to get per steradian to match mathis_isrf
            rad_field = sp_ext_wl.getValue(wavelength_list[w]) / PIx4;
        }

        // Calculate final emission
        sp_ext.setValue(w, wavelength_list[w], rad_field * dust->getForegroundExtinction(wavelength_list[w]));
    }
    sp_ext.createSpline();

    return true;
}

void CSourceISRF::createNextRay(photon_package * pp, ullong i_phot)
{
    double energy;
    StokesVector tmp_stokes_vector;
    pp->initRandomGenerator(i_phot);
    uint wID;

    if(pp->getDustWavelengthID() != MAX_UINT)
    {
        wID = pp->getDustWavelengthID();
        double pl = sp_ext.getValue(wavelength_list[wID]); //[W m^-2 m^-1 sr^-1]
        energy = pl * PI * PI * 3 * pow(radius * grid->getMaxLength(), 2) / nr_of_photons; //[W m^-1] energy per second and wavelength
        // if(g_zero > 0)
        //     energy *= PIx2;
    }
    else
    {
        wID = lam_pf.getXIndex(pp->getRND());
        energy = L / nr_of_photons;

        // Mol3D uses the upper value of the wavelength interval,
        // used for the selection of the emitting wavelengths from source!
        pp->setWavelength(wavelength_list[wID + 1], wID + 1);
    }

    tmp_stokes_vector = energy * StokesVector(1, c_q, c_u, c_v);

    // Get random direction for the postion (not the direction of travel) of the photon
    pp->calcRandomDirection();
    pp->initCoordSystem();

    // pos is center position of the sphere where the ISRF is emitted from
    pos = Vector3D(0,0,0);
    // set photon position to be on the surface of that sphere
    // negative sign, so that the direction of travel is inside of the cell
    pp->adjustPosition(pos, -sqrt(3) * radius * grid->getMaxLength() / 2);

    // the emission has to obey Lambert's cosine law
    // Thus, the square of cos(theta) of the direction of the photon
    // is given by a random number -> theta is only in (0,pi/2)
    double theta_direction = acos( sqrt( pp->getRND() ) );
    double phi_direction = PIx2 * pp->getRND();

    // the direction was calculated in the coord system of the photon
    // now we have to update the coord system acoordingly
    pp->updateCoordSystem(phi_direction, theta_direction);

    pp->setStokesVector(tmp_stokes_vector);
}

void CSourceISRF::createDirectRay(photon_package * pp, Vector3D dir_obs)
{
    // THIS IS PROBABLY STILL WRONG
    // actually, this function should never be used
    // because it puts all the energy of the isrf
    // into a single photon

    double energy;
    StokesVector tmp_stokes_vector;
    uint wID;

    if(pp->getDustWavelengthID() != MAX_UINT)
    {
        wID = pp->getDustWavelengthID();
        double pl = sp_ext.getValue(wavelength_list[wID]); //[W m^-2 m^-1 sr^-1]
        energy = pl * PI * 3 * pow(radius * grid->getMaxLength(), 2) / PIx2;
        if(g_zero > 0)
            energy *= PIx2; //[W m^-1] energy per second an wavelength
    }
    else
    {
        wID = lam_pf.getXIndex(pp->getRND());
        energy = L / PIx2;

        // Mol3D uses the upper value of the wavelength interval,
        // used for the selection of the emitting wavelengths from source!
        pp->setWavelength(wavelength_list[wID + 1], wID + 1);
    }

    tmp_stokes_vector = energy * StokesVector(1, c_q, c_u, c_v);

    e.rndDir(pp->getRND(), pp->getRND());
    l.setX(e.X() * sqrt(3) * radius * grid->getMaxLength() / 2);
    l.setY(e.Y() * sqrt(3) * radius * grid->getMaxLength() / 2);
    l.setZ(e.Z() * sqrt(3) * radius * grid->getMaxLength() / 2);

    // Set direction of the photon package to the observer
    if(dir_obs.length() > 0)
    {
        pp->setDirection(dir_obs);
        pp->initCoordSystem();
    }

    pos.set(l.X(), l.Y(), l.Z());

    pp->setPosition(pos);
    pp->setStokesVector(tmp_stokes_vector);
}

bool CSourceDust::initSource(uint id, uint max, bool use_energy_density)
{
    if(!use_energy_density)
    {
        cout << "\nERROR: The dust source for radiation field calculation "
             << "can only be used with energy density!" << endl;
        return false;
    }

    // Init variables
    ulong nr_of_cells = grid->getMaxDataCells();
    ulong nr_of_wavelengths = getNrOfWavelength();
    photon_package pp = photon_package();
    cell_prob = new prob_list[nr_of_wavelengths];
    total_energy = new double[nr_of_wavelengths];
    uint max_counter = nr_of_wavelengths * nr_of_cells;
    ullong per_counter = 0;
    float last_percentage = 0;

    // Show Initial message
    cout << CLR_LINE;
    cout << "-> Initiating dust grain emission          \r" << flush;

    for(uint w = 0; w < nr_of_wavelengths; w++)
    {
        // Init variables
        cell_prob[w].resize(nr_of_cells + 1);

        // Set wavelength of photon package
        pp.setWavelength(wavelength_list[w], w);

        // Set total energy to zero and starting value of prob_list
        total_energy[w] = 0;
        cell_prob[w].setValue(0, total_energy[w]);

// Causes problems. Find better solution!
//#pragma omp parallel for schedule(dynamic)
        for(long i_cell = 0; i_cell < long(nr_of_cells); i_cell++)
        {
            // Increase counter used to show progress
            per_counter++;

            // Calculate percentage of total progress
            float percentage = 100.0 * float(per_counter) / float(max_counter);

            // Show only new percentage number if it changed
            if((percentage - last_percentage) > PERCENTAGE_STEP)
            {
//#pragma omp critical
                {
                    cout << "-> Calculating prob. distribution for dust source: " << percentage << " [%]    \r"
                         << flush;
                    last_percentage = percentage;
                }
            }

            // Put photon package into current cell
            pp.setPositionCell(grid->getCellFromIndex(i_cell));

            // Get total energy of thermal emission
            total_energy[w] += dust->getTotalCellEmission(grid, pp);

            // Add energy to probability distribution
            cell_prob[w].setValue(i_cell + 1, total_energy[w]);
        }

        // Normalize probability distribution
        cell_prob[w].normalize(total_energy[w]);
    }

    // Show information
    cout << CLR_LINE;
    cout << "- Source (" << id + 1 << " of " << max
         << ") DUST: photons FOR EACH wavelength: " << nr_of_photons << "      " << endl;

    return true;
}

bool CSourceDust::initSource(uint w)
{
    // Init variables
    ulong nr_of_cells = grid->getMaxDataCells();
    photon_package pp = photon_package();
    uint max_counter = nr_of_cells;
    ullong per_counter = 0;
    float last_percentage = 0;

    // Init variables
    total_energy = new double[getNrOfWavelength()];
    cell_prob = new prob_list[getNrOfWavelength()];
    cell_prob[w].resize(nr_of_cells + 1);

    // Set wavelength of photon package
    pp.setWavelength(wavelength_list[w], w);

    // Set total energy to zero and starting value of prob_list
    total_energy[w] = 0;
    cell_prob[w].setValue(0, total_energy[w]);

    // Show Initial message
    cout << "-> Initiating dust grain emission          \r" << flush;

    #if (DUST_EMI_PROB)
// Causes problems. Find better solution!
//#pragma omp parallel for schedule(dynamic)
        for(long i_cell = 0; i_cell < long(nr_of_cells); i_cell++)
        {
            // Increase counter used to show progress
//#pragma omp atomic update
            per_counter++;

            // Calculate percentage of total progress
            float percentage = 100.0 * float(per_counter) / float(max_counter);

            // Show only new percentage number if it changed
            if((percentage - last_percentage) > PERCENTAGE_STEP)
            {
//#pragma omp critical
                {
                    cout << "-> Calculate prob. distribution for dust source: " << percentage << " [%]    \r"
                        << flush;
                    last_percentage = percentage;
                }
            }

            // Put photon package into current cell
            pp.setPositionCell(grid->getCellFromIndex(i_cell));

            // Get total energy of thermal emission
            total_energy[w] += dust->getTotalCellEmission(grid, pp);

            // Add energy to probability distribution
            cell_prob[w].setValue(i_cell + 1, total_energy[w]);
        }

        // Normalize probability distribution
        cell_prob[w].normalize(total_energy[w]);
    #else
        // reduce the total number of photons, so that each cell launches the same amount of photons,
        // i.e. (n_photon % n_cell) should be zero
        nr_of_photons -= nr_of_photons % grid->getMaxDataCells();
        nr_of_photons_per_cell = ullong(nr_of_photons / double(nr_of_cells));
    #endif

    return true;
}

void CSourceDust::createNextRay(photon_package * pp, ullong i_phot)
{
    // Init photon package and random direction
    pp->initRandomGenerator(i_phot);
    pp->calcRandomDirection();

    // Set wavelength of photon package
    uint w = pp->getDustWavelengthID();

    // Get index of current cell
    ulong i_cell = 0;
    #if (DUST_EMI_PROB)
        // Get random number
        double rnd = pp->getRND();

        i_cell = cell_prob[w].getIndex(rnd);
    #else
        i_cell = ulong(i_phot % grid->getMaxDataCells());
    #endif

    // Put photon package into current cell
    pp->setPositionCell(grid->getCellFromIndex(i_cell));

    // Set random position in cell
    grid->setRndPositionInCell(pp);

    // Set Stokes vector of photon package
    double energy = 0;
    #if (DUST_EMI_PROB)
        energy = total_energy[w] / double(nr_of_photons);
    #else
        energy = dust->getTotalCellEmission(grid, *pp) / double(nr_of_photons_per_cell);
    #endif

    // Set Stokes Vector
    pp->setStokesVector(StokesVector(energy, 0, 0, 0));

    // Init coordinate System for polarization
    pp->initCoordSystem();
}

void CSourceDust::createDirectRay(photon_package * pp, Vector3D dir_obs)
{
    // Set wavelength of photon package
    uint w = pp->getDustWavelengthID();

    // Get random number
    double rnd = pp->getRND();

    // Get index of current cell
    ulong i_cell = cell_prob[w].getIndex(rnd);

    // Put photon package into current cell
    pp->setPositionCell(grid->getCellFromIndex(i_cell));

    // Set random position in cell
    grid->setRndPositionInCell(pp);

    // Set Stokes vector of photon package
    double energy = total_energy[w] / PIx4;

    // Set Stokes Vector
    pp->setStokesVector(StokesVector(energy, 0, 0, 0));

    // Set direction of the photon package to the observer
    if(dir_obs.length() > 0)
    {
        pp->setDirection(dir_obs);
        pp->initCoordSystem();
    }
}

bool CSourceGas::initSource(uint id, uint max, bool use_energy_density)
{
    // Init variables
    ulong nr_of_cells = grid->getMaxDataCells();

    // Show information
    cout << "- Source GAS (lvl population):" << endl;
    cout << "    photons per cell: " << nr_of_photons << " (" << nr_of_cells << " cells)" << endl;

    return true;
}

void CSourceGas::createNextRayToCell(photon_package * pp, ullong i_phot, ulong i_cell, bool cell_as_border)
{
    // Init photon package and random direction
    pp->initRandomGenerator(i_phot);
    pp->calcRandomDirection();

    // Put photon package into current cell
    pp->setPositionCell(grid->getCellFromIndex(i_cell));

    // Set random position in cell
    grid->setRndPositionInCell(pp);

    // Save final position as position of last interaction
    pp->setBackupPosition(pp->getPosition());

    if(cell_as_border)
    {
        // Go one cell outward
        grid->next(pp);

        // Invert direction
        pp->multDirection(-1);

        // Go into cell again
        grid->next(pp);
    }
    else
    {
        // Move photon along the path outwards the grid
        pp->addPosition(pp->getDirection() * grid->maxLength());

        // Invert direction
        pp->multDirection(-1);
    }

    // // Init coordinate System for polarization
    pp->initCoordSystem();
}

bool CSourceLaser::initSource(uint id, uint max, bool use_energy_density)
{
    // Initial output
    cout << CLR_LINE << flush;
    cout << "-> Initiating source laser         \r" << flush;

    if(use_energy_density)
    {
        // For using energy density, only the photon number is required
        cout << "- Source (" << id + 1 << " of " << max << ") LASER: " << float(L)
             << " [W] (wavelength = " << float(wl) << " [m], FWHM = " << float(fwhm) << " [m])" << endl;
        cout << "    photons FOR EACH wavelength: " << nr_of_photons << "      " << endl;
        // cout << "\nERROR: Laser source cannot be used for dust temperature calculation!\n" << endl;
        // return false;
    }
    else
    {
        // Init variables
        dlist star_emi;
        double tmp_luminosity;

        for(uint w = 0; w < getNrOfWavelength(); w++)
        {
            double line_shape =
                1 / sqrt(2 * PI * sigma_sq) * exp(-pow(wavelength_list[w] - wl, 2) / (2 * sigma_sq));
            star_emi.push_back(L * line_shape);
        }

        tmp_luminosity = CMathFunctions::integ(wavelength_list, star_emi, 0, getNrOfWavelength() - 1);
        // L = tmp_luminosity;

        // For using energy density, only the photon number is required
        cout << "- Source (" << id + 1 << " of " << max << ") LASER: " << float(L)
             << " [W] (wavelength = " << float(wl) << " [m], FWHM = " << float(fwhm) << " [m])" << endl;
        cout << "    photons: " << nr_of_photons << "      " << endl;

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

void CSourceLaser::createNextRay(photon_package * pp, ullong i_phot)
{
    // Init variables
    StokesVector tmp_stokes_vector;
    uint wID = 0;

    pp->initRandomGenerator(i_phot);
    pp->setDirection(dir);
    pp->setPosition(pos);

    if(pp->getDustWavelengthID() != MAX_UINT)
        wID = pp->getDustWavelengthID();
    else
    {
        wID = lam_pf.getXIndex(pp->getRND());
        // Mol3D uses the upper value of the wavelength interval,
        // used for the selection of the emitting wavelengths from source!
        pp->setWavelength(wavelength_list[wID + 1], wID + 1);
    }

    double line_shape =
        1 / sqrt(2 * PI * sigma_sq) * exp(-pow(wavelength_list[wID] - wl, 2) / (2 * sigma_sq));
    tmp_stokes_vector = L / double(nr_of_photons) * line_shape * StokesVector(1.0, q, u, 0);
    pp->setStokesVector(tmp_stokes_vector);
    pp->initCoordSystem();
}

void CSourceLaser::createDirectRay(photon_package * pp, Vector3D dir_obs)
{
    // Init variables
    StokesVector tmp_stokes_vector;

    pp->setDirection(dir);
    pp->setPosition(pos);

    if(dir_obs == dir && pp->getDustWavelengthID() != MAX_UINT)
    {
        uint wID = pp->getDustWavelengthID();
        double line_shape =
            1 / sqrt(2 * PI * sigma_sq) * exp(-pow(wavelength_list[wID] - wl, 2) / (2 * sigma_sq));
        tmp_stokes_vector = L * line_shape * StokesVector(1.0, q, u, 0);
    }

    pp->setStokesVector(tmp_stokes_vector);

    // Set direction of the photon package to the observer
    if(dir_obs.length() > 0)
    {
        pp->initCoordSystem();
    }
}
