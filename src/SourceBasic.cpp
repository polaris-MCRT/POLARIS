/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "SourceBasic.hpp"
#include "CommandParser.hpp"
#include "GridBasic.hpp"
#include "MathFunctions.hpp"
#include "Parameters.hpp"
#include "Photon.hpp"

Vector3D CSourceBasic::getPosition()
{
    return pos;
}

double CSourceBasic::getRadius()
{
    return R;
}

double CSourceBasic::getLuminosity()
{
    return L;
}

double CSourceBasic::getTemperature()
{
    return T;
}

double CSourceBasic::getSublimationRadius()
{
    return r_sub;
}

void CSourceBasic::setParameter(parameters & param, CGridBasic * _grid, CDustMixture * _dust, uint p)
{
    setGrid(_grid);
    setDust(_dust);

    // Get global wavelength grid
    wavelength_list = _dust->getWavelengthList();

    setParameter(param, p);
}

void CSourceBasic::setParameter(parameters & param, uint p)
{}

bool CSourceBasic::setParameterFromFile(parameters & param, CGridBasic * _grid, CDustMixture * _dust, uint p)
{
    setGrid(_grid);
    setDust(_dust);

    // Get global wavelength grid
    wavelength_list = _dust->getWavelengthList();

    return setParameterFromFile(param, p);
}

bool CSourceBasic::setParameterFromFile(parameters & param, uint p)
{
    return false;
}

void CSourceBasic::setGrid(CGridBasic * _grid)
{
    grid = _grid;
}

void CSourceBasic::setDust(CDustMixture * _dust)
{
    dust = _dust;
}

bool CSourceBasic::initSource(uint id, uint max, bool use_energy_density)
{
    return true;
}

bool CSourceBasic::initSource(uint w)
{
    return true;
}

void CSourceBasic::clean()
{}

uint CSourceBasic::getID()
{
    return source_id;
}

void CSourceBasic::setNrOfPhotons(ullong val)
{
    nr_of_photons = val;
}

void CSourceBasic::updateNrOfPhotons(double val)
{
    nr_of_photons = ullong(nr_of_photons * val);
}

void CSourceBasic::setSideLength(double val)
{}

void CSourceBasic::createNextRay(photon_package * pp, CRandomGenerator * rand_gen)
{}

void CSourceBasic::createNextRayToCell(photon_package * pp,
                                       CRandomGenerator * rand_gen,
                                       ulong i_cell,
                                       bool cell_as_border)
{}

void CSourceBasic::createDirectRay(photon_package * pp, CRandomGenerator * rand_gen, Vector3D dir_obs)
{}

ullong CSourceBasic::getNrOfPhotons()
{
    return nr_of_photons;
}

uint CSourceBasic::getNrOfWavelength()
{
    return wavelength_list.size();
}

void CSourceBasic::setOrientation(Vector3D n1, Vector3D n2, double _theta, double _phi)
{}

StokesVector CSourceBasic::getStokesVector(photon_package * pp)
{
    return StokesVector(0, 0, 0, 0);
}

uint CSourceBasic::getBins()
{
    return 1;
}
