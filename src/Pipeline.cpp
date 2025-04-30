/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "Pipeline.hpp"
#include "CommandParser.hpp"
#include "GridCylindrical.hpp"
#include "GasSpecies.hpp"
#include "GridBasic.hpp"
#include "GridOcTree.hpp"
#include "RadiativeTransfer.hpp"
#include "SourceBasic.hpp"
#include "SourceStar.hpp"
#include "SourceStarField.hpp"
#include "SourceAGN.hpp"
#include "SourceISRF.hpp"
#include "SourceBackground.hpp"
#include "SourceDust.hpp"
#include "SourceGas.hpp"
#include "SourceLaser.hpp"
#include "GridSpherical.hpp"
#include "GridVoronoi.hpp"
#include "OPIATE.hpp"
#include "Detector.hpp"
#include "MathFunctions.hpp"

bool CPipeline::Init(int argc, char ** argv)
{
    end = 0, len = 0;
    begin = omp_get_wtime();
    srand(0);
    CMathFunctions mf;

    string out_string_prog = "*                      ";
    out_string_prog += PROG_ID;
    out_string_prog += "                       *";

    string out_string_vers = "*                      ";
    out_string_vers += VERS_ID;
    out_string_vers += "                       *";

    string out_string_copy = "*                      ";
    out_string_copy += COPY_ID;
    out_string_copy += "                       *";

    cout << SEP_LINE;

    cout << out_string_prog << endl;

    cout << CLR_LINE;

    cout << out_string_vers << endl;

    cout << CLR_LINE;

    cout << out_string_copy << endl;

    cout << SEP_LINE;

    /*if(argc != 2)
    {
        cout << ERROR_LINE << "Wrong amount of arguments!                     \n";
        cout << "       POLARIS requires only the path of a command file!            \n";
        Error();
        return false;
    }

    CCommandParser parser(argv[1]);*/
    
    string filename = "/mnt/c/Users/Stefan/Documents/NetBeansProjects/test_efrem/src/cmd_test";
    
    CCommandParser parser(filename);
    

    if(!parser.parse())
    {
        Error();
        return false;
    }

    param_list = parser.getParameterList();

    if(param_list.size() == 0)
    {
        cout << ERROR_LINE << "No tasks defined!" << endl;
        return false;
    }

    return true;
}

void CPipeline::Finish()
{
    end = omp_get_wtime();
    len = (long)(end - begin);

    h = len / 3600;
    m = (len - h * 3600) / 60;
    s = len - h * 3600 - m * 60;

    cout << SEP_LINE;
    printf("  Total time of processing: %lu h %02lu min %02lu sec  \n", h, m, s);

    cout << CLR_LINE;
    cout << SEP_LINE;
    cout << "* POLARIS SUCCESSFULLY FINISHED                                             "
            "        *"
         << endl;
    cout << SEP_LINE << endl;

    string opath = path_data + "time.txt";
    ofstream t_writer(opath.c_str());
    t_writer << "Total time:\r" << endl;
    t_writer << h << " h " << m << " min " << s << " sec\r\n\r" << endl;
    t_writer << len << " s\r" << endl;
    t_writer.close();

    /*#ifdef WINDOWS
     cout << "press any key ..." << endl;
     #pragma warning(suppress: 6031)
     getchar();
     #endif*/
}

void CPipeline::Run()
{
    bool result = false;
    uint size = uint(param_list.size());

    for(uint i = 0; i < size; i++)
    {
        parameters & param = param_list[i];

        omp_set_num_threads(param.getNrOfThreads());
        printParameters(param, size);

        switch(param.getCommand())
        {
            case CMD_TEMP:
                result = calcMonteCarloRadiationField(param);
                break;

            case CMD_TEMP_RAT:
                result = calcMonteCarloRadiationField(param);
                break;

            case CMD_RAT:
                result = calcMonteCarloRadiationField(param);
                break;

            case CMD_DUST_EMISSION:
                result = calcPolarizationMapsViaRayTracing(param);
                break;

            case CMD_DUST_SCATTERING:
                result = calcPolarizationMapsViaMC(param);
                break;

            case CMD_LINE_EMISSION:
                result = calcChMapsViaRayTracing(param);
                break;

            case CMD_FORCE:
                // Needs update!!!
                // result = calcRadPressure(param);
                break;

            case CMD_OPIATE:
                result = calcOpiateMapsViaRayTracing(param);
                break;

            case CMD_SYNCHROTRON:
                result = calcPolarizationMapsViaSynchrotron(param);
                break;

            default:
                cout << ERROR_LINE << "Command is unknown!" << endl;
        }

        cout << SEP_LINE << endl;
    }

    if(result)
        Finish();
    else
        Error();
}

void CPipeline::Error()
{
    end = omp_get_wtime();
    len = (long)(end - begin);

    h = len / 3600;
    m = (len - h * 3600) / 60;
    s = len - h * 3600 - m * 60;

    cout << SEP_LINE;
    printf("  Total time of processing: %lu h %02lu min %02lu sec  \n", h, m, s);

    cout << CLR_LINE;
    cout << SEP_LINE;
    cout << "* POLARIS ABORTED PROCESSING	                                         "
            "           *"
         << endl;
    cout << SEP_LINE << endl;
}

bool CPipeline::calcMonteCarloRadiationField(parameters & param)
{
    // Check if the energy density is used instead of launching photons with fixed energy
    // In case of (save radiation field), (calc RATs), and (calc stochastic heating
    // temperatures)
    bool use_energy_density = false;
    if(param.getSaveRadiationField() || param.isRatSimulation() || param.getStochasticHeatingMaxSize() > 0)
        use_energy_density = true;

    CGridBasic * grid = 0;
    CDustMixture * dust = new CDustMixture();

    if(!createOutputPaths(param.getPathOutput()))
        return false;

    if(!assignGridType(grid, param))
        return false;

    if(!createWavelengthList(param, dust, 0, 0))
        return false;

    if(!assignDustMixture(param, dust, grid))
        return false;

    grid->setSIConversionFactors(param);

    grid->setSpecLengthAsVector(use_energy_density);
    if(!grid->loadGridFromBinaryFile(param, use_energy_density ? 4 * WL_STEPS : WL_STEPS))
        return false;

    // Print helpfull information
    grid->createCellList();
    dust->printParameters(param, grid);
    grid->printParameters();

    if(!grid->writeMidplaneFits(path_data + "input_", param, param.getInpMidDataPoints(), true))
        return false;

    if(!grid->writePlotFiles(path_plot + "input_", param))
        return false;

    if(!grid->writeAMIRAFiles(path_plot + "input_", param, param.getInpAMIRAPoints()))
        return false;

    createSourceLists(param, dust, grid);
    if(sources_mc.size() == 0)
    {
        cout << ERROR_LINE << "No sources for Monte-Carlo simulations defined!" << endl;
        return false;
    }

    CRadiativeTransfer rad(param);

    rad.setGrid(grid);
    rad.setDust(dust);
    rad.setSourcesLists(sources_mc, sources_ray);
    rad.initiateRadFieldMC(param);

    if(param.isTemperatureSimulation())
    {
        if(param.getDustOffset())
            rad.convertTempInQB(param.getOffsetMinGasDensity(), false);
        else if(param.getDustGasCoupling())
            rad.convertTempInQB(param.getOffsetMinGasDensity(), true);
    }

    rad.calcMonteCarloRadiationField(param.getCommand(),
                                     use_energy_density,
                                     false); //(param.getCommand() == CMD_RAT));

    if(param.isTemperatureSimulation())
        rad.calcFinalTemperature(use_energy_density);

    if(param.isRatSimulation())
        rad.calcAlignedRadii();

    cout << SEP_LINE;

    if(!grid->writeMidplaneFits(path_data + "output_", param, param.getOutMidDataPoints()))
        return false;

    if(!grid->writePlotFiles(path_plot + "output_", param))
        return false;

    if(!grid->writeAMIRAFiles(path_plot + "output_", param, param.getOutAMIRAPoints()))
        return false;

    grid->writeSpecialLines(path_data);

    if(param.getSaveRadiationField())
        grid->saveRadiationField();
    if(param.isTemperatureSimulation())
        grid->saveBinaryGridFile(param.getPathOutput() + "grid_temp.dat");
    else if(param.getCommand() == CMD_RAT)
        grid->saveBinaryGridFile(param.getPathOutput() + "grid_rat.dat");

    delete grid;
    delete dust;
    deleteSourceLists();

    return true;
}

bool CPipeline::calcPolarizationMapsViaMC(parameters & param)
{
    CGridBasic * grid = 0;
    CDustMixture * dust = new CDustMixture();

    if(!createOutputPaths(param.getPathOutput()))
        return false;

    if(!assignGridType(grid, param))
        return false;

    if(!createWavelengthList(param, dust, 0, 0))
        return false;

    if(!assignDustMixture(param, dust, grid))
        return false;

    grid->setSIConversionFactors(param);

    if(!grid->loadGridFromBinaryFile(param, 0))
        return false;

    // Print helpfull information
    grid->createCellList();
    dust->printParameters(param, grid);
    grid->printParameters();

    if(!grid->writeMidplaneFits(path_data + "input_", param, param.getInpMidDataPoints(), true))
        return false;

    if(!grid->writePlotFiles(path_plot + "input_", param))
        return false;

    if(!grid->writeAMIRAFiles(path_plot + "input_", param, param.getInpAMIRAPoints()))
        return false;

    createSourceLists(param, dust, grid);
    if(sources_mc.size() == 0)
    {
        cout << ERROR_LINE << "No sources for Monte-Carlo simulations defined!" << endl;
        return false;
    }

    CDetector * detector = createDetectorList(param, dust, grid);

    if(detector == 0)
        return false;

    CRadiativeTransfer rad(param);

    rad.setGrid(grid);
    rad.setDust(dust);
    rad.setSourcesLists(sources_mc, sources_ray);
    rad.setDetectors(detector);

    rad.initiateDustMC(param);
    rad.calcPolMapsViaMC();

    cout << CLR_LINE;

    delete grid;
    delete dust;
    delete[] detector;
    deleteSourceLists();

    return true;
}

bool CPipeline::calcPolarizationMapsViaRayTracing(parameters & param)
{
    CGridBasic * grid = 0;
    CDustMixture * dust = new CDustMixture();

    if(!createOutputPaths(param.getPathOutput()))
        return false;

    if(!assignGridType(grid, param))
        return false;

    if(!createWavelengthList(param, dust,0, 0))
        return false;

    if(!assignDustMixture(param, dust, grid))
        return false;

    grid->setSIConversionFactors(param);

    if(!grid->loadGridFromBinaryFile(param, getNrOffsetEntriesRay(param, dust, grid)))
        return false;

    // Print helpfull information
    grid->createCellList();
    dust->printParameters(param, grid);
    grid->printParameters();

    if(!grid->writeMidplaneFits(path_data + "input_", param, param.getInpMidDataPoints(), true))
        return false;

    if(!grid->writePlotFiles(path_plot + "input_", param))
        return false;

    if(!grid->writeAMIRAFiles(path_plot + "input_", param, param.getInpAMIRAPoints()))
        return false;

    createSourceLists(param, dust, grid);
    if(sources_ray.size() == 0)
    {
        cout << ERROR_LINE << "No sources for raytracing simulations defined!" << endl;
        return false;
    }

    CRadiativeTransfer rad(param);

    rad.setGrid(grid);
    rad.setDust(dust);
    rad.setSourcesLists(sources_mc, sources_ray);

    if(!rad.initiateDustRaytrace(param))
        return false;

    if(param.getStochasticHeatingMaxSize() > 0)
        rad.calcStochasticHeating();

    // Calculate radiation field before raytracing (if sources defined and no radiation
    // field in grid)
    if(!grid->isRadiationFieldAvailable() && dust->getScatteringToRay() && !sources_mc.empty())
        rad.calcMonteCarloRadiationField(param.getCommand(), true, true);

    if(!rad.calcPolMapsViaRaytracing(param))
        return false;

    cout << CLR_LINE;

    if(!grid->writeMidplaneFits(path_data + "output_", param, param.getOutMidDataPoints()))
        return false;

    delete grid;
    delete dust;
    deleteSourceLists();

    param.setPathInput(path_data);
    param.setPathGrid("");
    param.resetNrOfDustComponents();

    return true;
}

bool CPipeline::calcChMapsViaRayTracing(parameters & param)
{
    CGridBasic * grid = 0;
    CDustMixture * dust = new CDustMixture();
    CGasMixture * gas = new CGasMixture();

    if(!createOutputPaths(param.getPathOutput()))
        return false;

    if(!assignGridType(grid, param))
        return false;

    if(!assignGasSpecies(param, gas, grid))
        return false;

    if(!createWavelengthList(param, dust, gas, 0))
        return false;

    if(!assignDustMixture(param, dust, grid))
        return false;

    grid->setSIConversionFactors(param);

    if(!grid->loadGridFromBinaryFile(param, gas->getNrOffsetEntries(grid, param)))
        return false;

    // Print helpfull information
    grid->createCellList();
    dust->printParameters(param, grid);
    gas->printParameters(param, grid);
    grid->printParameters();

    if(!grid->writeMidplaneFits(path_data + "input_", param, param.getInpMidDataPoints(), true))
        return false;

    if(!grid->writePlotFiles(path_plot + "input_", param))
        return false;

    if(!grid->writeAMIRAFiles(path_plot + "input_", param, param.getInpAMIRAPoints()))
        return false;

    createSourceLists(param, dust, grid);
    if(sources_ray.size() == 0)
    {
        cout << ERROR_LINE << "No sources for raytracing simulations defined!" << endl;
        return false;
    }

    CRadiativeTransfer rad(param);

    rad.setGrid(grid);
    rad.setDust(dust);
    rad.setGas(gas);
    rad.setSourcesLists(sources_mc, sources_ray);

    if(!rad.initiateLineRaytrace(param))
        return false;

    if(!rad.calcChMapsViaRaytracing(param))
        return false;

    delete grid;
    delete dust;
    delete gas;
    deleteSourceLists();

    return true;
}

bool CPipeline::calcOpiateMapsViaRayTracing(parameters & param)
{
    CGridBasic * grid = 0;
    CDustMixture * dust = new CDustMixture();
    //CGasMixture * gas = new CGasMixture();
    COpiateDataBase * op = new COpiateDataBase();

    if(!createOutputPaths(param.getPathOutput()))
        return false;

    if(!assignGridType(grid, param))
        return false;

    if(!op->readOpiateDataBase(param))
        return false;

    if(!createWavelengthList(param, dust, 0, op))
        return false;

    if(!assignDustMixture(param, dust, grid))
        return false;

    grid->setSIConversionFactors(param);

    if(!grid->loadGridFromBinaryFile(param,1))
        return false;

    grid->createCellList();

    // Print helpfull information
    dust->printParameters(param, grid);
    //gas->printParameters(param, grid);
    //op->printParameters(param,grid);
    grid->printParameters();

    if(!grid->writeMidplaneFits(path_data + "input_", param, param.getInpMidDataPoints(), true))
        return false;

    if(!grid->writePlotFiles(path_plot + "input_", param))
        return false;

    if(!grid->writeAMIRAFiles(path_plot + "input_", param, param.getInpAMIRAPoints()))
        return false;

    createSourceLists(param, dust, grid);
    if(sources_ray.size() == 0)
    {
        cout << ERROR_LINE << "No sources for raytracing simulations defined!" << endl;
        return false;
    }

    CRadiativeTransfer rad(param);

    rad.setGrid(grid);
    rad.setDust(dust);
    rad.setOpiateDataBase(op);
    rad.setSourcesLists(sources_mc, sources_ray);

    if(!rad.initiateOPIATERaytrace(param))
        return false;

    if(!rad.calcOPIATEMapsViaRaytracing(param))
        return false;

    delete grid;
    delete dust;
//    delete gas;
    delete op;
    deleteSourceLists();

    return true;
}


bool CPipeline::calcPolarizationMapsViaSynchrotron(parameters & param)
{
    CGridBasic * grid = 0;
    CDustMixture * dust = new CDustMixture();

    if(!createOutputPaths(param.getPathOutput()))
        return false;

    if(!assignGridType(grid, param))
        return false;

    if(!createWavelengthList(param, dust, 0, 0))
        return false;

    if(!assignDustMixture(param, dust, grid))
        return false;

    grid->setSIConversionFactors(param);

    if(!grid->loadGridFromBinaryFile(param, 0))
        return false;

    // Print helpfull information
    grid->createCellList();
    dust->printParameters(param, grid);
    grid->printParameters();

    if(!grid->writeMidplaneFits(path_data + "input_", param, param.getInpMidDataPoints(), true))
        return false;

    if(!grid->writePlotFiles(path_plot + "input_", param))
        return false;

    if(!grid->writeAMIRAFiles(path_plot + "input_", param, param.getInpAMIRAPoints()))
        return false;

    createSourceLists(param, dust, grid);
    if(sources_ray.size() == 0)
    {
        cout << ERROR_LINE << "No sources for raytracing simulations defined!" << endl;
        return false;
    }

    CRadiativeTransfer rad(param);

    rad.setGrid(grid);
    rad.setDust(dust);
    rad.setSourcesLists(sources_mc, sources_ray);

    if(!rad.initiateSyncRaytrace(param))
        return false;

    if(!rad.calcSyncMapsViaRaytracing(param))
        return false;

    delete grid;
    delete dust;
    deleteSourceLists();

    param.setPathInput(path_data);
    param.setPathGrid("");
    param.resetNrOfDustComponents();

    return true;
}

bool CPipeline::assignGridType(CGridBasic *& grid, parameters & param)
{
    string filename = param.getPathGrid();
    ifstream bin_reader(filename.c_str(), ios::in | ios::binary);
    ushort tmpID;

    if(bin_reader.fail())
    {
        cout << ERROR_LINE << "Cannot open binary grid file:" << endl;
        cout << filename << "\n" << endl;
        return false;
    }

    bin_reader.read((char *)&tmpID, 2);
    bin_reader.close();

    switch(tmpID)
    {
        case 0:
        case 1:
        case 6:
        case 7:
        case GRID_ID_OCT:
            grid = new CGridOcTree();
            break;

        case GRID_ID_SPH:
            grid = new CGridSpherical();
            break;

        case GRID_ID_CYL:
            grid = new CGridCylindrical();
            break;

        case GRID_ID_VOR:
            grid = new CGridVoronoi();
            break;

        default:
            cout << ERROR_LINE << "Grid type unknown!" << endl;
            return false;
            break;
    }

    return true;
}

CDetector * CPipeline::createDetectorList(parameters & param, CDustMixture * dust, CGridBasic * grid)
{
    CDetector * detector;
    dlist dust_mc_detectors = param.getDustMCDetectors();
    uint nr_mc_detectors = param.getNrOfDustMCDetectors();

    Vector3D axis1 = param.getAxis1();
    Vector3D axis2 = param.getAxis2();

    bool error = false;

    cout << CLR_LINE;
    cout << "-> Creating Monte-Carlo detector list           \r";

    detector = new CDetector[nr_mc_detectors];

    if(dust_mc_detectors.size() <= 0)
    {
        cout << ERROR_LINE << "No Monte-Carlo detector defined!" << endl;
        delete[] detector;
        detector = 0;
        return detector;
    }

    for(uint i = 0; i < dust_mc_detectors.size(); i += NR_OF_MC_DET)
    {
        uint pos = i / NR_OF_MC_DET;

        double lam_min = dust_mc_detectors[i + 0];
        double lam_max = dust_mc_detectors[i + 1];
        uint nr_spectral_bins = uint(dust_mc_detectors[i + 2]);

        double rot_angle_1 = PI * dust_mc_detectors[i + 3] / 180.0;
        double rot_angle_2 = PI * dust_mc_detectors[i + 4] / 180.0;

        double distance = dust_mc_detectors[i + 5];

        double sideLength_x = dust_mc_detectors[i + 6];
        double sideLength_y = dust_mc_detectors[i + 7];

        double map_shift_x = dust_mc_detectors[i + 8];
        double map_shift_y = dust_mc_detectors[i + 9];

        uint bins_x = uint(dust_mc_detectors[i + NR_OF_MC_DET - 2]);
        uint bins_y = uint(dust_mc_detectors[i + NR_OF_MC_DET - 1]);

        double max_length = grid->getMaxLength();

        if(sideLength_x <= 0)
            sideLength_x = max_length;
        if(sideLength_y <= 0)
            sideLength_y = max_length;

        if(error)
        {
            delete[] detector;
            detector = 0;
            return detector;
        }

        string path_out = path_data;

        #pragma warning(suppress : 6385)
        detector[pos].init(path_out,
                           bins_x,
                           bins_y,
                           sideLength_x,
                           sideLength_y,
                           map_shift_x,
                           map_shift_y,
                           distance,
                           lam_min,
                           lam_max,
                           nr_spectral_bins);
        detector[pos].setOrientation(axis1, axis2, rot_angle_1, rot_angle_2);
        /*if(param.getPeelOff() && param.getNrOfDustPhotons() != 0)
        {
            cout << INFO_LINE << "Peel-off technique disabled for self-scattering of dust grain
        emission!" << endl; param.setPeelOff(false);
        }*/
        if(param.getPeelOff() && param.getAcceptanceAngle() > 1.0)
            cout << INFO_LINE << "Peel-off technique needs no acceptance angle!" << endl;
        else
            detector[pos].setAcceptanceAngle(param.getAcceptanceAngle() * PI / 180.0);
    }

    // cout << "- Creating dust MC detectors    : done" << endl;
    return detector;
}

void CPipeline::createSourceLists(parameters & param, CDustMixture * dust, CGridBasic * grid)
{
    uint nr_ofSources = param.getNrOfSources();

    cout << CLR_LINE;
    cout << "-> Creating source list             \r" << flush;

    if(param.isRaytracingSimulation())
    {
        // Ray tracing simulations are only using background sources!
        if(param.getNrOfDiffuseSources() > 0)
        {
            cout << WARNING_LINE << "Diffuse sources cannot be considered in "
                 << "dust, line, or synchrotron emission!" << endl;
            nr_ofSources--;
        }

        // if(param.getISRFSource())
        // {
        //     cout << WARNING_LINE << "ISRF as radiation source cannot be considered in "
        //          << "dust, line, or synchrotron emission!" << endl;
        //     nr_ofSources--;
        // }

        if(param.getNrOfBackgroundSources() == 0)
        {
            cout << INFO_LINE << "No background source was defined!" << endl;
            cout << "- Default background source initiated." << endl;
            cout << SEP_LINE;

            sources_ray.clear();
            CSourceBasic * tmp_source = new CSourceBackground();
            tmp_source->setParameter(param, grid, dust, 0);
            sources_ray.push_back(tmp_source);

            nr_ofSources++;
        }
        else
        {
            for(uint s = 0; s < param.getBackgroundSources().size(); s += NR_OF_BG_SOURCES)
            {
                cout << "-> Creating background source list             \r" << flush;

                CSourceBasic * tmp_source = new CSourceBackground();
                string path = param.getBackgroundSourceString(s / NR_OF_BG_SOURCES);
                if(path.size() == 0)
                    tmp_source->setParameter(param, grid, dust, s);
                else
                {
                    if(!tmp_source->setParameterFromFile(param, grid, dust, s))
                    {
                        cout << ERROR_LINE << "Background source nr. " << s / NR_OF_BG_SOURCES + 1 << " undefined!"
                             << endl;
                        sources_ray.clear();
                    }
                }
                sources_ray.push_back(tmp_source);
            }
        }

        if(param.getISRFSource())
        {
            if(!dust->getScatteringToRay())
            {
                cout << WARNING_LINE << "Scattered radiation is disabled. ISRF source will be ignored." << endl;
                nr_ofSources--;
            }
            else if(param.getCommand() != CMD_DUST_EMISSION)
            {
                cout << WARNING_LINE << "ISRF source cannot be considered in line or synchrotron emission!" << endl;
            }
            else
            {
                CSourceBasic * tmp_source = new CSourceISRF();

                cout << "-> Creating ISRF source list             \r" << flush;

                if(param.getISRFPath() != "")
                {
                    if(!tmp_source->setParameterFromFile(param, grid, dust, 0))
                    {
                        cout << ERROR_LINE << "Interstellar radiation field undefined! \n" << flush;
                        sources_mc.clear();
                    }
                }
                else
                    tmp_source->setParameter(param, grid, dust, 0);

                sources_mc.push_back(tmp_source);
            }
        }

        if(param.getNrOfPointSources() > 0)
        {
            if(param.getCommand() != CMD_DUST_EMISSION)
            {
                cout << WARNING_LINE << "Stellar sources cannot be considered in line or synchrotron emission!" << endl;
            }
            else
            {
                if(!dust->getScatteringToRay())
                {
                    cout << WARNING_LINE << "Scattered radiation is disabled. Only direct stellar emission is added." << endl;
                }
                for(uint s = 0; s < param.getPointSources().size(); s += NR_OF_POINT_SOURCES)
                {
                    cout << "-> Creating star source list             \r" << flush;
                    string path = param.getPointSourceString(s / NR_OF_POINT_SOURCES);
                    CSourceBasic * tmp_source = new CSourceStar();

                    if(path.size() == 0)
                    {
                        tmp_source->setParameter(param, grid, dust, s);
                        // tmp_source->setNrOfPhotons(1);
                    }
                    else
                    {
                        if(!tmp_source->setParameterFromFile(param, grid, dust, s))
                        {
                            cout << ERROR_LINE << "Star source nr. " << s / NR_OF_POINT_SOURCES + 1
                                 << " undefined!" << endl;
                            sources_mc.clear();
                        }
                    }
                    sources_mc.push_back(tmp_source);
                }
            }
        }

        if(param.getNrOfLaserSources() > 0)
        {
            if(param.getCommand() != CMD_DUST_EMISSION)
            {
                cout << WARNING_LINE << "Laser sources cannot be considered in line or synchrotron emission!" << endl;
            }
            else
            {
                if(!dust->getScatteringToRay())
                {
                    cout << WARNING_LINE << "Scattered radiation is disabled. Only direct laser emission is added." << endl;
                }
                for(uint s = 0; s < param.getLaserSources().size(); s += NR_OF_LASER_SOURCES)
                {
                    cout << "-> Creating laser source list             \r" << flush;
                    CSourceBasic * tmp_source = new CSourceLaser();
                    tmp_source->setParameter(param, grid, dust, s);
                    sources_mc.push_back(tmp_source);
                }
            }
        }

        if(param.getDustSource())
        {
            if(!dust->getScatteringToRay() && grid->isRadiationFieldAvailable())
            {
                nr_ofSources--;
            }
            else if(param.getCommand() != CMD_DUST_EMISSION)
            {
                cout << WARNING_LINE << "Dust source cannot be considered in line or synchrotron emission!" << endl;
            }
            else
            {
                CSourceBasic * tmp_source = new CSourceDust();
                tmp_source->setParameter(param, grid, dust, 0);
                sources_mc.push_back(tmp_source);
            }
        }

        if(param.isGasSpeciesLevelPopMC())
        {
            if(param.getCommand() != CMD_LINE_EMISSION || !param.isGasSpeciesLevelPopMC())
            {
                nr_ofSources--;
            }
            else
            {
                CSourceBasic * tmp_source = new CSourceGas();
                tmp_source->setParameter(param, grid, dust, 0);
                sources_mc.push_back(tmp_source);
            }
        }
    }
    else
    {
        // Monte-Carlo simulations support various sources!
        for(uint s = 0; s < param.getPointSources().size(); s += NR_OF_POINT_SOURCES)
        {
            cout << "-> Creating star source list             \r" << flush;
            string path = param.getPointSourceString(s / NR_OF_POINT_SOURCES);
            CSourceBasic * tmp_source = new CSourceStar();

            if(path.size() == 0)
                tmp_source->setParameter(param, grid, dust, s);
            else
            {
                if(!tmp_source->setParameterFromFile(param, grid, dust, s))
                {
                    cout << ERROR_LINE << "Star source nr. " << s / NR_OF_POINT_SOURCES + 1 << " undefined!"
                         << endl;
                    sources_mc.clear();
                }
            }
            sources_mc.push_back(tmp_source);
        }

        for(uint s = 0; s < param.getDiffuseSources().size(); s += NR_OF_DIFF_SOURCES)
        {
            cout << "-> Creating starfield source list             \r" << flush;

            CSourceBasic * tmp_source = new CSourceStarField();
            string path = param.getDiffuseSourceString(s / NR_OF_DIFF_SOURCES);

            if(path.size() == 0)
                tmp_source->setParameter(param, grid, dust, s);
            else
            {
                if(!tmp_source->setParameterFromFile(param, grid, dust, s))
                {
                    cout << ERROR_LINE << "Sorce Starfield nr. " << s / NR_OF_DIFF_SOURCES + 1 << " undefined! \n"
                         << flush;
                    sources_mc.clear();
                }
            }
            sources_mc.push_back(tmp_source);
        }

        for(uint s = 0; s < param.getLaserSources().size(); s += NR_OF_LASER_SOURCES)
        {
            cout << "-> Creating laser source list             \r" << flush;
            CSourceBasic * tmp_source = new CSourceLaser();
            tmp_source->setParameter(param, grid, dust, s);
            sources_mc.push_back(tmp_source);
        }

        if(!param.getBackgroundSources().empty())
            cout << ERROR_LINE << "Background sources can only be used for raytracing "
                    "simulations!"
                 << endl;

        if(param.getISRFSource())
        {
            CSourceBasic * tmp_source = new CSourceISRF();

            cout << "-> Creating ISRF source list             \r" << flush;

            if(param.getISRFPath() != "")
            {
                if(!tmp_source->setParameterFromFile(param, grid, dust, 0))
                {
                    cout << ERROR_LINE << "Interstellar radiation field undefined! \n" << flush;
                    sources_mc.clear();
                }
            }
            else
                tmp_source->setParameter(param, grid, dust, 0);

            sources_mc.push_back(tmp_source);
        }

        if(param.getDustSource())
        {
            if(param.isTemperatureSimulation())
            {
                cout << ERROR_LINE << "Dust as radiation source cannot be considered in "
                     << "temperature calculations (use RAT to consider dust as a "
                        "separate source)!"
                     << endl;
                nr_ofSources--;
            }
            else
            {
                CSourceBasic * tmp_source = new CSourceDust();
                tmp_source->setParameter(param, grid, dust, 0);
                sources_mc.push_back(tmp_source);
            }
        }
    }

    if((sources_mc.size() + sources_ray.size()) != nr_ofSources)
    {
        cout << ERROR_LINE << "Not enough source(s) were initiated! \n" << flush;

        sources_mc.clear();
        sources_ray.clear();
    }
}

bool CPipeline::writeSources(parameters & param, CGridBasic * grid)
{
    dlist point_source = param.getPointSources();
    dlist diffuse_source = param.getDiffuseSources();
    dlist background_source = param.getBackgroundSources();

    stringstream str_header, str_plot, str_data;
    stringstream str_ps_data1, str_ps_data2;
    stringstream str_ds_data1, str_ds_data2;
    stringstream str_bs_data1, str_bs_data2;

    double max_len, half_length;

    cout << CLR_LINE;
    cout << "-> Plotting source list           \r" << flush;

    str_header.str("");
    str_plot.str("");
    str_data.str("");

    str_header << "reset" << endl;
    str_header << "set title \'Sources\'" << endl;
    str_header << "set ticslevel 0" << endl;
    str_header << "set size ratio -1" << endl;
    str_header << "set view 45,45" << endl;
    str_header << "set xlabel \'x [m]\'" << endl;
    str_header << "set ylabel \'y [m]\'" << endl;
    str_header << "set zlabel \'z [m]\'" << endl;
    str_header << "set border -1 front linetype -1 linewidth 1.000" << endl;
    str_header << "set ticslevel 0" << endl;
    str_header << "set xtics border" << endl;
    str_header << "set ytics border" << endl;
    str_header << "set ztics border" << endl;
    str_header << "set grid" << endl;
    str_header << "set nokey" << endl;
    str_header << "unset colorbox" << endl;

    if(grid != 0)
    {
        max_len = grid->getMaxLength();
        half_length = max_len / 2.0;
    }
    else
    {
        max_len = 1.0;
        half_length = max_len / 2.0;
    }

    str_header << "set xrange[" << -1.01 * half_length << ":" << 1.01 * half_length << "]" << endl;
    str_header << "set yrange[" << -1.01 * half_length << ":" << 1.01 * half_length << "]" << endl;
    str_header << "set zrange[" << -1.01 * half_length << ":" << 1.01 * half_length << "]" << endl;

    str_header << "set palette defined (0 0.5 0 0, 0.5 1 0 0, 2 1 1 0, 3 0.4 0.5 1)" << endl;
    str_header << "set style arrow 1 ls 1 lw 2 lc rgb 0x0000FF" << endl;
    str_header << "set style line 2 pt 7 ps variable lt palette" << endl;
    str_header << "set style line 3 pt 19 ps variable lt palette" << endl;
    str_header << "set cbrange[1000:10000]" << endl;
    str_header << "splot \'-\' w p ls 2 title \"\",\'-\' using 1:2:3:(sprintf(\"star ID %d, "
                  "%d\",$4,$5)) with labels left offset 2 notitle,\'-\' w p ls 3 title "
                  "\"\",\'-\' using 1:2:3:(sprintf(\"diff ID %d, %d\",$4,$5)) with labels left "
                  "offset 2 notitle,\'-\' with vectors as 1 title \"\",\'-\' using "
                  "1:2:3:(sprintf(\"bs ID %d\",$4)) with labels left offset 2 notitle"
               << endl;

    for(uint s = 0; s < point_source.size(); s += NR_OF_POINT_SOURCES)
    {
        uint i = s / NR_OF_POINT_SOURCES;
        double x = point_source[s + 0];
        double y = point_source[s + 1];
        double z = point_source[s + 2];

        double R = point_source[s + 3];
        double T = point_source[s + 4];
        ullong nr_of_photons = uint(point_source[s + 5]);

        if(R < 1)
            R = 1;
        if(R > 15)
            R = 15;

        if(T < 1000)
            T = 1000;
        if(T > 10000)
            T = 10000;

        str_ps_data1 << x << "\t" << y << "\t" << z << "\t" << R << "\t" << T << endl;
        str_ps_data2 << x << "\t" << y << "\t" << z << "\t" << i + 1 << "\t" << nr_of_photons << endl;
    }

    for(uint s = 0; s < diffuse_source.size(); s += 7)
    {
        uint i = s / 7;
        double R = diffuse_source[s + 3];
        double T = diffuse_source[s + 4];
        ullong nr_of_photons = uint(diffuse_source[s + 6]);

        if(R < 1)
            R = 1;
        if(R > 15)
            R = 10;

        if(T < 1000)
            T = 1000;
        if(T > 10000)
            T = 10000;

        str_ds_data1 << diffuse_source[s] << "\t" << diffuse_source[s + 1] << "\t" << diffuse_source[s + 2]
                     << "\t" << R << "\t" << T << endl;
        str_ds_data2 << diffuse_source[s] << "\t" << diffuse_source[s + 1] << "\t" << diffuse_source[s + 2]
                     << "\t" << i + 1 << "\t" << nr_of_photons << endl;
    }

    // for (uint s = 0; s < background_source.size(); s += 7) {
    //	uint i = s / 7;
    //}

    if(param.getDustSource())
    {
    }

    // path_plot = "E:\\gnutests\\";
    string plot_out = path_plot + "sources.plt";
    ofstream outStream(plot_out.c_str());

    if(outStream.fail())
    {
        cout << ERROR_LINE << "Can plot sources to:" << endl;
        cout << plot_out << endl;
        return false;
    }

    str_ps_data1 << "\ne" << endl;
    str_ps_data2 << "\ne" << endl;
    str_ds_data1 << "\ne" << endl;
    str_ds_data2 << "\ne" << endl;

    outStream << str_header.str() << endl;
    outStream << str_ps_data1.str() << endl;
    outStream << str_ps_data2.str() << endl;
    outStream << str_ds_data1.str() << endl;
    outStream << str_ds_data2.str() << endl;
    outStream << "\ne\ne" << endl;

    outStream.close();

    // cout << "- Plotting of sources           : done          \n" << flush;

    return true;
}

bool CPipeline::assignDustMixture(parameters & param, CDustMixture * dust, CGridBasic * grid)
{
    // Get the number of different mixtures of dust_components
    if(!dust->createDustMixtures(param, path_data, path_plot))
        return false;

    // Write plot files to show dust properties
    if(param.getWriteDustFiles())
        if(!dust->writeComponentPlot(path_plot))
            return false;
    if(!dust->writeComponentData(path_data))
        return false;

    dust->setGridRequirements(grid, param);

    // Check if either the radiation field is present or the radiation field can be
    // calculated otherwise, disable scattering added to the raytracing then
    if(param.getScatteringToRay())
        !grid->isRadiationFieldAvailable() && param.getNrOfPointSources() == 0 &&
                param.getNrOfDiffuseSources() == 0 && param.getNrOfLaserSources() == 0 &&
                !param.getISRFSource() && !param.getDustSource()
            ? dust->setScatteringToRay(false)
            : dust->setScatteringToRay(true);
    else
        dust->setScatteringToRay(false);
    return true;
}

bool CPipeline::assignGasSpecies(parameters & param, CGasMixture * gas, CGridBasic * grid)
{
    if(!gas->createGasSpecies(param))
        return false;

    return true;
}

void CPipeline::printParameters(parameters & param, uint max_id)
{
    cout << CLR_LINE;
    cout << "Input parameters (task " << param.getTaskID() << " of " << max_id << ")" << endl;
    cout << SEP_LINE;

    switch(param.getCommand())
    {
        case CMD_TEMP:
            cout << "- Command          : TEMPERATURE DISTRIBUTION" << endl;
            printPathParameters(param);
            printSourceParameters(param);
            printConversionParameters(param);
            printAdditionalParameters(param);
            printPlotParameters(param);
            break;

        case CMD_FORCE:
            cout << "- Command          : RADIATION FORCE" << endl;
            printPathParameters(param);
            printSourceParameters(param);
            printConversionParameters(param);
            printAdditionalParameters(param);
            printPlotParameters(param);
            break;

        case CMD_TEMP_RAT:
            cout << "- Command          : TEMPERATURE DISTRIBUTION and RAT ALIGNMENT" << endl;
            printPathParameters(param);
            printSourceParameters(param);
            printConversionParameters(param);
            printAdditionalParameters(param);
            printPlotParameters(param);
            break;

        case CMD_RAT:
            cout << "- Command          : RAT ALIGNMENT" << endl;
            printPathParameters(param);
            printSourceParameters(param, true);
            printConversionParameters(param);
            printAdditionalParameters(param);
            printPlotParameters(param);
            break;

        case CMD_DUST_EMISSION:
            cout << "- Command          : DUST EMISSION" << endl;
            printPathParameters(param);
            printSourceParameters(param, true);
            printConversionParameters(param);
            printAdditionalParameters(param);
            printAlignmentParameters(param);
            printDetectorParameters(param);
            printPlotParameters(param);
            break;

        case CMD_SYNCHROTRON:
            cout << "- Command          : SYNCHROTRON EMISSION" << endl;
            printPathParameters(param);
            printSourceParameters(param);
            printConversionParameters(param);
            printAdditionalParameters(param);
            printDetectorParameters(param);
            printPlotParameters(param);
            printSynchrotronParameters(param);
            break;

        case CMD_DUST_SCATTERING:
            cout << "- Command          : DUST SCATTERING (Monte-Carlo)" << endl;
            printPathParameters(param);
            printSourceParameters(param, true);
            printConversionParameters(param);
            printAdditionalParameters(param);
            printAlignmentParameters(param);
            printDetectorParameters(param, true);
            printPlotParameters(param);
            break;

        case CMD_OPIATE:
            cout << "todo: OPIATE parameter" << endl;

            break;
        case CMD_LINE_EMISSION:
            cout << "- Command          : SPECTRAL LINE EMISSION" << endl;
            printPathParameters(param);
            printSourceParameters(param);
            printConversionParameters(param);
            printAdditionalParameters(param);
            printDetectorParameters(param);
            printPlotParameters(param);
            break;

        default:
            cout << ERROR_LINE << "Command is unknown!" << endl;
            cout << "No parameters available! " << endl;
    };
    cout << SEP_LINE;
}

bool CPipeline::createWavelengthList(parameters & param, CDustMixture * dust, CGasMixture * gas, COpiateDataBase * op)
{
    dlist values;

    switch(param.getCommand())
    {
        case CMD_TEMP:
        case CMD_TEMP_RAT:
        case CMD_RAT:
            dust->addToWavelengthGrid(WL_MIN, WL_MAX, WL_STEPS);
            break;

        case CMD_DUST_EMISSION:
            // Add wavelength for stochastic heating
            if(param.getStochasticHeatingMaxSize() > 0)
                dust->addToWavelengthGrid(WL_MIN, WL_MAX, WL_STEPS, true);

            // Get detector parameters list
            values = param.getDustRayDetectors();

            // Check if a detector is defined
            if(values.empty())
            {
                cout << ERROR_LINE << "No dust raytracing detector defined (see "
                        "<detector_dust>)!"
                     << endl;
                return false;
            }

            // Add wavelength to global list of wavelength
            for(uint i = 0; i < values.size(); i += NR_OF_RAY_DET)
                dust->addToWavelengthGrid(values[i], values[i + 1], values[i + 2]);
            break;

        case CMD_DUST_SCATTERING:
            // Get detector parameters list
            values = param.getDustMCDetectors();

            // Check if a detector is defined
            if(values.empty())
            {
                cout << ERROR_LINE << "No dust Monte-Carlo detector defined (see "
                        "<detector_dust_mc>)!"
                     << endl;
                return false;
            }

            // Add wavelength to global list of wavelength
            for(uint i = 0; i < values.size(); i += NR_OF_MC_DET)
                dust->addToWavelengthGrid(values[i], values[i + 1], values[i + 2]);
            break;

        case CMD_LINE_EMISSION:
        {
            if(gas == 0)
                return false;

            if(param.isGasSpeciesLevelPopMC())
            {
                for(uint i_species = 0; i_species < gas->getNrOfSpecies(); i_species++)
                {
                    for(uint i_trans = 0; i_trans < gas->getNrOfTransitions(i_species); i_trans++)
                    {
                        // Calculate from frequency
                        double wavelength = con_c / gas->getTransitionFrequency(i_species, i_trans);

                        // Add wavelength to global list of wavelength
                        dust->addToWavelengthGrid(wavelength);
                    }
                }
            }
            else
            {
                // Get detector parameters list
                maplist line_ray_detector_list = param.getLineRayDetectors();
                maplist::iterator it;

                // Check if a detector is defined
                if(line_ray_detector_list.empty())
                {
                    cout << ERROR_LINE << "No spectral line detector of gas species defined!" << endl;
                    return false;
                }

                // Perform radiative transfer for each chosen gas species
                for(it = line_ray_detector_list.begin(); it != line_ray_detector_list.end(); ++it)
                {
                    // Get ID of the current gas species
                    uint i_species = it->first;

                    // Get number of spectral line transitions that have to be simulated
                    uint nr_of_spectral_lines = param.getNrOfSpectralLines(i_species);

                    // Perform radiative transfer for each chosen spectral line transition
                    for(uint i_line = 0; i_line < nr_of_spectral_lines; i_line++)
                    {
                        // Calculate from frequency
                        double wavelength = con_c / gas->getSpectralLineFrequency(i_species, i_line);

                        // Add wavelength to global list of wavelength
                        dust->addToWavelengthGrid(wavelength);
                    }
                }
            }
            break;
        }

        case CMD_OPIATE:

            if(op==0)
            {
                cout << ERROR_LINE << "No OPIATE database loaded!" << endl;
                return false;
            }

            // Get detector parameters list
            values = param.getOPIATERayDetectors();

            // Check if a detector is defined
            if(values.empty())
            {
                cout << ERROR_LINE << "No OPIATE detector defined (see <detector_opiate>)!" << endl;
                return false;
            }

            for(uint i=0;i<param.getNrOfOPIATESpecies();i++)
            {
                string spec_name=param.getOpiateSpec(i);
                if(op->findIndexByName(spec_name))
                {
                    double wavelength = con_c / op->getCurrentFrequency();
                    dust->addToWavelengthGrid(wavelength);
                }
                else
                {
                    cout << ERROR_LINE << "Label \"" << spec_name << "\" is not listed in OPIATE database! " << endl;
                    cout << "         Check the command \"<detector_opiate>\" in command file!" << endl;
                    return false;
                }
            }

            break;

        case CMD_SYNCHROTRON:
            // Get detector parameters list
            values = param.getSyncRayDetectors();

            // Check if a detector is defined
            if(values.empty())
            {
                cout << ERROR_LINE << "No synchrotron detector defined (see <detector_sync>)!" << endl;
                return false;
            }

            // Add wavelength to global list of wavelength
            for(uint i = 0; i < values.size(); i += NR_OF_RAY_DET)
                dust->addToWavelengthGrid(values[i], values[i + 1], values[i + 2]);
            break;

        default:
            break;
    }

    // Discard wavelengths that are duplicates
    dust->finalizeWavelengthList();

    return true;
}

bool CPipeline::solveNextVelocityField(uint itID)
{
    return true;
}

bool CPipeline::createOutputPaths(string path)
{
    const char * sep = SEP;
    string::size_type pos1 = 0, pos2 = 0;
    int len = -1;
    uint offset;

#ifdef WINDOWS
    offset = 0;
#else
    offset = 1;
#endif

    path_data = path + "data" + sep;
    path_plot = path + "plots" + sep;

    strlist folders;
    string folder;

    while(path.find(sep, offset) != string::npos)
    {
        pos1 = path.find(sep, offset);
        pos2 = path.find(sep, pos1 + 1);

        len = int(pos2 - pos1 + 1);

        if(len < 0)
            break;

        folder = path.substr(pos1 + 1, len - 2);
        path.erase(pos1 + 1, len - 1);
        folders.push_back(folder);
    }

    for(uint i = 0; i < folders.size(); i++)
    {
        folder = folders[i];
        path += folder + sep;

        if(!createPath(path))
        {
            cout << ERROR_LINE << "Failed to create output folder for data!" << endl;
            cout << path << std::endl;
            return false;
        }
    }

    if(!createPath(path))
    {
        cout << ERROR_LINE << "Failed to create output folder(s)!" << endl;
        cout << path << std::endl;
        return false;
    }

    if(!createPath(path_data))
    {
        cout << ERROR_LINE << "Failed to create output folder for data!" << endl;
        cout << path_data << std::endl;
        return false;
    }

    if(!createPath(path_plot))
    {
        cout << ERROR_LINE << "Failed to create output folder for plots!" << endl;
        cout << path_plot << std::endl;
        return false;
    }

    return true;
}

bool CPipeline::createPath(string path)
{
    const char * tmp_path = path.c_str();
#ifdef WINDOWS
    if(_access(tmp_path, 0) != 0)
        if(_mkdir(tmp_path))
            return false;
#else
    if(access(tmp_path, 0) != 0)
        if(mkdir(tmp_path, S_IRWXU | S_IRWXG | S_IRWXO))
            return false;
#endif

    return true;
}

void CPipeline::printConversionParameters(parameters & param)
{
    cout << "Conversion factors" << endl;
    cout << "- Conv. gas dens. in SI        : " << param.getSIConvDH() << endl;
    cout << "- Conv. length in SI           : " << param.getSIConvLength() << endl;
    cout << "- Conv. vel. field in SI.      : " << param.getSIConvVField() << endl;
    cout << "- Conv. mag. field in SI.      : " << param.getSIConvBField() << endl;
    if(param.getIndividualDustMassFractions())
        cout << "- Mass fraction (Mdust/Mgas)   : set by the dust components" << endl;
    else
        cout << "- Mass fraction (Mdust/Mgas)   : " << param.getDustMassFraction() << endl;
    cout << "- Relative molecular mass (mu) : " << param.getMu() << endl;
}

void CPipeline::printAdditionalParameters(parameters & param)
{
    cout << "Additional parameters" << endl;
    if(param.isRatSimulation())
    {
        cout << "- RAT efficiency ref. (Qref)   : " << param.getQref() << endl;
        cout << "- RAT efficiency exp. (alphaQ) : " << param.getAlphaQ() << endl;
    }
    else
        cout << "- None" << endl;
}

void CPipeline::printPathParameters(parameters & param)
{
    cout << "- Path grid file   : " << param.getPathGrid() << endl;
    cout << "- Path output      : " << param.getPathOutput() << endl;
    cout << "- Number of threads: " << param.getNrOfThreads() << endl;
}

void CPipeline::printPlotParameters(parameters & param, bool input_output)
{
    if(param.getNrOfPlotPoints() + param.getNrOfPlotVectors() + param.getInpMidDataPoints() +
            param.getOutMidDataPoints() + param.getInpAMIRAPoints() + param.getOutAMIRAPoints() >
        0)
        cout << "Plot parameters" << endl;

    if(param.getNrOfPlotPoints() + param.getNrOfPlotVectors() != 0)
    {
        cout << "- Raw data                     : " << param.getNrOfPlotPoints() << " points, ";
        cout << param.getNrOfPlotVectors() << " vectors, ";
        cout << param.getMaxPlotLines() << " lines" << endl;
    }

    if(param.getInpMidDataPoints() != 0 && param.getOutMidDataPoints() != 0)
    {
        cout << "- Midplane            (in,out) : " << param.getInpMidDataPoints() << ", "
                << param.getOutMidDataPoints() << " pixel" << endl;
    }
    else
    {
        if(param.getInpMidDataPoints() != 0)
            cout << "- Midplane                (in) : " << param.getInpMidDataPoints() << " pixel"
                    << endl;

        if(param.getOutMidDataPoints() != 0)
            cout << "- Midplane               (out) : " << param.getOutMidDataPoints() << " pixel"
                    << endl;
    }

    if(param.getInpAMIRAPoints() != 0 && param.getOutAMIRAPoints() != 0)
        cout << "- Amira points (in,out)        : " << param.getInpAMIRAPoints() << ", "
                << param.getOutAMIRAPoints() << " pixel" << endl;
    else
    {
        if(param.getInpAMIRAPoints() != 0)
            cout << "- Amira points            (in) : " << param.getInpAMIRAPoints() << " pixel" << endl;
        if(param.getOutAMIRAPoints() != 0)
            cout << "- Amira points           (out) : " << param.getOutAMIRAPoints() << " pixel" << endl;
    }
}

void CPipeline::printAlignmentParameters(parameters & param)
{
    cout << "Dust grain alignment" << endl;

    if(param.getAligRANDOM())
        cout << "- No alignment" << endl;
    else if(param.getAligPA())
        cout << "- Perfect alignment" << endl;
    else if(param.getAligNONPA())
    {
        cout << "- Non-perfect alignment" << endl;
        if(param.getRayleighReductionFactor() == -1)
            cout << "    Rayleigh reduction fact. (R) : Based on absolute B-field in grid cells" << endl;
        else
            cout << "    Rayleigh reduction fact. (R) : " << param.getRayleighReductionFactor() << endl;
    }
    else
    {
        if(param.getAligIDG())
            cout << "- IDG (paramagnetic alignment)" << endl;
        if(param.getAligRAT())
        {
            cout << "- RAT (radiative torques)" << endl;
            if(!param.getAligINTERNAL())
                cout << "    Rayleigh reduction fact. (R) : " << param.getRayleighReductionFactor() << endl;
        }
        if(param.getAligGOLD())
            cout << "- GOLD (mechanical alignment)" << endl;
        if(param.getAligINTERNAL())
        {
            cout << "- Internal alignment" << endl;
            if(param.getAligIDG() || param.getAligGOLD())
                cout << "    Correlation factor (f_c)     : " << param.getFcorr() << endl;
            if(param.getAligRAT())
                cout << "    Fraction at high-J (f_highJ) : " << param.getFHighJ() << endl;
        }
        cout << "    Align limit prefact. (larmF) : " << param.getLarmF() << endl;
    }
}

void CPipeline::printSourceParameters(parameters & param, bool show_dust)
{
    if(param.getNrOfSources() > 0 || param.isRaytracingSimulation())
    {
        cout << "Defined radiation source(s):" << endl;

        if(param.getNrOfPointSources() > 0 &&
            (param.isMonteCarloSimulation() || param.getCommand() == CMD_DUST_EMISSION || param.getScatteringToRay()))
        {
            dlist sources_list = param.getPointSources();
            for(uint s = 0; s < sources_list.size(); s += NR_OF_POINT_SOURCES)
            {
                // Index of current source
                uint pos = s / NR_OF_POINT_SOURCES;
                cout << "- Star " << pos + 1 << " :\n"
                        << "    Position    : "
                        << sources_list[s + 0] << ", "
                        << sources_list[s + 1] << ", "
                        << sources_list[s + 2] << " (x,y,z) [m]\n"
                        << "    Radius      : " << sources_list[s + 3] << " [R_sun]\n"
                        << "    Temperature : " << sources_list[s + 4] << " [K]\n"
                        << "    Stokes      : "
                        << sources_list[s + 5] << ", "
                        << sources_list[s + 6] << " (q,u)" << endl;
            }
        }
            // cout << "- Star(s)        : " << param.getNrOfPointSources() << endl;

        if(param.getNrOfDiffuseSources() > 0 && (!param.isRaytracingSimulation() || param.getScatteringToRay()))
        {
            dlist sources_list = param.getDiffuseSources();
            for(uint s = 0; s < sources_list.size(); s += NR_OF_DIFF_SOURCES)
            {
                // Index of current source
                uint pos = s / NR_OF_DIFF_SOURCES;
                cout << "- Starfield " << pos + 1 << " :\n"
                        << "    Position    : "
                        << sources_list[s + 0] << ", "
                        << sources_list[s + 1] << ", "
                        << sources_list[s + 2] << " (x,y,z) [m]\n"
                        << "    Radius      : " << sources_list[s + 3] << " [R_sun]\n"
                        << "    Temperature : " << sources_list[s + 4] << " [K]\n"
                        << "    Variance    : " << sources_list[s + 5] << " [m]\n"
                        << "    Stokes      : "
                        << sources_list[s + 6] << ", "
                        << sources_list[s + 7] << " (q,u)" << endl;
            }
        }
            // cout << "- Star field(s)  : " << param.getNrOfDiffuseSources() << endl;

        if(param.getNrOfLaserSources() > 0 && (!param.isRaytracingSimulation() || param.getScatteringToRay()))
        {
            dlist sources_list = param.getLaserSources();
            for(uint s = 0; s < sources_list.size(); s += NR_OF_LASER_SOURCES)
            {
                // Index of current source
                uint pos = s / NR_OF_LASER_SOURCES;
                cout << "- Laser " << pos + 1 << " :\n"
                        << "    Position           : "
                        << sources_list[s + 0] << ", "
                        << sources_list[s + 1] << ", "
                        << sources_list[s + 2] << " (x,y,z) [m]\n"
                        << "    Direction          : "
                        << sources_list[s + 3] << ", "
                        << sources_list[s + 4] << ", "
                        << sources_list[s + 5] << " (x,y,z) [m]\n"
                        << "    Total power        : " << sources_list[s + 6] << " [W]\n"
                        << "    Central wavelength : " << sources_list[s + 7] << " [m]\n"
                        << "    FWHM               : " << sources_list[s + 8] << " [m]\n"
                        << "    Stokes             : "
                        << sources_list[s + 9] << ", "
                        << sources_list[s + 10] << " (q,u)" << endl;
            }
        }
            // cout << "- Laser(s)  : " << param.getNrOfLaserSources() << endl;

        if(param.getNrOfBackgroundSources() > 0 || param.isRaytracingSimulation())
        {
            dlist sources_list = param.getBackgroundSources();
            for(uint s = 0; s < sources_list.size(); s += NR_OF_BG_SOURCES)
            {
                // Index of current source
                uint pos = s / NR_OF_BG_SOURCES;
                cout << "- Background " << pos + 1 << " :\n"
                        << "    Scaling factor  : " << sources_list[s + 0] << "\n"
                        << "    Temperature     : " << sources_list[s + 1] << " [K]\n"
                        << "    Stokes          : "
                        << sources_list[s + 2] << ", "
                        << sources_list[s + 3] << ", "
                        << sources_list[s + 4] << " (q,u,v)\n"
                        << "    Rotation angles : "
                        << sources_list[s + 5] << ", "
                        << sources_list[s + 6] << " (1,2) [°]" << endl;
            }
        }
            // cout << "- Background(s)  : " << max(uint(1), param.getNrOfBackgroundSources()) << endl;

        if(param.getDustSource() && show_dust)
            cout << "- Dust as sources: yes" << endl;

        if(param.getISRFSource() && (!param.isRaytracingSimulation() || param.getScatteringToRay()))
            cout << "- ISRF as sources: yes" << endl;
    }
}

void CPipeline::printDetectorParameters(parameters & param, bool monte_carlo)
{
    if(monte_carlo)
    {
        cout << "Monte-Carlo parameter" << endl;
        if(param.getNrOfDustMCDetectors() > 0)
            cout << "- Number of detectors   : " << param.getNrOfDustMCDetectors() << endl;
    }
    else
    {
        cout << "Raytrace parameter" << endl;
        if(param.getNrOfDustRayDetectors() > 0)
            cout << "- Number of detectors   : " << param.getNrOfDustRayDetectors() << endl;
        param.printRTGridDescription();
        cout << "- Start, Stop           : " << param.getStart() + 1 << ", " << param.getStop() + 1
                << endl;
    }
    cout << "- Observer distance [m] : " << param.getMinObserverDistance() << " - "
            << param.getMaxObserverDistance() << endl;
    cout << "- Pixel in x-dir/nsides : " << param.getMinDetectorPixelX() << " (min) - "
            << param.getMaxDetectorPixelX() << " (max)" << endl;
    if(param.getMinDetectorPixelY() != 0 || param.getMaxDetectorPixelY() != 0)
        cout << "- Pixel in y-direction  : " << param.getMinDetectorPixelY() << " (min) - "
                << param.getMaxDetectorPixelY() << " (max)" << endl;
    if(!monte_carlo)
        cout << "- Max subpixel level    : " << param.getMaxSubpixelLvl() << endl;

    if(param.getMaxSidelengthX() > 0)
    {
        cout << "- Sidelength X [m]      : " << param.getMinSidelengthX() << " (min) - "
                << param.getMaxSidelengthX() << " (max)" << endl;
        if(param.getUseGridSidelengthX())
            cout << "  (some detectors are using the grid size)" << endl;
    }
    else if(param.getUseGridSidelengthX())
        cout << "- Sidelength X [m]      : only the grid size is used" << endl;
    if(param.getMaxSidelengthY() > 0)
    {
        cout << "- Sidelength Y [m]      : " << param.getMinSidelengthY() << " (min) - "
                << param.getMaxSidelengthY() << " (max)" << endl;
        if(param.getUseGridSidelengthY())
            cout << "  (some detectors are using the grid size)" << endl;
    }
    else if(param.getUseGridSidelengthY())
        cout << "- Sidelength Y [m]      : only the grid size is used" << endl;

    if(!monte_carlo)
    {
        if(param.getMaxMapShiftX() != 0)
            cout << "- Map shift X [m]       : " << param.getMinMapShiftX() << " (min) - "
                    << param.getMaxMapShiftX() << " (max)" << endl;
        if(param.getMaxMapShiftY() != 0)
            cout << "- Map shift Y [m]       : " << param.getMinMapShiftY() << " (min) - "
                    << param.getMaxMapShiftY() << " (max)" << endl;
    }

    cout << "- Rotation axis         : n1 = ";
    param.getAxis1().printValues();
    cout << ", angle [°] = " << param.getMinDetectorAngle1() << " (min) - "
            << param.getMaxDetectorAngle1() << " (max)" << endl;
    cout << "- Rotation axis         : n2 = ";
    param.getAxis2().printValues();
    cout << ", angle [°] = " << param.getMinDetectorAngle2() << " (min) - "
            << param.getMaxDetectorAngle2() << " (max)" << endl;
}

void CPipeline::printSynchrotronParameters(parameters & param)
{
    cout << SEP_LINE;
    cout << "Synchotron parameter" << endl;
    cout << SEP_LINE;
    cout << "Observed wavelengths" << endl;
    dlist sync_ray_detectors = param.getSyncRayDetectors();
    for(uint i = 0; i < sync_ray_detectors.size(); i += NR_OF_RAY_DET)
    {
        uint pos = i / NR_OF_RAY_DET;
        double lam_min = sync_ray_detectors[i + 0];
        double lam_max = sync_ray_detectors[i + 1];
        uint lam_skip = uint(sync_ray_detectors[i + 2]);
        cout << "    - Synchotron emission detector " << (pos + 1) << ": from wl = " << lam_min
                << " [m] to wl = " << lam_max << " [m] with " << lam_skip << " step(s)" << endl;
    }
}

void CPipeline::deleteSourceLists()
{
    for(uint s = 0; s < sources_mc.size(); s++)
        delete sources_mc[s];
    sources_mc.clear();

    for(uint s = 0; s < sources_ray.size(); s++)
        delete sources_ray[s];
    sources_ray.clear();
}

uint CPipeline::getNrOffsetEntriesRay(parameters & param, CDustMixture * dust, CGridBasic * grid)
{
    uint nr_of_offset_entries = 0;
    if(param.getStochasticHeatingMaxSize() > 0)
    {
        // Add fields to store the stochastic heating propabilities for each cell
        for(uint i_mixture = 0; i_mixture < dust->getNrOfMixtures(); i_mixture++)
            nr_of_offset_entries +=
                dust->getNrOfStochasticSizes(i_mixture) * dust->getNrOfCalorimetryTemperatures(i_mixture);
    }
    else if(dust->getScatteringToRay())
    {
        // Set that the spectral length is saved as vector
        grid->setSpecLengthAsVector(true);

        // Get list of dust raytracing detector parameters
        dlist dust_ray_detectors = param.getDustRayDetectors();

        // Get number of dust raytracing detectors
        uint nr_ray_detectors = uint(dust_ray_detectors.size()) / NR_OF_RAY_DET;

        // Check if one of the detectors is using the healpix background grid
        bool detector_uses_healpix = false;
        for(uint i = 0; i < nr_ray_detectors; i++)
        {
            uint pos = i * NR_OF_RAY_DET;
            uint detector_id = uint(dust_ray_detectors[pos + NR_OF_RAY_DET - 3]);

            if(detector_id == DET_SPHER)
                detector_uses_healpix = true;
        }

        // Set up a list of all detector directions to use for the scattering via radiation field
        if(!detector_uses_healpix)
        {
            for(uint i_det = 0; i_det < nr_ray_detectors; i_det++)
            {
                uint nr_spectral_bins = uint(dust_ray_detectors[i_det * NR_OF_RAY_DET + 2]);

                // Add fields for the radiation field of each considered wavelength and detector
                nr_of_offset_entries += 4 * nr_spectral_bins;
            }
        }
        else
        {
            nr_of_offset_entries += 4 * dust->getNrOfWavelength();
        }
    }

    return nr_of_offset_entries;
}

double CPipeline::getNrOfRayDetector(parameters & param)
{
    dlist dust_ray_detectors = param.getDustRayDetectors();
    return uint(dust_ray_detectors.size()) / NR_OF_RAY_DET;
}

/*
bool CPipeline::calcRadPressure(parameter & param)
{
    CGridBasic * grid = 0;
    CDustMixture * dust = new CDustMixture();

    if(!createOutputPaths(param.getPathOutput()))
        return false;

    if(!assignGridType(grid, param))
        return false;

    if(!assignDustMixture(param, dust, grid))
        return false;

    grid->setSIConversionFactors(param);

    if(!grid->loadGridFromBinaryFile(param, 3 + 3 * WL_STEPS))
        return false;

    // Print helpfull information
    grid->createCellList();
    dust->printParameter(param, grid);
    grid->printParameter();

    if(!grid->writeMidplaneFits(path_data + "input_", param, param.getInpMidDataPoints(),
true)) return false;

    if(!grid->writePlotFiles(path_plot + "input_", param))
        return false;

    if(!grid->writeAMIRAFiles(path_plot + "input_", param, param.getInpAMIRAPoints()))
        return false;

    createSourceLists(param, dust, grid);
    if(sources_mc.size() == 0)
    {
        cout << ERROR_LINE << "No sources for Monte-Carlo simulations defined!" << endl;
        return false;
    }

    CRadiativeTransfer rad(param);

    rad.setGrid(grid);
    rad.setDust(dust);
    rad.setSourcesLists(sources_mc, sources_ray);
    rad.initiateRadFieldMC(param);

    //rad.calcRadiativePressure(param);
    preparePressureData(grid, dust, param, true, 0);
    //solveNextVelocityField(itID);

    cout << SEP_LINE;

    delete grid;
    delete dust;
    deleteSourceLists();

    return true;
}
*/

/*
bool CPipeline::preparePressureData(CGridBasic * grid, CDustMixture * dust, parameters &
param, bool plot, uint itID)
{
    ullong per_counter = 0;
    ulong max_cells = grid->getMaxDataCells();
    double mu = param.getMu();
    ulong data_off = grid->getDataOffset();
    bool res = true;
    uint nr_of_wavelength = dust->getNrOfWavelength();
    double Mstar = param.getStarMass(0);

    #pragma omp parallel for schedule(dynamic)
    for(long c_1 = 0; c_1 < long(max_cells - 1); c_1++)
    {
        cell_basic * cell_1 = (cell_basic*) grid->getCellFromIndex(c_1);
        double Md = dust->getAvgMass(grid, cell_1);
        double dens_1 = grid->getGasNumberDensity(cell_1);
        double vol_1 = grid->getVolume(cell_1);
        double M1 = dens_1 * vol_1 * mu * m_H;
        Vector3D pos_1;

        if(dens_1 == 0)
            continue;

        pos_1 = grid->getCenter(cell_1);

        Vector3D F1 = con_G * Mstar * Md / (pos_1.length() * pos_1.length()) *
pos_1.normalized();

        cell_1->updateData(data_off + 0, -F1.X());
        cell_1->updateData(data_off + 1, -F1.Y());
        cell_1->updateData(data_off + 2, -F1.Z());

        for(long c_2 = c_1 + 1; c_2 < long(max_cells - 1); c_2++)
        {
        cell_oc * cell_2 = (cell_oc*) grid->getCellFromIndex(c_2);
        double dens_2 = grid->getGasNumberDensity(cell_2);
        double vol_2 = grid->getVolume(cell_2);
        double M2 = dens_2 * vol_2 * mu * m_H;
        Vector3D pos_2;

        if(dens_2 ==0)
        continue;

        pos_2.setX(cell_2->x_min + (cell_2->x_max - cell_2->x_min) / 2.0);
        pos_2.setY(cell_2->y_min + (cell_2->y_max - cell_2->y_min) / 2.0);
        pos_2.setZ(cell_2->z_min + (cell_2->z_max - cell_2->z_min) / 2.0);

        Vector3D dist = pos_2 - pos_1;
        double R = dist.length();
        dist.normalize();

        Vector3D F1 = -con_G * M2*Md / (R * R) * dist;
        Vector3D F2 = -con_G * M1*Md / (R * R) * dist; // -con_G*M1*Mstar /
(pos_2.sq_length())*pos_2.normalized();

        cell_1->updateData(data_off + 0, -F1.X());
        cell_1->updateData(data_off + 1, -F1.Y());
        cell_1->updateData(data_off + 2, -F1.Z());

        cell_2->updateData(data_off + 0, F2.X());
        cell_2->updateData(data_off + 1, F2.Y());
        cell_2->updateData(data_off + 2, F2.Z());
        }

        double FradX = 0, FradY = 0, FradZ = 0;
        double * tmpX = new double[nr_of_wavelength];
        double * tmpY = new double[nr_of_wavelength];
        double * tmpZ = new double[nr_of_wavelength];

        for(uint w = 0; w < nr_of_wavelength; w++)
        {
            tmpX[w] = cell_1->getData(data_off + 3 + w * 3 + 0);
            tmpY[w] = cell_1->getData(data_off + 3 + w * 3 + 1);
            tmpZ[w] = cell_1->getData(data_off + 3 + w * 3 + 2);

            if(w > 0)
            {
                FradX += (wavelength_list[w] - wavelength_list[w - 1]) * tmpX[w - 1] + 0.5
* (wavelength_list[w] - wavelength_list[w - 1]) * (tmpX[w] - tmpX[w - 1]); FradY +=
(wavelength_list[w] - wavelength_list[w - 1]) * tmpY[w - 1] + 0.5 * (wavelength_list[w] -
wavelength_list[w - 1]) * (tmpY[w] - tmpY[w - 1]); FradZ += (wavelength_list[w] -
wavelength_list[w - 1]) * tmpZ[w - 1] + 0.5 * (wavelength_list[w] - wavelength_list[w -
1]) * (tmpZ[w] - tmpZ[w - 1]);
            }
        }

        delete[] tmpX;
        delete[] tmpY;
        delete[] tmpZ;

        cell_1->setData(data_off + 3 + 0, FradX);
        cell_1->setData(data_off + 3 + 1, FradY);
        cell_1->setData(data_off + 3 + 2, FradZ);
    }

    cout << CLR_LINE;
    cout << "- Calculation of final properties: done" << endl;

    uint bins = param.getInpMidDataPoints();
    double max_len = grid->getMaxLength();

    per_counter = 0;

    uint per_max = 3 * bins * bins;
    string prev[3] = {"_xy", "_xz", "_yz"};

    int b_limit = int(bins) / 2;
    double xyz_step = max_len / double(bins);

    double * v_pos = new double[100];
    double * v_Frad = new double[100];
    double * v_Fgra = new double[100];

    for(uint i = 0; i < 100; i++)
    {
        v_pos[i] = 0;
        v_Frad[i] = 0;
        v_Fgra[i] = 0;
    }

    double xy_min = param.getXYMin();
    double xy_max = param.getXYMax();
    double xy_steps = param.getXYSteps();
    string xy_label = param.getXYLabel();
    uchar vec_color[3];
    vec_color[0] = 255;
    vec_color[1] = 255;
    vec_color[2] = 255;

    double off_xyz = 0.5 * xyz_step;

    #pragma omp parallel for schedule(dynamic)
    for(int i = 1; i <= 3; i++)
    {
        photon_package pp = photon_package();

        string rad_filename = path_plot + "Frad" + prev[i - 1] + ".dat";
        string gra_filename = path_plot + "Fgrav" + prev[i - 1] + ".dat";

        ofstream gra_writer, rad_writer;

        char str_tmp[1024];
        char str_header[1024];

#ifdef WINDOWS
        strcpy_s(str_tmp, "# %s bins %03d x %03d   \n");
        sprintf_s(str_header, str_tmp, PROG_ID, bins, bins);
#else
        strcpy(str_tmp, "# %s bins %03d x %03d   \n");
        sprintf(str_header, str_tmp, PROG_ID, bins, bins);
#endif

        rad_writer.open(rad_filename.c_str(), ios::out);
        if(rad_writer.fail())
        {
            cout << ERROR_LINE << "Cannot write to:\n" << rad_filename
                    << endl;
            res = false;
            continue;
        }

        gra_writer.open(gra_filename.c_str(), ios::out);
        if(gra_writer.fail())
        {
            cout << ERROR_LINE << "Cannot write to:\n" << gra_filename
                    << endl;
            res = false;
            continue;
        }

        double sgx, sgy, sgz, tx, ty, tz;

        rad_writer << str_header;
        gra_writer << str_header;

        for(int j = -b_limit; j <= b_limit; j++)
        {
            if(j == 0)
                continue;

            for(int k = -b_limit; k <= b_limit; k++)
            {
                if(k == 0)
                    continue;

                switch(i)
                {
                    case PROJ_XY:
                        sgx = CMathFunctions::sgn(j);
                        sgy = CMathFunctions::sgn(k);
                        sgz = 0;
                        tx = double(j) * xyz_step - sgx * off_xyz;
                        ty = double(k) * xyz_step - sgy * off_xyz;
                        tz = 0.0;
                        break;

                    case PROJ_XZ:
                        sgx = CMathFunctions::sgn(j);
                        sgy = 0;
                        sgz = CMathFunctions::sgn(k);
                        tx = double(j) * xyz_step - sgx * off_xyz;
                        ty = 0.0;
                        tz = double(k) * xyz_step - sgz * off_xyz;
                        break;

                    default:
                        sgx = 0;
                        sgy = CMathFunctions::sgn(j);
                        sgz = CMathFunctions::sgn(k);
                        tx = 0.0;
                        ty = double(j) * xyz_step - sgy * off_xyz;
                        tz = double(k) * xyz_step - sgz * off_xyz;
                }

                pp->setPosition(Vector3D(tx, ty, tz));

                if(!grid->positionPhotonInGrid(pp))
                {
                    //rad_gr << Vector3D(1, 0, 0);
                    //gra_gr << Vector3D(1, 0, 0);
                    rad_writer << 1 << " " << 0 << " " << 0 << endl;
                    gra_writer << 1 << " " << 0 << " " << 0 << endl;
                }
                else
                {
                    double FgraX = pp->getPositionCell()->getData(data_off + 0);
                    double FgraY = pp->getPositionCell()->getData(data_off + 1);
                    double FgraZ = pp->getPositionCell()->getData(data_off + 2);

                    double FradX = pp->getPositionCell()->getData(data_off + 3 + 0);
                    double FradY = pp->getPositionCell()->getData(data_off + 3 + 1);
                    double FradZ = pp->getPositionCell()->getData(data_off + 3 + 2);

                    //if(k == 1 && j > 0)
                    //{
                    //v_Frad[j-1] += Vector3D(FradX, FradY, FradZ).length();
                    //v_Fgra[j-1] += Vector3D(FgraX, FgraY, FgraZ).length();
                    //}

                    //rad_gr << Vector3D(FradX, FradY, FradZ);
                    //gra_gr << Vector3D(FgraX, FgraY, FgraZ);
                    rad_writer << FradX << " " << FradY << " " << FradZ << endl;
                    gra_writer << FgraX << " " << FgraY << " " << FgraZ << endl;
                }

                per_counter++;
            }
        }

        //rad_gr.setPathOutput(path_plot + "Frad" + prev[i - 1] + ".svg");
        //gra_gr.setPathOutput(path_plot + "Fgra" + prev[i - 1] + ".svg");

        //rad_gr.createSVGFile();
        //gra_gr.createSVGFile();

        //rad_gr.rescaleVectorMap(30);
        //gra_gr.rescaleVectorMap(30);

        //rad_gr.createVectors(true, true, 1);
        //gra_gr.createVectors(true, true, 1);

        //rad_gr.writeSVGFile();
        //gra_gr.writeSVGFile();

        rad_writer.close();
        gra_writer.close();

    }

    photon_package pp = photon_package();
    uint pos_counter, avg_counter = 0;

    //0
    pp->setPosition(Vector3D(0, 0, 0));
    pp->setDirection(Vector3D(1.001, 1, 1).normalized());
    grid->positionPhotonInGrid(pp);

    pos_counter = 0;

    while(grid->next(pp))
    {
        if(pos_counter == 0)
            avg_counter++;

        double FgraX = pp->getPositionCell()->getData(data_off + 0);
        double FgraY = pp->getPositionCell()->getData(data_off + 1);
        double FgraZ = pp->getPositionCell()->getData(data_off + 2);

        double FradX = pp->getPositionCell()->getData(data_off + 3 + 0);
        double FradY = pp->getPositionCell()->getData(data_off + 3 + 1);
        double FradZ = pp->getPositionCell()->getData(data_off + 3 + 2);

        v_pos[pos_counter] = pp->getPosition().X();
        v_Frad[pos_counter] += Vector3D(FradX, FradY, FradZ).length();
        v_Fgra[pos_counter] += Vector3D(FgraX, FgraY, FgraZ).length();

        pos_counter++;
    }

    //1
    pp->setPosition(Vector3D(0, 0, 0));
    pp->setDirection(Vector3D(1.001, 1, -1).normalized());
    grid->positionPhotonInGrid(pp);

    pos_counter = 0;
    while(grid->next(pp))
    {
        if(pos_counter == 0)
            avg_counter++;

        double FgraX = pp->getPositionCell()->getData(data_off + 0);
        double FgraY = pp->getPositionCell()->getData(data_off + 1);
        double FgraZ = pp->getPositionCell()->getData(data_off + 2);

        double FradX = pp->getPositionCell()->getData(data_off + 3 + 0);
        double FradY = pp->getPositionCell()->getData(data_off + 3 + 1);
        double FradZ = pp->getPositionCell()->getData(data_off + 3 + 2);

        v_Frad[pos_counter] += Vector3D(FradX, FradY, FradZ).length();
        v_Fgra[pos_counter] += Vector3D(FgraX, FgraY, FgraZ).length();

        pos_counter++;
    }

    //2
    pp->setPosition(Vector3D(0, 0, 0));
    pp->setDirection(Vector3D(1.001, -1, 1).normalized());
    grid->positionPhotonInGrid(pp);

    pos_counter = 0;
    while(grid->next(pp))
    {
        if(pos_counter == 0)
            avg_counter++;

        double FgraX = pp->getPositionCell()->getData(data_off + 0);
        double FgraY = pp->getPositionCell()->getData(data_off + 1);
        double FgraZ = pp->getPositionCell()->getData(data_off + 2);

        double FradX = pp->getPositionCell()->getData(data_off + 3 + 0);
        double FradY = pp->getPositionCell()->getData(data_off + 3 + 1);
        double FradZ = pp->getPositionCell()->getData(data_off + 3 + 2);

        v_Frad[pos_counter] += Vector3D(FradX, FradY, FradZ).length();
        v_Fgra[pos_counter] += Vector3D(FgraX, FgraY, FgraZ).length();

        pos_counter++;
    }

    //3
    pp->setPosition(Vector3D(0, 0, 0));
    pp->setDirection(Vector3D(1.001, -1, -1).normalized());
    grid->positionPhotonInGrid(pp);

    pos_counter = 0;
    while(grid->next(pp))
    {
        if(pos_counter == 0)
            avg_counter++;

        double FgraX = pp->getPositionCell()->getData(data_off + 0);
        double FgraY = pp->getPositionCell()->getData(data_off + 1);
        double FgraZ = pp->getPositionCell()->getData(data_off + 2);

        double FradX = pp->getPositionCell()->getData(data_off + 3 + 0);
        double FradY = pp->getPositionCell()->getData(data_off + 3 + 1);
        double FradZ = pp->getPositionCell()->getData(data_off + 3 + 2);

        v_Frad[pos_counter] += Vector3D(FradX, FradY, FradZ).length();
        v_Fgra[pos_counter] += Vector3D(FgraX, FgraY, FgraZ).length();

        pos_counter++;
    }

    //4
    pp->setPosition(Vector3D(0, 0, 0));
    pp->setDirection(Vector3D(-1.001, 1, 1).normalized());
    grid->positionPhotonInGrid(pp);

    pos_counter = 0;
    while(grid->next(pp))
    {
        if(pos_counter == 0)
            avg_counter++;

        double FgraX = pp->getPositionCell()->getData(data_off + 0);
        double FgraY = pp->getPositionCell()->getData(data_off + 1);
        double FgraZ = pp->getPositionCell()->getData(data_off + 2);

        double FradX = pp->getPositionCell()->getData(data_off + 3 + 0);
        double FradY = pp->getPositionCell()->getData(data_off + 3 + 1);
        double FradZ = pp->getPositionCell()->getData(data_off + 3 + 2);

        v_Frad[pos_counter] += Vector3D(FradX, FradY, FradZ).length();
        v_Fgra[pos_counter] += Vector3D(FgraX, FgraY, FgraZ).length();

        pos_counter++;
    }

    //5
    pp->setPosition(Vector3D(0, 0, 0));
    pp->setDirection(Vector3D(-1.001, 1, -1).normalized());
    grid->positionPhotonInGrid(pp);

    pos_counter = 0;
    while(grid->next(pp))
    {
        if(pos_counter == 0)
            avg_counter++;

        double FgraX = pp->getPositionCell()->getData(data_off + 0);
        double FgraY = pp->getPositionCell()->getData(data_off + 1);
        double FgraZ = pp->getPositionCell()->getData(data_off + 2);

        double FradX = pp->getPositionCell()->getData(data_off + 3 + 0);
        double FradY = pp->getPositionCell()->getData(data_off + 3 + 1);
        double FradZ = pp->getPositionCell()->getData(data_off + 3 + 2);

        v_Frad[pos_counter] += Vector3D(FradX, FradY, FradZ).length();
        v_Fgra[pos_counter] += Vector3D(FgraX, FgraY, FgraZ).length();

        pos_counter++;
    }

    //6
    pp->setPosition(Vector3D(0, 0, 0));
    pp->setDirection(Vector3D(-1.001, -1, 1).normalized());
    grid->positionPhotonInGrid(pp);
    ;

    pos_counter = 0;
    while(grid->next(pp))
    {
        if(pos_counter == 0)
            avg_counter++;

        double FgraX = pp->getPositionCell()->getData(data_off + 0);
        double FgraY = pp->getPositionCell()->getData(data_off + 1);
        double FgraZ = pp->getPositionCell()->getData(data_off + 2);

        double FradX = pp->getPositionCell()->getData(data_off + 3 + 0);
        double FradY = pp->getPositionCell()->getData(data_off + 3 + 1);
        double FradZ = pp->getPositionCell()->getData(data_off + 3 + 2);

        v_Frad[pos_counter] += Vector3D(FradX, FradY, FradZ).length();
        v_Fgra[pos_counter] += Vector3D(FgraX, FgraY, FgraZ).length();

        pos_counter++;
    }

    //7
    pp->setPosition(Vector3D(0.1, 0.1, 0.1));
    pp->setDirection(Vector3D(-1.001, -1, -1).normalized());
    grid->positionPhotonInGrid(pp);

    pos_counter = 0;
    while(grid->next(pp))
    {
        if(pos_counter == 0)
            avg_counter++;

        double FgraX = pp->getPositionCell()->getData(data_off + 0);
        double FgraY = pp->getPositionCell()->getData(data_off + 1);
        double FgraZ = pp->getPositionCell()->getData(data_off + 2);

        double FradX = pp->getPositionCell()->getData(data_off + 3 + 0);
        double FradY = pp->getPositionCell()->getData(data_off + 3 + 1);
        double FradZ = pp->getPositionCell()->getData(data_off + 3 + 2);

        v_Frad[pos_counter] += Vector3D(FradX, FradY, FradZ).length();
        v_Fgra[pos_counter] += Vector3D(FgraX, FgraY, FgraZ).length();

        pos_counter++;
    }

    pos_counter--;

    string filename = path_plot;
    filename += "radial.txt";
    ofstream writer(filename.c_str());
    spline sFrad, sFgra;
    double frac = 1.0 / double(avg_counter);

    sFrad.resize(pos_counter);
    sFgra.resize(pos_counter);

    for(uint i = 0; i < pos_counter; i++)
    {
        v_Frad[i] *= frac;
        v_Fgra[i] *= frac;

        //cout << i << "\t" << v_pos[i] << endl;


        writer << v_pos[i] << "\t" << v_Frad[i] << "\t" << v_Fgra[i] << endl;

        sFrad.setValue(i, v_pos[i], v_Frad[i]);
        sFgra.setValue(i, v_pos[i], v_Fgra[i]);
    }

    sFrad.createSpline();
    sFgra.createSpline();

    double min_pos = v_pos[0];
    double max_pos = v_pos[pos_counter - 1];
    double step = (v_pos[pos_counter - 1] - v_pos[0]) / 101;

    //writer << "++++++" << endl;

    //for(double pos = min_pos; pos <= max_pos; pos += step)
    //    writer << pos << "\t" << abs(sFrad.getValue(pos)) << "\t" << sFgra.getValue(pos)
<< endl;

    writer.close();

    delete[] v_pos;
    delete[] v_Frad;
    delete[] v_Fgra;

    return true;
}
*/
