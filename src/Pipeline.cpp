#include "Pipeline.h"
#include "CommandParser.h"
#include "Cylindrical.h"
#include "Dust.h"
#include "GasSpecies.h"
#include "Grid.h"
#include "OcTree.h"
#include "RadiativeTransfer.h"
#include "Source.h"
#include "Spherical.h"
#include "Voronoi.h"

bool CPipeline::Init(int argc, char ** argv)
{
    end = 0, len = 0;
    begin = omp_get_wtime();
    srand(0);
    CMathFunctions mf;

    string out_string2 = "* ";
    out_string2 += PROG_ID;
    out_string2 += "                                                       *";

    uint str_len = uint(out_string2.size());
    string line_string = SEP_LINE;

    for(uint i = 0; i < str_len; i++)
    {
        cout << line_string.substr(0, i) << "| \r" << flush;
        mf.sleep(LINE_DELAY);
    }
    cout << SEP_LINE;

    for(uint i = 0; i <= str_len / 4; i++)
    {
        cout << out_string2.substr(0, i) << "| \r" << flush;
        mf.sleep(LINE_DELAY);
    }

    cout << out_string2 << endl;

    for(uint i = 0; i < str_len; i++)
    {
        cout << line_string.substr(0, i) << "| \r" << flush;
        mf.sleep(LINE_DELAY);
    }
    cout << SEP_LINE;

    if(argc != 2)
    {
        cout << "\nERROR: Wrong amount of arguments!                     \n";
        cout << "       POLARIS requires only the path of a command file!            \n";
        Error();
        return false;
    }

    CCommandParser parser(argv[1]); /**/

    // CCommandParser parser("/home/s0reissl/polaris projects/Francois/cmd_file");
    // CCommandParser parser("/home/s0reissl/polaris projects/Camilo/dustPolaris.cmd");

    if(!parser.parse())
    {
        Error();
        return false;
    }

    param_list = parser.getParameterList();

    if(param_list.size() == 0)
    {
        cout << "\nERROR: No tasks defined!" << endl;
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
    printf("  Total time of processing: %luh %02lumin. %02lusec.  \n", h, m, s);

    cout << CLR_LINE;
    cout << SEP_LINE;
    cout << "* POLARIS SUCCESSFULLY FINISHED                                             "
            "        *"
         << endl;
    cout << SEP_LINE << endl;

    string opath = path_data + "time.txt";
    ofstream t_writer(opath.c_str());
    t_writer << "Total time:\r" << endl;
    t_writer << h << "h " << m << "min " << s << "sec\r\n\r" << endl;
    t_writer << len << " ms\r" << endl;
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
                
            case CMD_DUST_TIME:
                result = calcMonteCarloTimeTransfer(param);
                break;

            case CMD_LINE_EMISSION:
                result = calcChMapsViaRayTracing(param);
                break;

            case CMD_FORCE:
                // Needs update!!!
                // result = calcRadPressure(param);
                break;

            case CMD_OPIATE:
                // TBD
                break;

            case CMD_SYNCHROTRON:
                result = calcPolarizationMapsViaSynchrotron(param);
                break;

            default:
                cout << "\nERROR: Command is unknown!" << endl;
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
    printf("  Total time of processing: %luh %02lumin. %02lusec.  \n", h, m, s);

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

    if(!createWavelengthList(param, dust))
        return false;

    if(!assignDustMixture(param, dust, grid))
        return false;

    grid->setSIConversionFactors(param);

    grid->setSpecLengthAsVector(use_energy_density);
    if(!grid->loadGridFromBinrayFile(param, use_energy_density ? 4 * WL_STEPS : WL_STEPS))
        return false;

    // Print helpfull information
    grid->createCellList();
    dust->printParameter(param, grid);
    grid->printParameters();

    if(!grid->writeMidplaneFits(path_data + "input_", param, param.getInpMidDataPoints(), true))
        return false;

    if(!grid->writeGNUPlotFiles(path_plot + "input_", param))
        return false;

    if(!grid->writeAMIRAFiles(path_plot + "input_", param, param.getInpAMIRAPoints()))
        return false;

    createSourceLists(param, dust, grid);
    if(sources_mc.size() == 0)
    {
        cout << "\nERROR: No sources for Monte-Carlo simulations defined!" << endl;
        return false;
    }

    CRadiativeTransfer rad(param);

    rad.setGrid(grid);
    rad.setDust(dust);
    rad.setSourcesLists(sources_mc, sources_ray);
    rad.initiateRadFieldMC(param);

    omp_set_num_threads(param.getNrOfThreads());

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

    if(!grid->writeGNUPlotFiles(path_plot + "output_", param))
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

bool CPipeline::calcMonteCarloTimeTransfer(parameters & param)
{   
    // Check same parameters as in MCRadiationField
    
    bool use_energy_density = false;
    if(param.getSaveRadiationField() || param.isRatSimulation() || param.getStochasticHeatingMaxSize() > 0)
        use_energy_density = true;

    CGridBasic * grid = 0;
    CDustMixture * dust = new CDustMixture();
 
    if(!createOutputPaths(param.getPathOutput()))
        return false;

    if(!assignGridType(grid, param))
        return false;

    if(!createWavelengthList(param, dust))
        return false;

    if(!assignDustMixture(param, dust, grid))
        return false; 
    
    grid->setSIConversionFactors(param);
    
    grid->setSpecLengthAsVector(use_energy_density);
    if(!grid->loadGridFromBinrayFile(param, use_energy_density ? 4 * WL_STEPS : WL_STEPS))
        return false;
    
    // Print helpfull information (Not working jet - Dust.cpp has to be modified)
    grid->createCellList();
    dust->printParameter(param, grid);
    grid->printParameters();
    
    if(!grid->writeMidplaneFits(path_data + "input_", param, param.getInpMidDataPoints(), true))
        return false;

    if(!grid->writeGNUPlotFiles(path_plot + "input_", param))
        return false;

    if(!grid->writeAMIRAFiles(path_plot + "input_", param, param.getInpAMIRAPoints()))
        return false;
    
    createSourceLists(param, dust, grid);
    if(sources_mc.size() == 0)
    {
        cout << "\nERROR: No sources for Monte-Carlo simulations defined!" << endl;
        return false;
    }
    
    // Detector similiar to calcPolarizationMapsViaMC 
    CDetector * detector = createDetectorList(param, dust, grid);
    
    if(detector == 0)
        return false;
    
    CRadiativeTransfer rad(param);
    
    rad.setGrid(grid);
    rad.setDust(dust);
    rad.setSourcesLists(sources_mc, sources_ray);
    rad.setDetectors(detector);
    
    // Both initiates should be included in initiateTimeMC
    rad.initiateDustMC(param);
    rad.initiateRadFieldMC(param);
    
    omp_set_num_threads(param.getNrOfThreads());
    
    // Here isTemperatureSimulation has to be adapated to TdMCRT
    if(param.isTemperatureSimulation())
    {
        if(param.getDustOffset())
            rad.convertTempInQB(param.getOffsetMinGasDensity(), false);
        else if(param.getDustGasCoupling())
            rad.convertTempInQB(param.getOffsetMinGasDensity(), true);
    }
    
    
    rad.calcMonteCarloTimeTransfer(param.getCommand(),
                                   use_energy_density,
                                   false);
    
    //Following is probably not needed due to no final temperature
    //if(param.isTemperatureSimulation())
    //    rad.calcFinalTemperature(use_energy_density);
    
    //Alignement theorys not jet included
    //if(param.isRatSimulation())
    //    rad.calcAlignedRadii();
    
    cout << SEP_LINE;
    
    
    //Here only final output, intermediate results will be saved somewhere else
    
    if(!grid->writeMidplaneFits(path_data + "output_", param, param.getOutMidDataPoints()))
        return false;

    if(!grid->writeGNUPlotFiles(path_plot + "output_", param))
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

    //Clean up
    delete grid;
    delete dust;
    delete[] detector;
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

    if(!createWavelengthList(param, dust))
        return false;

    if(!assignDustMixture(param, dust, grid))
        return false;

    grid->setSIConversionFactors(param);

    if(!grid->loadGridFromBinrayFile(param, 0))
        return false;

    // Print helpfull information
    grid->createCellList();
    dust->printParameter(param, grid);
    grid->printParameters();

    if(!grid->writeMidplaneFits(path_data + "input_", param, param.getInpMidDataPoints(), true))
        return false;

    if(!grid->writeGNUPlotFiles(path_plot + "input_", param))
        return false;

    if(!grid->writeAMIRAFiles(path_plot + "input_", param, param.getInpAMIRAPoints()))
        return false;

    createSourceLists(param, dust, grid);
    if(sources_mc.size() == 0)
    {
        cout << "\nERROR: No sources for Monte-Carlo simulations defined!" << endl;
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

    omp_set_num_threads(param.getNrOfThreads());

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

    if(!createWavelengthList(param, dust))
        return false;

    if(!assignDustMixture(param, dust, grid))
        return false;

    grid->setSIConversionFactors(param);

    if(!grid->loadGridFromBinrayFile(param, getNrOffsetEntriesRay(param, dust, grid)))
        return false;

    // Print helpfull information
    grid->createCellList();
    dust->printParameter(param, grid);
    grid->printParameters();

    if(!grid->writeMidplaneFits(path_data + "input_", param, param.getInpMidDataPoints(), true))
        return false;

    if(!grid->writeGNUPlotFiles(path_plot + "input_", param))
        return false;

    if(!grid->writeAMIRAFiles(path_plot + "input_", param, param.getInpAMIRAPoints()))
        return false;

    createSourceLists(param, dust, grid);
    if(sources_ray.size() == 0)
    {
        cout << "\nERROR: No sources for raytracing simulations defined!" << endl;
        return false;
    }

    CRadiativeTransfer rad(param);

    rad.setGrid(grid);
    rad.setDust(dust);
    rad.setSourcesLists(sources_mc, sources_ray);

    if(!rad.initiateDustRaytrace(param))
        return false;

    omp_set_num_threads(param.getNrOfThreads());

    if(param.getStochasticHeatingMaxSize() > 0)
        rad.calcStochasticHeating();

    // Calculate radiation field before raytracing (if sources defined and no radiation
    // field in grid)
    if(!grid->getRadiationFieldAvailable() && dust->getScatteringToRay() && !sources_mc.empty())
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

    if(!createWavelengthList(param, dust, gas))
        return false;

    if(!assignDustMixture(param, dust, grid))
        return false;

    grid->setSIConversionFactors(param);

    if(!grid->loadGridFromBinrayFile(param, 6 * param.getMaxNrOfLineTransitions()))
        return false;

    // Print helpfull information
    grid->createCellList();
    dust->printParameter(param, grid);
    gas->printParameter(param, grid);
    grid->printParameters();

    if(!grid->writeMidplaneFits(path_data + "input_", param, param.getInpMidDataPoints(), true))
        return false;

    if(!grid->writeGNUPlotFiles(path_plot + "input_", param))
        return false;

    if(!grid->writeAMIRAFiles(path_plot + "input_", param, param.getInpAMIRAPoints()))
        return false;

    createSourceLists(param, dust, grid);
    if(sources_ray.size() == 0)
    {
        cout << "\nERROR: No sources for raytracing simulations defined!" << endl;
        return false;
    }

    CRadiativeTransfer rad(param);

    rad.setGrid(grid);
    rad.setDust(dust);
    rad.setGas(gas);
    rad.setSourcesLists(sources_mc, sources_ray);

    if(!rad.initiateLineRaytrace(param))
        return false;

    omp_set_num_threads(param.getNrOfThreads());

    if(!rad.calcChMapsViaRaytracing(param))
        return false;

    delete grid;
    delete dust;
    delete gas;
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

    if(!createWavelengthList(param, dust))
        return false;

    if(!assignDustMixture(param, dust, grid))
        return false;

    grid->setSIConversionFactors(param);

    if(!grid->loadGridFromBinrayFile(param, 0))
        return false;

    // Print helpfull information
    grid->createCellList();
    dust->printParameter(param, grid);
    grid->printParameters();

    if(!grid->writeMidplaneFits(path_data + "input_", param, param.getInpMidDataPoints(), true))
        return false;

    if(!grid->writeGNUPlotFiles(path_plot + "input_", param))
        return false;

    if(!grid->writeAMIRAFiles(path_plot + "input_", param, param.getInpAMIRAPoints()))
        return false;

    createSourceLists(param, dust, grid);
    if(sources_ray.size() == 0)
    {
        cout << "\nERROR: No sources for raytracing simulations defined!" << endl;
        return false;
    }

    CRadiativeTransfer rad(param);

    rad.setGrid(grid);
    rad.setDust(dust);
    rad.setSourcesLists(sources_mc, sources_ray);

    if(!rad.initiateSyncRaytrace(param))
        return false;

    omp_set_num_threads(param.getNrOfThreads());
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
        cout << "\nERROR: Cannot open binary grid file:" << endl;
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
            cout << "\nERROR: Grid type unknown!" << endl;
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
        cout << "\nERROR: No Monte-Carlo detector defined!" << endl;
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
                           distance,
                           lam_min,
                           lam_max,
                           nr_spectral_bins);
        detector[pos].setOrientation(axis1, axis2, rot_angle_1, rot_angle_2);
        /*if(param.getPeelOff() && param.getNrOfDustPhotons() != 0)
        {
            cout << "\nHINT: Peel-off technique disabled for self-scattering of dust grain
        emission!" << endl; param.setPeelOff(false);
        }*/
        if(param.getPeelOff() && param.getAcceptanceAngle() > 1.0)
            cout << "\nHINT: Peel-off technique needs no acceptance angle!" << endl;
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

    if(param.isRaytracing())
    {
        // Raytracing simulations are only using background sources!
        if(param.getNrOfDiffuseSources() > 0)
        {
            cout << "\nWARNING: Diffuse sources cannot be considered in "
                 << "dust, line, or synchrotron emission!" << endl;
            nr_ofSources--;
        }

        // if(param.getISRFSource())
        // {
        //     cout << "\nWARNING: ISRF as radiation source cannot be considered in "
        //          << "dust, line, or synchrotron emission!" << endl;
        //     nr_ofSources--;
        // }

        if(param.getNrOfBackgroundSources() == 0)
        {
            cout << "\nHINT: No background source was defined!" << endl;
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
                        cout << "\nERROR: Background source nr. " << s / NR_OF_BG_SOURCES + 1 << " undefined!"
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
                nr_ofSources--;
            else if(param.getCommand() != CMD_DUST_EMISSION)
                cout << "\nWARNING: ISRF source cannot be considered in line or synchrotron emission!"
                     << endl;
            else
            {
                CSourceBasic * tmp_source = new CSourceISRF();

                cout << "-> Creating ISRF source list             \r" << flush;

                if(param.getISRFPath() != "")
                {
                    if(!tmp_source->setParameterFromFile(param, grid, dust, 0))
                    {
                        cout << "\nERROR: Interstellar radiation field undefined! \n" << flush;
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
            if(!dust->getScatteringToRay())
                nr_ofSources -= param.getNrOfPointSources();
            else if(param.getCommand() != CMD_DUST_EMISSION)
                cout << "\nWARNING: Point sources cannot be considered in line or "
                        "synchrotron emission!"
                     << endl;
            else
            {
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
                            cout << "\nERROR: Star source nr. " << s / NR_OF_POINT_SOURCES + 1
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
            if(!dust->getScatteringToRay())
                nr_ofSources -= param.getNrOfLaserSources();
            else if(param.getCommand() != CMD_DUST_EMISSION)
                cout << "\nWARNING: Laser sources cannot be considered in line or "
                        "synchrotron emission!"
                     << endl;
            else
            {
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
            if(!dust->getScatteringToRay() && grid->getRadiationFieldAvailable())
                nr_ofSources--;
            else if(param.getCommand() != CMD_DUST_EMISSION)
                cout << "\nWARNING: Dust source cannot be considered in line or "
                        "synchrotron emission!"
                     << endl;
            else
            {
                CSourceBasic * tmp_source = new CSourceDust();
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
                    cout << "\nERROR: Star source nr. " << s / NR_OF_POINT_SOURCES + 1 << " undefined!"
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
                    cout << "\nERROR: Sorce Starfield nr. " << s / NR_OF_DIFF_SOURCES + 1 << " undefined! \n"
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
            cout << "\nERROR: Background sources can only be used for raytracing "
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
                    cout << "\nERROR: Interstellar radiation field undefined! \n" << flush;
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
                cout << "\nERROR: Dust as radiation source cannot be considered in "
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
        cout << "\nERROR: Not enough source(s) were initiated! \n" << flush;

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
        cout << "\nERROR: Can plot sources to:" << endl;
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

    // Write Gnuplot files to show dust properties
    if(!dust->writeComponent(path_data, path_plot))
        return false;

    dust->setGridRequirements(grid, param);

    // Check if either the radiation field is present or the radiation field can be
    // calculated otherwise, disable scattering added to the raytracing then
    if(param.getScatteringToRay())
        !grid->getRadiationFieldAvailable() && param.getNrOfPointSources() == 0 &&
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
            printPlotParameters(param);
            break;

        case CMD_FORCE:
            cout << "- Command          : RADIATION FORCE" << endl;
            printPathParameters(param);
            printSourceParameters(param);
            printConversionParameters(param);
            printPlotParameters(param);
            break;

        case CMD_TEMP_RAT:
            cout << "- Command          : TEMPERATURE DISTRIBUTION and RAT ALIGNMENT" << endl;
            printPathParameters(param);
            printSourceParameters(param);
            printConversionParameters(param);
            printPlotParameters(param);
            break;

        case CMD_RAT:
            cout << "- Command          : RAT ALIGNMENT" << endl;
            printPathParameters(param);
            printSourceParameters(param, true);
            printConversionParameters(param);
            printPlotParameters(param);
            break;

        case CMD_DUST_EMISSION:
            cout << "- Command          : DUST EMISSION" << endl;
            printPathParameters(param);
            printSourceParameters(param, true);
            printConversionParameters(param);
            printAlignmentParameters(param);
            printDetectorParameters(param);
            printPlotParameters(param);
            break;

        case CMD_SYNCHROTRON:
            cout << "- Command          : SYNCHROTRON EMISSION" << endl;
            printPathParameters(param);
            printSourceParameters(param);
            printConversionParameters(param);
            printDetectorParameters(param);
            printPlotParameters(param);
            printSynchrotronParameters(param);
            break;

        case CMD_DUST_SCATTERING:
            cout << "- Command          : DUST SCATTERING (Monte-Carlo)" << endl;
            printPathParameters(param);
            printSourceParameters(param, true);
            printConversionParameters(param);
            printAlignmentParameters(param);
            printDetectorParameters(param, true);
            printPlotParameters(param);
            break;
            
        case CMD_DUST_TIME:
            cout << "- Command          : DUST TIME TRANSFER (Monte-Carlo)" << endl;
            // Parameter output tbd
            break;
            
        case CMD_OPIATE:
            cout << "todo: OPIATE parameter" << endl;
            break;
            
        case CMD_LINE_EMISSION:
            cout << "- Command          : SPECTRAL LINE EMISSION" << endl;
            printPathParameters(param);
            printSourceParameters(param);
            printConversionParameters(param);
            printDetectorParameters(param);
            printPlotParameters(param);
            break;

        default:
            cout << "\nERROR: Command is unknown!" << endl;
            cout << "No parameters available! " << endl;
    };
    cout << SEP_LINE;
}

bool CPipeline::createWavelengthList(parameters & param, CDustMixture * dust, CGasMixture * gas)
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
                cout << "\nERROR: No dust raytracing detector defined (see "
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
                cout << "\nERROR: No dust Monte-Carlo detector defined (see "
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

            // Get detector parameters list
            maplist line_ray_detector_list = param.getLineRayDetectors();
            maplist::iterator it;

            // Check if a detector is defined
            if(line_ray_detector_list.empty())
            {
                cout << "\nERROR: No spectral line detector of gas species defined!" << endl;
                return false;
            }

            // Perform radiative transfer for each chosen gas species
            for(it = line_ray_detector_list.begin(); it != line_ray_detector_list.end(); ++it)
            {
                // Get ID of the current gas species
                uint i_species = it->first;

                // Get number of spectral line transitions that have to be simulated
                uint nr_of_transitions = param.getNrOfGasSpeciesTransitions(i_species);

                // Perform radiative transfer for each chosen spectral line transition
                for(uint i_line = 0; i_line < nr_of_transitions; i_line++)
                {
                    // Calculate from frequency
                    double wavelength = con_c / gas->getTransitionFrequencyFromIndex(i_species, i_line);

                    // Add wavelength to global list of wavelength
                    dust->addToWavelengthGrid(wavelength);
                }
            }
            break;
        }

        case CMD_SYNCHROTRON:
            // Get detector parameters list
            values = param.getSyncRayDetectors();

            // Check if a detector is defined
            if(values.empty())
            {
                cout << "\nERROR: No synchrotron detector defined (see <detector_sync>)!" << endl;
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

    if(!grid->loadGridFromBinrayFile(param, 3 + 3 * WL_STEPS))
        return false;

    // Print helpfull information
    grid->createCellList();
    dust->printParameter(param, grid);
    grid->printParameter();

    if(!grid->writeMidplaneFits(path_data + "input_", param, param.getInpMidDataPoints(),
true)) return false;

    if(!grid->writeGNUPlotFiles(path_plot + "input_", param))
        return false;

    if(!grid->writeAMIRAFiles(path_plot + "input_", param, param.getInpAMIRAPoints()))
        return false;

    createSourceLists(param, dust, grid);
    if(sources_mc.size() == 0)
    {
        cout << "\nERROR: No sources for Monte-Carlo simulations defined!" << endl;
        return false;
    }

    CRadiativeTransfer rad(param);

    rad.setGrid(grid);
    rad.setDust(dust);
    rad.setSourcesLists(sources_mc, sources_ray);
    rad.initiateRadFieldMC(param);

    omp_set_num_threads(param.getNrOfThreads());

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

        /*for(long c_2 = c_1 + 1; c_2 < long(max_cells - 1); c_2++)
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

#pragma omp critical
        {
            per_counter++;
            if(per_counter % 100 == 0)
                cout << "-> Calculation of final properties: [ "
                    << 100.0 * float(per_counter) / float(max_cells)
                << " %]          \r";
        }
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
        photon_package * pp = new photon_package;

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
            cout << "\nERROR: Cannot write to:\n" << rad_filename
                    << endl;
            res = false;
            continue;
        }

        gra_writer.open(gra_filename.c_str(), ios::out);
        if(gra_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n" << gra_filename
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

#pragma omp critical
                {
                    if(per_counter % 220 == 0)
                        cout << " -> Writing midplane files: "
                            << 100.0 * float(per_counter) / float(per_max)
                        << " [%]             \r";
                }
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

        delete pp;
    }

    photon_package * pp = new photon_package;
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

    delete pp;

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
