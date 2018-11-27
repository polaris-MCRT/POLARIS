#include "GasSpecies.h"
#include "CommandParser.h"
#include "MathFunctions.h"
#include "typedefs.h"
#include "Grid.h"

#define TRANS_SIGMA_P  +1
#define TRANS_PI       0
#define TRANS_SIGMA_M  -1

// This function is based on
// Mol3d: 3D line and dust continuum radiative transfer code
// by Florian Ober 2015, Email: fober@astrophysik.uni-kiel.de

bool CGasSpecies::calcLTE(CGridBasic * grid)
{
    uint const nr_of_energy_levels = getNrEnergyLevels();
    uint const nr_of_transitions = getNrOfTransitions();
    float last_percentage = 0;
    long cell_count = 0;
    long max_cells = grid->getMaxDataCells();

    cout << CLR_LINE;

#pragma omp parallel for schedule(dynamic)
    for(long i_cell = 0; i_cell < long(max_cells); i_cell++)
    {
        cell_basic * cell = grid->getCellFromIndex(i_cell);
        // Calculate percentage of total progress per source
        float percentage = 100 * float(cell_count) / float(max_cells);

        // Show only new percentage number if it changed
        if((percentage - last_percentage) > PERCENTAGE_STEP)
        {
#pragma omp critical
            {
                cout << "-> Calculating LTE level population  : " << percentage << " [%]               \r";
                last_percentage = percentage;
            }
        }

        cell_count++;
        double * tmp_lvl_pop; // , tmp_lvl_pop_u, tmp_lvl_pop_l;
        tmp_lvl_pop = new double[nr_of_energy_levels];

        double temp_gas = grid->getGasTemperature(cell);
        if(temp_gas == 0)
        {
            tmp_lvl_pop[0] = 1;
            for(uint i = 1; i < nr_of_energy_levels; i++)
                tmp_lvl_pop[i] = 0;
        }
        else
        {
            double sum = 0;
            for(uint i = 0; i < nr_of_energy_levels; i++)
            {
                tmp_lvl_pop[i] = getGLevel(i) * exp(-con_h * getEnergylevel(i) * con_c * 100.0 / (con_kB * temp_gas));
                sum += tmp_lvl_pop[i];
            }

            for(uint i = 0; i < nr_of_energy_levels; i++)
                tmp_lvl_pop[i] /= sum;
        }

        for(uint i_line = 0; i_line < nr_of_transitions; i_line++)
        {
            grid->setLvlPopLower(cell, i_line, tmp_lvl_pop[getLowerTransition(getTransition(i_line))]);
            grid->setLvlPopUpper(cell, i_line, tmp_lvl_pop[getUpperTransition(getTransition(i_line))]);
        }

        delete[] tmp_lvl_pop;
    }
    return true;
}

// This function is based on
// Mol3d: 3D line and dust continuum radiative transfer code
// by Florian Ober 2015, Email: fober@astrophysik.uni-kiel.de

bool CGasSpecies::calcFEP(CGridBasic * grid)
{
    uint nr_of_energy_levels = getNrEnergyLevels();
    uint nr_of_total_transitions = getNrOfTotalTransitions();
    uint nr_of_transitions = getNrOfTransitions();
    float last_percentage = 0;
    long cell_count = 0;
    long max_cells = grid->getMaxDataCells();
    double * J_mid = new double[nr_of_total_transitions];
    bool no_error = true;

    cout << CLR_LINE;

    for(uint i_trans = 0; i_trans < nr_of_total_transitions; i_trans++)
        J_mid[i_trans] = CMathFunctions::planck_hz(getTransitionFrequency(i_trans), 2.75);

#pragma omp parallel for schedule(dynamic)
    for(long i_cell = 0; i_cell < long(max_cells); i_cell++)
    {
        cell_basic * cell = grid->getCellFromIndex(i_cell);
        double * tmp_lvl_pop = new double[nr_of_energy_levels];

        double temp_gas = grid->getGasTemperature(cell);
        double gas_number_density = grid->getGasNumberDensity(cell);

        // Calculate percentage of total progress per source
        float percentage = 100 * float(cell_count) / float(max_cells);

        // Show only new percentage number if it changed
        if((percentage - last_percentage) > PERCENTAGE_STEP)
        {
#pragma omp critical
            {
                cout << "-> Calculating FEP level population  : "
                        << last_percentage << " [%]                       \r";
                last_percentage = percentage;
            }
        }

        if(gas_number_density > 1e-200)
        {
            Matrix2D final_col_para, A;
            double * b = new double[nr_of_energy_levels];

            A.resize(nr_of_energy_levels, nr_of_energy_levels);

            cell_count++;
            final_col_para = calc_collision_parameter(grid, temp_gas, gas_number_density);
            createMatrix(J_mid, A, b, final_col_para);
            CMathFunctions::gauss(A, b, tmp_lvl_pop, nr_of_energy_levels);

            delete[] b;
            double sum_p = 0;
            for(uint i = 0; i < nr_of_energy_levels; i++)
            {
                if(tmp_lvl_pop[i] >= 0)
                    sum_p += tmp_lvl_pop[i];
                else
                {
                    cout << "WARNING: Level population element not greater than zero!" << endl;
                    no_error = false;
                }
            }
            for(uint i = 0; i < nr_of_energy_levels; i++)
                tmp_lvl_pop[i] /= sum_p;
        }
        else
        {
            for(uint i = 0; i < nr_of_energy_levels; i++)
                tmp_lvl_pop[i] = 0;
            tmp_lvl_pop[0] = 1;
        }

        for(uint i_line = 0; i_line < nr_of_transitions; i_line++)
        {
            grid->setLvlPopLower(cell, i_line, tmp_lvl_pop[getLowerTransition(getTransition(i_line))]);
            grid->setLvlPopUpper(cell, i_line, tmp_lvl_pop[getUpperTransition(getTransition(i_line))]);
        }

        delete[] tmp_lvl_pop;
    }

    delete[] J_mid;
    return no_error;
}

// This function is based on
// Mol3d: 3D line and dust continuum radiative transfer code
// by Florian Ober 2015, Email: fober@astrophysik.uni-kiel.de
// and based on Pavlyuchenkov et. al (2007)

bool CGasSpecies::calcLVG(CGridBasic * grid, double kepler_star_mass)
{
    CMathFunctions mf;
    uint nr_of_energy_levels = getNrEnergyLevels();
    uint nr_of_transitions = getNrOfTransitions();
    uint nr_of_total_transitions = getNrOfTotalTransitions();
    float last_percentage = 0;
    long cell_count = 0;
    long max_cells = grid->getMaxDataCells();
    bool no_error = true;

    cout << CLR_LINE;

    double * J_ext = new double[nr_of_total_transitions];
    for(uint i_trans = 0; i_trans < nr_of_total_transitions; i_trans++)
    {
        J_ext[i_trans] = mf.planck_hz(getTransitionFrequency(i_trans), 2.75);
    }

#pragma omp parallel for schedule(dynamic)
    for(long i_cell = 0; i_cell < long(max_cells); i_cell++)
    {
        cell_basic * cell = grid->getCellFromIndex(i_cell);
        Matrix2D final_col_para, A;
        double * J_mid = new double[nr_of_total_transitions];

        A.resize(nr_of_energy_levels, nr_of_energy_levels);
        double * b = new double[nr_of_energy_levels];
        double * new_pop = new double[nr_of_energy_levels];
        double * old_pop = new double[nr_of_energy_levels];

        double turbulent_velocity = grid->getTurbulentVelocity(cell);

        // Calculate percentage of total progress per source
        float percentage = 100 * float(cell_count) / float(max_cells);

        // Show only new percentage number if it changed
        if((percentage - last_percentage) > PERCENTAGE_STEP)
        {
#pragma omp critical
            {
                cout << "-> Calculating LVG level population : "
                        << last_percentage << " [%]                       \r";
                last_percentage = percentage;
            }
        }

        cell_count++;

        for(uint i = 1; i < nr_of_energy_levels; i++)
        {
            new_pop[i] = 0.0;
            old_pop[i] = 0.0;
        }
        old_pop[0] = 1.0;
        new_pop[0] = 1.0;

        double temp_gas = grid->getGasTemperature(cell);
        double gas_number_density = grid->getGasNumberDensity(cell);
        double dens_species = getNumberDensity(grid, cell);
        if(temp_gas == 0.0 || dens_species < 1.0e-200)
            continue;
        Vector3D pos_xyz_cell = grid->getCenter(cell);
        double abs_vel;
        if(kepler_star_mass > 0)
        {
            Vector3D velo = mf.calcKeplerianVelocity(pos_xyz_cell, kepler_star_mass);
            abs_vel = sqrt(pow(velo.X(), 2) + pow(velo.Y(), 2) + pow(velo.Z(), 2));
        }
        else if(grid->getDataID() >= 6)
        {
            Vector3D velo = grid->getVelocityField(cell);
            abs_vel = sqrt(pow(velo.X(), 2) + pow(velo.Y(), 2) + pow(velo.Z(), 2));
        }
        else
            abs_vel = 0.0;

        if(getGaussA(temp_gas, turbulent_velocity) * abs_vel < 1.0e-16)
            continue;

        double R_mid = sqrt(pow(pos_xyz_cell.X(), 2) + pow(pos_xyz_cell.Y(), 2));
        double L = R_mid * sqrt(2.0 / 3.0 / (getGaussA(temp_gas, turbulent_velocity) * abs_vel));
        uint i_iter = 0;

        final_col_para = calc_collision_parameter(grid, temp_gas, gas_number_density);

        for(i_iter = 0; i_iter < MAX_LVG_ITERATIONS; i_iter++)
        {
            for(uint i = 0; i < nr_of_energy_levels; i++)
                old_pop[i] = new_pop[i];

            for(uint i_trans = 0; i_trans < nr_of_total_transitions; i_trans++)
                J_mid[i_trans] = elem_LVG(grid, dens_species, new_pop,
                    getGaussA(temp_gas, turbulent_velocity), L, J_ext[i_trans], i_trans);

            createMatrix(J_mid, A, b, final_col_para);

            mf.gauss(A, b, new_pop, nr_of_energy_levels);

            double sum_p = 0;
            for(uint j = 0; j < nr_of_energy_levels; j++)
            {
                if(new_pop[j] >= 0)
                    sum_p += new_pop[j];
                else
                {
                    cout << "WARNING: Level population element not greater than zero!" << endl;
                    cout << new_pop[j] << endl;
                    no_error = false;
                }
            }
            for(uint j = 0; j < nr_of_energy_levels; j++)
                new_pop[j] /= sum_p;

            uint j = 0;
            for(uint i = 1; i < nr_of_energy_levels; i++)
                if(new_pop[i] > new_pop[j])
                    j = i;

            if(i_iter > 1)
            {
                if(abs(new_pop[j] - old_pop[j]) /
                        (old_pop[j] + numeric_limits<double>::epsilon()) < 1.0e-2)
                {
                    break;
                }
            }
        }
        if(i_iter == MAX_LVG_ITERATIONS)
            cout << "WARNING: Maximum iteration needed in cell: " << i_cell << endl;

        for(uint i_line = 0; i_line < nr_of_transitions; i_line++)
        {
            grid->setLvlPopUpper(cell, i_line, new_pop[getUpperTransition(getTransition(i_line))]);
            grid->setLvlPopLower(cell, i_line, new_pop[getLowerTransition(getTransition(i_line))]);
        }
        delete[] J_mid;
    }
    delete[] J_ext;
    return no_error;
}


// This function is based on
// Mol3d: 3D line and dust continuum radiative transfer code
// by Florian Ober 2015, Email: fober@astrophysik.uni-kiel.de

double CGasSpecies::elem_LVG(CGridBasic * grid, double dens_species, double * new_pop,
        double gauss_a, double L, double J_ext, uint i_trans)
{
    double lvl_pop_u = new_pop[getUpperTransition(i_trans)];
    double lvl_pop_l = new_pop[getLowerTransition(i_trans)];
    double Einst_A = getEinsteinA(i_trans);
    double Einst_B_u = getEinsteinBu(i_trans);
    double Einst_B_l = getEinsteinBl(i_trans);
    double j, alpha, tau, beta, J_mid, S;

    j = dens_species * lvl_pop_u * Einst_A * gauss_a * con_eps / sqrt(PI);
    alpha = dens_species * (lvl_pop_l * Einst_B_l - lvl_pop_u * Einst_B_u) *
            gauss_a * con_eps / sqrt(PI);

    if(alpha < 1e-20)
    {
        S = 0.0;
        alpha = 0.0;
    }
    else
        S = j / alpha;

    tau = alpha * L;

    if(tau < 1.0e-6)
        beta = 1.0 - 0.5 * tau;
    else
        beta = (1.0 - exp(-tau)) / tau;

    J_mid = (1.0 - beta) * S + beta * J_ext;
    return J_mid;
}

// This function is based on
// Mol3d: 3D line and dust continuum radiative transfer code
// by Florian Ober 2015, Email: fober@astrophysik.uni-kiel.de

void CGasSpecies::createMatrix(double * J_mid, Matrix2D & A, double * b, Matrix2D final_col_para)
{
    uint nr_of_energy_levels = getNrEnergyLevels();
    uint nr_of_collision_transition = getNrCollisionTransitions(0);
    uint k = 0;

    for(uint i = 0; i < nr_of_energy_levels; i++)
    {
        for(uint j = 0; j < nr_of_total_transitions; j++)
        {
            if(getUpperTransition(j) == i)
            {
                k = getLowerTransition(j);
                A(i, k) += getEinsteinBl(j) * J_mid[j];
                A(i, i) -= getEinsteinA(j) + getEinsteinBu(j) * J_mid[j];
            }
            else if(getLowerTransition(j) == i)
            {
                k = getUpperTransition(j);
                A(i, k) += getEinsteinA(j) + getEinsteinBu(j) * J_mid[j];
                A(i, i) -= getEinsteinBl(j) * J_mid[j];
            }
        }
        for(uint i_col_transition = 0; i_col_transition < nr_of_collision_transition; i_col_transition++)
        {
            if(getLowerCollision(0, i_col_transition) == i)
            {
                k = getUpperCollision(0, i_col_transition);
                A(i, k) += final_col_para(i_col_transition, 0);
                A(i, i) -= final_col_para(i_col_transition, 1);
            }
            else if(getUpperCollision(0, i_col_transition) == i)
            {
                k = getLowerCollision(0, i_col_transition);
                A(i, k) += final_col_para(i_col_transition, 1);
                A(i, i) -= final_col_para(i_col_transition, 0);
            }
        }
    }

    for(uint i = 0; i < nr_of_energy_levels; i++)
    {
        A(0, i) = 1.0;
        b[i] = 0;
    }
    b[0] = 1.0;
}

// This function is based on
// Mol3d: 3D line and dust continuum radiative transfer code
// by Florian Ober 2015, Email: fober@astrophysik.uni-kiel.de

Matrix2D CGasSpecies::calc_collision_parameter(CGridBasic * grid, double temp_gas, double gas_number_density)
{
    uint nr_col_trans = getNrCollisionTransitions(0);
    uint hi_i = 0;

    Matrix2D final_col_para;
    final_col_para.resize(nr_col_trans, 2);
    for(uint i_col_partner = 0; i_col_partner < getNrCollisionPartner(); i_col_partner++)
    {
        double dens = 0;
        dlist col_mtr_tmp_ul, col_mtr_tmp_lu;
        col_mtr_tmp_ul.resize(nr_col_trans);
        col_mtr_tmp_lu.resize(nr_col_trans);

        switch(getOrientation_H2(i_col_partner))
        {
            case (FULL):
                dens = gas_number_density;
                break;
            case (PARA):
                dens = gas_number_density * 0.25;
                break;
            case (ORTH):
                dens = gas_number_density * 0.75;
                break;
        }

        bool br = false;
        for(uint i_col_transition = 0; i_col_transition < nr_col_trans; i_col_transition++)
            if(getUpperCollision(i_col_partner, i_col_transition) == 0)
                br = true;

        if(br == false)
        {
            if(temp_gas < getCollisionTemp(i_col_partner, 0))
                hi_i = 0;
            else if(temp_gas > getCollisionTemp(i_col_partner, getNrCollisionTemps(i_col_partner) - 1))
                hi_i = getNrCollisionTemps(i_col_partner) - 2;
            else
            {
                for(uint i_col_temp = 1; i_col_temp < getNrCollisionTemps(i_col_partner); i_col_temp++)
                {
                    if(temp_gas < getCollisionTemp(i_col_partner, i_col_temp)
                            && temp_gas >= getCollisionTemp(i_col_partner, i_col_temp - 1))
                    {
                        hi_i = i_col_temp - 1;
                        break;
                    }
                }
            }

            for(uint i_col_transition = 0; i_col_transition < nr_col_trans; i_col_transition++)
            {
                double interp = CMathFunctions::interpolate(
                        getCollisionTemp(i_col_partner, hi_i),
                        getCollisionTemp(i_col_partner, hi_i + 1),
                        getCollisionMatrix(i_col_partner, i_col_transition, hi_i),
                        getCollisionMatrix(i_col_partner, i_col_transition, hi_i + 1),
                        temp_gas);

                if(interp > 0)
                {
                    col_mtr_tmp_ul[i_col_transition] = interp * dens;

                    col_mtr_tmp_lu[i_col_transition] = col_mtr_tmp_ul[i_col_transition] *
                            getGLevel(getUpperCollision(i_col_partner, i_col_transition)) /
                            getGLevel(getLowerCollision(i_col_partner, i_col_transition)) *
                            exp(-con_h * (getEnergylevel(getUpperCollision(i_col_partner, i_col_transition)) *
                            con_c * 100.0 - getEnergylevel(getLowerCollision(i_col_partner, i_col_transition)) *
                            con_c * 100.0) / (con_kB * temp_gas));

                    final_col_para(i_col_transition, 0) += col_mtr_tmp_ul[i_col_transition];
                    final_col_para(i_col_transition, 1) += col_mtr_tmp_lu[i_col_transition];
                }
            }
        }
    }
    return final_col_para;
}

StokesVector CGasSpecies::calcEmissivities(CGridBasic *grid, photon_package * pp, uint i_line)
{
    double gauss_a = grid->getGaussA(pp, i_line);
    uint i_trans = getTransition(i_line);
    // Calculate the optical depth of the gas particles and the dust grains
    double alpha_gas = (grid->getLvlPopLower(pp, i_line) * getEinsteinBl(i_trans)
        - grid->getLvlPopUpper(pp, i_line) * getEinsteinBu(i_trans)) * con_eps * gauss_a;

    return StokesVector(grid->getLvlPopUpper(pp, i_line) * getEinsteinA(i_trans) * con_eps * gauss_a,
        0, 0, 0, alpha_gas);
}

Matrix2D CGasSpecies::getGaussLineMatrix(CGridBasic * grid, photon_package * pp, uint i_line, double velocity)
{
    Matrix2D line_matrix(4, 4);
    double gauss_a = grid->getGaussA(pp, i_line);
    double line_amplitude = exp(-(pow(velocity, 2) * pow(gauss_a, 2))) / PIsq;
    for(uint i = 0; i < 4; i++)
        line_matrix(i, i) = line_amplitude;
    return line_matrix;
}

Matrix2D CGasSpecies::getZeemanSplittingMatrix(CGridBasic *grid, photon_package * pp, uint i_line, double velocity,
        Vector3D mag_field, double cos_theta, double sin_theta, double cos_2_phi, double sin_2_phi)
{
    // Init variables
    Matrix2D line_matrix(4, 4);
    double line_strength, freq_shift;
    double f_doppler, mult_A, mult_B;
    double frequency = getTransitionFrequencyFromIndex(i_line);
    uint i_pi = 0, i_sigma_p = 0, i_sigma_m = 0;

    double Gamma = grid->getGamma(pp, i_line);
    double doppler_width = grid->getDopplerWidth(pp, i_line);
    double voigt_a = grid->getVoigtA(pp, i_line);

    // Init the current value of the line function as a complex value
    complex<double> line_function;

    // Calculate the contribution of each allowed transition between Zeeman sublevels
    for(float i_sublvl_u = -getMaxMUpper(i_line); i_sublvl_u <= getMaxMUpper(i_line); i_sublvl_u++)
    {
        for(float i_sublvl_l = max(i_sublvl_u - 1, -getMaxMLower(i_line));
                i_sublvl_l <= min(i_sublvl_u + 1, getMaxMLower(i_line)); i_sublvl_l++)
        {
            // Calculate the frequency shift in relation to the unshifted line peak
            // Delta nu = (B * mu_Bohr) / h * (m' * g' - m'' * g'')
            freq_shift = mag_field.length() * con_mb / con_h *
                (i_sublvl_u * getLandeUpper(i_line) - i_sublvl_l * getLandeLower(i_line));

            // Calculate the frequency value of the current velocity channel in
            // relation to the peak of the line function
            f_doppler = CMathFunctions::Velo2Freq(velocity +
                    CMathFunctions::Freq2Velo(freq_shift, frequency), frequency) / doppler_width;

            // Calculate the line function value at the frequency
            // of the current velocity channel
            line_function = getLineShape_AB(f_doppler, voigt_a);

            // Multiply the line function value by PIsq for normalization
            mult_A = real(line_function) / PIsq;

            // Multiply the line function value by PIsq for normalization
            // and divide it by 2 to take the difference between the faddeeva
            // and Faraday-Voigt function into account
            mult_B = imag(line_function) / (PIsq * 2.0);

            // Use the correct propagation matrix and relative line strength that
            // depends on the current type of Zeeman transition (pi, sigma_-, sigma_+)
            switch(int(i_sublvl_l - i_sublvl_u))
            {
                case TRANS_SIGMA_P:
                    // Get the relative line strength from Zeeman file
                    line_strength = getLineStrengthSigmaP(i_line, i_sigma_p);

                    // Get the propagation matrix for extinction/emission
                    line_matrix += CMathFunctions::getPropMatrixASigmaP(cos_theta,
                            sin_theta, cos_2_phi, sin_2_phi, line_strength * mult_A);

                    // Get the propagation matrix for Faraday rotation
                    line_matrix += CMathFunctions::getPropMatrixBSigmaP(cos_theta,
                            sin_theta, cos_2_phi, sin_2_phi, line_strength * mult_B);

                    // Increase the sigma_+ counter to circle through the line strengths
                    i_sigma_p++;
                    break;
                case TRANS_PI:
                    // Get the relative line strength from Zeeman file
                    line_strength = getLineStrengthPi(i_line, i_pi);

                    // Get the propagation matrix for extinction/emission
                    line_matrix += CMathFunctions::getPropMatrixAPi(cos_theta,
                            sin_theta, cos_2_phi, sin_2_phi, line_strength * mult_A);

                    // Get the propagation matrix for Faraday rotation
                    line_matrix += CMathFunctions::getPropMatrixBPi(cos_theta,
                            sin_theta, cos_2_phi, sin_2_phi, line_strength * mult_B);

                    // Increase the pi counter to circle through the line strengths
                    i_pi++;
                    break;
                case TRANS_SIGMA_M:
                    // Get the relative line strength from Zeeman file
                    line_strength = getLineStrengthSigmaM(i_line, i_sigma_m);

                    // Get the propagation matrix for extinction/emission
                    line_matrix +=  CMathFunctions::getPropMatrixASigmaM(cos_theta,
                            sin_theta, cos_2_phi, sin_2_phi, line_strength * mult_A);

                    // Get the propagation matrix for Faraday rotation
                    line_matrix +=  CMathFunctions::getPropMatrixBSigmaM(cos_theta,
                            sin_theta, cos_2_phi, sin_2_phi, line_strength * mult_B);

                    // Increase the sigma_- counter to circle through the line strengths
                    i_sigma_m++;
                    break;
            }
        }
    }
    return line_matrix;
}

bool CGasSpecies::readGasParamaterFile(string _filename, uint id, uint max)
{
    uint line_counter, cmd_counter;
    uint pos_counter = 0;
    uint row_offset, i_col_partner, i_col_transition;
    CCommandParser ps;
    fstream reader(_filename.c_str());
    unsigned char ru[4] = {'|', '/', '-', '\\'};
    string line;
    dlist values;

    cout << CLR_LINE;

    if(reader.fail())
    {
        cout << "\nERROR: Cannot open gas_species catalog:" << endl;
        cout << _filename << endl;
        return false;
    }

    line_counter = 0;
    cmd_counter = 0;

    row_offset = 0;
    i_col_partner = 0;
    i_col_transition = 0;

    uint char_counter = 0;

    while(getline(reader, line))
    {
        line_counter++;

        if(line_counter % 20 == 0)
        {
            char_counter++;
            cout << "-> Reading gas species file nr. " << id << " of " << max << " : "
                << ru[(uint) char_counter % 4] << "           \r";
        }

        ps.formatLine(line);

        if(line.size() == 0)
            continue;

        if(cmd_counter != 0)
        {
            values = ps.parseValues(line);
            if(values.size() == 0)
                continue;
        }

        cmd_counter++;

        if(cmd_counter == 1)
            stringID = line;
        else if(cmd_counter == 2)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (gas species file)!" << endl;
                return false;
            }
            molecular_weight = values[0];
        }
        else if(cmd_counter == 3)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (gas species file)!" << endl;
                return false;
            }
            nr_of_energy_levels = uint(values[0]);

            energy_level = new double[nr_of_energy_levels];
            g_level = new double[nr_of_energy_levels];
            j_level = new double[nr_of_energy_levels];
        }
        else if(cmd_counter < 4 + nr_of_energy_levels && cmd_counter > 3)
        {
            if(values.size() < 4)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (gas species file)!" << endl;
                return false;
            }

            // ENERGIES(cm^-1)
            energy_level[pos_counter] = values[1];
            // WEIGHT
            g_level[pos_counter] = values[2];
            // J
            j_level[pos_counter] = values[3];

            pos_counter++;
        }
        else if(cmd_counter == 4 + nr_of_energy_levels)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (gas species file)!" << endl;
                return false;
            }
            nr_of_total_transitions = uint(values[0]);
            pos_counter = 0;

            trans_upper = new int[nr_of_total_transitions];
            trans_lower = new int[nr_of_total_transitions];
            trans_einstA = new double[nr_of_total_transitions];
            trans_freq = new double[nr_of_total_transitions];
            trans_inneregy = new double[nr_of_total_transitions];

            trans_einstB_u = new double[nr_of_total_transitions];
            trans_einstB_l = new double[nr_of_total_transitions];
        }
        else if(cmd_counter < 5 + nr_of_energy_levels + nr_of_total_transitions
                && cmd_counter > 4 + nr_of_energy_levels)
        {
            if(values.size() != 6)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (gas species file)!" << endl;
                return false;
            }

            for(uint i_line = 0; i_line < nr_of_transitions; i_line++)
            {
                if(transitions[i_line] == int(values[0] - 1))
                {
                    unique_transitions.push_back(i_line);
                    break;
                }
            }

            // UP (as index starting with 0)
            trans_upper[pos_counter] = int(values[1] - 1);
            // LOW (as index starting with 0)
            trans_lower[pos_counter] = int(values[2] - 1);
            // EINSTEINA(s^-1)
            trans_einstA[pos_counter] = values[3];
            // FREQ(GHz -> Hz)
            trans_freq[pos_counter] = values[4] * 1e9;
            // E_u(K)
            trans_inneregy[pos_counter] = values[5];

            trans_einstB_u[pos_counter] = values[3] * pow(con_c / (values[4] * 1e9), 2.0) /
                (2.0 * con_h * (values[4] * 1e9));

            trans_einstB_l[pos_counter] = g_level[int(values[1] - 1)] / g_level[int(values[2] - 1)]
                    * trans_einstB_u[pos_counter];

            pos_counter++;
        }
        else if(cmd_counter == 5 + nr_of_energy_levels + nr_of_total_transitions)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (gas species file)!" << endl;
                return false;
            }
            nr_of_col_partner = uint(values[0]);

            nr_of_collision_transition = new int[nr_of_col_partner];
            nr_of_col_temp = new int[nr_of_col_partner];
            orientation_H2 = new int[nr_of_col_partner];

            collision_temp = new double*[nr_of_col_partner];
            col_upper = new int*[nr_of_col_partner];
            col_lower = new int*[nr_of_col_partner];

            col_matrix = new double **[nr_of_col_partner];

            for(uint i = 0; i < nr_of_col_partner; i++)
            {
                collision_temp[i] = 0;
                col_upper[i] = 0;
                col_lower[i] = 0;

                nr_of_collision_transition[i] = 0;
                orientation_H2[i] = 0;
                nr_of_col_temp[i] = 0;
                col_matrix[i] = 0;
            }
        }
        else if(cmd_counter == 6 + nr_of_energy_levels + nr_of_total_transitions)
        {
            if(values[0] > 3 || values[0] < 1)
            {
                cout << "\nERROR: Line " << line_counter
                    << " wrong orientation of H2 collision partner (gas species file)!" << endl;
                return false;
            }
            orientation_H2[i_col_partner] = int(values[0]);
        }
        else if(cmd_counter == 7 + nr_of_energy_levels + nr_of_total_transitions + row_offset)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (gas species file)!" << endl;
                return false;
            }

            nr_of_collision_transition[i_col_partner] = int(values[0]);
            col_upper[i_col_partner] = new int[int(values[0])];
            col_lower[i_col_partner] = new int[int(values[0])];
            col_matrix[i_col_partner] = new double*[int(values[0])];

            for(uint i = 0; i < uint(values[0]); i++)
            {
                col_matrix[i_col_partner][i] = 0;
                col_upper[i_col_partner][i] = 0;
                col_lower[i_col_partner][i] = 0;
            }
        }
        else if(cmd_counter == 8 + nr_of_energy_levels + nr_of_total_transitions + row_offset)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (gas species file)!" << endl;
                return false;
            }

            nr_of_col_temp[i_col_partner] = int(values[0]);

            collision_temp[i_col_partner] = new double[int(values[0])];

            for(uint i = 0; i < uint(values[0]); i++)
                collision_temp[i_col_partner][i] = 0;
        }
        else if(cmd_counter == 9 + nr_of_energy_levels + nr_of_total_transitions + row_offset)
        {
            if(values.size() != uint(nr_of_col_temp[i_col_partner]))
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (gas species file)!" << endl;
                return false;
            }

            for(uint i = 0; i < uint(nr_of_col_temp[i_col_partner]); i++)
                collision_temp[i_col_partner][i] = values[i];
        }
        else if(cmd_counter < 10 + nr_of_energy_levels + nr_of_total_transitions +
                nr_of_collision_transition[i_col_partner] + row_offset &&
                cmd_counter > 9 + nr_of_energy_levels + row_offset)
        {
            if(values.size() != uint(nr_of_col_temp[i_col_partner] + 3))
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (gas species file)!" << endl;
                return false;
            }

            i_col_transition = cmd_counter - 10 - nr_of_energy_levels - nr_of_total_transitions - row_offset;
            col_upper[i_col_partner][i_col_transition] = int(values[1] - 1);
            col_lower[i_col_partner][i_col_transition] = int(values[2] - 1);

            col_matrix[i_col_partner][i_col_transition] = new double[nr_of_col_temp[i_col_partner]];

            for(uint i = 0; i < uint(nr_of_col_temp[i_col_partner]); i++)
                col_matrix[i_col_partner][i_col_transition][i] = values[3 + i] * 1.0e-6;
        }
        else if(i_col_partner < nr_of_col_partner)
        {
            if(values[0] > 3 || values[0] < 1)
            {
                cout << "\nERROR: Line " << line_counter
                    << " wrong orientation of H2 collision partner (gas species file)!" << endl;
                return false;
            }

            i_col_partner++;
            orientation_H2[i_col_partner] = int(values[0]);

            row_offset = cmd_counter - (6 + nr_of_energy_levels + nr_of_total_transitions);
        }
    }
    reader.close();

    return true;
}

bool CGasSpecies::readZeemanParamaterFile(string _filename)
{
    uint i_sigma_p_transition = 0, i_sigma_m_transition = 0;
    uint i_pi_transition = 0;
    bool splitting_exist = false;
    fstream reader(_filename.c_str());
    CCommandParser ps;
    string line;
    dlist values;

    uint sublevels_upper_nr = 0, sublevels_lower_nr = 0;
    uint offset_pi = 0, offset_sigma = 0;

    line_strength_pi = new double*[nr_of_transitions];
    line_strength_sigma_p = new double*[nr_of_transitions];
    line_strength_sigma_m = new double*[nr_of_transitions];

    nr_pi_transitions = new int[nr_of_transitions];
    nr_sigma_transitions = new int[nr_of_transitions];

    nr_sublevels_upper = new int[nr_of_transitions];
    nr_sublevels_lower = new int[nr_of_transitions];

    lande_upper = new double[nr_of_transitions];
    lande_lower = new double[nr_of_transitions];

    zeeman_transitions = new int[nr_of_transitions];

    for(uint i = 0; i < nr_of_transitions; i++)
    {
        line_strength_pi[i] = 0;
        line_strength_sigma_p[i] = 0;
        line_strength_sigma_m[i] = 0;

        nr_pi_transitions[i] = 0;
        nr_sigma_transitions[i] = 0;

        nr_sublevels_upper[i] = 0;
        nr_sublevels_lower[i] = 0;

        lande_upper[i] = 0;
        lande_lower[i] = 0;

        zeeman_transitions[i] = 0;
    }

    if(reader.fail())
    {
        cout << "\nERROR: Cannot open Zeeman splitting catalog:" << endl;
        cout << _filename << endl;
        return false;
    }

    uint line_counter = 0;
    uint cmd_counter = 0;

    while(getline(reader, line))
    {
        line_counter++;

        ps.formatLine(line);

        if(line.size() == 0)
            continue;

        if(cmd_counter != 0)
        {
            values = ps.parseValues(line);
            if(values.size() == 0)
                continue;
        }

        cmd_counter++;

        if(cmd_counter == 1)
        {
            if(line != stringID)
            {
                cout << "wrong Zeeman splitting catalog file chosen!" << endl;
                return false;
            }
        }
        else if(cmd_counter == 2)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (Zeeman file)!" << endl;
                return false;
            }
            gas_species_radius = values[0];
        }
        else if(cmd_counter == 3)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (Zeeman file)!" << endl;
                return false;
            }
            nr_zeeman_transitions = uint(values[0]);
        }
        else if(cmd_counter == 4)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (Zeeman file)!" << endl;
                return false;
            }
            bool not_unique = false;
            splitting_exist = false;
            for(uint i_line = 0; i_line < nr_of_transitions; i_line++)
            {
                if(transitions[i_line] == int(values[0] - 1))
                {
                    if(splitting_exist == false)
                    {
                        splitting_exist = true;
                        not_unique = false;
                    }
                    else
                        not_unique = true;
                }
                if(splitting_exist == true)
                    zeeman_transitions[i_line] = int(values[0] - 1);
            }
        }
        else if(cmd_counter == 5 && splitting_exist == true)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (Zeeman file)!" << endl;
                return false;
            }
            for(uint i_line = 0; i_line < nr_of_transitions; i_line++)
                if(transitions[i_line] == zeeman_transitions[i_line])
                    lande_upper[i_line] = values[0];
        }
        else if(cmd_counter == 6 && splitting_exist == true)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (Zeeman file)!" << endl;
                return false;
            }
            for(uint i_line = 0; i_line < nr_of_transitions; i_line++)
                if(transitions[i_line] == zeeman_transitions[i_line])
                    lande_lower[i_line] = values[0];
        }
        else if(cmd_counter == 7 && splitting_exist == true)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (Zeeman file)!" << endl;
                return false;
            }
            for(uint i_line = 0; i_line < nr_of_transitions; i_line++)
                if(transitions[i_line] == zeeman_transitions[i_line])
                    nr_sublevels_upper[i_line] = int(values[0]);
        }
        else if(cmd_counter == 7 && splitting_exist == false)
            sublevels_upper_nr = uint(values[0]);
        else if(cmd_counter == 8 && splitting_exist == true)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (Zeeman file)!" << endl;
                return false;
            }
            for(uint i_line = 0; i_line < nr_of_transitions; i_line++)
            {
                if(transitions[i_line] == zeeman_transitions[i_line])
                {
                    nr_sublevels_lower[i_line] = int(values[0]);
                    nr_pi_transitions[i_line] = min(nr_sublevels_upper[i_line], nr_sublevels_lower[i_line]);

                    if(nr_sublevels_upper[i_line] != nr_sublevels_lower[i_line])
                        nr_sigma_transitions[i_line] = min(nr_sublevels_upper[i_line], nr_sublevels_lower[i_line]);
                    else
                        nr_sigma_transitions[i_line] = nr_sublevels_upper[i_line] - 1;

                    offset_pi = nr_pi_transitions[i_line];
                    offset_sigma = nr_sigma_transitions[i_line];

                    line_strength_pi[i_line] = new double[nr_pi_transitions[i_line]];

                    for(int i = 0; i < nr_pi_transitions[i_line]; i++)
                        line_strength_pi[i_line][i] = 0;

                    line_strength_sigma_p[i_line] = new double[nr_sigma_transitions[i_line]];
                    line_strength_sigma_m[i_line] = new double[nr_sigma_transitions[i_line]];

                    for(int i = 0; i < nr_sigma_transitions[i_line]; i++)
                    {
                        line_strength_sigma_p[i_line][i] = 0;
                        line_strength_sigma_m[i_line][i] = 0;
                    }
                }
            }
            i_pi_transition = 0;
            i_sigma_p_transition = 0;
            i_sigma_m_transition = 0;
        }
        else if(cmd_counter == 8 && splitting_exist == false)
        {
            sublevels_lower_nr = uint(values[0]);
            offset_pi = min(sublevels_upper_nr, sublevels_lower_nr);
            if(sublevels_upper_nr != sublevels_lower_nr)
                offset_sigma = min(sublevels_upper_nr, sublevels_lower_nr);
            else
                offset_sigma = sublevels_upper_nr - 1;
        }
        else if(cmd_counter <= 8 + offset_pi && cmd_counter > 8 && splitting_exist == true)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (Zeeman file)!" << endl;
                return false;
            }
            for(uint i_line = 0; i_line < nr_of_transitions; i_line++)
                if(transitions[i_line] == zeeman_transitions[i_line])
                    line_strength_pi[i_line][i_pi_transition] = values[0];
            i_pi_transition++;
        }
        else if(cmd_counter <= 8 + offset_pi + offset_sigma && cmd_counter > 8 + offset_pi && splitting_exist == true)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (Zeeman file)!" << endl;
                return false;
            }
            for(uint i_line = 0; i_line < nr_of_transitions; i_line++)
                if(transitions[i_line] == zeeman_transitions[i_line])
                    line_strength_sigma_p[i_line][i_sigma_p_transition] = values[0];
            i_sigma_p_transition++;
        }
        else if(cmd_counter <= 8 + offset_pi + 2 * offset_sigma &&
                cmd_counter > 8 + offset_pi + offset_sigma && splitting_exist == true)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (Zeeman file)!" << endl;
                return false;
            }
            for(uint i_line = 0; i_line < nr_of_transitions; i_line++)
                if(transitions[i_line] == zeeman_transitions[i_line])
                    line_strength_sigma_m[i_line][i_sigma_m_transition] = values[0];
            i_sigma_m_transition++;
        }
        else if(cmd_counter == 9 + offset_pi + 2 * offset_sigma)
        {
            if(values.size() != 1)
            {
                cout << "\nERROR: Line " << line_counter << " wrong amount of numbers (Zeeman file)!" << endl;
                return false;
            }
            bool not_unique = false;
            splitting_exist = false;
            for(uint i_line = 0; i_line < nr_of_transitions; i_line++)
            {
                if(transitions[i_line] == int(values[0] - 1))
                {
                    if(splitting_exist == false)
                    {
                        splitting_exist = true;
                        not_unique = false;
                    }
                    else
                        not_unique = true;
                }
                if(splitting_exist == true)
                    zeeman_transitions[i_line] = int(values[0] - 1);
            }
            cmd_counter = 4;
        }
    }
    reader.close();

    for(uint i_line = 0; i_line < nr_of_transitions; i_line++)
        if(lande_upper[i_line] == 0)
        {
            cout << SEP_LINE;
            cout << "\nERROR: For transition number " << uint(getTransition(i_line) + 1)
                << " exists no Zeeman splitting data" << endl;
            return false;
        }

    zeeman_splitting = true;
    return true;
}

void CGasSpecies::calcLineBroadening(CGridBasic * grid)
{
    long max_cells = grid->getMaxDataCells();
#pragma omp parallel for schedule(dynamic)
    for(long i_cell = 0; i_cell < long(max_cells); i_cell++)
    {
        cell_basic * cell = grid->getCellFromIndex(i_cell);

        // Get necessary quantities from the current cell
        double temp_gas = grid->getGasTemperature(cell);
        double dens_gas = grid->getGasNumberDensity(cell);
        double dens_species = getNumberDensity(grid, cell);
        double turbulent_velocity = grid->getTurbulentVelocity(cell);

        for(uint i_line = 0; i_line < nr_of_transitions; i_line++)
        {
            uint i_trans = getTransition(i_line);
            double frequency = getTransitionFrequency(i_trans);

            double gauss_a = getGaussA(temp_gas, turbulent_velocity);
            double Gamma = 0, doppler_width = 0, voigt_a = 0;
            if(getZeemanSplitting())
            {
                Gamma = getGamma(i_trans, dens_gas, dens_species, temp_gas, turbulent_velocity);
                doppler_width = frequency / (con_c * gauss_a);
                voigt_a = Gamma / (4 * PI * doppler_width);
            }
            grid->setLineBroadening(cell, i_line, gauss_a, Gamma, doppler_width, voigt_a);
        }
    }
}

bool CGasMixture::createGasSpecies(parameters & param)
{
    nr_of_species = param.getNrOfGasSpecies();
    single_species = new CGasSpecies[nr_of_species];

    for(uint i_species = 0; i_species < nr_of_species; i_species++)
    {
        single_species[i_species].setAbundance(param.getGasSpeciesAbundance(i_species));
        single_species[i_species].setLevelPopType(param.getGasSpeciesLevelPopType(i_species));
        single_species[i_species].setNrTransitions(param.getNrOfGasSpeciesTransitions(i_species));
        single_species[i_species].setTransitions(param.getGasSpeciesTransitions(i_species));

        string path = param.getGasSpeciesCatalogPath(i_species);
        if(!single_species[i_species].readGasParamaterFile(path, i_species, nr_of_species))
            return false;

        if(param.getZeemanCatalog(i_species) != "")
            if(!single_species[i_species].readZeemanParamaterFile(param.getZeemanCatalog(i_species)))
                return false;

        if((single_species[i_species].getLevelPopType() == POP_FEP || 
                single_species[i_species].getLevelPopType() == POP_LVG) &&
                single_species[i_species].getNrCollisionPartner() == 0)
        {
            cout << "\nERROR: FEP and LVG level population approximations require a gas parameters file \n"
                "       with collisional data (e.g. from Leiden Atomic and Molecular Database)" << endl;
            return false;
        }
    }

    setKeplerStarMass(param.getKeplerStarMass());
    return true;
}

bool CGasMixture::calcLevelPopulation(CGridBasic * grid, uint i_species)
{
    uint lvl_pop_type = getLevelPopType(i_species);
    double kepler_star_mass = getKeplerStarMass();

    switch(lvl_pop_type)
    {
        case POP_LTE:
            if(!single_species[i_species].calcLTE(grid))
                return false;
            break;
        case POP_FEP:
            if(!single_species[i_species].calcFEP(grid))
                return false;
            break;
        case POP_LVG:
            if(!single_species[i_species].calcLVG(grid, kepler_star_mass))
                return false;
            break;
        default:
            return false;
            break;
    }
    return true;
}

void CGasMixture::printParameter(parameters & param, CGridBasic * grid)
{
    cout << CLR_LINE;
    cout << "Gas parameters                             " << endl;
    cout << SEP_LINE;
    cout << "- Velocity field                : ";
    if(getKeplerStarMass() > 0)
        cout << "kepler rotation, M_star: " << getKeplerStarMass()
            << " [M_sun]\n" << "    HINT: only available with one central star" << endl;
    else
        cout << "velocity field of the grid is used" << endl;
    cout << "- Turbulent Velocity            : ";
    if(param.getTurbulentVelocity() > 0)
        cout << param.getTurbulentVelocity() << " [m/s]" << endl;
    else
        cout << "turbulent velocity of the grid is used" << endl;

    for(uint i_species = 0; i_species < nr_of_species; i_species++)
    {
        cout << SEP_LINE;
        cout << "Gas species " << (i_species + 1) << " (" << getGasSpeciesName(i_species) << ")" << endl;

        if(single_species[i_species].getNrOfTransitions() == 0)
        {
            cout << "\nWARNING: No spectral lines selected!" << endl;
            return;
        }

        stringstream transition_str, vel_channels_str, max_vel_str;
        dlist line_ray_detectors = param.getLineRayDetector(i_species);
        for(uint i = 0; i < line_ray_detectors.size(); i += NR_OF_LINE_DET)
        {
            uint pos = i / NR_OF_LINE_DET;

            transition_str << uint(line_ray_detectors[i] + 1);
            max_vel_str << line_ray_detectors[i + 2];
            vel_channels_str << uint(line_ray_detectors[i + NR_OF_LINE_DET - 1]);
            if(i < line_ray_detectors.size() - NR_OF_LINE_DET)
            {
                transition_str << ", ";
                max_vel_str << ", ";
                vel_channels_str << ", ";
            }
        }
        cout << "- Line transition(s)            : " << transition_str.str() << endl;
        cout << "- Number of velocity channels      : " << vel_channels_str.str() << endl;
        cout << "- Velocity limit(s)             : " << max_vel_str.str() << " [m/s]" << endl;
        cout << "- Level population              : ";
        uint lvl_pop_type = getLevelPopType(i_species);
        switch(lvl_pop_type)
        {
            case POP_LTE:
                cout << "LTE" << endl;
                break;
            case POP_FEP:
                cout << "FEP" << endl;
                break;
            case POP_LVG:
                cout << "LVG" << endl;
                break;
            default:
                cout << "\nERROR: UNKNOWN!" << endl;
        }

        if(getZeemanSplitting(i_species) == true)
            cout << "- Particle radius (collisions)  : " << getCollisionRadius(i_species) << " [m]" << endl;

        cout << "- Molecular weight              : " << getMolecularWeight(i_species) << endl;
        double ab = getAbundance(i_species);
        if(ab > 0)
            cout << "- Abundance                     : " << ab << endl;
        else
        {
            double dens_species, min_ab = 1e200, max_ab = 0;
            for(long i_cell = 0; i_cell < grid->getMaxDataCells(); i_cell++)
            {
                cell_basic * cell = grid->getCellFromIndex(i_cell);
                double tmp_ab = grid->getCellAbundance(cell, uint(-ab - 1));
                if(tmp_ab < min_ab)
                    min_ab = tmp_ab;
                if(tmp_ab > max_ab)
                    max_ab = tmp_ab;
            }
            cout << "- Abundance from grid ID nr.    : " << int(-ab) << endl;
            cout << "                      (min,max) : [" << min_ab << ", " << max_ab << "] [m^-3]" << endl;
        }

        double total_species_mass = 0;
        for(long i_cell = 0; i_cell < grid->getMaxDataCells(); i_cell++)
        {
            cell_basic * cell = grid->getCellFromIndex(i_cell);
            total_species_mass += getMassDensity(grid, cell, i_species) * grid->getVolume(cell);
        }
        cout << "- Total mass                    : " << total_species_mass / M_sun << " [M_sun]" << endl;

        for(uint i = 0; i < getUniqueTransitions(i_species).size(); i++)
        {
            uint i_line = getUniqueTransitions(i_species, i);

            cout << SEP_LINE;
            cout << "Line transition " << uint(getTransition(i_species, i_line) + 1)
                << " (gas species " << (i_species + 1) << ")" << endl;
            cout << "- Transition frequency          : "
                << getTransitionFrequencyFromIndex(i_species, i_line) << " [Hz]" << endl;
            cout << "- Transition wavelength         : "
                    << (con_c / getTransitionFrequencyFromIndex(i_species, i_line)) << " [m]" << endl;
            if(getZeemanSplitting(i_species) == true)
            {
                cout << CLR_LINE;
                cout << "Zeeman splitting parameters                " << endl;
                cout << "- Lande factor of upper level   : " << getLandeUpper(i_species, i_line) << endl;
                cout << "- Lande factor of lower level   : " << getLandeLower(i_species, i_line) << endl;
                cout << "- Sublevels in upper level      : " << getNrSublevelsUpper(i_species, i_line) << endl;
                cout << "- Sublevels in lower level      : " << getNrSublevelsLower(i_species, i_line) << endl;

                uint i_pi = 0, i_sigma_p = 0, i_sigma_m = 0;

                cout << "- Sigma+ line strength\t\tm(upper)\t\tm(lower)" << endl;
                for(float i_sublvl_u = -getMaxMUpper(i_species, i_line);
                    i_sublvl_u <= getMaxMUpper(i_species, i_line); i_sublvl_u++)
                {
                    float i_sublvl_l = i_sublvl_u + 1;
                    if(abs(i_sublvl_l) <= getMaxMLower(i_species, i_line))
                    {
                        char LineStrengthTmp[16];
#ifdef WINDOWS
                        _snprintf_s(LineStrengthTmp, sizeof(LineStrengthTmp), "%.3f",
                            getLineStrengthSigmaP(i_species, i_line, i_pi));
#else
                        snprintf(LineStrengthTmp, sizeof(LineStrengthTmp), "%.3f",
                            getLineStrengthSigmaP(i_species, i_line, i_pi));
#endif

                        cout << "\t" << LineStrengthTmp << "\t\t\t   " << float(i_sublvl_u)
                            << "\t\t\t   " << float(i_sublvl_l) << endl;
                        i_sigma_p++;
                    }
                }
                cout << "- Pi line strength\t\tm(upper)\t\tm(lower)" << endl;
                for(float i_sublvl_u = -getMaxMUpper(i_species, i_line);
                    i_sublvl_u <= getMaxMUpper(i_species, i_line); i_sublvl_u++)
                {
                    float i_sublvl_l = i_sublvl_u;
                    if(abs(i_sublvl_l) <= getMaxMLower(i_species, i_line))
                    {
                        char LineStrengthTmp[16];
#ifdef WINDOWS
                        _snprintf_s(LineStrengthTmp, sizeof(LineStrengthTmp), "%.3f",
                            getLineStrengthPi(i_species, i_line, i_pi));
#else
                        snprintf(LineStrengthTmp, sizeof(LineStrengthTmp), "%.3f",
                            getLineStrengthPi(i_species, i_line, i_pi));
#endif
                        cout << "\t" << LineStrengthTmp << "\t\t\t   " << float(i_sublvl_u)
                            << "\t\t\t   " << float(i_sublvl_l) << endl;
                        i_pi++;
                    }
                }
                cout << "- Sigma- line strength\t\tm(upper)\t\tm(lower)" << endl;
                for(float i_sublvl_u = -getMaxMUpper(i_species, i_line);
                    i_sublvl_u <= getMaxMUpper(i_species, i_line); i_sublvl_u++)
                {
                    float i_sublvl_l = i_sublvl_u - 1;
                    if(abs(i_sublvl_l) <= getMaxMLower(i_species, i_line))
                    {
                        char LineStrengthTmp[16];
#ifdef WINDOWS
                        _snprintf_s(LineStrengthTmp, sizeof(LineStrengthTmp), "%.3f",
                            getLineStrengthSigmaM(i_species, i_line, i_pi));
#else
                        snprintf(LineStrengthTmp, sizeof(LineStrengthTmp), "%.3f",
                            getLineStrengthSigmaM(i_species, i_line, i_pi));
#endif

                        cout << "\t" << LineStrengthTmp << "\t\t\t   " << float(i_sublvl_u)
                            << "\t\t\t   " << float(i_sublvl_l) << endl;
                        i_sigma_m++;
                    }
                }
            }
        }
    }
    cout << SEP_LINE;
}
