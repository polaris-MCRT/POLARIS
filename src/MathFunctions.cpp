/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "MathFunctions.hpp"
#include "Faddeeva.hh"

bool CMathFunctions::isPowerOfTwo(int num)
{
    return ((num & (num - 1)) == 0);
}

Vector3D CMathFunctions::SolveEqSys(Matrix2D & inA, Vector3D b)
{
    Vector3D res;
    double sum;
    Matrix2D A = inA;
    uint M = A.get_m();
    double * a = new double[M];

    A.printMatrix();

    for(uint i = 1; i <= M; i++)
        a[i - 1] = 0;

    CholDec(A);

    A.printMatrix();

    for(uint i = 1; i <= M; i++)
    {
        sum = b(i - 1);
        for(uint k = i - 1; k >= 1; k--)
            sum -= A(i - 1, k - 1) * a[k - 1];

        a[i - 1] = sum / A(i - 1, i - 1);
    }

    for(uint i = M; i >= 1; i--)
    {
        sum = a[i - 1];
        for(uint k = i + 1; k <= M; k++)
            sum -= A(k - 1, i - 1) * a[k - 1];

        a[i - 1] = sum / A(i - 1, i - 1);
    }

    // res.setX(a[0]);
    // res.setY(a[1]);
    // res.setZ(a[2]);

    delete[] a;
    return res;
}

Vector3D CMathFunctions::gauss(Matrix2D inA, Vector3D inb)
{
    int n = 3;
    double x, sum;
    double * b = new double[n];

    b[0] = inb.X();
    b[1] = inb.Y();
    b[2] = inb.Z();

    for(int k = 0; k < n - 1; k++)
    {

        for(int i = k + 1; i < n; i++)
        {
            x = inA(i, k) / inA(k, k);

            for(int j = k + 1; j < n; j++)
            {
                inA(i, j) = inA(i, j) - inA(k, j) * x;
            }

            b[i] = b[i] - b[k] * x;
        }
    }

    // Resubstitution
    b[n - 1] = b[n - 1] / inA(n - 1, n - 1);
    for(int i = n - 2; i >= 0; i--)
    {
        sum = b[i];

        for(int j = i + 1; j < n; j++)
        {
            sum = sum - inA(i, j) * b[j];
        }

        b[i] = sum / inA(i, i);
    }

    Vector3D res(b[0], b[1], b[2]);
    delete[] b;

    return res;
}

void CMathFunctions::gauss(Matrix2D & A, double * b, int n)
{
    double x, sum;

    /*b[0] = inb.X();
    b[1] = inb.Y();
    b[2] = inb.Z();*/

    for(int k = 0; k < n - 1; k++)
    {

        for(int i = k + 1; i < n; i++)
        {
            x = A(i, k) / A(k, k);

            for(int j = k + 1; j < n; j++)
            {
                A(i, j) = A(i, j) - A(k, j) * x;
            }

            b[i] = b[i] - b[k] * x;
        }
    }

    // Resubstitution
    b[n - 1] = b[n - 1] / A(n - 1, n - 1);
    for(int i = n - 2; i >= 0; i--)
    {
        sum = b[i];

        for(int j = i + 1; j < n; j++)
        {
            sum = sum - A(i, j) * b[j];
        }

        b[i] = sum / A(i, i);
    }

    // return b;
}

double CMathFunctions::BE_Mass(double r, double rho0, double rc)
{
    double res = PIx4 * rho0 * rc * rc * (r - rc * atan(r / rc));
    return res;
}

void CMathFunctions::CholDec(Matrix2D & A)
{
    uint M = A.get_m();
    for(uint k = 1; k <= M; k++)
    {
        A(k - 1, k - 1) = sqrt(A(k - 1, k - 1));

        for(uint i = k + 1; i <= M; i++)
        {
            A(i - 1, k - 1) = A(i - 1, k - 1) / A(k - 1, k - 1);
            for(uint j = k + 1; j <= M; j++)
                A(i - 1, j - 1) = A(i - 1, j - 1) - A(i - 1, k - 1) * A(j - 1, k - 1);
        }
    }
}

void CMathFunctions::setHeader(uint pos,
                uint seqID,
                uint sourceID,
                double wavelength,
                double _rot_angle1,
                double _rot_angle2)
{
    stat_data(pos, 0) = double(seqID) + 1;
    stat_data(pos, 1) = double(sourceID) + 1;
    stat_data(pos, 2) = wavelength;
    stat_data(pos, 3) = _rot_angle1;
    stat_data(pos, 4) = _rot_angle2;
}

void CMathFunctions::updateSatistics(uint pos, StokesVector st)
{
    if(st.I() == 0)
        return;

    double Pl = sqrt(st.Q() * st.Q() + st.U() * st.U()) / st.I();
    double Pc = st.V() / st.I();

    stat_data(pos, 5)++;

    stat_data(pos, 7) += st.I();
    if(stat_data(pos, 8) > st.I())
        stat_data(pos, 8) = st.I();
    if(stat_data(pos, 9) < st.I())
        stat_data(pos, 9) = st.I();

    stat_data(pos, 10) += st.Q();
    if(stat_data(pos, 11) > st.Q())
        stat_data(pos, 11) = st.Q();
    if(stat_data(pos, 12) < st.Q())
        stat_data(pos, 12) = st.Q();

    stat_data(pos, 13) += st.U();
    if(stat_data(pos, 14) > st.U())
        stat_data(pos, 14) = st.U();
    if(stat_data(pos, 15) < st.U())
        stat_data(pos, 15) = st.U();

    stat_data(pos, 16) += st.V();
    if(stat_data(pos, 17) > st.V())
        stat_data(pos, 17) = st.V();
    if(stat_data(pos, 18) < st.V())
        stat_data(pos, 18) = st.V();

    stat_data(pos, 19) += st.T();
    if(stat_data(pos, 20) > st.T())
        stat_data(pos, 20) = st.T();
    if(stat_data(pos, 21) < st.T())
        stat_data(pos, 21) = st.T();

    stat_data(pos, 22) += Pl;
    if(stat_data(pos, 23) > Pl)
        stat_data(pos, 23) = Pl;
    if(stat_data(pos, 24) < Pl)
        stat_data(pos, 24) = Pl;

    stat_data(pos, 25) += Pc;
    if(stat_data(pos, 26) > Pc)
        stat_data(pos, 26) = Pc;
    if(stat_data(pos, 27) < Pc)
        stat_data(pos, 27) = Pc;
}

int CMathFunctions::insertInList(double val, dlist & list)
{
    uint N = uint(list.size());
    dlist::iterator it;

    if(N == 0)
    {
        list.push_back(val);
        return MAX_UINT;
    }

    if(N == 1)
    {
        if(list[0] == val)
            return 0;
    }

    if(N == 2)
    {
        if(val > list[0] && val < list[1])
        {
            it = list.begin();
            list.insert(it + 1, val);
            return MAX_UINT;
        }
    }

    if(val < list[0])
    {
        it = list.begin();
        list.insert(it, val);
        return MAX_UINT;
    }

    if(val > list[N - 1])
    {
        list.push_back(val);
        return MAX_UINT;
    }

    it = list.begin();
    uint index = biListIndexSearch(val, list);

    if(list[index + 1] == val)
        return index + 1;
    else
        list.insert(it + index + 1, val);

    return MAX_UINT;
}

uint CMathFunctions::biListIndexSearch(double val, const dlist & list)
{
    uint N = uint(list.size());

    if(val < list[0] || val > list[N - 1])
        return MAX_UINT;

    if(val == list[0])
        return 0;

    uint min = lower_bound(list.begin(), list.end(), val) - list.begin() - 1;

    return min;
}

uint CMathFunctions::biListIndexSearchRec(double val, const dlist & list)
{
    uint N = uint(list.size());

    if(val < list[0] || val > list[N - 1])
        return MAX_UINT;

    uint min = upper_bound(list.begin(), list.end(), val) - list.begin() - 1;

    return min;
}

uint CMathFunctions::biListIndexSearchRec(double val, const double * list, uint N)
{
    if(val < list[0] || val > list[N - 1])
        return MAX_UINT;

    uint min = upper_bound(list, list+N, val) - list - 1;

    return min;
}

uint CMathFunctions::biListIndexSearch(double val, const double * list, uint N)
{
    if(val < list[0] || val > list[N - 1])
        return MAX_UINT;

    if(val == list[0])
        return 0;

    uint min = lower_bound(list, list+N, val) - list - 1;

    return min;
}

double CMathFunctions::interpolate(double x_min, double x_max, double y_min, double y_max, double x_ipol)
{
    double ipol2_result;
    if(x_min == x_max)
        ipol2_result = y_min;
    else
        ipol2_result = (y_max - y_min) * (x_ipol - x_min) / (x_max - x_min) + y_min;

    return ipol2_result;
}

uint CMathFunctions::inSphereTest(Vector3D a, Vector3D b, Vector3D c, Vector3D d, Vector3D e)
{
    return 0;
}

uint CMathFunctions::orientationTest(Vector3D a, Vector3D b, Vector3D c, Vector3D d)
{
    return 0;
}

double CMathFunctions::det5x5(double mat[5][5], int N)
{
    double d = 0;

    if(N == 2)
        return ((mat[0][0] * mat[1][1]) - (mat[1][0] * mat[0][1]));
    else
    {
        double submat[5][5];

        for(int c = 0; c < N; c++)
        {
            int subi = 0;

            for(int i = 1; i < N; i++)
            {
                int subj = 0;
                for(int j = 0; j < N; j++)
                {
                    if(j == c)
                        continue;

                    submat[subi][subj] = mat[i][j];
                    subj++;
                }
                subi++;
            }

            d += pow(double(-1), double(c)) * mat[0][c] * det5x5(submat, N - 1);
        }
    }

    return d;
}

double CMathFunctions::det4x4(double mat[4][4], int N)
{
    double d = 0;

    if(N == 2)
        return ((mat[0][0] * mat[1][1]) - (mat[1][0] * mat[0][1]));
    else
    {
        double submat[4][4];

        for(int c = 0; c < N; c++)
        {
            int subi = 0;
            for(int i = 1; i < N; i++)
            {
                int subj = 0;
                for(int j = 0; j < N; j++)
                {
                    if(j == c)
                        continue;

                    submat[subi][subj] = mat[i][j];
                    subj++;
                }
                subi++;
            }

            d += pow(double(-1), double(c)) * mat[0][c] * det4x4(submat, N - 1);
        }
    }

    return d;
}

Vector3D CMathFunctions::calcKeplerianVelocity(Vector3D pos, double stellar_mass)
{
    Vector3D velo;

    double r = sqrt(pos.X() * pos.X() + pos.Y() * pos.Y());
    if(r > 0)
    {
        double kep_const = sqrt(con_G * stellar_mass * M_sun / r);
        velo.setX(-1.0 * pos.Y() / r * kep_const);
        velo.setY(pos.X() / r * kep_const);
    }

    return velo;
}

double CMathFunctions::getClosestLinePoint(Vector3D line_start, Vector3D line_end, Vector3D point)
{
    Vector3D line_diff = line_end - line_start;
    double line_length_sq = pow(line_diff.length(), 2);

    Vector3D line_to_point = point - line_start;

    double perc_along_line = (line_diff * line_to_point) / line_length_sq;

    if(perc_along_line < 0.0)
    {
        return 0.0;
    }
    else if(perc_along_line > 1.0)
    {
        return line_diff.length();
    }
    else
    {
        return perc_along_line * line_diff.length();
    }
}

double CMathFunctions::maxValue(double v1, double v2)
{
    if(v1 > v2)
        return v1;

    return v2;
}

uint CMathFunctions::maxValue(uint v1, uint v2)
{
    if(v1 > v2)
        return v1;

    return v2;
}

double CMathFunctions::grad2rad(double grad)
{
    return grad * PI / 180.0;
}

double CMathFunctions::rad2grad(double rad)
{
    return rad * 180.0 / PI;
}

bool CMathFunctions::areSame(double a, double b)
{
    if(fabs(a - b) < 1e-5 * fabs(a + b))
        return true;
    return false;
}

double CMathFunctions::generateGaussianNoise(const double & variance)
{
    bool hasSpare = false;
    double rand1, rand2;

    if(hasSpare)
    {
        hasSpare = false;
        return sqrt(variance * rand1) * sin(rand2);
    }

    hasSpare = true;

    rand1 = rand() / ((double)RAND_MAX);
    if(rand1 < 1e-100)
        rand1 = 1e-100;
    rand1 = -2 * log(rand1);
    rand2 = (rand() / ((double)RAND_MAX)) * PIx2;

    return sqrt(variance * rand1) * cos(rand2);
}

double CMathFunctions::diffSeries(double y, double upper_limit)
{
    double res = 0;
    double y_n;

    for(int n = 1; n <= upper_limit; n++)
    {
        y_n = pow(y, double(n));
        res += 2 * pow(double(-1), double(n + 1)) * y_n * y_n;
    }

    return res;
}

double CMathFunctions::getErfi(double x)
{
    double res = Faddeeva::erfi(x);
    if(res < 0)
        return 0;

    return res;
}

void CMathFunctions::sleep(int milliseconds)
{
#ifdef WINDOWS
    Sleep(milliseconds);
#else

    usleep(milliseconds * 1000);
#endif
}

double CMathFunctions::integ(const double * x, const double * y, uint xlow, uint xup)
{
    double res = 0;
    if(xlow != xup)
    {
        for(uint i = xlow + 1; i <= xup; i++)
            res += (x[i] - x[i - 1]) * y[i - 1] + 0.5 * (x[i] - x[i - 1]) * (y[i] - y[i - 1]);
    }

    return res;
}

double CMathFunctions::integ(const dlist & x, double * y, uint xlow, uint xup)
{
    double res = 0;
    if(xlow != xup)
        for(uint i = xlow + 1; i <= xup; i++)
            res += (x[i] - x[i - 1]) * y[i - 1] + 0.5 * (x[i] - x[i - 1]) * (y[i] - y[i - 1]);

    return res;
}

// for debugging reasons only
double CMathFunctions::integ1(dlist x, double * y, uint xlow, uint xup)
{
    double res = 0;
    if(xlow != xup)
        for(uint i = xlow + 1; i <= xup; i++)
        {
            res += (x[i] - x[i - 1]) * y[i - 1] + 0.5 * (x[i] - x[i - 1]) * (y[i] - y[i - 1]);

            cout << i << " " << res << " " << x[i] << " " << y[i] << endl;
        }

    return res;
}

double CMathFunctions::integ(const dlist & x, const dlist & y, uint xlow, uint xup)
{
    double res = 0;
    if(xlow != xup)
        for(uint i = xlow + 1; i <= xup; i++)
            res += (x[i] - x[i - 1]) * y[i - 1] + 0.5 * (x[i] - x[i - 1]) * (y[i] - y[i - 1]);

    return res;
}

double CMathFunctions::integ(const dlist & x, const uilist & y, uint xlow, uint xup)
{
    double res = 0;
    if(xlow != xup)
        for(uint i = xlow + 1; i <= xup; i++)
            res += (x[i] - x[i - 1]) * y[i - 1] + 0.5 * (x[i] - x[i - 1]) * (y[i] - y[i - 1]);

    return res;
}

double CMathFunctions::integ(const double * x, dlist & y, uint xlow, uint xup)
{
    double res = 0;
    if(xlow != xup)
        for(uint i = xlow + 1; i <= xup; i++)
            res += (x[i] - x[i - 1]) * y[i - 1] + 0.5 * (x[i] - x[i - 1]) * (y[i] - y[i - 1]);

    return res;
}

double CMathFunctions::integ(double * x, dlist & y, uint xlow, uint xup)
{
    double res = 0;
    if(xlow != xup)
        for(uint i = xlow + 1; i <= xup; i++)
            res += (x[i] - x[i - 1]) * y[i - 1] + 0.5 * (x[i] - x[i - 1]) * (y[i] - y[i - 1]);

    return res;
}

double CMathFunctions::integ_dust_size(const double * a_eff,
                                        const double * quantity,
                                        uint nr_of_dust_species,
                                        double a_min,
                                        double a_max)
{
    double res = 0;
    double dx, y0, m;

    // If maximum grain size is smaller than the minimum grain size -> return zero
    if(a_max < a_min)
        return res;

    // Return the quantity value if only one size is used
    if(nr_of_dust_species == 1)
        return quantity[0];

    for(uint i = 1; i < nr_of_dust_species; i++)
    {
        m = (quantity[i] - quantity[i - 1]) / (a_eff[i] - a_eff[i - 1]);
        if(a_eff[i - 1] >= a_min && a_eff[i] < a_max)
        {
            dx = a_eff[i] - a_eff[i - 1];
            y0 = quantity[i - 1];
        }
        else if(a_eff[i - 1] < a_min && a_eff[i] >= a_max)
        {
            dx = a_max - a_min;
            y0 = quantity[i - 1] + (a_min - a_eff[i - 1]) * m;
            if(dx == 0)
                return y0;
        }
        else if(a_eff[i] >= a_min && a_eff[i - 1] < a_min)
        {
            dx = a_eff[i] - a_min;
            y0 = quantity[i];
            m *= -1;
        }
        else if(a_eff[i - 1] < a_max && a_eff[i] >= a_max)
        {
            dx = a_max - a_eff[i - 1];
            y0 = quantity[i - 1];
        }
        else
            continue;

        res += y0 * dx + 0.5 * dx * dx * m;
    }

    return res;
}

StokesVector CMathFunctions::integ_dust_size(const double * a_eff,
                                            const StokesVector * stokes,
                                            uint nr_of_dust_species,
                                            double a_min,
                                            double a_max)
{
    StokesVector final_stokes;
    double m_I, m_Q, m_U, m_V;
    double y0_I, y0_Q, y0_U, y0_V;
    double dx;

    // If maximum grain size is smaller than the minimum grain size -> return zero
    if(a_max < a_min)
        return final_stokes;

    // Return the quantity value if only one size is used
    if(nr_of_dust_species == 1)
        return stokes[0];

    for(uint i = 1; i < nr_of_dust_species; i++)
    {
        m_I = (stokes[i].I() - stokes[i - 1].I()) / (a_eff[i] - a_eff[i - 1]);
        m_Q = (stokes[i].Q() - stokes[i - 1].Q()) / (a_eff[i] - a_eff[i - 1]);
        m_U = (stokes[i].U() - stokes[i - 1].U()) / (a_eff[i] - a_eff[i - 1]);
        m_V = (stokes[i].V() - stokes[i - 1].V()) / (a_eff[i] - a_eff[i - 1]);

        if(a_eff[i - 1] >= a_min && a_eff[i] < a_max)
        {
            dx = a_eff[i] - a_eff[i - 1];
            y0_I = stokes[i - 1].I();
            y0_Q = stokes[i - 1].Q();
            y0_U = stokes[i - 1].U();
            y0_V = stokes[i - 1].V();
        }
        else if(a_eff[i - 1] < a_min && a_eff[i] >= a_max)
        {
            dx = a_max - a_min;
            y0_I = stokes[i - 1].I() + (a_min - a_eff[i - 1]) * m_I;
            y0_Q = stokes[i - 1].Q() + (a_min - a_eff[i - 1]) * m_Q;
            y0_U = stokes[i - 1].U() + (a_min - a_eff[i - 1]) * m_U;
            y0_V = stokes[i - 1].V() + (a_min - a_eff[i - 1]) * m_V;
            if(dx == 0)
                return StokesVector(y0_I, y0_Q, y0_U, y0_V);
        }
        else if(a_eff[i] >= a_min && a_eff[i - 1] < a_min)
        {
            dx = a_eff[i] - a_min;
            y0_I = stokes[i].I();
            y0_Q = stokes[i].Q();
            y0_U = stokes[i].U();
            y0_V = stokes[i].V();
            m_I *= -1;
            m_Q *= -1;
            m_U *= -1;
            m_V *= -1;
        }
        else if(a_eff[i - 1] < a_max && a_eff[i] >= a_max)
        {
            dx = a_max - a_eff[i - 1];
            y0_I = stokes[i - 1].I();
            y0_Q = stokes[i - 1].Q();
            y0_U = stokes[i - 1].U();
            y0_V = stokes[i - 1].V();
        }
        else
            continue;

        final_stokes.addI(y0_I * dx + 0.5 * dx * dx * m_I);
        final_stokes.addQ(y0_Q * dx + 0.5 * dx * dx * m_Q);
        final_stokes.addU(y0_U * dx + 0.5 * dx * dx * m_U);
        final_stokes.addV(y0_V * dx + 0.5 * dx * dx * m_V);
    }

    return final_stokes;
}

double CMathFunctions::full_integ(dlist x, dlist y, double xlow, double xup)
{
    double res = 0;
    double dx, y0, m;

    if(xup < xlow)
        return res;

    for(uint i = 1; i < x.size(); i++)
    {
        m = (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
        if(x[i - 1] < xlow && x[i] >= xup)
        {
            dx = xup - xlow;
            y0 = y[i - 1] + (xlow - x[i - 1]) * m;
            if(dx == 0)
                return y0;
        }
        else if(x[i] >= xlow && (x[i - 1] < xlow))
        {
            dx = x[i] - xlow;
            y0 = y[i];
            m *= -1;
        }
        else if(x[i - 1] < xup && (x[i] >= xup))
        {
            dx = xup - x[i - 1];
            y0 = y[i - 1];
        }
        else if(x[i - 1] >= xlow && x[i] < xup)
        {
            dx = x[i] - x[i - 1];
            y0 = y[i - 1];
        }
        else
            continue;

        res += y0 * dx + 0.5 * dx * dx * m;
    }

    return res;
}

void CMathFunctions::probListInteg(const double * x, const double * y, prob_list & integ_spline)
{
    uint N = integ_spline.size();
    double res = 0;

    if(N == 1)
        integ_spline.setValue(0, 1);
    else
    {
        integ_spline.setValue(0, 0);
        for(uint i = 1; i < N; i++)
        {
            if(y[i - 1] > 0)
            {
                if(y[i] == 0)
                    res += 0.5 * (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1]) * (x[i] - x[i - 1]) * (x[i] - x[i - 1]);
                else
                    res += (x[i] - x[i - 1]) * y[i - 1] + 0.5 * (x[i] - x[i - 1]) * (y[i] - y[i - 1]);
            }
            integ_spline.setValue(i, res);
        }
    }

    integ_spline.normalize(res);
}

void CMathFunctions::probListInteg(dlist x, const double * y, prob_list & integ_spline)
{
    uint N = integ_spline.size();
    double res = 0;

    if(N == 1)
        integ_spline.setValue(0, 1);
    else
    {
        integ_spline.setValue(0, 0);
        for(uint i = 1; i < N; i++)
        {
            res += (x[i] - x[i - 1]) * y[i - 1] + 0.5 * (x[i] - x[i - 1]) * (y[i] - y[i - 1]);
            integ_spline.setValue(i, res);
        }
    }

    integ_spline.normalize(res);
}

double CMathFunctions::calc_delta(double B, double Td, double Tg, double ng)
{
    double den = ng * Td * sqrt(Tg);

    if(den == 0)
        return 1;

    return B * B / den;
}

double CMathFunctions::calc_mach(double vel, double Tg, double mu)
{
    return vel / sqrt(con_kB * Tg / (mu * m_H));
}

double CMathFunctions::calc_larm_limit(double B, double Td, double Tg, double ng, double s, double larm_f)
{
    double den = ng * Td * sqrt(Tg);

    if(den == 0)
        return 1;

    return s * s * B / den / larm_f;
}

double CMathFunctions::planck(double l, double T)
{
    return (2.0 * con_h * con_c * con_c) /
            ((l * l * l * l * l) * (exp((con_h * con_c) / (l * con_kB * T)) - 1));
}

double CMathFunctions::planck_hz(double f, double T)
{

    return (2.0 * con_h * f * f * f) / ((con_c * con_c) * (exp((con_h * f) / (con_kB * T)) - 1));
}

double CMathFunctions::dplanck_dT(double l, double T)
{
    double ex = exp((con_h * con_c) / (l * con_kB * T));
    ex = (2.0 * con_c * con_c * con_c * ex * con_h * con_h) /
            ((ex - 1) * (ex - 1) * con_kB * l * l * l * l * l * l * T * T);

    if(ex != ex)
        return 0;

    return ex;
}

double CMathFunctions::mathis_isrf(double wavelength)
{
    if(wavelength >= 0.0912e-6 && wavelength < 0.11e-6)
        return 3069.0 * pow((wavelength * 1e6), 3.4172);
    else if(wavelength >= 0.11e-6 && wavelength < 0.134e-6)
        return 1.627;
    else if(wavelength >= 0.134e-6 && wavelength < 0.25e-6)
        return 0.0566 * pow((wavelength * 1e6), -1.6678);
    else if(wavelength >= 0.25e-6)
        return 1e-14 * planck(wavelength, 7500.0) + 1e-13 * planck(wavelength, 4000.0) +
                4e-13 * planck(wavelength, 3000.0);

    return 0;
}

void CMathFunctions::SinList(double start, double stop, double * list, uint N, double f)
{
    if(f == 1)
    {
        double inter = (stop - start);
        double dang = PI / double(N - 1);

        double mid = (double(N) - 1.5) / 2.0;

        list[0] = start;

        for(uint i_x = 1; i_x < N - 1; i_x++)
        {
            if(double(i_x) <= mid)
                list[i_x] = start + inter * (0.5 * sin(i_x * dang));
            else
                list[i_x] = start + inter * (1 - 0.5 * sin(i_x * dang));
        }

        list[N - 1] = stop;
    }
    else
    {
        LinearList(start, stop, list, N);
    }
}

void CMathFunctions::LinearList(double start, double stop, double * list, uint N)
{
    list[0] = start;
    
    if(N > 1){
        double dx = (stop - start) / (N - 1);

        for(uint i_x = 1; i_x < N - 1; i_x++)
            list[i_x] = start + i_x * dx;

        list[N - 1] = stop;
    }
}

void CMathFunctions::LinearList(double start, double stop, dlist & list)
{
    uint N = list.size();
    list[0] = start;

    if(N > 1){
        double dx = (stop - start) / (N - 1);

        for(uint i_x = 1; i_x < N - 1; i_x++)
            list[i_x] = start + i_x * dx;

        list[N - 1] = stop;
    }
}

dlist CMathFunctions::LinearList(double start, double stop, uint N)
{
    dlist list(N);
    list[0] = start;
    
    if(N > 1){
        double dx = (stop - start) / (N - 1);

        for(uint i_x = 1; i_x < N - 1; i_x++)
            list[i_x] = start + i_x * dx;

        list[N - 1] = stop;
    }

    return list;
}

void CMathFunctions::ExpList(double start, double stop, double * list, uint N, double base)
{
    if(N == 1)
        list[0] = start;
    else if(base > 1)
    {
        N--;
        double dx = (stop - start) * (base - 1.0) / (pow(base, double(N)) - 1);

        for(uint i_x = 0; i_x <= N; i_x++)
        {
            list[i_x] = start + dx * (pow(base, double(i_x)) - 1) / (base - 1.0);
        }

        list[0] = start;
        list[N] = stop;
    }
    else
    {
        LinearList(start, stop, list, N);
    }
}

void CMathFunctions::ExpListSym(double start, double stop, double * list, uint N, double base)
{
    if(N == 1)
        list[0] = start;
    else if(base > 1)
    {
        uint tmpN;
        double tmp_mid, tmp_stop;

        tmp_mid = start + 0.5 * (stop - start);
        tmp_stop = stop;

        if(N % 2 == 0)
        {
            tmpN = N / 2 + 1;
            double * tmp_list = new double[tmpN];
            CMathFunctions::ExpList(tmp_mid, tmp_stop, tmp_list, tmpN, base);

            for(uint i_x = uint(N / 2); i_x < N; i_x++)
                list[i_x] = tmp_list[i_x - uint(N / 2) + 1];

            double diff = 0;

            for(uint i_x = 0; i_x < uint(N / 2); i_x++)
            {
                diff += tmp_list[i_x + 1] - tmp_list[i_x];
                list[uint(N / 2) - 1 - i_x] = tmp_mid - diff;
            }
            delete[] tmp_list;
        }
        else
        {
            tmpN = (N + 1) / 2;
            uint midN = uint((N - 1) / 2);
            double * tmp_list = new double[tmpN];
            CMathFunctions::ExpList(tmp_mid, tmp_stop, tmp_list, tmpN, base);

            for(uint i_x = midN; i_x < N; i_x++)
                list[i_x] = tmp_list[i_x - midN];

            double diff = 0;

            for(uint i_x = 0; i_x < midN; i_x++)
            {
                diff += tmp_list[i_x + 1] - tmp_list[i_x];
                list[midN - 1 - i_x] = tmp_mid - diff;
            }
            delete[] tmp_list;
        }
    }
    else
        LinearList(start, stop, list, N);

    list[0] = start;
    list[N - 1] = stop;
}

void CMathFunctions::LogList(double start, double stop, double * list, uint N, double base)
{
    if(N == 1)
        list[0] = start;
    else if(base > 1)
    {
        double log_start = log10(start) / log10(base);
        double log_stop = log10(stop) / log10(base);

        double dx = (log_stop - log_start) / (N - 1);

        list[0] = start;

        for(uint i_x = 1; i_x < N - 1; i_x++)
            list[i_x] = pow(base, log_start + i_x * dx);

        list[N - 1] = stop;
    }
    else
    {
        LinearList(start, stop, list, N);
    }
}

void CMathFunctions::LogList(double start, double stop, dlist & list, double base)
{
    uint N = list.size();

    if(N == 1)
        list[0] = start;
    else if(base > 1)
    {
        double log_start = log10(start) / log10(base);
        double log_stop = log10(stop) / log10(base);

        double dx = (log_stop - log_start) / (N - 1);

        list[0] = start;

        for(uint i_x = 1; i_x < N - 1; i_x++)
            list[i_x] = pow(base, log_start + i_x * dx);

        list[N - 1] = stop;
    }
}

void CMathFunctions::SymLogList(double start, double stop, double * list, uint N, double base)
{
    if(base > 1 && N > 3)
    {
        double inter = stop - start;
        double log_start = log10(start + 0.5 * inter) / log10(base);
        double log_stop = log10(stop) / log10(base);

        double dx = (log_stop - log_start) / (0.5 * double(N) - 1);

        list[0] = start;

        for(uint i_x = uint(0.5 * N); i_x < N - 1; i_x++)
        {
            uint pos1 = i_x;
            uint pos2 = uint(0.5 * N - i_x);
            list[pos1] = pow(base, log_start + i_x * dx);
            list[pos2] = pow(base, log_start + i_x * dx);
        }

        list[N - 1] = stop;
    }
    else
    {
        LinearList(start, stop, list, N);
    }
}

double CMathFunctions::sgn(double x)
{
    if(x < 0)
        return -1.0;
    else
        return 1.0;
}

double CMathFunctions::sgn(int x)
{
    if(x < 0)
        return -1.0;
    else
        return 1.0;
}

double CMathFunctions::sgn(int64_t x)
{
    if(x < 0)
        return -1.0;
    else
        return 1.0;
}

Matrix2D CMathFunctions::getRotationMatrix(double cos_phi,
                                            double sin_phi,
                                            double cos_theta,
                                            double sin_theta)
{
    Matrix2D D;
    D.resize(3, 3);

    D(0, 0) = cos_phi;
    D(1, 0) = sin_phi;
    D(2, 0) = 0.0;
    D(0, 1) = -sin_phi * cos_theta;
    D(1, 1) = cos_phi * cos_theta;
    D(2, 1) = -sin_theta;
    D(0, 2) = -sin_phi * sin_theta;
    D(1, 2) = cos_phi * sin_theta;
    D(2, 2) = cos_theta;

    return D;
}

double CMathFunctions::getRotationAngleObserver(const Vector3D & obs_ex,
                                                const Vector3D & photon_ex,
                                                const Vector3D & photon_ey)
{
    double cos_angle_1 = obs_ex * photon_ey;
    double cos_angle_2 = obs_ex * photon_ex;

    return Vector3D::atan3(cos_angle_1, cos_angle_2);
}

void CMathFunctions::getPropMatrixAPi(double cos_theta,
                                    double sin_theta,
                                    double cos_2_phi,
                                    double sin_2_phi,
                                    double mult,
                                    Matrix2D * propMatrix)
{
    propMatrix->addValue(0, 0, sin_theta * sin_theta * mult);
    propMatrix->addValue(0, 1, -cos_2_phi * sin_theta * sin_theta * mult);
    propMatrix->addValue(0, 2, -sin_2_phi * sin_theta * sin_theta * mult);
    propMatrix->addValue(1, 0, -cos_2_phi * sin_theta * sin_theta * mult);
    propMatrix->addValue(1, 1, sin_theta * sin_theta * mult);
    propMatrix->addValue(2, 0, -sin_2_phi * sin_theta * sin_theta * mult);
    propMatrix->addValue(2, 2, sin_theta * sin_theta * mult);
    propMatrix->addValue(3, 3, sin_theta * sin_theta * mult);
}

void CMathFunctions::getPropMatrixBSigmaP(double cos_theta,
                                        double sin_theta,
                                        double cos_2_phi,
                                        double sin_2_phi,
                                        double mult,
                                        Matrix2D * propMatrix)
{
    propMatrix->addValue(1, 2, -2 * cos_theta * mult);
    propMatrix->addValue(1, 3, sin_2_phi * sin_theta * sin_theta * mult);
    propMatrix->addValue(2, 1, +2 * cos_theta * mult);
    propMatrix->addValue(2, 3, -cos_2_phi * sin_theta * sin_theta * mult);
    propMatrix->addValue(3, 1, -sin_2_phi * sin_theta * sin_theta * mult);
    propMatrix->addValue(3, 2, cos_2_phi * sin_theta * sin_theta * mult);
}

void CMathFunctions::getPropMatrixASigmaP(double cos_theta,
                                        double sin_theta,
                                        double cos_2_phi,
                                        double sin_2_phi,
                                        double mult,
                                        Matrix2D * propMatrix)
{
    propMatrix->addValue(0, 0, (1 + cos_theta * cos_theta) * mult);
    propMatrix->addValue(0, 1, cos_2_phi * sin_theta * sin_theta * mult);
    propMatrix->addValue(0, 2, sin_2_phi * sin_theta * sin_theta * mult);
    propMatrix->addValue(0, 3, -2 * cos_theta * mult);
    propMatrix->addValue(1, 0, cos_2_phi * sin_theta * sin_theta * mult);
    propMatrix->addValue(1, 1, (1 + cos_theta * cos_theta) * mult);
    propMatrix->addValue(2, 0, sin_2_phi * sin_theta * sin_theta * mult);
    propMatrix->addValue(2, 2, (1 + cos_theta * cos_theta) * mult);
    propMatrix->addValue(3, 0, -2 * cos_theta * mult);
    propMatrix->addValue(3, 3, (1 + cos_theta * cos_theta) * mult);
}

void CMathFunctions::getPropMatrixASigmaM(double cos_theta,
                                        double sin_theta,
                                        double cos_2_phi,
                                        double sin_2_phi,
                                        double mult,
                                        Matrix2D * propMatrix)
{
    propMatrix->addValue(0, 0, (1 + cos_theta * cos_theta) * mult);
    propMatrix->addValue(0, 1, cos_2_phi * sin_theta * sin_theta * mult);
    propMatrix->addValue(0, 2, sin_2_phi * sin_theta * sin_theta * mult);
    propMatrix->addValue(0, 3, +2 * cos_theta * mult);
    propMatrix->addValue(1, 0, cos_2_phi * sin_theta * sin_theta * mult);
    propMatrix->addValue(1, 1, (1 + cos_theta * cos_theta) * mult);
    propMatrix->addValue(2, 0, sin_2_phi * sin_theta * sin_theta * mult);
    propMatrix->addValue(2, 2, (1 + cos_theta * cos_theta) * mult);
    propMatrix->addValue(3, 0, +2 * cos_theta * mult);
    propMatrix->addValue(3, 3, (1 + cos_theta * cos_theta) * mult);
}

void CMathFunctions::getPropMatrixBSigmaM(double cos_theta,
                                        double sin_theta,
                                        double cos_2_phi,
                                        double sin_2_phi,
                                        double mult,
                                        Matrix2D * propMatrix)
{
    propMatrix->addValue(1, 2, +2 * cos_theta * mult);
    propMatrix->addValue(1, 3, sin_2_phi * sin_theta * sin_theta * mult);
    propMatrix->addValue(2, 1, -2 * cos_theta * mult);
    propMatrix->addValue(2, 3, -cos_2_phi * sin_theta * sin_theta * mult);
    propMatrix->addValue(3, 1, -sin_2_phi * sin_theta * sin_theta * mult);
    propMatrix->addValue(3, 2, cos_2_phi * sin_theta * sin_theta * mult);
}

void CMathFunctions::getPropMatrixBPi(double cos_theta,
                                    double sin_theta,
                                    double cos_2_phi,
                                    double sin_2_phi,
                                    double mult,
                                    Matrix2D * propMatrix)
{
    propMatrix->addValue(1, 3, -sin_2_phi * sin_theta * sin_theta * mult);
    propMatrix->addValue(2, 3, cos_2_phi * sin_theta * sin_theta * mult);
    propMatrix->addValue(3, 1, sin_2_phi * sin_theta * sin_theta * mult);
    propMatrix->addValue(3, 2, -cos_2_phi * sin_theta * sin_theta * mult);
}

Vector3D CMathFunctions::interpVector(Vector3D y_1, Vector3D y_2, double x_1, double x_2)
{
    return (y_2 - y_1) / (x_2 - x_1);
}

double CMathFunctions::lum2Jy(double intensity, double wavelength, double distance)
{
    if(abs(intensity) < 1e-47)
        return 0;

    double freq = con_c / wavelength;             // [s^-1]
    double L = intensity * con_c / (freq * freq); // [W Hz^-1 sr^-1]
    return L * 1e+26 / (distance * distance);     // [Jy]
}

void CMathFunctions::lum2Jy(StokesVector * S, double wavelength, double distance)
{
    S->setI(lum2Jy(S->I(), wavelength, distance));
    S->setQ(lum2Jy(S->Q(), wavelength, distance));
    S->setU(lum2Jy(S->U(), wavelength, distance));
    S->setV(lum2Jy(S->V(), wavelength, distance));
}

uint CMathFunctions::findListIndex(double * list, uint l_limit, uint h_limit, double val)
{
    uint i;
    h_limit--;
    if(val <= list[l_limit])
        return l_limit;
    if(val >= list[h_limit])
        return h_limit;

    while(h_limit - l_limit > 1)
    {
        i = l_limit + (h_limit - l_limit) / 2;
        if(list[i] > val)
            h_limit = i;
        else
            l_limit = i;
    }
    return l_limit;
}

void CMathFunctions::initChiSq(uint _N, uint _M)
{
    N = _N;
    M = _M;
    if(b != 0)
        delete b;
    if(a != 0)
        delete a;
    if(d != 0)
        delete d;
    if(r != 0)
        delete r;
    if(y != 0)
        delete y;

    if(sigma != 0)
        delete sigma;

    a = new double[M];
    b = new double[M];
    d = new double[N];
    r = new double[N];
    y = new double[N];
    A.resize(M, M);
    C.resize(N, M);
}

double CMathFunctions::getPolarInterp(double t, double s, double Q1, double Q2, double Q3, double Q4)
{
    return Q1 * (1 - s) * (1 - t) + Q2 * s * (1 - t) + Q3 * (1 - s) * t + Q4 * s * t;
}

double CMathFunctions::phaseFunctionHG(double g, double theta)
{
    double g_sq = g * g;
    double res = 1.0 / PIx4;
    res *= (1.0 - g_sq);
    res /= pow(1.0 + g_sq - 2.0 * g * cos(theta), 1.5);

    return res;
}

double CMathFunctions::phaseFunctionDHG(double g, double alpha, double theta)
{
    double cos_th = cos(theta);
    double res = phaseFunctionHG(g, theta);
    res *= (1.0 + alpha * cos_th * cos_th);
    res /= (1.0 + alpha * (1.0 + 2.0 * g * g) / 3.0);

    return res;
}

double CMathFunctions::phaseFunctionTTHG(double g1, double g2, double weight, double theta)
{
    double res = weight * phaseFunctionHG(g1, theta);
    res += (1.0 - weight) * phaseFunctionHG(g2, theta);

    return res;
}

double CMathFunctions::getPhiIntegral(double phi, dlist args)
{
    double phipar = args[0];
    double rnd = args[1];

    double res = (phi - phipar * 0.5 * sin(2.0 * phi)) / PIx2;

    res -= rnd;
    if(abs(res) < 2.0 * EPS_DOUBLE) {
        res = 0.0;
    }
    return res;
}

double CMathFunctions::getDHGIntegral(double cos_theta, dlist args)
{
    double g = args[0];
    double alpha = args[1];
    double rnd = args[2];

    double g_sq = g * g;
    double res;

    if(abs(g) < 2.0 * EPS_DOUBLE) { // if g is close to zero, use integral of Rayleigh phase function
        res = alpha * (1.0 + cos_theta * cos_theta * cos_theta) / 3.0 + (1.0 + cos_theta);
        res /= (2.0 * (1.0 + alpha / 3.0));
    } else {
        res = 2.0 * alpha * (g_sq * g_sq - g_sq * g * cos_theta - 0.5 * g_sq * (cos_theta * cos_theta - 4.0) - g * cos_theta + 1.0) + (3.0 * g_sq);
        res /= (3.0 * g_sq * g * sqrt(1.0 + g_sq - 2.0 * g * cos_theta));
        res -= (2.0 * alpha * (g_sq * g_sq + g_sq * g + 1.5 * g_sq + g + 1.0) + 3.0 * g_sq) / (3.0 * g_sq * g * abs(1.0 + g));
        res *= 0.5 * (1.0 - g_sq) / (1.0 + alpha * (1.0 + 2.0 * g_sq) / 3.0);
    }

    res -= rnd;
    if(abs(res) < 2.0 * EPS_DOUBLE) {
        res = 0.0;
    }
    return res;
}

double CMathFunctions::getTTHGIntegral(double cos_theta, dlist args)
{
    double g1 = args[0];
    double g2 = args[1];
    double weight = args[2];
    double rnd = args[3];

    double g1_sq = g1 * g1;
    double g2_sq = g2 * g2;

    double res = 0.0;
    if(abs(g1) < 2.0 * EPS_DOUBLE) { // if g1 is close to zero, use integral of 1/2 (random isotropic direction)
        res += 0.5 * weight * (1.0 + cos_theta);
    } else {
        res += 0.5 * weight * ((1.0 - g1_sq) / sqrt(1.0 + g1_sq - 2.0 * g1 * cos_theta) - (1.0 - g1)) / g1;
    }
    if(abs(g2) < 2.0 * EPS_DOUBLE) { // if g2 is close to zero, use integral of 1/2 (random isotropic direction)
        res += 0.5 * (1.0 - weight) * (1.0 + cos_theta);
    } else {
        res += 0.5 * (1.0 - weight) * ((1.0 - g2_sq) / sqrt(1.0 + g2_sq - 2.0 * g2 * cos_theta) - (1.0 - g2)) / g2;
    }

    res -= rnd;
    if(abs(res) < 2.0 * EPS_DOUBLE) {
        res = 0.0;
    }
    return res;
}

double CMathFunctions::findRootBrent(double a, double b, double (*func)(double, dlist), dlist args)
{
    double t = 10.0 * EPS_DOUBLE;
    double sa, sb, c, d, e, fa, fb, fc, tol, m, p, q, r, s;
    uint maxiter = 100;
    uint iter = 0;

    sa = a;
    sb = b;
    fa = func(sa, args);
    fb = func(sb, args);

    if(fa * fb > 0.0) {
        cout << ERROR_LINE << "No root found. f(a) and f(b) must have different signs." << endl;
        cout << "f(a) = " << fa << ", f(b) = " << fb << endl;
        for(uint i = 0; i < args.size(); i++) {
            cout << "arg(" << i << ") = " << args[i] << endl;
        }
        return 0;
    }

    c = sa;
    fc = fa;
    e = sb - sa;
    d = e;

    while(iter < maxiter) {
        if(abs(fc) < abs(fb)) {
            sa = sb;
            sb = c;
            c = sa;
            fa = fb;
            fb = fc;
            fc = fa;
        }

        tol = 2.0 * EPS_DOUBLE * abs(sb) + t;
        m = 0.5 * (c - sb);

        if((abs(m) <= tol) || (fb == 0.0)) {
            break;
        }
        
        if((abs(e) >= tol) && (abs(fa) >= abs(fb))) { // see if a bisection is forced
            s = fb / fa;
            if(sa != c) { // inverse quadratic interpolation
                q = fa / fc;
                r = fb / fc;
                p = s * (2.0 * m * q * (q - r) - (sb - sa) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            } else { // linear interpolation
                p = 2.0 * m * s;
                q = 1.0 - s;
            }

            if(p <= 0.0) {
                p = -p;
            } else {
                q = -q;
            }

            s = e;
            e = d;

            if((2.0 * p >= 3.0 * m * q - abs(tol * q)) || (p >= abs(0.5 * s * q))) {
                e = m;
                d = e;
            } else {
                d = p / q;
            }
        } else {
            e = m;
            d = e;
        }

        sa = sb;
        fa = fb;

        if(abs(d) <= tol) {
            if(m <= 0.0) {
                sb = sb - tol;
            } else {
                sb = sb +  tol;
            }
        } else {
            sb = sb + d;
        }

        fb = func(sb, args);

        if(((fb > 0.0) && (fc > 0.0)) || ((fb <= 0.0) && (fc <= 0.0))) {
            c = sa;
            fc = fa;
            e = sb - sa;
            d = e;
        }

        iter++;
    }

    if(iter == maxiter)
    {
        cout << ERROR_LINE << "No root found. Maximum number of iterations reached." << endl;
        for(uint i = 0; i < args.size(); i++) {
            cout << "arg(" << i << ") = " << args[i] << endl;
        }
        return 0;
    }

    // sb   : root of the function
    // fb   : value of the function (approximately zero)
    // iter : number of iterations needed to find the root
    return sb;
}

double CMathFunctions::Freq2Velo(double _f, double _f0)
{
    return con_c * _f / _f0;
}

double CMathFunctions::Velo2Freq(double _v, double _f0)
{
    return _v / con_c * _f0;
}

double CMathFunctions::getWavelengthStepWidth(double lam_min, double lam_max, double nr_of_spectral_bins)
{
    double dw;
    if(nr_of_spectral_bins == 1)
        dw = 0;
    else
        dw = (lam_max - lam_min) / double(nr_of_spectral_bins - 1);
    return dw;
}

void CMathFunctions::gauss(Matrix2D & inA, double * b, double * res, uint n)
{
    // int n = inb.size();
    double x, sum;
    Matrix2D A = inA;
    // dlist res;

    /*b[0] = inb.X();
    b[1] = inb.Y();
    b[2] = inb.Z();*/

    for(int k = 0; k < int(n) - 1; k++)
    {

        for(int i = k + 1; i < int(n); i++)
        {
            x = A(i, k) / A(k, k);

            for(int j = k + 1; j < int(n); j++)
            {
                A(i, j) = A(i, j) - A(k, j) * x;
            }

            b[i] = b[i] - b[k] * x;
        }
    }

    // Resubstitution

    // res.resize(n);
    b[n - 1] = b[n - 1] / A(n - 1, n - 1);
    for(int i = int(n) - 2; i >= 0; i--)
    {
        sum = b[i];

        for(int j = i + 1; j < int(n); j++)
        {
            sum = sum - A(i, j) * b[j];
        }

        b[i] = sum / A(i, i);
        res[i] = b[i];
    }

    res[n - 1] = b[n - 1];
    // return res;
}

dlist CMathFunctions::gauss_solver(Matrix2D & A, uint n, dlist & b)
{
    /*Basic Gaussian elimination with pivoting*/

    /* DESCRIPTION:
        - Algorithm for solving system of n linear equations
        with n unknowns (x1,x2,...,xn)
        - Gaussian elimination is algorithm for solving
        system of linear equations. This method can be
        used to calculate determinant, inverse matrix or
        find inverse matrix.
        Author: Ervin B,
        May 2013
        */
    dlist x(n);
    for(uint i = 0; i < n; i++)
    {
        x[i] = 0;
    }

    for(uint i = 0; i < n; i++)
    {
        // Search for maximum in this column
        double maxEl = abs(A(i, i));
        int maxRow = i;
        for(uint k = i + 1; k < n; k++)
        {
            if(abs(A(k, i)) > maxEl)
            {
                maxEl = abs(A(k, i));
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        for(uint k = i; k < n; k++)
        {
            double tmp = A(maxRow, k);
            A(maxRow, k) = A(i, k);
            A(i, k) = tmp;
        }
        double tmp = b[maxRow];
        b[maxRow] = b[i];
        b[i] = tmp;

        // Make all rows below this one 0 in current column
        for(uint k = i + 1; k < n; k++)
        {
            double c = -A(k, i) / A(i, i);
            for(uint j = i; j < n; j++)
            {
                if(i == j)
                {
                    A(k, j) = 0;
                }
                else
                {
                    A(k, j) += c * A(i, j);
                }
                b[k] += c * b[i];
            }
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    for(int i = n - 1; i >= 0; i--)
    {
        x[i] = b[i] / A(i, i);
        for(int k = i - 1; k >= 0; k--)
        {

            b[k] -= A(k, i) * x[i];
        }
    }
    return x;
}

double CMathFunctions::fmodulo(double v1, double v2)
{
    if(v1 >= 0)
        return (v1 < v2) ? v1 : fmod(v1, v2);
    double tmp = fmod(v1, v2) + v2;
    return (tmp == v2) ? 0. : tmp;
}

int CMathFunctions::imodulo(int v1, int v2)
{
    int v = v1 % v2;
    return (v >= 0) ? v : v + v2;
}

dcomplex CMathFunctions::Csqrt(dcomplex z)
{
    dcomplex c;
    float x, y, w, r;

    if(real(z) == 0.0 && imag(z) == 0.0)
        c = dcomplex(0.0, 0.0);
    else
    {
        x = fabs(real(z));
        y = fabs(imag(z));
        if(x >= y)
        {
            r = y / x;
            w = sqrt(x) * sqrt(0.5 * (1.0 + sqrt(1.0 + r * r)));
        }
        else
        {
            r = x / y;
            w = sqrt(y) * sqrt(0.5 * (r + sqrt(1.0 + r * r)));
        }
        if(real(z) >= 0.0)
        {
            c.real(w);
            c.imag(imag(z) / (2.0 * w));
        }
        else
        {
            c.imag((imag(z) >= 0) ? w : -w);
            c.real(imag(z) / (2.0 * imag(c)));
        }
    }
    return c;
}

bool calcBHMie(double x,
                        dcomplex refractive_index,
                        double & qext,
                        double & qabs,
                        double & qsca,
                        double & gsca,
                        double * S11,
                        double * S12,
                        double * S33,
                        double * S34)
{
    dcomplex cxy = dcomplex(x, 0) * refractive_index;

    // Series expansion terminated after XSTOP terms
    float xstop = x + 4.0 * pow(x, 1.0 / 3.0) + 2.0;
    long nmx = fmax(xstop, abs(cxy)) + 15;

    if(nmx >= MAX_MIE_ITERATIONS)
    {
        cout << ERROR_LINE << "Failure in Mie-scattering calculation (NMX = " << nmx
                << " >= MAX_MIE_ITERATIONS = " << MAX_MIE_ITERATIONS << ")" << endl;
        return false;
    }

    float amu[NANG];
    float dang = 0.5 * PI / float(NANG - 1);
    for(int j = 0; j < NANG; j++)
        amu[j] = cos(float(j) * dang);

    // Logarithmic derivative D(J) calculated by downward recurrence
    // beginning with initial value (0., 0.) at J=NMX
    dcomplex cxd[MAX_MIE_ITERATIONS];
    cxd[nmx] = dcomplex(0, 0);

    dcomplex cxtemp;
    for(long n = 0; n < nmx - 1; n++)
    {
        float rn = float(nmx - n);
        cxd[nmx - (n + 1)] =
            dcomplex(rn, 0) / cxy - dcomplex(1.0, 0.0) / (cxd[nmx - n] + dcomplex(rn, 0) / cxy);
    }

    float pi[NANG], pi0[NANG], pi1[NANG];
    for(int j = 0; j < NANG; j++)
    {
        pi0[j] = 0.0;
        pi1[j] = 1.0;
    }

    dcomplex cxs1[2 * NANG - 1], cxs2[2 * NANG - 1];
    for(int j = 0; j < 2 * NANG - 1; j++)
    {
        cxs1[j] = dcomplex(0, 0);
        cxs2[j] = dcomplex(0, 0);
    }

    // Riccati-Bessel functions with real argument X calculated by upward recurrence
    double dn, psi;
    double psi0 = cos(x);
    double psi1 = sin(x);
    float rn, fn, apsi, chi;
    float chi0 = -sin(x);
    float chi1 = cos(x);
    float apsi0 = psi0;
    float apsi1 = psi1;
    float tau[NANG];
    qsca = 0.0;
    gsca = 0.0;
    dcomplex cxxi;
    dcomplex cxxi0 = dcomplex(apsi0, -chi0);
    dcomplex cxxi1 = dcomplex(apsi1, -chi1);
    dcomplex cxan, cxan1, cxbn, cxbn1;

    for(long n = 1; n <= long(xstop); n++)
    {
        dn = n;
        rn = n;
        fn = (2.0 * rn + 1.0) / (rn * (rn + 1.0));
        psi = (2.0 * dn - 1.0) * psi1 / x - psi0;
        apsi = psi;
        chi = (2.0 * rn - 1.0) * chi1 / x - chi0;
        cxxi = dcomplex(apsi, -chi);

        // Store previous values of AN and BN for use in computation of g = <cos(theta)>
        if(n > 1)
        {
            cxan1 = cxan;
            cxbn1 = cxbn;
        }

        // Compute AN and BN
        cxan = (cxd[n] / refractive_index + dcomplex(rn / x, 0)) * dcomplex(apsi, 0) - dcomplex(apsi1, 0);
        cxan = cxan / ((cxd[n] / refractive_index + dcomplex(rn / x, 0)) * cxxi - cxxi1);
        cxbn = (refractive_index * cxd[n] + dcomplex(rn / x, 0)) * dcomplex(apsi, 0) - dcomplex(apsi1, 0);
        cxbn = cxbn / ((refractive_index * cxd[n] + dcomplex(rn / x, 0)) * cxxi - cxxi1);

        // Augment sums for *qsca and g=<cos(theta)>
        qsca = qsca + (2.0 * rn + 1.0) * (abs(cxan) * abs(cxan) + abs(cxbn) * abs(cxbn));
        gsca = gsca + ((2.0 * rn + 1.0) / (rn * (rn + 1.0))) *
                            (real(cxan) * real(cxbn) + imag(cxan) * imag(cxbn));

        if(n > 1)
            gsca = gsca +
                    ((rn - 1.) * (rn + 1.0) / rn) * (real(cxan1) * real(cxan) + imag(cxan1) * imag(cxan) +
                                                    real(cxbn1) * real(cxbn) + imag(cxbn1) * imag(cxbn));

        for(int j = 0; j < NANG; j++)
        {
            int jj = 2 * NANG - 1 - j;
            pi[j] = pi1[j];
            tau[j] = rn * amu[j] * pi[j] - (rn + 1.0) * pi0[j];

            float p = pow(-1.0, n - 1);
            cxs1[j] =
                cxs1[j] + (dcomplex(fn, 0) * (cxan * dcomplex(pi[j], 0) + cxbn * dcomplex(tau[j], 0)));

            float t = pow(-1.0, n);
            cxs2[j] =
                cxs2[j] + dcomplex(fn, 0) * (cxan * dcomplex(tau[j], 0) + cxbn * dcomplex(pi[j], 0));

            if(j != jj)
            {
                cxs1[jj] = cxs1[jj] + dcomplex(fn, 0) * (cxan * dcomplex(pi[j] * p, 0) +
                                                            cxbn * dcomplex(tau[j] * t, 0));
                cxs2[jj] = cxs2[jj] + dcomplex(fn, 0) * (cxan * dcomplex(tau[j] * t, 0) +
                                                            cxbn * dcomplex(pi[j] * p, 0));
            }
        }

        psi0 = psi1;
        psi1 = psi;
        apsi1 = psi1;
        chi0 = chi1;
        chi1 = chi;
        cxxi1 = dcomplex(apsi1, -chi1);

        // For each angle J, compute pi_n+1 from PI = pi_n , PI0 = pi_n-1
        for(int j = 0; j < NANG; j++)
        {
            pi1[j] = ((2.0 * rn + 1.0) * amu[j] * pi[j] - (rn + 1.0) * pi0[j]) / rn;
            pi0[j] = pi[j];
        }
    }

    // Have summed sufficient terms. Now compute *qsca,*qext,*qback,and *gsca
    gsca = 2.0 * gsca / qsca;
    qsca = (2.0 / (x * x)) * qsca;
    qext = (4.0 / (x * x)) * real(cxs1[0]);
    qabs = qext - qsca;
    // qback = (4.0/(x * x)) * abs(cxs1[2 * NANG - 1]) * abs(cxs1[2 * NANG - 1]);

    for(int j = 0; j < 2 * NANG - 1; j++)
    {
        S11[j] = 0.5 * (abs(cxs2[j]) * abs(cxs2[j]) + abs(cxs1[j]) * abs(cxs1[j]));
        S12[j] = 0.5 * (abs(cxs2[j]) * abs(cxs2[j]) - abs(cxs1[j]) * abs(cxs1[j]));
        S33[j] = real(cxs2[j] * conj(cxs1[j]));
        S34[j] = imag(cxs2[j] * conj(cxs1[j]));
    }

    return true;
}

bool CMathFunctions::calcWVMie(double x,
                        dlist scat_angle,
                        dcomplex refractive_index,
                        double & qext,
                        double & qabs,
                        double & qsca,
                        double & gsca,
                        double * S11,
                        double * S12,
                        double * S33,
                        double * S34)
// Wolf & Voshchinnikov approximation of optical properties for spherical grains.
// see Wolf & Voshchinnikov (2004), Comput. Phys. Commun. 162, 113
{
    // Step width
    uint n_scat_angle = scat_angle.size();
    double factor = 1e250;

    if(x <= MIN_MIE_SIZE_PARAM) {
        cout << ERROR_LINE << "Mie scattering limit exceeded, current size parameter: " << x << endl;
        return false;
    }

    double ax = 1.0 / x;
    double b = 2.0 * ax * ax;
    dcomplex ss(0.0, 0.0);
    dcomplex s3(0.0, -1.0);
    double an = 3.0;

    // choice of number for subroutine aa [Loskutov (1971)]
    double y = abs(refractive_index) * x;
    uint num = uint(1.25 * y + 15.5);
    if(y < 1.0) {
        num = uint(7.5 * y + 9.0);
    } else if(y > 100.0 && y < 50000.0) {
        num = uint(1.0625 * y + 28.5);
    } else if(y >= 50000.0) {
        num = uint(1.005 * y + 50.5);
    }

    if(num > MAX_MIE_ITERATIONS) {
        cout << ERROR_LINE << "Maximum number of terms : " << MAX_MIE_ITERATIONS << ", number of terms required: " << num << endl;
        cout << "  increase default value of MAX_MIE_ITERATIONS in src/Typedefs.h" << endl;
        return false;
        // return calcGeometricOptics(x, refractive_index, qext, qabs,
        //    qabs, gsca, S11, S12, S33, S34);
    }

    // logarithmic derivative to Bessel function (complex argument)
    dcomplex *ru = new dcomplex[num];
    dcomplex s_tmp = ax / refractive_index;
    ru[num-1] = dcomplex(num + 1, 0.0) * s_tmp;
    for(uint i = num - 1; i >= 1; i--) {
        dcomplex s1 = double(i + 1) * s_tmp;
        ru[i-1] = s1 - dcomplex(1.0, 0.0) / (ru[i] + s1);
    }

    // Bessel functions
    double ass = 1.0 / sqrt(PI2 * ax);
    double w1 = invPI2 * ax;
    double Si = sin(x) / x;
    double Co = cos(x) / x;

    // n=0
    double besJ0 = Si * ass;
    double besY0 = -Co * ass;
    uint iu0 = 0;

    // n=1
    double besJ1 = (Si * ax - Co) * ass;
    double besY1 = (-Co * ax - Si) * ass;
    uint iu1 = 0;
    uint iu2 = 0;

    // Mie coefficients (first term)
    dcomplex s, s1, s2, ra0, rb0;

    // coefficient a_1
    s = ru[0] / refractive_index + ax;
    s1 = s * besJ1 - besJ0;
    s2 = s * besY1 - besY0;
    ra0 = s1 / (s1 - s3 * s2);

    // coefficient b_1
    s = ru[0] * refractive_index + ax;
    s1 = s * besJ1 - besJ0;
    s2 = s * besY1 - besY0;
    rb0 = s1 / (s1 - s3 * s2);

    // efficiency factors (first term)
    dcomplex r = -1.5 * (ra0 - rb0);
    qext = an * real(ra0 + rb0);
    qsca = an * (norm(ra0) + norm(rb0));

    // scattering amplitude functions
    double FN = 1.5;
    dlist dPI(n_scat_angle), dTAU(n_scat_angle);
    dlist dAMU(n_scat_angle), dPI0(n_scat_angle, 0), dPI1(n_scat_angle, 1);
    dcomplex SM1[n_scat_angle], SM2[n_scat_angle];
    for(uint i_scat_ang = 0; i_scat_ang < n_scat_angle; i_scat_ang++) {
        dAMU[i_scat_ang] = cos(scat_angle[i_scat_ang]);

        SM1[i_scat_ang] = dcomplex(0.0, 0.0);
        SM2[i_scat_ang] = dcomplex(0.0, 0.0);

        dTAU[i_scat_ang] = dAMU[i_scat_ang] * dPI1[i_scat_ang] - 2.0 * dPI0[i_scat_ang];

        SM1[i_scat_ang] = SM1[i_scat_ang] + FN * (ra0 * dPI1[i_scat_ang] + rb0 * dTAU[i_scat_ang]);
        SM2[i_scat_ang] = SM2[i_scat_ang] + FN * (ra0 * dTAU[i_scat_ang] + rb0 * dPI1[i_scat_ang]);

        dPI[i_scat_ang] = dPI1[i_scat_ang];
        dPI1[i_scat_ang] *= (dAMU[i_scat_ang] * 3.0);
        dPI1[i_scat_ang] -= (dPI0[i_scat_ang] * 2.0);
        dPI0[i_scat_ang] = dPI[i_scat_ang];
    }

    // 2., 3., ... num
    double z = -1.0, besY2, besJ2, an2, qq, r_iterm;
    dcomplex ra1, rb1, rr;
    for(uint iterm = 2; iterm <= num; iterm++) {
        an = an + 2.0;
        an2 = an - 2.0;

        // Bessel functions
        if(iu1 == iu0) {
            besY2 = an2 * ax * besY1 - besY0;
        } else {
            besY2 = an2 * ax * besY1 - besY0 / factor;
        }

        if(abs(besY2) > 1e200) {
            besY2 = besY2 / factor;
            iu2 = iu1 + 1;
        }

        // rbrunngraeber 10/14: Changed from besJ2 = (w1 + besY2 * besJ1) / besY1,
        // because besY2*besJ1 could become very large (1e300) for large grain sizes,
        // besY2/besY1 is about 1; suggested by fkirchschlager
        besJ2 = besY2 / besY1;
        besJ2 = w1 / besY1 + besJ2 * besJ1;

        // Mie coefficients
        r_iterm = double(iterm);

        s = ru[iterm-1] / refractive_index + r_iterm * ax;
        s1 = s * (besJ2 / factorial(iu2)) - besJ1 / factorial(iu1);
        s2 = s * (besY2 * factorial(iu2)) - besY1 * factorial(iu1);
        ra1 = s1 / (s1 - s3 * s2); // coefficient a_n, (n=iterm)

        s = ru[iterm-1] * refractive_index + r_iterm * ax;
        s1 = s * (besJ2 / factorial(iu2)) - besJ1 / factorial(iu1);
        s2 = s * (besY2 * factorial(iu2)) - besY1 * factorial(iu1);
        rb1 = s1 / (s1 - s3 * s2); // coefficient b_n, (n=iterm)

        // efficiency factors
        z = -z;
        rr = z * (r_iterm + 0.5) * (ra1 - rb1);
        r = r + rr;
        ss = ss + (r_iterm - 1.0) * (r_iterm + 1.0) / r_iterm * (ra0 * conj(ra1) + rb0 * conj(rb1)) +
                an2 / r_iterm / (r_iterm - 1.0) * (ra0 * conj(rb0));
        qq = an * real(ra1 + rb1);
        qext = qext + qq;
        qsca = qsca + an * (norm(ra1) + norm(rb1));

        // leaving-the-loop with error criterion
        if(isnan(qext)) {
            cout << ERROR_LINE << "Qext is not a number" << endl;
            return false;
        }

        // leaving-the-loop criterion
        if(abs(qq / qext) < MIE_ACCURACY) {
            break;
        }

        // Bessel functions
        besJ0 = besJ1;
        besJ1 = besJ2;
        besY0 = besY1;
        besY1 = besY2;
        iu0 = iu1;
        iu1 = iu2;
        ra0 = ra1;
        rb0 = rb1;

        // scattering amplitude functions
        FN = (2.0 * r_iterm + 1.0) / (r_iterm * (r_iterm + 1.0));
        for(uint i_scat_ang = 0; i_scat_ang < n_scat_angle; i_scat_ang++) {
            dTAU[i_scat_ang] = r_iterm * dAMU[i_scat_ang] * dPI1[i_scat_ang] - (r_iterm + 1.0) * dPI0[i_scat_ang];

            SM1[i_scat_ang] = SM1[i_scat_ang] + FN * (ra0 * dPI1[i_scat_ang] + rb0 * dTAU[i_scat_ang]);
            SM2[i_scat_ang] = SM2[i_scat_ang] + FN * (ra0 * dTAU[i_scat_ang] + rb0 * dPI1[i_scat_ang]);

            dPI[i_scat_ang] = dPI1[i_scat_ang];
            dPI1[i_scat_ang] *= (dAMU[i_scat_ang] * (2.0 + 1.0 / r_iterm));
            dPI1[i_scat_ang] -= (dPI0[i_scat_ang] * (1.0 + 1.0 / r_iterm));
            dPI0[i_scat_ang] = dPI[i_scat_ang];
        }
    }

    delete[] ru;

    // efficiency factors (final calculations)
    qext = b * qext;
    qsca = b * qsca;
    // double qbk = 2.0 * b * r * conj(r);
    double qpr = qext - 2.0 * b * real(ss);
    qabs = qext - qsca;
    gsca = (qext - qpr) / qsca;

    for(uint i_scat_ang = 0; i_scat_ang < n_scat_angle; i_scat_ang++) {
        S11[i_scat_ang] = 0.5 * (abs(SM2[i_scat_ang]) * abs(SM2[i_scat_ang]) + abs(SM1[i_scat_ang]) * abs(SM1[i_scat_ang]));
        S12[i_scat_ang] = 0.5 * (abs(SM2[i_scat_ang]) * abs(SM2[i_scat_ang]) - abs(SM1[i_scat_ang]) * abs(SM1[i_scat_ang]));
        S33[i_scat_ang] = real(SM2[i_scat_ang] * conj(SM1[i_scat_ang]));
        S34[i_scat_ang] = imag(SM2[i_scat_ang] * conj(SM1[i_scat_ang]));

        // if SM2 and SM1 get really large (if x >> 1 for instance)
        // then double precision may not be enough
        // to get S12/S34 = 0 for theta = 0/pi (cos(theta) = +-1)
        if(abs(dAMU[i_scat_ang]) == 1.0) {
            S12[i_scat_ang] = 0.0;
            S34[i_scat_ang] = 0.0;
        }
    }

    return true;
}

bool CMathFunctions::calcGeometricOptics(double x,
                                        dcomplex refractive_index,
                                        double & qext,
                                        double & qabs,
                                        double & qsca,
                                        double & gsca,
                                        double * S11,
                                        double * S12,
                                        double * S33,
                                        double * S34)
{
    // Efficiency for extinction is 2 in the limit of x>>1
    qext = 2.0;
    // Scattering Henyey-Greenstein g for Draine and Lee silicates
    gsca = 9.23e-1;

    // Set variables
    double res = 0;
    uint nr_angles = 5000;
    double d_ang = PI2 / double(nr_angles - 1);

    // Calculate from 0 to PI/2
    for(uint i = 0; i < nr_angles; i++)
    {
        if(i == 0 || i == nr_angles - 1)
            res = res + calcReflectionCoefficients(refractive_index, d_ang * double(i)) * d_ang;
        else
            res = res + 0.5 * calcReflectionCoefficients(refractive_index, d_ang * double(i)) * d_ang;
    }

    // Set scattering efficiency
    qsca = 1 + 2 * res;

    // Set absorption efficiency
    qabs = qext - qsca;

    for(int j = 0; j < 2 * NANG - 1; j++)
    {
        if(j == 0)
        {
            S11[j] = 1.0;
            S33[j] = 1.0;
        }
        else
        {
            S11[j] = 0.0;
            S33[j] = 0.0;
        }
        S12[j] = 0;
        S34[j] = 0;
    }

    return true;
}

double CMathFunctions::calcReflectionCoefficients(dcomplex refractive_index, double theta)
{
    dcomplex sin_theta = sin(theta) / refractive_index;
    dcomplex cos_theta = sqrt(dcomplex(1, 0) - (sin_theta * sin_theta));
    // r for E parallel to plane
    dcomplex rpll =
        (cos_theta - refractive_index * cos(theta)) / (cos_theta + refractive_index * cos(theta));

    // r for E perp. to plane
    dcomplex rper =
        (cos(theta) - refractive_index * cos_theta) / (cos(theta) + refractive_index * cos_theta);

    //  R = ½(|rpll|²+|rper|²)
    double res = (norm(rpll) + norm(rper)) / 2.0;
    return res * sin(theta) * cos(theta);
}

int CMathFunctions::factorial(int n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

void CMathFunctions::CholDec()
{
    for(uint k = 1; k <= M; k++)
    {
        A(k - 1, k - 1) = sqrt(A(k - 1, k - 1));

        for(uint i = k + 1; i <= M; i++)
        {
            A(i - 1, k - 1) = A(i - 1, k - 1) / A(k - 1, k - 1);

            for(uint j = k + 1; j <= M; j++)
                A(i - 1, j - 1) = A(i - 1, j - 1) - A(i - 1, k - 1) * A(j - 1, k - 1);
        }
    }
}

void CMathFunctions::SolveEqSys()
{
    double sum;
    for(uint i = 1; i <= M; i++)
    {
        sum = b[i - 1];
        for(uint k = i - 1; k >= 1; k--)
            sum -= A(i - 1, k - 1) * a[k - 1];

        a[i - 1] = sum / A(i - 1, i - 1);
    }

    for(uint i = M; i >= 1; i--)
    {
        sum = a[i - 1];

        for(uint k = i + 1; k <= M; k++)
            sum -= A(k - 1, i - 1) * a[k - 1];

        a[i - 1] = sum / A(i - 1, i - 1);
    }
}

double CMathFunctions::ChiSq()
{
    double sum = 0, Chi2 = 0;

    for(uint i = 0; i < N; i++)
        d[i] = -y[i] / sigma[i];

    for(uint k = 0; k < M; k++)
    {
        for(uint j = 0; j < M; j++)
        {
            sum = 0;
            for(uint i = 0; i < N; i++)
                sum += C(i, k) * C(i, j);

            A(k, j) = sum;
        }
    }

    for(uint k = 0; k < M; k++)
    {
        sum = 0;
        for(uint i = 0; i < N; i++)
            sum += y[i] * C(i, k);

        b[k] = sum;
    }

    CholDec();
    SolveEqSys();

    for(uint i = 0; i < N; i++)
    {
        for(uint j = 0; j < M; j++)
        {
            sum = 0;

            for(uint k = 0; k < M; k++)
                sum += C(i, k) * a[k];

            r[i] = sum + d[i];
        }
    }

    for(uint i = 0; i < N; i++)
        Chi2 += r[i] * r[i];

    return Chi2;
}
