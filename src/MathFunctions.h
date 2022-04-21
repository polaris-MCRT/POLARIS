#pragma once
#include "Faddeeva.hh"
#include "Matrix2D.h"
#include "Stokes.h"
#include "Typedefs.h"
#include <complex>
#include <limits>

#ifndef CMATHFUNCTIONS
#define CMATHFUNCTIONS

class spline
{
  public:
    spline()
    {
        N = 0;
        d = 0;
        u = 0;
        w = 0;
        p = 0;
        x = 0;
        y = 0;
    }

    spline(uint size)
    {
        N = size - 1;
        d = new double[size];
        u = new double[size];
        w = new double[size];
        p = new double[size];
        x = new double[size];
        y = new double[size];

        for(uint i = 0; i < size; i++)
        {
            d[i] = 0;
            u[i] = 0;
            w[i] = 0;
            p[i] = 0;
            x[i] = 0;
            y[i] = 0;
        }
    }

    ~spline()
    {
        if(d != 0)
            delete[] d;
        if(u != 0)
            delete[] u;
        if(w != 0)
            delete[] w;
        if(p != 0)
            delete[] p;
        if(x != 0)
            delete[] x;
        if(y != 0)
            delete[] y;
    }

    uint size() const
    {
        return N + 1;
    }

    void clear()
    {
        if(d != 0)
            delete[] d;
        d = 0;
        if(u != 0)
            delete[] u;
        u = 0;
        if(w != 0)
            delete[] w;
        w = 0;
        if(p != 0)
            delete[] p;
        p = 0;
        if(x != 0)
            delete[] x;
        x = 0;
        if(y != 0)
            delete[] y;
        y = 0;
    }

    void resize(uint size)
    {
        if(d != 0)
            delete[] d;
        if(u != 0)
            delete[] u;
        if(w != 0)
            delete[] w;
        if(p != 0)
            delete[] p;
        if(x != 0)
            delete[] x;
        if(y != 0)
            delete[] y;

        N = size - 1;
        d = new double[size];
        u = new double[size];
        w = new double[size];
        p = new double[size];
        x = new double[size];
        y = new double[size];

        for(uint i = 0; i < size; i++)
        {
            d[i] = 0;
            u[i] = 0;
            w[i] = 0;
            p[i] = 0;
            x[i] = 0;
            y[i] = 0;
        }
    }

    double getAverageY()
    {
        uint size = N + 1;
        double res = 0;

        for(uint i = 0; i < size; i++)
            res += y[i];

        res /= size;
        return res;
    }

    void setValue(uint pos, double _x, double _y)
    {
#ifdef DEBUG
        if(x == 0)
        {
            cout << "\nERROR: Spline was not initiated!" << endl;
            return;
        }
#endif
        x[pos] = _x;
        y[pos] = _y;
    }

    void addValue(uint pos, double _x, double _y)
    {
        x[pos] = _x;
        y[pos] += _y;
    }

    void setDynValue(double _x, double _y)
    {
        dx.push_back(_x);
        dy.push_back(_y);
    }

    /*
    void addYValueExt(uint pos, double _x, double _y)
    {
#ifdef DEBUG
        if(x == 0)
        {
            cout << "\nERROR: Spline was not initiated!" << endl;
            return;
        }
#endif

        x[pos] = _x;
        y[pos] += _y;
    }
    */

    void operator+=(spline spline2)
    {
        if(N + 1 == spline2.size())
            for(uint i = 0; i <= N; i++)
                y[i] += spline2.getValue(x[i]);
    }

    double f(double x) const
    {
        return x * x * x - x;
    }

    void printX()
    {
        cout << "spline" << endl;

        for(uint i = 0; i < N + 1; i++)
            cout << "   " << i << "  " << x[i] << endl; /**/
    }

    void createSpline()
    {
        if(N == 0)
            return;

        for(uint i = 1; i < N; i++)
            d[i] = 2 * (x[i + 1] - x[i - 1]);

        for(uint i = 0; i < N; i++)
        {
            // if((x[i + 1] - x[i]) == 0.0)
            //    cout << "\nERROR: Spline broken!" << endl;
            u[i] = x[i + 1] - x[i];
        }

        for(uint i = 1; i < N; i++)
        {
            w[i] = 6.0 * ((y[i + 1] - y[i]) / u[i] - (y[i] - y[i - 1]) / u[i - 1]);
        }

        for(uint i = 1; i < N - 1; i++)
        {
            w[i + 1] -= w[i] * u[i] / d[i];
            d[i + 1] -= u[i] * u[i] / d[i];
        }

        for(uint i = N - 1; i > 0; i--)
        {
            p[i] = (w[i] - u[i] * p[i + 1]) / d[i];
        }
    }

    void createDynSpline()
    {
        if(dx.size() <= 4)
            return;

        uint size = uint(dx.size());

        if(x != 0)
        {
            delete[] d;
            delete[] u;
            delete[] w;
            delete[] p;
            delete[] x;
            delete[] y;
        }

        N = size - 1;
        d = new double[size];
        u = new double[size];
        w = new double[size];
        p = new double[size];
        x = new double[size];
        y = new double[size];

        d[0] = 0.0;
        d[1] = 0.0;
        d[N] = 0.0;

        u[0] = dx[1] - dx[0];
        u[1] = 0.0;
        u[N] = 0.0;

        w[0] = 0.0;
        w[1] = 0.0;
        w[N] = 0.0;

        p[0] = 0.0;
        p[1] = 0.0;
        p[N] = 0.0;

        x[0] = dx[0];
        x[1] = dx[1];

        y[0] = dy[0];
        y[1] = dy[1];

        for(uint i = 1; i < N; i++)
        {
            x[i] = dx[i];
            y[i] = dy[i];

            x[i + 1] = dx[i + 1];
            y[i + 1] = dy[i + 1];

            d[i] = 2 * (x[i + 1] - x[i - 1]);

            u[i] = x[i + 1] - x[i];

            if(u[i] == 0)
            {
                cout << "\nERROR:  Identical values in x[" << i << "] and x[" << i + 1
                     << "]!\n                        ";
                cout << u[i]
                     << "        Try smaller step sizes in Spline. \n                    "
                        "    ";
                u[i] = 1;
            }

            w[i] = 6.0 * ((y[i + 1] - y[i]) / u[i] - (y[i] - y[i - 1]) / u[i - 1]);
            // cout << i << " x:" << x[i] << " y:" << y[i]<< " d:" << d[i] << " u:" <<
            // u[i] << " w:" << w[i] << endl;
        }

        for(uint i = 1; i < N - 1; i++)
        {
            w[i + 1] -= w[i] * u[i] / d[i];
            d[i + 1] -= u[i] * u[i] / d[i];

            // cout << i << " " << w[i] << " " << d[i] << endl;
        }

        for(uint i = N - 1; i > 0; i--)
            p[i] = (w[i] - u[i] * p[i + 1]) / d[i];
    }

    double getMinValue()
    {
        return y[0];
    }

    double getMaxValue()
    {
        return y[N];
    }

    double getLinear(uint i, double v) const
    {
        double t = v - x[i];
        return y[i] + t * (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
    }

    double getLogLinear(uint i, double v) const
    {
        double x1 = log10(x[i]);
        double x2 = log10(x[i+1]);

        double y1 = log10(y[i]);
        double y2 = log10(y[i+1]);

        v = log10(v);

        double res = y1 + ((y2-y1)/(x2-x1)) * (v-x1);
        res = pow(10.0, res);

        return res;
    }

    double getLinearValue(double v)
    {
        uint min = 0, max = N;

//        for(int i=0;i<=max;i++)
//            cout << x[i] << "\t" << y[i]<< endl;

        if(x == 0)
            return 0;

        if(v <= x[0])
           return getLinear(0, v);

        if(v >= x[N])
            return getLinear(N - 1, v);

        while(max - min > 1)
        {
            uint i = min + (max - min) / 2;
            if(x[i] >= v)
                max = i;
            else
                min = i;
        }

        double x1 = log10(x[min]);
        double x2 = log10(x[min+1]);

        double y1 = log10(y[min]);
        double y2 = log10(y[min+1]);

        v = log10(v);

        double res = y1 + ((y2-y1)/(x2-x1)) * (v-x1);
        res = pow(10.0, res);
        return res;
    }

    double getValue(double v, uint extrapolation = SPLINE) const
    {
        uint min = 0;

        if(x == 0)
            return 0;

        if(N == 0)
            return y[0];

        // Extrapolation
        if(v < x[0])
            switch(extrapolation)
            {
                case CONST:
                    return y[0];
                    break;

                case LOGLINEAR:
                    return getLogLinear(0, v);
                    break;

                case LINEAR:
                    return getLinear(0, v);
                    break;

                case SPLINE:
                    min = 0;
                    break;
            }
        else if(v > x[N])
            switch(extrapolation)
            {
                case CONST:
                    return y[N];
                    break;

                case LOGLINEAR:
                    return getLogLinear(N - 1, v);
                    break;

                case LINEAR:
                    return getLinear(N - 1, v);
                    break;

                case SPLINE:
                    min = N - 1;
                    break;
            }
        else if(v == x[0])
            min = 0;
        else
            min = lower_bound(x, x+N+1, v) - x - 1;

        double t = (v - x[min]) / (u[min]);

        return t * y[min + 1] + (1 - t) * y[min] +
               u[min] * u[min] * (f(t) * p[min + 1] + f(1 - t) * p[min]) / 6.0;
    }

    double getLimitedValue(double v, double limit)
    {
        double t; // int i=1;
        unsigned int min_ID = 0, max_ID = N;

        if(v < x[0] || v > x[N])
            return 0;

        while(max_ID - min_ID > 1)
        {
            uint i = min_ID + (max_ID - min_ID) / 2;
            if(x[i] > v)
                max_ID = i;
            else
                min_ID = i;
        }

        double min_val, max_val;

        if(min_ID == 0)
        {
            min_val = min(y[0], y[1]);
            max_val = max(y[0], y[1]);
        }
        else if(min_ID == N)
        {
            min_val = min(y[N - 1], y[N]);
            max_val = max(y[N - 1], y[N]);
        }
        else
        {
            min_val = min(y[min_ID], y[min_ID + 1]);
            max_val = max(y[min_ID], y[min_ID + 1]);
        }

        // double lim_min = limit*min_val;
        double lim_max = (1 + limit) * max_val;

        t = (v - x[min_ID]) / (u[min_ID]);

        double res = t * y[min_ID + 1] + (1 - t) * y[min_ID] +
                     u[min_ID] * u[min_ID] * (f(t) * p[min_ID + 1] + f(1 - t) * p[min_ID]) / 6.0;

        double diff, lim_f;

        if(res < min_val)
        {
            res = min_val;
        }

        if(res > max_val)
        {
            diff = res - max_val;
            lim_f = diff / lim_max;

            if(lim_f > 1)
                lim_f = 1;

            diff *= lim_f;

            res = max_val - diff;
        }

        return res;
    }

    double getX(uint pos)
    {
        return x[pos];
    }

    double getY(uint pos)
    {
        return y[pos];
    }

    uint getXIndex(double v)
    {
        if(v < x[0] || v > x[N] || N==1)
            return 0;

        uint min = upper_bound(x, x+N+1, v) - x - 1;

        return min;
    }

    uint getYIndex(double v) const
    {
        if(v < y[0] || v > y[N] || N==1)
            return 0;

        uint min = upper_bound(y, y+N+1, v) - y - 1;

        return min;
    }

    double getValue(uint i) const
    {
        if(i < 0 || i > N)
            return 0;

        return y[i];
    }

  private:
    uint N;
    double * d;
    double * u;
    double * w;
    double * p;
    double * x;
    double * y;

    dlist dx, dy;
};

class interp
{
  public:
    interp()
    {
        N = 0;
    }

    interp(uint size)
    {
        N = size - 1;
        x.resize(size);
        y.resize(size);
    }

    uint size() const
    {
        return N + 1;
    }

    void resize(uint size)
    {
        N = size - 1;
        x.resize(size);
        y.resize(size);
    }

    void setValue(uint pos, double _x, double _y)
    {
#ifdef DEBUG
        if(x == 0)
        {
            cout << "\nERROR: Linear interpolation was not initiated!" << endl;
            return;
        }
#endif
        x[pos] = _x;
        y[pos] = _y;
    }

    void addValue(double _x, double _y)
    {
        x.push_back(_x);
        y.push_back(_y);
    }

    double getLinear(uint i, double v) const
    {
        double t = v - x[i];
        return y[i] + t * (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
    }

    double getValue(double v, uint interpolation = LINEAR) const
    {
        if(N == 0)
            return y[0];

        if(v < x[0])
            switch(interpolation)
            {
                case CONST:
                    return y[0];
                    break;

                case LINEAR:
                    return getLinear(0, v);
                    break;
            }
        else if(v > x[N])
            switch(interpolation)
            {
                case CONST:
                    return y[N];
                    break;

                case LINEAR:
                    return getLinear(N - 1, v);
                    break;
            }
        else
        {
            uint min = 0;

            if(v != x[0])
                min = lower_bound(x.begin(), x.end(), v) - x.begin() - 1;

            switch(interpolation)
            {
                case CONST:
                    return y[min+1];
                    break;

                case LINEAR:
                    return getLinear(min, v);
                    break;
            }
        }

        return 0;
    }

  private:
    uint N;
    dlist x;
    dlist y;
};

class prob_list
{
  public:
    prob_list()
    {
        N = 0;
        x = 0;
    }

    prob_list(uint size)
    {
        N = size - 1;
        x = new double[size];

        for(uint i = 0; i < size; i++)
            x[i] = 0;
    }

    ~prob_list()
    {
        if(x != 0)
            delete[] x;
    }

    uint size()
    {
        return N + 1;
    }

    double X(uint i)
    {
        return x[i];
    }

    void resize(uint size)
    {
        if(x != 0)
            delete[] x;

        N = size - 1;
        x = new double[size];

        for(uint i = 0; i < size; i++)
            x[i] = 0;
    }

    void setValue(uint pos, double _x)
    {
#ifdef DEBUG
        if(x == 0)
        {
            cout << "\nERROR: Spline was not initiated!" << endl;
            return;
        }
#endif
        x[pos] = _x;
    }

    uint getIndex(double v)
    {
        uint min = 0;

        if(v < x[0] || v > x[N] || N == 1)
            return 0;

        min = upper_bound(x, x+N+1, v) - x - 1;

        return min;
    }

    void normalize(double integValue)
    {
        if(integValue > 0)
            for(uint i = 0; i <= N; i++)
                x[i] /= integValue;
        else
            resize(N + 1);
    }

  private:
    uint N;
    double * x;
};

namespace
{
inline prob_list operator-(prob_list & list1, prob_list & list2)
{
    if(list1.size() != list2.size())
        cout << "\nERROR: Probability lists have different lengths!";
    prob_list diff_list(list1.size());

    for(uint i = 0; i < list1.size(); i++)
        diff_list.setValue(i, list1.X(i) - list2.X(i));

    return diff_list;
}

inline spline operator*(double val, spline & spline2)
{
    spline res_spline(spline2.size());

    for(uint i = 0; i < spline2.size(); i++)
        res_spline.setValue(i, spline2.getX(i), spline2.getY(i) * val);
    return res_spline;
}
}

class CRandomGenerator
{
  public:
    CRandomGenerator()
    {
        // The standard seed values as proposed by George Marsaglia
        // These will be used is init is not called
        // x
        KISS_state[0] = 1234567890987654321ULL;
        // y
        KISS_state[1] = 362436362436362436ULL;
        // z
        KISS_state[2] = 1066149217761810ULL;
        // c
        KISS_state[3] = 123456123456123456ULL;
    }

    void init(ullong seed)
    {
        // KISS will work just fine even w/o the init process
        KISS_state[0] = seed;
        // The seed for CONG is given by the OMP thread number
        // call CONG several times to get a nice random seed for KISS
        for(int i = 0; i < 2000; i++)
            KISS_state[0] = CONG(KISS_state[0]);

        // fill the state array for KISS
        KISS_state[1] = CONG(KISS_state[0]);
        KISS_state[2] = CONG(KISS_state[1]);
        KISS_state[3] = CONG(KISS_state[2]);
    }

    ullong CONG(ullong current_state)
    {
        // CONG - very simple linear congruential generator
        return 6906969069ULL * current_state + 1234567;
    }

    double getRND()
    {
        // KISS (Keep it Simple Stupid) is a family of pseudorandom number generators
        // introduced by George Marsaglia.
        // Source: https://de.wikipedia.org/wiki/KISS_(Zufallszahlengenerator)

        // linear congruential generator
        KISS_state[2] = CONG(KISS_state[2]);

        // Xorshift
        KISS_state[1] ^= KISS_state[1] << 13;
        KISS_state[1] ^= KISS_state[1] >> 17;
        KISS_state[1] ^= KISS_state[1] << 43;

        // Multiply-with-carry
        ullong tmp = (KISS_state[0] << 58) + KISS_state[3];
        KISS_state[3] = (KISS_state[0] >> 6);
        KISS_state[0] += tmp;
        KISS_state[3] += (KISS_state[0] < tmp);

        // Return double between 0 and 1
        return double(KISS_state[0] + KISS_state[1] + KISS_state[2]) / 18446744073709551615ULL;
    }

    double getRNDnormal(double mu, double sigma)
    {
        double U1, U2, W, mult;
        double X1, X2;

        do
        {
            U1 = -1 + getRND() * 2;
            U2 = -1 + getRND() * 2;
            W = pow(U1, 2) + pow(U2, 2);
        } while(W >= 1 || W == 0);

        mult = sqrt((-2 * log(W)) / W);
        X1 = U1 * mult;
        X2 = U2 * mult;

        double res = mu + sigma * X1;

        if(res < 0)
            return getRNDnormal(mu, sigma);

        return res;
    }

  private:
    ullong KISS_state[4];
};

class CMathFunctions
{
  public:
    CMathFunctions(void)
    {
        b = 0;
        a = 0;
        r = 0;
        d = 0;
        y = 0;
        sigma = 0;

        nr_ofSeq = 0;
    }

    ~CMathFunctions(void)
    {
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
    }

    static inline bool isPowerOfTwo(int num)
    {
        return ((num & (num - 1)) == 0);
    }

    static inline Vector3D SolveEqSys(Matrix2D & inA, Vector3D b)
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

    static inline Vector3D gauss(Matrix2D inA, Vector3D inb)
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

    static inline void gauss(Matrix2D & A, double * b, int n)
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

    static inline double BE_Mass(double r, double rho0, double rc)
    {
        double res = PIx4 * rho0 * rc * rc * (r - rc * atan(r / rc));
        return res;
    }

    static inline void CholDec(Matrix2D & A)
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

    void setHeader(uint pos,
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

    void updateSatistics(uint pos, StokesVector st)
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

    /*     static inline double min(double * values, uint size)
     {
     double res_min=1e300;

     for(uint i=0; i<size; i++)
     {
     if(values[i]<res_min) res_min=values[i];
     }

     return res_min;
     }

     static inline double max(double * values, uint size)
     {
     double res_max=-1e300;

     for(uint i=0; i<size; i++)
     {
     if(values[i]>res_max) res_max=values[i];
     }

     return res_max;
     }

     static inline double mean(double * values, uint size)
     {
     double res_mean=0;

     for(uint i=0; i<size; i++)
     {
     res_mean+=values[i];
     }

     res_mean/=double(size);

     return res_mean;
     }*/

    /*static inline Vector3D sphere_dist(ulong N, ulong k)
    {
        Vector3D res;
        double inc = 3.1415 * (3.0 - sqrt(5.0));
        double off = 2.0 / double(N);
        double y = k * off - 1 + (off / 2);
        double r = sqrt(1 - y * y);
        double phi = k * inc;
        res.set(cos(phi) * r, y, sin(phi) * r);
        res.normalize();
        return res;
    }*/

    static inline int insertInList(double val, dlist & list)
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

    static inline uint biListIndexSearch(double val, const dlist & list)
    {
        uint N = uint(list.size());

        if(val < list[0] || val > list[N - 1])
            return MAX_UINT;

        if(val == list[0])
            return 0;

        uint min = lower_bound(list.begin(), list.end(), val) - list.begin() - 1;

        return min;
    }

    static inline uint biListIndexSearchRec(double val, const dlist & list)
    {
        uint N = uint(list.size());

        if(val < list[0] || val > list[N - 1])
            return MAX_UINT;

        uint min = upper_bound(list.begin(), list.end(), val) - list.begin() - 1;

        return min;
    }

    static inline uint biListIndexSearchRec(double val, const double * list, uint N)
    {
        if(val < list[0] || val > list[N - 1])
            return MAX_UINT;

        uint min = upper_bound(list, list+N, val) - list - 1;

        return min;
    }

    static inline uint biListIndexSearch(double val, const double * list, uint N)
    {
        if(val < list[0] || val > list[N - 1])
            return MAX_UINT;

        if(val == list[0])
            return 0;

        uint min = lower_bound(list, list+N, val) - list - 1;

        return min;
    }

    static inline double interpolate(double x_min, double x_max, double y_min, double y_max, double x_ipol)
    {
        double ipol2_result;
        if(x_min == x_max)
            ipol2_result = y_min;
        else
            ipol2_result = (y_max - y_min) * (x_ipol - x_min) / (x_max - x_min) + y_min;

        return ipol2_result;
    }

    static inline uint inSphereTest(Vector3D a, Vector3D b, Vector3D c, Vector3D d, Vector3D e)
    {

        return 0;
    }

    static inline uint orientationTest(Vector3D a, Vector3D b, Vector3D c, Vector3D d)
    {

        return 0;
    }

    static inline double det5x5(double mat[5][5], int N = 5)
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

    static inline double det4x4(double mat[4][4], int N = 4)
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

    static inline Vector3D calcKeplerianVelocity(Vector3D pos, double stellar_mass)
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

    static inline double getClosestLinePoint(Vector3D line_start, Vector3D line_end, Vector3D point)
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

    static inline double maxValue(double v1, double v2)
    {
        if(v1 > v2)
            return v1;

        return v2;
    }

    static inline uint maxValue(uint v1, uint v2)
    {
        if(v1 > v2)
            return v1;

        return v2;
    }

    static inline double grad2rad(double grad)
    {

        return grad * PI / 180.0;
    }

    static inline double rad2grad(double rad)
    {

        return rad * 180.0 / PI;
    }

    static inline bool areSame(double a, double b)
    {
        if(fabs(a - b) < 1e-5 * fabs(a + b))
            return true;
        return false;
    }

    static inline double generateGaussianNoise(const double & variance)
    {
        static bool hasSpare = false;
        static double rand1, rand2;

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

    static inline double diffSeries(double y, double upper_limit)
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

    void initErfi();

    static inline double getErfi(double x)
    {
        double res = Faddeeva::erfi(x);
        if(res < 0)

            return 0;
        return res;
    }

    void sleep(int milliseconds)
    {
#ifdef WINDOWS
        Sleep(milliseconds);
#else

        usleep(milliseconds * 1000);
#endif
    }

    static inline double integ(const double * x, const double * y, uint xlow, uint xup)
    {
        double res = 0;
        if(xlow != xup)
        {
            for(uint i = xlow + 1; i <= xup; i++)
                res += (x[i] - x[i - 1]) * y[i - 1] + 0.5 * (x[i] - x[i - 1]) * (y[i] - y[i - 1]);
        }

        return res;
    }

    static inline double integ(const dlist & x, double * y, uint xlow, uint xup)
    {
        double res = 0;
        if(xlow != xup)
            for(uint i = xlow + 1; i <= xup; i++)
                res += (x[i] - x[i - 1]) * y[i - 1] + 0.5 * (x[i] - x[i - 1]) * (y[i] - y[i - 1]);

        return res;
    }

    // for debugging reasons only
    static inline double integ1(dlist x, double * y, uint xlow, uint xup)
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

    static inline double integ(const dlist & x, const dlist & y, uint xlow, uint xup)
    {
        double res = 0;
        if(xlow != xup)
            for(uint i = xlow + 1; i <= xup; i++)
                res += (x[i] - x[i - 1]) * y[i - 1] + 0.5 * (x[i] - x[i - 1]) * (y[i] - y[i - 1]);

        return res;
    }

    static inline double integ(const dlist & x, const uilist & y, uint xlow, uint xup)
    {
        double res = 0;
        if(xlow != xup)
            for(uint i = xlow + 1; i <= xup; i++)
                res += (x[i] - x[i - 1]) * y[i - 1] + 0.5 * (x[i] - x[i - 1]) * (y[i] - y[i - 1]);

        return res;
    }

    static inline double integ(const double * x, dlist & y, uint xlow, uint xup)
    {
        double res = 0;
        if(xlow != xup)
            for(uint i = xlow + 1; i <= xup; i++)
                res += (x[i] - x[i - 1]) * y[i - 1] + 0.5 * (x[i] - x[i - 1]) * (y[i] - y[i - 1]);

        return res;
    }

    static inline double integ(double * x, dlist & y, uint xlow, uint xup)
    {
        double res = 0;
        if(xlow != xup)
            for(uint i = xlow + 1; i <= xup; i++)
                res += (x[i] - x[i - 1]) * y[i - 1] + 0.5 * (x[i] - x[i - 1]) * (y[i] - y[i - 1]);

        return res;
    }

    static inline double integ_dust_size(const double * a_eff,
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

    static inline StokesVector integ_dust_size(const double * a_eff,
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

    static inline double full_integ(dlist x, dlist y, double xlow, double xup)
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

    static inline void probListInteg(const double * x, const double * y, prob_list & integ_spline)
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

    static inline void probListInteg(dlist x, const double * y, prob_list & integ_spline)
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

    static inline double calc_delta(double B, double Td, double Tg, double ng)
    {
        double den = ng * Td * sqrt(Tg);

        if(den == 0)
            return 1;

        return 1.0 * B * B / den;
    }

    static inline double calc_mach(double vel, double Tg, double mu)
    {
        return vel / sqrt(con_kB * Tg / (mu * m_H));
    }

    static inline double calc_larm_limit(double B, double Td, double Tg, double ng, double s, double larm_f)
    {
        double den = 1.0 * ng * Td * sqrt(Tg);

        if(den == 0)
            return 1;

        return s * s * B / den / larm_f;
    }

    static inline double planck(double l, double T)
    {
        return (2.0 * con_h * con_c * con_c) /
               ((l * l * l * l * l) * (exp((con_h * con_c) / (l * con_kB * T)) - 1));
    }

    static inline double planck_hz(double f, double T)
    {

        return (2.0 * con_h * f * f * f) / ((con_c * con_c) * (exp((con_h * f) / (con_kB * T)) - 1));
    }

    static inline double dplanck_dT(double l, double T)
    {
        double ex = exp((con_h * con_c) / (l * con_kB * T));
        ex = (2.0 * con_c * con_c * con_c * ex * con_h * con_h) /
             ((ex - 1) * (ex - 1) * con_kB * l * l * l * l * l * l * T * T);

        if(ex != ex)
            return 0;

        return ex;
    }

    static inline double mathis_isrf(double wavelength)
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

    static inline void SinList(double start, double stop, double * list, uint N, double f)
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

    static inline void LinearList(double start, double stop, double * list, uint N)
    {
        double dx = (stop - start) / (N - 1);

        list[0] = start;

        for(uint i_x = 1; i_x < N - 1; i_x++)
            list[i_x] = start + i_x * dx;

        list[N - 1] = stop;
    }

    static inline dlist LinearList(double start, double stop, uint N)
    {
        dlist list(N);

        double dx = (stop - start) / (N - 1);

        list[0] = start;

        for(uint i_x = 1; i_x < N - 1; i_x++)
            list[i_x] = start + i_x * dx;

        list[N - 1] = stop;

        return list;
    }

    static inline void ExpList(double start, double stop, double * list, uint N, double base)
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

    static inline void ExpListSym(double start, double stop, double * list, uint N, double base)
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

    static inline void LogList(double start, double stop, double * list, uint N, double base)
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

    static inline void LogList(double start, double stop, dlist & list, double base)
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

    static inline void SymLogList(double start, double stop, double * list, uint N, double base)
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

    static inline double sgn(double x)
    {
        if(x < 0)
            return -1.0;

        else
            return 1.0;

        return 0;
    }

    static inline double sgn(int x)
    {
        if(x < 0)
            return -1.0;

        else
            return 1.0;
        return 0;
    }

    static inline Matrix2D getRotationMatrix(double cos_phi,
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

    static inline double getRotationAngleObserver(const Vector3D & obs_ex,
                                                  const Vector3D & photon_ex,
                                                  const Vector3D & photon_ey)
    {
        double cos_angle_1 = obs_ex * photon_ey;
        double cos_angle_2 = obs_ex * photon_ex;

        return atan3(cos_angle_1, cos_angle_2);
    }

    static inline void getPropMatrixAPi(double cos_theta,
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

    static inline void getPropMatrixBSigmaP(double cos_theta,
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

    static inline void getPropMatrixASigmaP(double cos_theta,
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

    static inline void getPropMatrixASigmaM(double cos_theta,
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

    static inline void getPropMatrixBSigmaM(double cos_theta,
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

    static inline void getPropMatrixBPi(double cos_theta,
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

    static inline Vector3D interpVector(Vector3D y_1, Vector3D y_2, double x_1, double x_2)
    {
        return (y_2 - y_1) / (x_2 - x_1);
    }

    static inline double lum2Jy(double intensity, double wavelength, double distance)
    {
        if(abs(intensity) < 1e-47)
            return 0;

        double freq = con_c / wavelength;             // [s^-1]
        double L = intensity * con_c / (freq * freq); // [W Hz^-1 sr^-1]
        return L * 1e+26 / (distance * distance);     // [Jy]
    }

    static inline void lum2Jy(StokesVector * S, double wavelength, double distance)
    {
        S->setI(lum2Jy(S->I(), wavelength, distance));
        S->setQ(lum2Jy(S->Q(), wavelength, distance));
        S->setU(lum2Jy(S->U(), wavelength, distance));
        S->setV(lum2Jy(S->V(), wavelength, distance));
    }

    static inline uint findListIndex(double * list, uint l_limit, uint h_limit, double val)
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

    void initChiSq(uint _N, uint _M)
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

    static inline double getPolarInterp(double t, double s, double Q1, double Q2, double Q3, double Q4)
    {
        return Q1 * (1 - s) * (1 - t) + Q2 * s * (1 - t) + Q3 * (1 - s) * t + Q4 * s * t;
    }

    static inline double phaseFunctionHG(double g, double theta)
    {
        double res = 1 / PIx4;
        res *= (1 - g * g);
        res /= pow(1 + g * g - 2 * g * cos(theta), 1.5);

        return res;
    }

    static inline double Freq2Velo(double _f, double _f0)
    {

        return con_c * _f / _f0;
    }

    static inline double Velo2Freq(double _v, double _f0)
    {

        return _v / con_c * _f0;
    }

    static inline double getWavelengthStepWidth(double lam_min, double lam_max, double nr_of_spectral_bins)
    {
        double dw;
        if(nr_of_spectral_bins == 1)
            dw = 0;
        else
            dw = (lam_max - lam_min) / double(nr_of_spectral_bins - 1);
        return dw;
    }

    static inline void gauss(Matrix2D & inA, double * b, double * res, uint n)
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

    static inline dlist gauss_solver(Matrix2D & A, uint n, dlist & b)
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

    static inline double fmodulo(double v1, double v2)
    {
        /*! Returns the remainder of the division \a v1/v2.
        The result is non-negative.
        \a v1 can be positive or negative; \a v2 must be positive. */
        if(v1 >= 0)
            return (v1 < v2) ? v1 : fmod(v1, v2);
        double tmp = fmod(v1, v2) + v2;
        return (tmp == v2) ? 0. : tmp;
    }

    static inline int imodulo(int v1, int v2)
    {
        /*! Returns the remainder of the division \a v1/v2.
        The result is non-negative.
        \a v1 can be positive or negative; \a v2 must be positive. */
        int v = v1 % v2;
        return (v >= 0) ? v : v + v2;
    }

    static inline dcomplex Csqrt(dcomplex z)
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

    /*static bool calcBHMie(double x,
                          dcomplex refractive_index,
                          double & qext,
                          double & qabs,
                          double & qsca,
                          double & gsca,
                          double * S11,
                          double * S12,
                          double * S33,
                          double * S34)
    // Subroutine BHMIE is the Bohren-Huffman Mie scattering subroutine
    // to calculate scattering and absorption by a homogenous isotropic
    // sphere.

    // Comment:
    //     NANG = number of angles between 0 and 90 degrees
    //             (will calculate 2 * NANG - 1 directions from 0 to 180 deg.)

    // Given:
    //     X = 2*pi*a/lambda
    //     REFREL = (complex refractive index of sphere) / (real index of medium)

    // Returns:
    //     S1(1 .. 2 * NANG - 1) =  (incident E perpendicular to scattering plane,
    //                               scattering E perpendicular to scattering plane)
    //     S2(1 .. 2 * NANG - 1) =  (incident E parallel to scattering plane,
    //                               scattering E parallel to scattering plane)
    //     QEXT = C_ext/pi*a**2 = efficiency factor for extinction
    //     QSCA = C_sca/pi*a**2 = efficiency factor for scattering
    //     QBACK = 4*pi*(dC_sca/domega)/pi*a**2
    //         = backscattering efficiency
    //     GSCA = <cos(theta)> for scattering

    // Original program taken from Bohren and Huffman (1983), Appendix A
    // Modified by B.T.Draine, Princeton Univ. Obs., 90/10/26
    // in order to compute <cos(theta)>
    {
        dcomplex cxy = dcomplex(x, 0) * refractive_index;

        // Series expansion terminated after XSTOP terms
        float xstop = x + 4.0 * pow(x, 1.0 / 3.0) + 2.0;
        long nmx = fmax(xstop, abs(cxy)) + 15;

        if(nmx >= MAX_MIE_ITERATIONS)
        {
            cout << "\nERROR: Failure in Mie-scattering calculation (NMX = " << nmx
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
    }*/

    static bool calcWVMie(double x,
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
    {
        // Step width
        uint n_scat_angle = scat_angle.size();
        double factor = 1e250;

        if(x <= MIN_MIE_SIZE_PARAM)
        {
            cout << "\nError: Mie scattering limit exceeded, current size parameter: " << x << "\n" << endl;
            return false;
        }

        double ax = 1 / x;
        double b = 2 * pow(ax, 2);
        dcomplex ss(0, 0);
        dcomplex s3(0, -1);
        double an = 3;

        // choice of number for subroutine aa [Loskutov (1971)]
        double y = abs(refractive_index) * x;
        uint num = uint(1.25 * y + 15.5);
        if(y < 1)
            num = uint(7.5 * y + 9.0);
        else if(y > 100 && y < 50000)
            num = uint(1.0625 * y + 28.5);
        else if(y >= 50000)
            num = uint(1.005 * y + 50.5);

        if(num >= MAX_MIE_ITERATIONS - 1)
        {
            cout << "\nError: Maximum number of terms : " << MAX_MIE_ITERATIONS << endl;
            cout << "       Number of terms required: " << num << endl;
            cout << "       Increase default value of the variable MAX_MIE_ITERATIONS in Typedefs.h\n" << endl;
            return false;
            // return calcGeometricOptics(x, refractive_index, qext, qabs,
            //    qabs, gsca, S11, S12, S33, S34);
        }

        // logarithmic derivative to Bessel function (complex argument)
        dcomplex *ru = new dcomplex[num + 1];
        dcomplex s_tmp = ax / refractive_index;
        ru[num] = dcomplex(num + 1, 0) * s_tmp;
        for(uint n = 1; n <= num - 1; n++)
        {
            double rn = double(num - n);
            dcomplex s1 = (rn + 1) * s_tmp;
            ru[num - n] = s1 - dcomplex(1, 0) / (ru[num - n + 1] + s1);
        }

        // initialize term counter
        uint iterm = 1;

        // Bessel functions
        double ass = sqrt(PI2 * ax);
        double w1 = invPI2 * ax;
        double Si = sin(x) / x;
        double Co = cos(x) / x;

        // n=0
        double besJ0 = Si / ass;
        double besY0 = -Co / ass;
        uint iu0 = 0;

        // n=1
        double besJ1 = (Si * ax - Co) / ass;
        double besY1 = (-Co * ax - Si) / ass;
        uint iu1 = 0;
        uint iu2 = 0;

        // Mie coefficients (first term)
        dcomplex s, s1, s2, ra0, rb0;

        // coefficient a_1
        s = ru[iterm] / refractive_index + ax;
        s1 = s * besJ1 - besJ0;
        s2 = s * besY1 - besY0;
        ra0 = s1 / (s1 - s3 * s2);

        // coefficient b_1
        s = ru[iterm] * refractive_index + ax;
        s1 = s * besJ1 - besJ0;
        s2 = s * besY1 - besY0;
        rb0 = s1 / (s1 - s3 * s2);

        // efficiency factors (first term)
        dcomplex r = -1.5 * (ra0 - rb0);
        qext = an * real(ra0 + rb0);
        qsca = an * (norm(ra0) + norm(rb0));

        // first term (iterm=1)
        double r_iterm = double(iterm);
        double FN = (2 * r_iterm + 1) / (r_iterm * (r_iterm + 1));

        dlist dPI(n_scat_angle), dTAU(n_scat_angle);
        dlist dAMU(n_scat_angle), dPI0(n_scat_angle, 0), dPI1(n_scat_angle, 1);

        dcomplex SM1[n_scat_angle], SM2[n_scat_angle];

        for(uint i_scat_ang = 0; i_scat_ang < n_scat_angle; i_scat_ang++)
        {
            dAMU[i_scat_ang] = cos(scat_angle[i_scat_ang]);

            SM1[i_scat_ang] = dcomplex(0, 0);
            SM2[i_scat_ang] = dcomplex(0, 0);

            dTAU[i_scat_ang] = r_iterm * dAMU[i_scat_ang] * dPI1[i_scat_ang] - (r_iterm + 1) * dPI0[i_scat_ang];

            SM1[i_scat_ang] = SM1[i_scat_ang] + FN * (ra0 * dPI1[i_scat_ang] + rb0 * dTAU[i_scat_ang]);
            SM2[i_scat_ang] = SM2[i_scat_ang] + FN * (ra0 * dTAU[i_scat_ang] + rb0 * dPI1[i_scat_ang]);

            dPI[i_scat_ang] = dPI1[i_scat_ang];
            dPI1[i_scat_ang] *= (2 + 1 / r_iterm) * dAMU[i_scat_ang];
            dPI1[i_scat_ang] -= dPI0[i_scat_ang] * (1 + 1 / r_iterm);
            dPI0[i_scat_ang] = dPI[i_scat_ang];
        }

        iterm++;
        r_iterm = double(iterm);

        double z = -1, besY2, besJ2, an2, qq;
        dcomplex ra1, rb1, rr;
        while(true)
        {
            // if(iterm % 1000 == 0)
            //    cout << x << TAB << iterm << endl;
            an = an + 2;
            an2 = an - 2;

            // Bessel functions
            if(iu1 == iu0)
                besY2 = an2 * ax * besY1 - besY0;
            else
                besY2 = an2 * ax * besY1 - besY0 / factor;

            if(abs(besY2) > 1e200)
            {
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

            dcomplex ru_tmp(0, 0);
            if(iterm <= num)
                ru_tmp = ru[iterm];

            s = ru_tmp / refractive_index + r_iterm * ax;
            s1 = s * (besJ2 / factorial(iu2)) - besJ1 / factorial(iu1);
            s2 = s * (besY2 * factorial(iu2)) - besY1 * factorial(iu1);
            ra1 = s1 / (s1 - s3 * s2); // coefficient a_n, (n=iterm)

            s = ru_tmp * refractive_index + r_iterm * ax;
            s1 = s * (besJ2 / factorial(iu2)) - besJ1 / factorial(iu1);
            s2 = s * (besY2 * factorial(iu2)) - besY1 * factorial(iu1);
            rb1 = s1 / (s1 - s3 * s2); // coefficient b_n, (n=iterm)

            // efficiency factors
            z = -z;
            rr = z * (r_iterm + 0.5) * (ra1 - rb1);
            r = r + rr;
            ss = ss + (r_iterm - 1) * (r_iterm + 1) / r_iterm * (ra0 * conj(ra1) + rb0 * conj(rb1)) +
                 an2 / r_iterm / (r_iterm - 1) * (ra0 * conj(rb0));
            qq = an * real(ra1 + rb1);
            qext = qext + qq;
            qsca = qsca + an * (norm(ra1) + norm(rb1));

            // leaving-the-loop with error criterion
            if(isnan(qext))
                return false;

            // leaving-the-loop criterion
            if(abs(qq / qext) < MIE_ACCURACY)
                break;

            // Bessel functions
            besJ0 = besJ1;
            besJ1 = besJ2;
            besY0 = besY1;
            besY1 = besY2;
            iu0 = iu1;
            iu1 = iu2;
            ra0 = ra1;
            rb0 = rb1;

            // terms iterm=2,...
            r_iterm = double(iterm);
            FN = (2 * r_iterm + 1) / (r_iterm * (r_iterm + 1));
            for(uint i_scat_ang = 0; i_scat_ang < n_scat_angle; i_scat_ang++)
            {
                dTAU[i_scat_ang] = r_iterm * dAMU[i_scat_ang] * dPI1[i_scat_ang] - (r_iterm + 1) * dPI0[i_scat_ang];

                SM1[i_scat_ang] = SM1[i_scat_ang] + FN * (ra0 * dPI1[i_scat_ang] + rb0 * dTAU[i_scat_ang]);
                SM2[i_scat_ang] = SM2[i_scat_ang] + FN * (ra0 * dTAU[i_scat_ang] + rb0 * dPI1[i_scat_ang]);

                dPI[i_scat_ang] = dPI1[i_scat_ang];
                dPI1[i_scat_ang] *= (2 + 1 / r_iterm) * dAMU[i_scat_ang];
                dPI1[i_scat_ang] -= dPI0[i_scat_ang] * (1 + 1 / r_iterm);
                dPI0[i_scat_ang] = dPI[i_scat_ang];
            }

            iterm++;
            r_iterm = double(iterm);

            if(iterm > num)
                break;
        }

        delete ru;

        // efficiency factors (final calculations)
        qext = b * qext;
        qsca = b * qsca;
        // qback = 2 * b * r * conj(r);
        double qpr = qext - 2 * b * real(ss);
        qabs = qext - qsca;
        gsca = (qext - qpr) / qsca;

        for(uint i_scat_ang = 0; i_scat_ang < n_scat_angle; i_scat_ang++)
        {
            S11[i_scat_ang] = 0.5 * (abs(SM2[i_scat_ang]) * abs(SM2[i_scat_ang]) + abs(SM1[i_scat_ang]) * abs(SM1[i_scat_ang]));
            S12[i_scat_ang] = 0.5 * (abs(SM2[i_scat_ang]) * abs(SM2[i_scat_ang]) - abs(SM1[i_scat_ang]) * abs(SM1[i_scat_ang]));
            S33[i_scat_ang] = real(SM2[i_scat_ang] * conj(SM1[i_scat_ang]));
            S34[i_scat_ang] = imag(SM2[i_scat_ang] * conj(SM1[i_scat_ang]));

            // if SM2 and SM1 get really large (if x >> 1 for instance)
            // then double precision may not be enough
            // to get S12/S34 = 0 for theta = 0/pi (cos(theta) = +-1)
            if(abs(dAMU[i_scat_ang]) == 1)
            {
                S12[i_scat_ang] = 0;
                S34[i_scat_ang] = 0;
            }
        }

        return true;
    }

    static inline bool calcGeometricOptics(double x,
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

    static inline double calcReflectionCoefficients(dcomplex refractive_index, double theta)
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

    static inline int factorial(int n)
    {
        return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
    }

    int IDUM;

  private:
    // statistics and SED

    Matrix2D stat_data;
    uint nr_ofSeq;

    // random number generator
    double RM;
    int IY, IFF;
    int IR[98];

    int Mr;
    int IA;
    int IC;

    // Chi squared

    double * b;
    double * a;
    double * r;
    double * d;
    double * y;
    double * sigma;
    Matrix2D A;
    Matrix2D C;
    uint N;
    uint M;

    ullong kiss_x;
    ullong kiss_y;
    ullong kiss_z;
    ullong kiss_c;

    void CholDec()
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

    void SolveEqSys()
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

    double ChiSq()
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
};

#endif
