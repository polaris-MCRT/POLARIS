#pragma once

#ifndef STOKES
#define STOKES
#include "Matrix2D.h"
#include "Typedefs.h"
#include "Vector.h"

class cross_sections
{
  public:
    cross_sections()
    {
        Cext = 0;
        Cpol = 0;
        Cabs = 0;
        Cpabs = 0;
        Ccirc = 0;
        Csca = 0;
    }

    void operator/=(double w)
    {
        Cext /= w;
        Cpol /= w;
        Cabs /= w;
        Cpabs /= w;
        Ccirc /= w;
        Csca /= w;
    }

    void operator*=(double w)
    {
        Cext *= w;
        Cpol *= w;
        Cabs *= w;
        Cpabs *= w;
        Ccirc *= w;
        Csca *= w;
    }

    cross_sections operator*(double w)
    {
        Cext *= w;
        Cpol *= w;
        Cabs *= w;
        Cpabs *= w;
        Ccirc *= w;
        Csca *= w;

        return *this;
    }

    void operator+=(const cross_sections & cs)
    {
        Cext += cs.Cext;
        Cpol += cs.Cpol;
        Cabs += cs.Cabs;
        Cpabs += cs.Cpabs;
        Ccirc += cs.Ccirc;
        Csca += cs.Csca;
    }

    void operator=(const double val)
    {
        Cext = val;
        Cpol = val;
        Cabs = val;
        Cpabs = val;
        Ccirc = val;
        Csca = val;
    }

    double Cext;
    double Cpol;
    double Cabs;
    double Cpabs;
    double Csca;
    double Ccirc;
};

class StokesVector
{
  public:
    StokesVector()
    {
        sI = 0;
        sQ = 0;
        sU = 0;
        sV = 0;
        sT = 0;
        sSp = 0;
    }

    StokesVector(double val)
    {
        sI = val;
        sQ = val;
        sU = val;
        sV = val;
        sT = val;
        sSp = val;
    }

    StokesVector(double I, double Q, double U, double V)
    {
        sI = I;
        sQ = Q;
        sU = U;
        sV = V;
        sT = 0;
        sSp = 0;
    }

    StokesVector(double I, double Q, double U, double V, double T)
    {
        sI = I;
        sQ = Q;
        sU = U;
        sV = V;
        sT = T;
        sSp = 0;
    }

    StokesVector(double I, double Q, double U, double V, double T, double Sp)
    {
        sI = I;
        sQ = Q;
        sU = U;
        sV = V;
        sT = T;
        sSp = Sp;
    }

    StokesVector(const StokesVector & st)
    {
        sI = st.I();
        sQ = st.Q();
        sU = st.U();
        sV = st.V();
        sT = st.T();
        sSp = st.Sp();
    }

    ~StokesVector(void)
    {}

    // Operators

    StokesVector & operator=(const StokesVector & ex)
    {
        sI = ex.I();
        sQ = ex.Q();
        sU = ex.U();
        sV = ex.V();
        sT = ex.T();
        sSp = ex.Sp();

        return *this;
    }

    StokesVector & operator=(double v)
    {
        sI = v;
        sQ = v;
        sU = v;
        sV = v;
        sT = v;
        sSp = v;
        return *this;
    }

    StokesVector operator+(const StokesVector & ex)
    {
        return StokesVector(sI + ex.I(), sQ + ex.Q(), sU + ex.U(), sV + ex.V(), sT + ex.T(), sSp + ex.Sp());
    }

    StokesVector operator-(const StokesVector & ex)
    {
        return StokesVector(sI - ex.I(), sQ - ex.Q(), sU - ex.U(), sV - ex.V(), sT - ex.T(), sSp - ex.Sp());
    }

    StokesVector & operator+=(const StokesVector & ex)
    {
        sI += ex.I();
        sQ += ex.Q();
        sU += ex.U();
        sV += ex.V();
        sT += ex.T();
        sSp += ex.Sp();
        return *this;
    }

    StokesVector & operator-=(const StokesVector & ex)
    {
        sI -= ex.I();
        sQ -= ex.Q();
        sU -= ex.U();
        sV -= ex.V();
        sT -= ex.T();
        sSp -= ex.Sp();
        return *this;
    }

    StokesVector & operator*=(double val)
    {
        sI *= val;
        sQ *= val;
        sU *= val;
        sV *= val;
        return *this;
    }

    StokesVector & operator/=(double val)
    {
        sI /= val;
        sQ /= val;
        sU /= val;
        sV /= val;
        return *this;
    }

    // linearly polarized intensity
    double iPol()
    {
        return sqrt(sU * sU + sQ * sQ);
    }

    // totaly polarized intensity
    double tPol()
    {
        return sqrt(sU * sU + sQ * sQ + sV * sV);
    }

    // degree of linear polarization
    double linPol()
    {
        if(sI != 0)
            return sqrt(sU * sU + sQ * sQ) / sI;

        return 0;
    }

    // degree of circular polarization
    double circPol()
    {
        if(sI != 0)
            return sV / sI;

        return 0;
    }

    // polarization angle
    double getAngle()
    {
        return 0.5 * angle(sQ, sU);
    }

    void setI(double _I)
    {
        sI = _I;
    }

    void setQ(double _Q)
    {
        sQ = _Q;
    }

    void setU(double _U)
    {
        sU = _U;
    }

    void setV(double _V)
    {
        sV = _V;
    }

    void setT(double _T)
    {
        sT = _T;
    }

    void setSp(double _Sp)
    {
        sSp = _Sp;
    }

    void set(double _I, double _Q, double _U, double _V, double _T)
    {
        sI = _I;
        sQ = _Q;
        sU = _U;
        sV = _V;
        sT = _T;
        sSp = 0;
    }

    void set(double _I, double _Q, double _U, double _V)
    {
        sI = _I;
        sQ = _Q;
        sU = _U;
        sV = _V;
        sT = 0;
        sSp = 0;
    }

    void set(StokesVector _S)
    {
        sI = _S.I();
        sQ = _S.Q();
        sU = _S.U();
        sV = _S.V();
        sT = _S.T();
        sSp = _S.Sp();
    }

    void addI(double _I)
    {
        sI += _I;
    }

    void addQ(double _Q)
    {
        sQ += _Q;
    }

    void addU(double _U)
    {
        sU += _U;
    }

    void addV(double _V)
    {
        sV += _V;
    }

    void addT(double _T)
    {
        sT += _T;
    }

    void addSp(double _Sp)
    {
        sSp += _Sp;
    }

    void addS(StokesVector _S)
    {
        sI += _S.I();
        sQ += _S.Q();
        sU += _S.U();
        sV += _S.V();
        sT += _S.T();
        sSp += _S.Sp();
    }

    void multI(double _I)
    {
        sI *= _I;
    }

    void multQ(double _Q)
    {
        sQ *= _Q;
    }

    void multU(double _U)
    {
        sU *= _U;
    }

    void multV(double _V)
    {
        sV *= _V;
    }

    void multT(double _T)
    {
        sT *= _T;
    }

    void multSp(double _Sp)
    {
        sSp *= _Sp;
    }

    void multS(double _S)
    {
        sI *= _S;
        sQ *= _S;
        sU *= _S;
        sV *= _S;
    }

    double I() const
    {
        return sI;
    }

    double Q() const
    {
        return sQ;
    }

    double U() const
    {
        return sU;
    }

    double V() const
    {
        return sV;
    }

    double T() const
    {
        return sT;
    }

    double Sp() const
    {
        return sSp;
    }

    void rot(double phi)
    {
        double tQ = sQ, tU = sU;
        double s = sin(2.0 * phi), c = cos(2.0 * phi);
        sQ = tQ * c - tU * s;
        sU = tQ * s + tU * c;
    }

    void rot(double sin_phi, double cos_phi)
    {
        double tQ = sQ, tU = sU;
        sQ = tQ * cos_phi - tU * sin_phi;
        sU = tQ * sin_phi + tU * cos_phi;
    }

    bool isConsistent()
    {
        if(sI * sI < sQ * sQ + sU * sU + sV * sV)
            return false;

        if(sI < 0)
            return false;

        if(sT < 0)
            return false;

        if(sSp < 0)
            return false;

        return true;
    }

    void normalize()
    {
        sQ /= sI;
        sU /= sI;
        sV /= sI;
        sI = 1;
    }

    void clear()
    {
        sI = 0;
        sQ = 0;
        sU = 0;
        sV = 0;
        sT = 0;
        sSp = 0;
    }

    void clearIntensity()
    {
        sI = 0;
        sQ = 0;
        sU = 0;
        sV = 0;
    }

  private:
    double sI, sQ, sU, sV, sT, sSp;
};

class MultiStokesVector
{
  public:
    MultiStokesVector()
    {
        s = 0;
    }

    MultiStokesVector(uint nr_stokes_vector)
    {
        s = new StokesVector[nr_stokes_vector];
    }

    ~MultiStokesVector(void)
    {
        if(s != 0)
            delete[] s;
    }

    void setS(StokesVector _S, uint _i)
    {
        s[_i] = _S;
    }
    StokesVector & S(uint _i)
    {
        return s[_i];
    }

  private:
    StokesVector * s;
};

namespace
{

inline ostream & operator<<(ostream & out, const StokesVector & ex)
{
    out << " I: " << ex.I() << " Q: " << ex.Q();
    out << " U: " << ex.U() << " V: " << ex.V() << " T: " << ex.T() << " ";
    return out;
}

inline StokesVector operator*(const Matrix2D & dM, const StokesVector & v)
{
    return StokesVector(v.I() * dM(0, 0) + v.Q() * dM(0, 1) + v.U() * dM(0, 2) + v.V() * dM(0, 3),
                        v.I() * dM(1, 0) + v.Q() * dM(1, 1) + v.U() * dM(1, 2) + v.V() * dM(1, 3),
                        v.I() * dM(2, 0) + v.Q() * dM(2, 1) + v.U() * dM(2, 2) + v.V() * dM(2, 3),
                        v.I() * dM(3, 0) + v.Q() * dM(3, 1) + v.U() * dM(3, 2) + v.V() * dM(3, 3));
}

inline StokesVector operator*(const StokesVector & v, const StokesVector & u)
{
    return StokesVector(
        v.I() * u.I(), v.Q() * u.Q(), v.U() * u.U(), v.V() * u.V(), v.T() * u.T(), v.Sp() * u.Sp());
}

inline StokesVector operator*(const StokesVector & v, double val)
{
    return StokesVector(v.I() * val, v.Q() * val, v.U() * val, v.V() * val, v.T(), v.Sp());
}

inline StokesVector operator/(const StokesVector & v, double val)
{
    return StokesVector(v.I() / val, v.Q() / val, v.U() / val, v.V() / val, v.T(), v.Sp());
}

}
#endif
